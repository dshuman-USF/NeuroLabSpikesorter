/* Copyright 2005-2020 Kendall F. Morris

    This file is part of the Spiksorter software suite.

    The Spikesorter software suite is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either
    version 3 of the License, or (at your option) any later version.

    The suite is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the suite.  If not, see <https://www.gnu.org/licenses/>.
*/


/* mpk_regions.c */


#define DEBUG 0

#include "nde.h"
#include "region.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <string.h>

#ifdef S_SPLINT_S
typedef unsigned int u_int;
#endif

/*@-namechecks @*/
#include "meschach/matrix.h"
#include "meschach/matrix2.h"
#include "meschach/matlab.h"
/*@=namechecks @*/

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <ctype.h>

#define STATIC static

typedef struct Pctr {
  double a, b;
} Pctr;

typedef struct Plane {
  int a, b;
  float va[SNAPLEN];
  float vb[SNAPLEN];
} Plane;

typedef struct PlaneHdr {
  int count;
  Plane *plane;
} PlaneHdr;

typedef struct Projection {
  float a, b;
} Projection;

typedef struct EDT EDT;

#define CTR 0
#define AMP(n) (fabs(buf[n]-mean))
#define BUFSAMPLES (8 MEG)
#define MEG *(1<<20)
#define ZAP_LEFT 10
#define PREPEAK ZAP_LEFT
#define ZAP_RIGHT 64
static double TAIL_THRESHOLD;
#define TAIL .000000000000001
#define MAX_UNITS 26
#define PXL_SET(x, y) pxl_done[(y*1600+x)/8] |= 1 << x%8
#define PXL_CLR(x, y) pxl_done[(y*1600+x)/8] &= ~(1 << x%8)
#define PXL_TST(x, y) (pxl_done[(y*1600+x)/8] & 1 << x%8)

int debug_region = -1;
int clip_hi, clip_lo;
static int file_size, bytes_read;
static enum {FT_UNKNOWN, FT_SHORT, FT_FLOAT} filetype = FT_SHORT;

static SnapList *global_cc;
static double cvec[SNAPLEN];
static char *pxl_done;
struct EDT {int id; int loc;};

FILE *spikecheck_file;
int *spikecheck;
int spikecheck_count;
static double lp_coeff_array[SNAPLEN+1];
double *lp_coeff = lp_coeff_array + SNAPLEN / 2;

static double **crrn[SNAPLEN+1];
char *file_name;
static char zapped[BUFSAMPLES];
double mean;
double unzapped_sd;
static float unzapped_min, unzapped_max;
float buf[BUFSAMPLES];
static float orig[BUFSAMPLES];
static float *origorig;
int bufsamples;
static int buf_end, buf_start;
static int unzapped_count;
static int zap_count;
static int zap_left = ZAP_LEFT;
static int zap_right = ZAP_RIGHT;
static double sd_sum;
static Snap mean_snap, uzmin_snap, uzmax_snap;
static Snap mean_snap0, uzmin_snap0, uzmax_snap0;
static SnapList hlines, hlines0;
double cov_rep_root[SNAPLEN][SNAPLEN];
static double icov[SNAPLEN][SNAPLEN];
static rfftw_plan fwd_snaplen, rev_snaplen;
//double radius = 7.95829;                         //1.00265 * sqrt(63) SD: median
static double radius = 9;
static int verbose;
static int chunk;
time_t start_time;
static float zero[SNAPLEN];
static int bad_start;
double detection_threshold;
static int *cluster_spikestart;
int starting_sample;
/* SpikeData *spikedata; */
/* float **spike_waveform; */
int spike_count, spikedata_alloc;
extern double cdf[SNAPLEN][CDFSZ];
static double dmap[SNAPLEN][CDFSZ];
static double lastgood[SNAPLEN];
static time_t last_time, now;
static int do_peaklocs;
int debug;

void
check_centers (SnapList *cc)
{
  int n, i;

  if (global_cc == 0)
    global_cc = cc;
  for (n = 0; n < cc->count; n++) {
    for (i = 0; i < SNAPLEN; i++) {
      if (fabs (cc->snap[n]->raw[i]) > 65536) {
        printf ("cluster %d sample %d: %f\n", n, i, cc->snap[n]->raw[i]);
        exit (DIE);
      }
    }
  }
}

STATIC double
map_radius_0 (double rin, int len)
{
  const int lastidx = 142;
  const int snaplenidx = SNAPLEN - 1;
  const double m = CDFMULT;
  
  int lenidx = len - 1;
  int ridx, hi, lo, mid;
  double frac, r;

  ridx = (int)floor (rin * m);
  frac = (ridx > CDFSZ - 2) ? 1.0 : cdf[lenidx][ridx] + (rin * m - ridx) * (cdf[lenidx][ridx+1] - cdf[lenidx][ridx]);
  if (frac >= 1)
    return (rin * m) / lastgood[lenidx] * lastgood[snaplenidx]/m;
  hi = lastidx;
  lo = 0;
  while (hi - lo > 1) {
    mid = (hi + lo) / 2;
    if (frac <= cdf[snaplenidx][mid])
      hi = mid;
    if (frac >= cdf[snaplenidx][mid])
      lo = mid;
  }
  r = ((hi == lo
        ? (double)lo
        : lo + (frac - cdf[snaplenidx][lo]) / (cdf[snaplenidx][hi] - cdf[snaplenidx][lo]))
       / m);
  if (r < 0)
    printf ("%f %d\n", rin, len),
      exit (DIE);
  return r;
}

double
map_radius (double rin, int len)
{
  const double m = CDFMULT;
  
  int lenidx = len - 1;
  int ridx;
  double r;

  ridx = (int)floor (rin * m);

  if (rin > 8)
    r = dmap[lenidx][(int)(8*m)] + (dmap[lenidx][(int)(8*m)] - dmap[lenidx][(int)(7*m)]) * (rin - 8);
  else
    r = dmap[lenidx][ridx] + (rin * m - ridx) * (dmap[lenidx][ridx+1] - dmap[lenidx][ridx]);
  return r;
}

STATIC void
init_cdf (void)
{
  int n, i;
  const double m = CDFMULT;

  for (n = 0; n < SNAPLEN; n++) {
    for (i = 0; i < CDFSZ - 1; i++)
      if (cdf[n][i+1] >= 1)
        break;
    lastgood[n] = (double)i;
  }
  for (n = 0; n < SNAPLEN; n++)
    for (i = 0; i < CDFSZ; i++) {
      if (i / m > 8)
        dmap[n][i] = dmap[n][(int)(8*m)] + (dmap[n][(int)(8*m)] - dmap[n][(int)(7*m)]) * (i/m - 8);
      else
        dmap[n][i] = map_radius_0 (i / m, n+1);
    }
}

STATIC SnapList *
snap_to_list (Snap *s)
{
  static SnapList sl;
  static Snap *snap;
  sl.count = 1;
  sl.snap = &snap;
  snap = s;
  return &sl;
}

void
add_snap (SnapList *sl, Snap *s)
{
  TREALLOC (sl->snap, sl->count+1);
  sl->snap[sl->count++] = s;  
  s->refcnt++;
}

STATIC long filesize (FILE *f)
{
  long bytes;

  fseek (f, 0, SEEK_END);
  bytes = ftell (f);
  rewind (f);
  return bytes;
}

void
shift_snap (float *v, double shift_amt)
{
  double real, imag, w, m, c, s;
  double tdat[SNAPLEN], fdat[SNAPLEN];
  int n;

  m = 2*M_PI*shift_amt / SNAPLEN;
  for (n = 0; n < SNAPLEN; n++)
    tdat[n] = v[n];
  rfftw_one(fwd_snaplen, tdat, fdat);
  for (n = 1; n <= SNAPLEN/2; n++) {
    w = m * n;
    real = fdat[n];
    imag = fdat[SNAPLEN-n];
    fdat[n] = real*(c=cos(w)) - imag*(s=sin(w));
    fdat[SNAPLEN-n] = real*s + imag*c;
  }
  rfftw_one(rev_snaplen, fdat, tdat);
  for (n = 0; n < SNAPLEN; n++)
    v[n] = tdat[n]/SNAPLEN;
}

STATIC void
lpfilter_coeff (void)
{
  fftw_real tdat[SNAPLEN], fdat[SNAPLEN];
  int n;

  rev_snaplen = rfftw_create_plan(SNAPLEN, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
  for (n = 0; n < SNAPLEN/2; n++)
    fdat[n] = 1.0;              /* all real parts 1 except for highest frequency (in SNAPLEN/2) */
  for (n = SNAPLEN / 2; n < SNAPLEN; n++)
    fdat[n] = 0.0;              /* imaginary parts and highest frequency (in SNAPLEN/2) set to 0 */
  for (n = 1; n < 6; n++)
    fdat[SNAPLEN/2-n] = 0.0;
  rfftw_one(rev_snaplen, fdat, tdat);
  for (n = 0; n < SNAPLEN / 2; n++)
    lp_coeff[n] = tdat[n] / SNAPLEN;
  for (n = 1; n < SNAPLEN / 2; n++)
    lp_coeff[n-SNAPLEN/2] = tdat[n+SNAPLEN/2] / SNAPLEN;
  lp_coeff[-SNAPLEN/2] = lp_coeff[SNAPLEN/2] = (tdat[SNAPLEN/2] / SNAPLEN) / 2;
}

STATIC void
whitener (void)
{
  fftw_real tdat[SNAPLEN], fdat[SNAPLEN];
  int n;

  fwd_snaplen = rfftw_create_plan(SNAPLEN, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  for (n = 0; n < SNAPLEN/2; n++)
    tdat[n] = cov_rep_root[n+SNAPLEN/2][SNAPLEN/2];
  for (n = 0; n < SNAPLEN/2; n++)
    tdat[n+SNAPLEN/2] = cov_rep_root[n][SNAPLEN/2];
  rfftw_one(fwd_snaplen, tdat, fdat);
}

static double
lp_filter_fold (float *v, int left, int right)
{
  int n, i;
  double sum = 0;

  for (n = -SNAPLEN / 2; n <= SNAPLEN / 2; n++) {
    if (n < left)
      i = left + (left - n);
    else if (n > right)
      i = right - (n - right);
    else
      i = n;
    sum += (double)v[i] * lp_coeff[n];
  }
  return sum;
}

double
lp_filter (float *v)
{
  int n;
  double sum = 0;

  for (n = -SNAPLEN / 2; n <= SNAPLEN / 2; n++)
    sum += (double)v[n] * lp_coeff[n];
  return sum;
}

double
lp_filter_dbl (double *v)
{
  int n;
  double sum = 0;

  for (n = -SNAPLEN / 2; n <= SNAPLEN / 2; n++)
    sum += (double)v[n] * lp_coeff[n];
  return sum;
}

double
distance (float *a, float *b)   /* debug */
{
  double sum = 0, tmp;
  int n;

  for (n = 0; n < SNAPLEN; n++)
    tmp = a[n] - b[n], sum += tmp*tmp;
  return sqrt(sum);
}

STATIC double
distance_full (Snap *a, Snap *b)
{
  double  sum = 0,  tmp;
  double sumr = 0, tmpr;
  int n;

  for (n = 0; n < SNAPLEN; n++)
    tmp = a->val[n] - b->val[n], sum += tmp*tmp,
      tmpr = a->raw[n] - b->raw[n], sumr += tmpr*tmpr;
#ifdef MIN_W_R
  return MIN (sqrt(sum), sqrt (sumr) / unzapped_sd);
#else
  return sqrt(sum);
#endif
}

STATIC void write_edt (char *file_name, SnapList *cc);

static void
write_status (char *file_name, SnapList *cc, char *msg)
{
  FILE *f;

  if (msg) {
    printf ("\n%s\n", msg);
    (f = fopen (change_filetype (file_name, ".msg"), "w")) || DIE;
    fprintf (f, "%s\n", msg);
    fclose (f);
  }
  write_edt (file_name, cc);
  (f = fopen (change_filetype (file_name, ".status"), "a")) || DIE;
  putc ('S', f);
  fclose (f);
  exit (0);
}

STATIC void find_mean (void)
{
  int n;
  double sum = 0;
  float max, min;

  max = min = buf[0];
  for (n = 0; n < bufsamples; n++) {
    sum += buf[n];
    if (buf[n] > max) max = buf[n];
    if (buf[n] < min) min = buf[n];
  }

  if (max - min < 100)
    write_status (file_name, 0, "flat signal - no spikes found");

  mean = sum / bufsamples;
  TMALLOC (mean_snap.val, SNAPLEN);
  TMALLOC (mean_snap0.val, SNAPLEN);
  TMALLOC (mean_snap.raw, SNAPLEN);
  TMALLOC (mean_snap0.raw, SNAPLEN);
  for (n = 0; n < SNAPLEN; n++)
    mean_snap0.raw[n] = mean_snap.raw[n] = mean_snap0.val[n] = mean_snap.val[n] = mean;
  printf ("max: %.0f, min: %.0f, mean: %18.16f \n", max, min, mean);
}

STATIC int spikedata_cmp (const SpikeData *a, const SpikeData *b)
{
  return a->sample < b->sample ? -1 : a->sample > b->sample; /* ascending */
}

STATIC double find_unzapped_sd (void)
{
  int count = 0, n;
  double sum = 0;
  double tmp;

  unzapped_max = unzapped_min = mean;
  for (n = 0; n < bufsamples; n++)
    if (!zapped[n]) {
      tmp = buf[n] - mean, sum += tmp*tmp, count++;
      if (buf[n] < unzapped_min)
        unzapped_min = buf[n];
      if (buf[n] > unzapped_max)
        unzapped_max = buf[n];
    }
  //  printf ("unzapped min: %f, unzapped max: %f\n", unzapped_min, unzapped_max);
  unzapped_count = count;

  sd_sum = sum;
  unzapped_sd = (double)(sqrt (sum / (count-1)));
  return unzapped_sd;
}

STATIC double find_overall_sd (void)
{
  int n;
  double sum = 0;
  double tmp;

  for (n = 0; n < bufsamples; n++)
    tmp = buf[n] - mean, sum += tmp*tmp;

  return (double)(sqrt (sum / (bufsamples-1)));
}

STATIC void
find_unzapped_minmax (void)
{
  int n;

  unzapped_max = unzapped_min = mean;
  for (n = 0; n < bufsamples; n++)
    if (!zapped[n]) {
      if (buf[n] < unzapped_min)
        unzapped_min = buf[n];
      if (buf[n] > unzapped_max)
        unzapped_max = buf[n];
    }
  if (verbose) {
    printf ("unzapped_min %f from mean\n", mean - unzapped_min);
    printf ("unzapped_max %f from mean\n", unzapped_max - mean);
  }
  TMALLOC (uzmin_snap.val, SNAPLEN);
  TMALLOC (uzmin_snap0.val, SNAPLEN);
  for (n = 0; n < SNAPLEN; n++)
    uzmin_snap0.val[n] = uzmin_snap.val[n] = unzapped_min;
  TMALLOC (uzmax_snap.val, SNAPLEN);
  TMALLOC (uzmax_snap0.val, SNAPLEN);
  for (n = 0; n < SNAPLEN; n++)
    uzmax_snap0.val[n] = uzmax_snap.val[n] = unzapped_max;
  add_snap (&hlines, &mean_snap);
  add_snap (&hlines, &uzmin_snap);
  add_snap (&hlines, &uzmax_snap);
  add_snap (&hlines0, &mean_snap0);
  add_snap (&hlines0, &uzmin_snap0);
  add_snap (&hlines0, &uzmax_snap0);
  if (verbose)
    printf ("unzapped min: %f, unzapped max: %f\n", unzapped_min, unzapped_max);
}

static void
filter_buf (float *v, int count)
{
  static float *t;
  static int alloc;
  int n;

  if (count > alloc)
    TREALLOC (t, alloc = count);
  for (n = 0; n < SNAPLEN/2; n++)
    t[n] = lp_filter_fold (v + n, -n, SNAPLEN/2);
  for (n = SNAPLEN/2; n  < count - SNAPLEN/2; n++)
    t[n] = lp_filter (v + n);
  for (n = count - SNAPLEN/2; n < count; n++)
    t[n] = lp_filter_fold (v + n, -SNAPLEN/2, count - n - 1);
  TMEMCPY (v, t, count);
}

static inline int
mid_cond (int prev, int a, int b, int next, int first)
{
  return first > prev && first <= a && first + SNAPLEN > b && first + SNAPLEN <= next;
}

static inline int
mult (int prev, int a, int b, int next)
{
  int min_first, max_first;

  if (0) {
    int n;
    (prev < a && a <= b && b < next) || DIE;
    
    prev == -1 || zapped[prev] || DIE;
    zapped[next] || DIE;
    if (1) {
      zapped[prev+1] == 0 || DIE;
      zapped[next-1] == 0 || DIE;
    }
    else
      for (n = prev + 1; n < next; n++) zapped[n] == 0 || DIE;
    next - prev > SNAPLEN || DIE;
    b - a < SNAPLEN || DIE;
  }

  min_first = MAX (prev + 1, b - SNAPLEN + 1);
  max_first = MIN (a, next - SNAPLEN);

  if (0) {
    int n;
    max_first >= min_first || DIE;
    if (0) {
      for (n = min_first; n <= max_first; n++)
        mid_cond (prev, a, b, next, n) || DIE;
    }
    else {
      mid_cond (prev, a, b, next, min_first) || DIE;
      mid_cond (prev, a, b, next, max_first) || DIE;
    }
    mid_cond (prev, a, b, next, min_first - 1) && DIE;
    mid_cond (prev, a, b, next, max_first + 1) && DIE;
  }
  return max_first - min_first + 1;
}

STATIC void
get_cvec (void)
{
  int count[SNAPLEN], a, b, d, unzapped_count;
  double sum[SNAPLEN];
  double unzapped_mean;
  int prev = -1, next;


  if (0)
  {
    float *tmp;
    printf ("%s line %d: filtering noise\n", __FILE__, __LINE__);

    filter_buf (buf, bufsamples);
    filter_buf (buf, bufsamples);
    filter_buf (buf, bufsamples);
    if (0) {
      TMALLOC (tmp, bufsamples);
      for (a = SNAPLEN/2; a + SNAPLEN/2 < bufsamples; a++)
        tmp[a] = lp_filter (buf + a);
      for (a = SNAPLEN/2; a + SNAPLEN/2 < bufsamples; a++)
        buf[a] = tmp[a];
      free (tmp);
    }

  }

  memset (count, 0, sizeof count);
  memset (sum, 0, sizeof sum);
  unzapped_mean = 0;
  unzapped_count = 0;
  for (a = 0; a < bufsamples; a++) {
    if (zapped[a])
      continue;
    unzapped_mean += buf[a];
    unzapped_count++;
  }
  unzapped_mean /= unzapped_count;

  for (a = 0; a < bufsamples; a++) {
    if (zapped[a]) {
      prev = a;
      continue;
    }
    for (next = a; next < bufsamples; next++)
      if (zapped[next])
        break;
    if (next - prev <= SNAPLEN)
      continue;
    {
      static time_t now, last_time;
      now = time (0);
      if (last_time == 0)
        last_time = now + 1;
      if (now > last_time) {
        last_time = now;
        printf ("get_cvec: sample %10d of %10d  \r", a, bufsamples);
        fflush (stdout);
      }
    }
    for (d = 0; d < SNAPLEN && (b = a + d) < next; d++) {
      int last_count, m;
      m = mult (prev, a, b, next);
      sum[d] += (buf[a]-unzapped_mean) * (buf[b]-unzapped_mean) * m;
      last_count = count[d];
      (count[d] += m) > last_count || DIE;;
    }
  }

  for (d = 0; d < SNAPLEN; d++)
    cvec[d] = sum[d] / (count[d]-1);

  printf ("unzapped_mean: %f\n", unzapped_mean);
}

double (*get_vc(void))[SNAPLEN];

STATIC double **
get_wmat_n (int len)
{
  int a, b, n;
  static MAT *cov, *tmp, *tmp2, *micov;
  static MAT *eigvec, *evdiag;
  static VEC *eigval;
  double **crr;
  double (*vc)[64];

  vc = get_vc ();

  m_resize_vars (len, len, &cov, &eigvec, &evdiag, 0);

  for (a = 0; a < len; a++)
    for (b = a; b < len; b++)
      cov->me[a][b] = cov->me[b][a] = vc[a][b];

  micov = m_inverse (cov, 0);
  for (a = 0; a < len; a++)
    for (b = 0; b < len; b++)
      icov[a][b] = micov->me[a][b];
  eigval = symmeig (cov, eigvec, eigval);
  for (n = 0; n < len; n++)
    if (eigval->ve[n] < 0) {
      note ("covariance matrix is not positive definite");
      exit (DIE);
    }

  (void)m_zero (evdiag);
  for (n = 0; n < len; n++)
    evdiag->me[n][n] = 1/sqrt(eigval->ve[n]);
  tmp = m_mlt (eigvec, evdiag, tmp);
  tmp2 = mmtr_mlt (tmp, eigvec, tmp2);
  TMALLOC (crr, len);
  for (a = 0; a < len; a++) {
    TMALLOC (crr[a], len);
    for (b = 0; b < len; b++)
      crr[a][b] = tmp2->me[a][b];
  }
  return crr;
}

STATIC void
get_all_wmat (void)
{
  int n, a, b;

  printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
  printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
  for (n = 64; n >= 1; n--)
    crrn[n] = get_wmat_n (n);
  printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
  if (0)
    for (a = 0; a < SNAPLEN; a++)
      for (b = 0; b < SNAPLEN; b++)
        crrn[64][a][b] == cov_rep_root[a][b] || DIE;
}

STATIC int
read_vc (double (*vc)[SNAPLEN])
{
  FILE *f;
  if (!((f = fopen (change_filetype (file_name, ".vc"), "rb"))
        && fread (vc, sizeof *vc, SNAPLEN, f) == SNAPLEN)) {
    return 0;
  }
  fclose (f);
  printf ("read vc from %s\n", change_filetype (file_name, ".vc"));
  return 1;
}

STATIC void
write_vc (double (*vc)[SNAPLEN])
{
  FILE *f;
  (f = fopen (change_filetype (file_name, ".vc"), "wb")) || DIE;
  fwrite (vc, sizeof *vc, SNAPLEN, f) == SNAPLEN || DIE;
  fclose (f);
}

double (*
get_vc(void))[SNAPLEN]
{
  int n, a, b, i, end, count;
  char *snap_zapped;
  float *bufp;
  long double sum[SNAPLEN], vcsum[SNAPLEN][SNAPLEN];
  double mean[SNAPLEN];
  int vccount[SNAPLEN][SNAPLEN];
  static double vc[SNAPLEN][SNAPLEN];
  static int calc_done;

  if (calc_done)
    return vc;
  if (read_vc (vc)) {
    calc_done = 1;
    return vc;
  }

  TCALLOC (snap_zapped, bufsamples);
  count = 0;
  for (n = 0; n < bufsamples; n++) {
    end = n + SNAPLEN;
    if (end > bufsamples) {
      snap_zapped[n] = 1;
      count++;
      continue;
    }
    for (a = n; a < end; a++)
      if (zapped[a]) {
        snap_zapped[n] = 1;
        count++;
        break;
      }
  }
  if (count == bufsamples) {
    printf ("spikes are everywhere\n");
    exit (1);
  }
  count = 0;
  memset (sum, 0, sizeof sum);
  for (n = 0; n < bufsamples; n++) {
    {
      static time_t now, last_time;
      now = time (0);
      if (last_time == 0)
        last_time = now + 1;
      if (now > last_time) {
        last_time = now;
        printf ("get_vc sum: sample %10d of %10d  \r", n, bufsamples); fflush (stdout);
      }
    }
    if (snap_zapped[n])
      continue;
    count++;
    bufp = buf + n;
    for (i = 0; i < SNAPLEN; i++)
      sum[i] += bufp[i];
  }
  for (i = 0; i < SNAPLEN; i++)
    mean[i] = sum[i] / count;

  memset (vccount, 0, sizeof vccount);
  memset (vcsum,   0, sizeof vcsum);

  for (n = 0; n < bufsamples; n++) {
    {
      static time_t now, last_time;
      now = time (0);
      if (last_time == 0)
        last_time = now + 1;
      if (now > last_time) {
        last_time = now;
        printf ("get_vc vcsum: sample %10d of %10d  \r", n, bufsamples); fflush (stdout);
      }
    }
    if (snap_zapped[n])
      continue;
    bufp = buf + n;
    for (a = 0; a < SNAPLEN; a++)
      for (b = a; b < SNAPLEN; b++) {
        vcsum[a][b] += (bufp[a] - mean[a]) * (bufp[b] - mean[b]);
        vccount[a][b]++;
      }
  }
  printf ("%50s\r", ""); fflush(stdout);
  for (a = 0; a < SNAPLEN; a++)
    for (b = a; b < SNAPLEN; b++)
      vc[a][b] = vcsum[a][b] / vccount[a][b];
  free (snap_zapped);
  calc_done = 1;
  write_vc (vc);
  return vc;
}

STATIC void find_unzapped_cov (void)
{
  int a, b, n;
  MAT *cov, *tmp, *tmp2, *micov;
  MAT *eigvec, *evdiag;
  VEC *eigval = 0;
  double (*vc)[SNAPLEN];

  cov = m_get (SNAPLEN, SNAPLEN);

  vc = get_vc ();
  for (a = 0; a < SNAPLEN; a++)
    for (b = a; b < SNAPLEN; b++)
      cov->me[a][b] = cov->me[b][a] = vc[a][b];

  micov = m_inverse (cov, 0);
  for (a = 0; a < SNAPLEN; a++)
    for (b = 0; b < SNAPLEN; b++)
      icov[a][b] = micov->me[a][b];
  eigvec = m_get (SNAPLEN, SNAPLEN);
  evdiag = m_get (SNAPLEN, SNAPLEN);
  eigval = symmeig (cov, eigvec, VNULL);

  for (n = 0; n < SNAPLEN; n++)
    if (eigval->ve[n] < 0) {
      printf ("covariance matrix is not positive definite\n");
      exit (1);
    }

  (void)m_zero (evdiag);
  for (n = 0; n < SNAPLEN; n++)
    evdiag->me[n][n] = 1/sqrt(eigval->ve[n]);
  tmp = m_mlt (eigvec, evdiag, MNULL);
  tmp2 = mmtr_mlt (tmp, eigvec, MNULL);
  for (a = 0; a < SNAPLEN; a++)
    for (b = 0; b < SNAPLEN; b++)
      cov_rep_root[a][b] = tmp2->me[a][b];

  if (0)
  {
    static double max, min;
    for (a = 0; a < SNAPLEN; a++)
      for (b = 0; b < SNAPLEN; b++) {
        if (cov_rep_root[a][b] < min)
          min = cov_rep_root[a][b];
        if (cov_rep_root[a][b] > max)
          max = cov_rep_root[a][b];
      }
    printf ("min: %.17f, max: %.17f\n", min, max);
  }
}

STATIC void find_white_cov (void)
{
  int a, b, n;
  MAT *cov, *tmp, *tmp2, *micov;
  MAT *eigvec, *evdiag;
  VEC *eigval = 0;
  double cov_rep_root_tmp[SNAPLEN][SNAPLEN];

  cov = m_get (SNAPLEN, SNAPLEN);
  for (a = 0; a < SNAPLEN; a++)
    for (b = 0; b < SNAPLEN; b++)
      cov->me[a][b] = a == b ? unzapped_sd * unzapped_sd : 0;
  micov = m_inverse (cov, 0);
  for (a = 0; a < SNAPLEN; a++)
    for (b = 0; b < SNAPLEN; b++)
      icov[a][b] = micov->me[a][b];
  eigvec = m_get (SNAPLEN, SNAPLEN);
  evdiag = m_get (SNAPLEN, SNAPLEN);
  eigval = symmeig (cov, eigvec, VNULL);
  (void)m_zero (evdiag);
  for (n = 0; n < SNAPLEN; n++)
    evdiag->me[n][n] = 1/sqrt(eigval->ve[n]);
  tmp = m_mlt (eigvec, evdiag, MNULL);
  tmp2 = mmtr_mlt (tmp, eigvec, MNULL);
  if (0) {
    for (a = 0; a < SNAPLEN; a++)
      for (b = 0; b < SNAPLEN; b++)
        cov_rep_root_tmp[a][b] = tmp2->me[a][b];
    printf ("cov_rep_root_tmp[0][0] = %.17f, cov_rep_root_tmp[0][1] = %.17f\n",
            cov_rep_root_tmp[0][0], cov_rep_root_tmp[0][1]);
  }
}

STATIC int
zap (int n)
{
  int samples_zapped = 0, start, end, z;
  double tmp;

  zap_count++;
  start = n - zap_left;
  end = n + zap_right + 1;
  if (start < 0)
    bad_start = 1,
      start = 0;
  if (end > bufsamples)
    end = bufsamples;
  for (z = start; z < end; z++)
    if (!zapped[z]) {
      tmp = buf[z] - mean;
      sd_sum -= tmp*tmp;
      zapped[z] = 1, samples_zapped++;
    }
  return samples_zapped;
}

STATIC size_t (*fileread)(void *, size_t, size_t, void *);

STATIC inline size_t
std_read (void *buf, size_t size, size_t count, void *p)
{
  size_t rcnt;
  unsigned short v = 1;

  rcnt = fread (buf, size, count, p);
  if (*(char *)&v == 0) {
    /* big endian */
    size_t scnt = size * rcnt / 2;
    unsigned short *sbuf = buf;
    int n;

    for (n = 0; n < scnt; n++)
      sbuf[n] = ((sbuf[n] >> 8) & 0xff) | ((sbuf[n] & 0xff) << 8);
  }
  return rcnt;
}

STATIC int
fillbuf2 (void *f)
{
  int n, i, readcnt;
  static long residue;

  if (f == 0) {
    chunk = 0;
    residue = 0;
    starting_sample = 0;
    bufsamples = 0;
    if (origorig)
      free (origorig);
  return 1;
  }

  chunk++;
  buf_start = buf_end + BEFORE;
  buf_end++;
  for (residue = n = 0; buf_end + n < bufsamples; n++)
    buf[residue++] = buf[buf_end+n];
  starting_sample += bufsamples - residue;
  buf_start -= bufsamples - residue;

  if (filetype == FT_FLOAT)
    readcnt = (int)fileread (buf+residue, sizeof *buf, (size_t)(BUFSAMPLES-residue), f);
  else
    readcnt = (int)fileread (buf+residue, sizeof *buf, (size_t)(BUFSAMPLES-residue)/2, f);

  if (filetype == FT_UNKNOWN) {
    for (n = 0; n < readcnt; n++) {
      if (buf[n] == 0)
        continue;
      if (buf[n] > 65536 || buf[n] < -65536 /*|| !finitef (buf[n])*/)
        break;
    }
    if (n < readcnt)
      filetype = FT_SHORT;
    else
      filetype = FT_FLOAT;
  }
  if (filetype == FT_SHORT) {
    short *sbuf = (short *) (buf+residue);
    float *fbuf =            buf+residue;
    if (chunk == 1 && verbose)
      printf ("data seem to be 16 bit integers\n");
    for (i = readcnt*2-1; i >= 0; i--)
      fbuf[i] = (float)sbuf[i];
    bufsamples = readcnt * 2 + residue; 
  }
  else {
    if (chunk == 1)
      printf ("data seem to be single precision floating point\n");
    if (readcnt == BUFSAMPLES/2)
      readcnt += fileread (buf+BUFSAMPLES/2, sizeof *buf, (size_t)BUFSAMPLES/2, f);
    bufsamples = readcnt + residue;
  }
  if (verbose)
    printf ("%ld samples read, residue was %ld\n", bufsamples-residue, residue);
  unzapped_count = bufsamples;
  TMEMCPY (orig, buf, bufsamples);
  return readcnt;
}

STATIC int
fillbuf (char *file_name)
{
  static void *f;

  if (file_name == 0) {
    rewind (f);
    return fillbuf2 (0);
    buf_start = 0;
  }
  if (f)
    bytes_read = ftell (f);
  printf ("%79s\r  chunk %d (%.0f%%): reading\r", "", chunk, bytes_read * 100.0 / file_size); fflush (stdout);
  if (f)
    return (fillbuf2 (f));
  if ((f = fopen (file_name, "rb"))) {
    file_size = filesize (f);
    fileread = std_read;
    return fillbuf2 (f);
  }

  if ((f = dmx_open (file_name, (size_t)BUFSAMPLES))) {
    fileread = dmx_read;
    return fillbuf2 (f);
  }

  exit (DIE);
  return 0;
}

STATIC void
free_spline (Snap *s)
{
  if (s->spline)
    ss_spline_free (s->spline);
  if (s->acc)
    gsl_interp_accel_free(s->acc);
  s->spline = 0;
  s->acc = 0;
}

STATIC void free_snap (Snap **sp)
{
  Snap *s = *sp;

  if (s == 0)
    return;
  if (s->refcnt <= 0)
    exit (DIE);
  if (--s->refcnt == 0) {
    free_spline (s);
    if (s->buf)
      free (s->buf);
    if (s->free_val)
      free (s->val);
    if (s->free_raw)
      free (s->raw);
    free (s);
    *sp = 0;
  }
}

void
free_snaplist (SnapList **slp)
{
  int n;
  SnapList *sl = *slp;
  if (sl == 0)
    return;
  for (n = 0; n < sl->count; n++)
    free_snap (&sl->snap[n]);
  free (sl->snap);
  free (sl);
  *slp = 0;
}

void
whiten_snap (float *v)          /* debug */
{
  int i, j;
  float tmp[SNAPLEN];
  double sum;
  float fmean = mean;

  for (i = 0; i < SNAPLEN; i++) {
    sum = 0;
    for (j = 0; j < SNAPLEN; j++)
      sum += (v[j]-fmean) * cov_rep_root[j][i];
    tmp[i] = (float) sum;
  }

  memcpy (v, tmp, SNAPLEN * sizeof *tmp);
}

STATIC void add_snap_alloc (SnapList *sl, float *v)
{
  Snap *s;

  TCALLOC (s, 1);
  s->count = 1;
  s->radius = radius;
  TMALLOC (s->val, SNAPLEN);
  TMALLOC (s->raw, SNAPLEN);
  memcpy (s->val, v, sizeof (SnapVals));
  memcpy (s->raw, v, sizeof (SnapVals));
  whiten_snap (s->val);
  add_snap (sl, s);
}

STATIC SnapList *find_noise_snaps (void)
{
  int n, count;
  SnapList *sl;

  TCALLOC (sl, 1);
  for (count = n = 0; n < bufsamples; n++) {
    if (!zapped[n]) {
      count++;
      if (count == SNAPLEN) {
        add_snap_alloc (sl, buf +  n - (SNAPLEN - 1));
        count = 0;
      }
    }
    else
      count = 0;
  }
  if (verbose)
    printf ("%d noise snapshots\n", sl->count);
  return sl;
}

void
whiten_samples (float *v, int count)
{
  int i, j;
  float tmp[count];
  double sum;
  float fmean = mean;

  for (i = 0; i < count; i++) {
    sum = 0;
    for (j = 0; j < count; j++)
      sum += (v[j]-fmean) * crrn[count][j][i];
    tmp[i] = (float) sum;
  }

  memcpy (v, tmp, count * sizeof *tmp);
}

#define OVER(n)  (v[n] > unzapped_max || v[n] < unzapped_min)
#define UNDER(n) (v[n] <= unzapped_max && v[n] >= unzapped_min)

//#define PARTITION_THRESHOLD 2.25
#define PARTITION_THRESHOLD 6

typedef struct PartitionCenter PartitionCenter;
struct PartitionCenter
{
  SnapValArray prev_v;
  SnapValArray best_v;
  double best_d;
};

STATIC double
distance_partsnap_debug (float *a, float *b, int count)
{
  double sum = 0, tmp;
  int n;

  for (n = 0; n < count; n++)
    tmp = a[n] - b[n], sum += tmp*tmp;
  tmp = sqrt(sum);
  printf ("map_radius (%18.16f, %d) = %18.16f\n", tmp, count, map_radius (tmp, count));
  return map_radius (tmp, count);
}

double
mdistance_at_shift_partsnap_debug (double shift, void *p)
{
  SnapList *sl = p;
  float v[SNAPLEN];
  float *fixed = sl->snap[1]->val;
  int ishift = (int)sl->snap[0]->shift;

  printf ("shift %f\n", shift);
  sl->count == 2 || DIE;
  memcpy (v, sl->snap[0]->raw, sizeof (SnapVals));
  shift_snap (v, shift);
  whiten_snap (v);
  printf ("count %d\n", SNAPLEN/2 + (ishift < 0) * abs (ishift));
  return distance_partsnap_debug (v, fixed, SNAPLEN/2 + (ishift < 0) * abs (ishift));
}

#undef OVER
#undef UNDER

STATIC void zap_to_threshold (double threshold)
{
  int n;

  printf ("%79s\r  chunk %d (%.0f%%): zapping\r", "", chunk, bytes_read * 100.0 / file_size); fflush (stdout);
  TMEMSET (zapped, 0, bufsamples);
  bad_start = zap_count = 0;
  for (n = 0; n < bufsamples; n++)
    if (fabs(buf[n]-mean) >= threshold)
      unzapped_count -= zap (n);
  if (verbose)
    printf ("bufsamples: %d, unzapped_count: %d\n", bufsamples, unzapped_count);
}

const int do_decompose_pdf = 0;

STATIC void
gen_cluster_spikepeak (SnapList *cc)
{
  int n, i, i_at_max = 0;
  double max, v;

  for (n = 0; n < cc->count; n++) {
    for (max = 0, i = 0; i < SNAPLEN; i++)
      if ((v = cc->snap[n]->raw[i] - mean) > max) {
        max = v;
        i_at_max = i;
      }
    cluster_spikestart[n] = i_at_max;
  }
}

STATIC void
gen_cluster_spikestart (SnapList *cc)
{
  int n, i;
  double spikestart_threshold;

  TREALLOC (cluster_spikestart, cc->count + 1);

  for (n = 0; n < cc->count; n++)
    cluster_spikestart[n] = PRESAMPLES;
  return;
  
  if (do_peaklocs) {
    gen_cluster_spikepeak (cc);
    return;
  }
    
  for (n = 0; n < cc->count; n++) {
    if (cc->snap[n]->nocenter)
      continue;
    //    spikestart_threshold = (.2 + .8/sqrt(cc->snap[n]->count)) * detection_threshold;
    spikestart_threshold = detection_threshold;
    if (verbose)
      printf ("cluster %d count %d detection_threshold %f spikestart_threshold %f\n", 
              n + 1, cc->snap[n]->count, detection_threshold, spikestart_threshold);
    for (i = 0; i < SNAPLEN; i++) {
      if (verbose)
        printf ("    sample %d val %f\n", i, fabs(cc->snap[n]->raw[i] - mean));
      if (fabs(cc->snap[n]->raw[i] - mean) >= spikestart_threshold) {
        cluster_spikestart[n] = i;
        break;
      }
    }
  }
}

char *
change_filetype (char *file_name, char *suffix)
{
  char *dot, *p;
  static char *new_file;

  TREALLOC (new_file, strlen (file_name) + strlen (suffix) + 1);
  strcpy (new_file, file_name);
  if ((dot = strrchr (new_file, '.')) == 0)
    strcat (new_file, suffix);
  else
    strcpy (dot, suffix);
  if ((p = strrchr (new_file, '/')) != 0
      || (p = strrchr (new_file, '\\')) != 0)
    return p + 1;
  else
    return new_file;
}

STATIC void
write_edt (char *file_name, SnapList *cc)
{
  FILE *f, *g;
  int n, ibuf[2];

  printf ("%79s\r  write edt\r", ""); fflush (stdout);
  qsort (spikedata, (size_t)spike_count, sizeof *spikedata, (int(*)(const void*,const void*))spikedata_cmp);

  (f = fopen (change_filetype (file_name, ".edt"), "w")) || DIE;
  (g = fopen (change_filetype (file_name, ".pos"), "wb")) || DIE;
  if (!do_peaklocs) {
    fprintf (f, "   33   3333333\n");
    fprintf (f, "   33   3333333\n");
  }
  ibuf[1] = 0;
  for (n = 0; n < spike_count; n++) {
    if (spikedata[n].partition == -1)
      continue;

    if (!cc->snap[spikedata[n].partition]->noise) {
      fwrite (&spikedata[n].sample, sizeof (double), 1, g) == 1 || DIE;
      ibuf[0] = spikedata[n].partition+1;
      fwrite (ibuf, sizeof ibuf, 1, g) == 1 || DIE;
    }
    if (! spikedata[n].leave_out) {
      if (!cc->snap[spikedata[n].partition]->noise) {
        if (do_peaklocs)
          fprintf (f, "%d %.0f\n", spikedata[n].partition+1, spikedata[n].sample);
        else
          fprintf (f, "%5d%10d\n", spikedata[n].partition+1, (int)floor(spikedata[n].sample/(SF/10000)+.5));
      }
    }
    else {
      if (do_peaklocs)
        fprintf (f, "%d %.0f\n", 0, spikedata[n].sample);
      else
        fprintf (f, "%5d%10d\n", 0, (int)floor(spikedata[n].sample/(SF/10000)+.5));
    }
  }
  fclose (f);
  fclose (g);
}

typedef struct VSum VSum;

struct VSum {
  double v[SNAPLEN];
  float mean[SNAPLEN];
  double cov[SNAPLEN][SNAPLEN];
  double *ve;
  int hist[20];
  int count;
  int n;
  int prin;
};
#define ROUND(x) floor ((x) * multiplier) / multiplier;

static inline int
ltr (int n)
{
  if (n >= 62)
    return '-';
  if (n < 10)
    return '0' + n;
  if (n < 36)
    return 'A' + n - 10;
  return 'a' + n - 36;
}

static inline int
nbr (int c)
{
  if (c >= '0' && c <= '9')
    return c - '0';
  if (c >= 'A' && c <= 'Z')
    return 10 + c - 'A';
  if (c >= 'A' + 128 && c <= 'Z' + 128)
    return 36 + c - 128 - 'A';
  return -1;
}

static inline void
maybe_flip (char *shw, size_t size, int c)
{
  int n;
  if ((n = nbr (c)) < 0 || (size_t)n >= size)
    return;
  shw[n] ^= 1;
}

#define VWRITE(v) count += sizeof (v) / sizeof (short), fwrite (&(v), sizeof (v), 1, f) == 1 || DIE
#define ORWRITE(v,m) count += sizeof buf / sizeof (short), buf = (unsigned short)(v | m), fwrite (&(buf), sizeof (buf), 1, f) == 1 || DIE
#define CWRITE(v) count += sizeof buf / sizeof (short), buf = (unsigned short)v, fwrite (&(buf), sizeof (buf), 1, f) == 1 || DIE

STATIC int count;
STATIC void
scatter (Projection *p, /*@unused@*/ Pctr *pctr, FILE *winp, int setlimits, FILE *f)
{
#define W 1600
#define H 1200
  double x, y;
  static double left, right, top, bottom;
  static int do_adjust;
  unsigned short xi, yi;

  if (p == 0) {
    top = right = -HUGE_VAL;
    bottom = left = HUGE_VAL;
    do_adjust = 1;
    return;
  }

  x = (double)p->a;
  y = (double)p->b;

  if (setlimits == 1) {
    if (x < left)
      left = x;
    if (x > right)
      right = x;
    if (y < bottom)
      bottom = y;
    if (y > top)
      top = y;
    return;
  }

  if (do_adjust) {
    double hpp, vpp, hd, vd, hc, vc;
    
    hd = right - left;
    hc = left + hd / 2;
    vd = top - bottom;
    vc = bottom + vd / 2;
    hpp = hd / W;
    vpp = vd / H;
    if (hpp > vpp) {
      left -= .03 * hd;
      right += .03 * hd;
      hd = right - left;
      bottom = vc - (hd / W * H) / 2;
      top = vc + (hd / W * H) / 2;
    }
    else {
      bottom -= .03 * vd;
      top += .03 * vd;
      vd = top - bottom;
      left = hc - (vd / H * W) / 2;
      right = hc + (vd / H * W) / 2;
    }
    do_adjust = 0;
  }
    
  if (x >= left && x <= right && y >= bottom && y <= top) {
    xi = (unsigned short)floor ((x - left) / (right - left) * (W-1) + .5);
    yi = (unsigned short)floor ((top - y) / (top - bottom) * (H-1) + .5);
    if (PXL_TST(xi, yi) == 0) {
      fprintf (winp, "%d %d d\n", (int)xi, (int)yi);
      VWRITE (xi);
      VWRITE (yi);
      PXL_SET (xi, yi);
    }
  }
}

STATIC Projection *
spike_projection (float *v, Plane *pln)
{
  static Projection proj;
  int i;
  double x, y;

  for (x = y = 0, i = 0; i < SNAPLEN; i++) {
    x += v[i] * pln->va[i];
    y += v[i] * pln->vb[i];
  }
  proj.a = x;
  proj.b = y;
  return &proj;
}

STATIC void
scatterplot (FILE *winp, SnapList *cc, Projection **noise_projections, int ncnt, Pctr **pctr_list, PlaneHdr *ph)
{
  static Pctr nctr;
  static Projection *proj;
  unsigned cc_color_list[] = {0xff000000U, 0xffdb0000U, 0x92ff0000U, 0x24b6ff00U,
                              0xdb00ff00U, 0xff00ff00U, 0x0000ff00U, 0x00ffff00U, 0x92490000U};
  int color_count = sizeof cc_color_list / sizeof cc_color_list[0];
  unsigned *cc_color;
  int i, pln, n, used;
  short a, b, ccn;
  unsigned short buf;
  char text[11];
  FILE *f;
  int limit, max_count_per_cluster, hit_limit = 0;

  if (ph->count == 0)
    return;
    
  (f = fopen (change_filetype (file_name, ".spl"), "wb")) || DIE;
  TMALLOC (pxl_done, 1600*1200/8);
  TMALLOC (cc_color, cc->count);
  for (i = n = 0; n < cc->count; n++) {
    if (!cc->snap[n]->noise && cc->snap[n]->count > 0) {
      cc_color[n] = cc_color_list[i];
      i = (i + 1) % color_count;
    }
    else cc_color[n] = 0;
  }
  for (used = 0, ccn = 0; ccn < cc->count; ccn++)
    if (!cc->snap[ccn]->noise && cc->snap[ccn]->count > 0)
      used++;
  used++;
  VWRITE (ph->count);
  VWRITE (used);
  max_count_per_cluster = 1000000000 / 2 / ph->count / used;
  for (used = 0, ccn = 0; ccn < cc->count; ccn++)
    if (!cc->snap[ccn]->noise && cc->snap[ccn]->count > 0) {
      float d = distance (cc->snap[ccn]->val, zero);
      VWRITE (d);
    }
      
  for (pln = 0; pln < ph->count; pln++) {
    printf ("%79s\r  scatterplot %d of %d\r", "", pln + 1, ph->count); fflush (stdout);
    scatter (0, 0, 0, 0, f);
    for (n = 0; n < spike_count; n++) {
      if (spikedata[n].partition == -1
          || spikedata[n].partition == cc->count - 1
          || spike_waveform[n] == 0
          || spikedata[n].nocenter
          )
        continue;
      if ((now = time (0)) != last_time) {
        printf ("%79s\r  scatterplot %d of %d: scale spike %d of %d\r",
                "", pln + 1, ph->count, n + 1, spike_count); fflush (stdout);
        last_time = now;
      }
      proj = spike_projection (spike_waveform[n], &ph->plane[pln]);
      scatter (proj, &pctr_list[spikedata[n].partition][pln], winp, 1, f);
    }
    for (n = 0; n < ncnt; n++)
      scatter (&noise_projections[n][pln], &nctr, winp, 1, f);
    start_pdf_page (winp);
    fprintf (winp, "6 setlinewidth 1 setlinecap\n");
    a = ph->plane[pln].a;
    b = ph->plane[pln].b;
    count = 0;
    VWRITE (a);
    VWRITE (b);
    {
      float d = 0;
      if (b != -1)
        d = distance (cc->snap[a]->val, cc->snap[b]->val);
      VWRITE (d);
    }
    setcolor (winp, cc_color[a]);
    snprintf (text, 11, "%d", a + 1);
    pdftext (winp, 10, 35, text);
    if (b != -1) {
      setcolor (winp, cc_color[b]);
      snprintf (text, 11, "%d", b + 1);
      pdftext (winp, 30, 35, text);
    }
    for (ccn = 0; ccn < cc->count; ccn++) {
      if (cc->snap[ccn]->noise || cc->snap[ccn]->count == 0)
        continue;
      VWRITE (ccn);
      memset (pxl_done, 0, 1600*1200/8);

      setcolor (winp, cc_color[ccn]);
      {
        float maxx = -FLT_MAX;
        int spike_at_max = 0;
        limit = count + max_count_per_cluster;
        for (n = 0; n < spike_count && count < limit; n++)
          if (spikedata[n].partition == ccn && !spikedata[n].leave_out && spike_waveform[n]) {
            if ((now = time (0)) != last_time) {
              printf ("%79s\r  scatterplot %d of %d, cluster %d of %d: spike %d of %d\r",
                      "", pln + 1, ph->count, ccn + 1, cc->count, n + 1, spike_count); fflush (stdout);
              last_time = now;
            }
            proj = spike_projection (spike_waveform[n], &ph->plane[pln]);
            if (proj->a > maxx)
              maxx = proj->a,
                spike_at_max = n;
            scatter (proj, &pctr_list[ccn][pln], winp, 0, f);
          }
        if (count >= limit) hit_limit = 1;
        if (0 && ccn == 10 && ((a == 10 && b == 1) || (a == 1 && b == 10)))
          printf ("rightmost spike center %d projection centers %d and %d: %d at %f \n",
                  ccn, a, b, spike_at_max, spikedata[spike_at_max].sample);
      }
      setcolor (winp, 0xa0a0a000U);
      if (0)
        for (n = 0; n < spike_count; n++)
          if (spikedata[n].partition == ccn && spikedata[n].leave_out && spike_waveform[n]) {
            if ((now = time (0)) != last_time) {
              printf ("%79s\r  scatterplot %d of %d, cluster %d of %d: spike %d of %d\r",
                      "", pln + 1, ph->count, ccn + 1, cc->count, n + 1, spike_count); fflush (stdout);
              last_time = now;
            }
            proj = spike_projection (spike_waveform[n], &ph->plane[pln]);
            scatter (proj, &pctr_list[ccn][pln], winp, 0, f);
          }

      setcolor (winp, 0xffffff00U);
      memset (pxl_done, 0, 1600*1200/8);
      proj = spike_projection (cc->snap[ccn]->val, &ph->plane[pln]);
      scatter (proj, &pctr_list[ccn][pln], winp, 0, f);
      CWRITE (0x8000);
    }
    printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
    setcolor (winp, 0);
    CWRITE (0xffff);
    memset (pxl_done, 0, 1600*1200/8);
    printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
    limit = count + max_count_per_cluster;
    for (n = 0; n < ncnt && count < limit; n++) {
      if ((now = time (0)) != last_time) {
        printf ("%79s\r  scatterplot %d of %d, noise: spike %d of %d\r",
                "", pln + 1, ph->count, n + 1, ncnt); fflush (stdout);
        last_time = now;
      }
      scatter (&noise_projections[n][pln], &nctr, winp, 0, f);
    }
    if (count >= limit) hit_limit = 1;
    setcolor (winp, 0xffffff00U);
    memset (pxl_done, 0, 1600*1200/8);
    proj = spike_projection (zero, &ph->plane[pln]);
    scatter (proj, &nctr, winp, 0, f);
    CWRITE (0x8000);
    fwrite (&count, sizeof count, 1, f) == 1 || DIE; /* not VWRITE so count is not incremented */
    end_pdf_page (winp);
    printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
  }
  if (hit_limit) {
    printf ("\n");
    printf ("%s line %d: limited scatterplots", __FILE__, __LINE__);
    printf ("\n");
    note ("%s line %d: limited scatterplots", __FILE__, __LINE__);
  }
  printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
  fclose (f);
  free (pxl_done);
  printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
}

#define SCALE(x) ((1200-1) - (int)(x/32768*(1200/2) + (1200/2)))
#define SCALE_LIMIT(x) x = x < SCALE (32767.0) ? SCALE (32767.0) : (x > SCALE (-32768.0) ? SCALE (-32768.0) : x)

STATIC void
put_wave_data (float *v, short ***wave_data)
{
  short ***p = wave_data;
  int n, a, b, i, new;

  for (n = 0; n + 1 < SNAPLEN; n++) {
    a = SCALE (v[n]);
    SCALE_LIMIT (a);
    b = SCALE (v[n+1]);
    SCALE_LIMIT (b);
    if (p[n][a]) {
      for (i = 1; i < p[n][a][0]; i++) {
        //      printf ("p[%d][%d][%d] = %d, %d?\n", n, a, i, p[n][a][i], b);
        if (p[n][a][i] == b)
          break;
      }
      if (i == p[n][a][0]) {
        new = p[n][a][0]++;
        TREALLOC (p[n][a], p[n][a][0]);
        //      printf ("setting p[%d][%d][%d] to %d\n", n, a, new, b);
        p[n][a][new] = b;
      }
    }
    else {
      TMALLOC (p[n][a], 2);
      p[n][a][0] = 2;
      p[n][a][1] = b;
      //      printf ("setting p[%d][%d][%d] to %d\n", n, a, 1, b);
    }
  }
}

STATIC void
put_scatter_data (float *v, PlaneHdr *ph, Projection *proj)
{
  int pln, i;
  double x, y;

  for (pln = 0; pln < ph->count; pln++) {
    for (x = y = 0, i = 0; i < SNAPLEN; i++) {
      x += v[i] * ph->plane[pln].va[i];
      y += v[i] * ph->plane[pln].vb[i];
    }
    proj[pln].a = x;
    proj[pln].b = y;
  }    
}

#if 1

STATIC int spikeloc_cmp (const int *a, const int *b)
{
  double sa, sb;
  return (sa = spikedata[*a].sample) < (sb = spikedata[*b].sample) ? -1 : (sa > sb); /* ascending */
}

static void
best_waveform (int spikeidx)
{
  float *ar = spike_waveform[spikeidx];


#ifdef MIN_W_R
  int n;
  double  sum = 0,  tmp;
  double sumr = 0, tmpr;
  float aw[SNAPLEN];
  float *br = cc->snap[spikedata[spikeidx].partition]->raw;
  float *bw = cc->snap[spikedata[spikeidx].partition]->val;

  SZMEMCPY (aw, ar);
  whiten_snap (aw);
  for (n = 0; n < SNAPLEN; n++)
    tmp = aw[n] - bw[n], sum += tmp*tmp,
      tmpr = ar[n] - br[n], sumr += tmpr*tmpr;
  if (sqrt (sum) < sqrt(sumr) / unzapped_sd)
    TMEMCPY (ar, aw, SNAPLEN);
  else
    for (n = 0; n < SNAPLEN; n++)
      ar[n] = aw[n] + (ar[n] - br[n]) / unzapped_sd;
#else
  whiten_snap (ar);
#endif

}

STATIC void
add_spike_waveforms_3 (SnapList *cc, short ****wave_data_list)
{
  int n, i;
  static int *idx_map;

  TREALLOC (idx_map, spike_count);
  for (i = 0; i < spike_count; i++)
    idx_map[i] = i;
  qsort (idx_map, (size_t)spike_count, sizeof *idx_map, (int(*)(const void*,const void*))spikeloc_cmp);

  for (i = 0; i < spike_count; i++) {
    n = idx_map[i];
    if (spike_waveform[n] == 0 || spikedata[n].partition < 0 || cc->snap[spikedata[n].partition]->noise)
      continue;
    put_wave_data (spike_waveform[n], wave_data_list[spikedata[n].partition]);
    //    whiten_snap (spike_waveform[n]);
    best_waveform (n);
  }
}

#endif

STATIC Projection **
get_noise_projections (SnapList *noise, PlaneHdr *ph)
{
  int n;
  static Projection **noise_projections;
  Pctr *nctr;

  TCALLOC (nctr, (size_t)ph->count);
  TMALLOC (noise_projections, noise->count);

  for (n = 0; n < noise->count; n++) {
    TMALLOC (noise_projections[n], ph->count);
    put_scatter_data (noise->snap[n]->val, ph, noise_projections[n]);
  }
  free (nctr);
  return noise_projections;
}

STATIC PlaneHdr *
get_planes (SnapList *cc)
{
  int a, b, ccnt, i, n, plane_count;
  PlaneHdr *ph;
  float *va, *vb;
  float v[SNAPLEN];
  double m;

  for (ccnt = n = 0; n < cc->count; n++)
    if (!cc->snap[n]->noise && cc->snap[n]->count > 0)
      ccnt++;
  plane_count = ccnt == 1 ? 1 : ccnt * (ccnt - 1) / 2;
  
  TMALLOC (ph, 1);
  ph->count = plane_count;
  TMALLOC (ph->plane, plane_count);
  for (n = a = 0; a < cc->count; a++) {
    if (cc->snap[a]->noise || cc->snap[a]->count == 0)
      continue;
    for (b = a + 1; b < cc->count || ccnt == 1; b++) {
      if (b < cc->count && (cc->snap[b]->noise || cc->snap[b]->count == 0))
        continue;

      ph->plane[n].a = a;
      ph->plane[n].b = ccnt == 1 ? -1 : b;
      va = ph->plane[n].va;
      vb = ph->plane[n].vb;
      n++;

      if (a < 0 || a >= cc->count)
        exit (DIE);
      memcpy (va, cc->snap[a]->val, sizeof (SnapVals));
      m = distance (va, zero);
      for (i = 0; i < SNAPLEN; i++)
        va[i] /= m;
      if (ccnt == 1) {
        float tmp;
        memcpy (vb, cc->snap[a]->val, sizeof (SnapVals));
        tmp = vb[0];
        vb[0] = vb[1];
        vb[1] = -tmp;
      }
      else
        memcpy (vb, cc->snap[b]->val, sizeof (SnapVals));
      for (m = 0, i = 0; i < SNAPLEN; i++)
        m += va[i] * vb[i];
      memcpy (v, va, sizeof (SnapVals));
      for (i = 0; i < SNAPLEN; i++)
        v[i] *= m;
      for (i = 0; i < SNAPLEN; i++)
        vb[i] -= v[i];
      m = distance (vb, zero);
      for (i = 0; i < SNAPLEN; i++)
        vb[i] /= m;
      if (ccnt == 1)
        break;
    }
  }
  return ph;
}

STATIC Pctr **
get_pctr_list (SnapList *cc, PlaneHdr *ph)
{
  int ccn, pln, i;
  Pctr **pctr;

  TMALLOC (pctr, cc->count);

  for (ccn = 0; ccn < cc->count; ccn++) {
    TCALLOC (pctr[ccn], (size_t)ph->count);
    for (pln = 0; pln < ph->count; pln++) {
      for (i = 0; i < SNAPLEN; i++) {
        pctr[ccn][pln].a += cc->snap[ccn]->val[i] * ph->plane[pln].va[i];
        pctr[ccn][pln].b += cc->snap[ccn]->val[i] * ph->plane[pln].vb[i];
      }
    }
  }
  return pctr;
}

STATIC short ****
alloc_wave_data_list (SnapList *cc)
{
  short ****wdl;
  int ccn, n;

  TMALLOC (wdl, cc->count + 1);
  for (ccn = 0; ccn <= cc->count; ccn++) {
    TMALLOC (wdl[ccn], SNAPLEN-1);
    for (n = 0; n < SNAPLEN-1; n++)
      TCALLOC (wdl[ccn][n], 1200);
  }
  return wdl;
}

STATIC inline int
isfull (char *p, int w)
{
  int n;

  for (n = 0; n <= w; n++)
    if (p[n] == 0)
      return 0;
  return 1;
}

STATIC void
find_block (short **lines, int w, int *startp, int *endp)
{
  char q[1200][27];
  int ys, y, n, x, y_at_max = 0, count, max;

  memset (q, 0, sizeof q);

  for (ys = 0; ys < 1200; ys++)
    if (lines[ys]) {
      for (n = 1; n < lines[ys][0]; n++)
        for (x = 0; x <= w; x++) {
          y = (int)floor (ys + (double)x / w * (lines[ys][n] - ys) + .5);
          if (x >= 27 || x < 0)
            exit (DIE);
          if (y < 0 || y >= 1200) {
            printf ("ys: %d, x: %d, w: %d, d: %d, y: %d\n", ys, x, w, lines[ys][n], y);
            exit (DIE);
          }
          q[y][x] = 1;
        }
    }

  for (count = max = y = 0; y < 1200; y++)
    if (isfull (q[y], w)) {
      if (++count > max)
        max = count,
          y_at_max = y;
    }
    else count = 0;
  *startp = y_at_max - max + 1;
  *endp = y_at_max;
}

#define scaleW(i) ((int)(((i) * ((1600-1.0)/63.0))))
STATIC void
write_wdt (char *file_name, SnapList *cc, short ****wave_data_list)
{
  FILE * f;
  unsigned short ccn, ccnt, i, left, right, y, ye, flagval, buf;
  const unsigned short endmark = 0xc000;
  int segcnt, start, end, n, count, segcnt2;
  short rectangle[(SNAPLEN-1)*2];
  short sv[SNAPLEN];

  printf ("%79s\r  write wdt\r", ""); fflush (stdout);
  (f = fopen (change_filetype (file_name, ".wdt"), "wb")) || DIE;

  VWRITE (mean);
  VWRITE (unzapped_min);
  VWRITE (unzapped_max);
  for (ccnt = ccn = 0; (int)ccn < cc->count; ccn++)
    ccnt += !cc->snap[ccn]->noise && cc->snap[ccn]->count > 0;
  printf ("\n%s line %d: %d non-noise centers of %d\n", __FILE__, __LINE__, ccnt, cc->count);
  VWRITE (ccnt);
  for (ccn = 0; (int)ccn < cc->count; ccn++) {
    if (cc->snap[ccn]->noise || cc->snap[ccn]->count == 0)
      continue;
    count = 0;
    VWRITE (ccn);
    for (segcnt = 0, i = 0; i < 63; i++) {
      left = (unsigned short)scaleW (i);
      right = (unsigned short)scaleW (i+1);
      find_block (wave_data_list[ccn][i], (int)(right - left), &start, &end);
      rectangle[i*2] = start;
      rectangle[i*2+1] = end;
      for (y = 0; y < 1200; y++)
        if (wave_data_list[ccn][i][y])
          for (n = 1; n < wave_data_list[ccn][i][y][0]; n++) {
            ye = (unsigned short)wave_data_list[ccn][i][y][n];
            if (! ((int)y >= start && (int)y <= end && (int)ye >= start && (int)ye <= end))
              segcnt++;
          }
    }
    ORWRITE (cluster_spikestart[ccn], 0);
    VWRITE (segcnt);
    VWRITE (rectangle);
    for (i = 0; i < SNAPLEN; i++) {
      sv[i] = SCALE (cc->snap[ccn]->raw[i]);
      SCALE_LIMIT (sv[i]);
    }
    VWRITE (sv);
    for (segcnt2 = 0, i = 0; i < 63; i++) {
      int newx = 1;
      start = rectangle[i*2];
      end = rectangle[i*2+1];
      for (y = 0; y < 1200; y++) {
        int newy = 1;
        if (wave_data_list[ccn][i][y]) {
          for (n = 1; n < wave_data_list[ccn][i][y][0]; n++) {
            ye = (unsigned short)wave_data_list[ccn][i][y][n];
            if (! ((int)y >= start && (int)y <= end && (int)ye >= start && (int)ye <= end)) {
              if (newy) {
                flagval = (unsigned short)(newx ? 0x8000 : 0x4000);
                ORWRITE (y, flagval);
                newx = newy = 0;
              }
              VWRITE (ye);
              segcnt2++;
            }
          }
        }
      }
      if (newx)
        CWRITE (0x8000);
    }
    if (segcnt2 != segcnt) {
      printf ("line %d: segcnt %d segcnt2 %d\n", __LINE__, segcnt, segcnt2);
      exit (DIE);
    }

    //    printf ("line %d: segcnt %d\n", __LINE__, segcnt);
    VWRITE (endmark);
    fwrite (&count, sizeof count, 1, f) == 1 || DIE;
  }
  fclose (f);
}
#undef VWRITE
#undef ORWRITE

double
white_distance_from_zero (float *a, int count)
{
  float v[SNAPLEN];
  double sum = 0;
  int n;

  memcpy (v, a, count * sizeof *v);
  whiten_samples (v, count);
  for (n = 0; n < count; n++)
    sum += v[n] * v[n];
  return map_radius (sqrt(sum), count);
}

static double
white_distance_from_zero_snaplen (float *a)
{
  float v[SNAPLEN];
  double sum = 0;
  int n;

  SZMEMCPY (v, a);
  whiten_snap (v);
  for (n = 0; n < SNAPLEN; n++)
    sum += v[n] * v[n];
  return sqrt(sum);
}

double
norm_distance_from_zero (float *a)
{
  float v[SNAPLEN];
  double sum = 0;
  double tmp;
  int n;
  float fmean = mean;

  memcpy (v, a, SNAPLEN * sizeof *v);
  for (n = 0; n < SNAPLEN; n++)
    tmp = v[n] - fmean,
      sum += tmp * tmp / cvec[0];
  return sqrt(sum);
}

#define D01_MAX 10
#define D12_MAX 25
#define D02_MAX (D01_MAX + D12_MAX)

STATIC inline int
do_mpk_4 (int a, int b, int c, int dir, RegionData *p, /*@only@*/ float *buf, int bufsamples, int starting_sample, int preexisting)
{
  int i, stop, lnear, rnear;
  int r0 = MPK_R0;
  static int r = MPK_R;
  double d;
  int b0 = b;
  int peakloc = starting_sample + b;
  //  int debug = (i = peakloc - 1207083) > 0 && i < SNAPLEN*2;
  //  int debug = peakloc > 255607615 + 63 && peakloc < 255607615 + 81;
  int debug = 0;
  int val = (int)(buf[b] - mean);

  if (fabs(buf[b]) > CLIPVAL)
    r0 += 2;

  if (dir * (buf[b] - mean) < 0) {
    if (debug) printf ("line %d: peak at %d %6d %2d doesn't cross mean\n", __LINE__, peakloc, val, dir);
    return 0;
  }

  stop = MAX (b - (2 * r0 - 1), 0);
  for (i = b - 1; i >= stop && buf[b] == buf[i]; i--);
  a = MAX (b - (2 * r0), 0);
  for ( ; i >= stop && dir * (buf[b] - buf[i]) >= 0; i--)
    if ((double)(dir * (buf[b] - buf[i])) >= detection_threshold) {
      a = i;
      break;
    }
  if (i < 0) {
    if (debug) printf ("line %d: peak at %d %6d %2d overlaps beginning of buffer\n", __LINE__, peakloc, val, dir);
    return 0;
  }

  stop = MIN (b + (2 * r0 - 1), bufsamples - 1);
  for (i = b + 1; i <= stop && buf[b] == buf[i]; i++);
  c = MIN (b + (2 * r0), bufsamples - 1);
  for ( ; i <= stop && dir * (buf[b] - buf[i]) > 0; i++)
    if ((double)(dir * (buf[b] - buf[i])) >= detection_threshold) {
      c = i;
      break;
    }
  if (i >= bufsamples) {
    if (debug) printf ("line %d: peak at %d %6d %2d overlaps end of buffer\n", __LINE__, peakloc, val, dir);
    return 0;
  }
  if ((double)(dir * (buf[b] - buf[c])) < detection_threshold || (double)(dir * (buf[b] - buf[a])) < detection_threshold || c - a > 2 * r0) {
    if (debug) printf ("line %d: peak at %d %6d %2d below detection threshold\n", __LINE__, peakloc, val, dir);
    return 0;
  }

  lnear = b0 - a;
  rnear = c - b0;

  if (b - a > r0 || c - b > r0)
    b = (a + c) / 2;

  stop = MAX (b - r, 0);
  for (i = a - 1; i >= stop && dir * (buf[b] - buf[i]) >= 0; i--)
    if (dir * (buf[a] - buf[i]) > 0)
      a = i;
  if (i < 0) {
    if (debug) printf ("line %d: peak at %d %6d %2d overlaps beginning of buffer\n", __LINE__, peakloc, val, dir);
    return 0;
  }

  stop = MIN (b + r, bufsamples - 1);
  for (i = c + 1; i <= stop && dir * (buf[b] - buf[i]) > 0; i++)
    if (dir * (buf[c] - buf[i]) > 0)
      c = i;
  if (i >= bufsamples) {
    if (debug) printf ("line %d: peak at %d %6d %2d overlaps end of buffer\n", __LINE__, peakloc, val, dir);
    return 0;
  }

#define LIM(a) MIN (a, bufsamples - SNAPLEN)
#ifdef BIGFIRST
  if ((d = white_distance_from_zero (buf + LIM (a), c - a + 1)) > 9) {
#else
  if ((d = white_distance_from_zero (buf + LIM (a), c - a + 1)) > (preexisting ? 11.5 : 9)) {
#endif
    {
      static int end_of_last_snapshot;
      double ds = 0;
      if (p->data_count == 0)
        end_of_last_snapshot = 0;
      if (peakloc < end_of_last_snapshot
          || (b0 >= PRESAMPLES && (ds = white_distance_from_zero (buf + LIM ((b0 - PRESAMPLES)), SNAPLEN)) > (preexisting == 1 ? 10.0 : 11.0))) {
        region_add_data (p, b0, b0 - a, lnear, rnear, c - b0, buf, bufsamples, starting_sample, preexisting);
        end_of_last_snapshot = peakloc - PRESAMPLES + SNAPLEN;
        if (debug) printf ("%s line %d: peak detected at %d\n", __FILE__, __LINE__, peakloc);
        b0 || DIE;
        return b0;
      }
      else {
        if (debug) printf ("%s line %d: distance %d from zero is only %f\n", __FILE__, __LINE__, peakloc, ds);
        return 0;
      }
    }
  }
  else {
    if (debug) printf ("%s line %d: distance %d %d %d from zero is only %f\n", __FILE__, __LINE__, b0-a, peakloc, c-b0, d);
    return 0;
  }
}

int
find_mpk_spikes (RegionData *rd, float *buf, int bufsamples, int starting_sample, int buf_start, int preexisting)
{
  int n, i, last, dir, p[3], count, done, pn;

  last = SGN (buf[buf_start + 1] - buf[buf_start]);
  i = 0;
  done = 0;
  count = 0;
  for (n = buf_start + 2; n < bufsamples; n++) {
    if ((dir = SGN (buf[n] - buf[n-1])) && dir != last) {
      p[i++] = n - 1;
      count++;
      last = dir;
      if (i == 2)
        break;
    }
  }
  for ( ; n < bufsamples; n++) {
    if ((dir = SGN (buf[n] - buf[n-1])) && dir != last) {
      p[i] = n - 1;
      count++;
      i = (i + 1) % 3;
      done = do_mpk_4 (p[i], p[(i + 1) % 3], p[(i + 2) % 3], dir, rd, buf, bufsamples, starting_sample, preexisting);

      if (bufsamples == SNAPLEN)
      {                         /* debug */
        //int sample = starting_sample + p[(i + 1) % 3];
        /* debug */
        //if (sample >= 44815200 && sample <= 44820000)
        printf ("peakloc %d %d %d, %2d: %d \n", p[i], p[(i + 1) % 3], p[(i + 2) % 3], dir, done);
      }

      last = dir;
    }
  }
  //  count + done > 3 || DIE;
  pn = done ? p[(i + 1) % 3] : p[i];
  (void)(MPK_R > 1 || DIE);
  return pn - BEFORE - 1;
}

static int pass_hist[100];

float *
optavg (SnapList *sl, int first_snap, int last_snap, float *p0, float *buf, int starting_sample, int bufsamples)
{
  static char map[(SNAPLEN+PAD*2)];
  static char hasv[5][4] = {{0,0,1,1},{0,1,1,1},{1,1,1,1},{1,1,1,0},{1,1,0,0}};
  int snapidx, i, n, j, k, e, snap_count, all_aligned;
  static int alloc;
  static int *bufidx;
  static double *x1, *x2, *x3;
  static double *xa[5], *xb[5], *xc[5], *x_[5], *xd[5], *xe[5], *xf[5];
  double peakloc, x4, x5, x6, x, ya, yb, yc, yd, ye, yf, m, b, change, tmp;
  double v[5];
  static float avg_array[(SNAPLEN+PAD*2)+6];
  float *avg = avg_array + 3;
  time_t long_time;
  int new, pass_count;

  //  printf ("%79s\roptavg %d %d\n", "", first_snap, last_snap);

  if (first_snap == last_snap) {
    double offset;
    float *bufp;

    peakloc = sl->snap[first_snap]->sample;
    offset = peakloc - ceil (peakloc) + 1;
    bufp = buf + ((int)ceil (peakloc) - PRESAMPLES - starting_sample - 2) - PAD;
    for (i = 0; i < SNAPLEN+PAD*2; i++)
      avg[i] = local_cubic (i+offset+1, bufp, SNAPLEN+PAD*2+3, mean);
    return avg;
  }

  if (map[1] == 0) {
    map[0] = 0, map[1] = 1, map[(SNAPLEN+PAD*2)-2] = 3, map[(SNAPLEN+PAD*2)-1] = 4;
    for (n = 2; n < (SNAPLEN+PAD*2)-2; n++)
      map[n] = 2;
  }
  if ((snap_count = last_snap - first_snap + 1) > alloc) {
    TREALLOC (bufidx, snap_count);
    TREALLOC (x1, snap_count);
    TREALLOC (x2, snap_count);
    TREALLOC (x3, snap_count);
    for (n = 0; n < 5; n++) TREALLOC (xa[n], snap_count);
    for (n = 0; n < 5; n++) TREALLOC (xb[n], snap_count);
    for (n = 0; n < 5; n++) TREALLOC (xc[n], snap_count);
    for (n = 0; n < 5; n++) TREALLOC (x_[n], snap_count);
    for (n = 0; n < 5; n++) TREALLOC (xd[n], snap_count);
    for (n = 0; n < 5; n++) TREALLOC (xe[n], snap_count);
    for (n = 0; n < 5; n++) TREALLOC (xf[n], snap_count);
    alloc = snap_count;
  }
  for (n = 0; n < (SNAPLEN+PAD*2)+6; n++)
    avg_array[n] = 0;
  if (p0)
    for (n = 0; n < SNAPLEN; n++)
      avg[PAD+n] = p0[n] - mean;

  all_aligned = 1;
  for (snapidx = first_snap; snapidx <= last_snap; snapidx++) {
    peakloc = sl->snap[snapidx]->sample;
    n = snapidx - first_snap;
    bufidx[n] = (int)(ceil(peakloc) - starting_sample - (PRESAMPLES+PAD) - 2);
    x1[n] = x = ceil(peakloc) - peakloc;
    if (x != 0) all_aligned = 0;
    x6 = (x5 = (x4 = (x3[n] = (x2[n] = x*x)*x)*x)*x)*x;

    xa[2][n] = -x6+3*x5-3*x4+x3[n];
    xb[2][n] = 6*x6-18*x5+15*x4-3*x2[n];
    xc[2][n] = -15*x6+45*x5-33*x4-9*x3[n]+12*x2[n];
    x_[2][n] = 20*x6-60*x5+42*x4+16*x3[n]-18*x2[n]+4;
    xd[2][n] = -15*x6+45*x5-33*x4-9*x3[n]+12*x2[n];
    xe[2][n] = 6*x6-18*x5+15*x4-3*x2[n];
    xf[2][n] = -x6+3*x5-3*x4+x3[n];

    xa[0][n] = 0;
    xb[0][n] = 0;
    xc[0][n] = -3*x6+11*x5-13*x4+3*x3[n]+4*x2[n]-2*x;
    x_[0][n] = 10*x6-34*x5+31*x4+8*x3[n]-19*x2[n]+4;
    xd[0][n] = -12*x6+38*x5-30*x4-8*x3[n]+12*x2[n];
    xe[0][n] = xe[2][n];
    xf[0][n] = xf[2][n];

    xa[1][n] = 0;
    xb[1][n] = 3*x6-10*x5+10*x4-2*x3[n]-x2[n];
    xc[1][n] = -12*x6+38*x5-30*x4-8*x3[n]+12*x2[n];
    x_[1][n] = 19*x6-58*x5+41*x4+16*x3[n]-18*x2[n]+4;
    xd[1][n] = xd[2][n];
    xe[1][n] = xe[2][n];
    xf[1][n] = xf[2][n];

    xa[3][n] = xa[2][n];
    xb[3][n] = xb[2][n];
    xc[3][n] = xc[2][n];
    x_[3][n] = 19*x6-56*x5+36*x4+20*x3[n]-19*x2[n]+4;
    xd[3][n] = -12*x6+34*x5-20*x4-12*x3[n]+8*x2[n]+2*x;
    xe[3][n] = 3*x6-8*x5+5*x4+2*x3[n]-2*x2[n];
    xf[3][n] = 0;

    xa[4][n] = xa[2][n];
    xb[4][n] = xb[2][n];
    xc[4][n] = -12*x6+34*x5-20*x4-12*x3[n]+8*x2[n]+2*x;
    x_[4][n] = 10*x6-26*x5+11*x4+8*x3[n]+x2[n];
    xd[4][n] = -3*x6+7*x5-3*x4-x3[n];
    xe[4][n] = 0;
    xf[4][n] = 0;

  }
  if (bufidx[snap_count - 1] + (SNAPLEN+PAD*2) - 1 + 3 >= bufsamples)
    return 0;
  long_time = time (0) + 2;
  new = 1;
  pass_count = 0;
  do {
    if (++pass_count > 1000)
      return 0;
    change = 0;
    for (i = 0; i < (SNAPLEN+PAD*2); i++) {
      e = map[i];
      ya = avg[i-3]; yb = avg[i-2]; yc = avg[i-1]; yd = avg[i+1]; ye = avg[i+2]; yf = avg[i+3];
      m = b = 0;
      for (n = 0; n < snap_count; n++) {
        for (k = 0; k < 4; k++) {
          v[k] = hasv[e][k] ? (buf[bufidx[n] + i + k]-mean) : 0;
        }
        for (j = 0; j < snap_count; j++)
          if (j != n)
            for (k = 0; k < 4; k++)
              if (hasv[e][k]) {
                int idx = i - 2 + k + bufidx[n] - bufidx[j];
                if (idx >= 0 && idx < (SNAPLEN+PAD*2)-1) {
                  v[k] -= local_cubic (1 + x1[j], avg + idx - 1, 4, 0);
                }
              }
        m += x_[e][n];
        b += (xa[e][n]*ya + xb[e][n]*yb + xc[e][n]*yc + xd[e][n]*yd + xe[e][n]*ye + xf[e][n]*yf
              +(2*v[3]-6*v[2]+6*v[1]-2*v[0])*x3[n] + (-4*v[3]+10*v[2]-8*v[1]+2*v[0])*x2[n] + (2*v[3]-2*v[1])*x1[n] - 4*v[2]);
      }
      if (m == 0) {
        (all_aligned && i == (SNAPLEN+PAD*2)-1) || DIE;
        avg[i] = buf[(int)sl->snap[last_snap]->sample - (PRESAMPLES+PAD) - starting_sample + (SNAPLEN+PAD*2) - 1] - mean;
        //printf ("\n%d %d %d %d %d\n", i, e, snap_count, x1[0] == 0, x1[1] == 0);
      }
      else {
        tmp = -b / m;
        //      printf ("b: %f, m: %f, -b/m: %f avg: %f\n", b, m, tmp, avg[i]);
        change += avg[i] - (float)tmp;
        avg[i] = tmp;
      }
    }

    if ((now = time(0)) > long_time) {
      if (new) {
        printf ("\nchange %.0f %.0f\n", sl->snap[first_snap]->sample, sl->snap[last_snap]->sample);
        new = 0;
      }
      printf ("%f %d\n", change, pass_count);
      long_time = now;
    }
  } while (change != 0);
  (void)(pass_count / 1000 < 100 && pass_hist[pass_count / 1000]++);

  for (n = 0; n < (SNAPLEN+PAD*2); n++)
    avg[n] += mean;
  return avg;
}


#define NSZ 128

STATIC void
add_mean_snap (SnapList *sl, int newcount)
{
  int i, ccn;
  float v[SNAPLEN];

  TREALLOC (cluster_spikestart, newcount);
  for (i = 0; i < SNAPLEN; i++)
    v[i] = mean;
  for (ccn = sl->count; ccn < newcount; ccn++) {
    add_snap_alloc (sl, v);
    cluster_spikestart[ccn] = PRESAMPLES;
    sl->snap[ccn]->nocenter = 1;
  }
}

STATIC int
read_white (void)
{
  FILE *f;
  if (!((f = fopen (change_filetype (file_name, ".wht"), "rb"))
        && fread (&mean, sizeof mean, 1, f) == 1
        && fread (cov_rep_root, sizeof cov_rep_root, 1, f) == 1
        && fread (cvec, sizeof cvec, 1, f) == 1)) {
    return 0;
  }
  fclose (f);
  printf ("read cov_rep_root from %s\n", change_filetype (file_name, ".wht"));
  return 1;
}

STATIC void
write_white (void)
{
  FILE *f;
  (f = fopen (change_filetype (file_name, ".wht"), "wb")) || DIE;
  fwrite (&mean, sizeof mean, 1, f) == 1 || DIE;
  fwrite (cov_rep_root, sizeof cov_rep_root, 1, f) == 1 || DIE;
  fwrite (cvec, sizeof cvec, 1, f) == 1 || DIE;
  fclose (f);
}

STATIC void
write_path (void)
{
  FILE *f;
  (f = fopen (change_filetype (file_name, ".path"), "wb")) || DIE;
  fwrite (file_name, strlen (file_name) + 1, 1, f) == 1 || DIE;
  fclose (f);
}

STATIC void
write_pdf (SnapList *cc,short ****wave_data_list, Projection **noise_projections, SnapList *noise, Pctr **pctr_list, PlaneHdr *ph)
{
  int h = 1200, w = 1600, scale = 32768, n, ccn, i, left, right, y;
  int start, end, ye;
  static Snap pdfsnap;
  char txt[80];
  FILE *winp;

  winp = create_pdf_file (change_filetype (file_name, ".ps"));
  printf ("%79s\r  writing waveform overlays\r", ""); fflush (stdout);
  for (ccn = 0; ccn < cc->count; ccn++) {
    if (cc->snap[ccn]->noise || cc->snap[ccn]->count == 0)
      continue;
    start_pdf_page (winp);
    pdflines (winp, w, h, 0, 0, 0xff000000U);

    pdfrender (&hlines0, winp, w, h, 0, 0, scale, 0xff000000U);

    setcolor (winp, 0);
    for (i = 0; i < 63; i++) {
      left = scaleW (i);
      right = scaleW (i+1);
      fprintf (winp, "/p {%d exch moveto %d exch rlineto stroke} bind def\n", left, right - left);
      find_block (wave_data_list[ccn][i], right - left, &start, &end);
      if (start <= end)
        fprintf (winp, "%d %d moveto %d %d lineto %d %d lineto %d %d lineto fill\n",
                 left, start, right, start, right, end, left, end);
      for (y = 0; y < 1200; y++)
        if (wave_data_list[ccn][i][y])
          for (n = 1; n < wave_data_list[ccn][i][y][0]; n++) {
            ye = wave_data_list[ccn][i][y][n];
            if (! (y >= start && y <= end && ye >= start && ye <= end))
              fprintf (winp, "%d %d p\n", wave_data_list[ccn][i][y][n] - y, y);
          }
    }
    pdfsnap.val = cc->snap[ccn]->raw;
    pdfrender (snap_to_list (&pdfsnap), winp, w, h, 0, 0, scale, 0xffff0000U);
    snprintf (txt, 80, "%d", ccn + 1);
    setcolor (winp, 0);
    pdftext (winp, 10, 35, txt);
    end_pdf_page (winp);
  }
  printf ("%79s\r  generating scatterplots\r", ""); fflush (stdout);
  scatterplot (winp, cc, noise_projections, noise->count, pctr_list, ph);

  if (0) {
    printf ("%79s\r  converting postscript to pdf\r", ""); fflush (stdout);
#if HAVE_SPAWNVP
    ps2pdf ();
#endif
    pdfdisplay (winp);
    rename ("psrender.pdf", change_filetype (file_name, ".pdf"));
  }
  else pdfdisplay (winp);

  printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
}

typedef struct {double weight, shift; size_t va, vb;} Edge;

static int edge_cmp (const Edge *a, const Edge *b)
{
  return a->weight < b->weight ? -1 : a->weight > b->weight; /* ascending */
}

#define SWAP(a, b, tmp) do {tmp = a; a = b; b = tmp;} while (0)

#if 0
static size_t *
spanning_tree (Edge *edge, size_t vcnt)
{
  size_t ecnt, v, e, tmp;
  char taken[vcnt];
  size_t *eorder;

  vcnt > 0 || DIE;
  ecnt = vcnt * (vcnt - 1) / 2;
  qsort (edge, ecnt, sizeof *edge, (int(*)(const void*,const void*))edge_cmp);
  TMALLOC (eorder, vcnt);
  memset (taken, 0, vcnt);
  taken[edge[0].va] = 1;
  printf ("vcnt %ld\n", (long)vcnt); fflush (stdout);
  for (v = 1; v < vcnt; v++) {
    for (e = 0; e < ecnt; e++)
      if (taken[edge[e].va] + taken[edge[e].vb] == 1)
        break;
    e < ecnt || DIE;
    if (taken[edge[e].vb]) {
      SWAP (edge[e].va, edge[e].vb, tmp);
      edge[e].shift *= -1;
    }
    taken[edge[e].vb] = 1;
    eorder[v-1] = e;
    printf ("eorder[%ld] = %ld\n", (long)v-1, (long)e); fflush (stdout);
  }
  return eorder;
}
#endif

static void
measure_centers (SnapList *cc)
{
  int n, i;
  float *v;
  double tmp, sum;

  for (n = 0; n < cc->count; n++) {
    v = cc->snap[n]->val;
    sum = 0;
    for (i = 0; i < SNAPLEN; i++) {
      tmp = v[i];
      sum += tmp * tmp;
    }
    cc->snap[n]->distance = sqrt (sum);
  }
}

typedef struct {int ccn; double sz;} CtrSize;

static int ctr_cmp_dn (const CtrSize *a, const CtrSize *b)
{
  return a->sz > b->sz ? -1 : (a->sz < b->sz); /* descending */
}

static int *
center_by_size (SnapList *cc)
{
  CtrSize *ctrsz;
  int *ctr;
  int n;

  TMALLOC (ctrsz, cc->count);
  TMALLOC (ctr, cc->count);
  measure_centers (cc);
  for (n = 0; n < cc->count; n++) {
    ctrsz[n].ccn = n;
    ctrsz[n].sz = cc->snap[n]->distance;
  }
  qsort (ctrsz, (size_t)cc->count, sizeof *ctrsz, (int(*)(const void*,const void*))ctr_cmp_dn);
  for (n = 0; n < cc->count; n++)
    ctr[n] = ctrsz[n].ccn;
  free (ctrsz);
  return ctr;
}

static size_t *
amplitude_tree (Edge *edge, size_t vcnt, int *ctr)
{
  size_t ecnt, v, e, tmp;
  char taken[vcnt];
  size_t *eorder;
  int cidx;

  vcnt > 0 || DIE;
  ecnt = vcnt * (vcnt - 1) / 2;
  qsort (edge, ecnt, sizeof *edge, (int(*)(const void*,const void*))edge_cmp);
  TMALLOC (eorder, vcnt);
  memset (taken, 0, vcnt);
  cidx = 0;
  taken[ctr[cidx++]] = 1;
  taken[ctr[cidx]] == 0 || DIE;
  printf ("vcnt %ld\n", (long)vcnt); fflush (stdout);
  for (v = 1; v < vcnt; v++) {
    for (e = 0; e < ecnt; e++) {
      if (taken[edge[e].va] && edge[e].vb == ctr[cidx])
        break;
      if (taken[edge[e].vb] && edge[e].va == ctr[cidx]) {
        SWAP (edge[e].va, edge[e].vb, tmp);
        edge[e].shift *= -1;
        break;
      }
    }
    e < ecnt || DIE;
    cidx++;
    taken[edge[e].vb] = 1;
    eorder[v-1] = e;
    printf ("eorder[%ld] = %ld\n", (long)v-1, (long)e); fflush (stdout);
  }
  return eorder;
}

void
show_centers (SnapList *avg)
{
  size_t vcnt = (size_t)avg->count;
  size_t n;
  static int c;
  char shw[vcnt];
  char msg[256];

  avg->count >= 0 || DIE;
  memset (shw, 0, sizeof shw);
  window ();
  while (c != 'c' && c != 'n') {
    maybe_flip (shw, vcnt, c);
    clear ();
    draw_snap (mean_snap0.raw, "black");
    for (n = 0; n < vcnt; n++)
      if (shw[n]) {
        //        draw_snap (avg->snap[n]->raw, "red");
        draw_snap (avg->snap[n]->raw, "black");
        if (n < 255)
          snprintf (msg + n, 256 - n, "%c", ltr ((int)n));
      }
      else snprintf (msg + n, 256 - n, "%c", ' ');
    draw_text (msg);
    c = show ();
    if (c == 'q')
      exit (0);
  }
  wclose ();
  if (c == 'n')
    c = 0;
}

static Snap *
cluster_avg_full (SnapList *sl, Snap *s)
{
  int n, i;
  double sum[SNAPLEN], sumr[SNAPLEN];
  Snap *avg;

  TCALLOC (avg, 1);
  TMALLOC (avg->val, SNAPLEN); avg->free_val = 1;
  TMALLOC (avg->raw, SNAPLEN); avg->free_raw = 1;
  memset (sum, 0, sizeof sum);
  memset (sumr, 0, sizeof sumr);

  for (n = 0; n < sl->count; n++) {
    sl->snap[n]->radius -= s->radius;
    if (sl->snap[n]->radius < radius) {
      sl->snap[n]->radius = distance_full (sl->snap[n], s);
      if (sl->snap[n]->radius < radius) {
        for (i = 0; i < SNAPLEN; i++)
          sum[i] += sl->snap[n]->val[i],
            sumr[i] += sl->snap[n]->raw[i];
        avg->count++;
      }
    }
  }
  for (i = 0; i < SNAPLEN; i++)
    avg->val[i] = sum[i] / avg->count,
      avg->raw[i] = sumr[i] / avg->count;
  avg->refcnt = avg->copies = 1;
  avg->radius = s->radius;
  return avg;
}

SnapList *
center_points_full (SnapList *sl, Snap *s)
{
  int i, loop_count;
  Snap *avg, *last_avg = 0;
  SnapList *ctr;
  static int max;

  TCALLOC (ctr, 1);

  for (i = 0; i < sl->count; i++)
    sl->snap[i]->radius = 0;
  s->radius = 0;
  avg = cluster_avg_full (sl, s);
  if (avg->count == 1) {
    free_snap (&avg);
    return ctr;
  }
  loop_count = 0;
  do {
    double d;
    if (last_avg)
      d = distance_full (avg, last_avg);
    else
      d = distance_full (avg, s);
    free_snap (&last_avg);
    last_avg = avg;
    last_avg->radius = d;
    avg = cluster_avg_full (sl, last_avg);
    loop_count++;
  } while (loop_count < 1000 && (avg->count != last_avg->count || memcmp (avg->val, last_avg->val, sizeof (SnapVals)) != 0));
  if (loop_count < 1000 && loop_count > max) {
    max = loop_count;
    //printf ("loop_count max %d\n", max); fflush (stdout);
  }
  if (loop_count >= 1000)
    printf ("loop_count = %d\n", loop_count);
  avg->radius = 0;            //because avg->val[n] == last_avg->val[n]
  free_snap (&last_avg);
                           
  for (i = 0; i < ctr->count; i++)
    if (distance_full (ctr->snap[i], avg) < .04) {
      ctr->snap[i]->copies++;
      break;
    }
  if (i == ctr->count && avg->count > 3)
    add_snap (ctr, avg);
  free_snap (&avg);

  return ctr;
}

void
region_filter (RegionData *p)
{
  int rn, i;
  time_t start;
  static int alloc;
  static double *x, *y;
  static float *a, *b;
  
  start = time (0);
  for (rn = 0; rn < p->region_count; rn++) {
    p->region[rn].count >= 0 || DIE;
    if (alloc < p->region[rn].count) {
      TREALLOC (x, p->region[rn].count);
      for (i = alloc; i < p->region[rn].count; i++)
        x[i] = (double)i;
      alloc = p->region[rn].count;
      TREALLOC (a, alloc);
      TREALLOC (b, alloc);
      TREALLOC (y, alloc);
    }
    for (i = 0; i < p->region[rn].count; i++)
      a[i] = ss_cspline_eval (p->regionm[rn].spline, (double)i, p->regionm[rn].acc);
    
    filter_buf (a, p->region[rn].count);
    filter_buf (a, p->region[rn].count);
    filter_buf (a, p->region[rn].count);
    for (i = 0; i < p->region[rn].count; i++)
      y[i] = a[i];
    gsl_interp_accel_free (p->regionm[rn].acc);
    ss_spline_free (p->regionm[rn].spline);
    p->regionm[rn].acc = gsl_interp_accel_alloc ();
    p->regionm[rn].spline = gsl_spline_alloc (gsl_interp_cspline, (size_t)p->region[rn].count);
    gsl_spline_init (p->regionm[rn].spline, x, y, (size_t)p->region[rn].count);

  }
  printf ("region_filter: %ld seconds\n", (long)(time(0) - start));
}

static double *a_dbl;

double
one_region_noise (float *v, int count, int *maxidx)
{
  static int alloc;
  static float *a, *b;
  double noise, max_noise;
  int i;

  if (alloc < count) {
    alloc = count;
    TREALLOC (a, alloc);
    TREALLOC (b, alloc);
  }
  TMEMCPY (a, v, count);
  TMEMCPY (b, v, count);

  //  lowpass (a, count);
  filter_buf (a, count);
  for (i = 0; i < count; i++)
    b[i] -= a[i];
  max_noise = 0;
  for (i = 0; i < count - SNAPLEN; i++)
    if ((noise = white_distance_from_zero (b + i, SNAPLEN)) > max_noise)
      max_noise = noise, *maxidx = i;
  return max_noise;
}

void
region_noise (RegionData *p)
{
  int rn, good, bad, good_clipped, bad_clipped, maxidx = 0, i;
  static time_t lasttime;
  double max_noise;
  static float *a;

  good_clipped = bad_clipped = good = bad = 0;
  for (rn = 0; rn < p->region_count; rn++) {
    TREALLOC (a, p->region_count);
    TREALLOC (a_dbl, p->region_count);
    for (i = 0; i < p->region[rn].count; i++)
      a[i] = a_dbl[i] = ss_cspline_eval (p->regionm[rn].spline, (double)i, p->regionm[rn].acc);
    max_noise = one_region_noise (a, p->region[rn].count, &maxidx);
    if (max_noise > 6.66) {
      static int count;
      if (count++ < 10)
        printf ("%s at %d: region %d at %d has extra noise %f\n", __FILE__, __LINE__, rn, maxidx, max_noise);
    }

    if (max_noise > 6.66)
    {
      int pn, j, showstart;
      float v[SNAPLEN];
      char msg[256];
      for (pn = 0; pn < p->peak_count; pn++)
        if (p->peak[pn].loc > p->region[rn].sample + maxidx)
          break;
      showstart = p->peak[pn].loc - p->region[rn].sample - SNAPLEN/2;
      snprintf (msg, 256, "region %d at %d for %d, %sclipped", rn, showstart, maxidx, p->region[rn].clipped ? "" : "not ");
      window ();
      clear ();
      for (j = 0; j < SNAPLEN; j++)
        v[j] = ss_cspline_eval (p->regionm[rn].spline, (double)showstart + j, p->regionm[rn].acc);
      draw_snap (v, "black");
      for (j = 0; j < SNAPLEN; j++)
        v[j] /= 2;
      draw_snap (v, "red");
      draw_text (msg);
      if (!p->region[rn].clipped)
        if (show () == 'q')
          exit (0);
    }
    
    if (max_noise <= 6.66) {
      good++;
      if (p->region[rn].clipped)
        good_clipped++;
    }
    else {
      bad++;
      if (p->region[rn].clipped)
        bad_clipped++;
    }
    if ((now = time (0)) > lasttime) {
      printf ("%d bad, %d good, %d total (of %d), %d%% bad, good_clipped: %d, bad_clipped: %d\n",
              bad, good, bad + good, p->region_count, (bad * 100) / (bad + good), good_clipped, bad_clipped);
      lasttime = now;
    }
  }
}

void
mem_use (char *file, int line)
{
/* removed from Mac OSX version */
#ifndef __APPLE__     
  struct mallinfo m;
  static struct mallinfo prev;
  m = mallinfo ();
  
  printf ("%s line %d: used: %10d, unused: %10d, mmapped: %10d, total: %10.0f, delta: %10.0f\n",
          file, line, m.uordblks, m.fordblks, m.hblkhd,
          (double)m.arena + m.hblkhd,
          ((double)m.arena + m.hblkhd) - ((double)prev.arena + prev.hblkhd));
  prev = m;
#endif
}

static int peakspan_threshold = 40;
static Center *global_center;

static int spread (Center *ctr, int *p0)
{
  int cz = ctr->peak_count - 1;
  int ctr_spread = ctr->peak[cz].loc - ctr->peak[0].loc;
  static int done[40];
  int n = ctr - global_center;
  int intraspike_peakspan_threshold = 14;
  
  *p0 = 0;
  if (ctr_spread > intraspike_peakspan_threshold) {
    float a = fabs (ctr->raw[ctr->peak[0].loc] - mean);
    float b = fabs (ctr->raw[ctr->peak[cz].loc] - mean);

    if (a > b) {
      cz--;
      if (0)
        if (!done[n])
          printf ("drop last peak of %2d\n", n);
    }
    else if (b > a) {
      *p0 = 1;
      if (0)
        if (!done[n])
          printf ("drop first peak of %2d\n", n);
    }
    ctr_spread = ctr->peak[cz].loc - ctr->peak[*p0].loc;
  }
  //  done[n] = 1;
  return ctr_spread;
}

static void
ctr_alignment_range (Center *ref, Center *ctr, int *minp, int *maxp)
{
  int c0, r0;
  int ctr_spread = spread (ctr, &c0);
  int ref_spread = spread (ref, &r0);
  int offset_to_align_peak_0 = -(ctr->peak[c0].loc - ref->peak[r0].loc);
  int threshold = 26;
  int extra, diff;


  if (ctr_spread > threshold || ref_spread > threshold)
    *minp = *maxp = 0;
  else if (ctr_spread >= ref_spread) {
    extra = threshold - ctr_spread;
    diff = ctr_spread - ref_spread;
    *minp = offset_to_align_peak_0 - diff - extra - 1;
    *maxp = offset_to_align_peak_0 + extra + 1;
  }
  else {
    extra = threshold - ref_spread;
    diff = ref_spread - ctr_spread;
    *minp = offset_to_align_peak_0 - extra - 1;
    *maxp = offset_to_align_peak_0 + diff + extra + 1;
  }
}

static double
distance_at_offset (double offset, void *p, int debug)
{
  Center **cp = (Center **)p;
  Center *ref = cp[0];
  Center *ctr = cp[1];
  float rw[SNAPLEN], cw[SNAPLEN];
  double sum, tmp;
  int i;

  TMEMCPY (rw, ref->raw, SNAPLEN);
  whiten_snap (rw);
  for (i = 0; i < SNAPLEN; i++)
    cw[i] = ss_cspline_eval (ctr->spline, i - offset, ctr->acc);
  whiten_snap (cw);
  sum = 0;
  for (i = 0; i < SNAPLEN; i++)
    tmp = cw[i] - rw[i], sum += tmp * tmp;
  return sqrt (sum);
}

static int
global_min (Center **ctr_list, int min, int max)
{
  int n, n_at_min = 0;
  double min_distance = HUGE_VAL, d;

  for (n = min; n <= max; n++)
    if ((d = distance_at_offset ((double)n, ctr_list, 0)) < min_distance) {
      min_distance = d;
      n_at_min = n;
    }
  return n_at_min;
}

static int
ctr_align_2 (Center *ref, Center *ctr, double *dp, double *offset, double *dref, double *dctr)
{
  Center *ctr_list[2];
  int min, max, status, i, startx;
  float v[SNAPLEN];
  double sum;
  int debug = ref - global_center == 13 && ctr - global_center == 29;
  //int debug = 0;
  
  ctr_list[0] = ref;
  ctr_list[1] = ctr;
  (ctr->peak_count >= 1 && ctr->peak_count <= 3) || DIE;
  //  printf ("%s line %d: %d %d\n", __FILE__, __LINE__, ctr->peak[0].loc, ctr->peak[ctr->peak_count-1].loc);
  ctr_alignment_range (ref, ctr, &min, &max);
  if (min == max) {
    if (debug)
      printf ("%ld %ld: min == max\n", (long)(ref - global_center), (long)(ctr - global_center));
    return 0;
  }
  startx = global_min (ctr_list, min, max);
  if (startx == min)
    status = -1;
  else if (startx == max)
    status = 1;
  else
    *offset = fmin2 (distance_at_offset, ctr_list, (double)startx, 1, (double)min, (double)max, &status);
  if (debug) {
    int n;
    printf ("%ld %ld: min %d, max %d, status = %d\n", (long)(ref - global_center), (long)(ctr - global_center), min, max, status);
    for (n = min; n <= max; n++)
      printf ("%3d: %f\n", n, distance_at_offset ((double)n, ctr_list, 0));
  }
  if (status)
    return 0;

  *dp = fmin_y ();

  TMEMCPY (v, ref->raw, SNAPLEN);
  whiten_snap (v);
  sum = 0;
  for (i = 0; i < SNAPLEN; i++)
    sum += v[i] * v[i];
  *dref = sqrt (sum);

  for (i = 0; i < SNAPLEN; i++)
    v[i] = ss_cspline_eval (ctr->spline, i - *offset, ctr->acc);
  whiten_snap (v);
  sum = 0;
  for (i = 0; i < SNAPLEN; i++)
    sum += v[i] * v[i];
  *dctr = sqrt (sum);
  return 1;
}

typedef struct {double left, right, d;} PeakSpan;

static PeakSpan
update_peak_span (PeakSpan *p, Center *c, double shift)
{
  PeakSpan ps;
  int c0;
  int cspread = spread (c, &c0);
  int left = c->peak[c0].loc;
  int right = left + cspread;

  if (p == 0) {
    ps.left = left + shift;
    ps.right = right + shift;
  }
  else {
    ps = *p;
    NEWMIN (ps.left, left + shift);
    NEWMAX (ps.right, right + shift);
  }
  ps.d = ps.right - ps.left;
  return ps;
}

static void
ctr_show_centers (DecomposeData *dd)
{
  size_t vcnt = (size_t)dd->center_count;
  int n;
  static int c;
  char shw[vcnt];
  char msg[256];
  static float vm[SNAPLEN];;
  
  dd->center_count >= 0 || DIE;
  if (vm[0] == 0)
    for (n = 0; n < SNAPLEN; n++)
      vm[n] = mean;

  memset (shw, 0, sizeof shw);
  window ();
  while (c != 'c' && c != 'n') {
    maybe_flip (shw, vcnt, c);
    clear ();
    draw_snap (vm, "black");
    for (n = 0; n < (int)vcnt; n++)
      if (shw[n]) {
        draw_snap (dd->center[n].raw, "black");
        if (n < 255)
          snprintf (msg + n, (size_t)(256 - n), "%c", ltr (n));
      }
      else if (n < 255) snprintf (msg + n, (size_t)(256 - n), "%c", ' ');
    draw_text (msg);
    c = show ();
    if (c == 'q')
      exit (0);
  }
  wclose ();
  if (c == 'n')
    c = 0;
}

double
next_largest (Edge *edge, size_t *eorder, int scnt, double weight_threshold)
{
  int n;
  double max;

  max = 0;
  for (n = 0; n < scnt; n++)
    if (edge[eorder[n]].weight < weight_threshold && edge[eorder[n]].weight > max)
      max = edge[eorder[n]].weight;
  return max;
}

void
write_dot_file (DecomposeData *dd, double *shift, Edge *edge, size_t *eorder, int scnt)
{
  FILE *f;
  int n;

  (f = fopen (change_filetype (file_name, ".dot"), "w")) || DIE;
  fprintf (f,
           "digraph align\n"
           "{\n"
           "    label=\"%s\";\n"
           "    rotate=90;\n"
           "    center=1;\n"
           "    size=\"8.5,11\";\n"
           "    page=\"8.5,11\";\n"
           "    rankdir=LR;\n"
           "    node [shape=circle];\n", change_filetype (file_name, ""));

  for (n = 0; n < scnt; n++)
    if (shift[edge[eorder[n]].vb] != 0)
      fprintf (f, "              %2d -> %2d[label = %2.0f];\n",
               dd->map[edge[eorder[n]].vb].old_to_new + 1,             
               dd->map[edge[eorder[n]].va].old_to_new + 1,
               floor (edge[eorder[n]].weight + .5));
  fprintf (f, "}\n");
  fclose (f);
}

static void
mfvec (char *name, float *data, int count)
{
  FILE *f;
  int n;

  (f = fopen ("vars.m", "a")) || DIE;
  fprintf (f, "%s = [", name);
  for (n = 0; n < count; n++)
    fprintf (f, "%s%.9e", n ? " " : "", data[n]);
  fprintf (f, "];\n");
  fclose (f);
}

static double *
ctr_pair_offset (DecomposeData *dd, SnapList *cc)
{
  size_t vcnt = (size_t)dd->center_count;
  size_t ecnt = vcnt * (vcnt - 1) / 2;
  size_t scnt = vcnt - 1;
  int psn;
  size_t a, b, e, n;
  double d, da, db;
  static Edge *edge;
  size_t *eorder = 0;
  static double *shift;
  static PeakSpan *peak_span;
  static int *psidx;
  double weight_threshold;

  dd->center_count > 0 || DIE;
  global_center = dd->center;
  TREALLOC (shift, vcnt);
  TREALLOC (peak_span, vcnt);
  TREALLOC (psidx, vcnt);
  TREALLOC (edge, ecnt);
  e = 0;
  for (a = 0; a < vcnt; a++) {
    for (b = a + 1; b < vcnt; b++) {
      double cosD;
      double offset = 0;
      if (ctr_align_2 (&dd->center[a], &dd->center[b], &d, &offset, &da, &db)) {
        cosD = (da*da + db*db - d*d) / (2*da*db);
        if (cosD < -1 || cosD > 1) {
          printf ("da: %f, db: %f, d: %f, cosD: %f\n", da, db, d, cosD);
          fflush (stdout);
          exit (DIE);
        }
        edge[e].weight = acos (cosD) / M_PI * 180;
        edge[e].weight = d;
        edge[e].shift = offset;
      }
      else {
        edge[e].weight = 90;
        edge[e].shift = 0;
      }
      edge[e].va = a;
      edge[e].vb = b;
      e++;
    }
  }
  e == ecnt || DIE;
  if (0)
    for (n = 0; n < vcnt; n++)
      cc->snap[n]->noise = white_distance_from_zero_snaplen (dd->center[n].raw) > norm_distance_from_zero (dd->center[n].raw);
  memset (shift, 0, vcnt * sizeof *shift);
  TMALLOC (dd->map, dd->center_count);

  if (vcnt > 1)
  {
    int *map, *newmap, newidx;
    double max_span;
    TMALLOC (map, vcnt);
    TMALLOC (newmap, vcnt);

    free (eorder);
    eorder = amplitude_tree (edge, vcnt, center_by_size (cc));
    //    eorder = spanning_tree (edge, vcnt);

    max_span = peakspan_threshold + 1.0;
    weight_threshold = DBL_MAX;
    while (max_span > peakspan_threshold) {
      memset (shift, 0, vcnt * sizeof *shift);
      a = edge[eorder[0]].va;
      (int)a < dd->center_count || DIE;
      psidx[a] = (int)a;
      peak_span[a] = update_peak_span (0, dd->center + a, 0);
      max_span = peak_span[a].d;

      map[0] = (int)a;
      for (n = 0; n < scnt; n++) {
        double s;
        PeakSpan pstmp;
        eorder[n] < ecnt || DIE;
        a = edge[eorder[n]].va;
        b = edge[eorder[n]].vb;
        (a < vcnt && b < vcnt) || DIE;
        map[n+1] = (int)b;
        s = shift[a] + edge[eorder[n]].shift;
        psn = psidx[a];
        pstmp = update_peak_span (&peak_span[psn], dd->center + b, s);
        if (edge[eorder[n]].weight < weight_threshold) {
          shift[b] = s;
          psidx[b] = psn;
          peak_span[psn] = pstmp;
        }
        else {
          shift[b] == 0 || DIE;
          psidx[b] = (int)b;
          pstmp = peak_span[b] = update_peak_span (0, dd->center + b, 0);
        }
        if (pstmp.d > max_span)
          max_span = pstmp.d;
      }
      weight_threshold = next_largest (edge, eorder, scnt, weight_threshold);
    }
    for (newidx = psn = 0; psn < (int)vcnt; psn++)
      for (n = 0; n < vcnt; n++)
        if (psidx[map[n]] == psn)
          newmap[newidx++] = map[n];
    newidx == (int)vcnt || DIE;
    for (n = 0; n < vcnt; n++) {
      dd->map[n].new_to_old = newmap[n];
      dd->map[newmap[n]].old_to_new = (int)n;
    }
    free (map);
    free (newmap);
  }
  if (vcnt > 1) {
    for (n = 0; n < vcnt; n++) {
      dd->map[dd->map[n].new_to_old].old_to_new == (int)n || DIE;
      dd->map[dd->map[n].old_to_new].new_to_old == (int)n || DIE;
    }
  }
  else {
    dd->map[0].old_to_new = 0;
    dd->map[0].new_to_old = 0;
  }
  if (vcnt > 1)
    write_dot_file (dd, shift, edge, eorder, scnt);
  free (eorder);
  if (vcnt > 1)
  {
    double min_left = (SNAPLEN - peakspan_threshold) / 2 - 1;
    for (n = 0; n < vcnt; n++)
      if (peak_span[psidx[n]].left < min_left)
        shift[n] += min_left - peak_span[psidx[n]].left;
  }

  for (n = 0; n < vcnt; n++) {
    Center *ctr;
    int i;
    double s = shift[n];
    ctr = dd->center + n;
    for (i = 0; i < SNAPLEN; i++)
      ctr->raw[i] = ss_cspline_eval (ctr->spline, i - s, ctr->acc);
    if (cc) {
      TMEMCPY (cc->snap[n]->raw, ctr->raw, SNAPLEN);
      TMEMCPY (cc->snap[n]->val, ctr->raw, SNAPLEN);
      whiten_snap (cc->snap[n]->val);
    }
    ss_spline_free (ctr->spline);
    gsl_interp_accel_free(ctr->acc);
    ctr->spline = 0;
    ctr->acc = 0;
    region_ctr_init (ctr);
  }

  fflush (stdout);
  if (0) {
    printf ("centers after ctr_pair_offset (%c)\n", ltr (dd->center_count - 1)); fflush (stdout);
    if (1)ctr_show_centers (dd);
  }
  return shift;
}

static void
shift_spikedata (DecomposeData *dd, double *shift)
{
  int sn, rn, i, ccn;
  RegionData *p = dd->p;
  double start;
      
  rn = 0;
  for (sn = 0; sn < spike_count; sn++) {
    if (spikedata[sn].sample < (double)p->region[rn].sample)
      rn = 0;
    while (spikedata[sn].sample >= (double)p->region[rn].sample + p->region[rn].count
           && rn < p->region_count)
      rn++;
    (rn >= 0 && rn < p->region_count) || DIE;

    p->regionm[rn].spline || DIE;
    ccn = spikedata[sn].partition;
    (ccn >= 0 && ccn < dd->center_count) || DIE;
    spikedata[sn].sample -= shift[ccn];
    start = (spikedata[sn].sample - PRESAMPLES) - p->region[rn].sample;
    for (i = 0; i < SNAPLEN; i++)
      spike_waveform[sn][i] = ss_cspline_eval (p->regionm[rn].spline, start + i, p->regionm[rn].acc);
  }
}

static void
sort_centers (DecomposeData *dd, SnapList **ccp)
{
  SnapList *cc = *ccp, *cc_new;
  Center *center;
  int n, ccn;

  TMALLOC (center, dd->center_count);
  TCALLOC (cc_new, 1);
  for (n = 0; n < dd->center_count; n++) {
    center[n] = dd->center[dd->map[n].new_to_old];
    add_snap (cc_new, cc->snap[dd->map[n].new_to_old]);
  }
  free_snaplist (&cc);
  *ccp = cc_new;
  free (dd->center);
  dd->center = center;
  for (n = 0; n < spike_count; n++)
    if ((ccn = spikedata[n].partition) >= 0 && ccn < dd->center_count)
      spikedata[n].partition = dd->map[ccn].old_to_new;
  for (n = 0; n < dd->p->region_count; n++)
    if ((ccn = dd->p->regionm[n].cluster) >= 0 && ccn < dd->center_count)
      dd->p->regionm[n].cluster  = dd->map[ccn].old_to_new;
}

static void
write_noise (char *file_name, SnapList *noise)
{
  FILE *f;
  int n;

  (f = fopen (file_name, "wb")) || DIE;
  fwrite (&noise->count, sizeof noise->count, 1, f) == 1 || DIE;
  for (n = 0; n < noise->count; n++)
    fwrite (noise->snap[n]->raw, sizeof *noise->snap[n]->raw, SNAPLEN, f) == SNAPLEN || DIE;
  fclose (f);
}

static void
histogram (float *buf, int count)
{
  int n, lo_first, hi_first;
  static int *hist;

  TCALLOC (hist, 65536);

  for (n = 0; n < count; n++)
    buf[n] == (int)buf[n] || DIE;

  for (n = 0; n < count; n++)
    hist[(int)buf[n] + 32768]++;

  n =     0; while (hist[n] == 0) n++; lo_first = n;
  n = 65535; while (hist[n] == 0) n--; hi_first = n;

  printf ("lo_first: %d (%d), hi_first: %d (%d)\n", lo_first, hist[lo_first], hi_first, hist[hi_first]);
  printf ("lo_first: %d (%d), hi_first: %d (%d)\n", lo_first+1, hist[lo_first+1], hi_first-1, hist[hi_first-1]);

  if (+  hist[hi_first] > 1 && hist[hi_first - 1] < hist[hi_first] / 2
      && hist[lo_first] > 1 && hist[lo_first + 1] < hist[hi_first] / 2) {
    clip_hi = hi_first - 32768;
    clip_lo = lo_first - 32768;
  }
  else {
    clip_hi = CLIPVAL;
    clip_lo = -CLIPVAL;
  }
  free (hist);
  printf ("clip_lo: %d, clip_hi = %d\n", clip_lo, clip_hi);
}

static RegionData *
duplicate_regiondata (RegionData *p)
{
  RegionData *q;

  TCALLOC (q, 1);
  q->data = p->data;
  q->data_count = p->data_count;
  q->data_alloc = p->data_alloc;
  TMALLOC (q->region, q->region_count = q->region_alloc = p->region_count);
  memcpy (q->region, p->region, q->region_count * sizeof (Region));
  TMALLOC (q->peak, q->peak_count = q->peak_alloc = p->peak_count);
  memcpy (q->peak, p->peak, q->peak_count * sizeof (Peak));
  return q;
}

static int dbl_cmp_dn (const double *a, const double *b)
{
  return *a > *b ? -1 : (*a < *b); /* descending */
}
static int rafc;

static void
adjust_peaks (Region *r, Peak *pk, int pkcnt, short *data, Peak **plp, int *pl_lenp)
{
  static int peak_alloc;
  static Peak *pl;
  static double *amp, *pamp;
  int n_at_max = 0, count, first;

  *plp = 0;
  *pl_lenp = 0;
  if (pkcnt == 1)
    return;
  if (peak_alloc < pkcnt) {
    TREALLOC (pl, peak_alloc = pkcnt);
    TREALLOC (amp, pkcnt);
    TREALLOC (pamp, pkcnt);
  }
  {
    double max = 0;
    int n;
    for (n = 0; n < pkcnt; n++) {
      double ampval = fabs (data[pk[n].loc - r->sample] - mean);
      if (ampval > max) {
        max = ampval;
        n_at_max = n;
      }
    }
    if (max / detection_threshold < 5)
      return;
  }
  first = n_at_max;
  count = 1;
  if (first > 0
      && SGN (data[pk[first].loc - r->sample]) == -SGN (data[pk[first-1].loc - r->sample])
      && (pk[first].loc - pk[first-1].loc) / SF < 0.5e-3) {
    first--;
    count++;
  }
  if (n_at_max + 1 < pkcnt
      && SGN (data[pk[n_at_max].loc - r->sample]) == -SGN (data[pk[n_at_max+1].loc - r->sample])
      && (pk[n_at_max+1].loc - pk[n_at_max].loc) / SF < 0.5e-3)
    count++;
  if (count == pkcnt)
    return;
  memcpy (pl, pk + first, count * sizeof (Peak));
  pl[0].loc < r->sample + r->count || DIE;
  *plp = pl;
  *pl_lenp = count;
}

 void
adjust_peaks_0 (Region *r, Peak *pk, int pkcnt, short *data, Peak **plp, int *pl_lenp)
{
  static int peak_alloc;
  static Peak *pl;
  static double *amp, *pamp;
  int n, bigcnt, first = 0;
  double threshold;
  
  *plp = 0;
  *pl_lenp = 0;
  if (pkcnt == 1)
    return;
  if (peak_alloc < pkcnt) {
    TREALLOC (pl, peak_alloc = pkcnt);
    TREALLOC (amp, pkcnt);
    TREALLOC (pamp, pkcnt);
  }
  for (n = 0; n < pkcnt; n++) {
    double ampval = data[pk[n].loc - r->sample] - mean;
    pamp[n] = amp[n] = fabs (ampval);
  }
  qsort (amp, (size_t)pkcnt, sizeof *amp, (int(*)(const void*,const void*))dbl_cmp_dn);
  if (DEBUG && amp[0] > 10000) {
    printf ("waveform.tcl *.chan %d %d\n", r->sample, r->sample + r->count);
    for (n = 0; n < pkcnt; n++)
      printf ("%6d %10d\n", data[pk[n].loc - r->sample], pk[n].loc);
  }
  for (n = 0; n < pkcnt - 1; n++)
    if (amp[n+1] == 0 || amp[n] / amp[n+1] > 5)
      break;
  if (n == pkcnt - 1) {
    if (DEBUG && amp[0] > 10000) printf ("adjust_peaks: no gap\n");
    return;
  }
  bigcnt = n + 1;
  threshold = (amp[n] + amp[n+1]) / 2;
  if (bigcnt > 3) {
    if (DEBUG && amp[0] > 10000) printf ("adjust_peaks: bigcnt > 3: %d\n", bigcnt);
    return;
  }

  for (n = 0; n < pkcnt; n++)
    if (pamp[n] > threshold)
      break;
  for (n = 1; n < bigcnt; n++)
    if (pamp[n] < threshold) {
      if (DEBUG && amp[0] > 10000) printf ("adjust_peaks: pamp[n] < threshold (%f < %f) %d\n",
                                           pamp[n], threshold, n);
      return;                   /* return if the big peaks are not consecutive */
    }
  if (DEBUG && amp[0] > 10000) printf ("adjust_peaks: bigcnt %d\n", bigcnt);
  memcpy (pl, pk + first, bigcnt * sizeof (Peak));
  pl[0].loc < r->sample + r->count || DIE;
  *plp = pl;
  *pl_lenp = bigcnt;
}

static void
adjust_regions (Region *r, Peak *pk, int pkcnt, Peak *pl, int pl_len, RegionData *q)
{
  int sample, count;
  Region *qr0 = 0, *qr1, *qr2 = 0;

  if (q->region_alloc < q->region_count + 3)
    TREALLOC (q->region, q->region_alloc = (q->region_count + 3) * 2);
  if (pl == 0) {
    memcpy (q->peak + q->peak_count, pk, pkcnt * sizeof (Peak));
    q->peak_count += pkcnt;
    q->region[q->region_count++] = *r;
    return;
  }
  sample = r->sample;
  count = r->count;
  qr0 = 0;
  if (pk[0].loc < sample + (pl[0].loc - pk[0].loc)) {
    qr0 = q->region + q->region_count;
    qr0->sample = r->sample;
    qr0->count = pl[0].loc - pk[0].loc;
    qr0->type = -3;
    qr0->clipped = 0;
    if (rafc) printf ("qr0 is region %d: %d to %d\n", q->region_count, qr0->sample, qr0->sample + qr0->count);
    q->region_count++;
    sample += qr0->count;
    count -= qr0->count;
  }
  qr1 = q->region + q->region_count;
  if (rafc) printf ("qr1 is region %d\n", q->region_count);
  qr1->sample = sample;
  if (qr0) qr1->sample == qr0->sample + qr0->count || DIE;
  qr1->type = 0;
  qr1->clipped = 0;
  q->region_count++;
  {
    int pe = pkcnt - 1;
    int ple = pl_len - 1;
    int tmpcnt = count - (pk[pe].loc - pl[ple].loc);
 
    if (pk[pe].loc >= qr1->sample + tmpcnt) {
      qr1->count = tmpcnt;
      sample += qr1->count;
      count  -= qr1->count;
      qr2 = q->region + q->region_count;
      if (rafc) printf ("qr2 is region %d\n", q->region_count);
      qr2->sample = sample;
      qr2->count = count;
      qr2->type = -3;
      qr2->clipped = 0;
      q->region_count++;
      qr1->sample + qr1->count == qr2->sample || DIE;
      qr2->sample + qr2->count == r->sample + r->count || DIE;
    }
    else {
      qr1->count = count;
      qr1->sample + qr1->count == r->sample + r->count || DIE;
    }
  }
  {
    int n = 0, e;
    while (pk[n].loc < qr1->sample) {
      if (rafc) printf ("peak at %d between %d and %d (qr0)\n",
                        pk[n].loc, qr0->sample, qr0->sample + qr0->count);
      n++;
    }
    memcpy (q->peak + q->peak_count, pk, n * sizeof (Peak));
    q->peak_count += n;

    if (rafc)
      for (n = 0; n < pl_len; n++)
        printf ("peak at %d between %d and %d (qr1)\n",
                pl[n].loc, qr1->sample, qr1->sample + qr1->count);
        
    memcpy (q->peak + q->peak_count, pl, pl_len * sizeof (Peak));
    q->peak_count += pl_len;
    n = 0;
    e = pkcnt - 1;
    sample = qr1->sample + qr1->count;
    while (pk[e - n].loc >= sample) {
      if (rafc) printf ("peak at %d between %d and %d (qr2)\n",
                        pk[e-n].loc, qr2->sample, qr2->sample + qr2->count);
      n++;
    }
    memcpy (q->peak + q->peak_count, pk + pkcnt - n, n * sizeof (Peak));
    q->peak_count += n;
  }
}

static void
region_adjust_for_clustering (RegionData *p, RegionData *q)
{
  int rn, pn, dn;

  q->peak_count = 0;
  q->region_count = 0;
  dn = pn = 0;
  for (rn = 0; rn < p->region_count; rn++) {
    while (pn < p->peak_count && p->peak[pn].loc < p->region[rn].sample)
      pn++;
    pn < p->peak_count || DIE;
    {
      Peak *pk = p->peak + pn;
      Region *r = p->region + rn;
      int next_region = r->sample + r->count;
      short *data = p->data + dn;
      int pkcnt = 0;
      Peak *pl;
      int pl_len, clipped = 0, n, end;

      while (pn + pkcnt < p->peak_count && pk[pkcnt].loc < next_region)
        pkcnt++;
      rafc = pn <= 49 && 49 < pn + pkcnt;
      end = r->count - AFTER;
      for (n = BEFORE; n <= end; n++)
        if (abs (data[n]) > CLIPVAL) {
          clipped = 1;
          break;
        }
      if (!clipped) adjust_peaks (r, pk, pkcnt, data, &pl, &pl_len); else pl = 0, pl_len = 0;
      pl_len < pkcnt || DIE;
      if (0)
      {
        int n;
        char *a, *b;
        a = strdup ("bigpeaks:");
        for (n = 0; n < pl_len; n++) {
          if (asprintf (&b, "%s %d", a, pl[n].loc) == -1) exit (1);
          free (a);
          a = b;
        }
        note (a);
        free (a);
      }
      q->peak_count + pkcnt <= p->peak_count || DIE;
      adjust_regions (r, pk, pkcnt, pl, pl_len, q);
      r->sample + r->count == q->region[q->region_count-1].sample + q->region[q->region_count-1].count || DIE;
      q->peak_count <= p->peak_count || DIE;
    }
    dn += p->region[rn].count;
  }
}

static void
free_q (RegionData *p)
{
  int n;

  if (p->regionm)
    for (n = 0; n < p->region_count; n++)
      ss_spline_free (p->regionm[n].spline);
  return;
  free (p->region);
  free (p->peak);
  free (p->regionm);
  free (p->peak_region);
}

static void
write_centers (SnapList *cc, char *file_name)
{
  FILE *f;
  int n;

  (f = fopen (change_filetype (file_name, ".ctr2"), "wb")) || DIE;
  fwrite (&cc->count, sizeof cc->count, 1, f) == 1 || DIE;
  for (n = 0; n < cc->count; n++) {
    fwrite (cc->snap[n]->raw, sizeof *cc->snap[n]->raw, SNAPLEN, f) == SNAPLEN || DIE;
    fwrite (cc->snap[n]->val, sizeof *cc->snap[n]->val, SNAPLEN, f) == SNAPLEN || DIE;
  }
  fclose (f);
}

SnapList *
read_centers (char *file_name)
{
  FILE *f;
  SnapList *cc;
  int ccn, count;

  (f = fopen (change_filetype (file_name, ".ctr2"), "rb")) || DIE;
  TCALLOC (cc, 1);
  fread (&count, sizeof count, 1, f) == 1 || DIE;

  for (ccn = 0; ccn < count; ccn++) {
    Snap *s;
    TCALLOC (s, 1);
    TMALLOC (s->raw, SNAPLEN);
    TMALLOC (s->val, SNAPLEN);
    fread (s->raw, sizeof *s->raw, SNAPLEN, f) == SNAPLEN || DIE;
    fread (s->val, sizeof *s->val, SNAPLEN, f) == SNAPLEN || DIE;
    add_snap (cc, s);
  }
  fclose (f);
  return cc;
}

SnapList *
read_noise (char *file_name)
{
  FILE *f;
  SnapList *cc;
  int ccn, count;

  (f = fopen (file_name, "rb")) || DIE;
  TCALLOC (cc, 1);
  fread (&count, sizeof count, 1, f) == 1 || DIE;
   
  for (ccn = 0; ccn < count; ccn++) {
    Snap *s;
    TCALLOC (s, 1);
    TMALLOC (s->raw, SNAPLEN);
    TMALLOC (s->val, SNAPLEN);
    fread (s->raw, sizeof *s->raw, SNAPLEN, f) == SNAPLEN || DIE;
    memcpy (s->val, s->raw, sizeof (SnapVals));
    whiten_snap (s->val);
    add_snap (cc, s);
  }
  fclose (f);
  return cc;
}

int
main (int argc, char **argv)
{
  SnapList *noise, *cc = 0;
  short ****wave_data_list;
  PlaneHdr *ph;
  Pctr **pctr_list;
  Projection **noise_projections;
  int detect = 0;

  if (argc > 1 && strcmp (argv[1], "version") == 0) {
    puts (PACKAGE_STRING);
    exit (0);
  }

  use_best_attempt_flag = false;

  if (strcmp (argv[argc-1], "detect") == 0) {
    detect = true;
    printf ("will quit after setting detection threshold\n");
    argc--;
  }

  if (strcmp (argv[argc-1], "use_best_attempt") == 0) {
    use_best_attempt_flag = true;
    printf ("using best attempt\n");
    argc--;
  }

  if (strcmp (argv[argc-1], "not_best_attempt") == 0) {
    use_best_attempt_flag = false;
    printf ("not using best attempt\n");
    argc--;
  }

  start_time = time (0);
  mem_use (__FILE__, __LINE__);
  lpfilter_coeff ();
  (void)(MPK_R < SNAPLEN / 2 || DIE);
  init_cdf();

  if (argc == 1) {
    printf ("usage: %s channel_file\n", argv[0]);
    exit (0);
  }
  file_name = argv[1];
  {
    FILE *f;
    (f = fopen (change_filetype (file_name, ".note"), "w")) || DIE;
    fclose (f);
  }
  note (PACKAGE_STRING);
  fillbuf (file_name);
  find_mean();
  note ("mean: %f\n", mean);
  histogram (buf, bufsamples);

  if (argc > 2)
    detection_threshold = atof (argv[2]);
  else 
    detection_threshold = find_detection_threshold (file_name);
  zap_to_threshold (detection_threshold);
  printf ("detection threshold: %.1f\n", detection_threshold);
  note ("detection threshold: %.1f\n", detection_threshold);

  note ("overall_sd : %.1f\n", find_overall_sd ());
  if (detect) return 0;
  printf ("zapped: %d, unzappped: %d, bufsamples: %d\n",
          bufsamples - unzapped_count, unzapped_count, bufsamples);
  if (unzapped_count > bufsamples - SNAPLEN * 3 && argc <= 2)
    write_status (file_name, cc, "signal doesn't cross threshold - no spikes found");

  find_unzapped_minmax();
  if (!read_white ()) {
    get_cvec ();
    find_unzapped_cov ();
  }
  (void)find_unzapped_sd ();
  note ("unzapped_sd: %.1f\n", unzapped_sd);
  find_white_cov ();
  whitener ();
  TAIL_THRESHOLD = sqrt (11.5*11.5 + (12731 / unzapped_sd) * (12731 / unzapped_sd));
  TAIL_THRESHOLD = 11.5;        /* debug */
  get_all_wmat ();
  write_white ();
  printf ("%79s\r", ""); fflush (stdout);
  //  printf ("SD: %.17f, cov_rep_root[0][0]: %.17f, cvec[0]: %.17f\n", unzapped_sd, cov_rep_root[0][0], cvec[0]);
  noise = find_noise_snaps();
  {
    RegionData *p = 0, *q = 0;
    DecomposeData *dd = 0;
    int rn = 0;

    if (argc == 4)
      rn = atoi (argv[3]);
    if (rn < 0) {
      debug_region = -rn;
      rn = 0;
    }
    if (rn)
      p = region_read_justone (change_filetype (file_name, ".regions"), rn);
    else
      p = region_read (change_filetype (file_name, ".regions"));

    if (p->region_count == 0) {
      int preexisting;
      TREALLOC (p->data, p->data_alloc = file_size / (int)sizeof *p->data);
      TREALLOC (p->peak, p->peak_alloc = file_size / (int)sizeof *p->peak);
      TREALLOC (p->region, p->region_alloc = file_size / (BEFORE + 1 + AFTER + 1));
      buf_end = find_mpk_spikes (p, buf, bufsamples, starting_sample, buf_start, preexisting = 0);
      while (fillbuf (file_name)) {
        zap_to_threshold (detection_threshold);
        printf ("%79s\r  chunk %d (%.0f%%): region_count %d, data_count %d\r",
                "", chunk, bytes_read * 100.0 / file_size, p->region_count, p->data_count);
        fflush (stdout);
        buf_end = find_mpk_spikes (p, buf, bufsamples, starting_sample, buf_start, preexisting = 0);
      }
      TREALLOC (p->data, p->data_count);
      TREALLOC (p->peak, p->peak_count);
      TREALLOC (p->region, p->region_count);
      region_classify (p);
      region_write (p, change_filetype (file_name, ".regions"));
      if (1) {
        region_adjust_for_clustering (p, q = duplicate_regiondata (p));
        if (DEBUG) exit (0);
        region_classify (q);
      }
    }
    if (p->region_count == 0)
      write_status (file_name, cc, "no spikes found");

    region_check (p) || DIE;
    if (q)
      region_check (q) || DIE;

    printf ("%d peaks, %d regions, %d samples\n", p->peak_count, p->region_count, p->data_count);

    region_spline_init (p);
    if (q) region_spline_init (q);

    TCALLOC (dd, 1);
    dd->p = p;

    if (q) dd->p = q;

    TREALLOC (spikedata, p->region_count);
    TREALLOC (spike_waveform, p->region_count);

    region_optavg (p);
    if (q) region_optavg (q);

    if (!(cc = region_read_centers (dd, file_name))) {
      cc = region_cluster (dd);
      if (dd->center_count == 0)
        write_status (file_name, cc, "no clusters found");

      if (q) spike_count = 0;
      region_write_centers (dd, file_name);
    }
    else if (argc != 4)
      read_spikedata (file_name, 0);

    dd->p = p;
    if (q) free_q (q);

    TREALLOC (spikedata, p->peak_count);
    TREALLOC (spike_waveform, p->peak_count);

    shift_spikedata (dd, ctr_pair_offset (dd, cc));
    sort_centers (dd, &cc);
    if (0)
    {
      int n;
      for (n = 0; n < dd->center_count; n++) {
        char *name;
        if (asprintf (&name, "center%d", n + 1) == -1) exit (1);
        mfvec (name, dd->center[n].raw, SNAPLEN);
        free (name);
      }
      exit (0);
    }

    dd->center_count == cc->count || DIE;
    dd->max_partition = cc->count - 1;
    decompose_unclassified_regions (dd);
    dd->max_partition++;
    dd->max_partition - cc->count <= CLIPPED_SPECIES || DIE;
    printf ("dd->max_partition: %d\n", dd->max_partition);
    write_centers (cc, file_name);
    write_spikedata (change_filetype (file_name, ".spd1"), 0);
    printf ("%d reassigned\n", revisit_assignments (dd, cc, noise));

    {
      int n, i;
      float max, v[SNAPLEN];
      for (n = 0; n < dd->center_count; n++) {
        max = 0;
        for (i = 0; i < SNAPLEN; i++) {
          v[i] = mean + (dd->center[n].raw[i] - mean) / 2;
          if (fabs (dd->center[n].raw[i] - mean) > max)
            max = fabs (dd->center[n].raw[i] - mean);
        }
        printf ("ccn %2d: %9.6f from zero, %12.6f max, %9.6f half, %d peaks %d %d\n", n,
                white_distance_from_zero_snaplen (dd->center[n].raw),
                max,
                white_distance_from_zero_snaplen (v),
                dd->center[n].peak_count,
                dd->center[n].peak[0].lfar,
                dd->center[n].peak[0].rfar
                );
      }                
      region_unclassified (p, cc->count);
    }

    {
      int n, ccn;
      for (ccn = 0; ccn < cc->count; ccn++)
        cc->snap[ccn]->count = 0;
      for (n = 0; n < spike_count; n++) {
        if (spikedata[n].partition == -2)
          spikedata[n].partition = dd->max_partition;
        if ((ccn = spikedata[n].partition) >= 0 && ccn < cc->count)
          cc->snap[ccn]->count++;
      }
    }    
    write_spikedata (change_filetype (file_name, ".spd2"), 1);
    write_noise (change_filetype (file_name, ".noise"), noise);

    gen_cluster_spikestart (cc);
    ph = get_planes (cc);
    pctr_list = get_pctr_list (cc, ph);
    noise_projections = get_noise_projections (noise, ph);
    add_mean_snap (cc, dd->max_partition + 1);
    {
      int ccn = cc->count - 1, n;
      cc->snap[ccn]->count = 0;
      for (n = 0; n < spike_count; n++)
        if (spikedata[n].partition == ccn)
          cc->snap[ccn]->count++;
    }
    wave_data_list = alloc_wave_data_list (cc);
    add_spike_waveforms_3 (cc, wave_data_list);
    write_pdf (cc, wave_data_list, noise_projections, noise, pctr_list, ph);
    write_wdt (file_name, cc, wave_data_list);
    write_path ();
    write_status (file_name, cc, "");
    exit (0);
  }
  return 0;
}
