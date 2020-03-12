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

/* region.c */


//#define JUSTONE 220982520
#define JUSTONER 38925
//#define SAMPLOC 162994168
#define SAMPLOC -1

extern int debug_region;
static char *reason;
extern int clip_hi, clip_lo;
 
static int mtotal[100];
int mtcnt;
#include "nde.h"
#include "region.h"
#include <limits.h>
#include <time.h>

#ifdef __APPLE__
#include <stdlib.h>
#else
#include <malloc.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_erf.h>
#include <string.h>
#include <stdarg.h>

#ifndef __APPLE__ 
extern size_t malloc_usable_size (void *);
#else
size_t
malloc_usable_size (void *p) 
{
  return 1;
}
#endif

#ifdef S_SPLINT_S
int asprintf (char **strp, const char *fmt, ...); 
#endif

int malloc_debug;

void
myfree (int line, void *p)
{
  if (malloc_debug && p)
    fprintf (stderr, "%d %ld %d\n", line, (long)p, 0);
  free (p);
}

//#define free(x) myfree(__LINE__, x)

#define STATIC static
#define BIGFIRST 1
static int do_scaled_distance = 0;

static FILE *spikedata_tmpfile;
static int still_unclassified;

static time_t last_time, now;
static int maxregion;

#ifdef JUSTONE
static int justone;
#endif

static int global_rn;
static double drstart;
static int *loose_matches;

void
note (char *fmt, ...)
{
  FILE *f;
  va_list args;

  (f = fopen (change_filetype (file_name, ".note"), "a")) || DIE;
  va_start (args, fmt);
  vfprintf (f, fmt, args);
  fprintf (f, "\n");
  va_end (args);
  fclose (f);
}

int
region_add_data (RegionData *p, int b, int lfar, int lnear, int rnear, int rfar, float *buf, int bufsamples, int starting_sample, int preexisting)
{
  int next, first, last, pn, rn, first_idx = 0;

  first = starting_sample + b - BEFORE;
  last = starting_sample + b + AFTER;
  
  if (preexisting == 0) {
    if (first + preexisting * 50 < starting_sample) {
      if (0)
        printf ("file %s line %d function %s: peak data overlaps beginning of buffer by %d\n",
                __FILE__, __LINE__, __FUNCTION__, starting_sample - first);
      return 0;
    }
    if (last - preexisting * 50 >= starting_sample + bufsamples) {
      if (0)
        printf ("file %s line %d function %s: peak data overlaps end of buffer by %d\n",
                __FILE__, __LINE__, __FUNCTION__, last - (starting_sample + bufsamples));
      return 0;
    }
  }

  if (!preexisting) {
    int n, no_clip_cnt;

    first_idx = first - starting_sample;
    no_clip_cnt = 0;

    for (n = b; n > 0; n--) {
      if (fabs (buf[n]) >= CLIPVAL)
        no_clip_cnt = 0;
      else
        no_clip_cnt++;
      if (no_clip_cnt >= BEFORE && n <= first_idx)
        break;
    }
    first_idx = n;
    first = starting_sample + first_idx;
    if (no_clip_cnt < BEFORE) {
      printf ("\n");
      printf ("%s line %d: no_clip_cnt %d < %d, sample %d", __FILE__, __LINE__, no_clip_cnt, BEFORE, starting_sample);
      printf ("\n");
      note ("%s line %d: no_clip_cnt %d < %d, sample %d", __FILE__, __LINE__, no_clip_cnt, BEFORE, starting_sample);
    }
  }

  rn = p->region_count - 1;
  if (rn < 0 || (first > (next = p->region[rn].sample + p->region[rn].count) && !preexisting))
  {
    rn = p->region_count++;
    if (p->region_count > p->region_alloc)
      TREALLOC (p->region, p->region_alloc += 1024*1024 / sizeof *p->region);
    next = p->region[rn].sample = first;
    p->region[rn].count = 0;
    if (preexisting) {
      p->region[rn].sample = starting_sample;
      p->region[rn].count = bufsamples;
      p->fdata = buf;
    }
  }

  pn = p->peak_count++;
  if (p->peak_count > p->peak_alloc)
    TREALLOC (p->peak, p->peak_alloc += 1024*1024 / sizeof *p->peak);
  p->peak[pn].loc = starting_sample + b;
  p->peak[pn].lfar = lfar;
  p->peak[pn].lnear = lnear;
  p->peak[pn].rfar = rfar;
  p->peak[pn].rnear = rnear;
            
  if (preexisting)
    return 1;
  {
    int overlap = next - first;
    int new_data, last_idx;
    int n, no_clip_cnt;
    float *from;
    short *to;

    no_clip_cnt = 0;
    last_idx = b + AFTER;
    for (n = b + 1; n < bufsamples; n++) {
      if (abs (buf[n]) >= CLIPVAL)
        no_clip_cnt = 0;
      else
        no_clip_cnt++;
      if (no_clip_cnt >= AFTER && n >= last_idx)
        break;
    }
    last_idx = n;
    if (last_idx >= bufsamples) {
      last_idx = bufsamples - 1;
      printf ("\n");
      printf ("%s line %d: no_clip_cnt %d < %d, sample %d", __FILE__, __LINE__, no_clip_cnt, AFTER, starting_sample + last_idx);
      printf ("\n");
      note ("%s line %d: no_clip_cnt %d < %d, sample %d", __FILE__, __LINE__, no_clip_cnt, AFTER, starting_sample + last_idx);
    }
    last = starting_sample + last_idx;
    from = buf + (first_idx + overlap);
    new_data = last - first + 1 - overlap;
    first_idx + overlap + new_data <= bufsamples || DIE;

    while (p->data_count + new_data > p->data_alloc)
      TREALLOC (p->data, p->data_alloc += 1024*1024 / sizeof *p->data);

    to = p->data + p->data_count;
    for (n = 0; n < new_data; n++) {
      to[n] = (short)from[n];
      (float)to[n] == from[n] || DIE;
    }
    //    memcpy (p->data + p->data_count, buf + (b - BEFORE + overlap), new_data * sizeof *p->data);

    p->data_count += new_data;
    p->region[rn].count += new_data;
    return 1;
  }
  exit (DIE);
  return 0;
}

int
region_check (RegionData *p)
{
  int rn, pn, data_count;
  int tcnt[9];

  memset (tcnt, 0, sizeof (tcnt));
  (p->region[0].type + 3 >= 0 && p->region[0].type + 3 < 9) || DIE;
  tcnt[p->region[0].type + 3]++;
  rn = 0;
  data_count = p->region[rn].count;
  p->region_count >= 0 || DIE;
  if (p->regionm == 0) {
    int n;
    TCALLOC (p->regionm, (size_t)p->region_count);
    for (n = 0; n < p->region_count; n++)
      p->regionm[n].cluster = -1;
    TREALLOC (p->peak_region, p->peak_count);
  }
  if (p->region[0].type >= 6 || p->region[0].type < -3) {
    printf ("file %s line %d function %s: region %d is type %d\n",
            __FILE__, __LINE__, __FUNCTION__, rn, p->region[rn].type);
    return 0;
  }
  if (p->peak_count && (p->peak[0].loc < 0 || p->peak[0].loc >= p->region[0].sample + p->region[0].count)) {
    printf ("file %s line %d function %s: peak %d is not in region %d\n",
            __FILE__, __LINE__, __FUNCTION__, 0, 0);
    return 0;
  }
  if (p->region[0].count > maxregion)
    maxregion = p->region[0].count;
  if (p->peak_count > 0)
    p->peak[p->peak_count-1].loc < p->region[p->region_count-1].sample + p->region[p->region_count-1].count || DIE;
  for (pn = 0; pn < p->peak_count; pn++) {
    p->peak_region[pn] = rn;
    while (p->peak[pn].loc >= p->region[rn].sample + p->region[rn].count) {
      rn++;
      p->peak_region[pn] = rn;
      if (p->region[rn].type >= 6 || p->region[rn].type < -3) {
        printf ("file %s line %d function %s: region %d is type %d\n",
              __FILE__, __LINE__, __FUNCTION__, rn, p->region[rn].type);
        return 0;
      }
      if (p->peak[pn].loc < p->region[rn].sample
          || (p->peak[pn].loc >= p->region[rn].sample + p->region[rn].count
              && p->region[rn].type != -3))
        {
          printf ("file %s line %d function %s: peak %d is not in region %d\n",
                  __FILE__, __LINE__, __FUNCTION__, pn, rn);
          printf ("pn: %d, p->peak[pn].loc: %d, p->region[rn].sample: %d, p->region[rn].count: %d\n",
                  pn, p->peak[pn].loc, p->region[rn].sample, p->region[rn].count);
          return 0;
        }
      tcnt[p->region[rn].type + 3]++;
      p->regionm[rn].pn = pn;

      //      printf ("%s line %d: region count %d, maxregion %d\n", __FILE__, __LINE__, p->region[rn].count, maxregion);
      if (p->region[rn].count > maxregion)
        maxregion = p->region[rn].count;
      data_count += p->region[rn].count;
    }
    if (rn >= p->region_count) {
      printf ("region %d, region_count %d, peak loc: %d, %d peaks\n",
              rn, p->region_count, p->peak[pn].loc, p->peak_count);
      printf ("file %s line %d function %s: peak %d is past last region\n",
              __FILE__, __LINE__, __FUNCTION__, pn);
      return 0;
    }
    if (p->fdata == 0 && p->region[rn].type != -3) {
      if (p->peak[pn].loc - BEFORE < p->region[rn].sample) {
        printf ("file %s line %d function %s: peak %d data starts before region %d\n",
                __FILE__, __LINE__, __FUNCTION__, pn, rn);
        return 0;
      }
      if (p->peak[pn].loc + AFTER >= p->region[rn].sample + p->region[rn].count && p->region[rn].type != -3) {
        printf ("file %s line %d function %s: peak %d data ends after region %d\n",
                __FILE__, __LINE__, __FUNCTION__, pn, rn);
        return 0;
      }
    }
  }
  if (p->fdata == 0) {
    if (rn != p->region_count - 1) {
      printf ("file %s line %d function %s: %d is last region used\n",
              __FILE__, __LINE__, __FUNCTION__, rn);
      return 0;
    }
    if (data_count != p->data_count) {
      printf ("file %s line %d function %s: sum of region data counts %d is not equal to total data count\n",
              __FILE__, __LINE__, __FUNCTION__, data_count);
      return 0;
    }
    if (0)
    {
      int n;
      for (n = 0; n < 6; n++)
        printf ("%5d type %2d\n", tcnt[n], n - 3);
    }
  }
  return 1;
}

void
region_write (RegionData *p, char *filename)
{
  FILE *f;

  (f = fopen (filename, "wb")) || DIE;
  fwrite (&p->peak_count, sizeof p->peak_count, 1, f) == 1 || DIE;
  p->peak_count >= 0 || DIE;
  fwrite (p->peak, sizeof *p->peak, (size_t)p->peak_count, f) == (size_t)p->peak_count || DIE;
  p->region_count >= 0 || DIE;
  fwrite (&p->region_count, sizeof p->region_count, 1, f) == 1 || DIE;
  fwrite (p->region, sizeof *p->region, (size_t)p->region_count, f) == (size_t)p->region_count || DIE;
  fwrite (&p->data_count, sizeof p->data_count, 1, f) == 1 || DIE;
  p->data_count >= 0 || DIE;
  fwrite (p->data, sizeof *p->data, (size_t)p->data_count, f) == (size_t)p->data_count || DIE;
  fclose (f);
}

static void
mclr (void)
{
  FILE *f;
  char *filename = "vars.m";

  if (access (filename, F_OK) == 0) {
    (f = fopen (filename, "w")) || DIE;
    fclose (f);
  }
}

static int justone_region;

RegionData *
region_read_justone (char *filename, int rn)
{
  FILE *f;
  RegionData *p;
  Peak *peak;
  Region *region;
  int peak_count, region_count, data_count, n, first_peak, peak_end, skip = 0;

  justone_region = rn;
  mclr ();
  TCALLOC (p, 1);
  if ((f = fopen (filename, "rb")) == 0)
    return p;

  fread (&peak_count, sizeof peak_count, 1, f) == 1 || DIE;
  peak_count >= 0 || DIE;
  TMALLOC (peak, peak_count);
  fread (peak, sizeof *peak, (size_t)peak_count, f) == (size_t)peak_count || DIE;

  fread (&region_count, sizeof region_count, 1, f) == 1 || DIE;
  (rn >= 0 && rn < region_count) || DIE;
  TMALLOC (region, rn);
  fread (region, sizeof *region, (size_t)rn, f) == (size_t)rn || DIE;
  for (n = 0; n < rn; n++)
    skip += region[n].count;
  TMALLOC (p->region, p->region_alloc = 1);
  fread (p->region, sizeof *p->region, 1, f) == 1 || DIE;
  fseek (f, (region_count - rn - 1) * (long)sizeof *region, SEEK_CUR);
  p->region_count = 1;

  for (n = 0; n < peak_count; n++)
    if (peak[n].loc >= p->region->sample)
      break;
  first_peak = n;
  for ( ; n < peak_count; n++)
    if (peak[n].loc >= p->region->sample + p->region->count)
      break;
  peak_end = n;
  p->peak_count = peak_end - first_peak;
  TMALLOC (p->peak, p->peak_alloc = p->peak_count);
  TMEMCPY (p->peak, peak + first_peak, p->peak_count);

  fread (&data_count, sizeof data_count, 1, f) == 1 || DIE;
  fseek (f, skip * (long)sizeof *p->data, SEEK_CUR);
  p->data_count = p->region->count;
  p->data_count >= 0 || DIE;
  TMALLOC (p->data, p->data_alloc = p->data_count);
  fread (p->data, sizeof *p->data, (size_t)p->data_count, f) == (size_t)p->data_count || DIE;
  fclose (f);
  free (region);
  free (peak);
  return p;
}

RegionData *
region_read (char *filename)
{
  FILE *f;
  RegionData *p;

  mclr ();
  TCALLOC (p, 1);
  if ((f = fopen (filename, "rb")) == 0)
    return p;

  fread (&p->peak_count, sizeof p->peak_count, 1, f) == 1 || DIE;
  p->peak_count >= 0 || DIE;
  TMALLOC (p->peak, p->peak_alloc = p->peak_count);
  fread (p->peak, sizeof *p->peak, (size_t)p->peak_count, f) == (size_t)p->peak_count || DIE;

  fread (&p->region_count, sizeof p->region_count, 1, f) == 1 || DIE;
  p->region_count >= 0 || DIE;
  TMALLOC (p->region, p->region_alloc = p->region_count);
  fread (p->region, sizeof *p->region, (size_t)p->region_count, f) == (size_t)p->region_count || DIE;

  fread (&p->data_count, sizeof p->data_count, 1, f) == 1 || DIE;
  p->data_count >= 0 || DIE;
  TMALLOC (p->data, p->data_alloc = p->data_count);
  fread (p->data, sizeof *p->data, (size_t)p->data_count, f) == (size_t)p->data_count || DIE;
  fclose (f);
  printf ("%s line %d: region 1 is %d samples starting at %d\n", __FILE__, __LINE__, p->region[1].count, p->region[1].sample);
  return p;
}

typedef struct {
  float *v;
  int vcnt;
  int *clip;
  int clipcnt;
} ClipParams;
  
double
clipval_noise (const gsl_vector *clipvals, void *clip_params)
{
  ClipParams *cp = (ClipParams *)clip_params;
  int maxidx, i;

  (cp->clipcnt >= 0 && clipvals->size == (size_t)cp->clipcnt) || DIE;
  for (i = 0; i < cp->clipcnt; i++)
    cp->v[cp->clip[i]] = gsl_vector_get (clipvals, (size_t)i);
  return one_region_noise (cp->v, cp->vcnt, &maxidx);
}

static void
pvec (char *lbl, gsl_vector_complex *a, FILE *f)
{
  size_t k;

  fprintf (f, "%s", lbl);
  for (k = 0; k < a->size; k++) {
    gsl_complex av = gsl_vector_complex_get (a, k);
    fprintf (f, "%s%.17e+%.17ei", k ? "; " : "", GSL_REAL (av), GSL_IMAG (av));
  }
}

static void
pmat (char *lbl, gsl_matrix_complex *a, FILE *f)
{
  size_t l, k;

  fprintf (f, "%s", lbl);
  for (l = 0; l < a->size1; l++) {
    for (k = 0; k < a->size2; k++) {
      gsl_complex av = gsl_matrix_complex_get (a, l, k);
      fprintf (f, "%s%.17e+%.17ei", k ? ", " : "", GSL_REAL (av), GSL_IMAG (av));
    }
    fprintf (f, ";\n");
  }
  fprintf (f, "\n");
}

void
to_octave (gsl_matrix_complex *T, gsl_vector_complex *b)
{
  FILE *f;

  (f = fopen ("A_k.m", "w")) || DIE;
  pmat ("A = [", T, f);
  fprintf (f, "];\n");
  pvec ("k = [", b, f);
  fprintf (f, "];\n");
  fclose (f);
}

static void
mvec (char *name, double *data, int count)
{
  FILE *f;
  int n;

  (f = fopen ("vars.m", "a")) || DIE;
  fprintf (f, "%s = [", name);
  for (n = 0; n < count; n++)
    fprintf (f, "%s%.16e", n ? " " : "", data[n]);
  fprintf (f, "];\n");
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

static void
mcvec (char *name, char *txt)
{
  FILE *f;

  (f = fopen ("vars.m", "a")) || DIE;
  fprintf (f, "%s = ", name);
  fprintf (f, "\"%s\"", txt);
  fprintf (f, ";\n");
  fclose (f);
}

static void
mivec (char *name, int *data, int count)
{
  FILE *f;
  int n;

  (f = fopen ("vars.m", "a")) || DIE;
  fprintf (f, "%s = [", name);
  for (n = 0; n < count; n++)
    fprintf (f, "%s%d", n ? " " : "", data[n]);
  fprintf (f, "];\n");
  fclose (f);
}

fftw_real *tdat0;

static double *
lowpass_f (Region *region, RegionM *regionm, double maxfreq)
{
  static fftw_real *tdat, *fdat;
  static int dat_alloc;
  int fftlen, pow2, n;
  static rfftw_plan fwd[32], rev[32];
  int passband, zerocount;

  (void)frexp ((double)region->count, &pow2);
  fftlen = (int)ldexp (1, pow2);
  (pow2 >= 0 && pow2 < 32) || DIE;
  if (!fwd[pow2]) {
    rev[pow2] = rfftw_create_plan (fftlen, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
    fwd[pow2] = rfftw_create_plan (fftlen, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
  }
  if (dat_alloc < fftlen) {
    TREALLOC (tdat, fftlen);
    TREALLOC (tdat0, fftlen);
    TREALLOC (fdat, fftlen);
    dat_alloc = fftlen;
  }
  passband = (int)((maxfreq * fftlen) / SF);
  region->count <= fftlen || DIE;
  for (n = 0; n <  region->count; n++)
    tdat[n] = ss_cspline_eval (regionm->spline, (double)n, regionm->acc);
  for (n = region->count; n < fftlen; n++)
    tdat[n] = mean;
  memcpy (tdat0, tdat, fftlen * sizeof *tdat);
  rfftw_one (fwd[pow2], tdat, fdat);
  zerocount = fftlen - (2 * passband + 1);
  (passband + 1 >= 0 && passband + 1 < fftlen && passband + zerocount >= 0 && passband + zerocount < fftlen) || DIE;
  for (n = passband + 1; n <= passband + zerocount; n++)
    fdat[n] = 0.0;
  rfftw_one (rev[pow2], fdat, tdat);
  for (n = 0; n <  region->count; n++)
    tdat[n] /= fftlen;
  return tdat;
}

static double *
lowpass (Region *region, RegionM *regionm)
{
  return lowpass_f (region, regionm, 10000);
}

static void
fdata_spline_init (RegionData *p)
{
  int count = p->region[0].count;

  (p->regionm[0].acc || p->regionm[0].spline) && DIE;
  count <= maxregion || DIE;
  p->regionm[0].spline = (gsl_spline *)cspline_new (p->fdata, 0, count);
}

static void
region_one_spline_free (RegionM *regionm)
{
  regionm->acc == 0 || DIE;
  regionm->spline->interp == 0 || DIE;
  ss_spline_free (regionm->spline);
  free (regionm->clip);
}

gsl_spline *
peak_spline_clip (float *v, int rcount, char *clip)
{
  static double *x, *y;
  int n, xyn;
  static gsl_spline spline;

  if (!x) {
    TMALLOC (x, maxregion);
    TMALLOC (y, maxregion);
  }
  
  for (xyn = n = 0; n < rcount; n++)
    if (!clip[n]) {
      x[xyn] = (double)n;
      y[xyn] = (double)v[n];
      xyn++;
    }
  spline.x = x;
  spline.y = y;
  spline.size = xyn;
  spline.interp = gsl_interp_alloc (gsl_interp_cspline, xyn);
  gsl_interp_init (spline.interp, x, y, xyn);

  for (n = 0; n < rcount; n++)
    if (clip[n])
      v[n] = gsl_spline_eval (&spline, n, 0);

  return &spline;
}

static int undershoot_debug;

static int
fix_undershoots (float *v, int rcount, char *clip)
{
  int n, i, cnt;

  cnt = n = 0;
  while (1) {
    for ( ; n < rcount; n++)
      if (clip[n] && v[n] < clip_hi && v[n] > clip_lo)
        break;
    if (n == rcount)
      return cnt;
    for (i = n + 1; i < rcount; i++)
      if (!(clip[i] && v[i] < clip_hi && v[i] > clip_lo))
        break;
    (n > 0 && i < rcount) || DIE;

    cnt += i - n;
    if (clip[n-1] && clip[i]) {
      printf ("clip: %d %d\n", clip[n-1], clip[i]);
      undershoot_debug = 1;
      return 0;
    }

    if (!clip[i]) {
      v[i-1] = (clip[i-1] > 0) ? clip_hi : clip_lo;
      clip[i-1] = 0;
    }

    if (clip[n] && !clip[n-1]) {
      v[n] = (clip[n] > 0) ? clip_hi : clip_lo;
      clip[n] = 0;
    }

//    if ((clip[n-1] && !clip[i]) || (!clip[n-1] && !clip[i] && fabs (v[i]) < fabs (v[n-1]))) {
//      v[i-1] = (v[n-1] > 0) ? clip_hi : clip_lo;
//      clip[i-1] = 0;
//    }
//    else if ((!clip[n-1] && clip[i]) || (!clip[n-1] && !clip[i] && fabs (v[i]) > fabs (v[n-1]))) {
//      v[n] = (v[i] > 0) ? clip_hi : clip_lo;
//      clip[n] = 0;
//    }

//    if (!clip[i]) {
//      v[i-1] = (v[i-1] > 0) ? clip_hi : clip_lo;
//      clip[i-1] = 0;
//    }
//    else if (!clip[n-1]) {
//      v[n] = (v[n] > 0) ? clip_hi : clip_lo;
//      clip[n] = 0;
//    }

    n = i + 1;
  }
}

typedef struct {int amp; int loc;} Amp;
typedef struct {Amp amplist[5]; int ampcount; int peakloc; int start;} AmpData;

static int next_spike (float *buf, int end, int bufidx, AmpData *ampdata, double *xtp);
#define clip_pol(x) (x >= CLIPVAL ? 1 : (x <= -CLIPVAL ? -1 : 0))

static void
region_one_spline_init (Region *region, RegionM *regionm, /*@observer@*/ short *y_short)
{
  float *y_float = 0;
  int n;
  static char *clip;
  static int pass_max;
  int passcnt;
  //  float fix[12][maxregion];

  if (!clip) TMALLOC (clip, maxregion);

  region->clipped = 0;
  memset (clip, 0, region->count);
  for (n = 0; n < region->count; n++)
    if ((clip[n] = clip_pol (y_short[n])) != 0)
      region->clipped = 1;
  if (region->clipped) {
    double *tdat;
    AmpData ampdata;
    double xt;

    TMALLOC (y_float, region->count);
    for (n = 0; n < region->count; n++)
      y_float[n] = y_short[n];
    ampdata.peakloc = -SNAPLEN;
    ampdata.ampcount = 0;
    if (next_spike (y_float, region->count, 0, &ampdata, &xt)) {
      note ("doubly clipped spike at %9d", ampdata.amplist[1].loc + region->sample);
      free(y_float);
      y_float = 0;
    }
    else {
      regionm->acc = 0;
      regionm->spline = peak_spline_clip (y_float, region->count, clip);
      //      memcpy (fix[0], y_float, region->count * 4);
      passcnt = 0;
      while (fix_undershoots (y_float, region->count, clip)) {
        gsl_interp_free (regionm->spline->interp);
        regionm->spline = peak_spline_clip (y_float, region->count, clip);
        passcnt++;
        //      memcpy (fix[passcnt], y_float, region->count * 4);
        if (passcnt > pass_max) {
          pass_max = passcnt;
          printf ("pass_max: %d\n", pass_max);
          if (undershoot_debug || passcnt > 10) {
            undershoot_debug = 1;
            break;
          }
        }
      }
      y_float || DIE;
      if (0 && undershoot_debug) {
        FILE *f;
        char name[10];
        (f = fopen ("undershoots.chan", "wb")) || DIE;
        fwrite (y_short, sizeof (short), region->count, f) == region->count || DIE;
        fclose (f);
        printf ("passcnt: %d\n", passcnt);
        for (n = 0; n <= passcnt; n++) {
          sprintf (name, "fix%d", n);
          (f = fopen (name, "wb")) || DIE;
          //      fwrite (fix[n], sizeof (float), region->count, f) == region->count || DIE;
          fclose (f);
        }
        exit (DIE);
      }

      tdat = lowpass (region, regionm);
      y_float || DIE;
      gsl_interp_free (regionm->spline->interp);
      y_float || DIE;
      for (n = 0; n < region->count; n++)
        y_float[n] = tdat[n];
      regionm->spline = 0;
    }
  }
  else {
    if (0) {
      static fftw_real *tdat;
      TMALLOC (y_float, region->count);
      for (n = 0; n < region->count; n++)
        y_float[n] = y_short[n];
      regionm->acc = 0;
      regionm->spline = peak_spline_clip (y_float, region->count, clip);
      tdat = lowpass (region, regionm);
      for (n = 0; n < region->count; n++)
        y_float[n] = tdat[n];
    }
  }
  regionm->spline = (gsl_spline *)cspline_new (y_float, y_short, region->count);
}

void
region_spline_init (RegionData *p)
{
  int rn, dn = 0;

  printf ("%79s\r region_spline_init\r", ""); fflush (stdout);

#ifdef JUSTONE
  if (JUSTONER)
    justone = JUSTONER;
  else {
    justone = -1;
    for (rn = 0; rn < p->region_count; rn++)
      if (p->region[rn].sample <= JUSTONE && JUSTONE < p->region[rn].sample + p->region[rn].count) {
        justone = rn;
        break;
      }
  }
  if (justone == -1) {
    for (rn = 0; rn < p->region_count; rn++)
      if (p->region[rn].sample > JUSTONE) {
        printf ("region %d: %d to %d\n", rn - 1, p->region[rn-1].sample, p->region[rn-1].sample + p->region[rn-1].count);
        printf ("justone: %d\n", JUSTONE);
        printf ("region %d: %d to %d\n", rn, p->region[rn].sample, p->region[rn].sample + p->region[rn].count);
        exit (DIE);
      }
  }
  for (rn = 0; rn <= justone; rn++) {
#else
  for (rn = 0; rn < p->region_count; rn++) {
#endif
#ifdef JUSTONE
    if (rn >= justone)          /* debug */
#endif
    region_one_spline_init (p->region + rn, p->regionm + rn, p->data + dn);
    dn += p->region[rn].count;
  }
  printf ("%79s\r region_spline_init done\r", ""); fflush (stdout);
}

SnapList *
region_snaplist (RegionData *p)
{
  int rn, pn;
  SnapList *sl;
  Snap *s;

  TCALLOC (sl, 1);
  rn = 0;
  for (pn = 0; pn < p->peak_count; pn++) {
    if (p->peak[pn].loc > p->region[rn].sample + p->region[rn].count)
      rn++;
    TCALLOC (s, 1);
    s->peak = p->peak[pn].loc;
    s->region = rn;
    add_snap (sl, s);
  }
  return sl;
}

static double
spline_peak_1 (gsl_spline *spline, gsl_interp_accel *acc, int max_x, int xa, int xc, double yb)
{
  double x0, x1, x2, x3, y0, y1, y2, y3, rad;
  double A, B, C, D, a, b, c;

  x0 = (double)xa;
  x1 = xa + 1.0/3 * (xc - xa);
  x2 = xa + 2.0/3 * (xc - xa);
  x3 = (double)xc;
  y0 = ss_cspline_eval (spline, x0, acc);
  y1 = ss_cspline_eval (spline, x1, acc);
  y2 = ss_cspline_eval (spline, x2, acc);
  y3 = ss_cspline_eval (spline, x3, acc);

  A = y0 / ((x0-x1)*(x0-x2)*(x0-x3));
  B = y1 / ((x1-x0)*(x1-x2)*(x1-x3));
  C = y2 / ((x2-x0)*(x2-x1)*(x2-x3));
  D = y3 / ((x3-x0)*(x3-x1)*(x3-x2));
  a = 3*(D+C+B+A);
  b = -2*(C*x3+B*x3+A*x3+D*x2+B*x2+A*x2+D*x1+C*x1+A*x1+D*x0+C*x0+B*x0);
  c = B*x2*x3+A*x2*x3+C*x1*x3+A*x1*x3+C*x0*x3+B*x0*x3+D*x1*x2+A*x1*x2+D*x0*x2+B*x0*x2+D*x0*x1+C*x0*x1;

  rad = -.5 * (b + SGN (b) * sqrt (b*b - 4*a*c));
  x1 = rad / a;
  x2 = c / rad;

  if (0) {
    if (x1 >= 0 && x1 <= (double)max_x)
      printf ("%d %d %f: %21.14f %21.14f %21.14f\n", xa, xc, x1,
              ss_cspline_eval (spline, x1-.000001, acc),
              ss_cspline_eval (spline, x1, acc),
              ss_cspline_eval (spline, x1+.000001, acc));
    else
      printf ("%d %d %f\n", xa, xc, x1);

    if (x2 >= 0 && x2 <= (double)max_x)
      printf ("%d %d %f: %21.14f %21.14f %21.14f\n", xa, xc, x2,
              ss_cspline_eval (spline, x2-.000001, acc),
              ss_cspline_eval (spline, x2, acc),
              ss_cspline_eval (spline, x2+.000001, acc));
    else
      printf ("%d %d %f\n", xa, xc, x1);
    printf ("%6.0f %6.0f %6.0f %6.0f\n",
            ss_cspline_eval (spline, (double)(xa - 1), acc),
            ss_cspline_eval (spline, (double)xa, acc),
            ss_cspline_eval (spline, (double)xc, acc),
            ss_cspline_eval (spline, (double)(xc + 1), acc));
  }

  
  {
    int sign = SGN (yb - mean);
    int x1in, x2in;

    if (x1 >= x0 && x1 <= x3)
      y1 = ss_cspline_eval (spline, x1, acc) - mean;
    if (x2 >= x0 && x2 <= x3)
      y2 = ss_cspline_eval (spline, x2, acc) - mean;
    x1in = x1 >= x0 && x1 <= x3 && SGN (y1) == sign;
    x2in = x2 >= x0 && x2 <= x3 && SGN (y2) == sign;

    if (x1in && !x2in)
      return x1;
    if (x2in && !x1in)
      return x2;
    if (!x1in && !x2in)
      return -1;
    if (x1in && x2in)
      return fabs (y1) > fabs (y2) ? x1 : x2;

  }
  return -1;
}


static double
spline_peak (gsl_spline *spline, gsl_interp_accel *acc, int max_x, int rpeak)
{
  double x0, x1, yb;

  yb = ss_cspline_eval (spline, (double)rpeak, acc);
  x0 = spline_peak_1 (spline, acc, max_x, rpeak - 1, rpeak, yb);
  x1 = spline_peak_1 (spline, acc, max_x, rpeak, rpeak + 1, yb);
  //  printf ("%d %f %f\n", rpeak, x0, x1);

  if (x0 >= 0 && (x1 == -1 || x1 == x0))
    return x0;
  if (x1 >= 0 && x0 == -1)
    return x1;
  if (x0 >= 0 && x1 >= 0) {
    double y0 = fabs (ss_cspline_eval (spline, x0, acc));
    double y1 = fabs (ss_cspline_eval (spline, x1, acc));
    if (y0 > y1)
      return x0;
    if (y1 > y0)
      return x1;
  }
  //printf ("x0: %.17f, x1: %.17f\n", x0, x1);
  {
    float y[SNAPLEN], m[SNAPLEN];
    int n, c;
    for (n = 0; n < SNAPLEN; n++) {
      y[n] = 90 * (ss_cspline_eval (spline, rpeak - 1 + n * 1.0/32, acc) + 10492);
      m[n] = mean;
    }

    window ();
    c = 0;
    while (c != 'c' && c != 'n') {
      clear ();
      draw_snap (m, "black");
      draw_snap (y, "black");
      c = show ();
      if (c == 'q')
        exit (0);
    }
    wclose ();
  }
  exit (DIE);
  return -1;
}

void
region_show_peaks (RegionData *p, double detection_threshold)
{
  float y[SNAPLEN], m[SNAPLEN], hi[SNAPLEN], lo[SNAPLEN], mh[SNAPLEN], ml[SNAPLEN];
  int rn, pn, pdn, c, n;
  short *rdata;
  //  char *cp_msg;

  rn = 0;
  rdata = p->data;
  window ();
  c = 0;
  for (n = 0; n < SNAPLEN; n++) {
    m[n] = mean;
    mh[n] = mean + detection_threshold;
    ml[n] = mean - detection_threshold;
    hi[n] = 25000;
    lo[n] = -25000;
  }
  for (pn = 0; pn < p->peak_count; pn++) {
    if (p->peak[pn].loc > p->region[rn].sample + p->region[rn].count)
      rdata += p->region[rn++].count;
    pdn = p->peak[pn].loc - p->region[rn].sample;
    //    if (fabs (rdata[pdn]) > 24000 && fabs (rdata[pdn]) < 25000 && !monotonic (rdata, pdn, detection_threshold)) {
    //    if (1) {
    {
      int a;
      printf ("peak %d\n", pn);

      for (a = -5; a <= +5; a++)
        printf (" %d", rdata[pdn+a]);
      printf ("\n");

    }
    //    if ((cp_msg = clipped_peak (rdata + pdn))) {
    if (fabs ((double)rdata[pdn]) > 25000) {
      char msg[256];
      for (n = 0; n < SNAPLEN; n++)
        //      y[n] = ss_cspline_eval (p->regionm[rn].spline, pdn - 16 + n / 2.0, p->regionm[rn].acc);
        y[n] = (float)rdata[pdn - 32 + n];
      {
        Peak *pk = p->peak + pn;
        float v[SNAPLEN], w[SNAPLEN];
        int count = pk->lnear + pk->rnear + 1;
        double tmp, tmpr, sum = 0, sumr = 0, mymean = 0;

        for (n = 0; n < count; n++)
          mymean += (v[n] = w[n] = (float)rdata[pdn - pk->lnear + n]);
        mymean /= count;
        whiten_samples (w, count);

        for (n = 0; n < count; n++)
          printf (" %f\n", w[n]);
        printf ("\n");

        for (n = 0; n < count; n++)
          tmp = w[n], sum += tmp*tmp,
            tmpr = v[n] - mymean, sumr += tmpr*tmpr;
        tmp = sqrt (sum);
        tmpr = sqrt (sumr);
        snprintf (msg, 256, "%d: %f %f", pn, map_radius (tmp, count), map_radius (tmpr / unzapped_sd, count));
      }
      clear ();
      draw_snap (m, "black");
      draw_snap (mh, "gray");
      draw_snap (ml, "gray");
      draw_snap (hi, "gray");
      draw_snap (lo, "gray");
      draw_snap (y, "black");
      draw_text (msg);
      c = show ();
      if (c == 'q')
        exit (0);
      if (c == 'c')
        break;
    }
  }
  wclose ();

}

static inline double
distance_samples (float *ar, float *aw, float *br, float *bw, int count)
{
  int n;
  double  sum = 0,  tmp;
#ifdef MIN_W_R
  double sumr = 0, tmpr;
#endif
  
  for (n = 0; n < count; n++) {
    tmp = aw[n] - bw[n], sum += tmp*tmp;
#ifdef MIN_W_R
    tmpr = ar[n] - br[n], sumr += tmpr*tmpr;
#endif
  }

#ifdef MIN_W_R
  return MIN (map_radius (sqrt(sum), count), map_radius (sqrt(sumr) / unzapped_sd, count));
#else
  return map_radius (sqrt(sum), count);
#endif


}

static int
doubly_clipped (short *y_short, int count)
{
  static int alloc;
  static float *y_float;
  AmpData ampdata;
  double xt;
  int n;

  if (count > alloc)
    TREALLOC (y_float, alloc = count);
  for (n = 0; n < count; n++)
    y_float[n] = y_short[n];
  ampdata.peakloc = -SNAPLEN;
  ampdata.ampcount = 0;
  return next_spike (y_float, count, 0, &ampdata, &xt);
}

void
region_classify (RegionData *p)
{
  int rn, pn, count, n, sign[3], i, sum;
  Peak *pk;
  short *data;
  float v[SNAPLEN], w[SNAPLEN], m[SNAPLEN], mw[SNAPLEN];
  int msgcnt = 0;
  int typecount[6];

  memset (typecount, 0, sizeof typecount);
  for (n = 0; n < SNAPLEN; n++)
    m[n] = mean, mw[n] = 0;
  pn = 0;
  data = p->data;
  for (rn = 0; rn < p->region_count; data += p->region[rn++].count) {

    count = 0;
    pk = p->peak + pn;
    while (pn < p->peak_count && p->peak[pn].loc < p->region[rn].sample + p->region[rn].count)
      pn++, count++;

    if (p->region[rn].type == -3)
      continue;
    p->region[rn].type = -1;

    if (0)
      if (rn < 10) printf ("%s line %d: region %d, %f to %f,  %d peaks\n",
                           __FILE__, __LINE__, rn, p->region[rn].sample / 25000.0,
                           (p->region[rn].sample + p->region[rn].count) / 25000.0, count);
    if (count == 0 || count > 3 || (count > 1 && doubly_clipped (data, p->region[rn].count)))
      continue;
    for (n = 0; n < count; n++)
      (sign[n] = SGN (data[pk[n].loc - p->region[rn].sample] - mean)) || DIE;
    for (n = 1; n < count; n++)
      if (sign[n] == sign[n-1] || pk[n-1].loc + pk[n-1].rfar < pk[n].loc - pk[n].lfar)
        break;
    //    if (rn < 10) printf ("%s line %d: region %d, %d peaks alternate sign\n", __FILE__, __LINE__, rn, n);
    if (n < count) {
      if (0)
        if (msgcnt++ < 10 && count == 3)
          printf ("region %d: only %d of %d peaks\n", rn, n, count);
      continue;
    }
    for (i = 0; i < SNAPLEN; i++)
      v[i] = (float)data[43 + i];
    memcpy (w, v, sizeof (SnapVals));
    whiten_snap (w);
    if (distance_samples (v, w, m, mw, SNAPLEN) < 11.5) {
      p->region[rn].type = -2;
      if (0) {
        if (rn < 10) printf ("%s line %d: region %d in the noise\n", __FILE__, __LINE__, rn);
        if (msgcnt++ < 10 && count == 3)
          printf ("region %d: distance is %f\n", rn, distance_samples (v, w, m, mw, SNAPLEN));
      }
      continue;
    }
    p->region[rn].type = (count - 1) * 2 + (sign[0] == 1);
    //    if (rn < 10) printf ("%s line %d: region %d type %d\n", __FILE__, __LINE__, rn, p->region[rn].type);
    (p->region[rn].type >= 0 && p->region[rn].type < 6) || DIE;
    typecount[(int)p->region[rn].type]++;
  }
  for (sum = n = 0; n < 6; n++) {
    sum += typecount[n];
    if (typecount[n])
      printf ("region type %d: %7d regions\n", n, typecount[n]);
  }
  printf ("%d typed regions\n", sum);
}


typedef struct {char count; struct {char s, c;} m[3];} RegionTableEntry;

/*@-initallelements -fullinitblock@*/
static RegionTableEntry region_table[6][6] =
  {{{1,{{0,0}}}, {0},         {1,{{0,0}}},       {1,{{0,1}}},       {0},                     {1,{{0,1}}}},
   {{0},         {1,{{0,0}}}, {1,{{0,1}}},       {1,{{0,0}}},       {1,{{0,1}}},             {0}},
   {{1,{{0,0}}}, {1,{{1,0}}}, {2,{{0,0},{1,1}}}, {0},               {2,{{0,0},{1,1}}},       {2,{{0,1},{1,2}}}},
   {{1,{{1,0}}}, {1,{{0,0}}}, {0},               {2,{{0,0},{1,1}}}, {2,{{0,1},{1,2}}},       {2,{{0,0},{1,1}}}},
   {{0},         {1,{{1,0}}}, {2,{{0,0},{1,1}}}, {2,{{1,0},{2,1}}}, {3,{{0,0},{1,1},{2,2}}}, {0}},
   {{1,{{1,0}}}, {0},         {2,{{1,0},{2,1}}}, {2,{{0,0},{1,1}}}, {0},                     {3,{{0,0},{1,1},{2,2}}}}
  };
/*@=initallelements =fullinitblock@*/

//static RegionTableEntry region_table[6][6] =
//  {{{1,{{0,0}}}, {0},         {0,{{0,0}}},       {0,{{0,1}}},       {0},                     {0,{{0,1}}}},
//   {{0},         {1,{{0,0}}}, {0,{{0,1}}},       {0,{{0,0}}},       {0,{{0,1}}},             {0}},
//   {{0,{{0,0}}}, {0,{{1,0}}}, {2,{{0,0},{1,1}}}, {0},               {0,{{0,0},{1,1}}},       {0,{{0,1},{1,2}}}},
//   {{0,{{1,0}}}, {0,{{0,0}}}, {0},               {2,{{0,0},{1,1}}}, {0,{{0,1},{1,2}}},       {0,{{0,0},{1,1}}}},
//   {{0},         {0,{{1,0}}}, {0,{{0,0},{1,1}}}, {0,{{1,0},{2,1}}}, {3,{{0,0},{1,1},{2,2}}}, {0}},
//   {{0,{{1,0}}}, {0},         {0,{{1,0},{2,1}}}, {0,{{0,0},{1,1}}}, {0},                     {3,{{0,0},{1,1},{2,2}}}}
//  };

typedef struct {
  int start, end;
} RegionRange;

static inline void
region_intersect (RegionRange *r1, RegionRange *r2)
{
  r1->start = MAX (r1->start, r2->start);
  r1->end = MIN (r1->end, r2->end);
}

static inline RegionRange *
region_range (Peak *s, Peak *c, Peak *first, Peak *last)
{
  static RegionRange r;
  int s0, s1, c0, c1;

  s0 = s->loc - s->lnear;
  s1 = s->loc + s->rnear;
  c0 = c->loc - c->lnear;
  c1 = c->loc + c->rnear;
  r.start = c0 + 1 - s1;
  r.end =   c1 - 1 - s0;

  if (last->loc + last->rnear >= r.start + SNAPLEN) /* keep all the peaks of the candidate in the window */
    r.start = last->loc + last->rnear - (SNAPLEN+1);
  if (first->loc - first->lnear < r.end)
    r.end = first->loc - first->lnear;

  return &r;
}

static double
region_distance_at_start (double start, void *p, int debug)
{
  Seed *s = (Seed *)p;
  RegionM *region = s->regionm;
  float v[SNAPLEN], w[SNAPLEN];
  int n;
  double  sum = 0,  tmp;
#ifdef MIN_W_R
  double sumr = 0, tmpr;
#endif

  for (n = 0; n < SNAPLEN; n++)
    w[n] = v[n] = ss_cspline_eval (region->spline, start + n, region->acc);
  whiten_snap (w);


#ifdef MIN_W_R
  for (n = 0; n < SNAPLEN; n++)
    tmp = s->val[n] - w[n], sum += tmp*tmp,
      tmpr = s->raw[n] - v[n], sumr += tmpr*tmpr;
  return MIN (sqrt(sum), sqrt (sumr) / unzapped_sd);
#else
  for (n = 0; n < SNAPLEN; n++)
    tmp = s->val[n] - w[n], sum += tmp*tmp;
  return sqrt(sum);
#endif
}

static SnapList *
region_aligned_snaps (RegionData *p, Seed *seed, int seed_region, int pass, int spike_count)
{
  int rn, m;
  RegionTableEntry *rte;
  RegionRange range;
  Peak *peak0, *s_peak, *peakn;
  SnapList *sl;
  RegionM *rm;

  TCALLOC (sl, 1);
  s_peak = seed->peak;
  for (rn = 0; rn < p->region_count; rn++) {
    p->regionm[rn].mark = 0;
    if (p->region[rn].type < 0 || p->regionm[rn].cluster >= 0)
      continue;

    if ((now = time (0)) > last_time) {
      if (spike_count != 0)
        printf ("%79s\r  region_cluster: %5d of %d pass %d, count %5d, region_aligned_snaps: %5d of %d\r",
                "", seed_region, p->region_count, pass, spike_count, rn, p->region_count);
      fflush (stdout);
      last_time = now;
    }

    if ((rte = &region_table[(int)seed->type][(int)p->region[rn].type])->count == 0)
      continue;
    range.start = -INT_MAX;
    range.end = INT_MAX;
    rm = p->regionm + rn;

    peak0 = rm->peak;
    peakn = rm->peak + rm->pkcnt - 1;

    for (m = 0; m < rte->count; m++) {
      peak0 + rte->m[m].c <= peakn || DIE;
      region_intersect (&range, region_range (s_peak + rte->m[m].s, peak0 + rte->m[m].c, peak0, peakn));
    }
    //printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);

    if (range.start > range.end) {
      int r0 = p->region[rn].sample;
      Peak *sp, *cp;
      //      printf ("region %5d seed type %d region type %d: no intersection\n", rn, seed->type, p->region[rn].type);
      range.start = -INT_MAX;
      range.end = INT_MAX;
      for (m = 0; m < rte->count; m++) {
        //      printf ("m = %d: %d %d\n", m, range.start - r0, range.end - r0);
        sp = s_peak + rte->m[m].s;
        cp = peak0 + rte->m[m].c;
        if (0) {
          printf ("%d: %d %d %d, %d: %d %d %d\n",
                  rte->m[m].s, sp->lnear, sp->loc, sp->rnear, rte->m[m].c, cp->lnear, cp->loc - r0, cp->rnear);
          printf ("m = %d: %d %d\n", m, range.start - r0, range.end - r0);
        }
      }
      continue;
    }
    if (range.start < p->region[rn].sample)
      range.start = p->region[rn].sample;
    if (!(range.start >= p->region[rn].sample && range.end <= p->region[rn].sample + 63))
      continue;
    range.start -= p->region[rn].sample;
    range.end -= p->region[rn].sample;
    {
      double start, newstart;
      int status;

      if (rte->m[0].s == 0)
        start = (double)((peak0[(int)rte->m[0].c].loc - p->region[rn].sample) - PRESAMPLES);
      else
        start = (range.start + range.end) / 2.0;

      seed->regionm = &p->regionm[rn];
      //printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
      newstart = fmin2 (&region_distance_at_start, seed, start, 1, (double)range.start, (double)range.end, &status);
      //printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
      if (status == 0) {
        int n, spike;
        Snap *s;

        if (fmin_y () < 11.5)
          p->regionm[rn].mark = 1;
        for (spike = 0; spike < rm->spike_count; spike++) {
          //printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
          TCALLOC (s, 1);
          TMALLOC (s->raw, SNAPLEN); s->free_raw = 1;
          TMALLOC (s->val, SNAPLEN); s->free_val = 1;
          for (n = 0; n < SNAPLEN; n++)
            s->raw[n] = ss_cspline_eval (p->regionm[rn].spline, newstart + n, p->regionm[rn].acc);
          memcpy (s->val, s->raw, sizeof (SnapVals));
          whiten_snap (s->val);
          s->sample = p->region[rn].sample + newstart + PRESAMPLES + rm->offset[spike];
          if ((s->distance = distance (s->val, seed->val)) < 11.5)
            s->close = p->regionm[rn].mark = 1;
          else
            s->close = 0;
          add_snap (sl, s);
        }
      }      
    }
  }
  return sl;
}

static inline Peak *
do_avg_peak (int b, int dir, float *v)
{
  int a, c, i, stop, lnear, rnear;
  static int r0 = MPK_R0;
  static int r = MPK_R;
  int b0 = b;
  int peakloc = starting_sample + b;
  int debug = 0;
  int val = (int)(v[b] - mean);
  static Peak pk;

  if (dir * (v[b] - mean) < 0) {
    if (debug) printf ("line %d: peak at %d %6d %2d doesn't cross mean\n", __LINE__, peakloc, val, dir);
    return 0;
  }

  stop = MAX (b - (2 * r0 - 1), 0);
  for (i = b - 1; i >= stop && v[b] == v[i]; i--);
  a = MAX (b - 2 * r0, 0);
  for ( ; i >= stop && dir * (v[b] - v[i]) > 0; i--)
    if ((double)(dir * (v[b] - v[i])) >= detection_threshold) {
      a = i;
      break;
    }
  if (i < 0) {
    if (debug) printf ("line %d: peak at %d %6d %2d overlaps beginning of buffer\n", __LINE__, peakloc, val, dir);
    return 0;
  }

  stop = MIN (b + (2 * r0 - 1), SNAPLEN - 1);
  for (i = b + 1; i <= stop && v[b] == v[i]; i++);
  c = MIN (b + 2 * r0, SNAPLEN - 1);
  for ( ; i <= stop && dir * (v[b] - v[i]) > 0; i++)
    if ((double)(dir * (v[b] - v[i])) >= detection_threshold) {
      c = i;
      break;
    }
  if (i >= SNAPLEN) {
    if (debug) printf ("line %d: peak at %d %6d %2d overlaps end of buffer\n", __LINE__, peakloc, val, dir);
    return 0;
  }
  if ((double)(dir * (v[b] - v[c])) < detection_threshold || (double)(dir * (v[b] - v[a])) < detection_threshold || c - a > 2 * r0) {
    if (debug) printf ("line %d: peak at %d %6d %2d below detection threshold\n", __LINE__, peakloc, val, dir);
    return 0;
  }

  lnear = b0 - a;
  rnear = c - b0;

  if (b - a > r0 || c - b > r0)
    b = (a + c) / 2;

  stop = MAX (b - r, 0);
  for (i = a - 1; i >= stop && dir * (v[b] - v[i]) > 0; i--)
    if (dir * (v[a] - v[i]) > 0)
      a = i;
  if (i < 0) {
    if (debug) printf ("line %d: peak at %d %6d %2d overlaps beginning of buffer\n", __LINE__, peakloc, val, dir);
    //    return 0;
  }

  stop = MIN (b + r, SNAPLEN - 1);
  for (i = c + 1; i <= stop && dir * (v[b] - v[i]) > 0; i++)
    if (dir * (v[c] - v[i]) > 0)
      c = i;
  if (i >= SNAPLEN) {
    if (debug) printf ("line %d: peak at %d %6d %2d overlaps end of buffer\n", __LINE__, peakloc, val, dir);
    //    return 0;
  }
  //  printf ("peak detected: %d %d %d\n", lnear, b0, rnear);
  pk.loc = b0;
  pk.lfar = b0 - a;
  pk.lnear = lnear;
  pk.rnear = rnear;
  pk.rfar = c - b0;
  return &pk;
}

static int
find_avg_peaks (float *v, Peak **pkp, int alloc)
{
  int n, i, last, dir, pkcnt, pn;
  Peak *pk;
  static Peak *s_peak;
  Peak *peak;

  peak = alloc ? 0 : s_peak;
  last = SGN (v[1] - v[0]);
  pkcnt = i = 0;
  for (n = 2; n < SNAPLEN; n++) {
    if ((dir = SGN (v[n] - v[n-1])) && dir != last) {
      pk = do_avg_peak (n - 1, -dir, v);
      if (pk) {
        pn = pkcnt++;
        TREALLOC (peak, pkcnt);
        peak[pn] = *pk;
      }
      last = dir;
    }
  }
  *pkp = peak;
  if (!alloc) s_peak = peak;
  return MIN (pkcnt, 3);
}

static double
extrap_val (float *v, double x)
{
  double x0 = 0, v0 = 0, wn, dx = 0;
  int n;
  wn = 2 * M_PI / 32;

  if (x <= 0) {
    x0 = v[0] - mean;
    v0 = v[0] - v[1];
    dx = -x;
  }
  else if (x >= SNAPLEN-1) {
    x0 = v[SNAPLEN-1] - mean;
    v0 = v[SNAPLEN-1] - v[SNAPLEN-2];
    dx = x - (SNAPLEN-1);
    return (x0 + (v0 + wn*x0) * dx) * exp (-wn * dx) + mean;
  }
  else if ((double)(n = (int)x) == x)
    return v[n];
  else exit (DIE);
  return (x0 + (v0 + wn*x0) * dx) * exp (-wn * dx) + mean;
}

static void
avg_spline_peak (float *v, Peak *peak)
{
  double x[SNAPLEN*2], y[SNAPLEN*2];
  double xp;
  int n;
  static gsl_interp_accel *acc;
  static gsl_spline *spline;

  if (!acc) {
    (acc || spline) && DIE;
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_cspline, SNAPLEN*2);
  }

  for (n = 0; n < SNAPLEN*2; n++)
    x[n] = (double)(n-PRESAMPLES), y[n] = extrap_val (v, (double)(n - PRESAMPLES));

  gsl_spline_init (spline, x, y, SNAPLEN*2);
  xp = spline_peak (spline, acc, SNAPLEN-1, peak->loc);
  for (n = 0; n < SNAPLEN; n++)
    v[n] = ss_cspline_eval (spline, xp - PRESAMPLES + n, acc);
}

#ifdef JUSTONE
static float dpeak[SNAPLEN], cpeak[SNAPLEN], pcount;
static double gscale;
#endif

static double
scaled_distance (float *to_be_scaled, float *fixed, int count, double *scalep)
{
  int n;
  double numer = 0, denom = 0, scale, tmp, sum;

  for (n = 0; n < count; n++) {
    numer += (double) to_be_scaled[n] * fixed[n];
    denom += (double) to_be_scaled[n] * to_be_scaled[n];
  }
  scale = numer / denom;

#ifdef JUSTONE
  gscale = scale;
#endif
  *scalep = scale;
  sum = 0;
  for (n = 0; n < count; n++) {
    tmp = scale * to_be_scaled[n] - fixed[n];
    sum += tmp * tmp;
  }

  return sqrt (sum);
}

static double
region_peak_distance_at_start (double start, void *p, int debug)
{
  Seed *s = (Seed *)p;
  RegionM *region = s->regionm;
  float v[SNAPLEN], w[SNAPLEN];
  int n, count;
  double scale, d;
  static double min_d;

  count = s->window_size;
  if (debug) {
    for (n = -1; n < count + 1; n++)
      printf ("%f\n", ss_cspline_eval (region->spline, start + n, region->acc));
    exit (DIE);
  }
  for (n = 0; n < count; n++)
    w[n] = v[n] = ss_cspline_eval (region->spline, start + n, region->acc);

  whiten_samples (w, count);
#ifdef JUSTONE
  TMEMCPY (dpeak, w, count);
  TMEMCPY (cpeak, s->val, count);
  pcount = count;
#endif
  if (do_scaled_distance) {
    d = scaled_distance (s->val, w, count, &scale);
    if (s->scale == 0 || d < min_d) {
      s->scale = scale;
      min_d = d;
    }
    s->scale || DIE;
    return d;
  }
  else
    return distance_samples (s->raw, s->val, v, w, count);
}

static void
snaplist_spikedata (SnapList *sl, int cluster)
{
  int n, i;
  Snap *s;

  for (n = 0; n < sl->count; n++)
    if ((s = sl->snap[n])->close) {
      i = spike_count++;
      memset (&spikedata[i], 0, sizeof spikedata[i]);
//      spikedata[i].leave_out = 0;
//      spikedata[i].nocenter = 0;
//      spikedata[i].scale = 0;
//      spikedata[i].overlapped = 0;
      spikedata[i].partition = cluster;
      spikedata[i].sample = s->sample;
      //    printf ("%s line %d: %f\n", __FILE__, __LINE__, spikedata[i].sample);
      spikedata[i].distance = s->distance;
      spike_waveform[i] = s->raw;
      s->free_raw = 0;
    }
}

static void
center_spline_init_1 (Center *ctr)
{
  double tol = .000001;
  int cnt = 10;
  static double *x, *y;
  static int alloc;
  double val;
  int len, left, last, start, n, close;

  maxregion || DIE;
  if (alloc < maxregion * 2 + SNAPLEN) {
    alloc = maxregion * 2 + SNAPLEN;
    TREALLOC (x, alloc);
    TREALLOC (y, alloc);
  }
    
  for (n = 0; n < SNAPLEN; n++)
    x[maxregion + n] = (double)n, y[maxregion + n] = ctr->raw[n];
  
  n = 0;
  val = fabs (ctr->raw[0] - mean);
  close = 0;
  while (close < cnt && n < maxregion) {
    n++;
    x[maxregion - n] = (double)-n;
    y[maxregion - n] = extrap_val (ctr->raw, (double)-n);
    val = fabs (y[maxregion - n] - mean);
    if (val < tol)
      close++;
    else
      close = 0;
  }
  if (n < maxregion) {
    n++;
    x[maxregion - n] = (double)-maxregion;
    y[maxregion - n] = mean;
  }
  left = n;
  
  n = 0;
  last = maxregion + SNAPLEN-1;
  val = fabs (ctr->raw[SNAPLEN-1] - mean);
  close = 0;
  while (close < cnt && n < maxregion) {
    n++;
    x[last + n] = (double)(SNAPLEN-1 + n);
    y[last + n] = extrap_val (ctr->raw, (double)(SNAPLEN-1 + n));
    val = fabs (y[last + n] - mean);
    if (val < tol)
      close++;
    else
      close = 0;
  }
  if (n < maxregion) {
    n++;
    x[last + n] = (double)(last + maxregion);
    y[last + n] = mean;
  }
  len = left + SNAPLEN + n;
  start = maxregion - left;

  ctr->acc = gsl_interp_accel_alloc ();
  len >= 0 || DIE;
  ctr->spline = gsl_spline_alloc (gsl_interp_cspline, (size_t)len);
  gsl_spline_init (ctr->spline, x+start, y+start, (size_t)len);
}

void
region_ctr_init (Center *ctr)
{
  int pcnt, n;
  Peak *peak;

  pcnt = find_avg_peaks (ctr->raw, &peak, 0);
  pcnt = pcnt > 3 ? 3 : pcnt;
  for (n = 0; n < pcnt; n++)
    ctr->peak[n] = peak[n];
  ctr->peak_count = pcnt;
  ctr->type = (pcnt - 1) * 2 + (SGN (ctr->raw[ctr->peak[0].loc] - mean) == 1);
  center_spline_init_1 (ctr);
}

SnapList *
region_cluster (DecomposeData *dd)
{
  RegionData *p = dd->p;
  RegionM *rm;
  int rn, n, i, start, pcnt, last_pcnt, ch = 0, cluster = 0;
  int pass, spike_count = 0;
  Seed seed;
  Peak peak[3];
  float v[SNAPLEN], w[SNAPLEN];
  SnapList *sl = 0, *ctr, *cc;
  Snap s;
  Peak *apk;
  float last[SNAPLEN];

  seed.peak = peak;
  seed.val = w;
  seed.raw = v;
  s.raw = v;
  s.val = w;
  TCALLOC (cc, 1);
  for (rn = 0; rn < p->region_count; rn++)
    p->regionm[rn].cluster = -1;
  for (rn = 0; rn < p->region_count; rn++) {
    if (p->region[rn].type < 0 || p->regionm[rn].cluster >= 0)
      continue;

    if ((now = time (0)) > last_time) {
      printf ("%79s\r  region_cluster: %5d of %d\r", "", rn, p->region_count);
      fflush (stdout);
      last_time = now;
    }

    rm = p->regionm + rn;
    seed.type = p->region[rn].type;

    start = BEFORE-PRESAMPLES;
    for (i = 0; i < SNAPLEN; i++)
      last[i] = w[i] = v[i] = ss_cspline_eval (p->regionm[rn].spline, (double)(start + i), p->regionm[rn].acc);
    whiten_snap (w);
    rm->pkcnt <= 3 || DIE;
    for (n = 0; n < rm->pkcnt; n++) {
      /*@-modobserver@*/
      peak[n] = rm->peak[n];
      /*@=modobserver@*/
      peak[n].loc -= p->region[rn].sample + start;
    }
    last_pcnt = rm->pkcnt;
    //    printf ("seed type %d\n", seed.type);
    if (0)
    {
      int pkcnt = seed.type / 2 + 1;
      int n;
      printf ("seed type %d: peaks at", seed.type);
      for (n = 0; n < pkcnt; n++)
        printf (" %d", peak[n].loc);
      printf ("\n");
    }
    pass = 0;
    ch = 'c';
    do {
      free_snaplist (&sl);
      if (cluster == 6)
        spike_count = 0;
      //      printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
      sl = region_aligned_snaps (p, &seed, rn, pass, spike_count);
      //      printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
      {
        spike_count = 0;
        for (n = 0; n < p->region_count; n++)
          spike_count += p->regionm[n].mark;
        //      printf ("pass %d: %2d spikes in cluster\n", pass, spike_count);
      }
      //      printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
      ctr = center_points_full (sl, &s);
      //      printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
      //      printf ("ctr->count: %d\n", ctr->count);
      pcnt = 0;
      if (ctr->count) {
        //      printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
        pcnt = find_avg_peaks (ctr->snap[0]->raw, &apk, 0);
        if (0 && cluster > 10) {
          char msg[256];
          float m[SNAPLEN];
          int n;
          snprintf (msg, 256, "cluster %d: %d peaks at %d %d %d before", cluster, pcnt, apk[0].loc, apk[1].loc, apk[2].loc);
          for (n = 0; n < SNAPLEN; n++)
            m[n] = mean;
          window ();
          clear ();
          draw_snap (m, "gray");
          draw_snap (ctr->snap[0]->raw, "red");
          draw_text (msg);
          show ();
        }
        if (pcnt)
          avg_spline_peak (ctr->snap[0]->raw, apk);
        pcnt = find_avg_peaks (ctr->snap[0]->raw, &apk, 0);
        if (0 && cluster > 10) {
          char msg[256];
          float m[SNAPLEN];
          int n;
          snprintf (msg, 256, "cluster %d: %d peaks at %d %d %d after", cluster, pcnt, apk[0].loc, apk[1].loc, apk[2].loc);
          for (n = 0; n < SNAPLEN; n++)
            m[n] = mean;
          window ();
          clear ();
          draw_snap (m, "gray");
          draw_snap (ctr->snap[0]->raw, "red");
          draw_text (msg);
          show ();
          wclose ();
        }
        
        if (pcnt && apk[0].loc != PRESAMPLES) {
          printf ("%s line %d: first peak not at %d\n",  __FILE__, __LINE__, PRESAMPLES);
          fflush (stdout);
          pcnt = 0;
        }
        //      printf ("peak count: %d\n", pcnt);
      }
      //      printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
      fflush (stdout);
      if (ch != 'c')
      {
        float m[SNAPLEN];
        int n;
        for (n = 0; n < SNAPLEN; n++)
          m[n] = mean;

        window ();
        clear ();
        draw_snap (m, "gray");
        draw_snap (v, "green");
        if (ctr->count)
          draw_snap (ctr->snap[0]->raw, "red");
        ch = show ();
        if (ch == 'q')
          exit (0);
        if (ch == 'c')
          wclose ();
      }
      if (ctr->count == 0)
        break;
      if (pcnt == 0)
        break;
      memcpy (v, ctr->snap[0]->raw, sizeof (SnapVals));
      memcpy (w, v, sizeof (SnapVals));
      whiten_snap (w);
      {
        int i;
        double max = 0, d;

        for (i = 0; i < SNAPLEN; i++)
          if ((d = fabs (v[i] - last[i])) > max)
            max = d;

        //      printf ("%s line %d: max change %f\n", __FILE__, __LINE__, max); fflush (stdout);
        if (pass++ == 40 || max < .5)
          ch = 'n';
        else {
          memcpy (last, v, sizeof (SnapVals));
          pcnt <= 3 || DIE;
          for (n = 0; n < pcnt; n++)
            peak[n] = apk[n];
          last_pcnt = pcnt;
        }
      }
    } while (ch != 'n');
    //    printf ("%d passes\n", pass);
    //printf ("%79s\r %s line %d\n", "", __FILE__, __LINE__); fflush (stdout);
    if (ctr->count && last_pcnt > 0 && pcnt > 0) {
      int classified = 0, notype = 0, unclassified = 0, spkcnt = 0, ccn;
      Snap *s;

      snaplist_spikedata (sl, cluster);
      ccn = dd->center_count++;
      TREALLOC (dd->center, dd->center_count);
      TMEMSET (&dd->center[ccn], 0, 1);
      SZMEMCPY (dd->center[ccn].peak, peak);
      dd->center[ccn].peak_count = last_pcnt;
      dd->center[ccn].type = (last_pcnt - 1) * 2 + (SGN (last[PRESAMPLES] - mean) == 1);
      SZMEMCPY (dd->center[ccn].raw, last);
      center_spline_init_1 (&dd->center[ccn]);
      TCALLOC (s, 1);
      TMALLOC (s->raw, SNAPLEN);
      TMALLOC (s->val, SNAPLEN);
      memcpy (s->raw, last, sizeof (SnapVals));
      memcpy (s->val, last, sizeof (SnapVals));
      whiten_snap (s->val);
      add_snap (cc, s);
      for (n = 0; n < p->region_count; n++)
        if (p->regionm[n].mark)
          p->regionm[n].cluster = cluster, classified++, spkcnt++, s->count++;
        else if (p->regionm[n].cluster >= 0)
          classified++;
        else if (p->region[n].type < 0)
          notype++;
        else
          unclassified++;
      cluster++;
      printf ("%79s\rcluster %2d: %4d spikes, %d peaks\n", "", cluster - 1, spkcnt, last_pcnt);
      //      printf ("%d classified, %d unclassified, %d no type\n", classified, unclassified, notype);
    }
  }
  //  centers_vs_regions (p, cc);
  //  show_centers (cc);
  return cc;
}

typedef struct {
  gsl_spline *spline;
  gsl_interp_accel *acc;
  double length;
  double *start;
  double *offset;
  int count;
} Composite;

void
lpfilter_snap (float *v, float *out)
{
  int n;
  float vx_array[SNAPLEN*2];
  float *vx = vx_array + SNAPLEN/2;

  for (n = -SNAPLEN / 2; n < SNAPLEN + SNAPLEN/2; n++)
    vx[n] = extrap_val (v, (double)n);
  for (n = 0; n < SNAPLEN; n++)
    out[n] = lp_filter (vx + n);
}

void
region_optavg (RegionData *p)
{
  int rn;
  RegionM *rm;

  for (rn = 0; rn < p->region_count; rn++) {
    rm = p->regionm + rn;
    if (p->region[rn].type == -1)
      if (p->region[rn].count > maxregion)
        maxregion = p->region[rn].count;

    /* optavg deleted for now */

    {
      static double offset;
      rm->spike_count = 1;
      rm->offset = &offset;
      rm->peak = p->peak + rm->pn;
      rm->pkcnt = p->region[rn].type / 2 + 1;
    }
  }
  if (0) printf ("%79s\rlargest region has %d samples\n", "", maxregion);
}

static RegionRange *
peak_range (Peak *c, Peak *s, int origin)
{
  static RegionRange r;
  int s0, s1, c0, c1, cloc;
  /* s is seed (center), and c is candidate (region peak) */

  cloc = c->loc - origin;
  s0 = s->loc - s->lnear;
  s1 = s->loc + s->rnear;
  c0 = cloc - c->lnear;
  c1 = cloc + c->rnear;
  r.start = c0 + 1 - s1;
  r.end   = c1 - 1 - s0;
  return &r;
}

static int
center_window_start (Seed *seed, int pnc0, int count)
{
  int diff, d, start0, window_size, s0, n;
  Peak *p0, *pn;
  
  p0 = seed->peak + pnc0;
  pn = p0 + (count - 1);

  if ((diff = 2 * WRADIUS0 - (d = p0->lnear + pn->rnear + (pn->loc - p0->loc))) > 0) {
    diff += diff & 1;
    s0 = diff / 2;
  }
  else {
    s0 = 0;
  }

  window_size = d + 1;
  start0 = p0->loc - p0->lnear - s0;
  if (start0 + window_size > SNAPLEN)
    window_size = SNAPLEN - start0;
  seed->raw = seed->cdata + start0;
  for (n = 0; n < window_size; n++)
    seed->val[n] = seed->raw[n];
  whiten_samples (seed->val, window_size);
  seed->window_size = window_size;
  return start0;
}

static void
range_trim (RegionRange *range, int start0, int window_size, int count)
{
  range->start += start0;
  range->end += start0;
  if (range->end + window_size > count) {
    range->end = count - window_size;
    if (range->start + window_size > count)
      range->start = count - window_size;
  }
}

#define SHOWPUSH 0
#define SHOWSUB 0

static int
region_center_peaks (RegionData *p, int pn0, Seed *seed, int pnc0, double *startp, int loose)
{
  int pcnt, region_end, rn, pn, pnc, start0, r0, status;
  double start, goodstart = 0, goodstart0 = 0, flstart;
  RegionRange ref0_range, refn_range, *r = &refn_range;
  double threshold;

  {
    double maxamp = 0, amp, tmp;
    int n;
    for (n = 0; n < SNAPLEN; n++)
      if ((amp = fabs (seed->cdata[n] - mean)) > maxamp)
        maxamp = amp;
    unzapped_sd || DIE;
    tmp = maxamp / unzapped_sd / 24;
    {
      double r = loose ? 20.2 : 11.8;
      threshold = sqrt (r * r + tmp * tmp * SNAPLEN);
    }
#ifdef JUSTONE
    printf ("threshold: %f, maxamp %f\n", threshold, maxamp);
#endif
  }

  rn = p->peak_region[pn0];
  seed->regionm = p->regionm + rn;

  r0 = p->region[rn].sample;
  region_end = r0 + p->region[rn].count;
  pcnt = seed->type / 2 + 1;
  ref0_range.start = 0;
  ref0_range.end = INT_MAX;
  for (pnc = pnc0, pn = pn0; pnc < pcnt; pnc++, pn++) {
    if (pn >= p->peak_count || p->peak[pn].loc >= region_end) {
#ifdef JUSTONE
      printf ("%s line %d: %d >= %d || %d >= %d\n", __FILE__, __LINE__, pn, p->peak_count, p->peak[pn].loc, region_end);
#endif
      break;
    }

    if (0)
    {
      Peak *r, *c;
      r = &p->peak[pn];
      c = &seed->peak[pnc];
      printf ("%d\n", p->region_count);
      printf ("pn: %d, pnc: %d\n", pn, pnc);
      printf ("%d %d, %d %d %d %d; %d,  %d %d %d %d\n",
              r0, r->loc - r0, r->lfar, r->lnear, r->rnear, r->rfar,
              c->loc, c->lfar, c->lnear, c->rnear, c->rfar);
      exit (0);
    }

    region_intersect (&ref0_range, peak_range (p->peak + pn, seed->peak + pnc, r0));
    start0 = center_window_start (seed, pnc0, pnc - pnc0 + 1);
    refn_range = ref0_range;
    range_trim (&refn_range, start0, seed->window_size, p->region[rn].count);
    if (refn_range.start > refn_range.end) {
#ifdef JUSTONE
      printf ("%s line %d: %d > %d\n",  __FILE__, __LINE__, refn_range.start, refn_range.end);
#endif
      break;
    }

    seed->scale = 0;
    start = fmin2 (&region_peak_distance_at_start, seed, (r->start + r->end) / 2.0, 1, (double)r->start, (double)r->end, &status);
    do_scaled_distance == 0 || seed->scale || DIE;
#ifdef JUSTONE
      printf ("%s line %d: start %f, fmin_y %f, count %d, val %f\n",  __FILE__, __LINE__, start, fmin_y (), seed->window_size, seed->raw[0]);
      if (start > 120.490601 && start < 120.490603) {
        int n;
        printf ("%s line %d\n",  __FILE__, __LINE__);
        mfvec ("dpeak", dpeak, pcount);
        if (0)
          for (n = 0; n < pcount; n++)
            cpeak[n] = gscale * (cpeak[n] - mean) + mean;
        mfvec ("cpeak", cpeak, pcount);
      }
#endif

    if (SHOWSUB) 
      printf ("%d %d: %d %f\n", pn, pnc, status, fmin_y ()); 

    if (p->region[0].sample == SAMPLOC)
      printf ("status %d, fmin_y %f vs %f\n", status, fmin_y (), threshold);
    if (status != 0 || fmin_y () > threshold)
      break;
    if (SHOWPUSH)
      printf ("center peak %d matches peak %d (%d samples at %f)\n", pnc, pn, seed->window_size, start);
    goodstart = start;
    goodstart0 = (double)start0;
  }
  if (pnc == pnc0)
    return 0;
  flstart = goodstart - goodstart0;
  if ((int)ceil (flstart) + SNAPLEN > p->region[rn].count)
    return 0;
  if (0)
    for (; pn < p->peak_count && (double)(p->peak[pn].loc - r0) + p->peak[pn].rnear < flstart + SNAPLEN; pn++) {
      seed->window_size = (int)((p->peak[pn].loc - r0) + p->peak[pn].rnear - goodstart);
      if (region_peak_distance_at_start (goodstart, seed, 0) > 11.5)
        break;
      if (SHOWPUSH)
        printf ("peak %d at %d rnear %d matches center starting at %f with window size %d\n",
                pn, p->peak[pn].loc - r0, p->peak[pn].rnear, goodstart, seed->window_size);
    }
  *startp = flstart;
#ifdef JUSTONE
  printf ("center %d: %d peaks matched\n", seed->ccn, pn - pn0);
#endif
  return pn - pn0;
}

static inline int
center_peak_polarity (DecomposeData *dd, int ccn, int pnc)
{
  int peak_sample;

  peak_sample = dd->center[ccn].peak[pnc].loc;
  return SGN (dd->center[ccn].raw[peak_sample] - mean);
}

static inline int
region_peak_polarity (DecomposeData *dd)
{
  int peak_sample;

  peak_sample = dd->p->peak[dd->pn].loc - dd->p->region[dd->rn].sample;
  return SGN (dd->p->fdata[peak_sample] - mean);
}

static inline double
region_peak_amplitude (DecomposeData *dd)
{
  int peak_sample;

  peak_sample = dd->p->peak[dd->pn].loc - dd->p->region[dd->rn].sample;
  return fabs (dd->p->fdata[peak_sample] - mean);
}

static void
free_region_data (/*@only@*/ RegionData *p)
{
  int n;
  CSpline *csp;

  if (p->data_alloc)
    free (p->data);
  else
    p->data == 0 || DIE;
  if (p->peak_alloc)
    free (p->peak);
  else
    p->peak == 0 || DIE;
  if (p->region_alloc) {
    free (p->region);
    if (p->regionm) {
      for (n = 0; n < p->region_count; n++) {
        gsl_interp_accel_free (p->regionm[n].acc);
        csp = (CSpline *)p->regionm[n].spline;
        if (csp->interp == 0 && csp->y_float == p->fdata)
          p->fdata = 0;
        ss_spline_free (p->regionm[n].spline);
      }
      free (p->regionm);
    }
  }
  else
    p->region == 0 || DIE;
  free (p->peak_region);
  free (p->fdata);
  free (p->done);
  free (p);
}

static RegionData *
new_region_data (/*@only@*/ float *buf, int count, int starting_sample)
{
  RegionData *p;
  TCALLOC (p, 1);
  TCALLOC (p->region, (size_t)(p->region_alloc = 1));
  p->region_count = 1;
  p->region[0].sample = starting_sample;
  p->region[0].count = count;
  p->fdata = buf;
  return p;
}

static double
ssdmp (float *v, int len, Center *cnp, double start)
{
  int istart = (int)floor (start + .5), left, right;
  int e = cnp->peak_count - 1, i;
  double sum, tmp;

  left = istart + cnp->peak[0].loc - cnp->peak[0].lnear;
  right = istart + cnp->peak[e].loc + cnp->peak[e].rnear;
  if (left < 0)
    left = 0;
  if (right >= len)
    right = len - 1;
  sum = 0;
  for (i = left; i <= right; i++) {
    tmp = v[i] - mean;
    if (justone_region)
      printf ("ssdmp %d: %f\n", i, tmp);
    sum += tmp * tmp;
  }
  tmp = unzapped_sd * unzapped_sd * (right - left + 1);
  if (justone_region)
    printf ("ssdmp: sum %f, tmp %f, diff %f\n", sum, tmp, sum - tmp);
  if (sum > tmp)
    return sum - tmp;
  return 1;
}

static double 
ssdm (float *v, int len, double start)
{
  int istart, count, n;
  double sum, tmp;

  istart = (int)floor (start + .5);
  count = SNAPLEN;
  if (istart < 0) {
    count += istart; 
    istart = 0;
  }
  if (istart + count > len)
    count = len - istart;
  for (sum = 0, n = 0; n < count; n++) {
    tmp = v[istart+n] - mean;
    sum += tmp * tmp;
  }
  return sum;
}

static double
get_amp_ratio (float *buf, Center *cnp, double start)
{
  int center_peakloc = cnp->peak[0].loc;
  int region_peakloc = (int)floor (center_peakloc + start + .5);
  double center_peakval = ss_cspline_eval (cnp->spline, region_peakloc - start, cnp->acc) - mean;
  double region_peakval = buf[region_peakloc] - mean;
  double amp_ratio = region_peakval / center_peakval;
  double frac = .50;

  if (amp_ratio > 1+frac)
    amp_ratio = 1+frac;
  else if (amp_ratio < 1-frac)
    amp_ratio = 1-frac;

  return amp_ratio;             /* debug */

  if (fabs (region_peakval) < CLIPVAL)
    return 1;
  return amp_ratio;
}

bool use_best_attempt_flag;
static Component *attempt_component;
static int attempt_component_count;
static int attempt_component_alloc;
static float *attempt_residue;
static double min_residue_sum;

static void
record_attempt (DecomposeData *dd, /*@unique@*/ float *buf, Component *cp)
{
  int e, rcnt, n;
  double sum;

  if (!use_best_attempt_flag || dd->decomposition_count > 0)
    return;
  rcnt = dd->p->region[0].count;
  for (sum = 0, n = 0; n < rcnt; n++)
    sum += (buf[n] - mean) * (buf[n] - mean);
  if (min_residue_sum != 0 && sum >= min_residue_sum)
    return;
  min_residue_sum = sum;

  if (attempt_residue == 0)
    TMALLOC (attempt_residue, maxregion);
  if ((e = dd->component_count) + 1 > attempt_component_alloc)
    TREALLOC (attempt_component, attempt_component_alloc = dd->component_count + 1);
  TMEMCPY (attempt_component, dd->component, e);
  attempt_component[e].ccn        = cp->ccn;
  attempt_component[e].peak_count = cp->peak_count;
  attempt_component[e].start      = cp->start;
  attempt_component[e].scale      = dd->seed.scale;
  attempt_component_count = e + 1;
  TMEMCPY (attempt_residue, buf, rcnt);
}

static void
use_best_attempt (DecomposeData *dd)
{
  if (dd->component_alloc < attempt_component_count)
    TREALLOC (dd->component, dd->component_alloc = attempt_component_count);
  TMEMCPY (dd->component, attempt_component, dd->component_count = attempt_component_count);
  dd->region_data_copy = attempt_residue;
}

static void
attempt_to_octave (DecomposeData *dd)
{
  int rcnt, n, ccn, cpn;
  char *name;
  static int alloc_count;
  static float **w;
  static int *ccns;
  Center *ctr;
  double start;

  if (justone_region == 0)
    return;
  rcnt = dd->p->region[0].count;
  if (alloc_count < dd->component_count) {
    TREALLOC (w, dd->component_count);
    for (n = alloc_count; n < dd->component_count; n++)
      TMALLOC (w[n], maxregion);
    alloc_count = dd->component_count;
    printf ("alloc_count = %d\n", alloc_count);
    TREALLOC (ccns, dd->component_count);
  }
  for (cpn = 0; cpn < dd->component_count; cpn++) {
    ccns[cpn] = ccn = dd->component[cpn].ccn;
    ctr = &dd->center[ccn];
    start = dd->component[cpn].start;
    printf ("%s line %d: component: %d, scale: %f, val: %f, val0: %f\n",
            __FILE__, __LINE__, cpn, dd->component[cpn].scale,
            ss_cspline_eval (ctr->spline, 0 - start, ctr->acc),
            ss_cspline_eval (ctr->spline, 0, ctr->acc));
    for (n = 0; n < rcnt; n++)
      w[cpn][n] = dd->component[cpn].scale * ss_cspline_eval (ctr->spline, n - start, ctr->acc);
  }
  mfvec ("residue", attempt_residue, rcnt);
  for (cpn = 0; cpn < dd->component_count; cpn++) {
    if (asprintf (&name, "ctr%d_%d", ccns[cpn], cpn) == -1) exit (1);
    mfvec (name, w[cpn], rcnt);
    free (name);
  }
  mivec ("ccns", ccns, attempt_component_count);
  mcvec ("reason", reason);
}

static RegionData *
subtract_center (DecomposeData *dd, int ccn, Component *cp)
{
  double start = cp->start;
  int rcnt, n, buf_start, preexisting;
  Center *cnp;
  float *buf;
  RegionData *p;
#ifdef JUSTONE
  char name[80];
#endif
  static int seq;
  static float *v;
  double ssdm_before, ssdm_after;
  double ssdmp_before, ssdmp_after;

  if (!v)
    TMALLOC (v, maxregion);

  seq++;
  rcnt = dd->p->region[0].count;
  TMALLOC (buf, rcnt);
  TMEMCPY (buf, dd->p->fdata, rcnt);
  cnp = &dd->center[ccn];

  if (justone_region)
  {
    int n;
    char name[80];
    snprintf (name, 80, "before%d", seq);
    mfvec (name, buf, rcnt);
    printf ("%d: before %d %d:", seq, ccn, do_scaled_distance);
    for (n = 0; n < dd->p->peak_count; n++)
      printf (" %d", dd->p->peak[n].loc - dd->p->region[0].sample + 1);
    printf ("\n");
    printf ("seq %d:", seq);
    for (n = 0; n < dd->component_count; n++)
      printf (" %d", dd->component[n].ccn);
    printf (", %d\n", ccn);
  }

  ssdm_before = ssdm (buf, rcnt, start);
  ssdmp_before = ssdmp (buf, rcnt, cnp, start);
  {
    double amp_ratio = get_amp_ratio (buf, cnp, start);

#ifdef JUSTONE
    printf ("%s line %d: amp_ratio %f vs. %f\n", __FILE__, __LINE__, amp_ratio, dd->seed.scale);
#endif
    if (do_scaled_distance && dd->seed.scale && amp_ratio / dd->seed.scale < 1.5) {
      //      if (dd->seed.scale > 0)
        amp_ratio = dd->seed.scale;
    }
    else amp_ratio = 1;

    dd->seed.scale = amp_ratio;
    {
      int istart = (int)ceil (start);
      int match = 1;
      float bw[SNAPLEN];
      if (dd->flat) {
        if (istart + SNAPLEN > rcnt)
          match = 0;
        else
          TMEMCPY (bw, buf + istart, SNAPLEN);
      }
      for (n = 0; n < rcnt; n++)
        buf[n] -= v[n] = (amp_ratio * (ss_cspline_eval (cnp->spline, n - start, cnp->acc) - mean));
      record_attempt (dd, buf, cp);
      if (dd->flat && !match)
        if (!reason) if (asprintf (&reason, "flat match center runs off end of region") == -1) exit (1);
      if (dd->flat && match) {
        double d;
        if ((d = white_distance_from_zero (buf + istart, SNAPLEN)) > 11.5) {
          if (!reason) if (asprintf (&reason, "flat match residue too far from noise: %f", d) == -1) exit (1);
          match = 0;
        }
        else {
          int i;
          double sum = 0, len, proj = 0;
          float cw[SNAPLEN];
          TMEMCPY (cw, v + istart, SNAPLEN);
          whiten_snap (cw);
          for (i = 0; i < SNAPLEN; i++)
            sum += cw[i] * cw[i];
          len = sqrt (sum);
          for (i = 0; i < SNAPLEN; i++)
            cw[i] /= len;
          whiten_snap (bw);
          for (i = 0; i < SNAPLEN; i++)
            proj += bw[i] * cw[i];
          if (proj < 0)
            match = 0;
          else if (proj < len) {
            double noise_expect = gsl_sf_erf_Q (proj);
            double cluster_expect = gsl_sf_erf_Q (len - proj);
            if (cluster_expect < 1000 * noise_expect)
              match = 0;
          }
          if (!reason && !match) if (asprintf (&reason, "flat match residue closer to noise than center") == -1) exit (1);
        }
      }
      if (!match) {
        free (buf);
        if (justone_region)
          printf ("seq %d subtract_center failed: %s\n", seq, reason);
        return 0;
      }
    }
  }
  ssdm_after = ssdm (buf, rcnt, start);
  ssdmp_after = ssdmp (buf, rcnt, cnp, start);

  if (justone_region) {
    char name[80];
    for (n = 0; n < rcnt; n++)
      v[n] += mean;
    printf ("seq %d center %d amp_ratio %f\n", seq, ccn, dd->seed.scale);
    snprintf (name, 80, "center%d", seq);
    mfvec (name, v, rcnt);
  }

  p = new_region_data (buf, rcnt, dd->p->region[0].sample);
  TMALLOC (p->done, rcnt);
  TMEMCPY (p->done, dd->p->done, rcnt);
  find_mpk_spikes (p, buf, rcnt, dd->p->region[0].sample, buf_start = 0, preexisting = 1);
  region_check (p);

  if (justone_region)
  {
    int n;
    char name[80];
    snprintf (name, 80, "after%d", seq);
    mfvec (name, buf, rcnt);
    printf ("%d: after %d:", seq, ccn);
    for (n = 0; n < p->peak_count; n++)
      printf (" %d", p->peak[n].loc - p->region[0].sample + 1);
    printf ("\n");
  }

  fdata_spline_init (p);
  if (SHOWSUB)
    printf ("center %d: before %d, after %d\n", ccn,  dd->p->peak_count, p->peak_count);

  if (ssdmp_after < ssdmp_before / 2) {
    if (justone_region)
      printf ("seq %d subtract_center succeeded\n", seq);
    return p;
  }

  if (justone_region) {
    printf ("center %d: no fewer peaks, ssdm: %f -> %f: %f\n", ccn, ssdm_before, ssdm_after, ssdm_before / ssdm_after);
    printf ("center %d: no fewer peaks, ssdmp: %f -> %f: %f\n", ccn, ssdmp_before, ssdmp_after, ssdmp_before/ ssdmp_after);
  }

  free_region_data (p);
  if (!reason) if (asprintf (&reason, "not enough progress") == -1) exit (1);
  if (justone_region) printf ("seq %d subtract_center failed: %s\n", seq, reason);
  return 0;
}

static int
max_peak (RegionData *p)
{
  int pn, max_pn = 0;
  double amp, max_amp = 0;

  for (pn = 0; pn < p->peak_count; pn++) {
    if ((amp = fabs (p->fdata[p->peak[pn].loc - p->region[0].sample] - mean)) > max_amp) {
#ifdef JUSTONE
      printf ("amp = %f\n", amp);
#endif
      max_amp = amp;
      max_pn = pn;
    }
  }
  return max_pn;
}

static int
flat_match_orig (DecomposeData *dd, Component *cp)
{
  int n0, match = 0;
  RegionData *p = dd->p;
  RegionData *p0 = dd->p0;
  int start = (int)ceil (cp->start) + p->region[0].sample;
  int end = start + SNAPLEN;
  if (dd->p->region[0].sample == SAMPLOC)
    printf ("start %d, end %d\n", start, end);
  for (n0 = 0; n0 < p0->peak_count; n0++) {
    if (dd->p->region[0].sample == SAMPLOC)
      printf ("peak %d at %d done %d\n", n0, p0->peak[n0].loc, p->done[n0]);
    if (p0->peak[n0].loc >= start && p0->peak[n0].loc < end) {
      if (p->done[n0] && dd->flat != 2)
        return 0;
      match = 1;
    }
  }
  return match;
}

static int
matches_remaining_orig_peaks (DecomposeData *dd, Component *cp)
{
  int pcnt = cp->peak_count;
  int n, n0, new_match, any_match;
  RegionData *p = dd->p;
  RegionData *p0 = dd->p0;

  if (dd->flat)
    return flat_match_orig (dd, cp);
  any_match = new_match = 0;
  for (n = dd->pn; n < dd->pn + pcnt; n++) {
    for (n0 = 0; n0 < p0->peak_count; n0++)
      if (p->peak[n].loc >= p0->peak[n0].loc - p0->peak[n0].lnear
          && p->peak[n].loc <= p0->peak[n0].loc + p0->peak[n0].rnear)
        break;
    if (n0 < p0->peak_count) {
      any_match = 1;
      if (p->done[n0])
        return 0;
      else
        new_match = 1;
    }
  }
  any_match = 0;
  for (n = 0; n < p->peak_count; n++) {
    for (n0 = 0; n0 < p0->peak_count; n0++)
      if (p->peak[n].loc + p->peak[n].rfar >= p0->peak[n0].loc - p0->peak[n0].lfar
          && p->peak[n].loc - p->peak[n].lfar <= p0->peak[n0].loc + p0->peak[n0].rfar)
        break;
    if (n0 < p0->peak_count) {
      any_match = 1;
      break;
    }
  }

#ifdef JUSTONE
  printf ("new match: %d, any_match: %d\n", new_match, any_match);
#endif
  return new_match || !any_match;
}

static void
flat_mark_peaks (DecomposeData *dd, double dstart, RegionData *newp)
{
  int n0;
  RegionData *p = dd->p;
  RegionData *p0 = dd->p0;
  int start = (int)ceil (dstart) + p->region[0].sample;
  int end = start + SNAPLEN;

  for (n0 = 0; n0 < p0->peak_count; n0++)
    if (p0->peak[n0].loc >= start && p0->peak[n0].loc < end)
      newp->done[n0] = 1;
}

static void
mark_peaks (DecomposeData *dd, Peak *peak, int pcnt, double start, RegionData *newp)
{
  int n, n0;
  RegionData *p = dd->p;
  RegionData *p0 = dd->p0;

  if (dd->flat) {
    flat_mark_peaks (dd, start, newp);
    return;
  }
  start += p->region[0].sample;
  for (n = 0; n < pcnt; n++) {
    for (n0 = 0; n0 < p0->peak_count; n0++)
      if (peak[n].loc + start >= (double)p0->peak[n0].loc - p0->peak[n0].lnear
          && peak[n].loc + start <= (double)p0->peak[n0].loc + p0->peak[n0].rnear)
        break;
    if (n0 < p0->peak_count) {
      if (justone_region)
        printf ("mark_peaks: number %d at %d done (of %d)\n",
                n0,
                p0->peak[n0].loc - p->region[0].sample,
                p0->peak_count);
      newp->done[n0] = 1;
    }
  }
}

static int
push_center (DecomposeData *dd, int ccn, int pnc0)
{
  Component cp;
  int cpn, loose;
  float v[SNAPLEN];
  RegionData *p;

#ifdef JUSTONE
  printf ("push_center: region peak %d at %d, center %d, center peak %d %s\n",
          dd->pn, dd->p->peak[dd->pn].loc - dd->p->region[0].sample, ccn, pnc0, (do_scaled_distance ? "scaled" : "unscaled"));
#endif

  if (center_peak_polarity (dd, ccn, pnc0) != region_peak_polarity (dd)) {
#ifdef JUSTONE
  printf ("push_center failed: polarity differs\n");
#endif
    return 0;
  }
  dd->seed.type = dd->center[ccn].type;
  dd->seed.peak = dd->center[ccn].peak;
  dd->seed.cdata = dd->center[ccn].raw;
  dd->seed.val = v;
  dd->seed.ccn = ccn;
  cp.ccn = ccn;
  if ((cp.peak_count = region_center_peaks (dd->p, dd->pn, &dd->seed, pnc0, &cp.start, loose = 0)) != 0
      && matches_remaining_orig_peaks (dd, &cp)
      && (p = subtract_center (dd, ccn, &cp))) {
    mark_peaks (dd, dd->seed.peak, dd->seed.type / 2 + 1, cp.start, p);
#ifdef BIGFIRST
    dd->pn = max_peak (p);
#else
    dd->pn = 0;
#endif
    dd->rn = 0;
    cpn = dd->component_count++;
    if (dd->component_count > dd->component_alloc) 
      TREALLOC (dd->component, dd->component_alloc = dd->component_count);
    cp.scale = dd->seed.scale;
    cp.overlapped = 0;
    cp.p = dd->p;
    dd->p = p;
    dd->region_data_copy = p->fdata;
    dd->component[cpn] = cp;
    if (justone_region || global_rn == debug_region) {
      static int count;
      printf ("pushpop: push center %d at %f succeeded\n", cp.ccn, cp.start);
      if (0) {
        printf ("%d push_center succeeded:  %d peaks\n", count++, cp.peak_count);
        for (cpn = 0; cpn < dd->component_count; cpn++)
          printf (" %d", dd->component[cpn].ccn);
        printf ("\n");
      }
    }
    dd->seed.val = 0;
    return 1;
  }
  else {
    if (0)
      if (justone_region == -1)
        printf ("push_center failed:  %s\n", cp.peak_count == 0 ? "no peaks match" : (matches_remaining_orig_peaks (dd, &cp) ? "subtract failed" : "doesn't match remaining orig peaks"));
  }
  dd->seed.val = 0;
  return 0;
}

static void
pop_center (DecomposeData *dd)
{
  int cpn;

  if (justone_region) printf ("pushpop: pop center\n");
  cpn = --dd->component_count;
#ifdef JUSTONE
  printf ("pop center %d\n", dd->component[cpn].ccn);
#endif
  free_region_data (dd->p);
  dd->p = dd->component[cpn].p;
  dd->region_data_copy = dd->component[cpn].p->fdata;
}

static double
component_distance (DecomposeData *dd, int cpn)
{
  int istart, n;
  float w[SNAPLEN], *v;
  static float zero[SNAPLEN], vmean[SNAPLEN];
  double d;

  if (vmean[0] == 0)
    for (n = 0; n < SNAPLEN; n++)
      vmean[n] = mean;
  istart = (int)floor (dd->component[cpn].start + .5);
  v = dd->region_data_copy + istart;
  SZMEMCPY (w, v);
  whiten_snap (w);
  d = distance_samples (v, w, vmean, zero, SNAPLEN);
#ifdef JUSTONE
  printf ("%d %3d %f\n", dd->component[cpn].ccn, istart, d);
#endif
  dd->component[cpn].distance = d;
  if (0)
  {
    static int c;
    char msg[256];
    window ();
    if (c != 'c') {
      clear ();
      draw_snap (dd->region_data_copy + istart, "black");
      snprintf (msg, 256, "region %d starting at %d from mean: %f", dd->rn, istart, d);
      draw_text (msg);
      c = show ();
      if (c == 'q')
        exit (0);
    }
    if (c == 'c')
      wclose ();
  }
  return d;
}

static void
residue_spline (DecomposeData *dd)
{
  static double *x, *y;
  static int alloc;
  int count, n;

  count = dd->p->region[dd->rn].count;
  if (alloc < count) {
    TREALLOC (x, count);
    TREALLOC (y, count);
    alloc = count;
  }

  if (dd->acc && dd->spline) {
    gsl_interp_accel_free (dd->acc);
    ss_spline_free (dd->spline);
  }
  else if (dd->acc == 0 && dd->spline == 0) {
  }
  else exit (DIE);

  dd->acc = gsl_interp_accel_alloc ();
  (count >= 0 && count <= maxregion) || DIE;
  dd->spline = gsl_spline_alloc (gsl_interp_cspline, (size_t)count);

  for (n = 0; n < count; n++) {
    x[n] = (double)n;
    y[n] = dd->region_data_copy[n];
  }

  gsl_spline_init (dd->spline, x, y, (size_t)count);
}

#define SNAPLEN2 (SNAPLEN*2)

void
do_shift_snap_128 (float *v, double shift_amt)
{
  static rfftw_plan fwd_snaplen2, rev_snaplen2;
  double real, imag, w, m, c, s;
  double tdat[SNAPLEN2], fdat[SNAPLEN2];
  int n;

  if (fwd_snaplen2 == 0) {
    fwd_snaplen2 = rfftw_create_plan(SNAPLEN2, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);
    rev_snaplen2 = rfftw_create_plan(SNAPLEN2, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
  }

  m = 2*M_PI*shift_amt / SNAPLEN2;
  for (n = 0; n < SNAPLEN2; n++)
    tdat[n] = v[n];
  rfftw_one(fwd_snaplen2, tdat, fdat);
  for (n = 1; n <= SNAPLEN2/2; n++) {
    w = m * n;
    real = fdat[n];
    imag = fdat[SNAPLEN2-n];
    fdat[n] = real*(c=cos(w)) - imag*(s=sin(w));
    fdat[SNAPLEN2-n] = real*s + imag*c;
  }
  rfftw_one(rev_snaplen2, fdat, tdat);
  for (n = 0; n < SNAPLEN2; n++)
    v[n] = tdat[n]/SNAPLEN2;
}

void
shift_snap_128 (float *v, double shift)
{
  int n;
  float vx_array[SNAPLEN*2];
  float *vx = vx_array + SNAPLEN/2;

  for (n = -SNAPLEN / 2; n < SNAPLEN + SNAPLEN/2; n++)
    vx[n] = extrap_val (v, (double)n);
  do_shift_snap_128 (vx_array, shift);
  for (n = 0; n < SNAPLEN; n++)
    v[n] = vx[n];
}

static void
get_component_raw (DecomposeData *dd, int cpn)
{
  int n, count;
  double start;
  Component *comp;
  float *raw, *craw;
  float v[SNAPLEN];
  float vs[SNAPLEN];
  int istart;
  double shift;

  comp = dd->component + cpn;
  start = comp->start;
  raw = comp->raw;
  craw = dd->center[comp->ccn].raw;
  comp->scale || DIE;
  count = dd->p->region[dd->rn].count;

  istart = floor (start + .5);
  shift = start - istart;
  (istart >= 0 && istart <= count - SNAPLEN) || DIE;
  for (n = 0; n < SNAPLEN; n++)
    vs[n] = dd->region_data_copy[istart + n];
  shift_snap_128 (vs, shift);

  for (n = 0; n < SNAPLEN; n++)
    if (start + n < 0 || start + n > (double)count - 1)
      raw[n] = comp->scale * (craw[n] - mean) + mean;
    else
      //      raw[n] = comp->scale * (craw[n] - mean) + (v[n] = ss_cspline_eval (dd->spline, start + n, dd->acc));
      raw[n] = comp->scale * (craw[n] - mean) + (v[n] = vs[n]);
  if (justone_region) {
    char *name;

    printf ("component %f from noise, residue is %f (%f) starting at %f\n",
            white_distance_from_zero (raw, SNAPLEN),
            white_distance_from_zero (v, SNAPLEN),
            white_distance_from_zero (vs, SNAPLEN),
            start);

    if (asprintf (&name, "calign%d", cpn) == -1) exit (1); mfvec (name, v, SNAPLEN); free (name);
    whiten_snap (v);
    if (asprintf (&name, "calign%dw", cpn) == -1) exit (1); mfvec (name, v, SNAPLEN); free (name);

    if (asprintf (&name, "calign%ds", cpn) == -1) exit (1); mfvec (name, vs, SNAPLEN); free (name);
    whiten_snap (vs);
    if (asprintf (&name, "calign%dws", cpn) == -1) exit (1); mfvec (name, vs, SNAPLEN); free (name);

    start = floor (start + .5);
    for (n = 0; n < SNAPLEN; n++)
      v[n] = ss_cspline_eval (dd->spline, start + n, dd->acc);

    printf ("residue is %f starting at %f\n",
            white_distance_from_zero (v, SNAPLEN), start);

    if (asprintf (&name, "ralign%d", cpn) == -1) exit (1); mfvec (name, v, SNAPLEN); free (name);
    whiten_snap (v);
    if (asprintf (&name, "ralign%dw", cpn) == -1) exit (1); mfvec (name, v, SNAPLEN); free (name);

  }
}

void
components_to_octave (DecomposeData *dd)
{
  char *name;
  int cpn, n;
  double start;
  Component *comp;
  Center *ctr;
  static float *buf;
  int count = dd->p->region[dd->rn].count;
  int dn = dd->decomposition_count;
  static int callcnt;

  callcnt++;

  TREALLOC (buf, count);
  if (asprintf (&name, "residue%d_%d", callcnt, dn) == -1) exit (1);
  mfvec (name, dd->region_data_copy, count);
  free (name);
  for (cpn = 0; cpn < dd->component_count; cpn++) {
    comp = dd->component + cpn;
    start = -comp->start;
    ctr = &dd->center[comp->ccn];
    for (n = 0; n < count; n++)
      buf[n] = comp->scale * (ss_cspline_eval (ctr->spline, start + n, ctr->acc) - mean);
    if (asprintf (&name, "component%d_%d_%d", callcnt, dn, cpn) == -1) exit (1);
    mfvec (name, buf, count);
    free (name);
  }
}

static double
snap_corr (float *a, float *b)
{
  int n;
  double sum = 0;

  for (n = 0; n < SNAPLEN; n++)
    sum += a[n] * b[n];
  return sum / (SNAPLEN - 1);
}

static double
component_residue_correlation (DecomposeData *dd, int cpn)
{
  int n, count;
  double start;
  Component *comp;
  float *craw;
  float v[SNAPLEN];
  float vs[SNAPLEN];
  int istart;
  double shift;

  comp = dd->component + cpn;
  start = comp->start;
  craw = dd->center[comp->ccn].raw;
  comp->scale || DIE;
  count = dd->p->region[dd->rn].count;

  istart = floor (start + .5);
  shift = start - istart;
  (istart >= 0 && istart <= count - SNAPLEN) || DIE;
  for (n = 0; n < SNAPLEN; n++)
    vs[n] = dd->region_data_copy[istart + n];
  shift_snap_128 (vs, shift);

  for (n = 0; n < SNAPLEN; n++)
    v[n] = mean + comp->scale * (craw[n] - mean);
  whiten_snap (v);
  whiten_snap (vs);
  return snap_corr (v, vs);
}

static int dbl_cmp_up (const double *a, const double *b)
{
  return *a < *b ? -1 : (*a > *b); /* ascending */
}

static int
too_close (DecomposeData *dd)
{
  static double *start;
  static int start_alloc;
  int cpn;

  if (dd->component_count < 3)
    return 0;
  if (start_alloc < dd->component_count)
    TREALLOC (start, start_alloc = dd->component_count);
  for (cpn = 0; cpn < dd->component_count; cpn++)
    start[cpn] = dd->component[cpn].start;
  qsort (start, (size_t)dd->component_count, sizeof *start, (int(*)(const void*,const void*))dbl_cmp_up);
  for (cpn = 0; cpn + 2 < dd->component_count; cpn++)
    if (start[cpn + 2] - start[cpn] <= SF / 2000.0)
      return 1;
  return 0;
}

static int
try_distance (DecomposeData *dd, int force)
{
  int cpn, dn, retval = 0;
  double d, max, corr, max_corr;

  dd->component_count > 0 || DIE;
  max = 0;
  max_corr = 0;
  for (cpn = 0; cpn < dd->component_count; cpn++) {
    d = component_distance (dd, cpn);
    if (SHOWSUB)
      printf (" %d", dd->component[cpn].ccn);
    if (!force) dd->component[cpn].overlapped = 1;
    if (d > max)
      max = d;
    if (!force) {
      corr = fabs (component_residue_correlation (dd, cpn));
      if (corr > max_corr)
        max_corr = corr;
    }
  }
  if (SHOWSUB)
    printf (": %f\n", max);

  if (force || (max_corr < .5 && max < 10 && !too_close (dd))) {
    //  if (1) {
    if (justone_region) {
      printf ("match: %f (%s scaling)\n", max, do_scaled_distance ? "with" : "without");
      components_to_octave (dd);
    }
      
    dn = dd->decomposition_count++;
    if (dn == 0 || max < dd->min_distance) {
      dd->min_distance = max;
      residue_spline (dd);
      for (cpn = 0; cpn < dd->component_count; cpn++)
        get_component_raw (dd, cpn);
    }
    if (dd->decomposition_count > dd->decomposition_alloc){
      TREALLOC (dd->decomposition, dd->decomposition_alloc = dd->decomposition_count);
      dd->decomposition[dn].component = 0;
    }
    TREALLOC (dd->decomposition[dn].component, dd->component_count);
    TMEMCPY (dd->decomposition[dn].component, dd->component, dd->component_count);
    dd->decomposition[dn].component_count = dd->component_count;
    dd->decomposition[dn].residue = max;
    retval = 1;
  }
  for (cpn = 0; cpn < dd->component_count; cpn++)
    dd->component[cpn].overlapped = 0;
  if (justone_region && retval == 0)
    printf ("try_distance failed: !(max_corr %f < .5 && max %f < 10 && !too_close)\n", max_corr, max);
  return retval;
}

#ifdef BIGFIRST

static int
look_back (DecomposeData *dd, int ccn, int pnc0)
{
  int pnc, ddpn, offset, r0;
  RegionRange range;

  pnc = pnc0;
  ddpn = dd->pn;
  offset = 1;
  r0 = dd->p->region[dd->rn].sample;
  range.start = 0;
  range.end = INT_MAX;

  while (1) {
    if (center_peak_polarity (dd, ccn, pnc) != region_peak_polarity (dd)) {
#ifdef JUSTONE
      printf ("look_back: center %d peak %d polarity %d, region peak %d polarity %d\n",
              ccn, pnc, center_peak_polarity (dd, ccn, pnc), dd->pn, region_peak_polarity (dd));
#endif
      break;
    }
    region_intersect (&range, peak_range (dd->p->peak + dd->pn, dd->center[ccn].peak + pnc, r0));
    if (range.start > range.end) {
#ifdef JUSTONE
      Peak *pr = dd->p->peak + dd->pn, *pc = dd->center[ccn].peak + pnc;
      printf ("look_back: %d %d %d, %d %d %d, %d\n", pr->lnear, pr->loc, pr->rnear, pc->lnear, pc->loc, pc->rnear, r0);
#endif
      break;
    }
    offset--;
    if (dd->pn == 0 || pnc == 0)
      break;
    dd->pn--;
    pnc--;
  }
  dd->pn = ddpn;
#ifdef JUSTONE
  printf ("look_back region_peak %d, center %d, center peak %d: offset %d\n", dd->pn, ccn, pnc0, offset);
#endif
  return offset;
}

static int
all_peaks_noise (DecomposeData *dd)
{
  RegionData *p = dd->p;
  int start, n, pn;
  float w[SNAPLEN], v[SNAPLEN];
  static float zero[SNAPLEN], vmean[SNAPLEN];
  double d;

  if (vmean[0] == 0)
    for (n = 0; n < SNAPLEN; n++)
      vmean[n] = mean;
  //  printf ("all_peaks_noise\n");
  for (pn = 0; pn < p->peak_count; pn++) {
    //    start = p->peak[pn].loc - SNAPLEN / 2;
    start = (p->peak[pn].loc - p->region[0].sample) - PRESAMPLES;
    if (start < 0)
      start = 0;
    if (start + SNAPLEN > p->region[0].count)
      start = p->region[0].count - SNAPLEN;
    TMEMCPY (v, p->fdata + start, SNAPLEN);
    SZMEMCPY (w, v);
    whiten_snap (w);
    d = distance_samples (v, w, vmean, zero, SNAPLEN);
#ifdef JUSTONE
    printf ("peak %d, window at %d: distance %f\n", pn, start, d);
#endif
    if (d > 11.5)
      return 0;
  }
#ifdef JUSTONE
  printf ("all remaining peaks noise\n");
#endif
  return 1;
}

static int
next_peak (DecomposeData *dd)
{
  RegionData *p = dd->p;
  int last_pn = dd->pn;
  int pn, max_pn = -1;
  double amp, max_amp = 0, last_amp, last_val, val;

  last_val = p->fdata[p->peak[last_pn].loc - p->region[0].sample];
  last_amp = fabs (last_val - mean);

  if (0 && do_scaled_distance) {
    for (pn = 0; pn < p->peak_count; pn++) {
      val = p->fdata[p->peak[pn].loc - p->region[0].sample];
      if (val > dd->good_lo && val < dd->good_hi)
        printf ("peak at %d\n", p->peak[pn].loc);
      if (val > dd->good_lo && val < dd->good_hi
          && (amp = fabs (val - mean)) > max_amp
          && (amp < last_amp || (amp == last_amp && pn > last_pn))) {
        max_amp = amp;
        max_pn = pn;
      }
    }
    if (max_pn != -1) {
      printf ("next_peak scaled at %d (last was at %d)\n", p->peak[max_pn].loc, p->peak[last_pn].loc);
      return max_pn;
    }
  }

  for (pn = 0; pn < p->peak_count; pn++) {
    if ((amp = fabs (p->fdata[p->peak[pn].loc - p->region[0].sample] - mean)) > max_amp
        && (amp < last_amp || (amp == last_amp && pn > last_pn))) {
#ifdef JUSTONE
      printf ("amp = %f\n", amp);
#endif
      max_amp = amp;
      max_pn = pn;
    }
  }
  if (justone_region) printf ("last_pn: %2d, last_amp: %11.6f at %d, max_pn: %2d, max_amp: %11.6f at %d\n",
                              last_pn, last_amp,
                              p->peak[last_pn].loc - p->region[0].sample + 1,
                              max_pn, max_amp,
                              p->peak[max_pn].loc - p->region[0].sample + 1
                              );
  return max_pn;
}
static int level;

typedef struct {int ccn; double start; double scale; } CtrLoc;

static int cl_cmp_up (const CtrLoc *a, const CtrLoc *b)
{
  return a->start < b->start ? -1 : (a->start > b->start); /* ascending */
}

static void
list_components (DecomposeData *dd)
{
  static int alloc;
  static CtrLoc *ctrloc;
  int n;

  if (do_scaled_distance == 0)
    return;
  if (alloc < dd->component_count)
    TREALLOC (ctrloc, alloc = dd->component_count);
  for (n = 0; n < dd->component_count; n++) {
    ctrloc[n].ccn = dd->component[n].ccn;
    ctrloc[n].start = dd->component[n].start;
    ctrloc[n].scale = dd->component[n].scale;
  }
  if (0) qsort (ctrloc, (size_t)dd->component_count, sizeof *ctrloc, (int(*)(const void*,const void*))cl_cmp_up);
  printf ("components %2d:", dd->component_count);
  for (n = 0; n < dd->component_count; n++) {
    //    printf (" %d %.1f %.1f", ctrloc[n].ccn, ctrloc[n].start, ctrloc[n].scale);
    int ccn = ctrloc[n].ccn;
    int pn;
    for (pn = 0; pn < dd->center[ccn].peak_count; pn++) {
      int ploc = dd->center[ccn].peak[pn].loc;
      printf (" %.0f", dd->p->region[0].sample + ctrloc[n].start + ploc);
    }
    printf (",");
  }
  printf ("\n");
}

static double
max_from_mean (float *v, int len)
{
  double max, d;
  int n;
  
  max = 0;
  for (n = 0; n < len; n++)
    if ((d = fabs (v[n] - mean)) > max)
      max = d;
  return max;
}

static int
match_at_windows (DecomposeData *dd)
{
  int cpn, ccn;
  double dr, dc, min, max, ratio;

  dr = max_from_mean (dd->region_data_copy, dd->p->region[0].count);
  min = DBL_MAX;
  max = 0;
  for (cpn = 0; cpn < dd->component_count; cpn++) {
    ccn = dd->component[cpn].ccn;
    dc = max_from_mean (dd->center[ccn].raw, SNAPLEN);
    ratio = dd->component[cpn].scale * dc / dr;
    if (ratio < min)
      min = ratio;
    if (ratio > max)
      max = ratio;
  }
  if (0)
    if (dd->p->region[0].sample == 29714478)
      printf ("%157s\rmatch_at_windows: %f %d %d %d\n",
              "", min, dd->component_count, do_scaled_distance, global_rn);
  if (justone_region)
    printf ("match_at_windows: %f %f %d %d %d\n",
            min, max, dd->component_count, do_scaled_distance, global_rn);
  if (max < 4.6)
    return 0;
  return dd->component_count;
}

static int
good_center (Center *c, double good_lo, double good_hi)
{
  int pidx;
  double peak_2val;

  for (pidx = 0; pidx < c->peak_count; pidx++) {
    peak_2val = c->raw[c->peak[pidx].loc] * 2;
    if (justone_region)
      printf ("good_center peak_2val %f\n", peak_2val);
    if (peak_2val >= good_lo && peak_2val <= good_hi)
      return c->peak[pidx].loc;
  }
  return 0;
}

static int
center_width (Center *c)
{
  Peak *p0, *pe;
  p0 = c->peak;
  pe = c->peak + c->peak_count - 1;
  return (pe->loc - p0->loc) + p0->lnear + pe->rnear;
}

static int
min_peak_spacing (DecomposeData *dd, double good_lo, double good_hi)
{
  int w, ccn, *cpkloc;
  Center *c;
  double max;

  TMALLOC (cpkloc, dd->center_count);
  w = 0;
  max = 0;
  for (ccn = 0; ccn < dd->center_count; ccn++) {
    c = dd->center + ccn;
    if (justone_region)
        printf ("good_center ccn %d\n", ccn);
    if ((cpkloc[ccn] = good_center (c, good_lo, good_hi))
        && (w = center_width (c)) > max)
      max = w;
  }
  if (justone_region)
    for (ccn = 0; ccn < dd->center_count; ccn++)
      printf ("good_center ccn %d, cpkloc %d (%g %g)\n", ccn, cpkloc[ccn], good_lo, good_hi);

  free (dd->cpkloc);
  dd->cpkloc = cpkloc;
  return w;
}

static int float_cmp_up (const float *a, const float *b)
{
  return *a < *b ? -1 : (*a > *b); /* ascending */
}

#define PK(t,n,c) (v[n] c t && v[n] c v[n-1] && (v[n] c v[n+1] || (v[n] == v[n+1] && v[n] c v[n+2])))

static int
find_clean_peaks (DecomposeData *dd)
{
  int n, pcnt, n_at_max = 0, pidx, pkmin, bpcnt, valley;
  static float *pval;
  double max, threshold, good_lo, good_hi, valley_lo, valley_hi;
  int *pkloc;
  float *v = dd->region_data_copy;
  int count = dd->p->region[0].count;

  dd->good_lo = 1;
  dd->good_hi = -1;
  free (dd->pkloc);
  dd->pkloc = 0;
  TREALLOC (pval, count);
  pcnt = 0;
  for (n = 1; n < count - 2; n++)
    if (PK(mean, n, >) || PK (mean, n, <))
      pval[pcnt++] = v[n];
  qsort (pval, (size_t)pcnt, sizeof *pval, (int(*)(const void*,const void*))float_cmp_up);
  max = 0;
  for (n = 1; n < pcnt; n++)
    if (pval[n] - pval[n-1] > max)
      max = pval[n] - pval[n-1],
        n_at_max = n;
  if (max < 3 * detection_threshold)
    return 0;
  if (justone_region)
    printf ("gap low: %f, gap high: %f\n", pval[n_at_max-1], pval[n_at_max]);
  if (pval[n_at_max] > mean) {
    threshold = pval[n_at_max] - 2 * detection_threshold;
    if (threshold < mean)
      return 0;
    good_lo = threshold;
    good_hi = DBL_MAX;
    valley_lo = -DBL_MAX;
    valley_hi = pval[n_at_max - 1];
    bpcnt = pcnt - n_at_max;
  }
  else {
    threshold = pval[n_at_max-1] + 2 * detection_threshold;
    if (threshold > mean)
      return 0;
    good_lo = -DBL_MAX;
    good_hi = threshold;
    valley_hi = DBL_MAX;
    valley_lo = pval[n_at_max];
    bpcnt = n_at_max;
  }
  TMALLOC (pkloc, bpcnt + 1);
  pkloc[0] = bpcnt;
  pidx = 1;

  valley = 0;
  for (n = 1; n < count - 2; n++) {
    if (valley == 0 && v[n] > valley_lo && v[n] < valley_hi)
      valley = 1;
    if ((threshold > mean && PK (threshold, n, >)) || (threshold < mean && PK (threshold, n, <))) {
      if (valley == 0 || pidx > bpcnt) {free (pkloc); return 0;}
      valley = 0;
      pkloc[pidx++] = n;
    }
  }
  pidx == bpcnt + 1 || DIE;
  
  pkmin = min_peak_spacing (dd, good_lo, good_hi);
  for (pidx = 1; pidx < bpcnt; pidx++)
    if (pkloc[pidx+1] - pkloc[pidx] < pkmin) {
      free (pkloc);
      if (justone_region) printf ("peak locs %d and %d: %d and %d, pkmin %d\n",
                                  pidx, pidx + 1, pkloc[pidx], pkloc[pidx+1], pkmin);
      return 0;
    }
  if (justone_region) {
    printf ("spike count: %d\n", spike_count);
    for (pidx = 1; pidx <= bpcnt; pidx++)
      printf ("clean peak at %d\n", dd->p->region[0].sample + pkloc[pidx]);
  }
  dd->pkloc = pkloc;
  dd->good_lo = good_lo;
  dd->good_hi = good_hi;
  return 1;
}

static int
match_any_component (DecomposeData *dd, int pkloc)
{
  int cpn, ccn;
  double d, dmin;
  Component *cmp;

  if (justone_region)
    printf ("match_any_component\n");
  dmin = maxregion;
  for (cpn = 0; cpn < dd->component_count; cpn++) {
    cmp = dd->component + cpn;
    ccn = cmp->ccn;
    if (justone_region)
      printf ("match_any_component: cpn %d, ccn %d, pkloc %d\n", cpn, ccn, dd->cpkloc[ccn]);
    d = fabs (cmp->start + dd->cpkloc[ccn] - pkloc);
    if (dd->cpkloc[ccn] && d < 1.75) {
      if (justone_region)
        printf ("match_any_component: cpn %d of %d, pkloc %d, dist %f\n",
                cpn, dd->component_count, dd->cpkloc[ccn], fabs (cmp->start + dd->cpkloc[ccn] - pkloc));
      return 1;
    }
    if (d < dmin)
      dmin = d;
  }
  if (justone_region)
    printf ("match_any_component: region pkloc %d, min dist %f\n", /*dd->p->region[0].sample + */pkloc, dmin);
  return 0;
}

static int
match_all_peaks (DecomposeData *dd)
{
  int bpn, bpcnt;

  if (justone_region)
    printf ("match_all_peaks dd->pkloc %ld, compcnt %d\n",
            (long)dd->pkloc, dd->component_count);
  if (dd->pkloc == 0)
    return 0;
  bpcnt = dd->pkloc[0];
  if (dd->component_count != bpcnt)
    return 0;
  for (bpn = 1; bpn <= bpcnt; bpn++)
    if (!match_any_component (dd, dd->pkloc[bpn])) {
      if (justone_region)
        printf ("match_all_peaks failed at %d of %d\n", bpn, bpcnt);
      return 0;
    }
  if (justone_region)
    printf ("match_all_peaks succeeded\n");
  return 1;
}

static int
decompose_region (DecomposeData *dd)
{
  int ccn, pn, ddpn, offset, td;
  int decomp_count = dd->decomposition_count;
  double now;
  int nopush;

  if (0 && justone_region) {
    int n;
    for (n = 0; n <= level; n++) putchar (' ');
    printf ("decompose_region level %d, dd->pn %d\n", level + 1, dd->pn);
  }
  now = (double)clock () / CLOCKS_PER_SEC;
  if (now > drstart + 1) {
    if (0 && justone_region) {
      int n;
      for (n = 0; n <= level; n++) putchar (' ');
      printf ("timeout level %d\n", level + 1);
    }
    return 0;
  }
  if (justone_region) list_components (dd);
  
  if (now < drstart)
    drstart = now;
  level++;
  if (0 && justone_region) {
    int n;
    for (n = 0; n < level; n++) putchar (' ');
    printf ("level %d, peak_count %d, drstart %f, now %f\n",
            level, dd->p->peak_count, drstart, now);
  }
  if (dd->p->peak_count == 0) {
    if (try_distance (dd, 0) && justone_region)
      printf ("decomposition at no peaks\n");
  }
  else if (all_peaks_noise (dd) && (td = try_distance (dd, 0))) {
    if (td && justone_region)
      printf ("decomposition at all peaks noise\n");
  }
  else if (0 && match_at_windows (dd) && try_distance (dd, 1)) {
  }
  else {
    nopush = 1;
    for (ccn = 0; ccn < dd->center_count; ccn++) {
      for (pn = 0; pn < dd->center[ccn].peak_count; pn++) {
        ddpn = dd->pn;
        for (offset = look_back (dd, ccn, pn); offset <= 0; offset++) {
          dd->pn = ddpn + offset;
          if (0 && justone_region) {
            int n;
            for (n = 0; n < level; n++) putchar (' ');
            printf ("level %d dd->pn %d, %d %d %d\n", level, dd->pn, ccn, pn, offset);
          }
          if (justone_region) {
            printf ("pushpop level %d: peak %d at %d, offset %d at %d, ccn %d peak %d\n", level,
                    ddpn, dd->p->peak[ddpn].loc,
                    dd->pn, dd->p->peak[dd->pn].loc, ccn, pn);
          }
          if (push_center (dd, ccn, pn + offset)) {
            if (justone_region) {
              int n;
              for (n = 0; n < level; n++) putchar (' ');
              printf ("level %d push ccn %d peak %d at %d succeeded\n",
                      level, ccn, pn + offset, dd->pn);
            }
            nopush = 0;
            if (decompose_region (dd) == 0 && dd->p->peak_count > 1 && dd->decomposition_count == 0)
              do {
                //              int last_ddpn = dd->pn;
                dd->pn = next_peak (dd);
                //              printf ("pushpop level %d: last peak %d, next peak %d\n", level, last_ddpn, dd->pn);
              } while (dd->pn != -1 && decompose_region (dd) == 0);
            if (justone_region) {
              int n;
              for (n = 0; n < level; n++) putchar (' ');
              printf ("level %d  pop ccn %d peak %d at %d\n",
                      level, ccn, pn + offset, ddpn + offset);
            }
            pop_center (dd);
            break;
          }
          else if (1) {
            if (justone_region) {
              int n;
              for (n = 0; n < level; n++) putchar (' ');
              printf ("level %d push ccn %d peak %d at %d + %d = %d failed\n",
                      level, ccn, pn, dd->pn, offset, pn + offset);
            }
          }
        }
        dd->pn = ddpn;
      }
    }
    if (nopush && 0 && justone_region && !reason) 
      {if (asprintf (&reason, "no match to peaks") == -1) exit (1);}
  }
  if (0 && justone_region) {
    int n;
    for (n = 0; n < level; n++) putchar (' ');
    printf ("return level %d: %d %d\n", level, dd->decomposition_count, decomp_count);
  }
  level--;
  if (dd->decomposition_count == decomp_count && do_scaled_distance && match_all_peaks (dd)) {
    if (justone_region)
      printf ("decomposition at all big peaks matched: %d components, scaled %d\n", 
              dd->component_count, do_scaled_distance);
    try_distance (dd, 1);
  }
  return dd->decomposition_count > decomp_count;
}

#else

static void
decompose_region (DecomposeData *dd)
{
  if (dd->p->peak_count == 0)
    try_distance (dd);
  else {
    int ccn;
    for (ccn = 0; ccn < dd->center_count; ccn++) {
      if (push_center (dd, ccn, 0)) {
        decompose_region (dd);
        pop_center (dd);
      }
    }
  }
}

#endif

typedef struct {SpikeData sd; float v[SNAPLEN];} SpikeDataVals;

static void
write_one_spikedata (SpikeData *sdp, /*@unique@*/ float *waveform, FILE *f)
{
  static SpikeDataVals *buf;
  static size_t bufidx;
  size_t bufcnt = 4096;

  if (sdp == 0) {
    fwrite (buf, sizeof *buf, bufidx, f) == bufidx || DIE;
    free (buf);
    buf = 0;
    return;
  }
  if (buf == 0)
    TMALLOC (buf, bufcnt);
  buf[bufidx].sd = *sdp;
  TMEMCPY (buf[bufidx].v, waveform, SNAPLEN);
  if (++bufidx == bufcnt) {
    fwrite (buf, sizeof *buf, bufcnt, f) == bufcnt || DIE;
    bufidx = 0;
  }
}

static void
read_all_spikedata (FILE *f, int extra)
{
  SpikeDataVals *buf;
  int file_size, n;
  size_t bufidx, bufcnt = 4096, rdcnt;

  fseek (f, 0, SEEK_END);
  file_size = ftell (f);
  rewind (f);
  spike_count = file_size / (int)sizeof *buf;
  spike_count * (int)sizeof *buf == file_size || DIE;
  spikedata_alloc = spike_count + extra;
  TMALLOC (spikedata, spikedata_alloc);
  TMALLOC (spike_waveform, spikedata_alloc);
  TMALLOC (buf, bufcnt);
  bufidx = bufcnt;
  for (n = 0; n < spike_count; n++) {
    if (bufidx == bufcnt) {
      rdcnt = fread (buf, sizeof *buf, bufcnt, f);
      rdcnt == bufcnt || ((int)rdcnt == spike_count - n && feof (f)) || DIE;
      bufidx = 0;
    }
    spikedata[n] = buf[bufidx].sd;
    TMALLOC (spike_waveform[n], SNAPLEN);
    TMEMCPY (spike_waveform[n], buf[bufidx].v, SNAPLEN);
    bufidx++;
  }
  memset (spikedata + spike_count, 0, extra * sizeof *spikedata);
  free (buf);
}

/*@dependent@*/static FILE *
write_all_spikedata (void)
{
  int n;
  FILE *f;

  (f = tmpfile ()) || DIE;
  for (n = 0; n < spike_count; n++) {
    write_one_spikedata (spikedata + n, spike_waveform[n], f);
    free (spike_waveform[n]);
  }
  free (spikedata);
  spikedata = 0;
  free (spike_waveform);
  spike_waveform = 0;
  spike_count = 0;
  spikedata_alloc = 0;
  return f;
}

static void
component_spikedata (DecomposeData *dd)
{
  double min;
  int min_n = 0, n, dn;
  Decomposition* decomp;
  Component *comp;
  static SpikeData sd;

  min = HUGE_VAL;
  for (n = 0; n < dd->decomposition_count; n++)
    if (dd->decomposition[n].residue < min) {
      min = dd->decomposition[n].residue;
      min_n = n;
    }
  decomp = dd->decomposition + min_n;
  if (justone_region || global_rn == debug_region)
    printf ("using decomposition %d\n", min_n);
  for (dn = 0; dn < decomp->component_count; dn++) {
    comp = decomp->component + dn;
    sd.partition = comp->ccn;
    sd.sample = dd->p->region[dd->rn].sample + comp->start + PRESAMPLES;

    //    printf ("%s line %d: %f %d %f %d\n", __FILE__, __LINE__, sd.sample, dd->p->region[dd->rn].sample, comp->start, PRESAMPLES);
    sd.distance = comp->distance;
    if (justone_region || global_rn == debug_region)
      printf ("component_spikedata: spike at %f, ccn %d, distance %f, scale: %f, from noise: %f\n",
              sd.sample, sd.partition, sd.distance, comp->scale,
              white_distance_from_zero (comp->raw, SNAPLEN));
    sd.leave_out = 0;
    sd.nocenter = 0;
    sd.scale = comp->scale;
    sd.overlapped = comp->overlapped;
    if (justone_region)
      printf ("spike at %f overlapped: %d\n", sd.sample, sd.overlapped);
    write_one_spikedata (&sd, comp->raw, spikedata_tmpfile);
    if (sd.partition == 0 && comp->raw[SNAPLEN-1] < -CLIPVAL) { /* debug */
      mfvec ("spike", comp->raw, SNAPLEN);
      printf ("\n");
      printf("%s line %d: sample %f distance %f", __FILE__, __LINE__, sd.sample, sd.distance);
      printf ("\n");
      note ("%s line %d: sample %f distance %f", __FILE__, __LINE__, sd.sample, sd.distance);
    }
  }
}

static void
free_region_data_partial (RegionData *p)
{
  int n;
  CSpline *csp;

  if (p->data_count)
    free (p->data);
  p->data = 0;
  p->data_count = 0;
  if (p->peak_count)
    free (p->peak);
  p->peak = 0;
  p->peak_count = 0;
  if (p->region_count) {
    if (p->regionm) {
      int start = 0, end = p->region_count;
#ifdef JUSTONE
      start = justone; end = justone + 1;
#endif
      for (n = start; n < end; n++) {
        gsl_interp_accel_free (p->regionm[n].acc);
        csp = (CSpline *)p->regionm[n].spline;
        if (csp->interp == 0 && csp->y_float == p->fdata)
          p->fdata = 0;
        ss_spline_free (p->regionm[n].spline);
      }
      free (p->regionm);
    }
  }
  free (p->peak_region);
  free (p->fdata);
  free (p->done);
}

/*@dependent@*/ static FILE *
write_unclassified_regions (RegionData *p)
{
  size_t bufcnt = (size_t)MAX (1024*1024, maxregion);
  short *buf;
  size_t space_left = bufcnt;
  size_t bufidx = 0;
  size_t data_left;
  int datidx, dn, rn, count;
  FILE *f;

  (f = tmpfile ()) || DIE;
  TMALLOC (buf, bufcnt);
  
  for (dn = rn = 0; rn < p->region_count; rn++) {
    if (p->regionm[rn].cluster == -1)
      p->region[rn].type = -1;
    if (p->region[rn].type == -1) {
      p->region[rn].count >= 0 || DIE;
      data_left = (size_t)p->region[rn].count;
      datidx = 0;
      while (data_left > 0) {
        count = (int)MIN (space_left, data_left);
        TMEMCPY (buf + bufidx, p->data + dn + datidx, count);
        space_left -= count;
        data_left -= count;
        bufidx += count;
        datidx += count;
        if (space_left == 0) {
          fwrite (buf, sizeof *buf, bufcnt, f) == bufcnt || DIE;
          space_left = bufcnt;
          bufidx = 0;
        }
      }
    }
    dn += p->region[rn].count;
  }
  if (bufidx > 0)
    fwrite (buf, sizeof *buf, bufidx, f) == bufidx || DIE;
  free (buf);
  return f;
}

static void
write_unclassified_region (FILE *f, int len, /*@unique@*/ short *data)
{
  size_t bufcnt = (size_t)MAX (1024*1024, maxregion);
  static size_t space_left, bufidx;
  static short *buf;

  if (len == 0) {
    fwrite (buf, sizeof *buf, bufidx, f) == bufidx || DIE;
    free (buf);
    buf = 0;
    return;
  }
  if (buf == 0) {
    TMALLOC (buf, bufcnt);
    space_left = bufcnt;
  }

  len >= 0 || DIE;
  if (space_left < (size_t)len) {
    TMEMCPY (buf + bufidx, data, space_left);
    fwrite (buf, sizeof *buf, bufcnt, f) == bufcnt || DIE;
    len - space_left < bufcnt || DIE;
    TMEMCPY (buf, data + space_left, len - space_left);
    bufidx = len - space_left;
    space_left = bufcnt - bufidx;
  }
  else {
    TMEMCPY (buf + bufidx, data, len);
    bufidx += len;
    space_left -= len;
  }
  still_unclassified++;
}

static short *
read_unclassified_region (FILE *f, int len)
{
  size_t bufcnt = (size_t)MAX (1024*1024, maxregion);
  size_t rdcnt;
  static size_t data_left, bufidx;
  static short *buf, *rbuf;
  short *data;

  if (len == -1) {
    free (buf);
    free (rbuf);
    buf = 0; rbuf = 0;
    return 0;
  }
  if (len == 0) {
    data_left = 0;
    rewind (f);
    return 0;
  }
  if (buf == 0) {
    TMALLOC (buf, bufcnt);
    TMALLOC (rbuf, maxregion);
  }

  if ((int)data_left < len) {
    feof (f) && DIE;
    TMEMCPY (rbuf, buf + bufidx, data_left);
    rdcnt = fread (buf, sizeof *buf, bufcnt, f);
    rdcnt == bufcnt || feof (f) || DIE;
    TMEMCPY (rbuf + data_left, buf, len - data_left);
    bufidx = len - data_left;
    data_left = rdcnt - bufidx;
    data = rbuf;
  }
  else {
    data = buf + bufidx;
    bufidx += len;
    data_left -= len;
  }
  return data;
}  

static int
flat_match_center (DecomposeData *dd, int ccn, double start, double *dp)
{
  int i;
  Center *cnp = &dd->center[ccn];
  float buf[SNAPLEN], ctr[SNAPLEN], dif[SNAPLEN];
  int istart = (int)ceil (start);
  double offset = istart - start;
  double sum, len, proj;
  double threshold = *dp;

  if (istart + SNAPLEN > dd->p->region[0].count)
    return 0;
  TMEMCPY (buf, dd->p->fdata + istart, SNAPLEN);
  for (i = 0; i < SNAPLEN; i++)
    ctr[i] = ss_cspline_eval (cnp->spline, i + offset, cnp->acc);
  for (i = 0; i < SNAPLEN; i++)
    dif[i] = buf[i] - (ctr[i] - mean);
  if ((*dp = white_distance_from_zero (dif, SNAPLEN)) > threshold)
    return 0;
  whiten_snap (ctr);
  sum = 0;
  for (i = 0; i < SNAPLEN; i++)
    sum += ctr[i] * ctr[i];
  len = sqrt (sum);
  for (i = 0; i < SNAPLEN; i++)
    ctr[i] /= len;
  proj = 0;
  whiten_snap (buf);
  for (i = 0; i < SNAPLEN; i++)
    proj += buf[i] * ctr[i];
  if (proj < 0)
    return 0;
  if (proj < len) {
    double noise_expect = gsl_sf_erf_Q (proj);
    double cluster_expect = gsl_sf_erf_Q (len - proj);
    if (cluster_expect < 1000 * noise_expect)
      return 0;
  }
  return 1;
}

static int
flat_try_center (DecomposeData *dd, int ccn, int pnc0, double *dp, Component **cpp)
{
  static Component cp;
  float v[SNAPLEN];
  int loose;
  int mrop = 0, fmc = 0;

  if (center_peak_polarity (dd, ccn, pnc0) != region_peak_polarity (dd)) {
    if (dd->p->peak[dd->pn].loc >= 27255524 && dd->p->peak[dd->pn].loc <= 27255550 && dd->flat == 2)
      printf ("center %d peak %d polarity differs from region peak %d\n", ccn, pnc0, dd->pn);
    return 0;
  }

  dd->seed.type = dd->center[ccn].type;
  dd->seed.peak = dd->center[ccn].peak;
  dd->seed.cdata = dd->center[ccn].raw;
  dd->seed.val = v;
  dd->seed.ccn = ccn;
  *dp == 11.5 || *dp == 20 || DIE;
  loose = *dp == 20;
  if ((cp.peak_count = region_center_peaks (dd->p, dd->pn, &dd->seed, pnc0, &cp.start, loose)) != 0
      && (mrop = matches_remaining_orig_peaks (dd, &cp))
      && (fmc = flat_match_center (dd, ccn, cp.start, dp))) {
    *cpp = &cp;
    dd->seed.val = 0;
    return 1;
  }
  dd->seed.val = 0;
  return 0;
}

static void
flat_add_center (DecomposeData *dd, Component *cp)
{
  int cpn;

  mark_peaks (dd, dd->seed.peak, dd->seed.type / 2 + 1, cp->start, dd->p);
  cpn = dd->component_count++;
  if (dd->component_count > dd->component_alloc) 
    TREALLOC (dd->component, dd->component_alloc = dd->component_count);
  cp->ccn = dd->seed.ccn;
  cp->scale = 1;
  cp->overlapped = 0;
  cp->p = dd->p;
  dd->component[cpn] = *cp;
}

static void
subtract_spikes (DecomposeData *dd)
{
  float *buf;
  int rcnt = dd->p->region[0].count;
  int n, cpn, buf_start, preexisting;
  Center *cnp;
  double start;
  RegionData *p;

  TMALLOC (buf, rcnt);
  TMEMCPY (buf, dd->p->fdata, rcnt);
  for (cpn = 0; cpn < dd->component_count; cpn++) {
    cnp = dd->center + dd->component[cpn].ccn;
    start = dd->component[cpn].start;
    for (n = 0; n < rcnt; n++)
      buf[n] -= ss_cspline_eval (cnp->spline, n - start, cnp->acc) - mean;
  }  
  p = new_region_data (buf, rcnt, dd->p->region[0].sample);
  TMALLOC (p->done, rcnt);
  TMEMCPY (p->done, dd->p->done, rcnt);
  find_mpk_spikes (p, buf, rcnt, dd->p->region[0].sample, buf_start = 0, preexisting = 2);
  region_check (p);
  fdata_spline_init (p);
  dd->p = p;
  dd->region_data_copy = p->fdata;
}

static DecomposeData *global_dd;

static int peak_cmp_dn (const int *a, const int *b)
{
  double amp_a, amp_b;
  global_dd->pn = *a; amp_a = region_peak_amplitude (global_dd);
  global_dd->pn = *b; amp_b = region_peak_amplitude (global_dd);
  return amp_a > amp_b ? -1 : (amp_a < amp_b); /* descending */
}

static int
flat_decompose_region_2 (DecomposeData *dd)
{
  int ccn, pn, match, ccn_at_best = 0, pn_at_best = 0, n;
  double d, best_d;
  Component *cp;
  static int peak_count;
  static int *peak_order;

  if (dd->p->peak_count > peak_count)
    TREALLOC (peak_order, peak_count = dd->p->peak_count);
  for (n = 0; n < dd->p->peak_count; n++)
    peak_order[n] = n;
  global_dd = dd;
  qsort (peak_order, (size_t)dd->p->peak_count, sizeof *peak_order, (int(*)(const void*,const void*))peak_cmp_dn);
  
  for (n = 0; n < dd->p->peak_count; n++) {
    dd->pn = peak_order[n];
    if (dd->p->done[dd->pn])
      continue;
    match = 0;
    best_d = DBL_MAX;
    for (ccn = 0; ccn < dd->center_count; ccn++)
      for (pn = 0; pn < dd->center[ccn].peak_count; pn++) {
        d = 11.5;
        if (flat_try_center (dd, ccn, pn, &d, &cp) && d < best_d) {
          best_d = d;
          ccn_at_best = ccn;
          pn_at_best = pn;
          match = 1;
        }
      }
    if (match) {
      d = 11.5;
      flat_try_center (dd, ccn_at_best, pn_at_best, &d, &cp);
      flat_add_center (dd, cp);
      if (justone_region || global_rn == debug_region)
        printf ("matched %d at %d, start %.1f\n", dd->pn, dd->p->peak[dd->pn].loc,
                dd->component[dd->component_count - 1].start);
    }
  }

  if (dd->component_count > 0) {
    subtract_spikes (dd);
    try_distance (dd, 1);
    dd->component_count = 0;
    return 1;
  }
  return 0;
}

static void
do_decompose_region (DecomposeData *dd)
{
  if (dd->p->peak_count == 0 || all_peaks_noise (dd))
    return;

#ifdef BIGFIRST
  dd->pn = max_peak (dd->p);
#else
  dd->pn = 0;
#endif
  drstart = (double)clock () / CLOCKS_PER_SEC;
  if (justone_region || global_rn == debug_region) {
    static char *name;
    static int passcnt;
    passcnt++;
    if (asprintf (&name, "orig%d", passcnt) == -1) exit (1);
    printf ("do_decompose_region at %d, %d samples - writing region data to vars.m as \"%s\"\n",
            dd->p->region[0].sample, dd->p->region[0].count, name);
    mfvec (name, dd->p->fdata, dd->p->region[0].count);
    free (name);
  }

  attempt_component_count = 0;
  min_residue_sum = 0;

  if (justone_region)
    free (reason), reason = 0;
  else if (reason == 0)
    {if (asprintf (&reason, "not used") == -1) exit (1);}

  decompose_region (dd);
  if (justone_region || global_rn == debug_region || dd->p->region[0].sample == 29714478)
    printf ("decompose_region found %d decompositions of region from %d %d\n",
            dd->decomposition_count,
            dd->p->region[0].sample,
            dd->p->region[0].sample + dd->p->region[0].count
            );
  if (1)
    if (dd->decomposition_count == 0) {
      if (justone_region || global_rn == debug_region || dd->p->region[0].sample == 29714478)
        printf ("trying scaled distance\n");
      do_scaled_distance = 1;
      drstart = (double)clock () / CLOCKS_PER_SEC;
      find_clean_peaks (dd);
      decompose_region (dd);
      if (justone_region || global_rn == debug_region || dd->p->region[0].sample == 29714478)
        printf ("decompose_region found %d scaled decompositions of region from %d %d\n",
                dd->decomposition_count,
                dd->p->region[0].sample,
                dd->p->region[0].sample + dd->p->region[0].count
                );
      do_scaled_distance = 0;
    }

  if (use_best_attempt_flag
      && dd->decomposition_count == 0
      && attempt_component_count > 0) {
    if (justone_region || global_rn == debug_region)
      printf ("trying best attempt\n");
    dd->component_count == 0 || DIE;
    use_best_attempt (dd);
    try_distance (dd, 1);
    attempt_to_octave (dd);
    dd->component_count = 0;
  }

  if (dd->decomposition_count)
    component_spikedata (dd);
}

static int uslevel;

static void
unclassified_spikedata (DecomposeData *dd, Peak *peak, int pcnt)
{
  static SpikeData sd;
  int start, n, max, n_at_max = 0;

  //  printf ("%d: %d %d\n", uslevel++, peak - dd->p->peak, pcnt);

  if (peak[pcnt-1].loc - peak[0].loc < SNAPLEN - PRESAMPLES) {
    sd.partition = -2;
    sd.sample = (double)peak[0].loc;
    sd.distance = -1;
    sd.leave_out = 0;
    sd.nocenter = 1;
    sd.scale = 0;
    sd.overlapped = 0;
    start = MAX (0, peak[0].loc - dd->p->region[0].sample - PRESAMPLES);
    if (start + SNAPLEN > dd->p->region[0].count)
      start = dd->p->region[0].count - SNAPLEN;
    if (dd->flag)
    {
      char msg[256];
      window ();
      clear ();
      draw_snap (dd->p->fdata + start, "black");
      snprintf (msg, 256, "unclassified at %d", dd->p->region[0].sample + start);
      draw_text (msg);
      if (show () == 'q')
        exit (0);
    }
    write_one_spikedata (&sd, dd->p->fdata + start, spikedata_tmpfile);
    if (justone_region) printf ("unclassified spike at %f\n", sd.sample);
  }
  else {
    max = 0;
    for (n = 1; n < pcnt; n++)
      if (peak[n].loc - peak[n-1].loc > max) {
        max = peak[n].loc - peak[n-1].loc;
        n_at_max = n;
      }
    unclassified_spikedata (dd, peak, n_at_max);
    unclassified_spikedata (dd, peak + n_at_max, pcnt - n_at_max);
  }
  uslevel--;
}

static int
flat_decompose_region_3 (DecomposeData *dd)
{
  int ccn, pn = 0;
  double d, min_d;
  Component cp_at_min, *cpp = 0;

  TMEMSET (dd->p->done, 0, dd->p->peak_count);
  dd->p0 = dd->p;
  for (dd->pn = 0; dd->pn < dd->p->peak_count; dd->pn++) {
    if (dd->p->done[dd->pn])
      continue;
    min_d = HUGE_VAL;
    for (ccn = 0; ccn < dd->center_count; ccn++) {
      for (pn = 0; pn < dd->center[ccn].peak_count; pn++) {
        d = 20;
        if (flat_try_center (dd, ccn, pn, &d, &cpp))
          break;
        if (dd->p->region[0].sample == SAMPLOC) {
          printf ("region peak %d, center %d peak %d: distance %f\n", dd->pn, ccn, pn, d);
        }
      }
      if (pn < dd->center[ccn].peak_count && d < min_d) {
        if (dd->p->region[0].sample == SAMPLOC)
          printf ("match d = %f\n", d);
        min_d = d;
        cp_at_min = *cpp;
        break;                  /* to take first match */
      }
    }
    if (min_d < HUGE_VAL)
      flat_add_center (dd, &cp_at_min);
      
  }
  dd->center_count >= 0 || DIE;
  if (dd->component_count > 0) {
    int cpn;
    if (!loose_matches)
      TCALLOC (loose_matches, (size_t)dd->center_count);
    for (cpn = 0; cpn < dd->component_count; cpn++) {
      int ccn = dd->component[cpn].ccn;
      (0 <= ccn && ccn < dd->center_count) || DIE;
      loose_matches[ccn]++;
    }
    subtract_spikes (dd);
    try_distance (dd, 1);
    dd->component_count = 0;
    return 1;
  }
  return 0;
}

#define CLIP_HWIDTH 4

static int
next_spike (float *buf, int end, int bufidx, AmpData *ampdata, double *xtp)
{
  int n, amp;
  Amp *amplist = ampdata->amplist;

  for (n = bufidx; n < end; n++) {
    amp = clip_pol (buf[n]);
    if (ampdata->ampcount == 0 || amp != amplist[4].amp) {
      memmove (amplist, amplist + 1, 4 * sizeof *amplist);
      amplist[4].loc = n;
      amplist[4].amp = amp;
      ampdata->ampcount++;
      if (ampdata->ampcount >= 5 && amplist[0].amp == 0 && amplist[2].amp == 0
          && amplist[4].amp == 0 && amplist[1].amp + amplist[3].amp == 0 && amplist[2].loc - amplist[1].loc < 40)
        ampdata->peakloc = (amplist[1].loc + amplist[4].loc) / 2, ampdata->start = 0;
      else if (ampdata->ampcount >= 5 && amplist[1].amp == 0 && amplist[4].amp == 0
               && amplist[2].amp + amplist[3].amp == 0 && amplist[3].loc - amplist[2].loc < 40)
        ampdata->peakloc = (amplist[2].loc + amplist[4].loc) / 2, ampdata->start = 1;
      else
        ampdata->peakloc = -SNAPLEN;
      //if it's not monotonic rail-to-rail, we don't want it:
      if (ampdata->peakloc > 0 && ampdata->start == 0) {
        int i = amplist[2].loc;
        int sign = SGN (buf[i + 1] - buf[i]);
        for (++i; i < amplist[3].loc; i++)
          if (SGN (buf[i + 1] - buf[i]) != sign)
            ampdata->peakloc = -SNAPLEN;
      }
    }
    if (n - ampdata->peakloc == 32 && n - amplist[ampdata->start].loc > SNAPLEN) {

      if (ampdata->start == 1) { /* amplist[3].loc and amplist[3].loc-1 at opposite rails */
        *xtp = amplist[3].loc - .5 - PRESAMPLES - CLIP_HWIDTH;
        if (floor (*xtp + .5) + SNAPLEN > end)
          return 0;
        return n + 1;
      }

      int plty, i = 0, type;
      double slope = 0, dx, dy, xt;

      plty = amplist[ampdata->start+1].amp;
      type = ampdata->start == 1 ? 0 : amplist[3].loc - amplist[2].loc;
      i = amplist[2].loc;
      if (type == 1) {
        if (plty * buf[i] > 0) slope = -plty * CLIPVAL - buf[i];
        else                   slope = -plty * CLIPVAL + buf[i];
      }
      else if (type < 6) {
        i = amplist[2].loc - 1;
        while ((buf[i] - mean) * (buf[i + 1] - mean) > 0)
          i++;
        slope = buf[i + 1] - buf[i];
      }
      else continue;
      dy = mean - buf[i];
      dx = dy / slope;
      xt = i + dx - PRESAMPLES - CLIP_HWIDTH;
      (i + dx > amplist[2].loc - 1 && i + dx < amplist[3].loc) || DIE;
      *xtp = xt;
      if (floor (*xtp + .5) + SNAPLEN > end)
        return 0;
      return n + 1;
    }
  }
  return 0;
}

static void
mark_clipped_peaks (DecomposeData *dd, AmpData *ampdata, double xt)
{
  Amp *amplist = ampdata->amplist;
  int pn, ss = dd->p0->region[0].sample, loc;
  Peak *pk = dd->p0->peak;
  float *data = dd->p->fdata;

  for (pn = 0; pn < dd->p0->peak_count; pn++) {
    loc = pk[pn].loc - ss;
    if (loc >= amplist[1].loc && loc < amplist[2].loc) {
      dd->p->done[pn] = 1;
      if (SGN (data[loc]) != amplist[1].amp)
        note ("peak at %d is clipped but has the wrong polarity", pk[pn].loc);
    }
    if (loc >= amplist[3].loc && loc < amplist[4].loc) {
      dd->p->done[pn] = 1;
      if (SGN (data[loc]) != amplist[3].amp)
        note ("peak at %d is clipped but has the wrong polarity", pk[pn].loc);
    }
  }
}

static void
clipped_spikedata (DecomposeData *dd)
{
  AmpData ampdata;
  int bufidx, start;
  double xt;
  static SpikeData sd;

  ampdata.peakloc = -SNAPLEN;
  ampdata.ampcount = 0;
  bufidx = 0;

  while ((bufidx = next_spike (dd->p->fdata, dd->p->region[0].count, bufidx, &ampdata, &xt)) != 0) {
    sd.partition = dd->center_count + (ampdata.amplist[1].amp > 0);
    if (sd.partition > dd->max_partition)
      dd->max_partition = sd.partition;
    start = floor (xt + .5);
    start >= 0 || DIE;
    start <= dd->p->region[0].count - SNAPLEN || DIE;
    sd.sample = start + PRESAMPLES + dd->p0->region[0].sample;
    //    printf ("spike at %15.6f\n", xt + dd->p0->region[0].sample);
    sd.distance = -1;
    sd.leave_out = 0;
    sd.nocenter = 1;
    sd.scale = 0;
    sd.overlapped = 0;
    mark_clipped_peaks (dd, &ampdata, xt);
    write_one_spikedata (&sd, dd->p->fdata + start, spikedata_tmpfile);
    if (justone_region) printf ("clipped spike at %f\n", sd.sample);
  }
}

static int
region_split_decompose (DecomposeData *dd)
{
  RegionData *fullp = dd->p;
  int segcnt;
  struct {int first, last;} *seg;
  int sn, pn, buf_start, preexisting, failcnt;

  TMALLOC (seg, fullp->peak_count);
  seg[sn=0].first = MAX (fullp->region[0].sample, fullp->peak[0].loc - BEFORE);
  for (pn = 1; pn < fullp->peak_count; pn++)
    if (fullp->peak[pn].loc - fullp->peak[pn-1].loc > BEFORE + AFTER + 1) {
      seg[sn++].last = fullp->peak[pn-1].loc + AFTER;
      seg[sn].first = fullp->peak[pn].loc - BEFORE;
    }
  seg[sn++].last = MIN (fullp->peak[pn-1].loc + AFTER, fullp->region[0].sample + fullp->region[0].count - 1);

  dd->p->region[0].count >= SNAPLEN || DIE;
  failcnt = 0;
  segcnt = sn;
  for (sn = 0; sn < segcnt; sn++) {
    RegionData *p = 0;
    float *buf;
    int count = seg[sn].last - seg[sn].first + 1;

    if (justone_region || global_rn == debug_region)
      printf ("segment %d of %d from %d to %d (%d)\n", sn, segcnt, seg[sn].first, seg[sn].last,
              dd->decomposition_count);

    count >= 0 || DIE;
    if (segcnt > 1) {
      TMALLOC (buf, count);
      TMEMCPY (buf, fullp->fdata + (seg[sn].first - fullp->region[0].sample), count);
      p = new_region_data (buf, count, seg[sn].first);
      TCALLOC (p->done, (size_t)count);
      find_mpk_spikes (p, buf, count, seg[sn].first, buf_start = 0, preexisting = 2);
      region_check (p);
      fdata_spline_init (p);
      dd->p0 = dd->p = p;
      dd->region_data_copy = p->fdata;
    }
    do_decompose_region (dd);
    if (dd->decomposition_count == 0) {
      dd->flat = 2;
      if (dd->p->region[0].sample == SAMPLOC) {
        int n;
        mfvec ("region", dd->p->fdata, dd->p->region[0].count);
        printf ("peak count: %d\n", dd->p->peak_count);
        for (n = 0; n < dd->p->peak_count; n++)
          printf ("%d: %d\n", n, dd->p->peak[n].loc - dd->p->region[0].sample + 1);
      }
      if (1)
        if (flat_decompose_region_3 (dd)) {
          dd->decomposition_count == 1 || DIE;
          component_spikedata (dd);
          dd->decomposition_count = 0;
        }
      dd->flat = 0;
      if (dd->p->peak_count > 0 && !all_peaks_noise (dd)) {
        int from = 0, to = 0;
        clipped_spikedata (dd);
        failcnt++;
        //      dd->flag = why_no_fit (dd);
        for (from = to = 0; from < dd->p0->peak_count; from++)
          if (!dd->p->done[from])
            dd->p0->peak[to++] = dd->p0->peak[from];
        dd->p0->peak_count = to;
        if (dd->p0->peak_count)
          unclassified_spikedata (dd, dd->p0->peak, dd->p0->peak_count);
      }
    }
    dd->decomposition_count = 0;

    if (dd->p != p && dd->p != fullp) {
      malloc_usable_size (dd->p) > 0 || DIE;
      free_region_data (dd->p);
    }
    if (segcnt > 1) {
      malloc_usable_size (p) > 0 || DIE;
      free_region_data (p);
    }
  }    
  dd->p = fullp;
  free (seg);
  return failcnt == 0;
}

static FILE *urfile2;

void
decompose_unclassified_regions (DecomposeData *dd)
{
  int rn, allnoise_flag, bigcnt = 0;
  int succeeded, failed, allnoise, nopeaks;
  time_t last_time = 0, now;
  double maxamp;
  int max_peak_count = 0;
  typedef struct {double min, sum, count, max;} Stats;
  static Stats *decomp_times;
  FILE *urfile1;
  int flat_spikes = 0;
  int flat_found, all_flat, splitcnt = 0, all_splits_decomposed = 0;

  fclose (fopen ("vars.m", "w"));

  nopeaks = allnoise = succeeded = failed = 0;
  dd->pn = 0;

  spikedata_tmpfile = write_all_spikedata ();

  urfile1 = write_unclassified_regions (dd->p);

  free_region_data_partial (dd->p);
  (void)read_unclassified_region (urfile1, 0);

  (urfile2 = tmpfile ()) || DIE;

  //  malloc_debug = 1;
#ifdef JUSTONE
  for (rn = justone; rn < dd->p->region_count; rn++) {
#else
  for (rn = 0; rn < dd->p->region_count; rn++) {
#endif
    global_rn = rn;
    allnoise_flag = 0;
    if (dd->p->region[rn].type == -1) {
      int buf_start, preexisting, count, i;
      float *buf;
      RegionData *p, *saved_p;
      RegionM rm;
      RegionM *rmp = &rm;
      short *y_short;
#ifdef JUSTONE
      printf ("%79s\rregion %d\n", "", rn);
#endif
      count = dd->p->region[rn].count;
      y_short = read_unclassified_region (urfile1, count);
      dd->region_data_ptr = y_short;
      TMEMSET (&rm, 0, 1);
      region_one_spline_init (dd->p->region + rn, rmp, y_short);
      
      TMALLOC (buf, count);
      for (i = 0; i < count; i++)
        buf[i] = ss_cspline_eval (rmp->spline, (double)i, rmp->acc);
      if (justone_region)
        mfvec ("buf", buf, count);
      p = new_region_data (buf, count, dd->p->region[rn].sample);
      dd->p0 = p;
      count >= 0 || DIE;
      TCALLOC (p->done, (size_t)count);
      find_mpk_spikes (p, buf, count, dd->p->region[rn].sample, buf_start = 0, preexisting = 1);
      if (justone_region) {
        int n;
        for (n = 0; n < p->peak_count; n++)
          printf ("peak at %d (%f)\n", p->peak[n].loc, (p->peak[n].loc)/2.5);
      }
      if (p->peak_count+1 > max_peak_count) {
        TREALLOC (decomp_times, p->peak_count+1);
        TMEMSET (decomp_times + max_peak_count, 0, (p->peak_count+1) - max_peak_count);
        max_peak_count = p->peak_count + 1;
      }
        
        
#ifdef JUSTONE
      {
        int n;
        printf ("%d peaks\n", p->peak_count);
        for (n = 0; n < p->peak_count; n++)
          printf ("  %9d - %9d = %3d\n", p->peak[n].loc, p->region[0].sample, p->peak[n].loc - p->region[0].sample);
        printf ("\n");
      }
#endif
      
      region_check (p);
      fdata_spline_init (p);
      dd->rn = 0;
      dd->decomposition_count = 0;
      saved_p = dd->p;
      dd->p = p;
      dd->region_data_copy = p->fdata;

      if (0)
      {
        double v[dd->p->region[0].count];
        int n;
        for (n = 0; n < dd->p->region[0].count; n++)
          v[n] = ss_cspline_eval (dd->p->regionm[0].spline, (double)n, dd->p->regionm[0].acc);
        mvec ("new", v, dd->p->region[0].count);
        exit (0);
      }

      {
        //      struct mallinfo m0, m1;
        //m0 = mallinfo ();
        //      mclr ();
        if (p->peak_count == 0)
          nopeaks++;
        all_flat = flat_found = 0;
        maxamp = 0;
        if (all_peaks_noise (dd)) {
          if (0) {
            if (p->peak_count > 0) {
              int n;
              printf ("region %d differs %d\n", rn, p->peak_count);
              printf ("count: %d, sample %d\n", p->region[0].count, p->region[0].sample);
              for (n = 0; n < p->peak_count; n++)
                printf ("peak %d at %d\n", n, p->peak[n].loc);
              exit (0);
            }
          }
          allnoise_flag = 1;
          allnoise++;
        }
        else {
          static int here_count;
          dd->flat = 1;
          if (justone_region || global_rn == debug_region)
            printf ("%d centers\n", dd->center_count);
          here_count++;
          if (flat_decompose_region_2 (dd)) {
            dd->decomposition_count == 1 || DIE;

            if (justone_region || global_rn == debug_region) {
              printf ("region %d: %d spikes by flat decomposition, %d samples\n", 
                      justone_region, dd->decomposition[0].component_count, p->region[0].count);
            }

            flat_spikes += dd->decomposition[0].component_count;
            component_spikedata (dd);
            dd->decomposition_count = 0;
            flat_found = 1;
          }
          if (justone_region) printf ("%s line %d: spike_count %d\n", __FILE__, __LINE__, spike_count);

          dd->flat = 0;

          all_splits_decomposed = 0;
          if (dd->p->peak_count > 0 && !all_peaks_noise (dd)) {
            all_splits_decomposed = region_split_decompose (dd);
            if (justone_region || global_rn == debug_region)
              printf ("remaining peaks were%s split decomposed\n", all_splits_decomposed ? "" : " not");
          }
          else {
            all_flat = flat_found;
            if (justone_region || global_rn == debug_region)
              printf ("no remaining non-noise peaks\n");
          }
        }

      }
      if (justone_region) printf ("%s line %d: spike_count %d\n", __FILE__, __LINE__, spike_count);
      
      if (dd->p != p) {
        malloc_usable_size (dd->p) > 0 || DIE;
        free_region_data (dd->p);
      }
          
      free_region_data (p);
      dd->p = saved_p;
      dd->rn = rn;
      if (all_splits_decomposed || all_flat) {
        dd->p->region[rn].type = -3;
        if (dd->decomposition_count)
          succeeded++;
      }
      else if (allnoise_flag)
        dd->p->region[rn].type = -2;
      else if (1)
        dd->p->region[rn].type = -4;
      else {
        write_unclassified_region (urfile2, count, y_short);
        if (!allnoise_flag) {
          if (maxamp > 16034) {
            bigcnt++;
            if (0)
              printf ("%79s\rbig %2d: region %5d, sample %9d, time %f\n",
                      "", bigcnt, rn, dd->p->region[rn].sample, dd->p->region[rn].sample / SF);
#ifdef JUSTONE //ifndef
            printf ("%79s\rdecompose region %d at %d failed\n", "", rn, dd->p->region[rn].sample);
            exit (0);
#endif
          }
          failed++;
        }
      }
      if (justone_region) printf ("%s line %d: spike_count %d\n", __FILE__, __LINE__, spike_count);
#ifdef JUSTONE
      exit (0);                 /* debug */
#endif
      if (0)
        if (justone_region) 
          exit (0);

      if ((now = time(0)) > last_time) {
        time_t projected_completion_time = start_time + (time_t)floor ((double)dd->p->region_count / (rn + 1) * (now - start_time) + .5);
        if (0) {
#ifndef __APPLE__
          struct mallinfo m;
          m = mallinfo ();
          printf ("arena: %d, mmaps: %d, mmspace: %d, used: %d, unused: %d in %d, top: %d\n",
                  m.arena, m.hblks, m.hblkhd, m.uordblks, m.fordblks, m.ordblks, m.keepcost);
          printf ("mtotal: %d\n", mtotal[0]);
#endif
        }

        if (!SHOWSUB) {
          char *t = ctime (&projected_completion_time);
          char * tstr =strdup (t);
          tstr[strlen (tstr) - 1] = 0;
          printf ("%157s\r  region %d of %d at %d: %d succeeded, %d failed, %d allnoise, %d splits, %d flat, %.1f so far, %.1f total, %s\r",
                  "", rn, dd->p->region_count, dd->p->region[rn].sample, succeeded, failed, allnoise, splitcnt, flat_spikes,
                  (double)(now - start_time) / 60, (double)(projected_completion_time - start_time) / 60, tstr);
          fflush (stdout);
          free (tstr);
        }

        fflush (stdout);
        last_time = now;
      }
      region_one_spline_free (rmp);
    }
  }
  if (justone_region) printf ("%s line %d: spike_count %d\n", __FILE__, __LINE__, spike_count);
  printf ("\nregion %d of %d: %d succeeded, %d failed, %d allnoise, %d nopeaks, %d big\n",
          rn, dd->p->region_count, succeeded, failed, allnoise, nopeaks, bigcnt);
  fclose (urfile1);
  write_unclassified_region (urfile2, 0, 0);
  write_one_spikedata (0, 0, spikedata_tmpfile);
  if (justone_region) printf ("%s line %d: spike_count %d\n", __FILE__, __LINE__, spike_count);
  read_all_spikedata (spikedata_tmpfile, still_unclassified);
  if (justone_region) printf ("%s line %d: spike_count %d\n", __FILE__, __LINE__, spike_count);
  fclose (spikedata_tmpfile);
  if (0)
    if (loose_matches) {
      int n;
      printf ("loose matches\n\ncenter count\n");
      for (n = 0; n < dd->center_count; n++)
        printf (" %2d: %5d\n", n, loose_matches[n]);
    }
}

void
read_spikedata (char *file_name, char *ext)
{
  FILE *f;
  int n;
  char buf[1];

  (f = fopen (change_filetype (file_name, ext ? ext : ".spd"), "rb")) || DIE;
  fread (&spike_count, sizeof spike_count, 1, f) == 1 || DIE;
  spike_count >= 0 || DIE;
  printf ("spd spike_count: %d\n", spike_count);
  if (!spikedata)
    TMALLOC (spikedata, spike_count);
  if (!spike_waveform)
    TMALLOC (spike_waveform, spike_count);
  fread (spikedata, sizeof *spikedata, (size_t)spike_count, f) == (size_t)spike_count || DIE;
  if (ext && strcmp (ext, ".spd1") != 0)
    for (n = 0; n < spike_count; n++) {
      TMALLOC (spike_waveform[n], SNAPLEN);
      fread (spike_waveform[n], sizeof **spike_waveform, SNAPLEN, f) == SNAPLEN || DIE;
    }
  fread (&buf, sizeof buf, 1, f) == 0 || DIE;
  feof (f) || DIE;
  fclose (f);
}

void
write_spikedata (char *file_name, int write_waveforms)
{
  FILE *f;
  int n;

  (f = fopen (file_name, "wb")) || DIE;
  spike_count >= 0 || DIE;
  fwrite (&spike_count, sizeof spike_count, 1, f) == 1 || DIE;
  fwrite (spikedata, sizeof *spikedata, (size_t)spike_count, f) == (size_t)spike_count || DIE;
  if (write_waveforms)
    for (n = 0; n < spike_count; n++)
      fwrite (spike_waveform[n], sizeof **spike_waveform, SNAPLEN, f) == SNAPLEN || DIE;
  fclose (f);
}

void
region_write_centers (DecomposeData *dd, char *file_name)
{
  FILE *f;
  char *cluster;
  int n;

  dd->center_count >= 0 || DIE;
  dd->p->region_count >= 0 || DIE;
  (f = fopen (change_filetype (file_name, ".ctr"), "wb")) || DIE;
  fwrite (&dd->center_count, sizeof dd->center_count, 1, f) == 1 || DIE;
  fwrite (dd->center, sizeof (Center), (size_t)dd->center_count, f) == (size_t)dd->center_count || DIE;
  TMALLOC (cluster, dd->p->region_count);
  for (n = 0; n < dd->p->region_count; n++)
    cluster[n] = dd->p->regionm[n].cluster;
  fwrite (cluster, sizeof cluster[0], (size_t)dd->p->region_count, f) == (size_t)dd->p->region_count || DIE;
  free (cluster);
  fclose (f);
  write_spikedata (change_filetype (file_name, ".spd"), 1);
}

SnapList *
region_read_centers (DecomposeData *dd, char *file_name)
{
  FILE *f;
  int ccn;
  SnapList *cc;

  TCALLOC (cc, 1);
  if ((f = fopen (change_filetype (file_name, ".ctr"), "rb")) == 0)
    return 0;
  if (fread (&dd->center_count, sizeof dd->center_count, 1, f) != 1)
    return 0;
  dd->center_count >= 0 || DIE;
  TMALLOC (dd->center, dd->center_count);
  if (fread (dd->center, sizeof (Center), (size_t)dd->center_count, f) != (size_t)dd->center_count) {
    free (dd->center);
    return 0;
  }
  dd->p->region_count >= 0 || DIE;
  for (ccn = 0; ccn < dd->center_count; ccn++) {
    Snap *s;
    dd->center[ccn].spline = 0;
    dd->center[ccn].acc = 0;
    if (0)
    {
#ifndef __APPLE__
      struct mallinfo m;
      m = mallinfo ();
      printf ("arena: %d, mmaps: %d, mmspace: %d, used: %d, unused: %d in %d, top: %d\n",
              m.arena, m.hblks, m.hblkhd, m.uordblks, m.fordblks, m.ordblks, m.keepcost);
#endif
    }
    center_spline_init_1 (&dd->center[ccn]);
    if (0)
    {
#ifndef __APPLE__
      struct mallinfo m;
      m = mallinfo ();
      printf ("arena: %d, mmaps: %d, mmspace: %d, used: %d, unused: %d in %d, top: %d\n",
              m.arena, m.hblks, m.hblkhd, m.uordblks, m.fordblks, m.ordblks, m.keepcost);
      printf ("\n");
#endif
    }
    TCALLOC (s, 1);
    TMALLOC (s->raw, SNAPLEN);
    TMALLOC (s->val, SNAPLEN);
    memcpy (s->raw, dd->center[ccn].raw, sizeof (SnapVals));
    memcpy (s->val, dd->center[ccn].raw, sizeof (SnapVals));
    whiten_snap (s->val);
    s->count = 1;
    add_snap (cc, s);
  }
  fclose (f);
  return cc;
}

void
region_unclassified (RegionData *p, int cluster)
{
  int rn, i, n, noise, unclassified, classified, decomposed, start;

  (void)read_unclassified_region (urfile2, 0);
  decomposed = classified = noise = unclassified = 0;
  for (rn = 0; rn < p->region_count; rn++) {
    if (p->region[rn].type == -3)
      decomposed++;
    else if (p->region[rn].type == -2)
      noise++;
    else if (p->region[rn].type == -4)
      unclassified++;
    else if (p->region[rn].type == -1) {
      exit (DIE);
      unclassified++;
      i = spike_count++;
      memset (&spikedata[i], 0, sizeof spikedata[i]);
//      spikedata[i].leave_out = 0;
//      spikedata[i].scale = 0;
      spikedata[i].partition = cluster;
      spikedata[i].sample = (double)p->region[rn].sample + BEFORE;
      //    printf ("%s line %d: %f\n", __FILE__, __LINE__, spikedata[i].sample);
      spikedata[i].distance = -1;
      TMALLOC (spike_waveform[i], SNAPLEN);
      start = BEFORE - PRESAMPLES;

      {
        RegionM rm;
        short *y_short;

        y_short = read_unclassified_region (urfile2, p->region[rn].count);
        TMEMSET (&rm, 0, 1);
        region_one_spline_init (p->region + rn, &rm, y_short);

        for (n = 0; n < SNAPLEN; n++)
          spike_waveform[i][n] = ss_cspline_eval (rm.spline, (double)start + n, rm.acc);
        region_one_spline_free (&rm);
      }
    }
    else if (p->region[rn].type < 0)
      exit (DIE);
    else classified++;
  }
  printf ("%d regions: %d classified, %d unclassified, %d decomposed, %d noise\n%d spikes\n",
          p->region_count, classified, unclassified, decomposed, noise, spike_count);
  (void)read_unclassified_region (urfile2, -1);
  fclose (urfile2);
}

void
revisit_assignments_0 (SnapList *cc)
{
  int n, ccn, ccn_at_min = 0, reassigned;
  double d, min;
  float w[SNAPLEN];
  static float zero[SNAPLEN];
  static int *count;

  TREALLOC (count, cc->count);
  TMEMSET (count, 0, cc->count);
  reassigned = 0;
  for (n = 0; n < spike_count; n++)
    if (spikedata[n].partition >= 0) {
      SZMEMCPY (w, spike_waveform[n]);
      whiten_snap (w);
      min = HUGE_VAL;
      for (ccn = 0; ccn < cc->count; ccn++) {
        d = distance_samples (spike_waveform[n], w, cc->snap[ccn]->raw, cc->snap[ccn]->val, SNAPLEN);
        if (ccn == 0) {
          if (d < 8)
            min = d;
          ccn_at_min = ccn;
        }
        else {
          if (d < min) {
            min = d;
            ccn_at_min = ccn;
          }
        }
      }
      min < HUGE_VAL || DIE;
      d = distance_samples (zero, w, zero, zero, SNAPLEN);
      if (d < min)
        min = d,
          ccn_at_min = -1;
      if (ccn_at_min != spikedata[n].partition) {
        reassigned++;
        spikedata[n].partition = ccn_at_min;
      }
      count[spikedata[n].partition]++;
    }
  for (ccn = 0; ccn < cc->count; ccn++)
    if (count[ccn] == 0)
      cc->snap[ccn]->noise = 1;
  printf ("%79s\r%d spikes reassigned\n", "", reassigned);
}

#define ROUND(x) floor ((x) * multiplier) / multiplier;

static void
region_count (double *v, int count, double b0, double region_size, int *lop, int *midp, int *hip)
{
  int n, lo, mid, hi;
  double b1 = b0 + region_size;
  double b2 = b1 + region_size;
  double b3 = b2 + region_size;

  lo = mid = hi = 0;
  for (n = 0; n < count; n++) {
    if (v[n] < b0)
      continue;
    if (v[n] > b3)
      break;
    if (v[n] < b1)
      lo++;
    else if (v[n] > b2)
      hi++;
    else
      mid++;
  }
  *lop = lo;
  *midp = mid;
  *hip = hi;
}

static double
get_multiplier (double *v, int count, double radius)
{
  double x;
  int ex;

  x = fabs (v[0]);
  if (fabs (v[count-1]) > x)
    x = fabs (v[count-1]);
  x += 2.1 * radius;
  (void)frexp (x, &ex);
  ex = 53 - ex;
  ex >= 0 || DIE;
  return (double)(1 << (unsigned)ex);
}

static double
set_1d_radius (double *v, int count, double *firstp, double *lastp)
{
  double b0, b0lo, b0hi, binc, b0end, first, last = 0;
  int lo, mid, hi, min;
  double region_size = 2;

  b0lo = v[0] - region_size;
  b0hi = v[count - 1] - 2 * region_size;
  binc = (b0hi - b0lo) / 100;
  b0end = b0hi + binc / 2;
  min = count / 20;
  //  printf ("count %d, min %d\n", count, min);
  for (first = HUGE_VAL, b0 = b0lo; b0 < b0end; b0 += binc) {
    region_count (v, count, b0, region_size, &lo, &mid, &hi);
    //    printf ("lo %d, mid %d, hi %d\n", lo, mid, hi);
    if (mid < lo && mid < hi && lo > min && hi > min) {
      if (first == HUGE_VAL)
        first = b0 + 1.5 * region_size;
      last = b0 + 1.5 * region_size;
    }
    //    printf ("%f: %d %4d\n", b0 + 1.5 * region_size, mid > lo || mid > hi, mid);
  }
  if (first == HUGE_VAL)
    return 0;
  //  printf ("line %d: first %f, last %f\n", __LINE__, first, last);
  *firstp = first;
  *lastp = last;
  return 1 + (last - first) / 2;
}

double
cluster_1d (double *v, int count, double radius)
{
  double pos, leading_edge, trailing_edge, multiplier, next_hi, next_lo;
  double sum, avg, movement, new_leading_edge, new_trailing_edge, valley_pos = 0;
  int lo, hi, n, dir, lastdir = 1, vlo, vhi = 0;
  int zlo, zmid, zhi, pcnt, peak_count, valley_count;
  struct {int count; int dir; double pos;} *p;
  double first, last;

  count >= 0 || DIE;
  TCALLOC (p, (size_t)count*2);
  qsort (v, (size_t)count, sizeof *v, (int(*)(const void*,const void*))dbl_cmp_up);
  radius = set_1d_radius (v, count, &first, &last);
  if (radius == 0) {
    printf ("%s line %d: radius == 0\n", __FILE__, __LINE__);
    return HUGE_VAL;
  }

  multiplier = get_multiplier (v, count, radius);
  for (n = 0; n < count; n++)
    v[n] = ROUND (v[n]);
  radius = ROUND (radius);
  lo = hi = 0;
  count > 0 || DIE;
  while (hi + 1 < count && v[hi+1] == v[hi])
    hi++;
  leading_edge = v[0];
  pos = leading_edge - radius;
  trailing_edge = pos - radius;
  pcnt = zlo = zmid = zhi = 0;
  printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
  while (1) {
    for (sum = 0, n = lo; n <= hi; n++)
      sum += v[n];
    avg = ROUND (sum / (hi - lo + 1));
    movement = avg - pos;
    new_leading_edge = leading_edge + movement;
    new_trailing_edge = trailing_edge + movement;
    if (avg > pos && ((hi + 1 < count && v[hi+1] <= new_leading_edge)
                      || v[lo] < new_trailing_edge))
      dir = 1;
    else if (avg < pos && ((lo > 0 && v[lo-1] >= new_trailing_edge)
                           || v[hi] > new_leading_edge))
      dir = -1;
    else dir = 0;

    if (dir == 0 && lastdir != 0) {
      if (pos < -1.3)
        zlo++;
      else if (pos > 1.3)
        zhi++;
      else zmid++;
    }
    lastdir = dir;
    p[pcnt].dir = dir;
    p[pcnt].count = hi - lo + 1;
    p[pcnt].pos = pos;
    pcnt++;

    //printf ("%9.6f %2d %3d\n", pos, dir, hi - lo + 1);
    /*  */
    if (hi + 1 < count)
      next_hi = v[hi+1] - leading_edge;
    else
      next_hi = -1;
    while (lo + 1 < count && v[lo+1] == v[lo])
      lo++;
    next_lo = v[lo] - trailing_edge + 1 / multiplier;

    if (next_hi == -1 && lo == hi)
      break;
    if (lo < hi && (next_hi == -1 || next_lo < next_hi)) {
      lo++;
      pos += next_lo;
    }
    else if (next_hi == -1 && lo > hi)
      exit (DIE);
    else if (next_hi == next_lo) {
      hi++;
      lo++;
      pos += next_hi;
    }
    else if (next_hi < next_lo || lo == hi) {
      hi++;
      pos += next_hi;
    }
    else {
      printf ("%d %d %f %f\n", lo, hi, next_lo, next_hi);
      exit (DIE);
    }
    while (hi + 1 < count && v[hi+1] == v[hi])
      hi++;
    leading_edge = pos + radius;
    trailing_edge = pos - radius;
  }
  printf ("%79s\r  line %d\r", "", __LINE__); fflush (stdout);
  //printf ("%.1f: %2d %2d %2d\n", radius, zlo, zmid, zhi);
  for (n = 0; n < pcnt; n++)
    sum += p[n].count;
  avg = sum / pcnt;
  lo = vlo = -1;

  for (n = 0; n < pcnt; n++) {
    if ((p[n].dir == 0 || (n > 0 && p[n].dir == 1 && p[n-1].dir == -1))
        && p[n].pos >= first && p[n].pos <= last
        ) {
      if (vlo == -1)
        vlo = n;
      vhi = n;
    }
  }
  valley_count = 1;
  valley_pos = (p[vlo].pos + p[vhi].pos)/2;

  if (0)
  for (peak_count = valley_count = n = 0; n < pcnt; n++) {
    if ((double)p[n].count > avg) {
      if (vlo != -1) {
        if (peak_count == 1) {
          valley_count == 0 || DIE;
          printf ("line %d: valley_pos: %f\n", __LINE__, valley_pos);
          valley_pos = (p[vlo].pos + p[vhi].pos)/2;
          valley_count++;
        }
        vlo = -1;
      }
      if (p[n].dir == 0) {
        if (lo == -1)
          lo = n;
        hi = n;
      }
    }
    else if ((double)p[n].count < avg) {
      if (lo != -1) {
        lo = -1;
        peak_count++;
      }
      if ((p[n].dir == 0 || (n > 0 && p[n].dir == 1 && p[n-1].dir == -1))) {
        if (vlo == -1)
          vlo = n;
        vhi = n;
      }
    }
  }
  free (p);
  printf ("line %d: valley_count: %d\n", __LINE__, valley_count);
  if (valley_count == 1)
    return valley_pos;
  else
    return HUGE_VAL;

  if (valley_count == 1 && peak_count == 2)
    return valley_pos;
  else if (valley_count == 0 && peak_count == 1)
    return HUGE_VAL;
  else {
    FILE *f;
    (f = fopen ("cluster_1d_error", "wb")) || DIE;
    fwrite (v, sizeof *v, (size_t)count, f) == (size_t)count || DIE;
    fclose (f);
    exit (DIE);
  }
  return 0;
}

static double
biggest_n_gap (double *v, int count, int n, int *ip)
{
  int i, i_at_max = 0;
  double max = 0, d;

  count >- 0 || DIE;
  qsort (v, (size_t)count, sizeof *v, (int(*)(const void*,const void*))dbl_cmp_up);
  for (i = 0; i < count - n; i++)
    if ((d = v[i+n] - v[i]) > max) {
      max = d;
      i_at_max = i;
    }
  *ip = i_at_max;
  return (v[i_at_max+n] + v[i_at_max]) / 2;
}

typedef struct {unsigned char first, second;} Two;

static Two *
closest_two (SnapList *cc)
{
  int n, ccn;
  double d, min, min2;
  float w[SNAPLEN];
  static float zero[SNAPLEN];
  Two *sd;

  cc->count <= UCHAR_MAX - 1 || DIE;
  TMALLOC (sd, spike_count);
  TMEMSET (sd, 255, spike_count);
  for (n = 0; n < spike_count; n++) {
    SZMEMCPY (w, spike_waveform[n]);
    whiten_snap (w);
    min2 = min = HUGE_VAL;
    for (ccn = 0; ccn <= cc->count; ccn++)
      if ((d = distance (w, ccn == cc->count ? zero : cc->snap[ccn]->val)) < min2 && d < 20) {
        if (d <= min) {
          min2 = min;
          min = d;
          sd[n].second = sd[n].first;
          sd[n].first = (unsigned char)ccn;
        }
        else {
          min2 = d;
          sd[n].second = (unsigned char)ccn;
        }
      }
    if (n == 87434)
      printf ("spike %d: %d and %d\n", n, sd[n].first, sd[n].second);
  }
  return sd;
}

static void
check_sd (Two *sd)
{
  static Two *sd0;
  int count, n;

  for (count = n = 0; n < spike_count; n++)
    if (spikedata[n].partition == -1)
      count++;
  printf ("check_sd: %d unclassified\n", count);

  if (sd0 == 0) {
    sd0 = sd;
    return;
  }
  for (count = n = 0; n < spike_count; n++)
    if (spikedata[n].partition != -1)
      if (sd0[n].first != sd[n].first || sd0[n].second != sd[n].second) {
        if (count < 10)
          printf ("spike %6d of %d: %d %d -> %d %d\n", n, spike_count, sd0[n].first, sd0[n].second, sd[n].first, sd[n].second);
        count++;
      }
  printf ("check_sd: %d changed\n", count);
}

int
revisit_assignments_1 (SnapList *cc, SnapList *noise)
{
  int an, bn, i, n, pcnt, ccn;
  float *aw, *bw;
  float v[SNAPLEN];
  float w[SNAPLEN];
  static float zero[SNAPLEN];
  double sum, len, splitval, proj_a, proj_b, min_noise;
  static double *proj;
  static double *proj_spike;
  static int *spike_idx;
  int change = 0, unclassified = 0;
  Two *sd;
  static int pass;

  pass++;

  printf ("revisit_assignments\n");

  TREALLOC (proj, spike_count + noise->count);
  TREALLOC (proj_spike, spike_count + noise->count);
  TREALLOC (spike_idx, spike_count + noise->count);


  /*
  n = 241393;
  printf ("before: spike %d, sample %f, distance %f, partition: %d\n",
          n, spikedata[n].sample, spikedata[n].distance, spikedata[n].partition);
  mfvec ("spike0", spike_waveform[n], SNAPLEN);
  */

  for (n = 0; n < spike_count; n++)
    if (spikedata[n].partition == 0 && spike_waveform[n][SNAPLEN-1] < -CLIPVAL) { /* debug */
      mfvec ("spike", spike_waveform[n], SNAPLEN);
      printf ("spike %d, sample %f, distance %f\n", n, spikedata[n].sample, spikedata[n].distance);
      exit (DIE);
    }
  sd = closest_two (cc);
  check_sd (sd);
  for (an = 0; an < cc->count; an++)
    for (bn = an + 1; bn <= cc->count; bn++) {
      aw = cc->snap[an]->val;
      if (bn == cc->count)
        bw = zero;
      else
        bw = cc->snap[bn]->val;
      for (i = 0; i < SNAPLEN; i++)
        v[i] = bw[i] - aw[i];
      sum = 0;
      for (i = 0; i < SNAPLEN; i++)
        sum += v[i] * v[i];
      len = sqrt (sum);
      for (i = 0; i < SNAPLEN; i++)
        v[i] /= len;
      pcnt = 0;

      proj_a = 0; for (i = 0; i < SNAPLEN; i++) proj_a += aw[i] * v[i];
      proj_b = 0; for (i = 0; i < SNAPLEN; i++) proj_b += bw[i] * v[i];

      for (n = 0; n < spike_count; n++) {
        if (((int)sd[n].first == an && (int)sd[n].second == bn) || ((int)sd[n].first == bn && (int)sd[n].second == an)) {
          //if ((ccn = spikedata[n].partition) == an || ccn == bn) {
          spike_idx[pcnt] = n;
          SZMEMCPY (w, spike_waveform[n]);
          whiten_snap (w);
          proj[pcnt] = 0;
          for (i = 0; i < SNAPLEN; i++)
            proj[pcnt] += w[i] * v[i];
          if (proj[pcnt] <= proj_a + 2.5) {
            if (n == 87434) printf ("%s line %d: spike %d set to %d from %d\n", __FILE__, __LINE__, n, an, spikedata[n].partition);
            spikedata[n].partition = an;
          }
          else if (proj[pcnt] >= proj_b - 2.5) {
            if (n == 87434) printf ("%s line %d: spike %d set to %d from %d\n",
                                    __FILE__, __LINE__, n, (bn == cc->count ? -1 : bn), spikedata[n].partition);
            spikedata[n].partition = (bn == cc->count ? -1 : bn);
          }
          else {
            if (n == 87434) printf ("%s line %d: spike %d set to proj # %d\n", __FILE__, __LINE__, n, pcnt);
            pcnt++;
          }
        }
      }
      min_noise = HUGE_VAL;
      if (bn == cc->count)
        for (n = 0; n < noise->count; n++) {
          spike_idx[pcnt] = -1;
          SZMEMCPY (w, noise->snap[n]->val);
          proj[pcnt] = 0;
          for (i = 0; i < SNAPLEN; i++)
            proj[pcnt] += w[i] * v[i];
          if (proj[pcnt] < min_noise)
            min_noise = proj[pcnt];
          if (proj[pcnt] < proj_b - 2.5)
            pcnt++;
        }

      TMEMCPY (proj_spike, proj, pcnt);
//      splitval = cluster_1d (proj, pcnt, 1.8);
      n = 0;
      n++; splitval = biggest_n_gap (proj, pcnt, n, &i);// printf ("%d %f %f %f\n", i, proj[i], splitval, proj[i+n]);
      if (min_noise < splitval)
        splitval = min_noise;
      printf ("centers %d and %d: %f %f %f (%d)\n", an, bn, proj_a, splitval, proj_b, pcnt);
      for (i = 0; i < pcnt; i++) {
        if ((n = spike_idx[i]) > -1) {
          int new;
          if (proj_spike[i] < splitval)
            new = an;
          else
            new = (bn == cc->count ? -1 : bn);
          if (n == 87434)
            printf ("spike %d: proj %f, new %d\n", n, proj_spike[i], new);
          if (new != spikedata[n].partition) {
            int gonna_be_minus_one = new >= 0 && new < cc->count && cc->snap[new]->noise;
            int is_minus_one = spikedata[n].partition == -1;
            unclassified += new == -1;
            if (!(gonna_be_minus_one && is_minus_one)) {
              if ((pass == 2 && change < 10))
                printf ("spike %d reassigned from %d to %d\n", n, spikedata[n].partition, new);
              
              change++;
            }
            if (n == 87434) printf ("%s line %d: spike %d set to %d from %d\n", __FILE__, __LINE__, n, new, spikedata[n].partition);
            spikedata[n].partition = new;
          }
        }         
      }
    }
  for (n = 0; n < spike_count; n++) {
    ccn = spikedata[n].partition;
    if (ccn >= 0 && ccn < cc->count && cc->snap[ccn]->noise) {
      if (n == 87434) printf ("%s line %d: spike %d set to %d from %d\n", __FILE__, __LINE__, n, -1, spikedata[n].partition);
      spikedata[n].partition = -1;
    }
  }

  if (0)
    for (n = 0; n < spike_count; n++)
      if (spikedata[n].partition == 0 && spike_waveform[n][SNAPLEN-1] < -CLIPVAL) { /* debug */
        mfvec ("spike", spike_waveform[n], SNAPLEN);
        printf ("spike %d, sample %f, distance %f\n", n, spikedata[n].sample, spikedata[n].distance);
        exit (DIE);
      }

  printf ("%d set to -1\n", unclassified);
  printf ("revisit_assignments done\n");
  //  free (sd); 
  printf ("%s line %d: spike %d set to %d\n", __FILE__, __LINE__, 87434, spikedata[87434].partition);
  return change;
}

typedef struct {int an, bn; float *v; double proj_a, proj_b;} Pair;

static Pair *
get_pairs (SnapList *cc)
{
  Pair *pair;
  int pairidx, an, bn, i;
  int paircnt = cc->count * (cc->count + 1) / 2;
  double proj_a, proj_b, sum, len;
  float *aw, *bw, *v;
  static float zero[SNAPLEN];

  TMALLOC (pair, paircnt);
  pairidx = 0;
  for (an = 0; an < cc->count; an++)
    for (bn = an + 1; bn <= cc->count; bn++) {
      pairidx < paircnt || DIE;
      TMALLOC (v, SNAPLEN);
      aw = cc->snap[an]->val;
      bw = (bn == cc->count) ? zero : cc->snap[bn]->val;
      for (i = 0; i < SNAPLEN; i++)
        v[i] = bw[i] - aw[i];
      sum = 0;
      for (i = 0; i < SNAPLEN; i++)
        sum += v[i] * v[i];
      len = sqrt (sum);
      for (i = 0; i < SNAPLEN; i++)
        v[i] /= len;
      proj_a = 0; for (i = 0; i < SNAPLEN; i++) proj_a += aw[i] * v[i];
      proj_b = 0; for (i = 0; i < SNAPLEN; i++) proj_b += bw[i] * v[i];
      pair[pairidx].an = an;
      pair[pairidx].bn = bn;
      pair[pairidx].v = v;
      pair[pairidx].proj_a = proj_a;
      pair[pairidx].proj_b = proj_b;
      pairidx++;
    }
  pairidx == paircnt || DIE;
  return pair;
}

typedef struct {double proj; int pair;} PairProj;

static PairProj *
get_pairproj (SnapList *cc, Pair *pair)
{
  int paircnt = cc->count * (cc->count + 1) / 2;
  static float zero[SNAPLEN];
  int n, ccn, ccn_at_min = 0, pairidx, i;
  float w[SNAPLEN];
  double min, d, proj, r, x, y;
  Pair *pr;
  PairProj *pairproj;

  TMALLOC (pairproj, spike_count);
  for (n = 0; n < spike_count; n++) {
    SZMEMCPY (w, spike_waveform[n]);
    whiten_snap (w);
    min = HUGE_VAL;
    for (ccn = 0; ccn <= cc->count; ccn++)
      if ((d = distance (w, ccn == cc->count ? zero : cc->snap[ccn]->val)) < min && d < 20) {
          min = d;
          ccn_at_min = ccn;
      }
    if (min == HUGE_VAL) {
      spikedata[n].partition = cc->count;
      pairproj[n].pair = -1;
      continue;
    }
    spikedata[n].partition = ccn_at_min;
    min = HUGE_VAL;
    for (pairidx = 0; pairidx < paircnt; pairidx++) {
      pr = pair + pairidx;
      if (ccn_at_min != pr->an && ccn_at_min != pr->bn)
        continue;
      proj = 0;
      for (i = 0; i < SNAPLEN; i++)
        proj += w[i] * pr->v[i];
      if (proj > pr->proj_a + 2.5 && proj < pr->proj_b - 2.5) {
        r = distance (w, cc->snap[pr->an]->val);
        x = proj - pr->proj_a;
        y = sqrt (r*r - x*x);
        if (y < min) {
          min = y;
          pairproj[n].proj = proj;
          pairproj[n].pair = pairidx;
        }
      }
    }
    if (min == HUGE_VAL) {
      spikedata[n].partition = ccn_at_min;
      pairproj[n].pair = -1;
      continue;
    }
  }
  return pairproj;
}

static PairProj *
get_noise_pairproj (SnapList *cc, Pair *pair, SnapList *noise)
{
  int paircnt = cc->count * (cc->count + 1) / 2;
  int n, pairidx, i;
  float *w;
  double min, proj, r, x, y;
  Pair *pr;
  PairProj *pairproj;

  TMALLOC (pairproj, noise->count);
  for (n = 0; n < noise->count; n++) {
    w = noise->snap[n]->val;
    min = HUGE_VAL;
    for (pairidx = 0; pairidx < paircnt; pairidx++) {
      pr = pair + pairidx;
      if (pr->bn != cc->count && pr->an != cc->count)
        continue;
      proj = 0;
      for (i = 0; i < SNAPLEN; i++)
        proj += w[i] * pr->v[i];
      if (proj > pr->proj_a + 2.5 && proj < pr->proj_b - 2.5) {
        r = distance (w, cc->snap[pr->an]->val);
        x = proj - pr->proj_a;
        y = sqrt (r*r - x*x);
        if (y < min) {
          min = y;
          pairproj[n].proj = proj;
          pairproj[n].pair = pairidx;
        }
      }
    }
    if (min == HUGE_VAL) {
      pairproj[n].pair = -1;
      continue;
    }
  }
  return pairproj;
}

static char **
get_pair_ok (SnapList *cc)
{
  int an, bn, cn, i;
  float *a, *b, *c;
  float v[SNAPLEN], w[SNAPLEN];
  static float zero[SNAPLEN];
  char **pair_ok;
  double lb, lc, bccos, cosine;
  double cosine_threshold = 0;

  cc->count >= 0 || DIE;
  TMALLOC (pair_ok, cc->count + 1);
  for (an = 0; an <= cc->count; an++) {
    a = an == cc->count ? zero : cc->snap[an]->val;
    TCALLOC (pair_ok[an], (size_t)cc->count + 1);
    for (bn = 0; bn <= cc->count; bn++) {
      if (bn == an)
        continue;
      b = bn == cc->count ? zero : cc->snap[bn]->val;
      for (i = 0; i < SNAPLEN; i++)
        v[i] = b[i] - a[i];
      lb = distance (v, zero);
      for (cn = 0; cn < cc->count; cn++) {
        if (cn == an || cn == bn)
          continue;
        c = cn == cc->count ? zero : cc->snap[cn]->val;
        for (i = 0; i < SNAPLEN; i++)
          w[i] = c[i] - a[i];
        lc = distance (w, zero);
        bccos = 0;
        for (i = 0; i < SNAPLEN; i++)
          bccos += v[i] * w[i];
        cosine = bccos / (lb * lc);
        if (cosine > cosine_threshold)
          break;
      }
      if (cn == cc->count)
        pair_ok[an][bn] = 1;
    }
  }
  return pair_ok;
}

static void
free_pair_ok (SnapList *cc, /*@only@*/ char **pair_ok)
{
  int n;
  for (n = 0; n <= cc->count; n++)
    free (pair_ok[n]);
  free (pair_ok);
}

int
revisit_assignments_2 (SnapList *cc, SnapList *noise)
{
  int paircnt = cc->count * (cc->count + 1) / 2;
  int pairidx, n, pcnt, i;
  double min_noise, splitval;
  Pair *pr;
  Pair *pair;
  PairProj *pairproj;
  PairProj *noise_pairproj;
  double *proj;
  double *proj_spike;
  int *spike_idx;
  int change = 0;
  static int pass;
  signed char *spike_ccn;
  char **pair_ok;
  int ok_count = 0;

  pass++;
  
  TMALLOC (spike_ccn, spike_count);
  for (n = 0; n < spike_count; n++)
    spike_ccn[n] = spikedata[n].partition;

  pair = get_pairs (cc);
  pairproj = get_pairproj (cc, pair);
  noise_pairproj = get_noise_pairproj (cc, pair, noise);

  TMALLOC (proj, spike_count + noise->count);
  TMALLOC (proj_spike, spike_count + noise->count);
  TMALLOC (spike_idx, spike_count + noise->count);

  pair_ok = get_pair_ok (cc);
  for (pairidx = 0; pairidx < paircnt; pairidx++) {
    pr = pair + pairidx;
    if (!pair_ok[pr->an][pr->bn] || !pair_ok[pr->bn][pr->an])
      continue;
    ok_count++;
    //    printf ("doing cluster separation %d to %d\n", pr->an, pr->bn);
    pcnt = 0;
    for (n = 0; n < spike_count; n++)
      if (pairproj[n].pair == pairidx) {
        spike_idx[pcnt] = n;
        proj[pcnt] = pairproj[n].proj;
        pcnt++;
      }
    min_noise = HUGE_VAL;
    for (n = 0; n < noise->count; n++)
      if (noise_pairproj[n].pair == pairidx) {
        spike_idx[pcnt] = -1;
        proj[pcnt] = noise_pairproj[n].proj;
        if (proj[pcnt] < min_noise)
          min_noise = proj[pcnt];
        pcnt++;
      }
    TMEMCPY (proj_spike, proj, pcnt);
    splitval = biggest_n_gap (proj, pcnt, 1, &i);
    if (min_noise < splitval)
      splitval = min_noise;
    for (i = 0; i < pcnt; i++) {
      if ((n = spike_idx[i]) > -1) {
        int new;
        if (proj_spike[i] < splitval)
          new = pr->an;
        else
          new = (pr->bn == cc->count ? -1 : pr->bn);
        if (new != spikedata[n].partition) {
          spikedata[n].partition = new;
          change++;
        }
      }   
    }
  }
  //  printf ("ok_count = %d of %d\n", ok_count, paircnt);
  free_pair_ok (cc, pair_ok);
  for (n = 0; n < spike_count; n++) {
    int ccn = spikedata[n].partition;
    if (ccn >= 0 && ccn < cc->count && cc->snap[ccn]->noise)
      spikedata[n].partition = -1;
  }

  change = 0;
  for (n = 0; n < spike_count; n++)
    if (spike_ccn[n] != spikedata[n].partition)
      change++;
  free (spike_ccn);

  printf ("revisit_assignments done\n");
  free (proj);
  free (proj_spike);
  free (spike_idx);
  free (pair);
  free (pairproj);
  free (noise_pairproj);
  return change;
}

static int pair3_4;

static Pair *
get_pairs_m1 (SnapList *cc)
{
  Pair *pair;
  int pairidx, an, bn, i;
  int paircnt = cc->count * (cc->count + 1) / 2;
  double proj_a, proj_b, sum, len;
  float *aw, *bw, *v;
  static float zero[SNAPLEN];

  TMALLOC (pair, paircnt);
  pairidx = 0;
  for (an = -1; an < cc->count; an++)
    for (bn = an + 1; bn < cc->count; bn++) {
      pairidx < paircnt || DIE;
      if (an == 3 && bn == 4)
        pair3_4 = pairidx;
      TMALLOC (v, SNAPLEN);
      aw = (an < 0) ? zero : cc->snap[an]->val;
      bw = cc->snap[bn]->val;
      for (i = 0; i < SNAPLEN; i++)
        v[i] = bw[i] - aw[i];
      sum = 0;
      for (i = 0; i < SNAPLEN; i++)
        sum += v[i] * v[i];
      len = sqrt (sum);
      for (i = 0; i < SNAPLEN; i++)
        v[i] /= len;
      proj_a = 0; for (i = 0; i < SNAPLEN; i++) proj_a += aw[i] * v[i];
      proj_b = 0; for (i = 0; i < SNAPLEN; i++) proj_b += bw[i] * v[i];
      pair[pairidx].an = an;
      pair[pairidx].bn = bn;
      pair[pairidx].v = v;
      pair[pairidx].proj_a = proj_a;
      pair[pairidx].proj_b = proj_b;
      pairidx++;
    }
  pairidx == paircnt || DIE;
  return pair;
}

static double*
gen_thresholds (SnapList *cc)
{
  static float zero[SNAPLEN];
  static double *threshold;
  int ccn;

  TREALLOC (threshold, cc->count);
  for (ccn = 0; ccn < cc->count; ccn++) {
    double m, d;
    d = distance (cc->snap[ccn]->val, zero);
    m = .5*d*(d > 33 ? 33 : d)/33;
    threshold[ccn] = sqrt (11.3*11.3 + m*m);
    //    printf ("ccn %2d: threshold %f\n", ccn, threshold[ccn]);
  }
  return threshold;
}

static int
check_distance (SnapList *cc, int n, int ccn_at_min)
{
  double scale, threshold, d;
  int ccn, i;
  float *craw, *mraw, v[SNAPLEN], diff[SNAPLEN];

  if (ccn_at_min < 0 || ccn_at_min >= cc->count)
    return ccn_at_min;
  if (spikedata[n].overlapped) {
    (scale = spikedata[n].scale) || DIE;
    ccn = spikedata[n].partition;
    craw = cc->snap[ccn]->raw;
    for (i = 0; i < SNAPLEN; i++)
      v[i] = mean + scale * (craw[i] - mean);
    mraw = cc->snap[ccn_at_min]->raw;
    d = white_distance_from_zero (mraw, SNAPLEN);
    threshold = .5*d*(d > 33 ? 33 : d)/33;
    for (i = 0; i < SNAPLEN; i++)
      diff[i] = v[i] - mraw[i] + mean;
    if (justone_region)
      printf ("spike at %f: ccn %d, ccn_at_min %d, threshold %f, distance %f\n",
              spikedata[n].sample, ccn, ccn_at_min, threshold, white_distance_from_zero (diff, SNAPLEN));
    if (white_distance_from_zero (diff, SNAPLEN) > threshold) {
      if (justone_region)
        printf ("check_distance no match: %f exceeds threshold of %f\n",
                white_distance_from_zero (diff, SNAPLEN), threshold);
      return -2;
    }
    if (justone_region) {
      static int count;
      char *name;
      if (asprintf (&name, "template%d", count) == -1) exit (1);
      mfvec (name, mraw, SNAPLEN);
      free (name);
      if (asprintf (&name, "match%d", count) == -1) exit (1);
      mfvec (name, v, SNAPLEN);
      free (name);
      count++;
    }
    
  }
  else if (justone_region)
    printf ("spike at %f: ccn %d, ccn_at_min %d not overlapped\n",
              spikedata[n].sample, spikedata[n].partition, ccn_at_min);

  if (justone_region)
    printf ("check distance match ccn %d\n", ccn_at_min);
  return ccn_at_min;
}

static void
assign_to_closest (SnapList *cc, double *threshold)
{
  int n, ccn, ccn_at_min = 0;
  static float zero[SNAPLEN];
  float w[SNAPLEN];
  double min, d;
  
  if (justone_region)
    for (ccn = 0; ccn < cc->count; ccn++)
      printf ("center %d is %f from the noise\n", ccn, distance (cc->snap[ccn]->val, zero));
  for (n = 0; n < spike_count; n++) {
    int close_count;
    SZMEMCPY (w, spike_waveform[n]);
    whiten_snap (w);
    min = distance (w, zero);
    ccn_at_min = -1;

    if (justone_region)
      printf ("spike at %f is %f from noise\n", spikedata[n].sample, min);
    close_count = 0;
    for (ccn = 0; ccn < cc->count; ccn++) {
      if ((d = distance (w, cc->snap[ccn]->val)) < min) {
          min = d;
          ccn_at_min = ccn;
      }
      close_count += d < threshold[ccn];
      if (justone_region) {
        if (d < threshold[ccn])
          printf ("spike at %f is %f from %d, within %f\n", spikedata[n].sample, d, ccn, threshold[ccn]);
        if (ccn == spikedata[n].partition)
          printf ("spike at %f: close to %d at distance %f\n", spikedata[n].sample, ccn, d);
      }
    }
    if (ccn_at_min != -1 && close_count == 0)
      min = HUGE_VAL;
    if (ccn_at_min == -1) {
      double min_ratio = 9;
      for (ccn = 0; ccn < cc->count; ccn++) {
        double scale, ratio;
        d = scaled_distance (cc->snap[ccn]->val, w, SNAPLEN, &scale);
        if (scale >= .75 && scale <= 2) {
          ratio = scale > 1 ? scale : 1 / scale;
          if (d < threshold[ccn] && (ratio < min_ratio || (ratio == min_ratio && d < min))) {
            min = d;
            ccn_at_min = ccn;
            min_ratio = ratio;
          }
        }
        if (justone_region) {
          if (d < threshold[ccn])
            printf ("spike at %f is %f from %d scaled %f, within %f\n",
                    spikedata[n].sample, d, ccn, scale, threshold[ccn]);
          if (ccn == spikedata[n].partition)
            printf ("spike at %f: close to %d at distance %f scaled %f\n", spikedata[n].sample, ccn, d, scale);
        }
      }
    }
    if (ccn_at_min == -1 && min >= 20) {
      if (justone_region)
        printf ("spike at %f matched no centers: ccn_at_min: %d, min: %f, threshold: 20\n",
                spikedata[n].sample, ccn_at_min, min);
      min = HUGE_VAL;
    }

    if (min == HUGE_VAL) {
      if (spikedata[n].partition < cc->count)
        spikedata[n].partition = -2;
    }
    else {
      spikedata[n].partition = check_distance (cc, n, ccn_at_min);
      spikedata[n].distance = min;
      if (justone_region)
        printf ("spike at %f: closest to %d at distance %f\n", spikedata[n].sample, ccn_at_min, min);
    }
  }
}

static double
overlap_index (DecomposeData *dd, int spikeidx0, int spikeidx1)
{
  int ccn0, ccn1, left_peak, left_edge, right_peak, right_edge;
  double d, margin;

  ccn0 = spikedata[spikeidx0].partition;
  ccn1 = spikedata[spikeidx1].partition;
  d = spikedata[spikeidx1].sample - spikedata[spikeidx0].sample;
  if (ccn0 < 0 || ccn1 < 0 || ccn0 >= dd->center_count || ccn1 >= dd->center_count) 
    return d <= SNAPLEN;
  left_peak = dd->center[ccn1].peak[0].loc;
  left_edge = SNAPLEN + left_peak - dd->center[ccn1].peak[0].lfar;
  right_peak = dd->center[ccn0].peak[dd->center[ccn0].peak_count - 1].loc;
  right_edge = right_peak + dd->center[ccn0].peak[dd->center[ccn0].peak_count - 1].rfar;
  margin = left_edge - right_edge;
  
  if (d <= SNAPLEN - margin)
    return 1;
  if (d >= SNAPLEN)
    return 0;
  return (SNAPLEN - d) / margin;
}


static int spikeloc_cmp_up (const int *a, const int *b)
{
  return (spikedata[*a].sample < spikedata[*b].sample
          ? -1
          : (spikedata[*a].sample > spikedata[*b].sample)); /* ascending */
}

static void
reject_loose_overlaps (DecomposeData *dd, double *threshold)
{
  int *sorted, n, i, j, ccn, chain_start = 0, bad_chain = 0;

  TMALLOC (sorted, spike_count);
  for (n = 0; n < spike_count; n++)
    sorted[n] = n;
  qsort (sorted, (size_t)spike_count, sizeof *sorted, (int(*)(const void*,const void*))spikeloc_cmp_up);
  for (i = 0; i < spike_count; i++) {
    double oi, oi_left, oi_right;
    n = sorted[i];
    if ((ccn = spikedata[n].partition) < 0 || ccn >= dd->center_count)
      spikedata[n].partition = -2;
    if (i == 0) oi_left = 0; else oi_left = overlap_index (dd, sorted[i - 1], n);
    if (oi_left != 0 && bad_chain)
      spikedata[n].partition = -2;
    if (spikedata[n].partition >= 0) {
      double new_threshold;
      if (i == spike_count - 1) oi_right = 0; else oi_right = overlap_index (dd, n, sorted[i + 1]);
      oi = MAX (oi_left, oi_right);
      new_threshold = oi * 10 + (1 - oi) * threshold[ccn];
      if (spikedata[n].distance > new_threshold)
        spikedata[n].partition = -2;
    }
    if (oi_left == 0) {
      chain_start = i;
      bad_chain = spikedata[n].partition < 0;
    }
    else if (spikedata[n].partition < 0) {
      bad_chain = 1;
      for (j = chain_start; j < i; j++)
        spikedata[sorted[j]].partition = -2;
    }
  }
}

static void
noise_correlation (SnapList *cc, SnapList *noise)
{
  int n, i, snapidx, ccn;
  float v[SNAPLEN], residue[SNAPLEN], *craw, underlying[SNAPLEN];
  double sum_of_corrs, corr, min, max, scale;

  return;
  if (!justone_region)
    return;
  for (n = 0; n < spike_count; n++) {
    memcpy (v, spike_waveform[n], sizeof v);
    whiten_snap (v);
    sum_of_corrs = 0;
    min = 1;
    max = -1;
    for (snapidx = 0; snapidx < noise->count; snapidx++) {
      corr = snap_corr (v, noise->snap[snapidx]->val);
      if (corr < min)
        min = corr;
      if (corr > max)
        max = corr;
      sum_of_corrs += fabs (corr);
      //      printf ("corr: %f\n", corr);
    }
    (scale = spikedata[n].scale) || DIE;
    ccn = spikedata[n].partition;
    craw = cc->snap[ccn]->raw;
    for (i = 0; i < SNAPLEN; i++)
      underlying[i] = (mean + scale * (craw[i] - mean));
    for (i = 0; i < SNAPLEN; i++)
      residue[i] = spike_waveform[n][i] - underlying[i];
    whiten_snap (residue);
    whiten_snap (underlying);
    corr = snap_corr (residue, underlying);

    printf ("spike %d, correlation with noise: min %f, avg %f, max %f; with residue %f\n",
            n, min, sum_of_corrs / noise->count, max, corr);
    
  }
}

STATIC int spiketime_cmp (const int *a, const int *b)
{
  return spikedata[*a].sample < spikedata[*b].sample 
    ? -1
    : spikedata[*a].sample > spikedata[*b].sample; /* ascending */
}

static void
no_doublets (void)
{
  int *s;
  int n;
  double half_ms = .5 * SF/1000;

  TMALLOC (s, spike_count);
  for (n = 0; n < spike_count; n++)
    s[n] = n;
  qsort (s, (size_t)spike_count, sizeof *s, (int(*)(const void*,const void*))spiketime_cmp);
  for (n = 1; n < spike_count; n++)
    if (spikedata[s[n]].sample - spikedata[s[n-1]].sample <= half_ms
        && spikedata[s[n]].partition == spikedata[s[n-1]].partition
        && spikedata[s[n]].partition != -1)
      spikedata[s[n]].partition = spikedata[s[n-1]].partition = -2;
  free (s);
}

int
revisit_assignments (DecomposeData *dd, SnapList *cc, SnapList *noise)
{
  int paircnt = cc->count * (cc->count + 1) / 2;
  int pairidx, n, pcnt, i, ccn;
  double max_noise, splitval;
  float w[SNAPLEN];
  Pair *pr;
  Pair *pair;
  static double *proj;
  static double *proj_spike;
  static int *spike_idx;
  struct {short ccn; short pair;} *adj;
  char *adjust;
  int change = 0;
  FILE *f = fopen ("reass", "w");
  double *threshold;

  noise_correlation (cc, noise);
  if (0) {
    FILE *f;
    int i;
    double sum, tmp;
    (f = fopen ("noise_distance", "w")) || DIE;
    for (n = 0; n < noise->count; n++) {
      sum = 0;
      for (i = 0; i < SNAPLEN; i++) {
        tmp = noise->snap[n]->val[i];
        sum += tmp * tmp;
      }
      fprintf (f, "%f\n", sqrt (sum));
    }
    fclose (f);
  }

  TMALLOC (proj, spike_count + noise->count + 2);
  TMALLOC (proj_spike, spike_count + noise->count + 2);
  TMALLOC (spike_idx, spike_count + noise->count + 2);
  TMALLOC (adjust, paircnt);
  TMEMSET (adjust, 1, paircnt);
  TMALLOC (adj, spike_count);
  for (n = 0; n < spike_count; n++)
    adj[n].ccn = -2;

  pair = get_pairs_m1 (cc);

  threshold = gen_thresholds (cc);
  assign_to_closest (cc, threshold);
  no_doublets ();
  return 0;
  if (0)
    reject_loose_overlaps (dd, threshold);

  for (pairidx = 0; pairidx < paircnt; pairidx++) {
    pr = pair + pairidx;
    pcnt = 0;
    for (n = 0; n < spike_count; n++)
      if ((ccn = spikedata[n].partition) == pr->an || ccn == pr->bn) {
        SZMEMCPY (w, spike_waveform[n]);
        whiten_snap (w);
        proj[pcnt] = 0;
        for (i = 0; i < SNAPLEN; i++)
          proj[pcnt] += w[i] * pr->v[i];
        spike_idx[pcnt] = n;
        if (proj[pcnt] > pr->proj_a + 2.5 && proj[pcnt] < pr->proj_b - 2.5)
          pcnt++;
      }
    max_noise = -HUGE_VAL;
    if (pr->an == -1)
      for (n = 0; n < noise->count; n++) {
        spike_idx[pcnt] = -1;
        SZMEMCPY (w, noise->snap[n]->val);
        proj[pcnt] = 0;
        for (i = 0; i < SNAPLEN; i++)
          proj[pcnt] += w[i] * pr->v[i];
        if (proj[pcnt] > max_noise)
          max_noise = proj[pcnt];
        if (proj[pcnt] > pr->proj_a + 2.5)
          pcnt++;
      }
    spike_idx[pcnt] = -1;
    proj[pcnt++] = pr->proj_a + 2.5;
    spike_idx[pcnt] = -1;
    proj[pcnt++] = pr->proj_b - 2.5;
    TMEMCPY (proj_spike, proj, pcnt);
    splitval = biggest_n_gap (proj, pcnt, 1, &i);
    if (max_noise > splitval)
        splitval = max_noise;
    for (i = 0; i < pcnt; i++) {
      if ((n = spike_idx[i]) > -1) {
        int new = proj_spike[i] < splitval ? pr->an : pr->bn;
        if (adj[n].ccn == -2 && new != spikedata[n].partition) {
          adj[n].ccn = new;
          adj[n].pair = pairidx;
        }
        else if (adj[n].ccn >= 0 && new != adj[n].ccn && new != spikedata[n].partition)
          adj[n].ccn = -3;
      }   
    }
  }
  for (pairidx = 0; pairidx < paircnt; pairidx++) {
    pr = pair + pairidx;
    pcnt = 0;
    for (n = 0; n < spike_count; n++)
      if ((ccn = spikedata[n].partition) == pr->an || ccn == pr->bn) {
        SZMEMCPY (w, spike_waveform[n]);
        whiten_snap (w);
        proj[pcnt] = 0;
        for (i = 0; i < SNAPLEN; i++)
          proj[pcnt] += w[i] * pr->v[i];
        spike_idx[pcnt] = n;
        if (proj[pcnt] > pr->proj_a + 2.5 && proj[pcnt] < pr->proj_b - 2.5)
          pcnt++;
      }
    max_noise = -HUGE_VAL;
    if (pr->an == -1)
      for (n = 0; n < noise->count; n++) {
        spike_idx[pcnt] = -1;
        SZMEMCPY (w, noise->snap[n]->val);
        proj[pcnt] = 0;
        for (i = 0; i < SNAPLEN; i++)
          proj[pcnt] += w[i] * pr->v[i];
        if (proj[pcnt] > max_noise)
          max_noise = proj[pcnt];
        if (proj[pcnt] > pr->proj_a + 2.5)
          pcnt++;
      }
    spike_idx[pcnt] = -1;
    proj[pcnt++] = pr->proj_a + 2.5;
    spike_idx[pcnt] = -1;
    proj[pcnt++] = pr->proj_b - 2.5;
    TMEMCPY (proj_spike, proj, pcnt);
    splitval = biggest_n_gap (proj, pcnt, 1, &i);
    if (max_noise > splitval)
        splitval = max_noise;
    for (i = 0; i < pcnt; i++) {
      if ((n = spike_idx[i]) > -1) {
        int new = proj_spike[i] < splitval ? pr->an : pr->bn;

        if (adj[n].ccn == -3) {
          if (spikedata[n].partition == 3 || spikedata[n].partition == 4)
            fprintf (f, "%2d %2d wants to reassign %6d from %2d to %2d\n", pr->an, pr->bn, n, spikedata[n].partition, new);
          adjust[pairidx] = 0;
        }
        else if (new != spikedata[n].partition)
          adj[n].ccn == -1 || adj[n].ccn == new || DIE;
        /*
        else if (new == spikedata[n].partition)
          adj[n].ccn == -2 || DIE;
        */
      }   
    }
  }

  {
    int count = 0;
    for (n = 0; n < paircnt; n++)
      if (adjust[n])
        count++;
    if (0) printf ("adjusting %d of %d pairs\n", count, paircnt);
  }
  for (n = 0; n < spike_count; n++)
    if ((adj[n].ccn >= 0 && adjust[adj[n].pair]) || adj[n].ccn == -1) {
      if (spikedata[n].partition != adj[n].ccn)
        change++;
      spikedata[n].partition = adj[n].ccn;
    }
  for (n = 0; n < spike_count; n++) {
    int ccn = spikedata[n].partition;
    if (ccn >= 0 && ccn < cc->count && cc->snap[ccn]->noise) {
      spikedata[n].partition = -1;
      change++;
    }
  }
  free (proj);
  free (proj_spike);
  free (spike_idx);
  free (pair);
  free (adjust);
  free (adj);
  //  printf ("revisit_assignments done\n");
  fclose (f);
  return change;
}
