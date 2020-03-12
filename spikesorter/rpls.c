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

//gcc -Wall --std=c99 -o rpulse rpulse.c -lgsl -lgslcblas
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <error.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <argp.h>
#include "config.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_cdf.h>

#define RC "200"                /* integration time constant in seconds */
#define MIN_I "110"             /* minimum E phase in ms */
#define MIN_E "75"              /* minimum I phase in ms */
#define EABS_THR "1.5"
#define IABS_THR "2"
#define IREL_THR "0.2"
#define MIN_I_RAW "250"         /* milliseconds */
#define MIN_E_RAW "100"         /* milliseconds */

static int interval = 60 * 25000;

#define SDSAMP 50

const char *argp_program_version = "rpulse, part of " PACKAGE_STRING;
const char *argp_program_bug_address = "<roconnor@health.usf.edu>";
static char args_doc[] = "FILE.edt";
static char doc[] = ("generate I and E pulses from FILE.chan and add them to FILE.edt"
                     "\vThresholds (eabs, iabs, irel) are specified as a multiple of the noise."
                     "  Above the iabs threshold is considered to be in I phase.  Below eabs is"
                     " normally E phase, but a prolonged increase by irel is an I phase."
                     "  The transition to I (E) is at the highest (lowest) slope of the integrated signal.");
struct arguments
{
  char *spikefile;
  char *chanfile;
  double tc;
  double eabs_thr, iabs_thr, irel_thr;
  int min_i_raw, min_e_raw;     /* in samples */
  int min_i, min_e;
};

static struct argp_option options[] = {
  {"time-constant",   't', "TC",  0, "Set integration time constant to TC ms (default " RC ")" },
  {"eabs",            'a', "THR", 0, "Set E threshold to THR (default " EABS_THR ")" },
  {"iabs",            'A', "THR", 0, "Set I absolute threshold to THR (default " IABS_THR ")" },
  {"irel",            'r', "THR", 0, "Set I relative threshold to THR (default " IREL_THR ")" },
  {"min-e-raw",       'E', "THR", 0, "Set min time below irel before E to THR ms (default " MIN_E_RAW ")" },
  {"min-i-raw",       'I', "THR", 0, "Set min time above irel before I to THR ms (default " MIN_I_RAW ")" },
  {"min_i",           'i', "THR", 0, "Set min I phase to THR ms (default " MIN_I ")" },
  {"min_e",           'e', "THR", 0, "Set min E phase to THR ms (default " MIN_E ")" },
  { 0 }
};

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 't':
      arguments->tc = atof (arg);
      break;
    case 'a':
      arguments->eabs_thr = atof (arg);
      break;
    case 'A':
      arguments->iabs_thr = atof (arg);
      break;
    case 'r':
      arguments->irel_thr = atof (arg);
      break;
    case 'E':
      arguments->min_e_raw = atof (arg) * 25;
      break;
    case 'I':
      arguments->min_i_raw = atof (arg) * 25;
      break;
    case 'e':
      arguments->min_e = atof (arg);
      break;
    case 'i':
      arguments->min_i = atof (arg);
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num >= 2)
        argp_usage (state);
      if (state->arg_num == 0) {
        char *p = malloc (strlen (arg) + 2);
        strcpy (p, arg);
        char *ext = strrchr (p, '.');
        if (ext == NULL || strcmp (ext, ".edt") != 0)
          argp_usage (state);
        strcpy (ext, ".chan");
        arguments->spikefile = arg;
        arguments->chanfile = p;
      }
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 1)
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

static char *
change_filetype_trim (const char *file_name, char *suffix, int trim)
{
  char *dot, *p;
  static char *new_file;

  new_file = realloc (new_file, strlen (file_name) + strlen (suffix) + 1);
  strcpy (new_file, file_name);
  if ((dot = strrchr (new_file, '.')) == 0)
    strcat (new_file, suffix);
  else
    strcpy (dot, suffix);
  if (trim
      && ((p = strrchr (new_file, '/')) != 0
	  || (p = strrchr (new_file, '\\')) != 0))
    return p + 1;
  else
    return new_file;
}

static void
plot (int *hist, int cnt, double mean, double sd, int N)
{
  FILE *f = fopen ("/tmp/mupplot", "w");
  for (int i = 0; i < cnt; i++)
    fprintf (f, "%d\n", hist[i]);
  fclose (f);
  f = popen ("gnuplot -persist", "w");
  if (mean > 0)
    fprintf (f, "plot \"/tmp/mupplot\" w his, %d/sqrt(2*pi*%g**2)*exp(-(x-%g)**2/(2*%g**2))\n", N, sd, mean, sd);
  else
    fprintf (f, "plot \"/tmp/mupplot\" w his\n");
  fclose (f);
  getchar ();
}

static int
sumsq (int *hist, int bincnt, double mean, double sd, int N)
{
  int minbin = nearbyint (mean - sd);
  int maxbin = nearbyint (mean + sd);
  double sum = 0;
  double last = gsl_cdf_gaussian_P (-mean - .5 + minbin, sd);
  for (int i = minbin; i <= maxbin; i++) {
    double p = gsl_cdf_gaussian_P (-mean + .5 + i, sd);
    sum += pow ((p - last) * N - hist[i], 2);
    last = p;
  }
  return sum;
}

typedef struct {int *hist; int bincnt;} Hist;

double
sumsq_wrap (const gsl_vector *v, void *params)
{
  Hist *p = (Hist *)params;

  double mean = gsl_vector_get(v, 0);
  double sd   = gsl_vector_get(v, 1);
  int    N    = gsl_vector_get(v, 2);

  return sumsq (p->hist, p->bincnt, mean, sd, N);
}

void
normfitls (int *hist, int bincnt, double *mean, double *sd, int *N)
{
  Hist h = {hist, bincnt};

  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc (3);
  gsl_vector_set (x, 0, *mean);
  gsl_vector_set (x, 1, *sd);
  gsl_vector_set (x, 2, *N / 2);

  /* Initial step sizes */
  ss = gsl_vector_alloc (3);
  gsl_vector_set (ss, 0, *mean / 10);
  gsl_vector_set (ss, 1, *sd / 10);
  gsl_vector_set (ss, 2, *N / 10);

  /* Initialize method and iterate */
  minex_func.n = 3;
  minex_func.f = sumsq_wrap;
  minex_func.params = &h;

  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (gsl_multimin_fminimizer_nmsimplex2, 3);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  size_t iter = 0, iter_limit = 200;
  int last_val = 0, last_siz = 0;
  int valcnt = 0, sizcnt = 0;
  do {
    status = gsl_multimin_fminimizer_iterate(s);

    if (status) {
      error_at_line (0, 0, __FILE__, __LINE__, "minimizer error %d", status);
      goto OUT;
    }

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1e-2);
    if (s->fval == last_val) valcnt++; else valcnt = 1, last_val = s->fval;
    int siz = nearbyint (size * 1000);
    if (   siz == last_siz) sizcnt++; else sizcnt = 1, last_siz =    siz;
    if (valcnt >= 4 && sizcnt >= 5) {
      status = GSL_SUCCESS;
      break;
    }
  }
  while (status == GSL_CONTINUE && ++iter < iter_limit);
  if (iter >= iter_limit) {
    error_at_line (0, 0, __FILE__, __LINE__, "minimizer hit limit");
    plot (hist, bincnt, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_vector_get (s->x, 2));
  }
  if (status == GSL_SUCCESS) {
    *mean = gsl_vector_get (s->x, 0);
    *sd   = gsl_vector_get (s->x, 1);
    *N    = gsl_vector_get (s->x, 2);
  }
 OUT:
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);
}

static bool debug = false;

static double
modal_sd (short *buf, int count, double mean)
{
  int max_hidx = 0;
  int hist[65536];
  
  memset (hist, 0, sizeof hist);

  for (short *bufp = buf; bufp + SDSAMP <= buf + count; bufp += SDSAMP) {
    double sd;
    int hidx;
    long long sqsum = 0;
    for (int n = 0; n < SDSAMP; n++)
      sqsum += pow (bufp[n] - mean, 2);
    sd = sqrt ((double)sqsum / SDSAMP);
    hidx = nearbyint (sd);
    if (hidx > max_hidx)
      max_hidx = hidx;
    hist[hidx]++;
  }
  {
    int max = 0, n_at_max = 0;
    for (int n = 1; n < max_hidx; n++)
      if (hist[n] > max)
	max = hist[n_at_max = n];
    double mean = n_at_max;
    double sd = n_at_max / 8;
    int N = nearbyint (count / (60 * 25000.) * 30000);
    if (debug) plot (hist, max_hidx, mean, sd, N);
    normfitls (hist, max_hidx, &mean, &sd, &N);
    if (debug) fprintf (stderr, "n_at_max: %d, mean: %g\n", n_at_max, mean);
    if (debug) plot (hist, max_hidx, mean, sd, N);
    return mean;
  }
}

static double *noise_sd;
static int interval_count;

static void
rectify (short *x, int count)
{
  interval_count = count / interval;
  noise_sd = malloc (interval_count * sizeof *noise_sd);
  int iidx = 0;
  int thiscount = 0;
  for (int start = 0; start < count; start += thiscount) {
    int end = start + interval;
    if (end + interval > count)
      end = count;
    thiscount = end - start;
    double sum = 0;
    for (int i = start; i < end; i++)
      sum += x[i];
    assert (iidx < interval_count);
    noise_sd[iidx++] = modal_sd (x + start, thiscount, sum / thiscount);
    short mean = nearbyint (sum / thiscount);
    for (int i = start; i < end; i++)
      x[i] =  abs (x[i] - mean);
  }
  double sum = 0, min = 65536, max = 0;
  for (int iidx = 0; iidx < interval_count; iidx++) {
    sum += noise_sd[iidx];
    if (noise_sd[iidx] < min)
      min = noise_sd[iidx];
    if (noise_sd[iidx] > max)
      max = noise_sd[iidx];
  }
}

static inline double
get_noise_sd (int i)
{
  int iidx = i / interval;
  if (iidx >= interval_count)
    iidx = interval_count - 1;
  double nsd = 0;
  int mid = iidx * interval + interval / 2;
  if (i <= interval / 2)
    nsd = noise_sd[0];
  else if (i >= interval_count * interval - interval / 2)
    nsd = noise_sd[interval_count - 1];
  else if (i == mid)
    nsd = noise_sd[iidx];
  else {
    int left, lidx;
    if (i < mid) {
      left = mid - interval;
      lidx = iidx - 1;
    } else {
      left = mid;
      lidx = iidx;
    }
    assert (lidx >= 0 && lidx + 1 < interval_count && i >= left && i - left <= interval);
    nsd = noise_sd[lidx] + (double) (i - left) / interval * (noise_sd[lidx + 1] - noise_sd[lidx]);
  }
  return nsd;
}

int
main (int argc, char **argv)
{
  struct arguments arg;
  arg.tc = atof (RC);
  arg.eabs_thr = atof (EABS_THR);
  arg.iabs_thr = atof (IABS_THR);
  arg.irel_thr = atof (IREL_THR);
  arg.min_i_raw = atof (MIN_I_RAW) * 25;
  arg.min_e_raw = atof (MIN_E_RAW) * 25;
  arg.min_i = atoi (MIN_I);
  arg.min_e = atoi (MIN_E);

  argp_parse (&argp, argc, argv, 0, 0, &arg);

  FILE *chan = fopen (arg.chanfile, "r");
  FILE *f[2], *outfile;
  if (chan == NULL)
    error (1, errno, "fopen %s", arg.chanfile);
  fseek (chan, 0, SEEK_END);
  int bytes = ftell (chan);
  rewind (chan);
  short *x = malloc (bytes);
  if (fread (x, 1, bytes, chan) != bytes)
    error (1, errno, "fread %s", arg.chanfile);
  fclose (chan);

  if ((f[0] = fopen (arg.spikefile, "r")) == NULL)
    error_at_line (1, errno, __FILE__, __LINE__, "Can't open %s for read, aborting", arg.spikefile);
  if ((f[1] = tmpfile ()) == NULL)
    error_at_line (1, errno, __FILE__, __LINE__, "Can't open temp file, aborting");
  
  int count = bytes / sizeof *x;
  rectify (x, count);
  double dt = 1. / 25000;
  double rc = arg.tc / 1000;    /* tc is in ms, rc is in seconds */
  long double alpha = dt / (rc + dt);
  long double s = 0;
#define CP (1 << 11)            /* checkpoint interval */
  long double sn[count / CP + 1];
  for (int i = count - 1; i >= 0; i--) {
    s = alpha * x[i] + (1 - alpha) * s;
    if (i % CP == 0)
      sn[i / CP] = s;
  }
  double left = 0, right = s;
  enum State {ST_START, ST_I, ST_E, ST_IE, ST_EI} state = ST_START;

  int min_slope_at = 0, min_slope = 0;
  int max_slope_at = 0, max_slope = 0;
  int last_time = -fmax (arg.min_i * 25, arg.min_e * 25);
  double Emin = DBL_MAX, Imax = 0;
  int Eraw = 0, Iraw = 0;
  int itime = 0;
  int last_tick = 0;
  bool large_I = false;
  for (int i = 0; i < count; i++) {
    if (i % CP == 0)
      right = sn[i / CP];

    double nsd = get_noise_sd (i);
    double Eabs = nsd * arg.eabs_thr;
    double Iabs = nsd * arg.iabs_thr;
    double Irel = nsd * arg.irel_thr;
    
    if (left < Emin && right < Emin) {
      Emin = fmax (left, right);
      max_slope = right - left;
      max_slope_at = i;
    }
    else if (left > Imax && right > Imax) {
      Imax = fmin (left, right);
      min_slope = right - left;
      min_slope_at = i;
    }

    if (left > Iabs && right > Iabs) {
      min_slope = right - left;
      min_slope_at = i;
    }

    if (state == ST_E) {
      if (left > Emin + Irel && right > Emin + Irel) Iraw++; else Iraw = 0;
      if (Iraw >= arg.min_i_raw || (large_I = left > Iabs && right > Iabs)) {
        if (max_slope_at - last_time > arg.min_e * 25) {
          itime = (int)nearbyint (max_slope_at / 2.5);
          last_time = max_slope_at;
        }
        state = ST_I;
        Emin = Iabs;
      }
    }
    else if (state == ST_I) {
      if (left < Imax - Irel && right < Imax - Irel) Eraw++; else Eraw = 0;
      if (!large_I && left > Iabs && right > Iabs) {
        large_I = true;
        itime = (int)nearbyint (max_slope_at / 2.5);
        last_time = max_slope_at;
      }
      if (left < Eabs && right < Eabs && (Imax > Iabs || Eraw >= arg.min_e_raw)) {
        if (min_slope_at - last_time > arg.min_i * 25 && itime) {
          if (itime >= last_tick) {fprintf (f[1], "%5d%10d\n", 97, itime); last_tick = itime;}
          int tick = (int)nearbyint (min_slope_at / 2.5);
          if (tick  >= last_tick) {fprintf (f[1], "%5d%10d\n", 98,  tick); last_tick =  tick;}
        }
        last_time = min_slope_at;
        itime = 0;
        state = ST_E;
        Imax = 0;
      }
    }
    else if (left < Eabs && right < Eabs) {
      state = ST_E;
    }
    else if (left > Iabs && right > Iabs) {
      state = ST_I;
    }

    int diff = right - left;
    if (diff > max_slope) {
      max_slope = diff;
      max_slope_at = i;
    }
    if (diff < min_slope) {
      min_slope = diff;
      min_slope_at = i;
    }

    left = alpha * x[i] + (1 - alpha) * left;
    right = (right - alpha * x[i]) / (1 - alpha);
  }
  if (itime && itime >= last_tick)
    fprintf (f[1], "%5d%10d\n", 97, itime);
      
  char *outname = change_filetype_trim (arg.spikefile, ".tmp", 0);

  if ((outfile = fopen (outname, "wb")) == NULL)
    error_at_line (1, errno, __FILE__, __LINE__, "Can't open %s for write, aborting", outname);

  rewind (f[1]);
  int valcode, time;
  if(fscanf (f[0], "%5d%10d", &valcode, &time));
  if (valcode == 0) 
    if (fscanf (f[0], "%5d%10d", &valcode, &time));
  if (fscanf (f[0], "%5d%10d", &valcode, &time));
  fprintf (outfile, "   33   3333333\n");
  fprintf (outfile, "   33   3333333\n");

  {
    struct {int val; int time;} look[2];
    int donecnt = 0, fidx;

    for (fidx = 0; fidx < 2; fidx++)
      if (fscanf (f[fidx], "%5d%10d", &look[fidx].val, &look[fidx].time) != 2) {
	look[fidx].time = INT_MAX;
	donecnt++;
      }
    while (donecnt < 2) {
      fidx = look[0].time < look[1].time ? 0 : 1;
      fprintf (outfile, "%5d%10d\n", look[fidx].val, look[fidx].time);
      if (fscanf (f[fidx] , "%5d%10d", &look[fidx].val, &look[fidx].time) != 2) {
	look[fidx].time = INT_MAX;
	donecnt++;
      }
    }
  }

  fclose (f[0]);
  fclose (f[1]) ;
  fclose (outfile);
  remove (arg.spikefile);
  rename (outname, arg.spikefile);

  return 0;
} 
