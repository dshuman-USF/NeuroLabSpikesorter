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

//gcc -Wall --std=c99 -O2 -o mup mup.c -lm -lgsl -lgslcblas
#define _GNU_SOURCE
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <error.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <argp.h>
#include <sys/types.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_multimin.h>
#define die abort
#include "tmalloc.h"
#include "config.h"

#define BUFSAMP (60*25000)
#define SPIKESAMP 50
#define MAX_NUM_EVENTS "500000"

#define SDSAMP SPIKESAMP

int malloc_debug;

static short buf[BUFSAMP];

const char *argp_program_version = "muph, part of " PACKAGE_STRING;
const char *argp_program_bug_address = "<roconnor@health.usf.edu>";
static char args_doc[] = "FILE.chan";
static char doc[] = ("generate phrenic spikes (code 89) from FILE.chan and add them to FILE.edt"
                     "\vOn successful completion, an 'R' will be appended to FILE.status.\n"
                     "\nThe maximize option may cause a spike on every tick at high amplitudes.\n"
                     "If the high amplitudes are brief glitches, this may be desirable.");
struct arguments
{
  char *chanfile;                /* ARG1 & ARG2 */
  bool maximize;
};

static struct argp_option options[] = {
  {"maximize",        'm',     0, 0, "Generate about " MAX_NUM_EVENTS " spikes" },
  { 0 }
};

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'm':
      arguments->maximize = true;
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num >= 1)
        argp_usage (state);
      arguments->chanfile = arg;
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

  TREALLOC (new_file, strlen (file_name) + strlen (suffix) + 1);
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

double noise_per_sample, max_sample_energy;

static int doplot;

static void
plot (int *hist, int cnt, double mean, double sd, int N)
{
  return;
  doplot = 0;
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

typedef struct {int *hist; int bincnt;} Hist;

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

static double
sumsq_wrap (const gsl_vector *v, void *params)
{
  Hist *p = (Hist *)params;

  double mean = gsl_vector_get(v, 0);
  double sd   = gsl_vector_get(v, 1);
  int    N    = gsl_vector_get(v, 2);

  return sumsq (p->hist, p->bincnt, mean, sd, N);
}

static void
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
  if (iter >= iter_limit)
    error_at_line (0, 0, __FILE__, __LINE__, "minimizer hit limit");
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

static int
modal_sd (double *meanp)
{
  int max_hidx = 0;
  long long sqsum = 0;
  double sum = 0, mean;
  int n;
  int hist[65536];
  
  memset (hist, 0, sizeof hist);

  for (n = 0; n < BUFSAMP; n++)
    sum += buf[n];
  *meanp = mean = sum / BUFSAMP;
  for (short *bufp = buf; bufp + SDSAMP <= buf + BUFSAMP; bufp += SDSAMP) {
    double sd;
    int hidx;
    sqsum = 0;
    for (n = 0; n < SDSAMP; n++)
      sqsum += pow (bufp[n] - mean, 2);
    sd = sqrt ((double)sqsum / SDSAMP);
    hidx = nearbyint (sd);
    if (hidx > max_hidx)
      max_hidx = hidx;
    hist[hidx]++;
  }
  if (doplot)
    plot (hist, max_hidx, 0, 0, 0);
  {
    int max = 0, n_at_max = 0;
    for (n = 1; n < max_hidx; n++)
      if (hist[n] > max)
	max = hist[n_at_max = n];
    plot (hist, max_hidx, n_at_max, n_at_max / 8, 30000);
    double mean = n_at_max;
    double sd = n_at_max / 8;
    int N = 30000;
    normfitls (hist, max_hidx, &mean, &sd, &N);
    plot (hist, max_hidx, mean, sd, N);

    return n_at_max;
  }
}

static double
stdev (double mean)
{
  double sqsum = 0;
  for (int i = 0; i < BUFSAMP; i++)
    sqsum += pow (buf[i] - mean, 2);
  return sqrt (sqsum / BUFSAMP);
}

int
main (int argc, char **argv)
{
  struct arguments arg;
  arg.maximize = false;
  argp_parse (&argp, argc, argv, 0, 0, &arg);
  int max_num_events = atoi (MAX_NUM_EVENTS);
  FILE *f, *of, *edt;
  unsigned long starting_sample, ecode, etime;
  char *edtname;
  char outname[] = "muph_XXXXXX";
  gsl_rng_env_setup ();
  gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);

  (f = fopen (arg.chanfile, "rb")) || DIE;
  edtname = change_filetype_trim (arg.chanfile, ".edt", 0);
  
  if ((edt = fopen (edtname, "r")) == NULL)
    error_at_line (1, errno, __FILE__, __LINE__, "Can't open %s for read", edtname);

  static double ebuf[SPIKESAMP];
  int ie = 0;
  double ebufsum = 0;
  double ebufsummax = 0;
  starting_sample = 0;
  int sampsread;
  double esum = 0;
  double modal_sum = 0;
  int modal_count = 0;
  double sdsum = 0;
  int ticksperspike = SPIKESAMP / 2.5;
  int loopcount = 0;
  double sdmax = 0;
  while ((sampsread = fread (buf, sizeof *buf, BUFSAMP, f)) > 0) {
    static double mean;
    static int sd;
    if (sampsread == BUFSAMP) {
      sd = modal_sd (&mean);
      modal_sum += sd;
      modal_count++;
    }
    double rsd = stdev (mean);
    if (sd > rsd)
      printf ("MODAL SD GREATER THAN OVERALL SD at %ld\n", ftell (f));
    double sdzero = sqrt (2) * sd;
    for (int i = 0; i < sampsread; i++) {
      double e = pow (buf[i] - mean, 2);
      esum += e;
      //      printf ("%.6f: %10.0f %10.0f\n", (starting_sample + i) / 25000., e, ebufsum / SPIKESAMP);
      ebufsum -= ebuf[ie];
      ebuf[ie] = e;
      ebufsum += e;
      loopcount++;
      if ((i + 1) % SPIKESAMP == 0) {
        double sd = sqrt (ebufsum / SPIKESAMP) - sdzero;
        if (sd > sdmax)
          sdmax = sd;
        if (sd > 0)
          sdsum += sd * ticksperspike;
        if (ebufsum > ebufsummax)
          ebufsummax = ebufsum;
      }
      ie = (ie + 1) % SPIKESAMP;
    }
    starting_sample += sampsread;
  }
  double eperspike = ebufsummax;
  //  printf ("sdsum: %g, loopcount: %d\n", sdsum, loopcount);

  eperspike = (modal_sum / modal_count) / SDSAMP * SPIKESAMP;

  double spikecnt = sdsum / sdmax;
  //  printf ("sdsum: %g, sdmax: %g, spikecount: %.1f\n", sdsum, sdmax, spikecnt);
  if (spikecnt > max_num_events || arg.maximize)
    sdmax = sdsum / max_num_events;
  spikecnt = sdsum / sdmax;
  //  printf ("sdsum: %g, sdmax: %g, spikecount: %.1f\n", sdsum, sdmax, spikecnt);

  int fd;
  (fd = mkstemp (outname)) >= 0 || DIE;
  (of = fdopen (fd, "w")) || DIE;

 LOOP:
  rewind (edt);
  for (int n = 0; n < 3; n++)
    if (fscanf (edt, "%5lu%10lu", &ecode, &etime) != 2) {
      etime = ULONG_MAX;
      break;
    }
  rewind (of);
  fprintf (of, "   33   3333333\n");
  fprintf (of, "   33   3333333\n");
  
  ie = 0;
  ebufsum = 0;
  starting_sample = 0;
  rewind (f);
  int spike_count = 0;
  memset (ebuf, 0, sizeof ebuf);
  while ((sampsread = fread (buf, sizeof *buf, BUFSAMP, f)) > 0) {
    static double mean;
    static int msd;
    if (sampsread == BUFSAMP)
      msd = modal_sd (&mean);
    double sdzero = sqrt (2) * msd;
    for (int i = 0; i < sampsread; i++) {
      double e = pow (buf[i] - mean, 2);
      ebufsum -= ebuf[ie];
      ebuf[ie] = e;
      ebufsum += e;
      if ((i + 1) % SPIKESAMP == 0) {
        double sd = sqrt (ebufsum / SPIKESAMP) - sdzero;
        if (sd > 0) {
          double p = sd / sdmax;
          for (int j = 0; j < ticksperspike; j++)
            if (gsl_rng_uniform (rng) < p) {
              int ptime = nearbyint ((starting_sample + (i + 1) - SPIKESAMP) / 2.5 + j);

              while (etime < ptime) {
                fprintf (of, "%5lu%10lu\n", ecode, etime);
                if (fscanf (edt, "%5lu%10lu", &ecode, &etime) != 2)
                  etime = ULONG_MAX;
              }
              fprintf (of, "%5d%10d\n", 89, ptime);
              spike_count++;
            }
        }
      }
      ie = (ie + 1) % SPIKESAMP;
    }
    starting_sample += sampsread;
  }
  //  printf ("%d spikes\n", spike_count);
  if (spike_count > max_num_events) {
    fflush (of);
    if (ftruncate (fd, 0) != 0)
      error_at_line (1, errno, __FILE__, __LINE__, "ftruncate");
    goto LOOP;
  }
  while (etime < ULONG_MAX) {
    fprintf (of, "%5lu%10lu\n", ecode, etime);
    if (fscanf (edt, "%5lu%10lu", &ecode, &etime) != 2)
      etime = ULONG_MAX;
  }
  fclose (f);
  fclose (of);
  fclose (edt);
  if (remove (edtname))
    error_at_line (1, errno, __FILE__, __LINE__, "error removing %s, aborting", edtname);
  if (rename (outname, edtname))
    error_at_line (1, errno, __FILE__, __LINE__, "error renaming %s to %s, aborting", outname, edtname);
  
  {
    FILE *f;
    (f = fopen (change_filetype_trim (arg.chanfile, ".status", 0), "a")) || DIE;
    putc ('R', f);
    fclose (f);
  }
  return 0;
}
