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

#include "nde.h"
#include <string.h>
#include <float.h>

typedef struct {double loc; int id; int pad;} Entry;
int malloc_debug;

static double
pfrac (double r, double t, int N, double *conf)
{
  double minp = 0, maxp = .5, p, f, err = .5;

  p = .5;
  f = pow ((1 - (1 - p) * r * t), (p * N)) * pow ((1 - p * r * t), ((1 - p) * N));
  if (f > err) {
    *conf = 1 - f;
    //    printf ("line %d: %f %f\n", __LINE__, p, f);
    return .5;
  }
  while (1) {
    p = (minp + maxp) / 2;
    if (maxp - minp < DBL_EPSILON * minp) {
      *conf = 1 - err;
      return p;
    }
    f = pow ((1 - (1 - p) * r * t), (p * N)) * pow ((1 - p * r * t), ((1 - p) * N));
    //    printf ("line %d: %f %f\n", __LINE__, p, f);
    if (f < err)
      maxp = p;
    else if (f > err)
      minp = p;
    else {
      *conf = 1 - err;
      return p;
    }
  }
}

static long filesize (FILE *f)
{
  long bytes;

  fseek (f, 0, SEEK_END);
  bytes = ftell (f);
  rewind (f);
  return bytes;
}

static char *
change_filetype_trim (char *file_name, char *suffix, int trim)
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

static Entry *
dofile (char *filename, int *countp)
{
  FILE *f;
  int count;
  Entry *ddt;

  (f = fopen (change_filetype_trim (filename, ".pos", 0), "r")) || DIE;
  count = filesize (f) / sizeof (Entry);
  filesize (f) == count *  sizeof (Entry) || DIE;

  TMALLOC (ddt, count);
  fread (ddt, sizeof *ddt, count, f) == count || DIE;
  *countp = count;
  return ddt;
}

static int
max_id (Entry *e, int count)
{
  int n, max;
  
  for (max = n = 0; n < count; n++)
    if (e[n].id > max)
      max = e[n].id;
  return max;
}

int main (int argc, char **argv)
{
  Entry *e;
  int count, n, ccnt, id;
  char *dounit;
  double isi, T;
  double last, min;
  int scnt, ltcnt;
  const double t = 2;
  double binstart, binend;

  for (n = 0; n < argc; n++)
    fprintf (stderr, "%s ", argv[n]);
  fprintf (stderr, "\n");

  if (argc < 4)
    return 1;
  e = dofile (argv[1], &count);
  ccnt = max_id (e, count);
  binstart = atof (argv[2]) * (SF / 1000);
  if ((binend = atof (argv[3]) * (SF / 1000)) == 0)
    binend = HUGE_VAL;
  TCALLOC (dounit, ccnt + 1);
  for (n = 4; n < argc; n++)
    if ((id = atoi (argv[n])) >= 0 && id <= ccnt)
      dounit[id] = 1;
  
  min = HUGE_VAL;
  for (last = 0, scnt = ltcnt = n = 0; n < count; n++) {
    if ((id = e[n].id) >= 0 && id <= ccnt && dounit[id] && e[n].loc >= binstart && e[n].loc <= binend) {
      scnt++;
      if (last) {
	isi = (e[n].loc - last) / (SF / 1000);
	if (isi < t)
	  ltcnt++;
	if (isi < min)
	  min = isi;
      }
      last = e[n].loc;
    }
  }
  T = (double)e[count-1].loc / (SF / 1000);

  {
    double c, f, r, p, conf;
    if (last == 0)
      return 0;
    printf ("%d spikes\n", scnt);
    f = (double)(ltcnt) / scnt;
    r = scnt / T;
    p = pfrac (r, min - .05, scnt, &conf);
    if (f > 1 - exp (-r * t))
      printf ("DIRTY\n%d < 2 ms\n%.3f EPR", ltcnt, f / (1 - exp (-r * t)));
    else if (min >= t) {
      printf ("CLEAN\n%.1f ms min ISI\n%.6f CCR", min, p);
    }
    else if (f <= .5 * r * t) {
      c = f / (2 * r * t);
      p = (1 - sqrt (1 - 4 * c)) / 2;
      printf ("MARGINAL\n%d < 2 ms\n%.6f MWR", ltcnt, p);
    }
    else {
      f = (double)ltcnt / scnt;
      n = 1 / (1 - f / (r * t));
      printf ("DIRTY\n%d < %.2f ms\n%d DNU", ltcnt, t, n);
    }
    printf ("\n");
  }
  return 0;
}
