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

#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "tmalloc.h"

#define ISAMP 50
#define ISM 51
#define UP .05
#define DN -.05
#define IHZ 200

int malloc_debug;

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

int
main (int argc, char **argv)
{
  double sm[ISM];
  struct {int val; int time;} ix[ISAMP];
  int itics, i, n, valcode, time, ichn, ichno, iflg, iminpt, imaxpt;
  double a, y, scale;
  FILE *f[2], *outfile;
  char *outname;

  if (argc != 3) {
    printf ("usage: %s edt_file channel\n", argv[0]);
    return 0;
  }

  itics = 10000 / IHZ;
  a = 0;

  for (i = 0; i < ISM; i++) {
    double rk = 4.0 * i / (ISM - 1);

    if      (rk <= 1) sm[i] = rk;
    else if (rk <= 3) sm[i] = 2 - rk;
    else              sm[i] = rk - 4;

    a += fabs (sm[i]);
  }

  for (i = 0; i < ISM; i++)
    sm[i] /= a;

  (f[0] = fopen (argv[1], "r")) || DIE;
  (f[1] = tmpfile ()) || DIE;

  ichn = atoi (argv[2]);
  ichno = ichn * 4096;

  iflg = n = 0;
  iminpt = imaxpt = -1;

  while (fscanf (f[0], "%5d%10d", &valcode, &time) == 2) {
    if (valcode / 4096 != ichn)
      continue;
    if ((ix[n].val = valcode - ichno) > 2047)
      ix[n].val -= 4096;
    ix[n++].time = time;
    if (n < ISAMP)
      continue;
    if (imaxpt < 0) {
      imaxpt = 0;
      for (i = 1; i < ISAMP; i++)
	if (ix[i].val > ix[imaxpt].val)
	  imaxpt = i;
    }
    else if (ix[n-1].val > ix[imaxpt].val)
      imaxpt = n - 1;
    if (iminpt < 0) {
      iminpt = 0;
      for (i = 1; i < ISAMP; i++)
	if (ix[i].val < ix[iminpt].val)
	  iminpt = i;
    }
    else if (ix[n-1].val < ix[iminpt].val)
      iminpt = n - 1;
    if (ix[ISAMP-1].time - ix[ISAMP-2].time > itics * 2) {
      iflg = n = 0;
      iminpt = imaxpt = -1;
      continue;
    }
    y = 0;
    if ((scale = (ix[imaxpt].val - ix[iminpt].val)) != 0)
      for (i = 0; i < ISAMP; i++)
	y += (double)ix[ISAMP - 1 - i].val / scale * sm[i];

    if (y > UP && iflg == 0) {
      iflg = 1;
      fprintf (f[1], "%5d%10d\n", 99, ix[ISAMP-ISM/2-1].time);
    }
    if (y < DN && iflg == 1)
      iflg = 0;

    for (i = 0; i < ISAMP - 1; i++)
      ix[i] = ix[i + 1];
    imaxpt--;
    iminpt--;
    n = ISAMP - 1;
  }
      
  outname = change_filetype_trim (argv[1], ".tmp", 0);

  (outfile = fopen (outname, "wb")) || DIE;

  rewind (f[0]);
  rewind (f[1]);
  if (fscanf (f[0], "%5d%10d", &valcode, &time));
  if (valcode == 0) 
    if (fscanf (f[0], "%5d%10d", &valcode, &time));
  if (fscanf (f[0], "%5d%10d", &valcode, &time));
  fprintf (outfile, "%5d%10d\n", valcode, time);
  fprintf (outfile, "%5d%10d\n", valcode, time);

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
  remove (argv[1]);
  rename (outname, argv[1]);

  {
    FILE *f;
    (f = fopen (change_filetype_trim (argv[1], ".status", 0), "a")) || DIE;
    putc ('C', f);
    fclose (f);
  }
  return 0;
}
