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

/* spiketime_to_region.c */


#include "nde.h"
#include <string.h>
#include "region.h"

int malloc_debug;

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

RegionData *
region_read (char *filename)
{
  FILE *f;
  RegionData *p;

  TCALLOC (p, 1);
  if ((f = fopen (filename, "rb")) == 0)
    return p;

  fread (&p->peak_count, sizeof p->peak_count, 1, f) == 1 || DIE;
  p->peak_count >= 0 || DIE;
  fseek (f, sizeof *p->peak * (size_t)p->peak_count, SEEK_CUR);
  fread (&p->region_count, sizeof p->region_count, 1, f) == 1 || DIE;
  p->region_count >= 0 || DIE;
  TMALLOC (p->region, p->region_alloc = p->region_count);
  fread (p->region, sizeof *p->region, (size_t)p->region_count, f) == (size_t)p->region_count || DIE;
  return p;
}

int
main (int argc, char **argv)
{
  int sample, rn;
  RegionData *p;

  if (argc < 3)
    return 0;
  p = region_read (change_filetype (argv[1], ".regions"));
  sample = (int)floor (atof (argv[2]) * 25000.0);
  for (rn = 0; rn < p->region_count; rn++)
    if (sample >= p->region[rn].sample && sample < p->region[rn].sample + p->region[rn].count) {
      printf ("%d\n", rn);
      exit (0);
    }
  return 0;
}
