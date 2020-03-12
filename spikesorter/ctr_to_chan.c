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
#include "region.h"
#include <string.h>

int
main (int argc, char **argv)
{
  FILE *f;
  int count, n, i, a;
  Center ctr;
  short val[SNAPLEN];
  static short zero[SNAPLEN];
  static short edge[SNAPLEN];

  argc > 1 || DIE;

  edge[1] = 32767;
  edge[2] = -32768;

  edge[SNAPLEN-2] = 32767;
  edge[SNAPLEN-3] = -32768;

  for (a = 1; a < argc; a++) {
    (f = fopen (argv[a], "rb")) || DIE;
    fread (&count, sizeof count, 1, f) == 1 || DIE;
    fprintf (stderr, "%s: %d centers\n", argv[a], count);
    for (n = 0; n < count; n++) {
      fread (&ctr, sizeof ctr, 1, f) == 1 || DIE;
      for (i = 0; i < SNAPLEN; i++) {
	if (ctr.raw[i] >= 32767)
	  val[i] = 32767;
	else if (ctr.raw[i] <= -32768)
	  val[i] = -32767;
	else
	  val[i] = floor (ctr.raw[i] + .5);
      }
      fwrite (val, sizeof val, 1, stdout) == 1 || DIE;
      fwrite (edge, sizeof edge[0], 4, stdout) == 4 || DIE;
    }
    fclose (f);
    fwrite (zero, sizeof zero, 1, stdout) == 1 || DIE;
  }
  return 0;
}
