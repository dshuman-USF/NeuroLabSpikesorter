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

/* residue.c */

#include "nde.h"
#include <string.h>

static int snap_count, point_count, rawsize;
static int *bufidx_array;
static double *x_array;
static float *raw, *buf0;
static int eval_count;

double
local_cubic (double x, float *y, int count, double edge)
{
  double a, b, c, d, y1, y2, y3, y4, f;
  int n;
  int n3;
  double x0;

  x0 = x;
  (x >= 0 && (int)ceil (x) <= count) || DIE;
  n = (int)floor (x);
  if ((double)n == x) {
    //    printf ("exact %f\n", y[n]);
    return y[n];
  }
  x -= n;
  if (edge == HUGE_VAL)
    y1 = n == 0 ? y[0] + (y[0] - y[1]) : y[n-1];
  else
    y1 = n == 0 ? edge : y[n-1];
  y2 = y[n++];
  n3 = n;
  y3 = y[n++];
  if (edge == HUGE_VAL)
    y4 = n == count ? y[n-1] + (y[n-1] - y[n-2]) : y[n];
  else
    y4 = n == count ? edge : y[n];
  a = (-y1 + 3 * y2 - 3 * y3 + y4) / 2;
  b = (2 * y1 - 5 * y2 + 4 * y3 - y4) / 2;
  c = (y3 - y1) / 2;
  d = y2;

  f = a * x * x * x + b * x * x + c * x + d;
  if (0)
    printf ("x: %f, y: %6.0f %6.0f %6.0f %6.0f, abcd: %8.1f %8.1f %8.1f %8.1f , f: %13.6f, x0: %f, n3: %d\n",
	    x, y1, y2, y3, y4, a, b, c, d, f, x0, n3);
  return f;
}

double
residue (float *avg)
{
  int snapidx, bufidx, n;
  double tmp, sum, x;

  rawsize >= 0 || DIE;
  memcpy (raw, buf0, (size_t)rawsize);
  for (snapidx = 0; snapidx < snap_count; snapidx++) {
    bufidx = bufidx_array[snapidx];
    x = x_array[snapidx];
    for (n = 0; n < SNAPLEN-1; n++)
      raw[bufidx++] -= local_cubic (n + x, avg, SNAPLEN, mean) - mean;
  }
  sum = 0;
  for (n = 0; n < point_count; n++)
    tmp = raw[n] - mean, sum += tmp * tmp;
  eval_count++;
  return sum + 1;
}

void
changed (float *avg, int print)
{
  int snapidx, bufidx, n;
  double x;
  static float *old;

  TREALLOC (old, point_count);
  rawsize >= 0 || DIE;
  memcpy (old, raw, (size_t)rawsize);
  memset (raw, 0, (size_t)rawsize);
  for (snapidx = 0; snapidx < snap_count; snapidx++) {
    bufidx = bufidx_array[snapidx];
    x = x_array[snapidx];
    for (n = 0; n < SNAPLEN-1; n++) {
      if (bufidx == 30)
	printf ("%3d: %f\n", bufidx, n+x);
      raw[bufidx++] += local_cubic (n + x, avg, SNAPLEN, mean) - mean;
    }
  }
  if (print) {
    printf ("changed:");
    for (n = 0; n < point_count; n++)
      if (raw[n] != old[n])
	printf ("%3d: %12.6f %5.0f %12.6f\n", n, raw[n], buf0[n], raw[n] - buf0[n]);
  }
}

float *
init_residue (SnapList *sl, int first_snap, int last_snap)
{
  int n, snapidx, first_sample;
  double peakloc, offset;
  float *bufptr;
  static float p0[SNAPLEN];

  snap_count = last_snap - first_snap + 1;
  TREALLOC (bufidx_array, snap_count);
  TREALLOC (x_array, snap_count);
  
  first_sample = (int)ceil (sl->snap[first_snap]->sample) - PRESAMPLES - starting_sample;
  for (snapidx = first_snap; snapidx <= last_snap; snapidx++) {
    peakloc = sl->snap[snapidx]->sample;
    n = snapidx - first_snap;
    x_array[n] = ceil (peakloc) - peakloc;
    bufidx_array[n] = (int)ceil (peakloc) - PRESAMPLES - first_sample;
  }
  point_count = bufidx_array[snap_count - 1] - bufidx_array[0] + SNAPLEN-1;
  TREALLOC (raw, point_count);
  buf0 = buf + first_sample;
  //  printf ("buf0: %d, points: %d\n", buf0 - buf, point_count);
  rawsize = point_count * (int)sizeof *raw;

  bufptr = buf + ((int)floor (sl->snap[first_snap]->sample) - PRESAMPLES - starting_sample);
  offset = sl->snap[first_snap]->sample - floor (sl->snap[first_snap]->sample);
  for (n = 0; n <= PRESAMPLES; n++)
    p0[n] = local_cubic (n + offset + 1, bufptr - 1, SNAPLEN+2, mean);
  bufptr = buf + ((int)floor (sl->snap[last_snap]->sample) + 1 - starting_sample);
  offset = sl->snap[last_snap]->sample - floor (sl->snap[last_snap]->sample);
  for (n = 0; n < SNAPLEN - PRESAMPLES - 1; n++)
    p0[n + PRESAMPLES + 1] = local_cubic (n + offset + 1, bufptr - 1, SNAPLEN+2, mean);

  if (0)
    for (n = 0; n < SNAPLEN; n++)
      p0[n] = mean;

  return p0;
}
