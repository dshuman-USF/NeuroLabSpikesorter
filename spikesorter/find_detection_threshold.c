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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <limits.h>

static int sample_count, noise_count;
double mean;
static short *mapbuf, *bufp;

#define SZ SNAPLEN
#define SQR(x) ((x)*(x))

static int hist[65536];

static void
map_file (char *name)
{
  int fd;

  fd = open (name, O_RDONLY);
  if (fd < 0) {
    fprintf (stderr, "Error opening %s: %s\n", name, strerror (errno));
    exit (1);
  }
  {
    struct stat s;
    int retval;
    ((retval = fstat (fd, &s)) >= 0 && s.st_size > 0) || DIE;
    sample_count = s.st_size / sizeof *mapbuf;
    mapbuf = mmap (0, s.st_size, PROT_READ , MAP_SHARED, fd, 0);
  }
}


static double
threshold_sd (int threshold) {
  int n, resume = 0, count = 0;
  double sqsum = 0;
  double sd;

  for (n = 0; n < BEFORE; n++)
    if (abs (mapbuf[n] - mean) >= threshold)
      resume = n + 1 + AFTER;

  for (n = 0; n < sample_count; n++) {
    if (n + BEFORE < sample_count && abs (mapbuf[n + BEFORE] - mean) >= threshold)
      resume = n + BEFORE + 1 + AFTER;
    if (n >= resume) {
      sqsum += SQR (mapbuf[n] - mean);
      count++;
    }
  }
  noise_count = count;
  sd = sqrt (sqsum / count);
  printf ("threshold %d: sd %f\n", threshold, sd);
  return sd;
}

static int
modal_sd ()
{
  int max_hidx = 0;
  long long sqsum = 0;
  long long sum = 0;
  int n;

  for (n = 0; n < sample_count; n++)
    sum += mapbuf[n];
  mean = (double) sum / sample_count;
  printf ("mean = %f\n", mean);
  for (bufp = mapbuf; bufp + SZ <= mapbuf + sample_count; bufp += SZ) {
    double sd;
    int hidx;
    sqsum = 0;
    for (n = 0; n < SZ; n++)
      sqsum += SQR (bufp[n] - mean);
    sd = sqrt ((double)sqsum/SZ);
    hidx = floor (sd + .5);
    if (hidx > max_hidx)
      max_hidx = hidx;
    hist[hidx]++;
  }
  {
    int max = 0, n_at_max = 0;
    for (n = 1; n < max_hidx; n++)
      if (hist[n] > max)
	max = hist[n_at_max = n];
    return n_at_max;
  }
}


static int
optimum_threshold (int target_sd)
{
  int thr = floor (3.227167 * target_sd + .5);
  double slope = 0.127339;
  double sd = threshold_sd (thr);

  double closest_d = fabs (sd - target_sd);
  int closest_thr = thr;
  double closest_sd = sd;
  int closest_noise_count = noise_count;

  int second_closest_thr = 0;
  double second_closest_sd = sd + 65536;
  double second_closest_d = 65536;

  if (sd == target_sd)
    return thr;

  while (1) {
    double d;

    thr = fabs (closest_thr + (target_sd - closest_sd) / slope + .5);
    if (thr == closest_thr)
      thr += target_sd > closest_sd ? 1 : -1;

    if ((sd = threshold_sd (thr)) == target_sd)
      return thr;
    d = fabs (sd - target_sd);
    if (d < closest_d || (d == closest_d && thr > closest_thr)) {
      second_closest_d = closest_d;
      second_closest_thr = closest_thr;
      second_closest_sd = closest_sd;
      closest_d = d;
      closest_thr = thr;
      closest_sd = sd;
      closest_noise_count = noise_count;
     }
    else if (d < second_closest_d) {
      second_closest_d = d;
      second_closest_thr = thr;
      second_closest_sd = sd;
    }
    if (abs (closest_thr - second_closest_thr) == 1) {
      noise_count = closest_noise_count;
      return closest_thr;
    }
    slope = (closest_sd - second_closest_sd) / (closest_thr - second_closest_thr);
  }
}

int
find_detection_threshold (char *file_name)
{
  int threshold;

  map_file (file_name);
  threshold = optimum_threshold (modal_sd ());
  munmap (mapbuf, sample_count * sizeof *mapbuf);
  return threshold;
}
