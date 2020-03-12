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

/* dmx_split.c */

#define _GNU_SOURCE
#include "lfs.h"
#include "tmalloc.h"
#include <string.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifdef _WIN32
#include <io.h>
#define mkdir(file,mode) mkdir (file)
#else
#include <unistd.h>
#endif

#include "dmx.h"

#undef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))

int malloc_debug;

static FILE **
open_files (char *fname, unsigned channel_count, char *enabled)
{
  char *fmt = "./split/%s_%02d.chan";
  char *namebuf, *inname, *s, *base_name;
  int n, overwrite = 0;
  FILE **file_list;

  TMALLOC (file_list, channel_count + 1);
  (inname = strdup (fname)) || DIE;
  for (n = 0; inname[n]; n++)
    if (inname[n] == '.')
      inname[n] = '_';
  base_name = inname;
  while ((s = strpbrk (base_name, "/\\")))
    base_name = s + 1;
  TMALLOC (namebuf, strlen (fname) + strlen (fmt) + 1);
  channel_count < 100 || DIE;	/* consider filename format & TMALLOC size before increasing this */
  for (n = 1; n <= channel_count; n++)
    if (enabled[n]) {
      sprintf (namebuf, fmt, base_name, n);
      if (access (namebuf, F_OK) == 0 && !overwrite) {
	printf ("overwrite existing output file(s)? (n/y) ");
	fflush (stdout);
	if (getchar () != 'y')
	  exit (0);
	while (getchar () != '\n')
	  ;
	overwrite = 1;
      }
      (file_list[n] = fopen (namebuf, "wb")) || DIE;
    }
    else file_list[n] = 0;

  free (inname);
  free (namebuf);
  return file_list;
}

int
main (int argc, char **argv)
{
  FILE *f, *of;
  GlobalHeader *ghdr;
  BoardHeader *bhdr;
  int n, i, sample_rate, channel_count, fpos;
  int *enabled_channels;
  char *enabled, *dmx_ofile;
  FILE **file_list;

  if (argc != 2) {
    printf ("usage: %s dmx_file_name\n", argv[0]);
    exit (0);
  }

  if ((f = fopen64 (argv[1], "rb")) == 0) {
    fprintf (stderr, "Can't open %s: %s\n", argv[1], strerror (errno));
    exit (1);
  }

  if (access ("./split", F_OK) != 0)
    if (mkdir("./split", 0777) != 0) {
      fprintf (stderr, "Can't mkdir ./split: %s\n", strerror (errno));
      exit (1);
    }
    
  if (access ("./clean", F_OK) != 0)
    if (mkdir("./clean", 0777) != 0) {
      fprintf (stderr, "Can't mkdir ./clean: %s\n", strerror (errno));
      exit (1);
    }

  TMALLOC (dmx_ofile, strlen(argv[1]) + 12);
  sprintf(dmx_ofile, "./clean/%s", argv[1]);
  if ((of = fopen64 (dmx_ofile, "wb")) == 0) {
    fprintf (stderr, "Can't open %s: %s\n", dmx_ofile, strerror (errno));
    exit (1);
  }

  TMALLOC (ghdr, 1);
  fread64 (ghdr, sizeof *ghdr, 1, f) == 1 || DIE;

  TMALLOC (bhdr, ghdr->board_count);
  TMALLOC (enabled_channels, ghdr->board_count);
  fread64 (bhdr, sizeof *bhdr, ghdr->board_count, f) == ghdr->board_count || DIE;
  sample_rate = 0;  
  channel_count = ghdr->board_count * 8;
  TMALLOC (enabled, channel_count + 1);
  for (n = 0; n < ghdr->board_count; n++) {
    enabled_channels[n] = 0;
    for (i = 0; i < 8; i++) {
      if ((enabled[n*8+i+1] = bhdr[n].chan[i].enable)) {
	enabled_channels[n]++;
	if (sample_rate == 0)
	  sample_rate = bhdr[n].chan[i].sample_rate;
	else
	  bhdr[n].chan[i].sample_rate == sample_rate || DIE;
	bhdr[n].chan[i].board == n || DIE;
	bhdr[n].chan[i].channel == i || DIE;
      }
    }
    enabled_channels[n] <= bhdr[n].channel_count || DIE;
  }
  file_list = open_files (argv[1], channel_count, enabled);

  fwrite (ghdr, sizeof *ghdr, 1, of) == 1 || DIE;
  fwrite (bhdr, sizeof *bhdr, ghdr->board_count, of) == ghdr->board_count || DIE;

  fpos = ftello64 (f);
  fpos == sizeof (*ghdr) + sizeof (*bhdr) * ghdr->board_count || DIE;

  fpos = (fpos / (1 << 16) + 1) * (1 << 16);
  fseeko64 (f, fpos, SEEK_SET);

  {
    int *swipe_size, swipe_size_total, bufsamples, j, ecb;
    int samples_per_100ms, board, chidx, chnum, inswipe, oswipe;
    short *buf, *bufp, *obufp, **obuf;
    time_t now, then;
    int swipes_per_buf = 5, swipecnt, swipe, oswipe_count;

    printf ("swipes_per_buf: %d\n", swipes_per_buf);
    TMALLOC (swipe_size, ghdr->board_count);
    samples_per_100ms = sample_rate / 10;
    for (swipe_size_total = n = 0; n < ghdr->board_count; n++)
      swipe_size_total += (swipe_size[n] = enabled_channels[n] * 2 * samples_per_100ms);
    bufsamples = swipe_size_total / sizeof *buf * swipes_per_buf;
    TMALLOC (buf, bufsamples);
    TMALLOC (obuf, channel_count + 1);
    for (chnum = 1; chnum <= channel_count; chnum++)
      if (enabled[chnum])
	TMALLOC (obuf[chnum], bufsamples);
      else
	obuf[chnum] = 0;
    then = 0;
    oswipe = 0;
    oswipe_count = bufsamples / samples_per_100ms;
    for (inswipe = 0; inswipe < ghdr->swipe_count; inswipe += swipecnt) {
      if (1)
	if ((now = time(0)) > then)
	  {
	    printf ("  %5d of %d swipes\r", inswipe, ghdr->swipe_count);
	    fflush (stdout);
	    then = now;
	  }
      swipecnt = MIN (swipes_per_buf, ghdr->swipe_count - inswipe);
      //      swipe_size_total * swipecnt / sizeof *buf <= bufsamples || DIE;
      fread64 (bufp = buf, swipe_size_total * swipecnt, 1, f) == 1 || DIE;
      for (swipe = 0; swipe < swipecnt; swipe++) {
	for (board = 0; board < ghdr->board_count; board++) {
	  ecb = enabled_channels[board];
	  chnum = board * 8 + 1;
	  for (chidx = 0; chidx < ecb; chidx++, chnum++) {
	    while (!enabled[chnum])
	      chnum++;
	    chnum < (board + 1) * 8 + 1 || DIE;
	    obufp = obuf[chnum] + oswipe * samples_per_100ms;
//	    obufp + samples_per_100ms - obuf[chnum] <= bufsamples || DIE;
//	    bufp + (samples_per_100ms - 1) * ecb - buf < bufsamples || DIE;
	    for (i = chidx, j = 0; j < samples_per_100ms; j++, i += ecb)
	      obufp[j] = bufp[i];
	  }
	  bufp += swipe_size[board] / sizeof *bufp;
	}
	if (++oswipe == oswipe_count || (swipe + 1 == swipecnt && inswipe + swipecnt == ghdr->swipe_count)) {
	  int obuf_count = oswipe * samples_per_100ms;
	  oswipe = 0;
	  obuf_count <= bufsamples || DIE;
	  for (chnum = 1; chnum <= channel_count; chnum++)
	    if (enabled[chnum])
	      fwrite (obuf[chnum], sizeof *obufp, obuf_count, file_list[chnum]) == obuf_count || DIE;
	}
      }
    }  
    printf ("  %5d of %d swipes\r", inswipe, ghdr->swipe_count);
    printf ("\n");
  }
  return 0;
}
