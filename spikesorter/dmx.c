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

/* dmx.c */

#include "nde.h"
#include <string.h>
#include <time.h>
#include "dmx.h"

int malloc_debug;
#define UI(x) *(unsigned int *)&(x)

# define SWAP32(x) \
        UI(x)=((((UI(x)) & 0xff000000) >> 24) | (((UI(x)) & 0x00ff0000) >>  8) |    \
         (((UI(x)) & 0x0000ff00) <<  8) | (((UI(x)) & 0x000000ff) << 24))

#define US(x) *(unsigned short *)&(x)

#define SWAP16(x) US(x) = ((US(x) >> 8) & 0xff) | ((US(x) & 0xff) << 8)

static void
swap_dinfo (DisplayInfo *d)
{
  SWAP32 (d->type);
  SWAP32 (d->board);
  SWAP32 (d->channel);
  SWAP32 (d->start);
  SWAP32 (d->end);
  SWAP32 (d->vscale);
  SWAP32 (d->voffset);
  SWAP32 (d->grid_enable);
  SWAP32 (d->eng_unit_enable);
  SWAP32 (d->color);
}

static void
swap_ghdr (GlobalHeader *g)
{
  int n;

  SWAP32 (g->version);
  SWAP32 (g->board_count);
  SWAP32 (g->channel_count);
  SWAP32 (g->swipe_count);
  SWAP32 (g->display_count);
  g->display_count <= 8 || DIE;
  for (n = 0; n < g->display_count; n++)
    swap_dinfo (&g->display_info[n]);
  SWAP16 (g->year);
  SWAP16 (g->month);
  SWAP16 (g->wday);
  SWAP16 (g->mday);
  SWAP16 (g->hour);
  SWAP16 (g->minute);
  SWAP16 (g->second);
  SWAP16 (g->msec);
  SWAP32 (g->note_size);

}

static void
swap_chan (ChannelInfo *c)
{
  SWAP32 (c->board);
  SWAP32 (c->channel);
  SWAP32 (c->bandwidth);
  SWAP32 (c->sample_rate);
  SWAP32 (c->scalar);
  SWAP32 (c->scale);
  SWAP32 (c->enable);
  SWAP32 (c->eng_unit_offset);
  SWAP32 (c->zero_enable);
  SWAP32 (c->zero_value);
}

static void
swap_bhdr (BoardHeader *b)
{
  int n;

  SWAP32 (b->serial);
  SWAP32 (b->channel_count);
  b->channel_count <= 8 || DIE;
  for (n = 0; n < b->channel_count; n++)
    swap_chan (&b->chan[n]);
}

int main (int argc, char **argv)
{
  FILE *f, *of;
  GlobalHeader *ghdr;
  BoardHeader *bhdr;
  char fid[9];
  int n, i, j, sample_rate;
  char *display_type[] = {"Wave", "Bar", "Digital"};
  char *dow[] = {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};
  int verbose = 0, labels = 0, split = 0, enabled = 0;
  int board_channel = 0, board = 0, channel = 0, chidx = -1, fpos;
  char *outname, *suffix;
  char *suffix_fmt = "%02d.chan";
  int *enabled_channels;
  unsigned short v = 1;

  argc > 1 || DIE;
  (f = fopen64 (argv[1], "rb")) || DIE;
  TMALLOC (ghdr, 1);
  fread64 (ghdr, sizeof *ghdr, 1, f) == 1 || DIE;

  if (*(char *)&v == 0)
    swap_ghdr (ghdr);		/* big endian */

  memcpy (fid, ghdr->fid, 8);
  fid[8] = 0;

  if (argc < 3)
    verbose = 1;
  else if (argc == 3 && strcmp (argv[2], "labels") == 0)
    labels = 1;
  else if (argc == 3 && strcmp (argv[2], "enabled") == 0)
    enabled = 1;
  else
    split = 1;

  if (labels || enabled)
    for (n = 0; n < ghdr->display_count; n++) {
      (ghdr->display_info[n].type < 3 && ghdr->display_info[n].type >= 0) || DIE;
    }
  if (verbose) {
    printf ("file identifier: %s\n", fid);
    printf ("version major: %d, minor %d\n", (ghdr->version >> 16) & 0xffff, ghdr->version & 0xffff);
    printf ("boards: %d\n", ghdr->board_count);
    printf ("channels: %d\n", ghdr->channel_count);
    printf ("swipes: %d\n", ghdr->swipe_count);
    printf ("displays: %d\n", ghdr->display_count);
    for (n = 0; n < ghdr->display_count; n++) {
      printf ("display %d\n", n + 1);
      (ghdr->display_info[n].type < 3 && ghdr->display_info[n].type >= 0) || DIE;
      printf ("\ttype: %s\n", display_type[ghdr->display_info[n].type]);
      printf ("\tboard: %d\n", ghdr->display_info[n].board);
      printf ("\tchannel: %d\n", ghdr->display_info[n].channel);
      printf ("\tstart: %d\n", ghdr->display_info[n].start);
      printf ("\tend: %d\n", ghdr->display_info[n].end);
      printf ("\tvscale: %d\n", ghdr->display_info[n].vscale);
      printf ("\tvoffset: %d\n", ghdr->display_info[n].voffset);
      printf ("\tgrid_enable: %d\n", ghdr->display_info[n].grid_enable);
      printf ("\teng_unit_enable: %d\n", ghdr->display_info[n].eng_unit_enable);
      printf ("\tcolor: %d\n", ghdr->display_info[n].color);
    }
    ghdr->wday < 7 || DIE;
    printf ("%s %d/%d/%d %02d:%02d:%02d.%03d\n", dow[ghdr->wday], ghdr->month, ghdr->mday, ghdr->year, ghdr->hour, ghdr->minute, ghdr->second, ghdr->msec);
    ghdr->note[ghdr->note_size] = 0;
    printf ("%s\n", ghdr->note);
  }
  if (split) {
    board_channel = atoi (argv[2]) - 1;
    (board_channel < ghdr->channel_count && board_channel >= 0) || DIE;
    board = board_channel / 8;
    channel = board_channel % 8;
  }
  TMALLOC (bhdr, ghdr->board_count);
  TMALLOC (enabled_channels, ghdr->board_count);
  fread64 (bhdr, sizeof *bhdr, ghdr->board_count, f) == ghdr->board_count || DIE;

  if (*(char *)&v == 0)
    for (n = 0; n < ghdr->board_count; n++) /* big endian */
      swap_bhdr (bhdr + n);

  for (sample_rate = n = 0; n < ghdr->board_count; n++) {
    enabled_channels[n] = 0;
    for (j = i = 0; i < 8; i++) {
      if (bhdr[n].chan[i].enable) {
	enabled_channels[n]++;
	if (sample_rate == 0)
	  sample_rate = bhdr[n].chan[i].sample_rate;
	else {
	  bhdr[n].chan[i].sample_rate == sample_rate || DIE;
	}
	bhdr[n].chan[i].board == n || DIE;
	bhdr[n].chan[i].channel == i || DIE;
	if (n == board && i == channel)
	  chidx = j;
	j++;
      }
    }
    enabled_channels[n] <= bhdr[n].channel_count || DIE;
    if (labels) {
      for (i = 0; i < bhdr[n].channel_count; i++) {
	char cbuf[81];

	memcpy (cbuf, bhdr[n].chan[i].label, 32);
	cbuf[32] = 0;
	printf ("%s\n", cbuf);
      }
    }
    if (enabled)
      for (i = 0; i < bhdr[n].channel_count; i++)
	printf ("%d\n", bhdr[n].chan[i].enable);
    if (verbose) {
      printf ("board: %d\n", n);
      printf ("\tserial: %d\n", bhdr[n].serial);
      printf ("\tchannel_count: %d\n", bhdr[n].channel_count);
      for (i = 0; i < bhdr[n].channel_count; i++) {
	char cbuf[81];
	printf ("\tchannel %d\n", i);
	printf ("\t\tboard: %d\n", bhdr[n].chan[i].board);
	printf ("\t\tchannel: %d\n", bhdr[n].chan[i].channel);
	printf ("\t\tbandwidth: %d\n", bhdr[n].chan[i].bandwidth);
	printf ("\t\tsample_rate: %d\n", bhdr[n].chan[i].sample_rate);
	
	{
	  char *txt[] = {"40V", "10V", "4V", "1V"};
	  (bhdr[n].chan[i].scalar < 4 && bhdr[n].chan[i].scalar >= 0) || DIE;
	  printf ("\t\tscalar: %s\n", txt[bhdr[n].chan[i].scalar]);
	}
	printf ("\t\tscale: %f\n", bhdr[n].chan[i].scale);

	memcpy (cbuf, bhdr[n].chan[i].unit, 8);
	cbuf[8] = 0;
	printf ("\t\tunit: %s\n", cbuf);

	memcpy (cbuf, bhdr[n].chan[i].label, 32);
	cbuf[32] = 0;
	printf ("\t\tlabel: %s\n", cbuf);

	memcpy (cbuf, bhdr[n].chan[i].comment, 80);
	cbuf[80] = 0;
	printf ("\t\tcomment: %s\n", cbuf);

	{
	  char *txt[] = {"Disable", "Enable"};
	  (bhdr[n].chan[i].enable < 2 && bhdr[n].chan[i].enable >= 0) || DIE;
	  printf ("\t\tenable: %s\n", txt[bhdr[n].chan[i].enable]);
	}

	printf ("\t\teng_unit_offset: %f\n", bhdr[n].chan[i].eng_unit_offset);

	{
	  char *txt[] = {"Disable", "Enable"};
	  (bhdr[n].chan[i].enable < 2 && bhdr[n].chan[i].enable >= 0) || DIE;
	  printf ("\t\tzero_enable: %s\n", txt[bhdr[n].chan[i].zero_enable]);
	}
	printf ("\t\tzero_value: %f\n", bhdr[n].chan[i].zero_value);
      }
    }
  }
  chidx >= 0 || DIE;

  if (!split)
    return 0;
  TMALLOC (outname, strlen (argv[1]) + strlen (suffix_fmt) + 1);

  {
    char *p;
    if ((p = strrchr (argv[1], '/')) == 0)
      if ((p = strrchr (argv[1], '\\')) == 0)
	p = argv[1];
    if (*p == '/' || *p == '\\')
      p++;
    strcpy (outname, p);
  }

  if ((suffix = strrchr (outname, '.')) == 0)
    suffix = outname + strlen (outname);
  sprintf (suffix, suffix_fmt, board_channel + 1);
  //  printf ("writing %s", outname);
  (of = fopen (outname, "wb")) || DIE;
  fpos = ftello64 (f);
  fpos == sizeof (*ghdr) + sizeof (*bhdr) * ghdr->board_count || DIE;

  fpos = (fpos / (1 << 16) + 1) * (1 << 16);
  fseeko64 (f, fpos, SEEK_SET);

  {
    int *swipe_size, swipe_size_total, j, ecb;
    int skip, samples_per_100ms;
    short *buf, *obuf;
    time_t now, then;

    TMALLOC (swipe_size, ghdr->board_count);
    samples_per_100ms = sample_rate / 10;
    for (swipe_size_total = n = 0; n < ghdr->board_count; n++)
      swipe_size_total += (swipe_size[n] = enabled_channels[n] * 2 * samples_per_100ms);
    for (skip = n = 0; n < board; n++)
      skip += swipe_size[n];
    fseeko64 (f, skip, SEEK_CUR);
    TMALLOC (buf,  swipe_size[board] / sizeof *buf);
    TMALLOC (obuf,  samples_per_100ms);
    skip = swipe_size_total - swipe_size[board];
    ecb = enabled_channels[board];
    then = time (0);
    //    printf ("swipe_size[board]: %d\n", swipe_size[board]);
    //    printf ("samples_per_100ms: %d\n", samples_per_100ms);
    //    for (n = 0; n < 3356; n++) {
    for (n = 0; n < ghdr->swipe_count; n++) {
      if (1)
	if ((now = time(0)) > then) {
	  printf ("  %5d of %d swipes\r", n, ghdr->swipe_count);
	  fflush (stdout);
	  then = now;
	}
      fread64 (buf, swipe_size[board], 1, f) == 1 || DIE;
      for (i = chidx, j = 0; j < samples_per_100ms; j++, i += ecb)
	obuf[j] = buf[i];
      fwrite (obuf, sizeof *obuf, samples_per_100ms, of) == samples_per_100ms || DIE;
      fseeko64 (f, skip, SEEK_CUR);
    }
    printf ("  %5d of %d swipes\r", n, ghdr->swipe_count);
    printf ("\n");
  }  
  return 0;
}
