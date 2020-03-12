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

/* dmx_th.c */


#include "nde.h"
#include <pthread.h>
#include <sched.h>
#include <semaphore.h>
#include <string.h>
#include <ctype.h>
#include "dmx.h"

#ifdef S_SPLINT_S
typedef long long off64_t;
int fseeko64 (void *stream, off64_t offset, int whence);
void *fopen64 (const char *filename, const char *mode);
off64_t ftello64 (void *stream);
#endif

extern time_t start_time;

typedef struct {
  short *bufptr;
  size_t bufsize, bufidx, total_samples, samples_read;
  int channel;
  int chidx;
  int swipe_size_board;
  int enabled_channels_board;
  int swipe_count;
  long skip;
  FILE *file;
  sem_t empty;
  sem_t full;
  pthread_t thread;
} DMX;

static void *
dmx_background (void *arg)
{
  int count, ibuf_samples, n, i, j, ecb;
  short *ibuf, *obuf;
  DMX *dmx = arg;
  time_t now, then;
  int sec;

  //  printf ("line %d: dmx_background\n", __LINE__);
  ecb = dmx->enabled_channels_board;
  ibuf_samples = dmx->swipe_size_board / (int)sizeof *ibuf;
  TMALLOC (ibuf, ibuf_samples);
  then = time (0);
  sec = 0;
  //  printf ("\nline %d: %5ld waiting on empty\n", __LINE__, time (0) - start_time);
  sem_wait (&dmx->empty);
  //  printf ("\nline %d: %5ld filling buffer\n", __LINE__, time (0) - start_time);
  for (n = 0; n < dmx->swipe_count; n++) {
    dmx->swipe_size_board >= 0 || DIE;
    fread64 (ibuf, (size_t)dmx->swipe_size_board, 1, dmx->file) == 1 || DIE;
    //    printf ("%d: %5d of %d\n", sec++, n, dmx->swipe_count);
    if (0)
    if ((now = time(0)) > then) {
      printf ("%d: %5d of %d\n", sec++, n, dmx->swipe_count);
      then = now;
    }
    for (i = dmx->chidx; i < ibuf_samples; ) {
      count = MIN ((ibuf_samples - (i - dmx->chidx)) / ecb, (int)(dmx->bufsize - dmx->bufidx));
      obuf = dmx->bufptr + dmx->bufidx;
      for (j = 0; j < count; j++, i += ecb)
	obuf[j] = ibuf[i];
      dmx->bufidx += count;
      //      printf ("line %d %d %d %d %d %d\n", __LINE__, dmx->bufidx, dmx->bufsize, count, ibuf_samples, i);
      if ((dmx->bufidx == dmx->bufsize) || (n + 1 == dmx->swipe_count && i == ibuf_samples)) {
	sem_post (&dmx->full);
	//	printf ("\nline %d: %5ld waiting on empty\n", __LINE__, time (0) - start_time);
	sem_wait (&dmx->empty);
	//	printf ("\nline %d: %5ld filling buffer\n", __LINE__, time (0) - start_time);
	sec = 0;
      }
    }
    if (n + 1 < dmx->swipe_count)
      /*@+voidabstract@*/ fseeko64 (dmx->file, dmx->skip, SEEK_CUR); /*@=voidabstract@*/
  }
  fclose64 (dmx->file) == 0 || DIE;
  while (dmx->bufidx > 0) {
    sem_post (&dmx->full);
    //    printf ("line %d: waiting on empty\n", __LINE__);
    sem_wait (&dmx->empty);
    //    printf ("line %d: \"filling\" buffer\n", __LINE__);
  } 
  free (ibuf);
  return 0;
}

void *
dmx_open (char *file_name, size_t bufsize)
{
  char *fname, *dot, *chan_txt;
  int swipe_size_total, skip, samples_per_swipe_time, bytes_per_swipe_time;
  int n, i, j, sample_rate, chan, board, channel, chidx = -1;
  off64_t fpos;
  DMX *dmx;
  static GlobalHeader *ghdr;
  static BoardHeader *bhdr;
  static int *enabled_channels;
  static int *swipe_size;
  FILE *f;

  (fname = strdup (file_name)) || DIE;
  if (!(dot = strrchr (fname, '.')))
    dot = fname + strlen (fname);
  if (dot - fname < 3)
    goto fail;
  chan_txt = dot - 2;
  if (!isdigit ((int)chan_txt[0]) || !isdigit ((int)chan_txt[1]))
    goto fail;
  if ((chan = atoi (chan_txt)) == 0)
    goto fail;
  chan--;
  memmove (chan_txt, dot, strlen (dot) + 1);
  f = fopen64 (fname, "rb");
  free (fname);
  if (f == 0)
    goto fail;

  TREALLOC (ghdr, 1);
  if (fread64 (ghdr, sizeof *ghdr, 1, f) != 1
      || strncmp (ghdr->fid, "DATAMAX ", 8) != 0
      || chan >= ghdr->channel_count) {
    fclose (f);
    return 0;
  }
  board = chan / 8;
  channel = chan % 8;
  ghdr->board_count >= 0 || DIE;
  TREALLOC (bhdr, ghdr->board_count);
  TREALLOC (enabled_channels, ghdr->board_count);
  fread64 (bhdr, sizeof *bhdr, (size_t)ghdr->board_count, f) == (size_t)ghdr->board_count || DIE;
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
  }
  chidx >= 0 || DIE;

  /*@+voidabstract@*/ 
  fpos = ftello64 (f);
  /*@=voidabstract@*/

  fpos == (off64_t)(sizeof (*ghdr) + sizeof (*bhdr) * ghdr->board_count) || DIE;
  fpos = (fpos / (1 << 16) + 1) * (1 << 16);
  /*@+voidabstract@*/ 
  fseeko64 (f, fpos, SEEK_SET);
  /*@=voidabstract@*/

  TREALLOC (swipe_size, ghdr->board_count);
  sample_rate >= 0 || DIE;
  samples_per_swipe_time = sample_rate / 10;
  bytes_per_swipe_time = samples_per_swipe_time * (int)sizeof (short);
  for (swipe_size_total = n = 0; n < ghdr->board_count; n++)
    swipe_size_total += (swipe_size[n] = enabled_channels[n] * bytes_per_swipe_time);
  for (skip = n = 0; n < board; n++)
    skip += swipe_size[n];
  /*@+voidabstract@*/ 
  fseeko64 (f, skip, SEEK_CUR);
  /*@=voidabstract@*/
  TCALLOC (dmx, 1);
  dmx->chidx = chidx;
  dmx->skip = swipe_size_total - swipe_size[board];
  dmx->bufsize = bufsize;
  TMALLOC (dmx->bufptr, bufsize);
  dmx->file = f;
  dmx->swipe_size_board = swipe_size[board];
  dmx->enabled_channels_board = enabled_channels[board];
  ghdr->swipe_count >= 0 || DIE;
  dmx->swipe_count = ghdr->swipe_count;
  dmx->total_samples = (size_t)(dmx->swipe_count * samples_per_swipe_time);
  dmx->samples_read = 0;
  sem_init (&dmx->empty, 0, 1);
  sem_init (&dmx->full, 0, 0);

  {
    struct sched_param param;
    int policy;

    pthread_create (&dmx->thread, 0, dmx_background, dmx) == 0 || DIE;
    param.sched_priority = sched_get_priority_max (SCHED_FIFO);
    pthread_setschedparam (dmx->thread, SCHED_FIFO, &param);

    if (0) {
      pthread_getschedparam (dmx->thread, &policy, &param);
      printf ("policy: %d, priority: %d\n", policy, param.sched_priority);
      printf ("policy %d: min %d, max %d\n", SCHED_FIFO, sched_get_priority_min (SCHED_FIFO), sched_get_priority_max (SCHED_FIFO));
      printf ("policy %d: min %d, max %d\n", SCHED_RR, sched_get_priority_min (SCHED_RR), sched_get_priority_max (SCHED_RR));
      printf ("policy %d: min %d, max %d\n", SCHED_OTHER, sched_get_priority_min (SCHED_OTHER), sched_get_priority_max (SCHED_OTHER));
    }
  }
  return dmx;
 fail:
  return 0;
}

size_t
dmx_read (void *buf, size_t size, size_t count, void *p)
{
  size_t samplecount, samples_left;
  DMX *dmx = p;
  static int herecount;

  if (count == 0 || dmx->samples_read == dmx->total_samples)
    return 0;
  samplecount = (count * size) / sizeof *dmx->bufptr;
  //  printf ("	line %d: %5ld waiting on full\n", __LINE__, time (0) - start_time);
  sem_wait (&dmx->full);
  herecount++;
  //  printf ("	line %d: %5ld emptying buffer %d\n", __LINE__, time (0) - start_time, herecount);
  if (dmx->bufidx < samplecount)
    samplecount = dmx->bufidx;
  memcpy (buf, dmx->bufptr, samplecount * sizeof *dmx->bufptr);
  samples_left = dmx->bufidx - samplecount;
  memmove (dmx->bufptr, dmx->bufptr + samplecount, samples_left * sizeof *dmx->bufptr);
  dmx->bufidx = samples_left;
  dmx->samples_read += samplecount;
  sem_post (&dmx->empty);
  //  printf ("	line %d: %5ld doing spikesort\n", __LINE__, time (0) - start_time);
  return (samplecount * sizeof *dmx->bufptr) / size;
}
