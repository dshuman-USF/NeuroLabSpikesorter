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


/* chan2hdt.c */


#include <config.h>
#define _FILE_OFFSET_BITS 64
#include "nde.h"
#include <inttypes.h>
#include <limits.h>

int malloc_debug;

typedef struct
{
  char *filename;
  FILE *file;
  int16_t *buf;
  short min, max;
  int bufidx;
  int samples_left;
} Chan;

static inline void
swapshort (char *bp)
{
  char c;
  char *sp;

  c = *bp;
  *bp = *(sp = bp + 1);
  *sp = c;
}
                                                                                
static inline void
swaplong (char *bp)
{
  char c;
  char *sp;
                                                                                
  c = *(sp = bp + 3);
  *sp = *bp;
  *bp++ = c;
  c = *(sp = bp + 1);
  *sp = *bp;
  *bp = c;
}

#if defined i386 && defined __GNUC__
#define SWAP16(x) __asm__("xchgb %b0,%h0" : "=q" (x) :  "0" (x))
#else
#define SWAP16(x) swapshort ((char *)&x)
#endif

#define SWAP32(x) swaplong ((char *)&x)

#ifdef WORDS_BIGENDIAN
# define l_e_to_host_16(x) SWAP16(x)
# define l_e_to_host_32(x) SWAP32(x)
# define host_to_b_e_16(x)
# define host_to_b_e_32(x)
#else
# define l_e_to_host_16(x)
# define l_e_to_host_32(x)
# define host_to_b_e_16(x) SWAP16(x)
# define host_to_b_e_32(x) SWAP32(x)
#endif

static int64_t ticks_per_second = 2000;
static int64_t samples_per_second = SF;

#define BUFSAMPLE_TARGET (10*1024*1024)
int bufsamples;
static int min_samples;

static char *hdt_filename;
static FILE *hdt_file;
static int hdt_hdr_size;
static int16_t *hdt_buf;

static uint16_t sample_rate;
static uint16_t channel_count;
static int16_t min;
static int16_t max;
static uint32_t block_count;
static uint32_t first;
static uint32_t last;

static int
get_chan (Chan *chan, int chancnt, int chanidx)
{
  int count = bufsamples / chancnt;
  int min_samples_per_chan = min_samples / chancnt;
  int n;
  int16_t *ibuf;
  int16_t *obufp;

  chan->bufidx + chan->samples_left > bufsamples && DIE;
  if (feof (chan->file) && chan->samples_left == 0)
    return 0;
  chan->bufidx + chan->samples_left == bufsamples || feof (chan->file) || ftell (chan->file) == 0 || DIE;
  chan->bufidx % count == 0 || DIE;
  chan->samples_left % min_samples_per_chan == 0 || DIE;
  if (chan->samples_left == 0) {
    int readcnt = fread (chan->buf, sizeof *chan->buf, bufsamples, chan->file);
    if (readcnt < bufsamples)
      feof (chan->file) || DIE;
    chan->samples_left = readcnt / min_samples_per_chan * min_samples_per_chan;
    chan->bufidx = 0;
  }
  if (chan->samples_left < count)
    count = chan->samples_left;
  ibuf = chan->buf + chan->bufidx;
  obufp = hdt_buf + chanidx;
  for (n = 0; n < count; n++) {
    int16_t val = ibuf[n];
    l_e_to_host_16 (val);
    if (val > chan->max) chan->max = val;
    if (val < chan->min) chan->min = val;
    host_to_b_e_16 (val);
    *obufp = val;
    obufp += chancnt;
  }
  chan->samples_left -= count;
  chan->bufidx += count;
  return count;
}

static int
fill_hdt_buf (Chan *chan, int chancnt)
{
  int n, chan_readcnt, chan0_readcnt;
  
  chan0_readcnt = get_chan (chan, chancnt, 0);
  for (n = 1; n < chancnt; n++)
    if ((chan_readcnt = get_chan (chan + n, chancnt, n)) != chan0_readcnt) {
      fprintf (stderr, "FATAL ERROR: %s is a different size than the preceding .chan files\n", chan[n].filename);
      exit (1);
    }
  return chan0_readcnt * chancnt;
}

static FILE *
file_open (char *filename, char *mode)
{
  FILE *f;
  if ((f = fopen (filename, mode)) == 0) {
    char *emsg;
    if (asprintf (&emsg, "Can't open %s %s", filename, mode) == -1) exit (1);
    perror (emsg);
    exit (1);
  }
  return f;
}

static int
gcd (int a, int b)
{
  while (a != b)
    if (a < b) b -= a;
    else       a -= b;
  return a;
}

#define WEIRDNESS(x) x = (x << 1 & 0xffff0000) | (x & 0x7fff); host_to_b_e_32 (x)

int
main (int argc, char **argv)
{
  int chancnt, n, samples_written, readcnt;
  Chan *chan;

  if (argc < 3) {
    printf ("this program is part of %s\n", PACKAGE_STRING);
    printf ("usage: %s chan_file_prefix chan# ...\n", argv[0]);
    exit (0);
  }
  chancnt = argc - 2;

  TCALLOC (chan, chancnt);
  {
    int gcdval = gcd (ticks_per_second, samples_per_second);

    min_samples = samples_per_second / gcdval * chancnt;
    bufsamples = (BUFSAMPLE_TARGET / min_samples + 1) * min_samples;
    bufsamples % chancnt == 0 || DIE;
  }
  for (n = 0; n < chancnt; n++) {
    int channum = atoi (argv[n + 2]);

    if (asprintf (&chan[n].filename, "%s_%02d.chan", argv[1], channum) == -1) exit (1);
    chan[n].file = file_open (chan[n].filename, "rb");
    TMALLOC (chan[n].buf, bufsamples);
    chan[n].min = SHRT_MAX;
    chan[n].max = SHRT_MIN;
  }
  if (asprintf (&hdt_filename, "%s.hdt", argv[1]) == -1) exit (1);
  hdt_file = file_open (hdt_filename, "wb");
  TMALLOC (hdt_buf, bufsamples);

  hdt_hdr_size = (sizeof sample_rate
		  + sizeof channel_count
		  + chancnt * (sizeof min + sizeof max)
		  + sizeof block_count
		  + 1 * (sizeof first + sizeof last));
  bufsamples * sizeof *hdt_buf > hdt_hdr_size || DIE;
  fwrite (hdt_buf, 1, hdt_hdr_size, hdt_file) == hdt_hdr_size || DIE;

  samples_written = 0;
  while ((readcnt = fill_hdt_buf (chan, chancnt)) > 0) {
    fwrite (hdt_buf, sizeof *hdt_buf, readcnt, hdt_file) == readcnt || DIE;
    samples_written += readcnt;
  }
  samples_written % chancnt == 0 || DIE;
  {
    int64_t samples_per_channel = samples_written / chancnt;
    int64_t ticks_per_channel = samples_per_channel * ticks_per_second / samples_per_second;
    samples_per_channel * ticks_per_second % samples_per_second == 0 || DIE;
    ticks_per_channel <= 0xffffffff || DIE;
    last = ticks_per_channel;
  }
  rewind (hdt_file);

  sample_rate   = (int)SF; host_to_b_e_16 (sample_rate  ); fwrite (&sample_rate  , sizeof sample_rate  , 1, hdt_file) == 1 || DIE;
  channel_count = chancnt; host_to_b_e_16 (channel_count); fwrite (&channel_count, sizeof channel_count, 1, hdt_file) == 1 || DIE;
  for (n = 0; n < chancnt; n++) {
    min = chan[n].min; host_to_b_e_16 (min); fwrite (&min, sizeof min, 1, hdt_file) == 1 || DIE;
    max = chan[n].max; host_to_b_e_16 (max); fwrite (&max, sizeof max, 1, hdt_file) == 1 || DIE;
  }
  block_count   = 1      ; WEIRDNESS      (block_count  ); fwrite (&block_count  , sizeof block_count  , 1, hdt_file) == 1 || DIE;
  first         = 0      ; WEIRDNESS      (first        ); fwrite (&first        , sizeof first        , 1, hdt_file) == 1 || DIE;
  last          = last   ; WEIRDNESS      (last         ); fwrite (&last         , sizeof last         , 1, hdt_file) == 1 || DIE;
  ftell (hdt_file) == hdt_hdr_size || DIE;
  return 0;
}
