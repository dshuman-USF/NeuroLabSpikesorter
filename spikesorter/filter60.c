#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <error.h>
#include <errno.h>
#include <complex.h>
#include <fftw3.h>

static void
flatten (int freq, double fbin, complex *cplx)
{
  double widthfrac = .002;
  double fixfrac = .0001;
  int widthbins = nearbyint (widthfrac * freq / fbin);
  if (widthbins % 2 == 0) widthbins++;
  int centerbin = nearbyint (freq / fbin);
  int leftstart = nearbyint (centerbin - (widthbins - 1) / 2 - widthbins);
  int midstart = leftstart + widthbins;
  int rightstart = midstart + widthbins;
  int rightend = rightstart + widthbins;
  double sum = 0;
  for (int i = leftstart; i < midstart; i++)
    sum += cabs (cplx[i]);
  for (int i = rightstart; i < rightend; i++)
    sum += cabs (cplx[i]);
  double mean = sum / (widthbins * 2);
  int fixbins = nearbyint (fixfrac * freq / fbin);
  if (fixbins < 1) fixbins = 1;
  if (fixbins > widthbins) fixbins = widthbins;
  for (int i = midstart; i < rightstart; i += fixbins) {
    sum = 0;
    for (int j = 0; j < fixbins; j++)
      sum += cabs (cplx[i + j]);
    double fixmean = sum / fixbins;
    if (fixmean > mean) {
      double ratio = mean / fixmean;
      for (int j = 0; j < fixbins; j++)
        cplx[i + j] *= ratio;
    }
  }
}

int
main (int argc, char **argv)
{
  if (argc < 2) {
    printf ("usage: %s filename.chan\n", argv[0]);
    return 0;
  }
  FILE *f = fopen (argv[1], "r");
  if (f == NULL)
    error (1, errno, "can't open %s for read", argv[1]);
  fseek (f, 0, SEEK_END);
  int file_bytes = ftell (f);
  rewind (f);
  int fftlen = 640000;
  int fftbuflen = 2*(fftlen/2+1);

  double *real = (double*) fftw_malloc (sizeof(double) * fftbuflen);
  fftw_complex *cplx = (fftw_complex*) real;

  int buf_bytes = fftlen * sizeof (short);
  short *buf = malloc (buf_bytes);
  double fbin = 25000. / fftlen;
  printf ("creating plans\n");
  fftw_plan fwd = fftw_plan_dft_r2c_1d (fftlen, real, cplx, FFTW_ESTIMATE);
  fftw_plan rev = fftw_plan_dft_c2r_1d (fftlen, cplx, real, FFTW_ESTIMATE);
  printf ("plans done\n");

  char *ofilename;
  if (asprintf (&ofilename, "filtered_%s", argv[1]) == -1) exit (1);
  FILE *g = fopen (ofilename, "w");
  if (g == NULL)
    error (1, errno, "can't open %s for write", ofilename);

  int remaining_bytes = file_bytes;
  int bufs_written = 0;
  while (remaining_bytes > 0) {
    if (remaining_bytes < buf_bytes)
      fseek (f, -buf_bytes, SEEK_END);
    if (fread (buf, 1, buf_bytes, f));
    for (int i = 0; i < fftlen; i++)
      real[i] = buf[i];
    fftw_execute (fwd);
    for (int harmonic = 1; harmonic <= 26; harmonic++)
      flatten (harmonic * 60, fbin, cplx);
    fftw_execute (rev);
    for (int i = 0; i < fftlen; i++)
      buf[i] = nearbyint (real[i] / fftlen);
    if (remaining_bytes < buf_bytes) {
      fwrite (buf + (buf_bytes - remaining_bytes) / 2, 1, remaining_bytes, g);
      printf ("%d, %d\n", remaining_bytes, remaining_bytes + bufs_written * buf_bytes);
    }
    else {
      fwrite (buf, 1, buf_bytes, g);
      bufs_written++;
    }
    remaining_bytes -= buf_bytes;
  }
  fclose (g);
  fftw_destroy_plan (fwd);
  fftw_destroy_plan (rev);
  fftw_free(real);
  return 0;
}
