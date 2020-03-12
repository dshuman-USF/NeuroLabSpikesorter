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

/* ksone.c */


#include <math.h>
#include <stdlib.h>

static int cmp (const float *a, const float *b)
{
  return *a < *b ? -1 : (*a > *b); /* ascending */
}

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(float)(a),maxarg2=(float)(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


#define EPS1 0.001
#define EPS2 1.0e-8
float probks(float alam)
{
  int j;
  float a2,fac=2.0,sum=0.0,term,termbf=0.0;
  a2 = (float)(-2.0*alam*alam);
  for (j=1;j<=100;j++) {
    term=(float)(fac*exp(a2*j*j));
    sum += term;
    if (fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
    fac = -fac;
    termbf=(float)fabs(term);
  }
  return 1.0;
}

static void
ksone_nosort(float data[], unsigned long n, float (*func)(float), float *d, float *prob)
{
  float probks(float alam);
  unsigned long j;
  float dt,en,ff,fn,fo=0.0;
  en=(float)n;
  *d=0.0;
  for (j=1;j<=n;j++) {
    fn=(float)j/en;
    ff=(*func)(data[j]);
    dt=FMAX(fabs(fo-ff),fabs(fn-ff));
    if (dt > *d) *d=dt;
    fo=fn;
  }
  en=(float)sqrt(en);
  *prob=probks((float)(en+0.12+0.11/en)*(*d));
}

void ksone(float data[], unsigned long n, float (*func)(float), float *d, float *prob)
{
  qsort (data+1, n, sizeof *data, (int(*)(const void*,const void*))cmp);
  ksone_nosort(data, n, func, d, prob);
}
