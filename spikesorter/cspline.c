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


/* cspline.c */


#include "nde.h"

//#define free(x) myfree(__LINE__, x)

#define PARABLEND(x, y)				\
{						\
  double x2 = x * x;				\
  double x3 = x2 * x;				\
						\
  return (( (   x3 -   x2)     * y[3]		\
	  + (-3*x3 + 4*x2 + x) * y[2]		\
	  + ( 3*x3 - 5*x2 + 2) * y[1]		\
	  + (  -x3 + 2*x2 - x) * y[0])		\
	  / 2);					\
}

#define CSPLINE_EVAL						\
{								\
  if (x == (double)xi && x >= 0 && x < (double)count)		\
    return (double)y[xi];					\
  else if (x < 1) {						\
    yt[0] = (double)y[0] - (y[1] - y[0]);			\
    for (n = 0; n < 3; n++)					\
      yt[n+1] = (double)y[n];					\
    PARABLEND (x, yt);						\
  }								\
  else if (x > (double)count - 2) {				\
    for (n = 0; n < 3; n++)					\
      yt[n] = (double)y[count-3+n];				\
    yt[3] = (double)y[count-1] + (y[count-1] - y[count-2]);	\
    x -= count - 2;						\
    PARABLEND (x, yt);						\
  }								\
  else {							\
    x -= xi;							\
    y += xi - 1;						\
    PARABLEND (x, y);						\
  }								\
}

static double
cspline_eval (const CSpline *csp, double x)
{
  int count = csp->size;
  int xi = (int)floor (x);
  double yt[4];
  int n;

  if (csp->y_float) {
    float *y = csp->y_float;
    CSPLINE_EVAL;
  }
  else {
    short *y = csp->y_short;
    CSPLINE_EVAL;
  }
  exit (DIE); return 0;
}

CSpline *
cspline_new (float *y_float, short *y_short, int size)
{
  CSpline *csp;

  TCALLOC (csp, 1);
  csp->y_float = y_float;
  csp->y_short = y_short;
  csp->size = size;
  return csp;
}

double
ss_cspline_eval (const gsl_spline *spline, double x, gsl_interp_accel *a)
{
  if (spline->interp)
    return gsl_interp_eval (spline->interp, spline->x, spline->y, x, a);
  return cspline_eval ((CSpline *)spline, x);
}

void
/*@-incondefs@*/
ss_spline_free (/*@only@*/ gsl_spline * spline)
/*@=incondefs@*/
{
  if (spline->interp) {
    gsl_interp_free (spline->interp);
    free (spline->x);
    free (spline->y);
    free (spline);
  }
  else {
    CSpline *csp = (CSpline *)spline;
    free (csp->y_float);
    free (csp->c);
    free (csp);
  }
}
