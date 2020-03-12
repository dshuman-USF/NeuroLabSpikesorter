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

/* min.c */


#include "nde.h"
#include <time.h>

#define M 4
#define M2 4000000
#define N 2

static int print;

#define PARAFIT												\
do {													\
  int pull;												\
  (x1 > x0 && x2 > x1) || DIE;										\
  x10 = x1 - x0;											\
  x21 = x2 - x1;											\
  x20 = x2 - x0;											\
  m01 = (y1 - y0)/x10;											\
  m12 = (y2 - y1)/x21;											\
  m12 > m01 || DIE;											\
  m02 = (y2 - y0)/x20;											\
  den = (m12 - m01)/x20;										\
  para = 1;												\
  x = .5 * (x0 + x2 - m02/den);										\
  if (count < -10)											\
    printf ("x: %f\n", x);										\
  pull = 0;												\
  if (x2 - x > M * (x - x0)) {										\
    while (m * (x - x0) > (x2 - x) / 2)									\
      m /= N;												\
    x += m * (x - x0);											\
    pull = 1;												\
  }													\
  else if (x - x0 > M * (x2 - x)) {									\
    while (m * (x2 - x) > (x - x0) / 2)									\
      m /= N;												\
    x -= m * (x2 - x);											\
    pull = 1;												\
  }													\
  else if (1) {}											\
  else if (y2 - y1 > M2 * (y0 - y1)) {									\
    printf ("%.17f %.17f %.17f, %.17f %.17f %.17f\n", x0, x1, x2, y0, y1, y2);				\
    printf ("%.17f\n", x);										\
    x = (MAX (x, x1) + x2) / 2;										\
    printf ("%.17f\n", x);										\
    print = 1;												\
  }													\
  else if (y0 - y1 > M2 * (y2 - y1)) {									\
    printf ("%.17f %.17f %.17f, %.17f %.17f %.17f\n", x0, x1, x2, y0, y1, y2);				\
    printf ("%.17f\n", x);										\
    x = (MIN (x, x1) + x0) / 2;										\
    printf ("%.17f\n", x);										\
    print = 1;												\
  }													\
  if (print) printf ("%.17f %.17f %.17f, %.17f %.17f %.17f, %17f\n", x0, x1, x2, y0, y1, y2, x);	\
  y = func (x, p, 0);											\
  if (pull && ((((x > x1 && x21 < x10) || (x < x1 && x10 < x21)) && y > y1)				\
	       || (((x > x1 && x21 > x10) || (x < x1 && x10 > x21)) && y < y1)				\
	       || x == x1))										\
    m *= N;												\
  if (count < -10)											\
    printf ("pull: %d, y: %f, y1: %f, m: %f\n", pull, y, y1, m);					\
  hi = lo = x1;												\
} while (0)

static double fmin_y_val;

static union
{
  double d;
  struct
  {
    unsigned int wx;
    unsigned int wy;
  } w;
} inA, outA, inB, outB;

#define SCRAPE(val, a) (in##a.d = val, out##a.w.wx = in##a.w.wx, out##a.w.wy = in##a.w.wy, out##a.d)
#define GT(a,b) (SCRAPE (a,A)  > SCRAPE (b,B))
#define LT(a,b) (SCRAPE (a,A)  < SCRAPE (b,B))
#define EQ(a,b) (SCRAPE (a,A) == SCRAPE (b,B))

#define CHECK(new)				\
if (new < min || new > max)			\
{						\
      fmin_y_val = y1;				\
      *status = new < min ? -1 : 1;		\
      return x1;				\
}

double
fmin2 (double (*func)(double x, void *param, int debug), void *p, double startx, double radius, double min, double max, int *status)
{
  double x10, x21, x20, m01, m12, m02, dx;
  double x0, x1, x2, y0, y1, y2, x, y, den;
  int count, para, lastpara;
  static int max1, max2;
  static int invocations, evaluations;
  double m, hi, lo, xt, yt;
  static int debug;

  print = 0;
  radius > 0 || DIE;
  x0 = startx - radius;
  x1 = startx;
  x2 = startx + radius;
  y0 = func (x0, p, 0);
  y1 = func (x1, p, 0);
  y2 = func (x2, p, 0);
  
  debug = 0;
  count = 0;
  while (y1 >= y0 || y1 >= y2) {
    if (++count > max1) {
      max1 = count;
      if (0)
	if (max1 > 27)
	  printf ("\nmax1: %d max2: %d\n", max1, max2);
    }
    if (count > 100) {
      *status = 2;
      fmin_y_val = y1;
      return x1;
      debug = 1;
    }

    if (debug) printf ("line %d: %f %f %f => %f %f %f\n", __LINE__, x0, x1, x2, y0, y1, y2);

    if (y0 == y1 && y1 <= y2) {
      double midx = (x0+x1)/2;
      double midy = func (midx, p, 0);
      if (midy < y0) {
	CHECK (midx);
	x2 = x1, y2 = y1;
	x1 = midx, y1 = midy;
	if (debug) printf ("line %d: %f %f => %f %f\n", __LINE__, x1, x2, y1, y2);
	break;
      }
    }
    if (y0 > y1 && y1 == y2) {
      double midx = (x1+x2)/2;
      double midy = func (midx, p, 0);
      if (midy < y1) {
	CHECK (midx);
	x0 = x1, y0 = y1;
	x1 = midx, y1 = midy;
	if (debug) printf ("line %d: %f %f => %f %f\n", __LINE__, x0, x1, y0, y1);
	break;
      }
    }
    if (y0 < y2) {
      x2 = x1, y2 = y1;
      CHECK (x0);
      x1 = x0, y1 = y0;
      x0 = x1 - radius;
      y0 = func (x0, p, 0);
      if (debug) printf ("line %d: %f %f %f => %f %f %f\n", __LINE__, x0, x1, x2, y0, y1, y2);
    }
    else {
      x0 = x1, y0 = y1;
      CHECK (x2);
      x1 = x2, y1 = y2;
      x2 = x1 + radius;
      y2 = func (x2, p, 0);
      if (debug) printf ("line %d: %f %f %f => %f %f %f\n", __LINE__, x0, x1, x2, y0, y1, y2);
    }
//    printf ("%f: %f\n", x0, y0);
//    printf ("%f: %f\n", x1, y1);
//    printf ("%f: %f\n", x2, y2);

    if (debug) printf ("\n");
  }
  m = 1;
  {
    x10 = x1 - x0;
    x21 = x2 - x1;
    m01 = (y1 - y0)/x10;
    m12 = (y2 - y1)/x21;
    if (!(x1 > x0 && x2 > x1)
	|| !(m12 > m01)) {
      if (0)
	printf ("\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e %.16e\n",
		x0, x1, x2, y0, y1, y2, x10, x21, m01, m12);
      *status = 3;
      fmin_y_val = y1;
      return x1;
      exit (DIE);
    }
  }
  PARAFIT;

  count = 0;
  lastpara = 1;
  //  printf ("\ncount %d: %f %f %f, %f %f %f, %f %f (%d)\n", count, x0, x1, x2, y0, y1, y2, x, y, lastpara);

  //  if (invocations % 1000 == 0)
  if (0)
    if (invocations == 200000)
      printf ("%79s\r%s line %d: %d invocations, %f evaluations per invocation\n",
	      "", __FILE__, __LINE__, invocations, (double) evaluations / invocations);

  invocations++;
  evaluations++;

  while (1) {
    lastpara = 0;		/* debug */

    if (++count > max2) {
      max2 = count;
      //      if (max2 > 46301)
      if (0)
	printf ("\nmax1: %d max2: %d\n", max1, max2);
    }

    if (count > 1000) {
      printf ("\ncount %d: %.16f %.16f %.16f, %.15f %.15f %.15f, %.16f %.15f (%.0f)\n", count, x0, x1, x2, y0, y1, y2, x, y, m);
      if (count > 1050)
	func (x1, p, 1);
    }

    if (count < 0)
      printf ("\ncount %d: %.16f %.16f %.16f, %.15f %.15f %.15f, %.16f %.15f (%.0f)\n", count, x0, x1, x2, y0, y1, y2, x, y, m);

    para = 0;
    if (0)
    if (((SnapList *)p)->snap[0]->val[0] == 1577.618408203125
	&& ((SnapList *)p)->snap[1]->val[0] == 0.81232863664627075)
      printf ("count %d: %20.17f %20.17f %20.17f\n         %.17f %.17f %.17f, %.17f %.17f (%.0f)\n", count, x0, x1, x2, y0, y1, y2, x, y, m), fflush (stdout);
    if (print) printf ("%.17f %.17f %.17f, %.17f %.17f %.17f, %.17f %.17f\n", x0, x1, x2, y0, y1, y2, x, y);
    if (count > 100)
      debug = 1;
    if (0)
    if (count > 100) {
      printf ("%.17f %.17f\n", ((SnapList *)p)->snap[0]->val[0], ((SnapList *)p)->snap[1]->val[0]);
      exit (DIE);
    }
    if (LT(x, x1)) {
      dx = x1 - x;
      if (LT(y, y1)) {
	x2 = x1, y2 = y1;
	CHECK (x);
	x1 = x, y1 = y;
	if (lastpara && LT(dx, (x1 - x0) / 2))
	  y = func (x = x1 - dx, p, 0);
	else
	  PARAFIT;
      }
      else if (GT(y, y1)) {
	x0 = x, y0 = y;
	if (lastpara && LT(dx, (x2 - x1) / 2))
	  y = func (x = x1 + dx, p, 0);
	else
	  PARAFIT;
      }
      else {//y == y1
	xt = (x + x1) / 2;
	yt = func (xt, p, 0);
	if (LT (yt, y)) {
	  x0 = x;  y0 = y;
	  x2 = x1; y2 = y;
	  x1 = xt; y1 = yt;
	  PARAFIT;
	}
	else {
	  lo = x;
	  if (lo - x0 > x2 - hi)
	    x = (x0 + lo) / 2;
	  else
	    x = (hi + x2) / 2;
	  y = func (x, p, 0);
	}
      }
    }
    else if (GT(x, x1)) {
      dx = x - x1;
      if (LT(y, y1)) {
	x0 = x1, y0 = y1;
	CHECK (x);
	x1 = x, y1 = y;
	if (lastpara && LT(dx, (x2 - x1) / 2))
	  y = func (x = x1 + dx, p, 0);
	else
	  PARAFIT;
      }
      else if (GT(y, y1)) {
	x2 = x, y2 = y;
	if (lastpara && LT(dx, (x1 - x0) / 2))
	  y = func (x = x1 - dx, p, 0);
	else
	  PARAFIT;
      }
      else {//y == y1
	xt = (x + x1) / 2;
	yt = func (xt, p, 0);
	if (LT (yt, y)) {
	  x0 = x1; y0 = y;
	  x1 = xt; y1 = yt;
	  x2 = x;  y2 = y;
	  PARAFIT;
	}
	else {
	  hi = x;
	  if (lo - x0 > x2 - hi)
	    x = (x0 + lo) / 2;
	  else
	    x = (hi + x2) / 2;
	  y = func (x, p, 0);
	}
      }
    }
    else {
      if (y2 - y1 > y0 - y1) {
	if ((x = (x1 + x2) / 2) == x1 || x == x2)
	  x = nextafter (x2, x1);
      }
      else {
	if ((x = (x1 + x0) / 2) == x1 || x == x0)
	  x = nextafter (x0, x1);
      }
      y = func (x, p, 0);
      if (0)
	printf ("%.17f %.17f %.17f, %.17f %.17f %.17f, %.17f %.17f\n", x0, x1, x2, y0, y1, y2, x, y);
    }
    if (print) printf ("%.17f %.17f %.17f, %.17f %.17f %.17f, %.17f %.17f\n", x0, x1, x2, y0, y1, y2, x, y);
    //    printf ("\ncount %d: %f %f %f, %f %f %f, %f %f (%d)\n", count, x0, x1, x2, y0, y1, y2, x, y, lastpara);

    evaluations++;

    //    if (EQ(y0, y1) && EQ(y1, y2))
    if (y0 - y1 < .01 && y2 - y1 < .01) {
      *status = 0;
      if (y < y1) {
	CHECK (x);
	fmin_y_val = y;
	return x;
      }
      else {
	fmin_y_val = y1;
	return x1;
      }
    }
    lastpara = para;
  }
}

double
fmin1 (double (*func)(double x, void *param, int debug), void *p, double startx, double radius)
{
  double x10, x21, x20, m01, m12, m02, dx;
  double x0, x1, x2, y0, y1, y2, x, y, den;
  int para, count = 0;
  double m = 1, hi, lo;

  if (0) printf ("%d %g", para, hi);

  radius > 0 || DIE;
  x0 = startx - radius;
  x1 = startx;
  x2 = startx + radius;
  y0 = func (x0, p, 0);
  y1 = func (x1, p, 0);
  y2 = func (x2, p, 0);
  
  while (y1 >= y0 || y1 >= y2) {
    static int debug;

    if (debug) printf ("line %d: %f %f %f => %f %f %f\n", __LINE__, x0, x1, x2, y0, y1, y2);

    if (y0 == y1 && y1 <= y2) {
      double midx = (x0+x1)/2;
      double midy = func (midx, p, 0);
      if (debug) printf ("line %d: %f %f\n", __LINE__, midx, midy);
      if (midy < y0) {
	x2 = x1, y2 = y1;
	x1 = midx, y1 = midy;
	if (debug) printf ("line %d: %f %f => %f %f\n", __LINE__, x1, x2, y1, y2);
	break;
      }
    }
    if (y0 > y1 && y1 == y2) {
      double midx = (x1+x2)/2;
      double midy = func (midx, p, 0);
      if (debug) printf ("line %d: %f %f\n", __LINE__, midx, midy);
      if (midy < y1) {
	x0 = x1, y0 = y1;
	x1 = midx, y1 = midy;
	if (debug) printf ("line %d: %f %f => %f %f\n", __LINE__, x0, x1, y0, y1);
	break;
      }
    }
    if (y0 < y2) {
      x2 = x1, y2 = y1;
      x1 = x0, y1 = y0;
      x0 = x1 - radius;
      y0 = func (x0, p, 0);
      if (debug) printf ("line %d: %f %f %f => %f %f %f\n", __LINE__, x0, x1, x2, y0, y1, y2);
    }
    else {
      x0 = x1, y0 = y1;
      x1 = x2, y1 = y2;
      x2 = x1 + radius;
      y2 = func (x2, p, 0);
      if (debug) printf ("line %d: %f %f %f => %f %f %f\n", __LINE__, x0, x1, x2, y0, y1, y2);
    }
    if (debug) printf ("\n");
  }

  (x1 > x0 && x2 > x1) || DIE;
  x10 = x1 - x0;
  x21 = x2 - x1;
  x20 = x2 - x0;
  m01 = (y1 - y0)/x10;
  m12 = (y2 - y1)/x21;
  if (m12 <= m01) {
    double mdistance_at_shift_partsnap_debug (double shift, void *p);

    printf ("line %d: %f %f %f => %f %f %f\n", __LINE__, x0, x1, x2, y0, y1, y2),
      printf ("line %d: %f %f %f => %f %f\n", __LINE__, x10, x21, x20, m01, m12);
    (void)mdistance_at_shift_partsnap_debug (x0, p);
    (void)mdistance_at_shift_partsnap_debug (x1, p);
    (void)mdistance_at_shift_partsnap_debug (x2, p);
  }

  m12 > m01 || DIE;
  m02 = (y2 - y0)/x20;
  den = (m12 - m01)/x20;
  x = .5 * (x0 + x2 - m02/den);
  y = func (x, p, 0);

  //  PARAFIT;
  while (1) {
    if (LT(x, x1)) {
      dx = x1 - x;
      if (LT(y, y1)) {
	x2 = x1, y2 = y1;
	x1 = x, y1 = y;
	if (LT(dx, (x1 - x0) / 2))
	  y = func (x = x1 - dx, p, 0);
	else
	  PARAFIT;
      }
      else if (GT(y, y1)) {
	x0 = x, y0 = y;
	if (LT(dx, (x2 - x1) / 2))
	  y = func (x = x1 + dx, p, 0);
	else
	  PARAFIT;
      }
      else {//y == y1
	x0 = x;
	x2 = x1;
	y0 = y2 = y;
	y1 = func (x1 = (x0 + x2) / 2, p, 0);
	if (GT(y1, y0)) {
	  double mdistance_at_shift_partsnap_debug (double shift, void *p);
	  printf ("line %d: %f %f %f => %f %f %f\n", __LINE__, x0, x1, x2, y0, y1, y2);
	  (void)mdistance_at_shift_partsnap_debug (x0, p);
	  (void)mdistance_at_shift_partsnap_debug (x1, p);
	  (void)mdistance_at_shift_partsnap_debug (x2, p);
	}
	GT(y1, y0) && DIE;
      }
    }
    else if (GT(x, x1)) {
      dx = x - x1;
      if (LT(y, y1)) {
	x0 = x1, y0 = y1;
	x1 = x, y1 = y;
	if (LT(dx, (x2 - x1) / 2))
	  y = func (x = x1 + dx, p, 0);
	else
	  PARAFIT;
      }
      else if (GT(y, y1)) {
	x2 = x, y2 = y;
	if (LT(dx, (x1 - x0) / 2))
	  y = func (x = x1 - dx, p, 0);
	else
	  PARAFIT;
      }
      else {//y == y1
	x0 = x1;
	x2 = x;
	y0 = y2 = y;
	y1 = func (x1 = (x0 + x2) / 2, p, 0);
	GT(y1, y0) && DIE;
      }
    }
    else {
      x = x1 + sqrt ((nextafter (y1, y2) - y1) / den);
      (void)(LT(x, (x1 + x2) / 2) || (x = (x1 + x2) / 2));
      y = func (x, p, 0);
    }
    //    if (EQ(y0, y1) && EQ(y1, y2))
    if (y0 - y1 < .001 && y2 - y1 < .001) {
      fmin_y_val = y1;
      return x1;
    }
  }
}

double
fmin_y (void)
{
  return fmin_y_val;
}
