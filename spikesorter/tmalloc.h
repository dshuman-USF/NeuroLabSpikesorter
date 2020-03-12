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

/* tmalloc.h */

#define __NO_MINGW_LFS 1

#include <stdio.h>
#include <stdlib.h>

extern int malloc_debug;
void die (void);

#define DIE (fprintf (stderr, "fatal error in %s at line %d \n", __FILE__, __LINE__), die(), 0)

#define TMALLOC(buf, n) ((buf = malloc ((n) * sizeof *(buf)))								\
			 && (malloc_debug == 0 || fprintf (stderr, "%d %d %ld\n", __LINE__, 0, (long)buf))) || DIE

#define TCALLOC(buf, n) ((buf = calloc ((n), sizeof *(buf)))								\
			 && (malloc_debug == 0 || fprintf (stderr, "%d %d %ld\n", __LINE__, 0, (long)buf))) || DIE

#define TREALLOC(buf, n) ((malloc_debug == 0 || fprintf (stderr, "%d %ld", __LINE__, (long)buf))	\
			  && (buf = realloc (buf, (n) * sizeof *(buf)))			\
			  && (malloc_debug == 0 || fprintf (stderr, " %ld\n", (long)buf))) || DIE
