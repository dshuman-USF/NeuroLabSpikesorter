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

/* die.c */


#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#ifdef _WIN32
#define sleep _sleep
#else
#include <unistd.h>
#endif

static void
die_break (void)
{
  exit (EXIT_FAILURE);
}

void die (void)
{
  if (errno)
    perror ("last system error");
  fflush (stdout);
  fflush (stderr);
  sleep (1);
  die_break ();
}
