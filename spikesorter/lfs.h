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

/* lfs.h */


#if defined(_WIN32) || defined(__CYGWIN__)

#include <stddef.h>

typedef long long off64_t;

void * fopen64 (const char *filename, const char *mode);

int fclose64 (void *h);

size_t fread64 (void * buf, size_t size, size_t count, void *h);

int fseeko64 (void *h, off64_t offset, int whence);

off64_t ftello64 (void *h);

#else

#define _LARGEFILE64_SOURCE
#define fread64 fread

#endif
