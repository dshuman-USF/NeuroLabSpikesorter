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


/* lfs.c */


#if defined(_WIN32) || defined(__CYGWIN__)

#define __NO_MINGW_LFS 1

#ifndef __CYGWIN__
#include <windows.h>
#else
#include <w32api/windows.h>
#include <sys/cygwin.h>
#endif

#include "tmalloc.h"
#include "lfs.h"
#include <io.h>
#include <string.h>
#include <limits.h>

void *
fopen64 (const char *filename, const char *mode)
{
  HANDLE h;
  char win32_full_path[PATH_MAX];

  strcmp (mode, "rb") == 0 || DIE;

#ifdef __CYGWIN__
  cygwin32_conv_to_full_win32_path (filename, win32_full_path);
#else
  strcpy (win32_full_path, filename);
#endif

  h = CreateFile (win32_full_path, GENERIC_READ, FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
  if (h == INVALID_HANDLE_VALUE)
    return 0;

  return h;
}

int
fclose64 (void *h)
{
  if (CloseHandle (h) == 0) {
    fprintf (stderr, "win32 error %ld\n", GetLastError ());
    DIE;
  }
  return 0;
}

size_t
fread64 (void * buf, size_t size, size_t count, void *h)
{
  DWORD bytes_read;
  int items_read;

  if (ReadFile (h, buf, size * count, &bytes_read, 0) == 0) {
    LPVOID lpMsgBuf;
    FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		  FORMAT_MESSAGE_FROM_SYSTEM | 
		  FORMAT_MESSAGE_IGNORE_INSERTS,
		  NULL,
		  GetLastError(),
		  MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
		  (LPTSTR) &lpMsgBuf,
		  0,
		  NULL);
    fprintf (stderr, "win32 error reading %d bytes: %s\n", size * count, (char *)lpMsgBuf);
    LocalFree( lpMsgBuf );
    
    //    fprintf (stderr, "win32 error %ld\n", GetLastError ());
    DIE;
  }
  items_read = bytes_read / size;
  items_read * size == bytes_read || DIE;
  return items_read;
}

int
fseeko64 (void *h, off64_t offset, int whence)
{
  LARGE_INTEGER li;
  DWORD origin = 0;

  li.QuadPart = offset;

  if      (whence == SEEK_SET) origin = FILE_BEGIN;
  else if (whence == SEEK_CUR) origin = FILE_CURRENT;
  else if (whence == SEEK_END) origin = FILE_END;
  else DIE;

  if (SetFilePointer (h, li.LowPart, &li.HighPart, origin) == 0xFFFFFFFF && GetLastError () != NO_ERROR)
    return -1;
  else
    return 0;
}

off64_t
ftello64 (void *h)
{
  LARGE_INTEGER li;

  li.HighPart = 0;
  li.LowPart = SetFilePointer(h, 0, &li.HighPart, FILE_CURRENT);
  return li.QuadPart;
}

#endif
