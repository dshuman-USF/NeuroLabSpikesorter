AC_DEFUN([SS_NO_CYGWIN],
[case $host_os in
  *cygwin*) NO_CYGWIN=-mno-cygwin ;;
  *) NO_CYGWIN= ;;
esac
AC_SUBST(NO_CYGWIN)])

AC_DEFUN([SS_TCLTK_TEST],
[AC_MSG_CHECKING([for tcltk library >= 8.4.4])

AC_DEFUN([SS_CHECK_TCLTK],
ss_save_LIBS=$LIBS
LIBS="$LIBS -ltcl8.5 -ltk8.5 $X_LIBS -lX11 $X_EXTRA_LIBS"
[AC_RUN_IFELSE([AC_LANG_SOURCE([
#include <tcl.h>
#include <tk.h>
#define MAJOR 8
#define MINOR 4
#define PATCH 4
int
Tcl_AppInit(Tcl_Interp *interp)
{
  return TCL_OK;
}
int
main (int argc, char **argv)
{
  int major, minor, patch, type;
  Tcl_GetVersion (&major, &minor, &patch, &type);

  if (type != TCL_FINAL_RELEASE)
    return 1;
  if (major < MAJOR || (major == MAJOR && minor < MINOR) || (major == MAJOR && minor == MINOR && patch < PATCH))
    return 1;

  if (TK_MAJOR_VERSION < MAJOR
      || (TK_MAJOR_VERSION == MAJOR && TK_MINOR_VERSION < MINOR)
      || (TK_MAJOR_VERSION == MAJOR && TK_MINOR_VERSION == MINOR && TK_RELEASE_SERIAL < PATCH))
    return 1;

  Tk_GetNumMainWindows ();

  return 0;
}
])], [ss_pass=yes], [ss_pass=no])
LIBS=$ss_save_LIBS
])

ss_old_LDFLAGS=$LDFLAGS
ss_old_CPPFLAGS=$CPPFLAGS
ss_done=no
ss_pass=no

if test $ss_done = no && test -f /usr/local/spikesorter/lib/libtk8.4.a ; then
    LDFLAGS="-L/usr/local/spikesorter/lib $LDFLAGS"
    CPPFLAGS="-I/usr/local/spikesorter/include $CPPFLAGS"
    SS_CHECK_TCLTK
fi

if test "$ss_pass" = yes ; then
    TCLTK_LDFLAGS=-L/usr/local/spikesorter/lib
    TCLTK_CPPFLAGS=-I/usr/local/spikesorter/include
    AC_MSG_RESULT([yes, in /usr/local/spikesorter])
    ss_done=yes
else
    LDFLAGS=$ss_old_LDFLAGS
    CPPFLAGS="-I/usr/include/tcl8.5 $ss_old_CPPFLAGS"
    SS_CHECK_TCLTK
fi

if test $ss_done = no && test $ss_pass = yes ; then
    unset TCLTK_LDFLAGS
    TCLTK_CPPFLAGS=-I/usr/include/tcl8.5
    AC_MSG_RESULT([yes])
    ss_done=yes
fi

if test $ss_done = no ; then
    LDFLAGS="-L/usr/local/lib $ss_old_LDFLAGS"
    CPPFLAGS="-I/usr/local/include $ss_old_CPPFLAGS"
    SS_CHECK_TCLTK
fi

if test $ss_done = no && test $ss_pass = yes ; then
    TCLTK_LDFLAGS=-L/usr/local/lib
    TCLTK_CPPFLAGS=-I/usr/local/include
    AC_MSG_RESULT([yes, in /usr/local])
    ss_done=yes
fi

if test $ss_pass = no && test -f tk8.4.4-src.tar.gz ; then
    TCLTK_LDFLAGS=-L/usr/local/spikesorter/lib
    TCLTK_CPPFLAGS=-I/usr/local/spikesorter/include
    TCLTK_DEP=tk
    rm -f tk
    AC_MSG_RESULT([no, will install in /usr/local/spikesorter])
    ss_done=yes
fi

if test $ss_done = no ; then
    AC_MSG_FAILURE([no tcltk library found])
fi

LDFLAGS=$ss_old_LDFLAGS
CPPFLAGS=$ss_old_CPPFLAGS

AC_SUBST(TCLTK_DEP)
AC_SUBST(TCLTK_LDFLAGS)
AC_SUBST(TCLTK_CPPFLAGS)

])


AC_DEFUN([SS_FFTW_TEST],
[AC_MSG_CHECKING([for fftw library >= 2.1.3])

AC_DEFUN([SS_CHECK_FFTW],
ss_save_LIBS=$LIBS
LIBS="$LIBS -lrfftw -lfftw"
[AC_RUN_IFELSE([AC_LANG_SOURCE([
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fftw.h>
#define MAJOR 2
#define MINOR 1
#define PATCH 3
int
main (void)
{
  int major, minor, patch = 0;
  char *p, *q;
  if (strncmp (fftw_version, "FFTW V", 6) != 0)
    return 1;
  major = strtol ((p = (char *)(fftw_version + 6)), &q, 10);
  if (q == p || *q++ != '.')
    return 1;
  p = q;
  minor = strtol (p, &q, 10);
  if (q == p)
    return 1;
  if (*q++ == '.')
    patch = strtol (q, 0, 10);
  if (major < MAJOR || (major == MAJOR && minor < MINOR) || (major == MAJOR && minor == MINOR && patch < PATCH))
    return 1;
  return 0;
}
])], [ss_pass=yes], [ss_pass=no])
LIBS=$ss_save_LIBS
])

ss_old_LDFLAGS=$LDFLAGS
ss_old_CPPFLAGS=$CPPFLAGS
ss_done=no
ss_pass=no

if test -f /usr/local/spikesorter/lib/libfftw.a ; then
    LDFLAGS="-L/usr/local/spikesorter/lib $LDFLAGS"
    CPPFLAGS="-I/usr/local/spikesorter/include $CPPFLAGS"
    SS_CHECK_FFTW
fi

if test "$ss_pass" = yes ; then
    FFTW_LDFLAGS=-L/usr/local/spikesorter/lib
    FFTW_CPPFLAGS=-I/usr/local/spikesorter/include
    AC_MSG_RESULT([yes, in /usr/local/spikesorter])
    ss_done=yes
else
    LDFLAGS=$ss_old_LDFLAGS
    CPPFLAGS=$ss_old_CPPFLAGS
    SS_CHECK_FFTW
fi

if test $ss_done = no && test $ss_pass = yes ; then
    unset FFTW_LDFLAGS
    unset FFTW_CPPFLAGS
    AC_MSG_RESULT([yes])
    ss_done=yes
fi

if test $ss_done = no ; then
    LDFLAGS="-L/usr/local/lib $LDFLAGS"
    CPPFLAGS="-I/usr/local/include $CPPFLAGS"
    SS_CHECK_FFTW
fi

if test $ss_done = no && test $ss_pass = yes ; then
    FFTW_LDFLAGS=-L/usr/local/lib
    FFTW_CPPFLAGS=-I/usr/local/include
    AC_MSG_RESULT([yes, in /usr/local])
    ss_done=yes
fi

if test $ss_done = no && test -f fftw-2_1_3_tar.gz ; then
    FFTW_LDFLAGS=-L/usr/local/spikesorter/lib
    FFTW_CPPFLAGS=-I/usr/local/spikesorter/include
    FFTW_DEP=fftw
    rm -f fftw
    AC_MSG_RESULT([no, will install in /usr/local/spikesorter])
    ss_done=yes
fi

if test $ss_done = no ; then
    AC_MSG_FAILURE([no fftw library found])
fi

LDFLAGS=$ss_old_LDFLAGS
CPPFLAGS=$ss_old_CPPFLAGS

AC_SUBST(FFTW_DEP)
AC_SUBST(FFTW_LDFLAGS)
AC_SUBST(FFTW_CPPFLAGS)

])


AC_DEFUN([SS_MESCHACH_TEST],
[AC_MSG_CHECKING([for meschach library >= 1.2b])
AC_DEFUN([SS_CHECK_MESCHACH],
ss_save_LIBS=$LIBS
LIBS="$LIBS -lmeschach"
[rm -f meschach_test.out
AC_RUN_IFELSE([AC_LANG_SOURCE([
#include <stdio.h>
#include <meschach/matrix.h>
int
main (void)
{
  freopen ("meschach_test.out", "w", stdout) || (exit (1), 1);
  m_version ();
  return 0;
}
])])
LIBS=$ss_save_LIBS
[ss_mvers=`head -1 meschach_test.out 2> /dev/null | sed -n 's/.* \([0-9]\.[0-9][a-zA-Z]*\)$/\1/p'`]
ss_least=`(echo $ss_mvers; echo 1.2b) | sort | head -1`
if test "$ss_least" = 1.2b ; then ss_pass=yes; else ss_pass=no; fi
rm -f meschach_test.out
])

ss_old_LDFLAGS=$LDFLAGS
ss_old_CPPFLAGS=$CPPFLAGS
ss_done=no
ss_pass=no

if test -f /usr/local/spikesorter/lib/libmeschach.a ; then
    LDFLAGS="-L/usr/local/spikesorter/lib $LDFLAGS"
    CPPFLAGS="-I/usr/local/spikesorter/include $CPPFLAGS"
    SS_CHECK_MESCHACH
fi

if test "$ss_pass" = yes ; then
    MESCHACH_LDFLAGS=-L/usr/local/spikesorter/lib
    MESCHACH_CPPFLAGS=-I/usr/local/spikesorter/include
    AC_MSG_RESULT([yes, in /usr/local/spikesorter])
    ss_done=yes
else
    LDFLAGS=$ss_old_LDFLAGS
    CPPFLAGS=$ss_old_CPPFLAGS
    SS_CHECK_MESCHACH
fi

if test $ss_done = no && test $ss_pass = yes ; then
    unset MESCHACH_LDFLAGS
    unset MESCHACH_CPPFLAGS
    AC_MSG_RESULT([yes])
    ss_done=yes
fi

if test $ss_done = no ; then
    LDFLAGS="-L/usr/local/lib $LDFLAGS"
    CPPFLAGS="-I/usr/local/include $CPPFLAGS"
    SS_CHECK_MESCHACH
fi

if test $ss_done = no && test $ss_pass = yes ; then
    MESCHACH_LDFLAGS=-L/usr/local/lib
    MESCHACH_CPPFLAGS=-I/usr/local/include
    AC_MSG_RESULT([yes, in /usr/local])
    ss_done=yes
fi

if test $ss_done = no && test -f mesch12b.tar.gz ; then
    MESCHACH_LDFLAGS=-L/usr/local/spikesorter/lib
    MESCHACH_CPPFLAGS=-I/usr/local/spikesorter/include
    MESCHACH_DEP=meschach_stamp
    rm -f meschach_stamp
    AC_MSG_RESULT([no, will install in /usr/local/spikesorter])
    ss_done=yes
fi

if test $ss_done = no ; then
    AC_MSG_FAILURE([no meschach library found])
fi

LDFLAGS=$ss_old_LDFLAGS
CPPFLAGS=$ss_old_CPPFLAGS

AC_SUBST(MESCHACH_DEP)
AC_SUBST(MESCHACH_LDFLAGS)
AC_SUBST(MESCHACH_CPPFLAGS)

])

# Configure path for the GNU Scientific Library
# Christopher R. Gabriel <cgabriel@linux.it>, April 2000
AC_DEFUN([SS_GSL_TEST],
[
AC_ARG_WITH(gsl-prefix,[  --with-gsl-prefix=PFX   Prefix where GSL is installed (optional)],
            gsl_prefix="$withval", gsl_prefix="")
AC_ARG_WITH(gsl-exec-prefix,[  --with-gsl-exec-prefix=PFX Exec prefix where GSL is installed (optional)],
            gsl_exec_prefix="$withval", gsl_exec_prefix="")
AC_ARG_ENABLE(gsltest, [  --disable-gsltest       Do not try to compile and run a test GSL program],
		    , enable_gsltest=yes)

  if test "x${GSL_CONFIG+set}" != xset ; then
     if test "x$gsl_prefix" != x ; then
         GSL_CONFIG="$gsl_prefix/bin/gsl-config"
     fi
     if test "x$gsl_exec_prefix" != x ; then
        GSL_CONFIG="$gsl_exec_prefix/bin/gsl-config"
     fi
  fi

  AC_PATH_PROG(GSL_CONFIG, gsl-config, no, [/usr/local/spikesorter/bin:$PATH])
  min_gsl_version=1.4
  AC_MSG_CHECKING(for GSL - version >= $min_gsl_version)
  no_gsl=""
  if test "$GSL_CONFIG" = "no" ; then
    no_gsl=yes
  else
    if test ! -x $GSL_CONFIG ; then
      GSL_CONFIG="/bin/sh $GSL_CONFIG"
    fi
    GSL_CFLAGS=`$GSL_CONFIG --cflags`
    GSL_LIBS=`$GSL_CONFIG --libs`

    gsl_major_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\).*/\1/'`
    if test "x${gsl_major_version}" = "x" ; then
       gsl_major_version=0
    fi

    gsl_minor_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\2/'`
    if test "x${gsl_minor_version}" = "x" ; then
       gsl_minor_version=0
    fi

    gsl_micro_version=`$GSL_CONFIG --version | \
           sed 's/^\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\)\.\{0,1\}\([[0-9]]*\).*/\3/'`
    if test "x${gsl_micro_version}" = "x" ; then
       gsl_micro_version=0
    fi

    if test "x$enable_gsltest" = "xyes" ; then
      ac_save_CFLAGS="$CFLAGS"
      ac_save_LIBS="$LIBS"
      CFLAGS="$CFLAGS $GSL_CFLAGS"
      LIBS="$LIBS $GSL_LIBS"

      rm -f conf.gsltest
      AC_TRY_RUN([
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* my_strdup (const char *str);

char*
my_strdup (const char *str)
{
  char *new_str;
  
  if (str)
    {
      new_str = (char *)malloc ((strlen (str) + 1) * sizeof(char));
      strcpy (new_str, str);
    }
  else
    new_str = NULL;
  
  return new_str;
}

int main (void)
{
  int major = 0, minor = 0, micro = 0;
  int n;
  char *tmp_version;

  system ("touch conf.gsltest");

  /* HP/UX 9 (%@#!) writes to sscanf strings */
  tmp_version = my_strdup("$min_gsl_version");

  n = sscanf(tmp_version, "%d.%d.%d", &major, &minor, &micro) ;

  if (n != 2 && n != 3) {
     printf("%s, bad version string\n", "$min_gsl_version");
     exit(1);
   }

   if (($gsl_major_version > major) ||
      (($gsl_major_version == major) && ($gsl_minor_version > minor)) ||
      (($gsl_major_version == major) && ($gsl_minor_version == minor) && ($gsl_micro_version >= micro)))
    {
      exit(0);
    }
  else
    {
      exit(1);
    }
}

],, no_gsl=yes,[echo $ac_n "cross compiling; assumed OK... $ac_c"])
       CFLAGS="$ac_save_CFLAGS"
       LIBS="$ac_save_LIBS"
     fi
  fi
  if test "x$no_gsl" = x ; then
     AC_MSG_RESULT(yes)
  else
      if test -f gsl-1.4.tar.gz ; then
	  GSL_DEP=gsl
	  rm -f gsl
	  GSL_CFLAGS=-I/usr/local/spikesorter/include
	  GSL_LIBS="-L/usr/local/spikesorter/lib -lgsl -lgslcblas"
	  AC_MSG_RESULT([no, will install in /usr/local/spikesorter])
      else
	  AC_MSG_FAILURE([no gsl library found])
      fi
  fi
  AC_SUBST(GSL_CFLAGS)
  AC_SUBST(GSL_LIBS)
  AC_SUBST(GSL_DEP)
  rm -f conf.gsltest
])
