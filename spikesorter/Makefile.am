AM_CXXFLAGS = -Wall -std=gnu++11

bin_PROGRAMS = spikesort dmx isi_pos split_merge_abf			\
	       spiketime_to_region tkss dmx_split integrate rpls cpls	\
	       trachpls waveform muph edt_merge ctr_to_chan chan2hdt	\
	       nfs_lock edt2pos dup_ids filter60 crossings2pos		\
	       crossings2edt align

align_SOURCES = align.cc
align_CPPFLAGS = $(GTK_CFLAGS) $(GSL_CFLAGS)
align_LDADD = $(GTK_LIBS) $(GSL_LIBS)

crossings2pos_SOURCES = crossings2pos.c util.c util.h

crossings2edt_SOURCES = crossings2edt.c util.c util.h

filter60_SOURCES = filter60.c
filter60_CPPFLAGS = $(FFTW3_CFLAGS)
filter60_LDADD = $(FFTW3_LIBS)

edt2pos_SOURCES = edt2pos.c die.c util.h

dup_ids_SOURCES = dup_ids.c

nfs_lock_SOURCES = nfs_lock.c

bin_PROGRAMS += requester

requester_SOURCES = requester.c old_requester.c

man_MANS = spikesort_control_panel.1

ctr_to_chan_SOURCES = ctr_to_chan.c nde.h region.h die.c

chan2hdt_SOURCES = chan2hdt.c nde.h die.c

bin_SCRIPTS = tkss.tcl spikesort_control_panel waveform.tcl	\
info_file.pl alert
EXTRA_DIST = waveform.tcl info_file.pl $(man_MANS) debian	\
user_guide.tex utilities.tex userguide-img001.png		\
userguide-img002.png raw.png diagram.png results.png		\
userguide-img006.png userguide-img007.png userguide-img008.png	\
tools.jpg user_guide.bbl


spikesort_SOURCES = cdf.c cspline.c dmx.h ksone.c min.c nde.h	\
qselect.c region.h window.c die.c dmx_th.c lfs.c mpk_regions.c	\
pdfrender.c region.c residue.c find_detection_threshold.c

FFTW_CFLAGS = 
FFTW_LIBS = -lrfftw -lfftw

MESCHACH_CFLAGS = -I/usr/include/meschach
MESCHACH_LIBS = -lmeschach

TCLTK_CFLAGS = -I/usr/include/tcl -I/usr/include/tk
TCLTK_LIBS = -ltcl -ltk

READLINE_CFLAGS = -I/usr/include/readline
READLINE_LIBS =

AM_CPPFLAGS = $(FFTW_CFLAGS) $(GSL_CFLAGS)
AM_CFLAGS = -Wall --std=c99
LDADD = $(LIBOBJS)

spikesort_LDADD = $(MESCHACH_LIBS) $(GSL_LIBS) $(X_LIBS) $(FFTW_LIBS)	\
	 $(LIBOBJS)

spikesort_CPPFLAGS = $(AM_CPPFLAGS) $(GSL_CFLAGS) $(MESCHACH_CFLAGS)	\
$(X_CFLAGS)

tkss_SOURCES = tkss.c nde.h lfs.c dmx.h die.c

tkss_LDADD = $(TCLTK_LIBS) $(X_LIBS) $(LIBOBJS) $(GSL_LIBS)
tkss_CPPFLAGS = $(AM_CPPFLAGS) $(TCLTK_CFLAGS) $(X_CFLAGS) $(GSL_CFLAGS)

waveform_SOURCES = waveform.c tmalloc.h lfs.c dmx.h die.c

waveform_LDADD = $(TCLTK_LIBS) $(X_LIBS) $(LIBOBJS)
waveform_CPPFLAGS = $(AM_CPPFLAGS) $(TCLTK_CFLAGS) $(X_CFLAGS)

dmx_SOURCES = dmx.c lfs.c nde.h dmx.h die.c
dmx_CFLAGS = -fno-strict-aliasing

dmx_split_SOURCES = dmx_split.c lfs.c tmalloc.h dmx.h lfs.h die.c

isi_pos_SOURCES = isi_pos.c die.c nde.h

split_merge_abf_SOURCES = split_merge_abf.c

spiketime_to_region_SOURCES = die.c nde.h region.h spiketime_to_region.c

integrate_SOURCES = integrate.cc

muph_SOURCES = muph.c die.c
muph_CFLAGS = -fno-strict-aliasing --std=c99
muph_LDADD = $(GSL_LIBS)

rpls_SOURCES = rpls.c die.c
rpls_LDADD = $(GSL_LIBS)

cpls_SOURCES = cpls.c die.c

trachpls_SOURCES = trachpls.c die.c

edt_merge_SOURCES = edt_merge.c getline.c getline.h
edt_merge_CPPFLAGS = $(READLINE_CFLAGS)

dist_doc_DATA = README utilities.pdf user_guide.pdf

install-exec-hook:
	cd $(DESTDIR)/usr/local/bin; rm -f digitize; $(LN_S) integrate digitize

SRC_FILES = $(ctr_to_chan_SOURCES) $(chan2hdt_SOURCES)					\
            $(spikesort_SOURCES) $(tkss_SOURCES) $(waveform_SOURCES)			\
            $(dmx_SOURCES) $(dmx_split_SOURCES) $(isi_pos_SOURCES)			\
            $(split_merge_abf_SOURCES) $(spiketime_to_region_SOURCES)			\
            $(integrate_SOURCES) $(muph_SOURCES) $(rpls_SOURCES) $(cpls_SOURCES)	\
            $(trachpls_SOURCES) $(edt_merge_SOURCES) $(bin_SCRIPTS)

srcfiles: $(SRC_FILES) Makefile.am
	ls $(SRC_FILES) | sort -u > srcfiles

user_guide.pdf: user_guide.tex
	latexmk -f -gg -pdf user_guide.tex

utilities.pdf: utilities.tex
	latexmk -f -gg -pdf utilities.tex



checkin_files = $(SRC_FILES) Makefile.am configure.ac ChangeLog


add_head:
	for n in $(checkin_files) ; do \
		touch -r $$n oldtime ; \
		vi $$n ; \
		touch -r oldtime $$n ;\
	done

deb:
	@echo 'Making debian packages'
	make distdir &&\
	cd $(distdir) &&\
	dpkg-buildpackage -Zgzip -uc -us && \
	cd .. && \
	rm -rf $(distdir)



