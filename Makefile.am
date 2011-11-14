LIBDIR=@LIBDIR@

ACLOCAL_AMFLAGS = -I m4

SUBDIRS=doc src
RJVERABBR=$(PACKAGE_VERSION)-$(build)$(MPISEP)$(MPIINFIX)
RJVERABBREXE=$(RJVERABBR)$(EXEEXT)
RJBPVA="rjungle-bin-pkg-$(RJVERABBR)"

eintest:
	echo $(MAKE_FOR_WINDOWS)
	echo $(RJVERABBREXE)

dist-hook:
	cp $(top_srcdir)/VERSION $(distdir)/
	cp $(top_srcdir)/THANKS $(distdir)/
	mkdir $(distdir)/tests
	cp $(top_srcdir)/tests/data_from_r.dat $(distdir)/tests/
	cp $(top_srcdir)/tests/test.ped $(distdir)/tests/
	cp $(top_srcdir)/tests/test.dat $(distdir)/tests/
	cp $(top_srcdir)/tests/check.* $(distdir)/tests/
	cp $(top_srcdir)/src/$(LIBDIR)/error_codes.txt $(distdir)/src/$(LIBDIR)/

dist-internet: dist dist-internet-doc dist-internet-static
	$(top_srcdir)/doc/putfilestoftp.sh randomjungle-$(PACKAGE_VERSION).tar.gz

dist-internet-doc: dist
	cd $(top_srcdir)/doc && ./putfilestoftp.sh randomjungle.pdf
	cd $(top_srcdir)/doc && ./putfilestoftp.sh randomjungle.html
	cd $(top_srcdir)/doc && cp ../src/$(LIBDIR)/error_codes.txt . && ./putfilestoftp.sh error_codes.txt && rm error_codes.txt
	cd $(top_srcdir)/doc && cp ../NEWS . && ./putfilestoftp.sh NEWS && rm NEWS

dist-internet-static:
	mkdir -p $(top_srcdir)/$(RJBPVA)
	mkdir -p $(distdir)/bin-pkg
	cp -Rf $(top_srcdir)/pkgskeleton$(MPIINFIX)/* $(top_srcdir)/$(RJBPVA)/
	find $(top_srcdir)/$(RJBPVA)/ -name .svn | xargs -IFILES rm -Rf FILES
	if test "$(MAKE_FOR_WINDOWS)" == 'yes' ; then rm -f $(top_srcdir)/$(RJBPVA)/*.sh; else rm -f $(top_srcdir)/$(RJBPVA)/*.bat; rm -f $(top_srcdir)/$(RJBPVA)/bin/*.dll; fi
	cp -Rf $(top_srcdir)/$(RJBPVA)/* $(distdir)/bin-pkg/
	if test "$(MAKE_FOR_WINDOWS)" == 'yes' ; then cp $(top_srcdir)/src/rjunglesparse_static$(EXEEXT) $(top_srcdir)/$(RJBPVA)/bin/rjunglesparse$(EXEEXT) && cp $(top_srcdir)/src/rjungle_static$(EXEEXT) $(top_srcdir)/$(RJBPVA)/bin/rjungle$(EXEEXT); zip -r $(RJBPVA).zip $(top_srcdir)/$(RJBPVA); $(top_srcdir)/doc/putfilestoftp.sh $(RJBPVA).zip; else cp $(top_srcdir)/src/rjunglesparse_static$(EXEEXT) $(top_srcdir)/$(RJBPVA)/bin/rjunglesparse-$(RJVERABBREXE) && cp $(top_srcdir)/src/rjungle_static$(EXEEXT) $(top_srcdir)/$(RJBPVA)/bin/rjungle-$(RJVERABBREXE); ln -s rjunglesparse-$(RJVERABBREXE) $(top_srcdir)/$(RJBPVA)/bin/rjunglesparse && ln -s rjungle-$(RJVERABBREXE) $(top_srcdir)/$(RJBPVA)/bin/rjungle; tar -zcf $(RJBPVA).tar.gz $(top_srcdir)/$(RJBPVA); $(top_srcdir)/doc/putfilestoftp.sh $(RJBPVA).tar.gz; fi

check-am:
	$(top_srcdir)/src/rjungle$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.ped -p -o $(top_srcdir)/tests/testa1
	$(top_srcdir)/src/rjunglesparse$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.ped -p -o $(top_srcdir)/tests/testsparsea1
	test "`diff $(top_srcdir)/tests/testa1.importance $(top_srcdir)/tests/testsparsea1.importance | wc | tr -s ' '`" = " 0 0 0"
	$(top_srcdir)/src/rjungle$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.dat -D y -o $(top_srcdir)/tests/testa2
	$(top_srcdir)/src/rjungle$(EXEEXT) --testlib
	$(top_srcdir)/src/rjunglesparse$(EXEEXT) --testlib
	$(top_srcdir)/src/rjungle$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.dat -D y -t1 -w2 -U1 -o $(top_srcdir)/tests/testjungle
	test "`diff $(top_srcdir)/tests/check.importance $(top_srcdir)/tests/testjungle.importance | wc | tr -s ' '`" = " 0 0 0"
	test "`diff $(top_srcdir)/tests/check.jungle.xml $(top_srcdir)/tests/testjungle.jungle.xml | wc | tr -s ' '`" = " 0 0 0"
	$(top_srcdir)/src/rjungle$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.dat -D y -t1 -w2 -U1 -o $(top_srcdir)/tests/testpred -P $(top_srcdir)/tests/testjungle.jungle.xml
	test "`diff $(top_srcdir)/tests/check.prediction $(top_srcdir)/tests/testpred.prediction | wc | tr -s ' '`" = " 0 0 0"
	$(top_srcdir)/src/rjungle_static$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.ped -p -o $(top_srcdir)/tests/testb1
	$(top_srcdir)/src/rjunglesparse_static$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.ped -p -o $(top_srcdir)/tests/testsparseb1
	test "`diff $(top_srcdir)/tests/testb1.importance $(top_srcdir)/tests/testsparseb1.importance | wc | tr -s ' '`" = " 0 0 0"
	$(top_srcdir)/src/rjungle_static$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.dat -D y -o $(top_srcdir)/tests/testb2
	$(top_srcdir)/src/rjungle_static$(EXEEXT) --testlib
	$(top_srcdir)/src/rjunglesparse_static$(EXEEXT) --testlib
	$(top_srcdir)/src/rjungle_static$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.dat -D y -t1 -w2 -U1 -o $(top_srcdir)/tests/testjungle
	test "`diff $(top_srcdir)/tests/check.importance $(top_srcdir)/tests/testjungle.importance | wc | tr -s ' '`" = " 0 0 0"
	test "`diff $(top_srcdir)/tests/check.jungle.xml $(top_srcdir)/tests/testjungle.jungle.xml | wc | tr -s ' '`" = " 0 0 0"
	$(top_srcdir)/src/rjungle_static$(EXEEXT) -z 123 -f $(top_srcdir)/tests/test.dat -D y -t1 -w2 -U1 -o $(top_srcdir)/tests/testpred -P $(top_srcdir)/tests/testjungle.jungle.xml
	test "`diff $(top_srcdir)/tests/check.prediction $(top_srcdir)/tests/testpred.prediction | wc | tr -s ' '`" = " 0 0 0"


clean-am:
	rm -f $(top_srcdir)/tests/test1.*
	rm -f $(top_srcdir)/tests/test2.*
	rm -Rf $(top_srcdir)/$(RJBPVA)
