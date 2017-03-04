SUBDIRS = tests

buildbot:
	buildbot try --connect=pb --master=131.112.5.69:443 --username=exafmm --passwd=exafmm --vc=git
cleanbin:
	@find . -name "*.o" -o -name "*.out*" | xargs rm -rf
cleandat:
	@find . -name "*.dat" | xargs rm -f
cleanlib:
	@find . -name "*.a" -o -name "*.so" | xargs rm -f
cleanreg:
	@find . -name "*.reg" | xargs rm -f
cleanall:
	make cleanbin
	make cleandat
	make cleanlib
	make cleanreg
