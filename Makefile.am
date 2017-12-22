SUBDIRS = 3dp

buildbot:
	buildbot try --connect=pb --master=131.112.5.69:443 --username=exafmm --passwd=exafmm --vc=git
cleanbin:
	@find . -name "*.o" -o -name "*.out*" | xargs rm -rf
cleandat:
	@find . -name "*.dat" | xargs rm -rf
cleanlib:
	@find . -name "*.a" -o -name "*.so" | xargs rm -f
cleanlog:
	@find . -name "*.log" -o -name "*.p" | xargs rm -f
cleanall:
	make cleanbin
	make cleandat
	make cleanlib
	make cleanlog
