SUBDIRS = 2d 3d 3dp 3dp_gpu 3dp_gpu_mpi 3dp_simd 3dp_simd_mpi

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
