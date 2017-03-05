#ifndef args_h
#define args_h
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include "namespace.h"
#include <stdint.h>
#include "types.h"

namespace EXAFMM_NAMESPACE {
  static struct option long_options[] = {
    {"accuracy",     no_argument,       0, 'a'},
    {"ncrit",        required_argument, 0, 'c'},
    {"cutoff",       required_argument, 0, 'C'},
    {"distribution", required_argument, 0, 'd'},
    {"dual",         no_argument,       0, 'D'},
    {"graft",        no_argument,       0, 'g'},
    {"getMatrix",    no_argument,       0, 'G'},
    {"help",         no_argument,       0, 'h'},
    {"images",       required_argument, 0, 'i'},
    {"IneJ",         no_argument,       0, 'j'},
    {"numBodies",    required_argument, 0, 'n'},
    {"path",         required_argument, 0, 'p'},
    {"P",            required_argument, 0, 'P'},
    {"repeat",       required_argument, 0, 'r'},
    {"nspawn",       required_argument, 0, 's'},
    {"theta",        required_argument, 0, 't'},
    {"threads",      required_argument, 0, 'T'},
    {"verbose",      no_argument,       0, 'v'},
    {"write",        no_argument,       0, 'w'},
    {0, 0, 0, 0}
  };

  class Args {
  public:
    int accuracy;
    int ncrit;
    double cutoff;
    const char * distribution;
    int dual;
    int graft;
    int getMatrix;
    int images;
    int IneJ;
    int numBodies;
    const char * path;
    int P;
    int repeat;
    int nspawn;
    double theta;
    int threads;
    int verbose;
    int write;

  private:
    void usage(char * name) {
      fprintf(stderr,
	      "Usage: %s [options]\n"
	      "Long option (short option)       : Description (Default value)\n"
	      " --accuracy (-a)                 : Regression for accuracy only (%d)\n"
	      " --ncrit (-c)                    : Number of bodies per leaf cell (%d)\n"
	      " --cutoff (-C)                   : Cutoff distance of interaction (%f)\n"
	      " --distribution (-d) [c/l/o/p/s] : lattice, cube, sphere, octant, plummer (%s)\n"
	      " --dual (-D)                     : Use dual tree traversal (%d)\n"
	      " --graft (-g)                    : Graft remote trees to global tree (%d)\n"
	      " --getMatrix (-G)                : Write G matrix to file (%d)\n"
	      " --help (-h)                     : Show this help document\n"
	      " --images (-i)                   : Number of periodic image levels (%d)\n"
	      " --IneJ (-j)                     : Use different sources & targets (%d)\n"
	      " --numBodies (-n)                : Number of bodies (%d)\n"
	      " --path (-p)                     : Path to save files (%s)\n"
	      " --P (-P) not working            : Order of expansion (%d)\n"
	      " --repeat (-r)                   : Number of iteration loops (%d)\n"
	      " --nspawn (-s)                   : Threshold for stopping task creation during recursion (%d)\n"
	      " --theta (-t)                    : Multipole acceptance criterion (%f)\n"
	      " --threads (-T)                  : Number of threads (%d)\n"
	      " --verbose (-v)                  : Print information to screen (%d)\n"
	      " --write (-w)                    : Write timings to file (%d)\n",
	      name,
              accuracy,
	      ncrit,
	      cutoff,
	      distribution,
	      dual,
	      graft,
	      getMatrix,
	      images,
	      IneJ,
	      numBodies,
              path,
	      P,
	      repeat,
	      nspawn,
	      theta,
	      threads,
	      verbose,
	      write);
    }

    const char * parseDistribution(const char * arg) {
      switch (arg[0]) {
      case 'c':
	return "cube";
      case 'l':
	return "lattice";
      case 'o':
	return "octant";
      case 'p':
	return "plummer";
      case 's':
	return "sphere";
      default:
	fprintf(stderr, "invalid distribution %s\n", arg);
	abort();
      }
      return "";
    }

    uint64_t getBasisNum(const char * _basis) {
      switch (_basis[0]) {
      case 'C':
        return 0;
      case 'S':
        return 1;
      case 'P':
        return 2;
      case 'E':
        return 3;
      default:
        fprintf(stderr, "invalid basis %s\n", _basis);
        abort();
      }
      return 0;
    }

    uint64_t getDistNum(const char * _distribution) {
      switch (_distribution[0]) {
      case 'c':
        return 0;
      case 'l':
        return 1;
      case 'o':
        return 2;
      case 'p':
        return 3;
      case 's':
        return 4;
      case 'e':
        return 7;
      default:
        fprintf(stderr, "invalid distribution %s\n", _distribution);
        abort();
      }
      return 0;
    }

    uint64_t getEqNum() {
#if EXAFMM_LAPLACE
      return 0;
#elif EXAFMM_HELMHOLTZ
      return 1;
#elif EXAFMM_STOKES
      return 2;
#elif EXAFMM_BIOTSAVART
      return 3;
#endif
    }

    uint64_t getConfigNum() {
      uint64_t key = 0;
#if EXAFMM_SINGLE
      key |= 1;
#endif
#if EXAFMM_USE_SIMD
      key |= 2;
#endif
#if EXAFMM_USE_KAHAN
      key |= 4;
#endif
      return key;
    }

  public:
    Args(int argc=0, char ** argv=NULL) :
      accuracy(0),
      ncrit(64),
      cutoff(.0),
      distribution("cube"),
      dual(0),
      graft(0),
      getMatrix(0),
      images(0),
      IneJ(0),
      numBodies(1000000),
      path("./"),
      P(10),
      repeat(1),
      nspawn(5000),
      theta(.4),
      threads(16),
      verbose(0),
      write(0) {
      while (1) {
	int option_index;
	int c = getopt_long(argc, argv, "ab:c:d:De:gGhi:jmMn:op:P:r:s:t:T:vwx", long_options, &option_index);
	if (c == -1) break;
	switch (c) {
	case 'a':
	  accuracy = 1;
	  break;
	case 'c':
	  ncrit = atoi(optarg);
	  break;
	case 'C':
	  cutoff = atof(optarg);
	  break;
	case 'd':
	  distribution = parseDistribution(optarg);
	  break;
	case 'D':
	  dual = 1;
	  break;
	case 'g':
	  graft = 1;
	  break;
	case 'G':
	  getMatrix = 1;
	  break;
	case 'h':
	  usage(argv[0]);
	  abort();
	case 'i':
	  images = atoi(optarg);
	  break;
	case 'j':
	  IneJ = 1;
	  break;
	case 'n':
	  numBodies = atoi(optarg);
	  break;
        case 'p':
          path = optarg;
          break;
	case 'P':
	  P = atoi(optarg);
	  break;
	case 'r':
	  repeat = atoi(optarg);
	  break;
	case 's':
	  nspawn = atoi(optarg);
	  break;
	case 't':
	  theta = atof(optarg);
	  break;
	case 'T':
	  threads = atoi(optarg);
	  break;
	case 'v':
	  verbose = 1;
	  break;
	case 'w':
	  write = 1;
	  break;
	default:
	  usage(argv[0]);
	  abort();
	}
      }
    }

    uint64_t getKey(int mpisize=1) {
      uint64_t key = 0;
      key |= uint64_t(round(log(ncrit)/log(2)));
      key |= getDistNum(distribution) << 4;
      key |= dual << 7;
      key |= getEqNum() << 8;
      key |= graft << 11;
      key |= images << 12;
      key |= IneJ << 14;
      key |= uint64_t(round(log(numBodies)/log(10))) << 15;
      key |= P << 19;
      key |= uint64_t(round(log(nspawn)/log(10))) << 25;
      key |= uint64_t(theta*14) << 28;
      key |= uint64_t(round(log(threads)/log(2))) << 31;
      key |= getConfigNum() << 34;
      key |= uint64_t(round(log(mpisize)/log(2))) << 37;
      assert( uint64_t(round(log(mpisize)/log(2))) < 26 );
      return key;
    }

    void print(int stringLength) {
      if (verbose) {
	std::cout << std::setw(stringLength) << std::fixed << std::left
		  << "accuracy" << " : " << accuracy << std::endl
		  << std::setw(stringLength)
		  << "ncrit" << " : " << ncrit << std::endl
		  << std::setw(stringLength)
		  << "cutoff" << " : " << cutoff << std::endl
		  << std::setw(stringLength)
		  << "distribution" << " : " << distribution << std::endl
		  << std::setw(stringLength)
		  << "dual" << " : " << dual << std::endl
		  << std::setw(stringLength)
		  << "graft" << " : " << graft << std::endl
		  << std::setw(stringLength)
		  << "getMatrix" << " : " << getMatrix << std::endl
		  << std::setw(stringLength)
		  << "images" << " : " << images << std::endl
		  << std::setw(stringLength)
		  << "IneJ" << " : " << IneJ << std::endl
		  << std::setw(stringLength)
		  << "numBodies" << " : " << numBodies << std::endl
		  << std::setw(stringLength)
		  << "path" << " : " << path << std::endl
		  << std::setw(stringLength)
		  << "P" << " : " << P << std::endl
		  << std::setw(stringLength)
		  << "repeat" << " : " << repeat << std::endl
		  << std::setw(stringLength)
		  << "nspawn" << " : " << nspawn << std::endl
		  << std::setw(stringLength)
		  << "theta" << " : " << theta << std::endl
		  << std::setw(stringLength)
		  << "threads" << " : " << threads << std::endl
		  << std::setw(stringLength)
		  << "verbose" << " : " << verbose << std::endl
		  << std::setw(stringLength)
		  << "write" << " : " << write << std::endl;
      }
    }
  };
}
#endif
