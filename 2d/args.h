#ifndef args_h
#define args_h
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <iomanip>

namespace exafmm {
  static struct option long_options[] = {
    {"ncrit",        required_argument, 0, 'c'},
    {"distribution", required_argument, 0, 'd'},
    {"help",         no_argument,       0, 'h'},
    {"numBodies",    required_argument, 0, 'n'},
    {"P",            required_argument, 0, 'P'},
    {"theta",        required_argument, 0, 't'},
    {"verbose",      no_argument,       0, 'v'},
    {0, 0, 0, 0}
  };

  class Args {
  public:
    int ncrit;
    const char * distribution;
    int numBodies;
    int P;
    double theta;
    int verbose;

  private:
    void usage(char * name) {
      fprintf(stderr,
              "Usage: %s [options]\n"
              "Long option (short option)       : Description (Default value)\n"
              " --ncrit (-c)                    : Number of bodies per leaf cell (%d)\n"
              " --distribution (-d) [c/l/o/p/s] : lattice, cube, sphere, octant, plummer (%s)\n"
              " --help (-h)                     : Show this help document\n"
              " --numBodies (-n)                : Number of bodies (%d)\n"
              " --P (-P) not working            : Order of expansion (%d)\n"
              " --theta (-t)                    : Multipole acceptance criterion (%f)\n"
              " --verbose (-v)                  : Print information to screen (%d)\n",
              name,
              ncrit,
              distribution,
              numBodies,
              P,
              theta,
              verbose);
    }

    const char * parseDistribution(const char * arg) {
      switch (arg[0]) {
        case 'c': return "cube";
        case 'l': return "lattice";
        case 'o': return "octant";
        case 'p': return "plummer";
        case 's': return "sphere";
        default:
          fprintf(stderr, "invalid distribution %s\n", arg);
          abort();
      }
      return "";
    }

  public:
    Args(int argc=0, char ** argv=NULL)
      : ncrit(64), 
        distribution("cube"),
        numBodies(1000000),
        P(10),
        theta(.4),
        verbose(0) 
    {
      while (1) {
        int option_index;
        int c = getopt_long(argc, argv, "c:d:hn:P:t:v",
                            long_options, &option_index);
        if (c == -1) break;
        switch (c) {
          case 'c':
            ncrit = atoi(optarg);
            break;
          case 'd':
            distribution = parseDistribution(optarg);
            break;
          case 'h':
            usage(argv[0]);
            abort();
          case 'n':
            numBodies = atoi(optarg);
            break;
          case 'P':
            P = atoi(optarg);
            break;
          case 't':
            theta = atof(optarg);
            break;
          case 'v':
            verbose = 1;
            break;
          default:
            usage(argv[0]);
            abort();
        }
      }
    }

    void print(int stringLength) {
      if (verbose) {
        std::cout << std::setw(stringLength) << std::fixed << std::left
                  << "ncrit" << " : " << ncrit << std::endl
                  << std::setw(stringLength)
                  << "distribution" << " : " << distribution << std::endl
                  << std::setw(stringLength)
                  << "numBodies" << " : " << numBodies << std::endl
                  << std::setw(stringLength)
                  << "P" << " : " << P << std::endl
                  << std::setw(stringLength)
                  << "theta" << " : " << theta << std::endl
                  << std::setw(stringLength)
                  << "verbose" << " : " << verbose << std::endl
                  << std::setw(stringLength);
      }
    }
  };
}
#endif
