#ifndef verify_h
#define verify_h
#include <fstream>
#include <map>
#include "exafmm.h"

namespace exafmm {
  class Verify {
    typedef std::map<uint64_t,double> Record;
    typedef Record::iterator R_iter;

  private:
    const char * path;

  public:
    bool verbose;
    double average, average2;

    Verify() : path("./"), average(0), average2(0) {}
    Verify(const char * _path) : path(_path), average(0), average2(0) {}

    double getSumScalar(const Bodies & bodies) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += bodies[b].p * bodies[b].q;
      }
      return v;
    }

    double getNrmScalar(const Bodies & bodies) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += std::abs(bodies[b].p * bodies[b].p);
      }
      return v;
    }

    double getDifScalar(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += std::abs((bodies[b].p - bodies2[b].p) *
                      (bodies[b].p - bodies2[b].p));
      }
      return v;
    }

    double getRelScalar(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += std::abs(((bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p))
                      / (bodies2[b].p * bodies2[b].p));
      }
      return v;
    }

    double getNrmVector(const Bodies & bodies) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += std::abs(norm(bodies[b].F));
      }
      return v;
    }

    double getDifVector(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += std::abs((bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
                      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]));
      }
      return v;
    }

    double getRelVector(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        v += std::abs(((bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
                       (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]))
                      / norm(bodies2[b].F));
      }
      return v;
    }

    bool accuracyRegression(const uint64_t key, double value, double value2) {
      bool pass = false;
      Record record, record2;
      std::stringstream name;
      name << path << "accuracy.reg";
      std::ifstream ifile(name.str().c_str());
      if (verbose) std::cout << "opening: " << name.str() << std::endl;
      int numKeys;
      if (ifile.good()) {
        if (verbose) std::cout << "file exists" << std::endl;
        ifile >> numKeys;
        for (int i=0; i<numKeys; i++) {
          uint64_t readKey;
          ifile >> readKey;
          ifile >> record[readKey];
          ifile >> record2[readKey];
        }
      }

      if (record[key] != 0 && verbose) std::cout << "entry exists" << std::endl;
      double threshold = 1 + 1e-8 + 0.01;
      if ((record[key] == 0 || value <= threshold*record[key]) && (average < 5e-4) &&
        (record2[key] == 0 || value2 <= threshold*record2[key]) && (average < 5e-3)) {
        pass = true;
        record[key] = value;
        record2[key] = value2;
      }

      if (!pass) {
        std::cout << "Accuracy regression failed: " <<
              std::scientific << value << " / " << record[key] << std::endl;
        std::cout << "                            " <<
              value2 << " / " << record2[key] << std::endl;
      } else {
        std::ofstream ofile(name.str().c_str());
        ofile << record.size() << std::endl;
        R_iter R2 = record2.begin();
        for (R_iter R=record.begin(); R!=record.end(); R++,R2++) {
          ofile << R->first << " " << R->second;
          ofile << " " << R2->second;
          ofile << std::endl;
        }
      }
      return pass;
    }


    bool timeRegression(const uint64_t key, double value) {
      bool pass = false;
      Record record, record2;
      std::stringstream name;
      name << path << "time.reg";
      std::ifstream ifile(name.str().c_str());
      if (verbose) std::cout << "opening: " << name.str() << std::endl;
      int numKeys;
      if (ifile.good()) {
        if (verbose) std::cout << "file exists" << std::endl;
        ifile >> numKeys;
        for (int i=0; i<numKeys; i++) {
          uint64_t readKey;
          ifile >> readKey;
          ifile >> record[readKey];
          ifile >> record2[readKey];
        }
      }
    }

    if (record[key] != 0 && verbose) std::cout << "entry exists" << std::endl;
    double threshold = 1 + 1e-8 + 0.01;
    if (record[key] == 0) {
      pass = true;
      record[key] = value;
      record2[key] = 1;
    } else if (value <= threshold*record[key]) {
      pass = true;
      record[key] = (record[key] * record2[key] + value) / (record2[key] + 1);
      record2[key] += 1;
    }

    if (!pass) {
      std::cout << "Time regression failed: " <<
        value << " / " << record[key] << std::endl;
    } else {
      std::ofstream ofile(name.str().c_str());
      ofile << record.size() << std::endl;
      R_iter R2 = record2.begin();
      for (R_iter R=record.begin(); R!=record.end(); R++,R2++) {
        ofile << R->first << " " << R->second;
        ofile << " " << (int)(R2->second);
        ofile << std::endl;
      }
    }
    return pass;
  };
}
#endif
