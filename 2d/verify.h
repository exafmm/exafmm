#ifndef verify_h
#define verify_h
#include <fstream>
#include <map>
#include "exafmm.h"

namespace exafmm {
  //! Verify results
  class Verify {
    typedef std::map<uint64_t,double> Record;                   //!< Map of regression key value pair
    typedef Record::iterator R_iter;                            //!< Iterator of regression map

  private:
    const char * path;                                          //!< Path to save files

  public:
    bool verbose;                                               //!< Print to screen
    double average, average2;                                   //!< Average for regression

    //! Constructor
    Verify() : path("./"), average(0), average2(0) {}

    //! Constructor with argument
    Verify(const char * _path) : path(_path), average(0), average2(0) {}

    //! Get sum of scalar component of a vector of target bodies
    double getSumScalar(const Bodies & bodies) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies
        v += bodies[b].p * bodies[b].q;                         //  Sum of scalar component for Laplace
      }                                                         // End loop over bodies
      return v;                                                 // Return difference
    }

    //! Get norm of scalar component of a vector of target bodies
    double getNrmScalar(const Bodies & bodies) {
      double v = 0;                                             // Initialize norm
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies
        v += std::abs(bodies[b].p * bodies[b].p);               //  Norm of scalar component
      }                                                         // End loop over bodies
      return v;                                                 // Return norm
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getDifScalar(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies & bodies2
        v += std::abs((bodies[b].p - bodies2[b].p) *            //  Difference of scalar component
                      (bodies[b].p - bodies2[b].p));
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get relative difference between scalar component of two vectors of target bodies
    double getRelScalar(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies & bodies2
        v += std::abs(((bodies[b].p - bodies2[b].p) * (bodies[b].p - bodies2[b].p))
		                / (bodies2[b].p * bodies2[b].p));           //  Difference of scalar component
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get norm of scalar component of a vector of target bodies
    double getNrmVector(const Bodies & bodies) {
      double v = 0;                                             // Initialize norm
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies
        v += std::abs(norm(bodies[b].F));                       //  Norm of vector component
      }                                                         // End loop over bodies
      return v;                                                 // Return norm
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getDifVector(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies & bodies2
        v += std::abs((bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
                      (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]));
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getRelVector(const Bodies & bodies, const Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      for (size_t b=0; b<bodies.size(); b++) {                  // Loop over bodies & bodies2
        v += std::abs(((bodies[b].F[0] - bodies2[b].F[0]) * (bodies[b].F[0] - bodies2[b].F[0]) +
		                   (bodies[b].F[1] - bodies2[b].F[1]) * (bodies[b].F[1] - bodies2[b].F[1]))
		                  / norm(bodies2[b].F));
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Perform regression test on accuracy
    bool accuracyRegression(const uint64_t key, double value, double value2) {
      bool pass = false;                                        // Flag for regression test
      Record record, record2;                                   // Map for regression value
      std::stringstream name;                                   // File name for regression
      name << path << "accuracy.reg";                           // Append path and file name
      
      {                                                         // Read record into regression maps
        std::ifstream ifile(name.str().c_str());                //  Input file stream
        if (verbose) std::cout << "opening: " << name.str() << std::endl;  // Print file name
        int numKeys;                                            //  Number of keys stored in file
        if (ifile.good()) {                                     //  If file exists
          if (verbose) std::cout << "file exists" << std::endl; //   Print if file exists
          ifile >> numKeys;                                     //   Read number of keys
          for (int i=0; i<numKeys; i++) {                       //   Loop over regression values
            uint64_t readKey;                                   //    Read key buffer
            ifile >> readKey;                                   //    Read key
            ifile >> record[readKey];                           //    Read value
            ifile >> record2[readKey];                          //    Read value2
          }                                                     //   End loop over regression values
        }                                                       //  End if file exists
      }                                                         // End read record

      if (record[key] != 0 && verbose) std::cout << "entry exists" << std::endl;       // Print if entry exits
      double threshold = 1 + 1e-8 + 0.01;                       // Accuracy threshold
      if ((record[key] == 0 || value <= threshold*record[key]) && (average < 5e-4) &&  // If new record or pass
        (record2[key] == 0 || value2 <= threshold*record2[key]) && (average < 5e-3)) {
        pass = true;                                            //  Change flag to pass
        record[key] = value;                                    //  Add key value pair
        record2[key] = value2;                                  //  Add key value2 pair
      }                                                         // End if for better value

      if (!pass) {                                              // If test does not pass
        std::cout << "Accuracy regression failed: " <<          //  Print value vs. record
              std::scientific << value << " / " << record[key] << std::endl;
        std::cout << "                            " <<          //  Print value2 vs. record2
              value2 << " / " << record2[key] << std::endl;
      } else {                                                  // else test passes
        std::ofstream ofile(name.str().c_str());                // Output file stream
        ofile << record.size() << std::endl;                    // Write number of keys
        R_iter R2 = record2.begin();                            // Iterator for record2
        for (R_iter R=record.begin(); R!=record.end(); R++,R2++) { // Loop over regression values
          ofile << R->first << " " << R->second;                //  Write key value pair
          ofile << " " << R2->second;                           //  Write second value
          ofile << std::endl;                                   //  End line
        }                                                       // End loop over regression values
      }
      return pass;                                              // return regression result
    }

    //! Perform regression test on performance
    bool timeRegression(const uint64_t key, double value) {
      bool pass = false;                                        // Flag for regression test
      Record record, record2;                                   // Map for regression value
      std::stringstream name;                                   // File name for regression
      name << path << "time.reg";                               // Append path and file name

      {                                                         // Read record into regression maps
        std::ifstream ifile(name.str().c_str());                //  Input file stream
        if (verbose) std::cout << "opening: " << name.str() << std::endl;//  Print file name
        int numKeys;                                            //  Number of keys stored in file
        if (ifile.good()) {                                     //  If file exists
          if (verbose) std::cout << "file exists" << std::endl; //   Print if file exists
          ifile >> numKeys;                                     //   Read number of keys
          for (int i=0; i<numKeys; i++) {                       //   Loop over regression values
            uint64_t readKey;                                   //    Read key buffer
            ifile >> readKey;                                   //    Read key
            ifile >> record[readKey];                           //    Read value
            ifile >> record2[readKey];                          //    Read value2
          }                                                     //  End loop over regression values
        } 
      }

      if (record[key] != 0 && verbose) std::cout << "entry exists" << std::endl;// Print if entry exits
      double threshold = 1 + 1e-8 + 0.01;                       // Accuracy threshold
      if (record[key] == 0) {                                   // If new record
        pass = true;                                            //  Change flag to pass
        record[key] = value;                                    //  Add key value pair
        record2[key] = 1;                                       //  Add key count pair
      } else if (value <= threshold*record[key]) {              // Else if passes
        pass = true;                                            //  Change flag to pass
        record[key] = (record[key] * record2[key] + value) / (record2[key] + 1);  // Update timing
        record2[key] += 1;                                      //  increment count
      }

      if (!pass) {                                              // If regression failed
        std::cout << "Time regression failed: " <<              //  Print message for time regression
                    value << " / " << record[key] << std::endl; //  Print value and record 
      } else {                                                  // If regression passed
        std::ofstream ofile(name.str().c_str());                //  Output file stream
        ofile << record.size() << std::endl;                    //   Write number of keys
        R_iter R2 = record2.begin();                            //   Iterator for record2
        for (R_iter R=record.begin(); R!=record.end(); R++,R2++) { //   Loop over regression values
          ofile << R->first << " " << R->second;                   //    Write key timing pair
          ofile << " " << (int)(R2->second);                       //    Write key count pair
          ofile << std::endl;                                      //    End line
        }                                                         //    End loop over regression values
      }
      return pass;
    }
  };
}
#endif
