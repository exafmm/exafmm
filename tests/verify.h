#ifndef verify_h
#define verify_h
#include "logger.h"
#include "namespace.h"
#include "types.h"

namespace EXAFMM_NAMESPACE {
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
    double getSumScalar(Bodies & bodies) {
      double v = 0;                                             // Initialize difference
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
#if EXAFMM_LAPLACE
	v += B->TRG[0] * B->SRC;                                //  Sum of scalar component for Laplace
#elif EXAFMM_HELMHOLTZ
	v += std::abs(B->TRG[0] * B->SRC);                      //  Sum of scalar component for Helmholtz
#elif EXAFMM_BIOTSAVART
	v += B->TRG[0];                                         //  Sum of x component for Biot-Savart
#endif
      }                                                         // End loop over bodies
      return v;                                                 // Return difference
    }

    //! Get norm of scalar component of a vector of target bodies
    double getNrmScalar(Bodies & bodies) {
      double v = 0;                                             // Initialize norm
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	v += std::abs(B->TRG[0] * B->TRG[0]);                   //  Norm of scalar component
      }                                                         // End loop over bodies
      return v;                                                 // Return norm
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getDifScalar(Bodies & bodies, Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      B_iter B2 = bodies2.begin();                              // Set iterator of bodies2
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
	v += std::abs((B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0])); //  Difference of scalar component
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getRelScalar(Bodies & bodies, Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      B_iter B2 = bodies2.begin();                              // Set iterator of bodies2
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
	v += std::abs(((B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]))
		      / (B2->TRG[0] * B2->TRG[0]));             //  Difference of scalar component
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get norm of scalar component of a vector of target bodies
    double getNrmVector(Bodies & bodies) {
      double v = 0;                                             // Initialize norm
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over bodies
	v += std::abs(B->TRG[1] * B->TRG[1] +                   //  Norm of vector x component
		      B->TRG[2] * B->TRG[2] +                   //  Norm of vector y component
		      B->TRG[3] * B->TRG[3]);                   //  Norm of vector z component
      }                                                         // End loop over bodies
      return v;                                                 // Return norm
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getDifVector(Bodies & bodies, Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      B_iter B2 = bodies2.begin();                              // Set iterator of bodies2
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
	v += std::abs((B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]) + //  Difference of vector x component
		      (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]) + //  Difference of vector y component
		      (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3])); //  Difference of vector z component
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Get difference between scalar component of two vectors of target bodies
    double getRelVector(Bodies & bodies, Bodies & bodies2) {
      double v = 0;                                             // Initialize difference
      B_iter B2 = bodies2.begin();                              // Set iterator of bodies2
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
	v += std::abs(((B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]) +//  Difference of vector x component
		       (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]) +//  Difference of vector y component
		       (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]))//   Difference of vector z component
		      / (B2->TRG[1] * B2->TRG[1] +              //  Norm of vector x component
			 B2->TRG[2] * B2->TRG[2] +              //  Norm of vector y component
			 B2->TRG[3] * B2->TRG[3]));             //  Norm of vector z component
      }                                                         // End loop over bodies & bodies2
      return v;                                                 // Return difference
    }

    //! Print relative L2 norm scalar error
    void print(std::string title, double v) {
      if (logger::verbose) {                                    // If verbose flag is true
	std::cout << std::setw(logger::stringLength) << std::left //  Set format
		  << title << " : " << std::setprecision(logger::decimal) << std::scientific // Set title
		  << v << std::endl;                            //  Print potential error
      }                                                         // End if for verbose flag
    }

    //! Compare data for regression
    bool regression(uint64_t key, bool time, int iteration, double value, double value2=0) {
      bool pass = false;                                        // Flag for regression test
      bool secondValue = false;                                 // Whether there is a second value
      if (value2 != 0) secondValue = true;                      // If there is a second value
      Record record, record2;                                   // Map for regression value
      const char * host = getenv("WORKERNAME");                 // Get workername
      std::stringstream name;                                   // File name for regression
      name << path;                                             // Append path to file name
      if (time) name << "time_" << host << ".reg";              // If time regression
      else name << "accuracy.reg";                              // Else if accuracy regression
      std::fstream file;                                        // File id for regression
      file.open(name.str().c_str(),std::fstream::in);           //  Open regression file
      if (verbose) std::cout << "opening: " << name.str() << std::endl;// Print file name
      int numKeys;                                              // Number of keys stored in file
      if (file.good()) {                                        // If file exists
        if (verbose) std::cout << "file exists" << std::endl;   //  Print if file exists
        file >> numKeys;                                        //  Read number of keys
        for (int i=0; i<numKeys; i++) {                         //  Loop over regression values
          uint64_t readKey;                                     //   Read key buffer
          file >> readKey;                                      //   Read key
          file >> record[readKey];                              //   Read value
          if (secondValue) file >> record2[readKey];            //   Read value2
        }                                                       //  End loop over regression values
      }                                                         // End if for file existence
      file.close();                                             // Close regression file
      average = (average * iteration + value) / (iteration + 1);// Average of all iterations
      if (secondValue) average2 = (average2 * iteration + value2) / (iteration + 1);// Average of all iterations
      if (record[key] != 0 && verbose) std::cout << "entry exists" << std::endl;// Print if entry exits
      double threshold = record[key]*(1+EPS+iteration*.01);     // Accuracy threshold
      if ((record[key] == 0 || average <= threshold) && (average < 5e-4 || time)) { // If new record or pass
        pass = true;                                            //  Change flag to pass
        record[key] = average;                                  //  Add key value pair
      }                                                         // Endif for better value
      if (secondValue) {                                        // If second value
        threshold = record2[key]*(1+EPS+iteration*.01);         // Accuracy threshold
        if ((record2[key] == 0 || average2 <= threshold) && (average2 < 5e-3 || time)) { // If new record2 or pass
          pass &= true;                                         //  Change flag to pass
          record2[key] = average2;                              //  Add key value pair
        } else {                                                // Else if fail
          pass = false;                                         //  Change flag to fail
        }                                                       // Endif for second value
      }                                                         // Endif for better value
      file.open(name.str().c_str(),std::fstream::out);          // Open regression file
      file << record.size() << std::endl;                       // Write number of keys
      R_iter R2 = record.begin();                               // Dummy iterator for record2
      if (secondValue) R2 = record2.begin();                    // Iterator for record2
      for (R_iter R=record.begin(); R!=record.end(); R++,R2++) {// Loop over regression values
        file << R->first << " " << R->second;                   // Write key value pair
        if (secondValue) file << " " << R2->second;             // Write second value
        file << std::endl;                                      // End line
      }                                                         // End loop over regression values
      file.close();                                             // Close regression file
      if (!pass && verbose) {                                   // If regression failed
        if (time) std::cout << "Time regression failed: " <<    //  Print message for time regression
                    average << " / " << record[key] << std::endl;//  Print value and record
        else {                                                  // If accuracy regression
          std::cout << "Accuracy regression failed: " <<        //  Print message for accuracy regression
            average << " / " << record[key] << std::endl;       //  Print value and record
          if (secondValue) std::cout << "                            " << average2 //  Print value2
                                     << " / " << record2[key];  //  Print record2
          std::cout << std::endl;                               //  End line
        }                                                       // Endif accuracy regression
      }                                                         // Endif for failed regression
      return pass;                                              // Return flag for regression test
    }
  };
}
#endif
