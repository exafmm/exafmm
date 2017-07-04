#ifndef buildtree_h
#define buildtree_h
#include "exafmm.h"

namespace exafmm {
  //! Get bounding box of bodies
  void getBounds(Bodies & bodies, real_t & R0, vec3 & X0) {
    vec3 Xmin = bodies[0].X;                                    // Initialize Xmin
    vec3 Xmax = bodies[0].X;                                    // Initialize Xmax
    for (size_t b=0; b<bodies.size(); b++) {                    // Loop over range of bodies
      Xmin = min(bodies[b].X, Xmin);                            //  Update Xmin
      Xmax = max(bodies[b].X, Xmax);                            //  Update Xmax
    }                                                           // End loop over range of bodies
    X0 = (Xmax + Xmin) / 2;                                     // Calculate center of domain
    R0 = fmax(max(X0-Xmin), max(Xmax-X0));                      // Calculate max distance from center
    R0 *= 1.00001;                                              // Add some leeway to radius
  }

  //! Build cells of tree adaptively using a top-down approach based on recursion
  void buildCells(Body * bodies, Body * buffer, int begin, int end, Cell * cell, Cells & cells,
                  const vec3 & X, real_t R, int level=0, bool direction=false) {
    // Create a tree cell
    cell->BODY = bodies + begin;                                // Pointer of first body in cell
    if(direction) cell->BODY = buffer + begin;                  // Pointer of first body in cell
    cell->NBODY = end - begin;                                  // Number of bodies in cell
    cell->NCHILD = 0;                                           // Initialize counter for child cells
    cell->X = X;                                                // Center position of cell
    cell->R = R / (1 << level);                                 // Cell radius
    // If cell is a leaf
    if (end - begin <= ncrit) {                                 // If number of bodies is less than threshold
      if (direction) {                                          //  If direction of data is from bodies to buffer
        for (int i=begin; i<end; i++) {                         //   Loop over bodies in cell
          buffer[i].X = bodies[i].X;                            //    Copy bodies coordinates to buffer
          buffer[i].q = bodies[i].q;                            //    Copy bodies source to buffer
        }                                                       //   End loop over bodies in cell
      }                                                         //  End if for direction of data
      return;                                                   //  Return without recursion
    }                                                           // End if for number of bodies
    // Count number of bodies in each octant
    int size[8] = {0,0,0,0,0,0,0,0};
    vec3 x;                                                     // Coordinates of bodies
    for (int i=begin; i<end; i++) {                             // Loop over bodies in cell
      x = bodies[i].X;                                          //  Position of body
      int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to
      size[octant]++;                                           //  Increment body count in octant
    }                                                           // End loop over bodies in cell
    // Exclusive scan to get offsets
    int offset = begin;                                         // Offset of first octant
    int offsets[8], counter[8];                                 // Offsets and counter for each octant
    for (int i=0; i<8; i++) {                                   // Loop over elements
      offsets[i] = offset;                                      //  Set value
      offset += size[i];                                        //  Increment offset
      if (size[i]) cell->NCHILD++;                              //  Increment child cell counter
    }                                                           // End loop over elements
    // Sort bodies by octant
    for (int i=0; i<8; i++) counter[i] = offsets[i];            // Copy offsets to counter
    for (int i=begin; i<end; i++) {                             // Loop over bodies
      x = bodies[i].X;                                          //  Position of body
      int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);// Which octant body belongs to`
      buffer[counter[octant]].X = bodies[i].X;                  // Permute bodies coordinates out-of-place according to octant
      buffer[counter[octant]].q = bodies[i].q;                  //  Permute bodies sources out-of-place according to octant
      counter[octant]++;                                        //  Increment body count in octant
    }                                                           // End loop over bodies
    // Loop over children and recurse
    vec3 Xchild;                                                // Coordinates of children
    cells.resize(cells.size()+cell->NCHILD);                    // Resize cell vector
    Cell * child = &cells.back() - cell->NCHILD + 1;            // Pointer for first child cell
    cell->CHILD = child;                                        // Point to first child cell
    int c = 0;                                                  // Counter for child cells
    for (int i=0; i<8; i++) {                                   // Loop over children
      Xchild = X;                                               //  Initialize center position of child cell
      real_t r = R / (1 << (level + 1));                        //  Radius of cells for child's level
      for (int d=0; d<3; d++) {                                 //  Loop over dimensions
        Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);         //   Shift center position to that of child cell
      }                                                         //  End loop over dimensions
      if (size[i]) {                                            //  If child exists
        buildCells(buffer, bodies, offsets[i], offsets[i] + size[i],// Recursive call for each child
                   &child[c], cells, Xchild, R, level+1, !direction);
        c++;                                                    //   Increment child cell counter
      }                                                         //  End if for child
    }                                                           // End loop over children
  }

  Cells buildTree(Bodies & bodies) {
    real_t R0;                                                  // Radius of root cell
    vec3 X0;                                                    // Center of root cell
    getBounds(bodies, R0, X0);                                  // Get bounding box from bodies
    Bodies buffer = bodies;                                     // Copy bodies to buffer
    Cells cells(1);                                             // Vector of cells
    cells.reserve(bodies.size());                               // Reserve memory space
    buildCells(&bodies[0], &buffer[0], 0, bodies.size(), &cells[0], cells, X0, R0);// Build tree recursively
    return cells;                                               // Return pointer of root cell
  }
}
#endif
