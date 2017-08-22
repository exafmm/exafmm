#ifndef buildtree_h
#define buildtree_h
#include "exafmm.h"

namespace exafmm {
  //! Get bounding box of bodies
  void getBounds(Bodies & bodies, real_t & R0, vec2 & X0) {
    vec2 Xmin = bodies[0].X;
    vec2 Xmax = bodies[0].X;
    for (size_t b=0; b<bodies.size(); b++) {
      Xmin = min(bodies[b].X, Xmin);
      Xmax = max(bodies[b].X, Xmax);
    }
    X0 = (Xmax + Xmin) / 2;
    R0 = fmax(max(X0-Xmin), max(Xmax-X0));
    R0 *= 1.00001;
  }

  //! Build cells of tree adaptively using a top-down approach based on recursion
  void buildCells(Body * bodies, Body * buffer, int begin, int end, Cell * cell, Cells & cells,
                  const vec2 & X, real_t R, int level=0, bool direction=false) {
    //! Create a tree cell
    cell->body = bodies + begin;
    if(direction) cell->body = buffer + begin;
    cell->numBodies = end - begin;
    cell->numChilds = 0;
    cell->X = X;
    cell->R = R / (1 << level);
    //! If cell is a leaf
    if (end - begin <= NCRIT) {
      if (direction) {
        for (int i=begin; i<end; i++) {
          buffer[i].X = bodies[i].X;
          buffer[i].q = bodies[i].q;
        }
      }
      return;
    }
    //! Count number of bodies in each quadrant
    int size[4] = {0,0,0,0};
    vec2 x;
    for (int i=begin; i<end; i++) {
      x = bodies[i].X;
      int quadrant = (x[0] > X[0]) + ((x[1] > X[1]) << 1);
      size[quadrant]++;
    }
    //! Exclusive scan to get offsets
    int offset = begin;
    int offsets[4], counter[4];
    for (int i=0; i<4; i++) {
      offsets[i] = offset;
      offset += size[i];
      if (size[i]) cell->numChilds++;
    }
    //! Sort bodies by quadrant
    for (int i=0; i<4; i++) counter[i] = offsets[i];
    for (int i=begin; i<end; i++) {
      x = bodies[i].X;
      int quadrant = (x[0] > X[0]) + ((x[1] > X[1]) << 1);
      buffer[counter[quadrant]].X = bodies[i].X;
      buffer[counter[quadrant]].q = bodies[i].q;
      counter[quadrant]++;
    }
    //! Loop over children and recurse
    vec2 Xchild;
    cells.resize(cells.size()+cell->numChilds);
    Cell * child = &cells.back() - cell->numChilds + 1;
    cell->child = child;
    int c = 0;
    for (int i=0; i<4; i++) {
      Xchild = X;
      real_t r = R / (1 << (level + 1));
      for (int d=0; d<2; d++) {
        Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);
      }
      if (size[i]) {
        buildCells(buffer, bodies, offsets[i], offsets[i] + size[i],
                   &child[c], cells, Xchild, R, level+1, !direction);
        c++;
      }
    }
  }

  Cells buildTree(Bodies & bodies) {
    real_t R0;
    vec2 X0;
    getBounds(bodies, R0, X0);
    Bodies buffer = bodies;
    Cells cells(1);
    cells.reserve(bodies.size());
    buildCells(&bodies[0], &buffer[0], 0, bodies.size(), &cells[0], cells, X0, R0);
    return cells;
  }
}
#endif
