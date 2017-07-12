## Development Guide

**Flow of development:**

1. Experiment with design decisions in 2d/

2. Propagate to other directories when decision is reached

3. Write incremental tests and leave breadcrumbs

4. Document each new code

5. Unit test for each component

6. Integration/regression test for each commit


**Variable name convention:**

Use all capital for
 - Macro variables
 - Environtment values
 - Makefile variables
 - Global variables

Start with capital for
 - Class & Struct names
 - User defined types

Start with small case for
 - everything else

Use camel case for variable names, underscore for file names
Avoid abbreviation as much as possible
 - Short scope, short name
 - Use plural for vectors (e.g. vector<cell> cells;)
 - For equations use the same characters in the equations

Common characters (prefix)
 - d* : Difference
 - i* : Loop counter (don't just use i)
 - k* : Waves
 - n* : Size (for inner kernels)
 - X* : Coordinates

Common characters (postfix)
 - *i : target
 - *j : source
 - *p : parent
 - *c : child
 - *0 : first element
 - *N : last element

Common three letter keywords (prefix)
 - max* : Maximum
 - num* : Size (use "n" for inner kerenls)
 - old* : Previous value
 - sum* : Sum
 - tmp* : Temporary value
 - vec* : Short vector type

Common keywords (postfix)
 - *begin & *end (following STL)
 - *send & *recv (following MPI)
 - *displ & *count (following MPI)
 - *global & *local (for distributed memory)
 - *min & *max
 - *key & *color

GPU memory keywords (postfix)
 - *glob : Global memory
 - *shrd : Shared memory
