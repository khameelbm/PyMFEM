%module plinearform
%{
#include <mpi.h>
#define MFEM_USE_MPI        
#include "fem/plinearform.hpp"
#include "numpy/arrayobject.h"
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%init %{
import_array();
%}

%import "cpointer.i"
%import "linearform.i"
%import "pfespace.i"
%import "hypre.i"

%pointer_class(int, intp);

#define MFEM_USE_MPI
%include "fem/plinearform.hpp"
