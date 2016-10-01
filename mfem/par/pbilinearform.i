%module pbilinearform
%{
#include <mpi.h>
#define MFEM_USE_MPI    
#include "fem/pbilinearform.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"           

%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);


%init %{
import_array();
%}

%import "cpointer.i"
%import "bilinearform.i"
%import "pfespace.i"
%import "hypre.i"

%pointer_class(int, intp);

#define MFEM_USE_MPI  
%include "fem/pbilinearform.hpp"
