%module pnonlinearform
%{
#include <mpi.h>
#define MFEM_USE_MPI
#include "fem/pnonlinearform.hpp"
#include "fem/linearform.hpp"
#include "pyoperator.hpp"           
%}

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
/*
%init %{
import_array();
%}
*/

%import "cpointer.i"
%import "nonlinearform.i"
%import "pfespace.i"
%import "pgridfunc.i"

%pointer_class(int, intp);

#define MFEM_USE_MPI
%include "fem/pnonlinearform.hpp"
