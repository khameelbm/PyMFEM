%module pncmesh
%{
#include <mpi.h>
#define MFEM_USE_MPI    
#include "mesh/mesh_headers.hpp"
#include "mpi4py/mpi4py.h"
#include "numpy/arrayobject.h"  
%}
%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);
#define MFEM_USE_MPI
%init %{
import_array();
%}

%init %{
import_array();
%}

%import "cpointer.i"
%import "mesh.i"
%import "ncmesh.i"
%import "communication.i"

%pointer_class(int, intp);

#define MFEM_USE_MPI  
%include "mesh/pncmesh.hpp"
