%module pmesh
%{
#include <mpi.h>
#define MFEM_USE_MPI
#include "mesh/mesh.hpp"  
#include "mesh/pmesh.hpp"
#include "general/communication.hpp"  
#include "numpy/arrayobject.h"
  //#include "mpi4py/mpi4py.h"  
%}

%init %{
import_array();
%}

%include mpi4py/mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);

%import "cpointer.i"
%import "mesh.i"
%import "pncmesh.i"
%import "communication.i"

%pointer_class(int, intp);

%immutable face_nbr_elements;
%immutable face_nbr_vertices;
%immutable gtopo;


#define MFEM_USE_MPI  
%include "mesh/pmesh.hpp"
