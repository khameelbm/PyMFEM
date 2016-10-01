%module mesh_operators

%{
#define MFEM_USE_MPI      
#include "numpy/arrayobject.h"
#include "mesh/mesh_operators.hpp"
#include "mesh/mesh_operators.hpp"    
  %}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}

%import "array.i"
%import "vector.i"
%import "mesh.i"
%import "estimators.i"

%include "mesh/mesh_operators.hpp"  
