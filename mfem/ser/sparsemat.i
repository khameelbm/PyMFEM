%module sparsemat

%{
#include "linalg/sparsemat.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"
#include <fstream>  
#include <iostream>  
  %}
// initialization required to return numpy array from SWIG
%init %{
import_array();
%}
%import "array.i"
%import "vector.i"
%import "operators.i"
%import "matrix.i"
%import "densemat.i"
%ignore Walk;

%pythonappend mfem::SparseMatrix::operator*= %{
    val.thisown = self.thisown
%}
%pythonappend mfem::SparseMatrix::operator+= %{
    val.thisown = self.thisown
%}

// RAP_P, RAP_R replaces RAP, since RAP has two definition one accept
// pointer and the other accept reference. From Python, two
// can not be distingished..
%inline %{
  mfem::SparseMatrix *RAP_P (const mfem::SparseMatrix &A,
			     const mfem::SparseMatrix &R,
                             mfem::SparseMatrix *ORAP){
    return RAP(A, R, ORAP);
  }    
    
  mfem::SparseMatrix *RAP_R(const mfem::SparseMatrix &Rt,
			    const mfem::SparseMatrix &A,
                            const mfem::SparseMatrix &P){
    return RAP(Rt, A, P);
  }
%}
%rename(add_sparse) mfem::Add;
%include "linalg/sparsemat.hpp"

%extend mfem::SparseMatrix {
PyObject* GetIArray(void) const{
  int * I = self->GetI();
  int L = self->Size();

  npy_intp dims[] = { L+1 };
  return  PyArray_SimpleNewFromData(1, dims, NPY_INT, I);
  }
PyObject* GetJArray(void) const{
  int * I = self->GetI();
  int * J = self->GetJ();
  int L = self->Size();
  npy_intp dims[] = { I[L]};
  return  PyArray_SimpleNewFromData(1, dims, NPY_INT, J);
  }
PyObject* GetDataArray(void) const{
  int * I = self->GetI();
  double * A = self->GetData();    
  int L = self->Size();
  npy_intp dims[] = {I[L]};
  return  PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, A);
  }
void Print(const char *file){
  std::ofstream mat_stream(file);
  if (!mat_stream) {
     std::cerr << "\nCan not open output file: " << file << '\n' << std::endl;
     return;
     }
  self->Print(mat_stream);
  } 
};
