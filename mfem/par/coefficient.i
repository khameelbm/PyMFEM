/*
   Coefficient.i
   SWIG interface file for coefficient.hpp

   Features
   1) function callback for VectorFunctionCoefficent


   Date: 2016. 2. 18
   Author: S. Shiraiwa (MIT)
 */
%module(directors="1")  coefficient
/*%module  coefficient*/
%{
#define MFEM_USE_MPI    
#include "fem/fem.hpp"
#include "fem/fe_coll.hpp"
#include "fem/fespace.hpp"
#include "fem/eltrans.hpp"
#include "fem/intrules.hpp"
#include "linalg/vector.hpp"
#include "mesh/pmesh.hpp"      
#include "fem/coefficient.hpp"
#include "linalg/densemat.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include <cstring>
#include <ctime>
#include "numpy/arrayobject.h"  
#include "pycoefficient.hpp"
#define MFEM_USE_MPI  
%}
%init %{
import_array();
%}
#define MFEM_USE_MPI    
//%import "general/array.hpp"
%import "array.i"
%import "matrix.i"
%import "intrules.i"
%import "sparsemat.i"
%import "densemat.i"
%import "vector.i"
%import "eltrans.i"
//%import "pmesh.i"
%ignore Function;
%ignore DeltaCoefficient;
%feature("notabstract") VectorFunctionCoefficient;
%feature("notabstract") VectorConstantCoefficient;

namespace mfem { 
%pythonprepend DeltaCoefficient::SetWeight %{
    w.thisown=0 
%}
%pythonprepend VectorArrayCoefficient::Set %{ 
    c.thisown=0 
%}
%pythonprepend MatrixArrayCoefficient::Set %{ 
    c.thisown=0 
%}
%pythonappend VectorRestrictedCoefficient::VectorRestrictedCoefficient %{
    self._ref_to_vc = vc
%}
%pythonappend RestrictedCoefficient::RestrictedCoefficient %{
    self._ref_to_c = _c
%}
%pythonappend MatrixRestrictedCoefficient::MatrixRestrictedCoefficient %{
    self._ref_to_mc = mc
%}
}

%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }    
    //catch (...){
    //  SWIG_fail;
    //}
    //    catch (Swig::DirectorMethodException &e) { SWIG_fail; }
    //    catch (std::exception &e) { SWIG_fail; }    
}
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}

%typemap(in) const mfem::IntegrationRule *irs[]{
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    //std::cout << std::to_string(size) << "\n";
    int i = 0;
    $1 = (mfem::IntegrationRule **) malloc((size+1)*sizeof(mfem::IntegrationRule *));
    for (i = 0; i < size; i++) {
       //std::cout << std::to_string(i) << "\n";      
       PyObject *o = PyList_GetItem($input,i);
       void *temp;
       if (SWIG_ConvertPtr(o, &temp,
	   $descriptor(mfem::IntegrationRule *),SWIG_POINTER_EXCEPTION) == -1){
   	   //std::cout << "Failed to pointer conversion" << "\n";      	 
           return NULL;
       }
       $1[i] = reinterpret_cast<mfem::IntegrationRule *>(temp);
     }
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}
%typemap(typecheck) const mfem::IntegrationRule *irs[]{
   $1 = PyList_Check($input) ? 1 : 0;
}

%include "fem/coefficient.hpp"

%feature("director") mfem::VectorPyCoefficientBase;
%feature("director") mfem::PyCoefficientBase;
%feature("director") mfem::MatrixPyCoefficientBase;

%inline %{
double fake_func(const mfem::Vector &x)
{
  return 0.0;
}
void fake_func_vec(const mfem::Vector &x, mfem::Vector &Ht)
{
     Ht(0) = 0.0;
     Ht(1) = 0.0;
     Ht(2) = 0.0;
}
void fake_func_mat(const mfem::Vector &x, mfem::DenseMatrix &Kt)
{
  Kt(0,0) = 1.0;
  Kt(1,0) = 0.0;
  Kt(2,0) = 0.0;
  Kt(0,1) = 0.0;
  Kt(1,1) = 1.0;
  Kt(2,1) = 0.0;
  Kt(0,2) = 0.0;
  Kt(1,2) = 0.0;
  Kt(2,2) = 1.0;
}

namespace mfem{ 
double PyCoefficientBase::Eval(ElementTransformation &T,
                               const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);

   if (isTimeDependent)
   {
     return _EvalPyT(transip, GetTime());
   }
   else
   {
     return _EvalPy(transip);
   }
}
void VectorPyCoefficientBase::Eval(Vector &V, ElementTransformation &T,
                                     const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);

   V.SetSize(vdim);
   if (isTimeDependent)
   {
      _EvalPyT(transip, GetTime(),  V);          
   }
   else
   {
      _EvalPy(transip, V);
   }
}

void VectorPyCoefficientBase::Eval(DenseMatrix &M, ElementTransformation &T,
                                  const IntegrationRule &ir)

{
   Vector Mi;
   M.SetSize(vdim, ir.GetNPoints());
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      M.GetColumnReference(i, Mi);
      const IntegrationPoint &ip = ir.IntPoint(i);
      T.SetIntPoint(&ip);
      Eval(Mi, T, ip);
   }
}
void MatrixPyCoefficientBase::Eval(DenseMatrix &K, ElementTransformation &T,
                                     const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);

   T.Transform(ip, transip);
   K.SetSize(vdim);
   if (isTimeDependent)
   {
      _EvalPyT(transip, GetTime(),  K);          
   }
   else
   {
      _EvalPy(transip, K);
   }
}

}  /* end of name space*/
%}

%include "pycoefficient.hpp"

%pythoncode %{
class PyCoefficient(PyCoefficientBase):
   def __init__(self):
       PyCoefficientBase.__init__(self, 0)
   def _EvalPy(self, x):
       return self.EvalValue(x.GetDataArray())
   def EvalValue(self, x):
       return 0.0
  
class PyCoefficientT(PyCoefficientBase):
   def __init__(self):
       PyCoefficientBase.__init__(self, 1)
   def _EvalPyT(self, x, t, V):
       return self.EvalValue(x.GetDataArray(), t)
   def EvalValue(self, x, t):
       return 0.0
	 
class VectorPyCoefficient(VectorPyCoefficientBase):
   def __init__(self, dim):
       self.sdim = dim
       VectorPyCoefficientBase.__init__(self, dim, 0)
   def _EvalPy(self, x, V):
       v = self.EvalValue(x.GetDataArray())
       for i in range(self.sdim):
           V[i] = v[i]
   def _EvalPyT(self, x, t, V):
       v = self.EvalValue(x.GetDataArray())  
       for i in range(self.sdim):
           V[i] = v[i]
   def EvalValue(self, x):
       return [0,0,0]
  
class VectorPyCoefficientT(VectorPyCoefficientBase):
   def __init__(self, dim):
       self.sdim = dim  
       VectorPyCoefficientBase.__init__(self, dim, 1)
   def _EvalPy(self, x, V):
       v = self.EvalValue(x.GetDataArray(), 0)
       for i in range(self.sdim):
           V[i] = v[i]
   def _EvalPyT(self, x, t, V):
       v = self.EvalValue(x.GetDataArray(), t)  
       for i in range(self.sdim):
           V[i] = v[i]
   def EvalValue(self, x, t):
       return [0,0,0]

class MatrixPyCoefficient(MatrixPyCoefficientBase):
   def __init__(self, dim):
       self.sdim = dim
       MatrixPyCoefficientBase.__init__(self, dim, 0)
   def _EvalPy(self, x, K):
       k = self.EvalValue(x.GetDataArray())
       for i in range(self.sdim):
           for j in range(self.sdim):
               K[i, j] = k[i, j]
   def EvalValue(self, x):
       return np.array([[0,0,0], [0,0,0] [0,0,0]])
  
class MatrixPyCoefficientT(MatrixPyCoefficientBase):
   def __init__(self, dim):
       self.sdim = dim  
       MatrixPyCoefficientBase.__init__(self, dim, 1)
   def _EvalPyT(self, x, t, K):
       k = self.EvalValue(x.GetDataArray(), t)  
       for i in range(self.sdim):
           for j in range(self.sdim):
               K[i, j] = k[i, j]
   def EvalValue(self, x, t):
       return np.array([[0,0,0], [0,0,0] [0,0,0]])
	 
%}

