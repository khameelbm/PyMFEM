%module(directors="1")  bilinearform
%{
#include "fem/bilinearform.hpp"
#include "numpy/arrayobject.h"
#include "pyoperator.hpp"         
%}

%init %{
import_array();
%}

%import "array.i"
%import "fespace.i"
%import "fe_coll.i"
%import "intrules.i"
%import "matrix.i"
%import "vector.i"
%import "densemat.i"
%import "sparsemat.i"
%import "lininteg.i"
%import "eltrans.i"
%import "bilininteg.i"
%import "linearform.i"
%import "gridfunc.i"

%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }    
}

 //%include "fem/coefficient.hpp"
namespace mfem { 
%pythonprepend BilinearForm::AddDomainIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0 
   %}
%pythonprepend BilinearForm::AddBoundaryIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0 
   %} 
%pythonprepend BilinearForm::AddBdrFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0 
   %} 
%pythonprepend BilinearForm::AddInteriorFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0 
   %}
%pythonprepend MixedBilinearForm::AddDomainIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0 
   %}
%pythonprepend MixedBilinearForm::AddBoundaryIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0 
   %} 
%pythonprepend MixedBilinearForm::AddTraceFaceIntegrator %{
    if not hasattr(self, "_integrators"): self._integrators = []
    self._integrators.append(bfi)
    bfi.thisown=0 
   %} 
%pythonappend MixedBilinearForm::SpMat %{
    if not hasattr(self, "_spmat"): self._spmat = []
    self._spmat.append(val)
    val.thisown=0 
   %}
%pythonappend BilinearForm::SpMat %{
    if not hasattr(self, "_spmat"): self._spmat = []
    self._spmat.append(val)
    val.thisown=0 
   %}
%pythonappend BilinearForm::EnableHybridization %{
    constr_integ.thisown = 0
   %} 
} 




%include "fem/bilinearform.hpp"
