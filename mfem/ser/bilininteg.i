%module(directors="1")  bilininteg
%{
#include "fem/gridfunc.hpp"  
#include "fem/linearform.hpp"
#include "fem/bilininteg.hpp"
#include "pycoefficient.hpp"
#include "numpy/arrayobject.h"      
%}

%init %{
import_array();
%}

%import "array.i"
%import "coefficient.i"
%import "matrix.i"
%import "vector.i"
%import "gridfunc.i"
%import "fespace.i"
%import "fe_coll.i"
%import "intrules.i"
%import "densemat.i"
%import "sparsemat.i"
%import "lininteg.i"
%import "eltrans.i"
%import "linearform.i"
%import "fe.i"
 //%template(IntegrationPointArray) mfem::Array<mfem::IntegrationPoint>;

%exception {
    try { $action }
    catch (Swig::DirectorException &e) { SWIG_fail; }    
}
namespace mfem {
%pythonprepend TransposeIntegrator::TransposeIntegrator %{
    if _own_bfi == 1:  _bfi.thisown = 0
%}
%pythonprepend InverseIntegrator::InverseIntegrator %{
    if own_integ == 1:  integ.thisown = 0
%}
%pythonprepend SumIntegrator::AddIntegrator %{
    integ.thisown = 0
%}
%pythonappend CurlCurlIntegrator::CurlCurlIntegrator %{
    self._coeff = args[0]
%}
%pythonappend VectorFEMassIntegrator::VectorFEMassIntegrator %{
    self._coeff = args[0]
%}
}

%include "fem/bilininteg.hpp"
%inline %{
namespace mfem{
  
class GradientVectorFEIntegrator : public BilinearFormIntegrator
{
private:
   Coefficient *Q;
   VectorCoefficient *VQ;
   MatrixCoefficient *MQ;

   DenseMatrix Jinv;
   DenseMatrix dshape;
   DenseMatrix gshape;
   DenseMatrix elmat;
#ifndef MFEM_THREAD_SAFE
   Vector shape;
   Vector D;
   DenseMatrix K;
   DenseMatrix test_vshape;
   DenseMatrix trial_vshape;
#endif
   

public:
   GradientVectorFEIntegrator() { Q = NULL; VQ = NULL; MQ = NULL;}
   GradientVectorFEIntegrator(Coefficient &q) { Q = &q; VQ = NULL; MQ = NULL;}
   virtual ~GradientVectorFEIntegrator() { }
   virtual void AssembleElementMatrix2(const FiniteElement &trial_el,
                                       const FiniteElement &test_el,
                                      ElementTransformation &Trans,
				       DenseMatrix &elmat){
   // test  H1
   // trial ND
   int dim = test_el.GetDim();
   int trial_dof = trial_el.GetDof();
   int test_dof = test_el.GetDof();
   int spaceDim = Trans.GetSpaceDim();   
   double w;
   
   //std::cout << "dim :" << std::to_string(dim) << "\n";
   //std::cout << "space_dim :" << std::to_string(spaceDim) << "\n";         
   //std::cout << "test_dof:" << std::to_string(test_dof) << "\n";
   //std::cout << "trial_dof:" << std::to_string(trial_dof) << "\n";   

   elmat.SetSize(test_dof, trial_dof);
   Jinv.  SetSize(dim);
   dshape.SetSize(test_dof, dim);
   gshape.SetSize(test_dof, dim);
#ifdef MFEM_THREAD_SAFE
   Vector D(VQ ? VQ->GetVDim() : 0);
   DenseMatrix trial_vshape(trial_dof, spaceDim);
   DenseMatrix K(MQ ? MQ->GetVDim() : 0, MQ ? MQ->GetVDim() : 0);
#else
   trial_vshape.SetSize(trial_dof, spaceDim);
   D.SetSize(VQ ? VQ->GetVDim() : 0);
   K.SetSize(MQ ? MQ->GetVDim() : 0, MQ ? MQ->GetVDim() : 0);
#endif

   //std::cout << "trial_vshape H:" << std::to_string(trial_vshape.Height()) << "\n";
   //std::cout << "K width:" << std::to_string(K.Width()) << "\n";   
   DenseMatrix tmp(trial_vshape.Height(), K.Width());
   DenseMatrix tmp2(test_dof, trial_dof);
   
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // integrand is rational function if det(J) is not constant
      int order = Trans.OrderGrad(&test_el) + trial_el.GetOrder();
      if (test_el.Space() == FunctionSpace::rQk)
      {
         ir = &RefinedIntRules.Get(test_el.GetGeomType(), order);
      }
      else
      {
         ir = &IntRules.Get(test_el.GetGeomType(), order);
      }
   }
   for (int i = 0; i < ir -> GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      test_el.CalcDShape (ip, dshape);
      Trans.SetIntPoint (&ip);     
      CalcInverse (Trans.Jacobian(), Jinv);
      Mult (dshape, Jinv, gshape);
      Trans.SetIntPoint (&ip);

      trial_el.CalcVShape(Trans, trial_vshape);
      w = ip.weight * Trans.Weight();
      if (MQ)
      {
         MQ->Eval(K, Trans, ip);
         K *= w;
         Mult(gshape, K, tmp);
         AddMultABt(tmp,trial_vshape,elmat);
      }
      else if (VQ)
      {
         VQ->Eval(D, Trans, ip);
         D *= w;
         MultADBt(gshape, D, trial_vshape, tmp2);	 
         elmat += tmp2;
      }
      else
      {
         if (Q)
         {
            w *= Q -> Eval (Trans, ip);
         }
         gshape *= w;
         AddMultABt (gshape, trial_vshape, elmat);
      }
   }
   };
};   
}
%}
