/**
 Integrator for (grad u, Qv) 
    u : nodal element
    v : nd element
**/
#include "fem/bilininteg.hpp"

class GradientVectorFEIntegrator : public BilinearFormIntegrator
{
private:
   Coefficient *Q;

   DenseMatrix Jinv;
   DenseMatrix dshape;
   DenseMatrix gshape;
   DenseMatrix pelmat;

public:
   GradientVectorFEIntegrator() { Q = NULL; }
   GradientVecotrFEIntegrator(Coefficient &q) { Q = &q; }

   virtual void AssembleElementMatrix2(const FiniteElement &trial_el,
                                       const FiniteElement &test_fe,
                                      ElementTransformation &Trans,
				       DenseMatrix &elmat){
   // test  H1
   // trial ND
   int dim = test_el.GetDim();
   int trial_dof = trial_el.GetDof();
   int test_dof = test_el.GetDof();
   double norm;

   elmat.SetSize (test_dof, trial_dof);

   Jinv.  SetSize (dim);
   dshape.SetSize (test_dof, dim);
   gshape.SetSize (test_dof, dim);
     
   const IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // integrand is rational function if det(J) is not constant
      int order = Trans.OrderGrad(&test_el) + trail_fe.GetOrder()
      if (el.Space() == FunctionSpace::rQk)
      {
         ir = &RefinedIntRules.Get(el.GetGeomType(), order);
      }
      else
      {
         ir = &IntRules.Get(el.GetGeomType(), order);
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
         AddMultADBt(gshape, D, trial_vshape, elmat);
      }
      else
      {
         if (Q)
         {
            w *= Q -> Eval (Trans, ip);
         }
         gshape *= w;
         AddMult_a_ABt (gshape, trial_vshape, elmat);
      }
   }
};
