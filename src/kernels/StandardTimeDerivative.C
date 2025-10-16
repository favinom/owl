//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StandardTimeDerivative.h"

registerMooseObject("OwlApp", StandardTimeDerivative);

InputParameters
StandardTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
    params.addRequiredCoupledVar("J_k", "determinante di K, parametro d'ordine per questo modello");
    params.addRequiredCoupledVar("mu_k", "potenziale chimico");
  return params;
}

StandardTimeDerivative::StandardTimeDerivative(const InputParameters & parameters) : 
TimeKernel(parameters),
_id_Jk(coupled("J_k")),
_id_muk(coupled("mu_k")),
_Jk(coupledValue("J_k")),
_Jk_old(coupledValueOld("J_k")),
_muk(coupledValue("mu_k"))
{}

Real StandardTimeDerivative::computeQpResidual()
{
  return  (_Jk[_qp]-_Jk_old[_qp])/_dt*_test[_i][_qp];   //_grad_u[_qp] * _grad_test[_i][_qp]+ _u[_qp]*_test[_i][_qp];
}

Real StandardTimeDerivative::computeQpJacobian()
{
  return _phi[_j][_qp]/_dt * _test[_i][_qp];
}


Real StandardTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar){

return 0.0;

}
