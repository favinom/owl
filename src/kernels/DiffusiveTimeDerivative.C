//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusiveTimeDerivative.h"

registerMooseObject("OwlApp", DiffusiveTimeDerivative);

InputParameters
DiffusiveTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
    params.addParam<Real>("m", 1.0, "First Lamé constant");
  params.addParam<Real>("kv", 1.0, "Second Lamé constant");
    params.addRequiredCoupledVar("J_k", "determinante di K, parametro d'ordine per questo modello");
        params.addRequiredCoupledVar("mu_k", "potenziale chimico");
  return params;
}

DiffusiveTimeDerivative::DiffusiveTimeDerivative(const InputParameters & parameters) : 
TimeKernel(parameters),
_m(getParam<Real>("m")),
_kv(getParam<Real>("kv")),
_id_Jk(coupled("J_k")),
_id_muk(coupled("mu_k")),
_Jk(coupledValue("J_k")),
_Jk_old(coupledValueOld("J_k")),
_grad_Jk(coupledGradient("J_k")),
_grad_Jk_old(coupledGradientOld("J_k")),
_muk(coupledValue("mu_k"))
{}

Real DiffusiveTimeDerivative::computeQpResidual()
{
  return  _m*_kv*(_grad_Jk[_qp]/_Jk[_qp] - _grad_Jk_old[_qp]/_Jk_old[_qp] )/_dt*_grad_test[_i][_qp];   //_grad_u[_qp] * _grad_test[_i][_qp]+ _u[_qp]*_test[_i][_qp];
}

Real DiffusiveTimeDerivative::computeQpJacobian()
{
  return  _m*_kv*(_Jk[_qp]*_grad_phi[_j][_qp] - _phi[_j][_qp]*_grad_Jk[_qp])/(_dt*_Jk[_qp]*_Jk[_qp])* _grad_test[_i][_qp];
}


Real DiffusiveTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar){

return 0.0;

}
