//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MyDiffusion.h"

registerMooseObject("OwlApp", MyDiffusion);

InputParameters
MyDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  	params.addParam<Real>("diff_coeff", 1.0, "Coefficient to multiply by the body force term");
  	params.addParam<Real>("react_coeff", 0.0, "Coefficient to multiply by the body force term");
  return params;
}

MyDiffusion::MyDiffusion(const InputParameters & parameters) : 
Kernel(parameters),
    	_D(getParam<Real>("diff_coeff")),
    	_R(getParam<Real>("react_coeff"))
         {}

Real
MyDiffusion::computeQpResidual()
{
  return _D*_grad_u[_qp] * _grad_test[_i][_qp]+ _R*_u[_qp]*_test[_i][_qp];
}

Real
MyDiffusion::computeQpJacobian()
{
  return _D*_grad_phi[_j][_qp] * _grad_test[_i][_qp] + _R*_phi[_j][_qp]*_test[_i][_qp];
}
