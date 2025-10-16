//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MyTimeDerivative.h"

registerMooseObject("OwlApp", MyTimeDerivative);

InputParameters
MyTimeDerivative::validParams()
{
  InputParameters params = TimeKernel::validParams();
  return params;
}

MyTimeDerivative::MyTimeDerivative(const InputParameters & parameters) : TimeKernel(parameters) {}

Real
MyTimeDerivative::computeQpResidual()
{
  return  _u_dot[_qp]*_test[_i][_qp];   //_grad_u[_qp] * _grad_test[_i][_qp]+ _u[_qp]*_test[_i][_qp];
}

Real
MyTimeDerivative::computeQpJacobian()
{
  return _phi[_j][_qp]/_dt * _test[_i][_qp];
}
