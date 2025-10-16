//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SAC_Diffusion.h"

registerMooseObject("OwlApp", SAC_Diffusion);

InputParameters
SAC_Diffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  	params.addParam<Real>("alpha", 1.0, "Coefficient to multiply by the body force term");
    	params.addParam<Real>("beta", 1.0, "Coefficient to multiply by the body force term");
    	
  return params;
}

SAC_Diffusion::SAC_Diffusion(const InputParameters & parameters) : 
Kernel(parameters),
    	_alpha(getParam<Real>("alpha")),
        _beta(getParam<Real>("beta"))
         {
         _min1 = 0.2;
         _max = 0.3;
         _min2 = 0.8;
         _A(0,0) = 1;
         _A(0,1) = 0;
         _A(1,1) = 1; 
         _A(1,0) = 0;
         }
//# (u-0.2)^2*(u-0.8)^2


Real
SAC_Diffusion::computeQpResidual()
{
Real r = (_u[_qp]-_min1)*(_u[_qp]-_max)*(_u[_qp]-_min2);

  return _alpha*_grad_u[_qp] * _grad_test[_i][_qp]+_beta*r*_test[_i][_qp];
}

Real
SAC_Diffusion::computeQpJacobian()
{
Real r1 = (_u[_qp]-_min1)*(_u[_qp]-_max);
Real r2 = (_u[_qp]-_max)*(_u[_qp]-_min2);
Real r3 = (_u[_qp]-_min1)*(_u[_qp]-_min2);
Real r = r1 + r2 + r3;
  return _alpha*_grad_phi[_j][_qp] * _grad_test[_i][_qp]+_beta*r*_phi[_j][_qp]*_test[_i][_qp];
}
