//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearElasticityMaterial.h"

registerMooseObject("OwlApp", LinearElasticityMaterial);

InputParameters
LinearElasticityMaterial::validParams()
{
  InputParameters params = ElasticMaterial::validParams();
  params.addParam<Real>("mu", 1.0, "First Lamé constant");
  params.addParam<Real>("lambda", 1.0, "Second Lamé constant");
  return params;
}

LinearElasticityMaterial::LinearElasticityMaterial(const InputParameters & parameters) : 
ElasticMaterial(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda"))
{
}


void
LinearElasticityMaterial::computeQpProperties()
{
_elasticMaterial[_qp] = this;
RealTensorValue U;
RealTensorValue eps;
RealTensorValue sigma;
Real treps;
U(0,0) = _grad_disp_x[_qp](0); U(0,1) = _grad_disp_x[_qp](1); U(0,2) = _grad_disp_x[_qp](2);
U(1,0) = _grad_disp_y[_qp](0); U(1,1) = _grad_disp_y[_qp](1); U(1,2) = _grad_disp_y[_qp](2);
U(2,0) = _grad_disp_z[_qp](0); U(2,1) = _grad_disp_z[_qp](1); U(2,2) = _grad_disp_z[_qp](2);
eps = (U+U.transpose())/2.0;

treps = eps.tr();

	_stress[_qp] = 2*_mu*eps + _lambda*treps*_I;
}





RealTensorValue LinearElasticityMaterial::evaluateJac(RealTensorValue const & H, int const & qp){
RealTensorValue epsH;
RealTensorValue sigmaH;
Real trepsH;

epsH = (H+H.transpose())/2.0;

trepsH = epsH.tr();

sigmaH = 2*_mu*epsH + _lambda*trepsH*_I;


	return sigmaH;
}





