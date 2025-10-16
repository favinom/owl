//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DeSaintVenant.h"

registerMooseObject("OwlApp", DeSaintVenant);

InputParameters
DeSaintVenant::validParams()
{
  InputParameters params = ElasticMaterial::validParams();
  params.addParam<Real>("mu", 1.0, "First Lamé constant");
  params.addParam<Real>("lambda", 1.0, "Second Lamé constant");
  return params;
}

DeSaintVenant::DeSaintVenant(const InputParameters & parameters) : 
ElasticMaterial(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda"))
{
}



void DeSaintVenant::computeQpProperties()
{
if(_Kmat){
	_K = _Kmat[_qp];
} 
else{
	_K = _I;
}
RealTensorValue C_K;
RealTensorValue B_K;

C_K = _K.transpose()*_K;
B_K = C_K.inverse();

_elasticMaterial[_qp] = this;
RealTensorValue U;
RealTensorValue F;
RealTensorValue C;
RealTensorValue E;
Real trE;
U(0,0) = _grad_disp_x[_qp](0); U(0,1) = _grad_disp_x[_qp](1); U(0,2) = _grad_disp_x[_qp](2);
U(1,0) = _grad_disp_y[_qp](0); U(1,1) = _grad_disp_y[_qp](1); U(1,2) = _grad_disp_y[_qp](2);
U(2,0) = _grad_disp_z[_qp](0); U(2,1) = _grad_disp_z[_qp](1); U(2,2) = _grad_disp_z[_qp](2);

F = _I+U;

C = F.transpose()*F;

E = (C - _I)/2.0;

trE = E.tr();

	//_stress[_qp] = 2*_mu*F*E + _lambda*trE*F;
	_stress[_qp] = _K.det()*F*(2*_mu*(B_K*C*B_K - B_K) +_lambda*((B_K*C).tr()-3)*B_K);

}





RealTensorValue DeSaintVenant::evaluateJac(RealTensorValue const & H, int const & qp){
if(_Kmat){
	_K = _Kmat[qp];
} 
else{
	_K = _I;
}
RealTensorValue C_K;
RealTensorValue B_K;

C_K = _K.transpose()*_K;
B_K = C_K.inverse();

RealTensorValue U;
RealTensorValue F;
RealTensorValue C;
RealTensorValue S;

U(0,0) = _grad_disp_x[qp](0); U(0,1) = _grad_disp_x[qp](1); U(0,2) = _grad_disp_x[qp](2);
U(1,0) = _grad_disp_y[qp](0); U(1,1) = _grad_disp_y[qp](1); U(1,2) = _grad_disp_y[qp](2);
U(2,0) = _grad_disp_z[qp](0); U(2,1) = _grad_disp_z[qp](1); U(2,2) = _grad_disp_z[qp](2);

F = _I+U;

C = F.transpose()*F;

S = 2*_K.det()*(_mu*(B_K*C*B_K - B_K) +_lambda/2.0*((B_K*C).tr()-3)*B_K);


RealTensorValue C_H;
RealTensorValue S_H;

C_H = F.transpose()*H + H.transpose()*F;

S_H = 2*_K.det()*(_mu*B_K*C_H*B_K+_lambda/2.0*(C_H*B_K).tr()*B_K);

	return H*S + F*S_H;
}









