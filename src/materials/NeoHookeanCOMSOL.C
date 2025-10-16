//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NeoHookeanCOMSOL.h"

registerMooseObject("OwlApp", NeoHookeanCOMSOL);

InputParameters
NeoHookeanCOMSOL::validParams()
{
  InputParameters params = ElasticMaterial::validParams();
  params.addParam<Real>("mu", 1.0, "First Lamé constant");
  params.addParam<Real>("lambda", 1.0, "Second Lamé constant");
  return params;
}

NeoHookeanCOMSOL::NeoHookeanCOMSOL(const InputParameters & parameters) : 
ElasticMaterial(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda"))
{
}



void NeoHookeanCOMSOL::computeQpProperties()
{
if(_Kmat){
	_K = _Kmat[_qp];
} 
else{
	_K = _I;
}
RealTensorValue C_K;
RealTensorValue B_K;
Real J_K;

C_K = _K.transpose()*_K;
B_K = C_K.inverse();
J_K = _K.det();

_elasticMaterial[_qp] = this;
RealTensorValue U;
RealTensorValue F;
RealTensorValue C;
RealTensorValue FinvT;
Real trC;
Real J;
Real c1;
Real d1;

U(0,0) = _grad_disp_x[_qp](0); U(0,1) = _grad_disp_x[_qp](1); U(0,2) = _grad_disp_x[_qp](2);
U(1,0) = _grad_disp_y[_qp](0); U(1,1) = _grad_disp_y[_qp](1); U(1,2) = _grad_disp_y[_qp](2);
U(2,0) = _grad_disp_z[_qp](0); U(2,1) = _grad_disp_z[_qp](1); U(2,2) = _grad_disp_z[_qp](2);

F = _I+U;

J = F.det();

FinvT = F.inverse().transpose();

C = F.transpose()*F;

trC = C.tr();

c1 = _mu/2.0;
d1 = _lambda/2.0;

	_stress[_qp] = _mu*J_K*F*B_K  + (_lambda*J*(log(J) - log(J_K)) - _mu*J_K )*FinvT;
}




RealTensorValue NeoHookeanCOMSOL::evaluateJac(RealTensorValue const & H, int const & qp){
if(_Kmat){
	_K = _Kmat[qp];
} 
else{
	_K = _I;
}
RealTensorValue C_K;
RealTensorValue B_K;
Real J_K;

C_K = _K.transpose()*_K;
B_K = C_K.inverse();
J_K = _K.det();

RealTensorValue U;
RealTensorValue F;
RealTensorValue Finv;
RealTensorValue FinvT;
RealTensorValue segnaposto;
RealTensorValue C;
RealTensorValue S;
Real J;
Real trCB_K;
Real c1;
Real d1;

U(0,0) = _grad_disp_x[qp](0); U(0,1) = _grad_disp_x[qp](1); U(0,2) = _grad_disp_x[qp](2);
U(1,0) = _grad_disp_y[qp](0); U(1,1) = _grad_disp_y[qp](1); U(1,2) = _grad_disp_y[qp](2);
U(2,0) = _grad_disp_z[qp](0); U(2,1) = _grad_disp_z[qp](1); U(2,2) = _grad_disp_z[qp](2);

F = _I+U;

Finv = F.inverse();

FinvT = Finv.transpose();

J = F.det();

C = F.transpose()*F;

trCB_K = (C*B_K).tr();

segnaposto = (Finv*H*Finv).transpose();


c1 = _mu/2.0;
d1 = _lambda/2.0;

// DA CALCOLARE PER BENE, I CALCOLI DOVREBBERO ESSERE SEMPLICI

	return _mu*J_K*H*B_K + _lambda*J*(log(J)-log(J_K)+1)*(FinvT.contract(H))*FinvT - (_lambda*J*(log(J)-log(J_K)) - _mu*J_K)*segnaposto;
	
}



