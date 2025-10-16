//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NeoHookean.h"

registerMooseObject("OwlApp", NeoHookean);

InputParameters
NeoHookean::validParams()
{
  InputParameters params = ElasticMaterial::validParams();
  params.addParam<Real>("mu", 1.0, "First Lamé constant");
  params.addParam<Real>("lambda", 1.0, "Second Lamé constant");
  return params;
}

NeoHookean::NeoHookean(const InputParameters & parameters) : 
ElasticMaterial(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda"))
{
}



void NeoHookean::computeQpProperties()
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

	_stress[_qp] = 2*c1*J_K*F*B_K + (2*d1/J_K*J*(J-J_K) - 2*c1*J_K)*FinvT;
	//_mu*J_K*(-(C*B_K).tr()/3.0*FinvT + F*B_K) + (_lambda+2.0/3.0*_mu)*J*(J-J_K)/J_K*FinvT;
	//_mu*(-trC/3.0*FinvT+F) + (_lambda+2.0/3.0*_mu)*J*(J-1)*FinvT;
}




RealTensorValue NeoHookean::evaluateJac(RealTensorValue const & H, int const & qp){
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


	return 2*c1*J_K*H*B_K + 2*d1/J_K*(2*J-J_K)*J*(FinvT.contract(H))*FinvT - (2*d1/J_K*J*(J-J_K) - 2*c1*J_K)*segnaposto;

	//return 2J_K*_mu*(H*B_K - 2/3*((F*B_K).contract(H))*FinvT + trCB_K/3*segnaposto) + (_lambda+2/3*_mu)/J_K*( J*(2*J-J_K)*(FinvT.contract(H))*FinvT - J*(J-J_K)*segnaposto);



	//return J_K*_mu*((-2*(F*B_K*H).tr()*FinvT + (C*B_K).tr()*segnaposto)/3.0 + H*B_K) +(_lambda+2.0/3.0*_mu)/J_K*(J*(2*J-J_K)*FinvT.contract(H)*FinvT - J*(J-J_K)*segnaposto);
	
	
}



