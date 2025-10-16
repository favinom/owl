//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HolmesMow.h"

registerMooseObject("OwlApp", HolmesMow);

InputParameters
HolmesMow::validParams()
{
  InputParameters params = ElasticMaterial::validParams();
  params.addParam<Real>("mu", 1.0, "First Lamé constant");
  params.addParam<Real>("lambda", 1.0, "Second Lamé constant");
  return params;
}

HolmesMow::HolmesMow(const InputParameters & parameters) : 
ElasticMaterial(parameters),
_mu(getParam<Real>("mu")),
_lambda(getParam<Real>("lambda"))
{
}





void HolmesMow::computeQpProperties()
{
if(_Kmat){
	_K = _Kmat[_qp];
} 
else{
	_K = _I;
}
// Quantità dipendenti da K
RealTensorValue C_K;
RealTensorValue B_K;
Real J_K;

C_K = _K.transpose()*_K;
B_K = C_K.inverse();
J_K = _K.det();

// Proprietà di F
_elasticMaterial[_qp] = this;
RealTensorValue U;
RealTensorValue F;
RealTensorValue FinvT;
RealTensorValue C;
RealTensorValue CB_K;
Real trC;
Real J;

U(0,0) = _grad_disp_x[_qp](0); U(0,1) = _grad_disp_x[_qp](1); U(0,2) = _grad_disp_x[_qp](2);
U(1,0) = _grad_disp_y[_qp](0); U(1,1) = _grad_disp_y[_qp](1); U(1,2) = _grad_disp_y[_qp](2);
U(2,0) = _grad_disp_z[_qp](0); U(2,1) = _grad_disp_z[_qp](1); U(2,2) = _grad_disp_z[_qp](2);

F = _I+U;

J = F.det();

FinvT = F.inverse().transpose();

C = F.transpose()*F;

trC = C.tr();

CB_K = C*B_K;

// Coefficienti del materiale
Real alpha0;
Real alpha1;
Real alpha2;
Real alpha3;

alpha0 = (2*_mu+_lambda)/4;
alpha1 = (2*_mu-_lambda)/(2*_mu+_lambda);
alpha2 = _lambda/(2*_mu+_lambda);
alpha3 = 1;


// Invarianti
Real I_1e;
Real I_2e;
Real I_3e;
Real Psi;

I_1e = CB_K.tr();
I_2e = 0.5*(I_1e*I_1e - (CB_K*CB_K).tr());
I_3e = (J*J)/(J_K*J_K);

Psi = exp(alpha1*(I_1e -3) + alpha2*(I_2e-3) - alpha3*log(I_3e));


	_stress[_qp] = 2*alpha0*J_K*Psi*( (alpha1+alpha2*I_1e)*F*B_K - alpha2*F*B_K*CB_K - alpha3*FinvT);

// std::cout << (alpha1+alpha2*I_1e)*F*B_K - alpha2*F*B_K*CB_K;

}







RealTensorValue HolmesMow::evaluateJac(RealTensorValue const & H, int const & qp){
if(_Kmat){
	_K = _Kmat[qp];
} 
else{
	_K = _I;
}

// Quantità dipendenti da K
RealTensorValue C_K;
RealTensorValue B_K;
Real J_K;

C_K = _K.transpose()*_K;
B_K = C_K.inverse();
J_K = _K.det();

// Proprietà di F
RealTensorValue U;
RealTensorValue F;
RealTensorValue FinvT;
RealTensorValue C;
RealTensorValue Cinv;
RealTensorValue CB_K;
RealTensorValue FB_K;
//RealTensorValue W_FTH;
RealTensorValue D_FTH;
Real trC;
Real J;

U(0,0) = _grad_disp_x[qp](0); U(0,1) = _grad_disp_x[qp](1); U(0,2) = _grad_disp_x[qp](2);
U(1,0) = _grad_disp_y[qp](0); U(1,1) = _grad_disp_y[qp](1); U(1,2) = _grad_disp_y[qp](2);
U(2,0) = _grad_disp_z[qp](0); U(2,1) = _grad_disp_z[qp](1); U(2,2) = _grad_disp_z[qp](2);

F = _I+U;

J = F.det();

FinvT = (F.inverse()).transpose();

C = F.transpose()*F;

Cinv = C.inverse();

trC = C.tr();

CB_K = C*B_K;

FB_K = F*B_K;

D_FTH = ((F.transpose()*H) + H.transpose()*F)/2.0;

// Coefficienti del materiale
Real alpha0;
Real alpha1;
Real alpha2;
Real alpha3;

alpha0 = (2*_mu+_lambda)/4;
alpha1 = (2*_mu-_lambda)/(2*_mu+_lambda);
alpha2 = _lambda/(2*_mu+_lambda);
alpha3 = 1;


// Invarianti
Real I_1e;
Real I_2e;
Real I_3e;
Real ePsi;
Real e_Psi;
Real Psi;

I_1e = CB_K.tr();
I_2e = 0.5*(I_1e*I_1e - (CB_K*CB_K).tr());
I_3e = (J*J)/(J_K*J_K);
Psi = alpha1*(I_1e -3) + alpha2*(I_2e-3) - alpha3*log(I_3e);
ePsi = exp(Psi);
e_Psi = exp(-Psi);

// Secondo Piola-Kirchhoff
RealTensorValue S;
Real DerivataPsi;
S = 2*alpha0*J_K*ePsi*((alpha1+alpha2*I_1e)*B_K - alpha2*B_K*CB_K - alpha3*Cinv);
//segnaposto = 2*(alpha1*FB_K +alpha2*I_1e*FB_K - alpha2*FB_K*CB_K - alpha3*FinvT );
DerivataPsi = (2*alpha1+2*alpha2*I_1e)*FB_K.contract(H)-2*alpha2*(CB_K*D_FTH*B_K).tr()-2*alpha3*FinvT.contract(H);

//std::cout << DerivataPsi;

	//return (H + DerivataPsi*F)*S + 2*alpha0*J_K*ePsi*( 2*alpha2*(D_FTH.tr())*FB_K - 2*alpha2*FB_K*D_FTH*B_K + alpha3*(H*F.inverse() + FinvT*H.transpose())*FinvT);
	
	
	return (H + DerivataPsi*F)*S + 2*alpha0*J_K*ePsi*( 2*alpha2*(FB_K.contract(H))*FB_K - 2*alpha2*FB_K*D_FTH*B_K + alpha3*(H*F.inverse() + FinvT*H.transpose())*FinvT);


	//return (H + segnaposto.contract(H)*F)*S + 2*alpha0*J_K*ePsi*( 2*alpha2*FB_K.contract(H)*FB_K - 2*alpha2*FB_K*D_FTH*B_K + alpha3*H*Cinv + alpha3*FinvT*H.transpose()*F.transpose());

	//return (H + DerivataPsi*F)*S + 2*alpha0*J_K*ePsi*( 2*alpha2*(D_FTH.tr())*FB_K - 2*alpha2*FB_K*D_FTH*B_K + alpha3*H*Cinv + alpha3*FinvT*H.transpose()*F.transpose());




//std::cout << 8*alpha0*J_K*Psi*(0.5*FB_K.contract(H)*FB_K - 0.5*alpha2*FB_K*D_FTH*B_K + alpha3*FinvT*D_FTH*Cinv);



	//return (2*((alpha1+alpha2*I_1e)*(FB_K).contract(H) - alpha2*(FB_K*CB_K).contract(H) - alpha3*FinvT.contract(H))*F+H)*S + 8*alpha0*J_K*Psi*(0.5*FB_K.contract(H)*FB_K - 0.5*alpha2*FB_K*D_FTH*B_K + alpha3*FinvT*D_FTH*Cinv); 
	
	//return _I;

}
