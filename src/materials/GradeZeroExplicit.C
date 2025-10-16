//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GradeZeroExplicit.h"

registerMooseObject("OwlApp", GradeZeroExplicit);

InputParameters
GradeZeroExplicit::validParams()
{
  InputParameters params = Microstructure::validParams();

  return params;
}

GradeZeroExplicit::GradeZeroExplicit(const InputParameters & parameters) : 
Microstructure(parameters),
_detK(declareProperty<Real>("determinanteK")),
_B_K_old(getMaterialPropertyOld<RealTensorValue>("microstructureBK"))
{
}


void GradeZeroExplicit::initQpStatefulProperties(){
_B_K[_qp] = _I;
_detK[_qp] = 1;
}

void GradeZeroExplicit::computeQpProperties()
{
RealTensorValue F;
RealTensorValue P;


// mettere if

F = _F_old[_qp];
P = _P_old[_qp];

_sigma = (P*F.transpose())/(F.det());
_SIGMA = (F.transpose())*P;

RealTensorValue devsigma;
RealTensorValue devSIGMA;

devsigma = _sigma - (_sigma.tr()/3.0)*_I;
devSIGMA = _SIGMA - (_SIGMA.tr()/3.0)*_I;

Real normdevsigma;
normdevsigma = std::sqrt(devsigma.contract(devsigma));

Real gamma_p;
if(normdevsigma > 0.000001){
	gamma_p = _lambda_p/normdevsigma*(normdevsigma-std::sqrt(2.0/3.0)*_sigma_y + std::fabs(normdevsigma-std::sqrt(2.0/3.0)*_sigma_y ));
} 
else{
	gamma_p = 0;
}

RealTensorValue A;
A = _I + gamma_p*_dt*devSIGMA;
_B_K[_qp] = _B_K_old[_qp]*A.inverse();


Real toll;
Real J_K;
toll = 0.01;
J_K = 1/(std::sqrt(_B_K[_qp].det()));
_detK[_qp] = J_K;

//if(abs(J_K-1) > toll){std::cout << J_K << std::endl;}




}
