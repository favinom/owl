//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DoubleWellCH.h"

registerMooseObject("OwlApp", DoubleWellCH);

InputParameters
DoubleWellCH::validParams()
{
  InputParameters params = Kernel::validParams();
      	params.addParam<Real>("coeff", 2.0, "Coefficient to multiply by the body force term");
  	params.addParam<Real>("min1", 1.0, "Coefficient to multiply by the body force term");
    	params.addParam<Real>("max", 1.5, "Coefficient to multiply by the body force term");
    	params.addParam<Real>("min2", 2.0, "Coefficient to multiply by the body force term");
    	params.addRequiredCoupledVar("J_k", "determinante di K, parametro d'ordine per questo modello");
    	params.addRequiredCoupledVar("mu_k", "potenziale chimico");
  return params;
}

DoubleWellCH::DoubleWellCH(const InputParameters & parameters) : 
Kernel(parameters),
    	_coeff(getParam<Real>("coeff")),
    	_min1(getParam<Real>("min1")),
        _max(getParam<Real>("max")),
        _min2(getParam<Real>("min2")),
        _id_Jk(coupled("J_k")),
	_id_muk(coupled("mu_k")),
	_Jk(coupledValue("J_k")),
	_muk(coupledValue("mu_k"))
         {
         }


Real DoubleWellCH::computeQpResidual(){

Real r = (_Jk[_qp]-_min1)*(_Jk[_qp]-_max)*(_Jk[_qp]-_min2);

  return _coeff*r*_test[_i][_qp];
}



Real DoubleWellCH::computeQpJacobian(){

	return 0.0;
}



Real DoubleWellCH::computeQpOffDiagJacobian(unsigned int jvar){
if(jvar == _id_Jk){
	Real r1 = (_Jk[_qp]-_min1)*(_Jk[_qp]-_max);
	Real r2 = (_Jk[_qp]-_max)*(_Jk[_qp]-_min2);
	Real r3 = (_Jk[_qp]-_min1)*(_Jk[_qp]-_min2);
	Real r = r1 + r2 + r3;
  	return _coeff*r*_phi[_j][_qp]*_test[_i][_qp];

} else if(jvar == _id_muk){

	return 0.0;

}

return 0.0;

}












