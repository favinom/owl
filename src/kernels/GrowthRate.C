//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrowthRate.h"

registerMooseObject("OwlApp", GrowthRate);

InputParameters
GrowthRate::validParams()
{
  InputParameters params = Kernel::validParams();
        params.addParam<Real>("Phi_natural", 0.5, "First Lamé constant");
        params.addParam<Real>("alpha", 1.0, "First Lamé constant");
        params.addParam<Real>("beta", 1.0, "First Lamé constant");
    	params.addRequiredCoupledVar("J_k", "determinante di K, parametro d'ordine per questo modello");
    	params.addRequiredCoupledVar("mu_k", "parte energetica potenziale chimico");
        params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  return params;
}

GrowthRate::GrowthRate(const InputParameters & parameters) : 
Kernel(parameters),
_phi_nat(getParam<Real>("Phi_natural")),
_alpha(getParam<Real>("alpha")),
_beta(getParam<Real>("beta")),
_dim(_mesh.dimension()),
_id_Jk(coupled("J_k")),
_id_muk(coupled("mu_k")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_Jk(coupledValue("J_k")),
//_grad_Jk(coupledGradient("J_k")),
//_muk(coupledValue("mu_k"))
_J(getMaterialProperty<Real>("deformationdeterminant")),
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
//_F_old(getMaterialPropertyOld<RealTensorValue>("deformationgradient")),
_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff")),
_elastoplasticmaterial(getMaterialProperty<ElastoPlasticMaterial *>("elastoplasticmaterial"))
{}





Real GrowthRate::computeQpResidual(){
  return  -_Jk[_qp]*_phi_nat*(_alpha*(_J[_qp]-_Jk[_qp]*_phi_nat) - _beta)*_test[_i][_qp];
}



Real GrowthRate::computeQpJacobian(){
  return (-_phi[_j][_qp]*_phi_nat*(_alpha*(_J[_qp]-_Jk[_qp]*_phi_nat)-_beta)+_Jk[_qp]*_phi_nat*_phi_nat*_alpha*_phi[_j][_qp])*_test[_i][_qp]  ;
}



Real GrowthRate::computeQpOffDiagJacobian(unsigned int jvar){
RealTensorValue H;
//RealTensorValue Jac;
int component_H;

if(jvar == _id_Jk){
	return 0.0;

} else if(jvar == _id_muk){
	return 0.0;

} else if(jvar==_id_x){
	component_H = 0;
}
else if(jvar==_id_y){
	component_H = 1;
}
else if(jvar==_id_z){
	component_H = 2;
} else{
	return 0.0;
}

H(component_H,0) = _grad_phi[_j][_qp](0); H(component_H,1) = _grad_phi[_j][_qp](1); H(component_H,2) = _grad_phi[_j][_qp](2);

//Jac = _elastoplasticMaterial[_qp][0].evaluateJac(H,_qp);

return -_Jk[_qp]*_phi_nat*_alpha*_J[_qp]*((_F[_qp].inverse()).transpose()).contract(H)*_test[_i][_qp];
}













