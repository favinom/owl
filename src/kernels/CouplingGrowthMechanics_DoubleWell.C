//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CouplingGrowthMechanics_DoubleWell.h"

registerMooseObject("OwlApp", CouplingGrowthMechanics_DoubleWell);

InputParameters
CouplingGrowthMechanics_DoubleWell::validParams()
{
  InputParameters params = Kernel::validParams();
  	params.addRequiredCoupledVar("J_k", "determinante di K, parametro d'ordine per questo modello");
        params.addRequiredCoupledVar("mu_k", "parte energetica potenziale chimico");
        params.addRequiredCoupledVar("disp_x", "Displacement along x");
  	params.addRequiredCoupledVar("disp_y", "Displacement along y");
  	params.addCoupledVar("disp_z","Displacement along z");
  return params;
}

CouplingGrowthMechanics_DoubleWell::CouplingGrowthMechanics_DoubleWell(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_id_Jk(coupled("J_k")),
_id_muk(coupled("mu_k")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_Jk_old(coupledValueOld("J_k")),
//_grad_Jk(coupledGradient("J_k")),
//_muk(coupledValue("mu_k"))
_F(getMaterialProperty<RealTensorValue>("deformationgradient")),
//_F_old(getMaterialPropertyOld<RealTensorValue>("deformationgradient")),
_P(getMaterialProperty<RealTensorValue>("firstpiolakirchhoff")),
_Psi(getMaterialProperty<Real>("referenceenergy")),
_elastoplasticMaterial(getMaterialProperty<ElastoPlasticMaterial *>("elastoplasticmaterial"))
{}

Real CouplingGrowthMechanics_DoubleWell::computeQpResidual(){
Real min1 = 1;
Real max = 1.01;
Real min2 = 1.05;
Real r = 1*(_Jk_old[_qp]-min1)*(_Jk_old[_qp]-max)*(_Jk_old[_qp]-min2);

  return 1/_Jk_old[_qp]*(_Psi[_qp]-_P[_qp].contract(_F[_qp])/3.0 )*_test[_i][_qp] + r*_test[_i][_qp];
}

Real CouplingGrowthMechanics_DoubleWell::computeQpJacobian(){
  return 0.0;
}


Real CouplingGrowthMechanics_DoubleWell::computeQpOffDiagJacobian(unsigned int jvar){
RealTensorValue H;
RealTensorValue Jac;
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

Jac = _elastoplasticMaterial[_qp][0].evaluateJac(H,_qp);

  return (2.0/3.0*_P[_qp].contract(H)-Jac.contract(_F[_qp])/3.0)*_test[_i][_qp]/_Jk_old[_qp];
}


