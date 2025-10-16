//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemicalPotentialReaction.h"

registerMooseObject("OwlApp", ChemicalPotentialReaction);

InputParameters
ChemicalPotentialReaction::validParams()
{
  InputParameters params = Kernel::validParams();
  	params.addParam<Real>("d", 1.0, "Second Lamé constant");
  	params.addRequiredCoupledVar("J_k", "determinante di K, parametro d'ordine per questo modello");
        params.addRequiredCoupledVar("mu_k", "parte energetica potenziale chimico");
  return params;
}

ChemicalPotentialReaction::ChemicalPotentialReaction(const InputParameters & parameters) : 
Kernel(parameters),
_d(getParam<Real>("d")),
_id_Jk(coupled("J_k")),
_id_muk(coupled("mu_k")),
_Jk(coupledValue("J_k")),
_grad_Jk(coupledGradient("J_k")),
_muk(coupledValue("mu_k"))
{}

Real ChemicalPotentialReaction::computeQpResidual(){
  return _muk[_qp]*_test[_i][_qp] + _d *_grad_Jk[_qp] * _grad_test[_i][_qp];
}

Real ChemicalPotentialReaction::computeQpJacobian(){
  return _phi[_j][_qp]*_test[_i][_qp];
}


Real ChemicalPotentialReaction::computeQpOffDiagJacobian(unsigned int jvar){
if(jvar == _id_Jk){

	return +_d *_grad_phi[_j][_qp] * _grad_test[_i][_qp];

} else if(jvar == _id_muk){

	return 0.0;

}

return 0.0;

}


