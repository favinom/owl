//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChemicalPotentialDiffusion.h"

registerMooseObject("OwlApp", ChemicalPotentialDiffusion);

InputParameters
ChemicalPotentialDiffusion::validParams()
{
  InputParameters params = Kernel::validParams();
    	params.addParam<Real>("m", 1.0, "First Lamé constant");
  	params.addRequiredCoupledVar("J_k", "determinante di K, parametro d'ordine per questo modello");
        params.addRequiredCoupledVar("mu_k", "parte energetica potenziale chimico");
  return params;
}

ChemicalPotentialDiffusion::ChemicalPotentialDiffusion(const InputParameters & parameters) : 
Kernel(parameters),
_m(getParam<Real>("m")),
_id_Jk(coupled("J_k")),
_id_muk(coupled("mu_k")),
_Jk(coupledValue("J_k")),
_muk(coupledValue("mu_k")),
_grad_muk(coupledGradient("mu_k"))
{}

Real ChemicalPotentialDiffusion::computeQpResidual(){
  return -_m*_grad_muk[_qp] * _grad_test[_i][_qp];
}

Real ChemicalPotentialDiffusion::computeQpJacobian(){
  return 0.0;
}


Real ChemicalPotentialDiffusion::computeQpOffDiagJacobian(unsigned int jvar){
if(jvar == _id_Jk){

	return 0.0;

} else if(jvar == _id_muk){

	return -_m*_grad_phi[_j][_qp] * _grad_test[_i][_qp];

}

return 0.0;

}


