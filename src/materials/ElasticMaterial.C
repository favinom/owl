//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElasticMaterial.h"

//registerMooseObject("OwlApp", ElasticMaterial);

InputParameters
ElasticMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  return params;
}

ElasticMaterial::ElasticMaterial(const InputParameters & parameters) : 
Material(parameters),
_dim(_mesh.dimension()),
_grad_disp_x(coupledGradient("disp_x")),
_grad_disp_y(coupledGradient("disp_y")),
_grad_disp_z(_dim>2 ? coupledGradient("disp_z") : _grad_zero),
_stress(declareProperty<RealTensorValue>("stress")),
_Kmat(getOptionalMaterialProperty<RealTensorValue>("microstructure")),
_elasticMaterial(declareProperty<ElasticMaterial *>("elasticmaterial"))
{
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
	_I(i,j) = (i==j);
	}
}
//if(_dim == 2){_I(2,2) = 0.0;}

}
