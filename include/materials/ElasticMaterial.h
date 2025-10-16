//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

class ElasticMaterial : public Material
{
public:
  ElasticMaterial(const InputParameters & parameters);

  static InputParameters validParams();
  
  virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp) = 0;

protected:
  // virtual void computeQpProperties() override;
    	int const _dim;
    	VariableGradient const & _grad_disp_x;
    	VariableGradient const & _grad_disp_y;
    	VariableGradient const & _grad_disp_z;
  	RealTensorValue _I;
  	RealTensorValue _K;
  /// member variable to hold the computed diffusivity coefficient
  	MaterialProperty<RealTensorValue> & _stress;
  	
  	OptionalMaterialProperty<RealTensorValue> const & _Kmat;
  	
  	MaterialProperty<ElasticMaterial *> & _elasticMaterial;

	

};
