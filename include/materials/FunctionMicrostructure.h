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
#include "Function.h"

class FunctionMicrostructure : public Material
{
public:
  FunctionMicrostructure(const InputParameters & parameters);

  static InputParameters validParams();
  
 // virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp);
  
  int const _dim;
  const Function & _set_K00;
  const Function & _set_K01;
  const Function & _set_K02;
  const Function & _set_K10;
  const Function & _set_K11;
  const Function & _set_K12;
  const Function & _set_K20;
  const Function & _set_K21;
  const Function & _set_K22;  
  
  
protected:
   virtual void computeQpProperties() override;
   
  	RealTensorValue _I;
  	
  /// member variable to hold the computed diffusivity coefficient
  	MaterialProperty<RealTensorValue> & _K;
  	
  	MaterialProperty<RealTensorValue> & _B_K;
  	

	

};
