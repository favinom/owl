//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElasticMaterial.h"

class LinearElasticityMaterial : public ElasticMaterial
{
public:
  LinearElasticityMaterial(const InputParameters & parameters);

  static InputParameters validParams();
  
  virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp);

protected:
  virtual void computeQpProperties() override;
	Real const _mu;
	Real const _lambda;
  };
