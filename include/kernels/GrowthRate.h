//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"
#include "ElastoPlasticMaterial.h"

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class GrowthRate : public Kernel
{
public:
  static InputParameters validParams();

  GrowthRate(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  	
  	Real const _phi_nat;
	Real const _alpha;
	Real const _beta;
  	
  	int const _dim;
	int const _id_Jk;
	int const _id_muk;
    	int const _id_x;
    	int const _id_y;
    	int const _id_z;

  	VariableValue const & _Jk;
  	
	MaterialProperty<Real> const & _J;

	MaterialProperty<RealTensorValue> const & _F;
    	
  	MaterialProperty<RealTensorValue> const & _P;
  	
  	MaterialProperty<ElastoPlasticMaterial *> const & _elastoplasticmaterial;
};
