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
class PoroElastoPlasticityMassBalance : public Kernel
{
public:
  static InputParameters validParams();

  PoroElastoPlasticityMassBalance(const InputParameters & parameters);

protected:
	int const _dim;
	Real const _k;
    	VariableGradient const & _grad_p;
    	int const _id_x;
    	int const _id_y;
    	int const _id_z;
    	int const _id_p;
    	
    	
    	MaterialProperty<RealTensorValue> const & _F;
    	MaterialProperty<Real> const & _J;
    	MaterialProperty<Real> const & _J_old;
    	MaterialProperty<RealTensorValue> const & _P;
    
 // MaterialProperty<RealTensorValue> const & _stress;
     
 // MaterialProperty<ElasticMaterial *> const & _elasticmaterial;
  
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
};


