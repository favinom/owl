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
class PoroElastoPlasticity : public Kernel
{
public:
  static InputParameters validParams();

  PoroElastoPlasticity(const InputParameters & parameters);

protected:
	int const _dim;
	int const _component;
	VariableValue const & _p;
    	int const _id_x;
    	int const _id_y;
    	int const _id_z;
	int const _id_p;
    
    
        MaterialProperty<RealTensorValue> const & _F;
    	MaterialProperty<Real> const & _J;
    	MaterialProperty<RealTensorValue> const & _P;
    	
  //mMaterialProperty<RealTensorValue> const & _stress;
     
  MaterialProperty<ElastoPlasticMaterial *> const & _elastoplasticmaterial;
  
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
};


