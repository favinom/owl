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
class ElastoPlasticity : public Kernel
{
public:
  static InputParameters validParams();

  ElastoPlasticity(const InputParameters & parameters);

protected:
	int const _dim;
	int const _component;
   	int const _id_x;
    	int const _id_y;
    	int const _id_z;
    
  MaterialProperty<RealTensorValue> const & _P;
     
  MaterialProperty<ElastoPlasticMaterial *> const & _elastoplasticMaterial;
  
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
};


