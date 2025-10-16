//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once
// keyword per evitare loop negli include

#include "Kernel.h"
#include "ElastoPlasticMaterial.h"
// include della classe madre

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class CouplingGrowthMechanics : public Kernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  CouplingGrowthMechanics(const InputParameters & parameters); // costruttore della classe

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  	int const _dim;
	int const _id_Jk;
    	int const _id_muk;
    	int const _id_x;
    	int const _id_y;
    	int const _id_z;
  	VariableValue const & _Jk_old;
//  	VariableGradient const & _grad_Jk;
//  	VariableValue const & _muk;  

	MaterialProperty<RealTensorValue> const & _F;
    	
//    	MaterialProperty<RealTensorValue> const & _F_old;
  
  	MaterialProperty<RealTensorValue> const & _P;
  	
  	MaterialProperty<Real> const & _Psi;
  	
  	MaterialProperty<ElastoPlasticMaterial *> const & _elastoplasticMaterial;
  
};
