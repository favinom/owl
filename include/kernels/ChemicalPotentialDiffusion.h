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
// include della classe madre

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class ChemicalPotentialDiffusion : public Kernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  ChemicalPotentialDiffusion(const InputParameters & parameters); // costruttore della classe

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  	Real const _m;
	int const _id_Jk;
    	int const _id_muk;
  	VariableValue const & _Jk;
  	VariableValue const & _muk;
  	VariableGradient const & _grad_muk;
  
  
  
};
