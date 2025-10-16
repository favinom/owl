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

#include "TimeKernel.h"
// include della classe madre

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class DiffusiveTimeDerivative : public TimeKernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  DiffusiveTimeDerivative(const InputParameters & parameters); // costruttore della classe

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
  	Real const _m;
	Real const _kv;
  	int const _id_Jk;
    	int const _id_muk;
  	VariableValue const & _Jk;
  	VariableValue const & _Jk_old;
  	VariableGradient const & _grad_Jk;
  	VariableGradient const & _grad_Jk_old;
  	VariableValue const & _muk;
  	
};
