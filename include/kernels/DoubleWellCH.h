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
class DoubleWellCH : public Kernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  DoubleWellCH(const InputParameters & parameters); // costruttore della classe

protected:
	Real const _coeff;
	Real const _min1;
	Real const _max;
	Real const _min2;
	
	int const _id_Jk;
    	int const _id_muk;
  	VariableValue const & _Jk;
  	VariableValue const & _muk;

  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
  
    virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
};
