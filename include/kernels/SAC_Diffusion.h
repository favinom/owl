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
class SAC_Diffusion : public Kernel
{
public:
  static InputParameters validParams(); //funzione che legge parametri del Kernel

  SAC_Diffusion(const InputParameters & parameters); // costruttore della classe

protected:
	Real _alpha;
	Real _beta;
	RealTensorValue _A;
	Real _min1;
	Real _min2;
	Real _max;

  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
};
