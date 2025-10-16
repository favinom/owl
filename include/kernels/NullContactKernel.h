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
#include "ObstacleFunctions.h"
/**
 *
 */
class NullContactKernel : public Kernel
{
public:
  static InputParameters validParams();

  NullContactKernel(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// filler value to put on the on-diagonal Jacobian
  const Real _jacobian_fill;
  Real const _eps;
  
  MaterialProperty<Real> const & _g;
};
