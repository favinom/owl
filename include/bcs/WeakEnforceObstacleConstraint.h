//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "IntegratedBC.h"
#include "ObstacleFunctions.h"

/**
 * A FluxBC which is consistent with the boundary terms arising from
 * the Diffusion Kernel. The residual contribution is:
 *
 * \f$ F(u) = - \int_{\Gamma} \nabla u * \hat n * \phi d\Gamma \f$
 *
 * This class is essentially identical to the DiffusionFluxBC, but it
 * is not a part of the FluxBC hierarchy. It does not actually impose
 * any boundary condition, instead it computes the residual
 * contribution due to the boundary term arising from integration by
 * parts of the Diffusion Kernel.
 */
class WeakEnforceObstacleConstraint : public IntegratedBC
{
public:
  static InputParameters validParams();

  WeakEnforceObstacleConstraint(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  
    	int const _dim;
    	Real const _eps;
    	int const _activateNanson;
    	
    	VariableGradient const & _grad_disp_x;
    	VariableGradient const & _grad_disp_y;
    	VariableGradient const & _grad_disp_z;
	
	// identifiers
	int const _id_x;
    	int const _id_y;
    	int const _id_z;
    	
    	 MaterialProperty<Real> const & _g;
    	 MaterialProperty<RealVectorValue> const & _gradg;
    	 MaterialProperty<RealTensorValue> const & _Hg;
};
