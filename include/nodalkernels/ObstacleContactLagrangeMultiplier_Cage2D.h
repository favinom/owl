//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NodalKernel.h"

/**
 * Class used to enforce a lower bound on a coupled variable
 */
class ObstacleContactLagrangeMultiplier_Cage2D : public NodalKernel
{
public:
  static InputParameters validParams();

  ObstacleContactLagrangeMultiplier_Cage2D(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

private:
    	int const _dim;
    	int const _component;
  /// The value of the coupled variable
  	const VariableValue & _disp_x;
    	const VariableValue & _disp_y;
	const VariableValue & _disp_z;
	const VariableValue & _lambda;
	
	// identifiers
	int const _id_x;
    	int const _id_y;
    	int const _id_z;
    	int const _id_lambda;


  /// Boundaries on which we should not execute this object
  //std::set<BoundaryID> _bnd_ids;

  /// The lower bound on the coupled variable
  //const Real _lower_bound;
};
