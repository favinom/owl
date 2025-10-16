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
//#include "ObstacleFunctions.h"

/**
 * Class used to enforce a lower bound on a coupled variable
 */
class IdentityNoContact : public NodalKernel
{
public:
  static InputParameters validParams();

  IdentityNoContact(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

    	int const _dim;
    	Real const _eps;
  /// The value of the coupled variable
  	const VariableValue & _disp_x;
    	const VariableValue & _disp_y;
	const VariableValue & _disp_z;	
	
	
	
    	//MaterialProperty<Real> const & _g;

  /// Boundaries on which we should not execute this object
  //std::set<BoundaryID> _bnd_ids;

  /// The lower bound on the coupled variable
  //const Real _lower_bound;
};
