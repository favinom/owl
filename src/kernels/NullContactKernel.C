//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NullContactKernel.h"

registerMooseObject("OwlApp", NullContactKernel);

InputParameters
NullContactKernel::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Kernel that sets a zero residual.");
  params.addParam<Real>(
      "jacobian_fill",
      1,
      "On diagonal Jacobian fill term to retain an invertible matrix for the preconditioner");
  return params;
}

NullContactKernel::NullContactKernel(const InputParameters & parameters)
  : Kernel(parameters), 
  _jacobian_fill(getParam<Real>("jacobian_fill")),
  _eps(1e-10),
  _g(getMaterialProperty<Real>("obstacleshape"))
{
}

Real NullContactKernel::computeQpResidual(){
  return 0.0;
}




Real NullContactKernel::computeQpJacobian(){  
if(_g[_qp] >= -_eps && _u[_qp] <= _eps){
  	return _jacobian_fill;
} 

return 0.0;
}
  
  



Real NullContactKernel::computeQpOffDiagJacobian(unsigned int jvar){

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

return 0.0;

}


