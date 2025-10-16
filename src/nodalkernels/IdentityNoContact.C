//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IdentityNoContact.h"

Real g2( Real x, Real y, Real z, Real t){ 
return -(x*x+y*y-x-2*y*(-3+1)-2*3+1) -0.1*t;
//return (4*(x-0.5)*(x-0.5)-0.1*t+1-y); 
};

//x*x+y*y-x-2*y*(-3+1)-2*3+1 -0.1*t 

registerMooseObject("MooseApp", IdentityNoContact);

InputParameters
IdentityNoContact::validParams()
{
  InputParameters params = NodalKernel::validParams();
      	params.addRequiredCoupledVar("disp_x", "displacement along x");
      	params.addRequiredCoupledVar("disp_y", "displacement along y");
        params.addCoupledVar("disp_z", "displacement along z");
  return params;
}

// La variabile locale _u è il nome del moltiplicatore di Lagrange associato al contatto nella classe che stai utilizzando

IdentityNoContact::IdentityNoContact(const InputParameters & parameters) : 
NodalKernel(parameters),
	_dim(_mesh.dimension()),
	_eps(1e-10),
	_disp_x(coupledValue("disp_x")),
    	_disp_y(coupledValue("disp_y")),
	_disp_z(_dim>2 ? coupledValue("disp_z") : _zero)
	//_g(getMaterialProperty<Real>("obstacleshape"))
{
}





Real IdentityNoContact::computeQpResidual(){
Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

Real xx=(*_current_node)(0);
Real yy=(*_current_node)(1);

//if ( 0.4-_eps<xx && xx < 0.4+_eps && yy > 1-_eps ){
	//std::cout << (*_current_node) << "  " <<  _u[_qp] << "   " << g2(chi_x,chi_y,chi_z,_t) << std::endl;
//}


// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

// mettere tolleranza per moltiplicatore
if(g2(chi_x,chi_y,chi_z,_t) < -_eps || _u[_qp] > _eps){
	return 0.0;
} else if (_u[_qp] < 0.0) { 
	// PROVA PER ESSERE SICURI CHE LAMBDA SIA NON NEGATIVO
	//std::cout << (*_current_node) << "  " <<  _u[_qp] << std::endl;
	return _u[_qp];
	//return 0.0;
}

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

return 0.0;
}





Real IdentityNoContact::computeQpJacobian(){
Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

// mettere tolleranza per moltiplicatore
if(g2(chi_x,chi_y,chi_z,_t) >= -_eps && _u[_qp] <= _eps){
	return 1.0;
} 

return 0.0;
}






Real IdentityNoContact::computeQpOffDiagJacobian(unsigned int jvar){

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

return 0.0;

}
