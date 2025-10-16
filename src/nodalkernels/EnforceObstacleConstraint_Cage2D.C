//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EnforceObstacleConstraint_Cage2D.h"

/*Real g( Real x, Real y, Real z, Real t){ return (1-0.1*t-y)*(x>0.4)*(x<0.6); };
Real gx( Real x, Real y, Real z, Real t){ return 0; };
Real gy( Real x, Real y, Real z, Real t){ return -(x>0.4)*(x<0.6); };
Real gz( Real x, Real y, Real z, Real t){ return 0; };*/
/*Real g( Real x, Real y, Real z, Real t){ return 4*(x-0.5)*(x-0.5) -0.1*t*(t<3.95) -0.39*(t>3.95)*(t<5.95)+(0.1*(t-5.9)-0.39)*(t>5.95) +1-y; };
Real gx( Real x, Real y, Real z, Real t){ return 8*(x-0.5); };
Real gy( Real x, Real y, Real z, Real t){ return -1; };
Real gz( Real x, Real y, Real z, Real t){ return 0; };*/
/*Real g( Real x, Real y, Real z, Real t){ return 1-y-0.06*t*(t<2.95)-0.18*(t>2.95)*(t<34.95)+(0.06*(t-35)-0.18)*(t>34.95);};
Real gx( Real x, Real y, Real z, Real t){ return 0; };
Real gy( Real x, Real y, Real z, Real t){ return -1; };
Real gz( Real x, Real y, Real z, Real t){ return 0; };*/
// Test per verificare dipendenza dalla normale e non dal gradiente di g
/*Real g( Real x, Real y, Real z, Real t){ return 1000*(4*(x-0.5)*(x-0.5) -0.1*t*(t<3.95) -0.39*(t>3.95)*(t<5.95)+(0.1*(t-5.9)-0.39)*(t>5.95) +1-y); };
Real gx( Real x, Real y, Real z, Real t){ return 1000*8*(x-0.5); };
Real gy( Real x, Real y, Real z, Real t){ return -1000; };
Real gz( Real x, Real y, Real z, Real t){ return 0; };*/

Real cage( Real x, Real y, Real z, Real t){ return (1.2-y)*(1.2-x)*(0.2+x);  };
Real cagex( Real x, Real y, Real z, Real t){ return (1.2-y)*(-2*x+1); };
Real cagey( Real x, Real y, Real z, Real t){ return -(1.2-x)*(0.2+x);; };
Real cagez( Real x, Real y, Real z, Real t){ return 0; };

Real eps_cage2D = 1e-10;

registerMooseObject("MooseApp", EnforceObstacleConstraint_Cage2D);

InputParameters
EnforceObstacleConstraint_Cage2D::validParams()
{
  InputParameters params = NodalKernel::validParams();
    	params.addRequiredCoupledVar("disp_x", "displacement along x");
      	params.addRequiredCoupledVar("disp_y", "displacement along y");
        params.addCoupledVar("disp_z", "displacement along z");
  return params;
}

// La variabile locale _u è il nome del moltiplicatore di Lagrange associato al contatto nella classe che stai utilizzando

EnforceObstacleConstraint_Cage2D::EnforceObstacleConstraint_Cage2D(const InputParameters & parameters) : 
	NodalKernel(parameters),
	_dim(_mesh.dimension()),

    	_disp_x(coupledValue("disp_x")),
    	_disp_y(coupledValue("disp_y")),
	_disp_z(_dim>2 ? coupledValue("disp_z") : _zero),
    
        _id_x(coupled("disp_x")),
        _id_y(coupled("disp_y")),
	_id_z(_dim>2 ? coupled("disp_z") : -999999)
        

{
  /*if (_var.number() == _v_var)
    mooseError("Coupled variable 'v' needs to be different from 'variable' with "
               "LowerBoundNodalKernel");

  const auto & bnd_names = getParam<std::vector<BoundaryName>>("exclude_boundaries");
  for (const auto & bnd_name : bnd_names)
    _bnd_ids.insert(_mesh.getBoundaryID(bnd_name));*/
}

Real
EnforceObstacleConstraint_Cage2D::computeQpResidual()
{
Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

// mettere tolleranza per moltiplicatore
if(cage(chi_x,chi_y,chi_z,_t) < -eps_cage2D || _u[_qp] > eps_cage2D){
	return cage(chi_x,chi_y,chi_z,_t);
} else{
	// PROVA PER ESSERE SICURI CHE LAMBDA SIA NON NEGATIVO
	return _u[_qp];
	//return 0.0;
}

  /*for (auto bnd_id : _bnd_ids)
    if (_mesh.isBoundaryNode(_current_node->id(), bnd_id))
      return _u[_qp];

  return std::min(_u[_qp], _v[_qp] - _lower_bound);*/
}





Real EnforceObstacleConstraint_Cage2D::computeQpJacobian()
{
Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

// mettere tolleranza per moltiplicatore
if(cage(chi_x,chi_y,chi_z,_t) >= -eps_cage2D && _u[_qp] <= eps_cage2D){
	return 1.0;
} 

return 0.0;
}

Real
EnforceObstacleConstraint_Cage2D::computeQpOffDiagJacobian(unsigned int jvar){

Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO


if(cage(chi_x,chi_y,chi_z,_t) < -eps_cage2D || _u[_qp] > eps_cage2D){
	if(jvar==_id_x){
		return cagex(chi_x,chi_y,chi_z,_t);
	}
	else if(jvar==_id_y){
		return cagey(chi_x,chi_y,chi_z,_t);
	}
	else if(jvar==_id_z){
		return cagez(chi_x,chi_y,chi_z,_t);
	}
}

return 0.0;

}
