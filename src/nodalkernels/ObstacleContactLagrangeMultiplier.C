//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ObstacleContactLagrangeMultiplier.h"

/*Real g1( Real x, Real y, Real z, Real t){ return (1-0.1*t-y)*(x>0.4)*(x<0.6); };
Real gx1( Real x, Real y, Real z, Real t){ return 0; };
Real gy1( Real x, Real y, Real z, Real t){ return -(x>0.4)*(x<0.6); };
Real gz1( Real x, Real y, Real z, Real t){ return 0; };*/
/*Real g1( Real x, Real y, Real z, Real t){ return 4*(x-0.5)*(x-0.5) -0.1*t*(t<3.95) -0.39*(t>3.95)*(t<5.95)+(0.1*(t-5.9)-0.39)*(t>5.95) +1-y; };
Real gx1( Real x, Real y, Real z, Real t){ return 8*(x-0.5); };
Real gy1( Real x, Real y, Real z, Real t){ return -1; };
Real gz1( Real x, Real y, Real z, Real t){ return 0; };
Real gxx1( Real x, Real y, Real z, Real t){ return 8; };
Real gxy1( Real x, Real y, Real z, Real t){ return 0; };
Real gxz1( Real x, Real y, Real z, Real t){ return 0; };
Real gyy1( Real x, Real y, Real z, Real t){ return 0; };
Real gyz1( Real x, Real y, Real z, Real t){ return 0; };
Real gzz1( Real x, Real y, Real z, Real t){ return 0; };*/
// Test per verificare dipendenza dalla normale e non dal gradiente di g
/*Real g1( Real x, Real y, Real z, Real t){ return 1000*(4*(x-0.5)*(x-0.5) -0.1*t*(t<3.95) -0.39*(t>3.95)*(t<5.95)+(0.1*(t-5.9)-0.39)*(t>5.95) +1-y); };
Real gx1( Real x, Real y, Real z, Real t){ return 1000*8*(x-0.5); };
Real gy1( Real x, Real y, Real z, Real t){ return -1000; };
Real gz1( Real x, Real y, Real z, Real t){ return 0; };
Real gxx1( Real x, Real y, Real z, Real t){ return 8*1000; };
Real gxy1( Real x, Real y, Real z, Real t){ return 0; };
Real gxz1( Real x, Real y, Real z, Real t){ return 0; };
Real gyy1( Real x, Real y, Real z, Real t){ return 0; };
Real gyz1( Real x, Real y, Real z, Real t){ return 0; };
Real gzz1( Real x, Real y, Real z, Real t){ return 0; };*/

/*Real g1( Real x, Real y, Real z, Real t){ return 1-y-0.06*t*(t<2.95)-0.18*(t>2.95)*(t<34.95)+(0.06*(t-35)-0.18)*(t>34.95); };
Real gx1( Real x, Real y, Real z, Real t){ return 0; };
Real gy1( Real x, Real y, Real z, Real t){ return -1; };
Real gz1( Real x, Real y, Real z, Real t){ return 0; };
Real gxx1( Real x, Real y, Real z, Real t){ return 0; };
Real gxy1( Real x, Real y, Real z, Real t){ return 0; };
Real gxz1( Real x, Real y, Real z, Real t){ return 0; };
Real gyy1( Real x, Real y, Real z, Real t){ return 0; };
Real gyz1( Real x, Real y, Real z, Real t){ return 0; };
Real gzz1( Real x, Real y, Real z, Real t){ return 0; };*/

Real g1( Real x, Real y, Real z, Real t){ return 4*(x-0.5)*(x-0.5) -0.05*t*(t<6.95)-0.35*(t>6.95)*(t<24.95)+(0.05*(t-25)-0.35)*(t>24.95)+1-y; };
Real gx1( Real x, Real y, Real z, Real t){ return 8*(x-0.5); };
Real gy1( Real x, Real y, Real z, Real t){ return -1; };
Real gz1( Real x, Real y, Real z, Real t){ return 0; };
Real gxx1( Real x, Real y, Real z, Real t){ return 8; };
Real gxy1( Real x, Real y, Real z, Real t){ return 0; };
Real gxz1( Real x, Real y, Real z, Real t){ return 0; };
Real gyy1( Real x, Real y, Real z, Real t){ return 0; };
Real gyz1( Real x, Real y, Real z, Real t){ return 0; };
Real gzz1( Real x, Real y, Real z, Real t){ return 0; };


Real eps1 = 1e-10;


registerMooseObject("MooseApp", ObstacleContactLagrangeMultiplier);

InputParameters
ObstacleContactLagrangeMultiplier::validParams()
{
  InputParameters params = NodalKernel::validParams();
   	params.addRequiredParam<int>("component", "Component");
    	params.addRequiredCoupledVar("disp_x", "displacement along x");
      	params.addRequiredCoupledVar("disp_y", "displacement along y");
        params.addCoupledVar("disp_z", "displacement along z");
      	params.addRequiredCoupledVar("lambda", "Lagrange multiplier");
  return params;
}

ObstacleContactLagrangeMultiplier::ObstacleContactLagrangeMultiplier(const InputParameters & parameters) : 
	NodalKernel(parameters),
	_dim(_mesh.dimension()),
	_component(getParam<int>("component")),

    	_disp_x(coupledValue("disp_x")),
    	_disp_y(coupledValue("disp_y")),
	_disp_z(_dim>2 ? coupledValue("disp_z") : _zero),
	_lambda(coupledValue("lambda")),
    
        _id_x(coupled("disp_x")),
        _id_y(coupled("disp_y")),
	_id_z(_dim>2 ? coupled("disp_z") : -999999),
	_id_lambda(coupled("lambda"))
        
{
  /*if (_var.number() == _v_var)
    mooseError("Coupled variable 'v' needs to be different from 'variable' with "
               "LowerBoundNodalKernel");

  const auto & bnd_names = getParam<std::vector<BoundaryName>>("exclude_boundaries");
  for (const auto & bnd_name : bnd_names)
    _bnd_ids.insert(_mesh.getBoundaryID(bnd_name));*/
}










Real
ObstacleContactLagrangeMultiplier::computeQpResidual()
{
Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO


// mettere tolleranza per moltiplicatore
if(g1(chi_x,chi_y,chi_z,_t) < -eps1 || _lambda[_qp] > eps1){
	if(_component ==_id_x){
		return -_lambda[_qp]*gx1(chi_x,chi_y,chi_z,_t);
	}
	else if(_component ==_id_y){
		return -_lambda[_qp]*gy1(chi_x,chi_y,chi_z,_t);
	}
	else if(_component ==_id_z){
		return -_lambda[_qp]*gz1(chi_x,chi_y,chi_z,_t);
	}
} 

return 0.0;

}












Real ObstacleContactLagrangeMultiplier::computeQpJacobian(){

Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// Hessiana fatta ancora a mano
RealTensorValue Hg;
Hg(0,0) = gxx1(chi_x,chi_y,chi_z,_t);
Hg(0,1) = gxy1(chi_x,chi_y,chi_z,_t);
Hg(0,2) = gxz1(chi_x,chi_y,chi_z,_t);
Hg(1,1) = gyy1(chi_x,chi_y,chi_z,_t);
Hg(1,2) = gyz1(chi_x,chi_y,chi_z,_t);
Hg(2,2) = gzz1(chi_x,chi_y,chi_z,_t);

Hg(1,0) = Hg(0,1);
Hg(2,0) = Hg(0,2);
Hg(2,1) = Hg(1,2);

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

// mettere tolleranza per moltiplicatore
if(g1(chi_x,chi_y,chi_z,_t) < -eps1 || _lambda[_qp] > eps1){
	if(_component ==_id_x){
		return -_lambda[_qp]*Hg(0,0);
	}
	else if(_component ==_id_y){
		return -_lambda[_qp]*Hg(1,1);
	}
	else if(_component ==_id_z){
		return -_lambda[_qp]*Hg(2,2);
	}
} 

// SE LA FUNZIONE CHE DESCRIVE L'OSTACOLO E' LINEARE NELLE VARIABILI SPAZIALI IL CONTRIBUTO AL JACOBIANO, CHE DIPENDE DALLA MATRICE HESSIANA, E' NULLO.




return 0.0;

}













Real
ObstacleContactLagrangeMultiplier::computeQpOffDiagJacobian(unsigned int jvar){

	Real chi_x = (*_current_node)(0) + _disp_x[_qp];
	Real chi_y = (*_current_node)(1) + _disp_y[_qp];
	Real chi_z = (*_current_node)(2) + _disp_z[_qp];
	
if(g1(chi_x,chi_y,chi_z,_t) < -eps1 || _lambda[_qp] > eps1){


	
	// Hessiana fatta ancora a mano
	RealTensorValue Hg;
	Hg(0,0) = gxx1(chi_x,chi_y,chi_z,_t);
	Hg(0,1) = gxy1(chi_x,chi_y,chi_z,_t);
	Hg(0,2) = gxz1(chi_x,chi_y,chi_z,_t);
	Hg(1,1) = gyy1(chi_x,chi_y,chi_z,_t);
	Hg(1,2) = gyz1(chi_x,chi_y,chi_z,_t);
	Hg(2,2) = gzz1(chi_x,chi_y,chi_z,_t);

	Hg(1,0) = Hg(0,1);
	Hg(2,0) = Hg(0,2);
	Hg(2,1) = Hg(1,2);



	if(jvar == _id_x){
		if(_component==_id_x){
			return 0.0;
		}
		else if(_component==_id_y){
			return -_lambda[_qp]*Hg(1,0);
		}
		else if(_component==_id_z){
			return -_lambda[_qp]*Hg(2,0);
		}
	}
	else if(jvar == _id_y){
		if(_component==_id_x){
			return -_lambda[_qp]*Hg(0,1);
		}
		else if(_component==_id_y){
			return 0.0;
		}
		else if(_component==_id_z){
			return -_lambda[_qp]*Hg(2,1);
		}
	}
	else if(jvar == _id_z){
		if(_component==_id_x){
			return -_lambda[_qp]*Hg(0,2);
		}
		else if(_component==_id_y){
			return -_lambda[_qp]*Hg(1,2);
		}
		else if(_component==_id_z){
			return 0.0;
		}
	}
	else if(jvar == _id_lambda){
		if(_component==_id_x){
			return -gx1(chi_x,chi_y,chi_z,_t);
		}
		else if(_component==_id_y){
			return -gy1(chi_x,chi_y,chi_z,_t);
		}
		else if(_component==_id_z){
			return -gz1(chi_x,chi_y,chi_z,_t);
		}
	}
	
}

return 0.0;




/*if(jvar == _id_lambda){

Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO


if(g1(chi_x,chi_y,chi_z,_t) < -eps1 || _lambda[_qp] > eps1){
	if(_component==_id_x){
		return -gx1(chi_x,chi_y,chi_z,_t);
	}
	else if(_component==_id_y){
		return -gy1(chi_x,chi_y,chi_z,_t);
	}
	else if(_component==_id_z){
		return -gz1(chi_x,chi_y,chi_z,_t);
	}
}
}

return 0.0;*/

}
