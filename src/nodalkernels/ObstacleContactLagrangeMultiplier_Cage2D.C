//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ObstacleContactLagrangeMultiplier_Cage2D.h"

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

Real cage1( Real x, Real y, Real z, Real t){ return (1.2-y)*(1.2-x)*(0.2+x); };
Real cagex1( Real x, Real y, Real z, Real t){ return (1.2-y)*(-2*x+1); };
Real cagey1( Real x, Real y, Real z, Real t){ return -(1.2-x)*(0.2+x); };
Real cagez1( Real x, Real y, Real z, Real t){ return 0; };
Real cagexx1( Real x, Real y, Real z, Real t){ return -2*(1.2-y); };
Real cagexy1( Real x, Real y, Real z, Real t){ return 2*x-1; };
Real cagexz1( Real x, Real y, Real z, Real t){ return 0; };
Real cageyy1( Real x, Real y, Real z, Real t){ return 0; };
Real cageyz1( Real x, Real y, Real z, Real t){ return 0; };
Real cagezz1( Real x, Real y, Real z, Real t){ return 0; };


Real eps1_cage2D = 1e-10;


registerMooseObject("MooseApp", ObstacleContactLagrangeMultiplier_Cage2D);

InputParameters
ObstacleContactLagrangeMultiplier_Cage2D::validParams()
{
  InputParameters params = NodalKernel::validParams();
   	params.addRequiredParam<int>("component", "Component");
    	params.addRequiredCoupledVar("disp_x", "displacement along x");
      	params.addRequiredCoupledVar("disp_y", "displacement along y");
        params.addCoupledVar("disp_z", "displacement along z");
      	params.addRequiredCoupledVar("lambda", "Lagrange multiplier");
  return params;
}

ObstacleContactLagrangeMultiplier_Cage2D::ObstacleContactLagrangeMultiplier_Cage2D(const InputParameters & parameters) : 
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
ObstacleContactLagrangeMultiplier_Cage2D::computeQpResidual()
{
Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO


// mettere tolleranza per moltiplicatore
if(cage1(chi_x,chi_y,chi_z,_t) < -eps1_cage2D || _lambda[_qp] > eps1_cage2D){
	if(_component ==_id_x){
		return -_lambda[_qp]*cagex1(chi_x,chi_y,chi_z,_t);
	}
	else if(_component ==_id_y){
		return -_lambda[_qp]*cagey1(chi_x,chi_y,chi_z,_t);
	}
	else if(_component ==_id_z){
		return -_lambda[_qp]*cagez1(chi_x,chi_y,chi_z,_t);
	}
} 

return 0.0;

}












Real ObstacleContactLagrangeMultiplier_Cage2D::computeQpJacobian(){

Real chi_x = (*_current_node)(0) + _disp_x[_qp];
Real chi_y = (*_current_node)(1) + _disp_y[_qp];
Real chi_z = (*_current_node)(2) + _disp_z[_qp];

// Hessiana fatta ancora a mano
RealTensorValue Hg;
Hg(0,0) = cagexx1(chi_x,chi_y,chi_z,_t);
Hg(0,1) = cagexy1(chi_x,chi_y,chi_z,_t);
Hg(0,2) = cagexz1(chi_x,chi_y,chi_z,_t);
Hg(1,1) = cageyy1(chi_x,chi_y,chi_z,_t);
Hg(1,2) = cageyz1(chi_x,chi_y,chi_z,_t);
Hg(2,2) = cagezz1(chi_x,chi_y,chi_z,_t);

Hg(1,0) = Hg(0,1);
Hg(2,0) = Hg(0,2);
Hg(2,1) = Hg(1,2);

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

// mettere tolleranza per moltiplicatore
if(cage1(chi_x,chi_y,chi_z,_t) < -eps1_cage2D || _lambda[_qp] > eps1_cage2D){
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
ObstacleContactLagrangeMultiplier_Cage2D::computeQpOffDiagJacobian(unsigned int jvar){

	Real chi_x = (*_current_node)(0) + _disp_x[_qp];
	Real chi_y = (*_current_node)(1) + _disp_y[_qp];
	Real chi_z = (*_current_node)(2) + _disp_z[_qp];
	
if(cage1(chi_x,chi_y,chi_z,_t) < -eps1_cage2D || _lambda[_qp] > eps1_cage2D){


	
	// Hessiana fatta ancora a mano
	RealTensorValue Hg;
	Hg(0,0) = cagexx1(chi_x,chi_y,chi_z,_t);
	Hg(0,1) = cagexy1(chi_x,chi_y,chi_z,_t);
	Hg(0,2) = cagexz1(chi_x,chi_y,chi_z,_t);
	Hg(1,1) = cageyy1(chi_x,chi_y,chi_z,_t);
	Hg(1,2) = cageyz1(chi_x,chi_y,chi_z,_t);
	Hg(2,2) = cagezz1(chi_x,chi_y,chi_z,_t);

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
			return -cagex1(chi_x,chi_y,chi_z,_t);
		}
		else if(_component==_id_y){
			return -cagey1(chi_x,chi_y,chi_z,_t);
		}
		else if(_component==_id_z){
			return -cagez1(chi_x,chi_y,chi_z,_t);
		}
	}
	
}

return 0.0;

}
