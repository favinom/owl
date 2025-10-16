//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PoroLinearElasticity.h"

registerMooseObject("OwlApp", PoroLinearElasticity);

InputParameters
PoroLinearElasticity::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredParam<int>("component", "Component");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  params.addRequiredCoupledVar("p", "pressure");
  //params.addCoupledVar("density",
  //                     "The name of the temperature variable used in the "
  //                     "ComputeThermalExpansionEigenstrain.  (Not required for "
  //                     "simulations without temperature coupling.)");
  return params;
}

PoroLinearElasticity::PoroLinearElasticity(const InputParameters & parameters) : 
Kernel(parameters),
_dim(_mesh.dimension()),
_component(getParam<int>("component")),
_p(coupledValue("p")),
_id_x(coupled("disp_x")),
_id_y(coupled("disp_y")),
_id_z(_dim>2 ? coupled("disp_z") : -999999),
_id_p(coupled("p")),
_stress(getMaterialProperty<RealTensorValue>("stress")),
_elasticmaterial(getMaterialProperty<ElasticMaterial *>("elasticmaterial"))
 {
 //std::cout << _dim << std::endl;
 }



Real PoroLinearElasticity::computeQpResidual(){
RealTensorValue I;
RealTensorValue V;

for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
		V(i,j) = 0;
		I(i,j) = (i==j);
	}
}

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

  return (_stress[_qp]-_p[_qp]*I).contract(V);
}





Real
PoroLinearElasticity::computeQpJacobian(){
RealTensorValue H;
RealTensorValue Jac;
RealTensorValue V;

H(_component,0) = _grad_phi[_j][_qp](0); H(_component,1) = _grad_phi[_j][_qp](1); H(_component,2) = _grad_phi[_j][_qp](2);

Jac = _elasticmaterial[_qp][0].evaluateJac(H,_qp);

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

  return Jac.contract(V);
}




Real
PoroLinearElasticity::computeQpOffDiagJacobian(unsigned int jvar) 
{
RealTensorValue H;
RealTensorValue Jac;
RealTensorValue V;
int component_H;

V(_component,0) = _grad_test[_i][_qp](0); V(_component,1) = _grad_test[_i][_qp](1); V(_component,2) = _grad_test[_i][_qp](2);

if(jvar == _id_p){
	return -_phi[_j][_qp]*V(_component,_component);
}
else if(jvar==_id_x){
	component_H = 0;
}
else if(jvar==_id_y){
	component_H = 1;
}
else if(jvar==_id_z){
	component_H = 2;
}
else{
	return 0.0;
}

H(component_H,0) = _grad_phi[_j][_qp](0); H(component_H,1) = _grad_phi[_j][_qp](1); H(component_H,2) = _grad_phi[_j][_qp](2);

Jac = _elasticmaterial[_qp][0].evaluateJac(H,_qp);

  return Jac.contract(V);
}


