//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeakEnforceObstacleConstraint.h"

registerMooseObject("MooseApp", WeakEnforceObstacleConstraint);

InputParameters
WeakEnforceObstacleConstraint::validParams()
{
  InputParameters params = IntegratedBC::validParams();
    	params.addRequiredCoupledVar("disp_x", "displacement along x");
      	params.addRequiredCoupledVar("disp_y", "displacement along y");
        params.addCoupledVar("disp_z", "displacement along z");
  return params;
}




WeakEnforceObstacleConstraint::WeakEnforceObstacleConstraint(const InputParameters & parameters) : 
IntegratedBC(parameters),
	_dim(_mesh.dimension()),
	_eps(1e-10),
	_activateNanson(0),
	
	_grad_disp_x(coupledGradient("disp_x")),
	_grad_disp_y(coupledGradient("disp_y")),
	_grad_disp_z(_dim>2 ? coupledGradient("disp_z") : _grad_zero),
    
        _id_x(coupled("disp_x")),
        _id_y(coupled("disp_y")),
	_id_z(_dim>2 ? coupled("disp_z") : -999999),
	
	_g(getMaterialProperty<Real>("obstacleshape")),
	_gradg(getMaterialProperty<RealVectorValue>("obstaclegradient")),
	_Hg(getMaterialProperty<RealTensorValue>("obstaclehessian"))
{
}




Real WeakEnforceObstacleConstraint::computeQpResidual(){

  
if(_g[_qp] < -_eps || _u[_qp] > _eps){
	RealTensorValue F;
	RealVectorValue temp;
	Real da;
	for  (int j=0; j<3; ++j){
		F(0,j) = _grad_disp_x[_qp](j)+(j==0); 
		
		F(1,j) = _grad_disp_y[_qp](j)+(j==1);
		F(2,j) = _grad_disp_z[_qp](j)+(j==2);
	}
	temp = ((F.inverse()).transpose())*_normals[_qp];
	//da = (_activateNanson == 0)*F.det()*std::sqrt(temp*temp)+(_activateNanson == 1);
	da = F.det()*std::sqrt(temp*temp);



	return _g[_qp]*_test[_i][_qp]*da;
}

return 0.0;

// _normals da usare per formula di Nanson
  //return (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp];
}





Real WeakEnforceObstacleConstraint::computeQpJacobian(){

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

// I TERMINI DIAGONALI DEVONO ESSERE DATI UNA SOLA VOLTA: VEDERE IL NODAL KERNEL DI RIFERIMENTO
//if(_g[_qp] >= -_eps && _u[_qp] <= _eps){
//	return 1.0;
//} 

return 0.0;
}





Real WeakEnforceObstacleConstraint::computeQpOffDiagJacobian(unsigned int jvar){
	//std::cout<<"offdiagonal\n" << _g[_qp] << "   " << _u[_qp] << std::endl;
if(_g[_qp] < -_eps || _u[_qp] > _eps){
	//std::cout<<"sto calcolando off diagonal\n";
	RealTensorValue F;
	RealTensorValue Finv;
	RealVectorValue temp;
	RealTensorValue H;
	Real da;
	Real daprime;
	for  (int j=0; j<3; ++j){
		F(0,j) = _grad_disp_x[_qp](j)+(j==0); 
		F(1,j) = _grad_disp_y[_qp](j)+(j==1);
		F(2,j) = _grad_disp_z[_qp](j)+(j==2);
		H(0,j) = 0; H(1,j) = 0; H(2,j) = 0;
	}
	Finv = F.inverse();
	temp = (Finv.transpose())*_normals[_qp];
	
	int componentH;
	if(jvar == _id_x){
		componentH = 0;
	} else if(jvar == _id_y){
		componentH = 1;
	}else if(jvar == _id_z){
		componentH = 2;
	} else{
		std::cout<<"qui non dovrei mai esserci id\n";
		return 0.0;
	}
	H(componentH,0) = _grad_phi[_j][_qp](0); H(componentH,1) = _grad_phi[_j][_qp](1); H(componentH,2) = _grad_phi[_j][_qp](2);
	
	
	//da = (_activateNanson == 0)*F.det()*std::sqrt(temp*temp)+(_activateNanson == 1);
	da = F.det()*std::sqrt(temp*temp);
	daprime = (da*(Finv.transpose()).contract(H) - F.det()/(std::sqrt(temp*temp)) * ((H*Finv)*temp)*temp  );
	
	//std::cout<<"offdiagonal\n" << _gradg[_qp](componentH)*da << "   " << daprime << std::endl;
	
	RealVectorValue gradgmat = _gradg[_qp];
	bool _gradgF = false;
	if(_gradgF){  gradgmat = F.transpose()*_gradg[_qp];}
	Real grad_comp = _gradg[_qp](componentH);
	
	
	return grad_comp*_phi[_j][_qp]*_test[_i][_qp]*da + _g[_qp]*_test[_i][_qp]*daprime;
	
	
	

	/*if(jvar == _id_x){
		H(0,0) = _grad_phi[_j][_qp](0); H(0,1) = _grad_phi[_j][_qp](1); H(0,2) = _grad_phi[_j][_qp](2);
		return _gradg[_qp](0)*_phi[_j][_qp]*_test[_i][_qp];
		
		//return _gradg[_qp]*(F.transpose()).row(0)*_phi[_j][_qp]*_test[_i][_qp];
		
		//return _gradg[_qp]*(F.transpose()).row(0)*da*_phi[_j][_qp]*_test[_i][_qp] + _g[_qp]*_test[_i][_qp]*(da*(Finv.transpose()).contract(H) - F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp));
	}
	else if(jvar == _id_y){
		H(1,0) = _grad_phi[_j][_qp](0); H(1,1) = _grad_phi[_j][_qp](1); H(1,2) = _grad_phi[_j][_qp](2);
		return _gradg[_qp](1)*_phi[_j][_qp]*_test[_i][_qp];		
		
		//return _gradg[_qp]*(F.transpose()).row(1)*_phi[_j][_qp]*_test[_i][_qp];
		
		//return _gradg[_qp]*(F.transpose()).row(1)*da*_phi[_j][_qp]*_test[_i][_qp] + _g[_qp]*_test[_i][_qp]*(da*(Finv.transpose()).contract(H) - F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp));

	}
	else if(jvar == _id_z){
		H(2,0) = _grad_phi[_j][_qp](0); H(2,1) = _grad_phi[_j][_qp](1); H(2,2) = _grad_phi[_j][_qp](2);
		return _gradg[_qp](2)*_phi[_j][_qp]*_test[_i][_qp];		
		
		//return _gradg[_qp]*(F.transpose()).row(2)*_phi[_j][_qp]*_test[_i][_qp];

		//return _gradg[_qp]*(F.transpose()).row(2)*da*_phi[_j][_qp]*_test[_i][_qp] + _g[_qp]*_test[_i][_qp]*(da*(Finv.transpose()).contract(H) - F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp));
	}*/
}

return 0.0;

}
















