//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WeakObstacleContactLagMult.h"

registerMooseObject("MooseApp", WeakObstacleContactLagMult);

InputParameters
WeakObstacleContactLagMult::validParams()
{
  InputParameters params = IntegratedBC::validParams();
   	params.addRequiredParam<int>("component", "Component");
   	params.addRequiredParam<bool>("activate_Nanson", "flag to activate Nanson");
   	params.addRequiredParam<bool>("scoperta_del_secolo", "scoperta del secolo");
    	params.addRequiredCoupledVar("disp_x", "displacement along x");
      	params.addRequiredCoupledVar("disp_y", "displacement along y");
        params.addCoupledVar("disp_z", "displacement along z");
      	params.addRequiredCoupledVar("lambda", "Lagrange multiplier");
  return params;
}




WeakObstacleContactLagMult::WeakObstacleContactLagMult(const InputParameters & parameters) : 
IntegratedBC(parameters),
	_dim(_mesh.dimension()),
	_eps(1e-10),
	_activateNanson(getParam<bool>("activate_Nanson")),
	_HFbool(getParam<bool>("scoperta_del_secolo")),
	_component(getParam<int>("component")),
	
	_grad_disp_x(coupledGradient("disp_x")),
	_grad_disp_y(coupledGradient("disp_y")),
	_grad_disp_z(_dim>2 ? coupledGradient("disp_z") : _grad_zero),
	_lambda(coupledValue("lambda")),
    
        _id_x(coupled("disp_x")),
        _id_y(coupled("disp_y")),
	_id_z(_dim>2 ? coupled("disp_z") : -999999),
	_id_lambda(coupled("lambda")), 
	
	_g(getMaterialProperty<Real>("obstacleshape")),
	_gradg(getMaterialProperty<RealVectorValue>("obstaclegradient")),
	_normgradg(getMaterialProperty<Real>("obstaclenormgradient")),
	_Hg(getMaterialProperty<RealTensorValue>("obstaclehessian"))
{
}




Real WeakObstacleContactLagMult::computeQpResidual(){
  
if(_g[_qp] < -_eps || _lambda[_qp] > _eps){
	RealTensorValue F;
	RealVectorValue temp;
	Real da;
	for  (int j=0; j<3; ++j){
		F(0,j) = _grad_disp_x[_qp](j)+(j==0);
		F(1,j) = _grad_disp_y[_qp](j)+(j==1);
		F(2,j) = _grad_disp_z[_qp](j)+(j==2);
	}
	temp = ((F.inverse()).transpose())*_normals[_qp];
	da = 1;
	if(_activateNanson){da = F.det()*std::sqrt(temp.contract(temp));}
	
	return -_lambda[_qp]*_gradg[_qp](_component)*da*_test[_i][_qp];	
	
	/*if(_component ==_id_x){
		return -_lambda[_qp]*_gradg[_qp](0)*da*_test[_i][_qp];	
		//return -_lambda[_qp]*_gradg[_qp](0)/_normgradg[_qp]*da*_test[_i][_qp]; //forma completa con la norma del gradiente
		//return -_lambda[_qp]*_gradg[_qp](0)*_test[_i][_qp];
	}
	else if(_component ==_id_y){
		return -_lambda[_qp]*_gradg[_qp](1)*da*_test[_i][_qp];	
		//return -_lambda[_qp]*_gradg[_qp](1)/_normgradg[_qp]*da*_test[_i][_qp];  //forma completa con la norma del gradiente
		//return -_lambda[_qp]*_gradg[_qp](1)*_test[_i][_qp];
	}
	else if(_component ==_id_z){
		return -_lambda[_qp]*_gradg[_qp](2)*da*_test[_i][_qp];	
		//return -_lambda[_qp]*_gradg[_qp](2)/_normgradg[_qp]*da*_test[_i][_qp]; //forma completa con la norma del gradiente
		//return -_lambda[_qp]*_gradg[_qp](2)*_test[_i][_qp];
	}*/
} 

return 0.0;

// _normals da usare per formula di Nanson
  //return (_grad_u[_qp] * _normals[_qp]) * _test[_i][_qp];
}





Real WeakObstacleContactLagMult::computeQpJacobian(){

// STANDARD: SONO IN CONTATTO SE INEQUALITY = 0 E MOLTIPLICATORE POSITIVO
// STANDARD: HO VIOLATO LA CONDIZIONE DI CONTATTO SE INEQUALITY < 0 
// STANDARD: SONO NEL SET AMMISSIBILE (LIBERO DI MUOVERMI NELLO SPAZIO) SE INEQUALITY >= 0 E MOLTIPLICATORE NULLO

if(_g[_qp] < -_eps || _lambda[_qp] > _eps){
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
	
	da = 1;
	daprime = 0;
	H(_component,0) = _grad_phi[_j][_qp](0); H(_component,1) = _grad_phi[_j][_qp](1); H(_component,2) = _grad_phi[_j][_qp](2);
	if(_activateNanson){
		da = F.det()*std::sqrt(temp.contract(temp));
		daprime = ( da*(Finv.transpose()).contract(H) - F.det()/(std::sqrt(temp*temp)) * ((H*Finv)*temp)*temp  );
	}
	

	
	RealTensorValue Hmat=_Hg[_qp];
	
	if(_HFbool){	Hmat = _Hg[_qp]*F;}
	
	Real Hess_comp = Hmat(_component, _component);
	
	//std::cout << _Hg[_qp] << " " << Hmat << std::endl;
	//std::cout << -_lambda[_qp]*Hess_comp*da*_phi[_j][_qp]*_test[_i][_qp]  << std::endl;
	
	//Real daprime =  ( da*(Finv.transpose()).contract(H) - F.det()/(std::sqrt(temp*temp)) * ((H*Finv)*temp)*temp  );
	
	return -_lambda[_qp]*Hess_comp*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)*_test[_i][_qp]*daprime;
	
	
	/*if(_component ==_id_x){

		
		
		return 
		
		
		
		
		-_lambda[_qp]/_normgradg[_qp]*(_Hg[_qp].row(_component)*(F.transpose()).row(_component) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(_component))/(_normgradg[_qp]*_normgradg[_qp]) )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
		
		//return -_lambda[_qp]/_normgradg[_qp]*(_Hg[_qp].row(_component)*(F.transpose()).row(_component) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(_component))/(_normgradg[_qp]*_normgradg[_qp]) )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
		
		//return -_lambda[_qp]*(_Hg[_qp](0,0)*F(0,0) + _Hg[_qp](0,1)*F(1,0) + _Hg[_qp](0,2)*F(2,0)  )*_phi[_j][_qp]*_test[_i][_qp];
		return -_lambda[_qp]*_Hg[_qp](0,0)*_phi[_j][_qp]*_test[_i][_qp];
	}
	else if(_component ==_id_y){
		
		//return -_lambda[_qp]/_normgradg[_qp]*(_Hg[_qp].row(_component)*(F.transpose()).row(_component) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(_component))/(_normgradg[_qp]*_normgradg[_qp]) )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
		
		//return -_lambda[_qp]*(_Hg[_qp](1,0)*F(0,1) + _Hg[_qp](1,1)*F(1,1) + _Hg[_qp](1,2)*F(2,1) )*_phi[_j][_qp]*_test[_i][_qp];
		return -_lambda[_qp]*_Hg[_qp](1,1)*_phi[_j][_qp]*_test[_i][_qp];
	}
	else if(_component ==_id_z){
		
		//return -_lambda[_qp]/_normgradg[_qp]*(_Hg[_qp].row(_component)*(F.transpose()).row(_component) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(_component))/(_normgradg[_qp]*_normgradg[_qp]) )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
		
		//return -_lambda[_qp]*(_Hg[_qp](2,0)*F(0,2) + _Hg[_qp](2,1)*F(1,2) + _Hg[_qp](2,2)*F(2,2) )*_phi[_j][_qp]*_test[_i][_qp];
		return -_lambda[_qp]*_Hg[_qp](2,2)*_phi[_j][_qp]*_test[_i][_qp];
	}*/
} 

return 0.0;
}













Real WeakObstacleContactLagMult::computeQpOffDiagJacobian(unsigned int jvar){
if(_g[_qp] < -_eps || _lambda[_qp] > _eps){
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
	
	//INIZIALIZZAZIONE
	int componentH;
	if(jvar == _id_x){
		componentH = 0;
	} else if(jvar == _id_y){
		componentH = 1;
	}else if(jvar == _id_z){
		componentH = 2;
	}
	H(componentH,0) = _grad_phi[_j][_qp](0); H(componentH,1) = _grad_phi[_j][_qp](1); H(componentH,2) = _grad_phi[_j][_qp](2);
	
	da = 1;
	daprime = 0;
	if(_activateNanson){
		da = F.det()*std::sqrt(temp.contract(temp));
		daprime = ( da*(Finv.transpose()).contract(H) - F.det()/(std::sqrt(temp*temp)) * ((H*Finv)*temp)*temp  );
	}
	
	
	if(jvar == _id_lambda){
		//std::cout<<"e invece\n";
		return -_phi[_j][_qp]*_gradg[_qp](_component)*da*_test[_i][_qp];	
	}
	

	RealTensorValue Hmat=_Hg[_qp];
	
	if(_HFbool){	Hmat = _Hg[_qp]*F;}
	
	Real Hess_comp = Hmat( _component,componentH);
	
	
	return -_lambda[_qp]*Hess_comp*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)*_test[_i][_qp]*daprime;
	
	//return -_lambda[_qp]*Hess_comp*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)*_test[_i][_qp]*( da*(Finv.transpose()).contract(H) - F.det()/(std::sqrt(temp*temp)) * ((H*Finv)*temp)*temp  );
	
	
	/*
	if(jvar == _id_x){
		if(_component==_id_x){
		std::cout<<"qui non dovrei mai esserci x\n";
			return 0.0;
		}
		else if(_component==_id_y){
			H(jvar,0) = _grad_phi[_j][_qp](0); H(jvar,1) = _grad_phi[_j][_qp](1); H(jvar,2) = _grad_phi[_j][_qp](2);
		
			//return -_lambda[_qp]/_normgradg[_qp]*(  _Hg[_qp].row(_component)*(F.transpose()).row(jvar) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(jvar))/(_normgradg[_qp]*_normgradg[_qp])  )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
		

		
			return -_lambda[_qp]*_Hg[_qp](1,0)*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_z){
			H(jvar,0) = _grad_phi[_j][_qp](0); H(jvar,1) = _grad_phi[_j][_qp](1); H(jvar,2) = _grad_phi[_j][_qp](2);
		
			//return -_lambda[_qp]/_normgradg[_qp]*(  _Hg[_qp].row(_component)*(F.transpose()).row(jvar) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(jvar))/(_normgradg[_qp]*_normgradg[_qp])  )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
		
		
			return -_lambda[_qp]*_Hg[_qp](2,0)*_test[_i][_qp]*_phi[_j][_qp];
		}
	}
	else if(jvar == _id_y){
		if(_component==_id_x){
			H(jvar,0) = _grad_phi[_j][_qp](0); H(jvar,1) = _grad_phi[_j][_qp](1); H(jvar,2) = _grad_phi[_j][_qp](2);
		
			//return -_lambda[_qp]/_normgradg[_qp]*(  _Hg[_qp].row(_component)*(F.transpose()).row(jvar) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(jvar))/(_normgradg[_qp]*_normgradg[_qp])  )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
		
		
			return -_lambda[_qp]*_Hg[_qp](0,1)*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_y){
			std::cout<<"qui non dovrei mai esserci y\n";
			return 0.0;
		}
		else if(_component==_id_z){
			H(jvar,0) = _grad_phi[_j][_qp](0); H(jvar,1) = _grad_phi[_j][_qp](1); H(jvar,2) = _grad_phi[_j][_qp](2);
		
			//return -_lambda[_qp]/_normgradg[_qp]*(  _Hg[_qp].row(_component)*(F.transpose()).row(jvar) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(jvar))/(_normgradg[_qp]*_normgradg[_qp])  )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
			
			return -_lambda[_qp]*_Hg[_qp](2,1)*_test[_i][_qp]*_phi[_j][_qp];
		}
	}
	else if(jvar == _id_z){
		if(_component==_id_x){
			H(jvar,0) = _grad_phi[_j][_qp](0); H(jvar,1) = _grad_phi[_j][_qp](1); H(jvar,2) = _grad_phi[_j][_qp](2);
		
			//return -_lambda[_qp]/_normgradg[_qp]*(  _Hg[_qp].row(_component)*(F.transpose()).row(jvar) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(jvar))/(_normgradg[_qp]*_normgradg[_qp])  )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
			
			
			
			return -_lambda[_qp]*_Hg[_qp](0,2)*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_y){
			H(jvar,0) = _grad_phi[_j][_qp](0); H(jvar,1) = _grad_phi[_j][_qp](1); H(jvar,2) = _grad_phi[_j][_qp](2);
		
			//return -_lambda[_qp]/_normgradg[_qp]*(  _Hg[_qp].row(_component)*(F.transpose()).row(jvar) -_gradg[_qp](_component)*_gradg[_qp]*(_Hg[_qp]*(F.transpose()).row(jvar))/(_normgradg[_qp]*_normgradg[_qp])  )*da*_phi[_j][_qp]*_test[_i][_qp] - _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*da*_test[_i][_qp] *(Finv.transpose()).contract(H) +  _lambda[_qp]*_gradg[_qp](_component)/_normgradg[_qp]*_test[_i][_qp]*F.det()/(std::sqrt(temp*temp))*temp*((H*Finv)*temp);  
			
			
			
			
			return -_lambda[_qp]*_Hg[_qp](1,2)*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_z){
			std::cout<<"qui non dovrei mai esserci z\n";
			return 0.0;
		}
	}
	else if(jvar == _id_lambda){
		if(_component==_id_x){
			return -_gradg[_qp](0)*_test[_i][_qp]*_phi[_j][_qp];
			//return -_gradg[_qp](0)/_normgradg[_qp]*da*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_y){
			return -_gradg[_qp](1)*_test[_i][_qp]*_phi[_j][_qp];			
			//return -_gradg[_qp](1)/_normgradg[_qp]*da*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_z){
			return -_gradg[_qp](2)*_test[_i][_qp]*_phi[_j][_qp];		
			//return -_gradg[_qp](2)/_normgradg[_qp]*da*_test[_i][_qp]*_phi[_j][_qp];
		}
	}*/
	
}

return 0.0;





// vecchia implementazione
/*if(_g[_qp] < -_eps || _lambda[_qp] > _eps){

	if(jvar == _id_x){
		if(_component==_id_x){
			return 0.0;
		}
		else if(_component==_id_y){
			return -_lambda[_qp]*_Hg[_qp](1,0)*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_z){
			return -_lambda[_qp]*_Hg[_qp](2,0)*_test[_i][_qp]*_phi[_j][_qp];
		}
	}
	else if(jvar == _id_y){
		if(_component==_id_x){
			return -_lambda[_qp]*_Hg[_qp](0,1)*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_y){
			return 0.0;
		}
		else if(_component==_id_z){
			return -_lambda[_qp]*_Hg[_qp](2,1)*_test[_i][_qp]*_phi[_j][_qp];
		}
	}
	else if(jvar == _id_z){
		if(_component==_id_x){
			return -_lambda[_qp]*_Hg[_qp](0,2)*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_y){
			return -_lambda[_qp]*_Hg[_qp](1,2)*_test[_i][_qp]*_phi[_j][_qp];
		}
		else if(_component==_id_z){
			return 0.0;
		}
	}
	else if(jvar == _id_lambda){
		if(_component==_id_x){
			return -_gradg[_qp](0)*_test[_i][_qp]*_phi[_j][_qp];;
		}
		else if(_component==_id_y){
			return -_gradg[_qp](1)*_test[_i][_qp]*_phi[_j][_qp];;
		}
		else if(_component==_id_z){
			return -_gradg[_qp](2)*_test[_i][_qp]*_phi[_j][_qp];;
		}
	}
	
}

return 0.0;*/

}
















