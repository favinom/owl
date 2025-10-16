//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FunctionMicrostructure.h"

registerMooseObject("OwlApp", FunctionMicrostructure);

InputParameters
FunctionMicrostructure::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<FunctionName>("set_K00", "1", "description");
  params.addParam<FunctionName>("set_K01", "0", "description");
  params.addParam<FunctionName>("set_K02", "0", "description");
  params.addParam<FunctionName>("set_K10", "0", "description");
  params.addParam<FunctionName>("set_K11", "1", "description");
  params.addParam<FunctionName>("set_K12", "0", "description");
  params.addParam<FunctionName>("set_K20", "0", "description");
  params.addParam<FunctionName>("set_K21", "0", "description");
  params.addParam<FunctionName>("set_K22", "1", "description");
  return params;
}

FunctionMicrostructure::FunctionMicrostructure(const InputParameters & parameters) : 
Material(parameters),
_dim(_mesh.dimension()),
_set_K00(getFunction("set_K00")),
_set_K01(getFunction("set_K01")),
_set_K02(getFunction("set_K02")),
_set_K10(getFunction("set_K10")),
_set_K11(getFunction("set_K11")),
_set_K12(getFunction("set_K12")),
_set_K20(getFunction("set_K20")),
_set_K21(getFunction("set_K21")),
_set_K22(getFunction("set_K22")),
_K(declareProperty<RealTensorValue>("microstructure")),
_B_K(declareProperty<RealTensorValue>("microstructureBK"))
{
for(int i=0;i<3;++i){
	for (int j=0;j<3;++j){
	_I(i,j) = (i==j);
	}
}
//if(_dim == 2){_I(2,2) = 0.0;}

}


void FunctionMicrostructure::computeQpProperties()
{
_K[_qp](0,0) = _set_K00.value(_t,_q_point[_qp]);
_K[_qp](0,1) = _set_K01.value(_t,_q_point[_qp]);
_K[_qp](0,2) = _set_K02.value(_t,_q_point[_qp]);
_K[_qp](1,0) = _set_K10.value(_t,_q_point[_qp]);
_K[_qp](1,1) = _set_K11.value(_t,_q_point[_qp]);
_K[_qp](1,2) = _set_K12.value(_t,_q_point[_qp]);
_K[_qp](2,0) = _set_K20.value(_t,_q_point[_qp]);
_K[_qp](2,1) = _set_K21.value(_t,_q_point[_qp]);
_K[_qp](2,2) = _set_K22.value(_t,_q_point[_qp]);

_B_K[_qp] = _K[_qp].transpose()*_K[_qp];
}
