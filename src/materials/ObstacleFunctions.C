//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ObstacleFunctions.h"

registerMooseObject("OwlApp", ObstacleFunctions);

InputParameters
ObstacleFunctions::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredParam<FunctionName>("set_g", "description");
  params.addRequiredParam<FunctionName>("set_gx", "description");
  params.addRequiredParam<FunctionName>("set_gy", "description");
  params.addParam<FunctionName>("set_gz", "0", "description");
  params.addRequiredParam<FunctionName>("set_gxx", "description");
  params.addRequiredParam<FunctionName>("set_gxy", "description");
  params.addParam<FunctionName>("set_gxz", "0", "description");
  params.addRequiredParam<FunctionName>("set_gyy", "description");
  params.addParam<FunctionName>("set_gyz", "0", "description");
  params.addParam<FunctionName>("set_gzz", "0", "description");
  params.addRequiredCoupledVar("disp_x", "Displacement along x");
  params.addRequiredCoupledVar("disp_y", "Displacement along y");
  params.addCoupledVar("disp_z","Displacement along z");
  return params;
}

ObstacleFunctions::ObstacleFunctions(const InputParameters & parameters) : 
Material(parameters),
_dim(_mesh.dimension()),
_g(getFunction("set_g")),
_gx(getFunction("set_gx")),
_gy(getFunction("set_gy")),
_gz(getFunction("set_gz")),
_gxx(getFunction("set_gxx")),
_gxy(getFunction("set_gxy")),
_gxz(getFunction("set_gxz")),
_gyy(getFunction("set_gyy")),
_gyz(getFunction("set_gyz")),
_gzz(getFunction("set_gzz")),
_geval(declareProperty<Real>("obstacleshape")),
_gradgeval(declareProperty<RealVectorValue>("obstaclegradient")),
_normgradgeval(declareProperty<Real>("obstaclenormgradient")),
_Hgeval(declareProperty<RealTensorValue>("obstaclehessian")),
_disp_x(coupledValue("disp_x")),
_disp_y(coupledValue("disp_y")),
_disp_z(_dim>2 ? coupledValue("disp_z") : _zero)
{
}
//if(_dim == 2){_I(2,2) = 0.0;}




void ObstacleFunctions::computeQpProperties()
{
RealVectorValue chi;
chi(0) = (_q_point[_qp])(0) + _disp_x[_qp];
chi(1) = (_q_point[_qp])(1) + _disp_y[_qp];
chi(2) = (_q_point[_qp])(2) + _disp_z[_qp];

_geval[_qp] = _g.value(_t,chi);
_gradgeval[_qp](0) = _gx.value(_t,chi);
_gradgeval[_qp](1) = _gy.value(_t,chi);
_gradgeval[_qp](2) = _gz.value(_t,chi);
_normgradgeval[_qp] = std::sqrt(_gradgeval[_qp]*_gradgeval[_qp]);
_Hgeval[_qp](0,0) = _gxx.value(_t,chi);
_Hgeval[_qp](0,1) = _gxy.value(_t,chi);
_Hgeval[_qp](0,2) = _gxz.value(_t,chi);
_Hgeval[_qp](1,0) = _gxy.value(_t,chi);
_Hgeval[_qp](1,1) = _gyy.value(_t,chi);
_Hgeval[_qp](1,2) = _gyz.value(_t,chi);
_Hgeval[_qp](2,0) = _gxz.value(_t,chi);
_Hgeval[_qp](2,1) = _gyz.value(_t,chi);
_Hgeval[_qp](2,2) = _gzz.value(_t,chi);
}








