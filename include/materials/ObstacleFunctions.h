//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "Function.h"

class ObstacleFunctions : public Material
{
public:
  ObstacleFunctions(const InputParameters & parameters);

  static InputParameters validParams();
  
 // virtual RealTensorValue evaluateJac(RealTensorValue const & H, int const & qp);
  
  int const _dim;
  const Function & _g;
  const Function & _gx;
  const Function & _gy;
  const Function & _gz;
  const Function & _gxx;
  const Function & _gxy;
  const Function & _gxz;
  const Function & _gyy;
  const Function & _gyz;
  const Function & _gzz;  
  
  
protected:
   virtual void computeQpProperties() override;
   
   // devo definire dei metodi che restituiscano degli scalari
   
   // virtual Real computeActiveSet(Real const x, Real const y, Real const z, Real const t);
   
   // virtual Real computeDerivativeObstacle(int const i, Real const x, Real const y, Real const z, Real const t);
   
   // virtual RealTensorValue computeHessianObstacle(Real const x, Real const y, Real const z, Real const t);
   
	
	  // valutazione della 
    MaterialProperty<Real> & _geval;
  
    MaterialProperty<RealVectorValue> & _gradgeval;
    
    MaterialProperty<Real> & _normgradgeval;
  	
    MaterialProperty<RealTensorValue> & _Hgeval;
  	
  	
    const VariableValue & _disp_x;
    const VariableValue & _disp_y;
    const VariableValue & _disp_z;

	

};
