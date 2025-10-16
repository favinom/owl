//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GradeOneGrowth.h"

registerMooseObject("OwlApp", GradeOneGrowth);

InputParameters
GradeOneGrowth::validParams()
{
  InputParameters params = Microstructure::validParams();
	params.addParam<bool>("is_explicit" , 1 , "aaaaaaaaaaa");
  	params.addRequiredCoupledVar("J_k", "determinante di K, parametro d'ordine per questo modello");
  return params;
}

GradeOneGrowth::GradeOneGrowth(const InputParameters & parameters) : 
Microstructure(parameters),
_is_explicit(getParam<bool>("is_explicit")),
_Jk(coupledValue("J_k")),
_Jk_old(coupledValueOld("J_k"))
{
}


void GradeOneGrowth::initQpStatefulProperties(){
_B_K[_qp] = std::pow(_Jk[_qp],-2.0/3.0)*_I;
}

void GradeOneGrowth::computeQpProperties()
{


if(_is_explicit){
	_B_K[_qp] = std::pow(_Jk_old[_qp],-2.0/3.0)*_I;
} else{
	_B_K[_qp] = std::pow(_Jk[_qp],-2.0/3.0)*_I;
}

}
