/***************************************************************************
 *            function.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <cstring>

#include "base/stlio.h"
#include "base/exceptions.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/exceptions.h"
#include "function/affine_variable.h"
#include "function/taylor_derivative.h"
#include "function/virtual_machine.h"
#include "input/modelica.h"

#include "function/virtual_machine.template.h"

#include "interpreted_function.h"



namespace Ariadne {


std::ostream& 
Function::operator<<(std::ostream& os, const Variable& var) 
{
  os << " Real";
  if(var.array_flag==true) { os << "[" << var.size << "]"; }
  os << " " << var.name;
  return os;
}

std::ostream& 
Function::operator<<(std::ostream& os, const FunctionVariable& var) 
{
  switch(var.type) {
    case FunctionVariable::OUTPUT: os << "output "; break;
    case FunctionVariable::INPUT: os << "input "; break;
    case FunctionVariable::INTERMEDIATE: os << ""; break;
    case FunctionVariable::CONSTANT: os << "input "; break;
  }
  os << "Real";
  if(var.array_flag==true) { os << "[" << var.size << "]"; }
  os << " " << var.name;
  return os;
}


template<class R>
void
Function::InterpretedFunction<R>::_initialise()
{
  this->_argument_size=0;
  this->_result_size=0;
  this->_intermediates_size=0;
  this->_constants_size=0;

  for(size_type i=0; i!=this->_variables.size(); ++i) {
    const FunctionVariable& var=this->_variables[i];
    
    switch(var.type) {
    case FunctionVariable::OUTPUT: 
      this->_result_size+=var.size; 
      break;
    case FunctionVariable::INPUT: 
      this->_argument_size+=var.size; 
      break;
    case FunctionVariable::INTERMEDIATE: 
      this->_intermediates_size+=var.size; 
      break;
    case FunctionVariable::CONSTANT: 
      _constants_size+=var.size;
      break;
    }
  }
  assert(this->_constants.size()==this->_constants_size);
}


template<class R>
Function::InterpretedFunction<R>::InterpretedFunction()
{
  this->_initialise();
}


template<class R>
Function::InterpretedFunction<R>::InterpretedFunction(const std::string& str)
{
  std::stringstream ss(str);
  this->read(ss);
}


template<class R>
Function::InterpretedFunction<R>::InterpretedFunction(std::istream& is)
{
  this->read(is);
}





template<class R>
LinearAlgebra::Vector<typename Function::InterpretedFunction<R>::A>
Function::InterpretedFunction<R>::evaluate(const LinearAlgebra::Vector<A>& x) const
{
  ARIADNE_CHECK_ARGUMENT_SIZE(*this,x,"Function::image(Vector x)");
  LinearAlgebra::Vector<A> y(this->result_size());
  array<A> t(this->_intermediates_size);
  const array<A>& c=this->_constants;

  A* args[4]={ y.begin(),const_cast<A*>(x.begin()),t.begin(),const_cast<A*>(c.begin()) };
  VirtualMachine().evaluate(this->_operations,args);
  return y;
}




template<class R>
typename Function::TaylorDerivative<typename Function::InterpretedFunction<R>::A>
Function::InterpretedFunction<R>::derivative(const LinearAlgebra::Vector<A>& x, const smoothness_type& s) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
LinearAlgebra::Matrix<typename Function::InterpretedFunction<R>::A>
Function::InterpretedFunction<R>::jacobian(const LinearAlgebra::Vector<A>& x) const
{
  typedef AffineVariable<A> FD;
  size_type m=this->result_size();
  size_type n=this->argument_size();
  FD dx=FD::constant(n,0.0);
  array<FD> y(m,dx);
  array<FD> a(n,dx);
  array<FD> t(this->_intermediates_size);
  array<FD> c(this->_constants_size);
  for(size_type i=0; i!=n; ++i) {
    a[i]=FD::variable(n,x[i],i);
  }
  for(size_type i=0; i!=c.size(); ++i) {
    c[i]=FD::constant(n,this->_constants[i]);
  }
  ARIADNE_LOG(7,"y="<<y<<"  a="<<a<<"  t="<<t<<"  c="<<c);
  FD* args[4]={ y.begin(), a.begin(), t.begin(), c.begin() };
  VirtualMachine machine;
  machine.evaluate(this->_operations,args);
  ARIADNE_LOG(7,"y="<<y<<"  t="<<t);
  LinearAlgebra::Matrix<A> result(m,n);
  for(size_type i=0; i!=m; ++i) {
    for(size_type j=0; j!=n; ++j) {
      result(i,j)=y[i].derivative(j);
    }
  }
  return result;
}


template<class R>
Function::InterpretedFunction<R>*
Function::InterpretedFunction<R>::clone() const 
{
  return new InterpretedFunction<R>(*this);
}


template<class R>
std::string
Function::InterpretedFunction<R>::name() const
{
  return this->_name;
}


template<class R>
void
Function::InterpretedFunction<R>::read(const std::string& filename)
{
  std::ifstream ifs(filename.c_str());
  if(ifs) {
    this->read(ifs);
  } else {
    ARIADNE_THROW(InvalidInput,"Function::read(string filename)",": Could not open file \""<<filename<<"\"");
  }
}


template<class R>
std::istream&
Function::InterpretedFunction<R>::read(std::istream& is)
{
  Input::ModelicaParser parser(is);

  this->_name=parser.function_name();
  this->_constants.assign(parser.constants().begin(),parser.constants().end());
  this->_variables.assign(parser.variables().begin(),parser.variables().end());
  this->_operations.assign(parser.operations().begin(),parser.operations().end());
  this->_initialise();
  return is;
}


template<class R>
std::ostream&
Function::InterpretedFunction<R>::write(std::ostream& os) const
{
  VirtualMachine::Index index;
  int value;

  os << "InterpretedFunction(\n" << std::flush;
  os << "  result_size=" << this->result_size() << ",\n" << std::flush;
  os << "  argument_size=" << this->argument_size() << ",\n" << std::flush;
  os << "  variables=" << this->_variables << ",\n" << std::flush;
  if(this->_constants_size>0) { os << "  constants=" << this->_constants << ",\n" << std::flush; }
  os << "  program=[ ";
  for(array<VirtualMachine::ByteCode>::const_iterator stack_iter=this->_operations.begin();
      stack_iter!=this->_operations.end(); ++stack_iter)
  {
    switch(stack_iter->op) {
    case VirtualMachine::PUSH:
      index=(++stack_iter)->ind;
      os << "PUSH[" << index[0] << "," << index[1] << "] ";
      break;
    case VirtualMachine::PULL:
      index=(++stack_iter)->ind;
      os << "PULL[" << index[0] << "," << index[1] << "] ";
      break;
    case VirtualMachine::CONST:
      value=(++stack_iter)->val;
      os << "CONST(" << value << ") ";
      break;
    case VirtualMachine::POS:
      os << "POS " << std::flush;
      break;
    case VirtualMachine::NEG:
      os << "NEG " << std::flush;
      break;
    case VirtualMachine::ADD:
      os << "ADD " << std::flush;
      break;
    case VirtualMachine::SUB:
      os << "SUB " << std::flush;
      break;
    case VirtualMachine::MUL:
      os << "MUL " << std::flush;
      break;
    case VirtualMachine::DIV:
      os << "DIV " << std::flush;
      break;
    case VirtualMachine::POW:
      value=(++stack_iter)->val;
      os << "POW(" << value << ") " << std::flush;
      break;
    case VirtualMachine::MIN:
      os << "MIN " << std::flush;
      break;
    case VirtualMachine::MAX:
      os << "MAX " << std::flush;
      break;
    case VirtualMachine::ABS:
      os << "ABS " << std::flush;
      break;
    case VirtualMachine::EXP:
      os << "EXP " << std::flush;
      break;
    case VirtualMachine::LOG:
      os << "LOG " << std::flush;
      break;
    case VirtualMachine::SIN:
      os << "SIN " << std::flush;
      break;
    case VirtualMachine::COS:
      os << "COS " << std::flush;
      break;
    case VirtualMachine::TAN:
      os << "TAN " << std::flush;
      break;
    case VirtualMachine::ASIN:
      os << "ASIN " << std::flush;
      break;
    case VirtualMachine::ACOS:
      os << "ACOS " << std::flush;
      break;
    case VirtualMachine::ATAN:
      os << "ATAN " << std::flush;
      break;
    }
  }
  os << "]\n)";
  return os;
}



template<class R>
smoothness_type
Function::InterpretedFunction<R>::smoothness() const
{
  //return std::numerical_limits<smoothness_type>::max();
  return 1; 
}


template<class R>
size_type
Function::InterpretedFunction<R>::argument_size() const
{
  return this->_argument_size;
}


template<class R>
size_type
Function::InterpretedFunction<R>::result_size() const
{
  return this->_result_size;
}







}















