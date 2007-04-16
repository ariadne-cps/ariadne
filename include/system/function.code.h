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
#include <cassert>
#include <cstring>

#include "../base/exceptions.h"
#include "../numeric/rational.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne {

namespace System {
std::ostream& 
operator<<(std::ostream& os, const Variable& var) 
{
  return os << var.name << "(" << var.num << ")[" << var.ind << "]";
}
std::ostream& 
operator<<(std::ostream& os, const TokenType& type) 
{
  switch(type) {
  case STR: os << "STR"; break;
  case NUM: os << "NUM"; break;
  case INT: os << "INT"; break;
  case REAL: os << "REAL"; break;
  case END: os << "END"; break;
  case ASSIGN: os << "ASSIGN"; break;
  case LP: os << "LP"; break;
  case RP: os << "RP"; break;
  case LSP: os << "LSP"; break;
  case RSP: os << "RSP"; break;
  case PLUS: os << "PLUS"; break;
  case MINUS: os << "MINUS"; break;
  case TIMES: os << "TIMES"; break;
  case DIVIDES: os << "DIVIDES"; break;
  }
  return os << "(" << int(type) << ")";
}
std::ostream& 
operator<<(std::ostream& os, const Token& tok) 
{
  os << "Token( type=" << tok.type;
  if(tok.type==STR) { os << ", str=\"" << tok.str<<"\""; }
  if(tok.type==INT) { os << ", int=" << tok.val; }
  if(tok.type==REAL) { os << ", rl=" << tok.val; }
  os << " )";
  return os;
}
}

template<class R>
LinearAlgebra::Vector<typename System::GeneralFunction<R>::A>
System::GeneralFunction<R>::image(const LinearAlgebra::Vector<A>& x) const
{
  LinearAlgebra::Vector<A> result(this->result_dimension());
  A* res_ptr=result.begin();
  A* arg_ptr=const_cast<A*>(x.begin());
  A* args[2] = { res_ptr, arg_ptr };
  this->_evaluate(args);
  return result;
}


template<class R>
typename System::GeneralFunction<R>::A
System::GeneralFunction<R>::derivative(const LinearAlgebra::Vector<A>& x, const size_type& i, const LinearAlgebra::MultiIndex& j) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
 LinearAlgebra::Matrix<typename System::GeneralFunction<R>::A>
System::GeneralFunction<R>::jacobian(const LinearAlgebra::Vector<A>& x) const
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
System::GeneralFunction<R>*
System::GeneralFunction<R>::clone() const 
{
  return new GeneralFunction<R>(*this);
}


template<class R>
System::Token
System::GeneralFunction<R>::peek_token(std::istream& is)
{
  Token result;
  char c;
  is>>c;
  is.putback(c);
  //std::cerr << "peek_token: c='" << c << "'" << std::endl;
  if(isalpha(c)) {
    result=Token(STR);
  } else if(isdigit(c)) {
    result=Token(NUM);
  } else {
    switch(c) {
    case '(': result=Token(LP); break;
    case ')': result=Token(RP); break;
    case '[': result=Token(LSP); break;
    case ']': result=Token(RSP); break;
    case '=': result=Token(ASSIGN); break;
    case ';': result=Token(END); break;
    case '+': result=Token(PLUS); break;
    case '-': result=Token(MINUS); break;
    case '*': result=Token(TIMES); break;
    case '/': result=Token(DIVIDES); break;
    default: assert(false);
    }
  }
  std::cerr << "Peek" << result << std::endl;
  return result;
}



template<class R>
System::Token
System::GeneralFunction<R>::get_token(std::istream& is)
{
  Token result;
  char c;
  is >> c;
  //std::cerr << "get_token: c='" << c << "'" << std::endl;
  if(is.eof()) {
    std::cerr << "end_of_file reached" << std::endl;
    result=Token(END);
  } else if(isalpha(c)) {
    std::string str;
    while(isalpha(c)) {
      str+=c;
      c=is.get();
    }
    is.putback(c);
    result=Token(str);
  } else if(isdigit(c) || c=='.') {
    int n;
    is.putback(c);
    is >> n;
    c=is.peek();
    if(c=='.') {
      double x;
      is >> x;
      x+=n;
      result=Token(x);
    } else {
      result=Token(n);
    }
  } else {
    switch(c) {
    case '(': result=Token(LP); break;
    case ')': result=Token(RP); break;
    case '[': result=Token(LSP); break;
    case ']': result=Token(RSP); break;
    case '=': result=Token(ASSIGN); break;
    case ';': result=Token(END); break;
    case '+': result=Token(PLUS); break;
    case '-': result=Token(MINUS); break;
    case '*': result=Token(TIMES); break;
    case '/': result=Token(DIVIDES); break;
    default: assert(false);
    }
  }
  std::cerr << result << std::endl;
  return result;
}



template<class R>
std::string
System::GeneralFunction<R>::read_string(std::istream& is)
{
  std::string result;
  char c=is.peek();
  while(isalpha(c)) {
    c=is.get();
    result+=c;
    c=is.peek();
  }
  return result;
}

template<class R>
int
System::GeneralFunction<R>::read_int(std::istream& is)
{
  int n;
  is >> n;
  return n;
}


template<class R>
System::Variable
System::GeneralFunction<R>::get_variable(std::istream& is)
{
  Variable res;
  Token tok;
  res.num=0;
  res.ind=0;
  std::string name;
  char c;
  tok=this->get_token(is);
  res.name=tok.str;
  for(std::vector<std::string>::const_iterator names_iter=this->_variable_names.begin();
      names_iter!=this->_variable_names.end(); ++names_iter)
  {
    if(*names_iter==res.name) {
      break;
    }
    ++res.num;
  }
  if(res.num==this->_variable_names.size()) {
    ARIADNE_THROW(SyntaxError,"read_variable","unknown variable: "<<name);
  }
  if(is.peek()=='[') {
    tok=this->get_token(is);
    is >> res.ind;
    is >> c;
    if(c!=']') {
      ARIADNE_THROW(SyntaxError,"read_variable","expected ']' after '" << name << '[' << res.ind);
    }
  }
  std::cerr << res << "\n";
  //std::cerr << "Variable: " << res << std::endl;
  return res;
}



template<class R>
void
System::GeneralFunction<R>::read_assignment(std::istream& is)
{
  char c;
  is >> c;
  if(c!='=') {
    ARIADNE_THROW(SyntaxError,"read_assignement","expected '=' assignment");
  }
}
  

template<class R>
uint
System::GeneralFunction<R>::variable_number(const std::string& name) const
{
  uint num=0;
  for(std::vector<std::string>::const_iterator names_iter=this->_variable_names.begin();
      names_iter!=this->_variable_names.end(); ++names_iter)
  {
    if(*names_iter==name) {
      break;
    }
    ++num;
  }
  if(num==this->_variable_names.size()) {
    ARIADNE_THROW(SyntaxError,"read_variable","unknown variable: "<<name);
  }
  return num;
}




template<class R>
void
System::GeneralFunction<R>::read_statement(std::istream& is)
{
  ByteCode code;
  Variable lhs=this->get_variable(is);
  
  this->read_assignment(is);
  
  this->read_expression(is);
  
  code.op=PULL;
  this->_operations.push_back(code);
  code.var=lhs.num;
  this->_operations.push_back(code);
  code.ind=lhs.ind;
  this->_operations.push_back(code);

}


template<class R>
void
System::GeneralFunction<R>::read_expression(std::istream& is)
{
  // FIXME: Syntax should be 
  // expression: 
  //   term ( [ PLUS | MINUS ] term )*

  // term: 
  //   term
  //   expression+term
  //   expression-term

  std::cerr << "read_expression"<<std::endl;
  ByteCode code;
  Token tok;
  read_term(is);
  tok=this->peek_token(is);
  switch(tok.type) {
  case PLUS:
    tok=this->get_token(is);
    code.op=ADD;
    read_expression(is);
    this->_operations.push_back(code);
    break;
  case MINUS:
    tok=this->get_token(is);
    code.op=SUB;
    read_expression(is);
    this->_operations.push_back(code);
    break;
  default: 
    { }
  }
}

template<class R>
void
System::GeneralFunction<R>::read_term(std::istream& is)
{
  // FIXME: Syntax should be 
  //   term: primary ( [ TIMES | DIVIDES ] primary ) *

  // term: 
  //   primary*term
  //   primary/term
  //   primary

  std::cerr << "read_term"<<std::endl;
  ByteCode code;
  
  Token tok;
  read_primary(is);
  tok=this->peek_token(is);
  switch(tok.type) {
  case TIMES:
    tok=this->get_token(is);
    code.op=MUL;
    read_term(is);
    this->_operations.push_back(code);
    break;
  case DIVIDES:
    tok=this->get_token(is);
    code.op=DIV;
    read_term(is);
    this->_operations.push_back(code);
  default:
    { }
  }
}

template<class R>
void
System::GeneralFunction<R>::read_factor(std::istream& is)
{
  std::cerr << "read_factor"<<std::endl;
} 



template<class R>
void
System::GeneralFunction<R>::read_primary(std::istream& is)
{
  std::cerr << "read_primary"<<std::endl;
  ByteCode code;
  Token tok;
  Variable var;
  int n;
  
  tok=this->get_token(is);
  switch(tok.type) {
  case STR: 
    var.name=tok.str;
    var.num=this->variable_number(tok.str);
    tok=this->peek_token(is);
    if(tok.type==LSP) {
      tok=this->get_token(is);
      tok=this->get_token(is);
      var.ind=tok.val;
      tok=this->get_token(is);
      assert(tok.type==RSP);
    } else {
      var.ind=0;
    }
    code.op=PUSH;
    this->_operations.push_back(code);
    code.var=var.num;
    this->_operations.push_back(code);
    code.ind=var.ind;
    this->_operations.push_back(code);
    break;
  case NUM: 
    n=tok.val;
    code.op=CONST;
    this->_operations.push_back(code);
    code.val=n;
    this->_operations.push_back(code);
    break;
  case LP:
    this->read_expression(is);
    tok=this->get_token(is);
    if(tok.type!=RP) {
      ARIADNE_THROW(SyntaxError,"read_expression","expected ')'");
    }
  default:
    ;//ARIADNE_THROW(SyntaxError,"read_primary","tok=" << tok << "; expected variable | constant | (expression)");
  }
}


template<class R>
std::istream&
System::GeneralFunction<R>::read(std::istream& is)
{
  this->_variable_names.push_back("y");
  this->_variable_names.push_back("x");
  this->_variable_names.push_back("p");
  this->_variable_names.push_back("c");

  this->read_statement(is);

  return is;
}




template<class R>
std::ostream&
System::GeneralFunction<R>::write(std::ostream& os) const
{
  os << "Function( " << std::flush;
  os << "program_size=" << this->_operations.size() << ",   " << std::flush;
  for(std::vector<ByteCode>::const_iterator stack_iter=this->_operations.begin();
      stack_iter!=this->_operations.end(); ++stack_iter)
  {
    switch(stack_iter->op) {
    case PUSH:
      os << "PUSH(" << std::flush;
      os << this->_variable_names[(++stack_iter)->var];
      os << "[" << (++stack_iter)->ind << "]) ";
      break;
    case PULL:
      os << "PULL(" << std::flush;
      os << this->_variable_names[(++stack_iter)->var];
      os << "[" << (++stack_iter)->ind << "]) ";
      break;
    case CONST:
      os << "CONST " << std::flush;
      os << stack_iter->val;
      break;
    case NEG:
      os << "NEG " << std::flush;
      break;
    case ADD:
      os << "ADD " << std::flush;
      break;
    case SUB:
      os << "SUB " << std::flush;
      break;
    case MUL:
      os << "MUL " << std::flush;
      break;
    case DIV:
      os << "DIV " << std::flush;
      break;
    }
  }
  os << ")";
  return os;
}


template<class R>
System::GeneralFunction<R>::GeneralFunction(const std::string& str)
{
  std::stringstream ss(str);
  this->read(ss);
}

template<class R>
void 
System::GeneralFunction<R>::_evaluate(A** args) const
{
  
  std::vector<A> stack;
  int v;
  int i;
  for(std::vector<ByteCode>::const_iterator op_iter=this->_operations.begin();
      op_iter!=this->_operations.end(); ++op_iter)
  {
    switch(op_iter->op) {
    case PUSH:
      v=(++op_iter)->var;
      i=(++op_iter)->ind;
      stack.push_back(args[v][i]); 
      break;
    case PULL:
      v=(++op_iter)->var;
      i=(++op_iter)->ind;
      args[v][i]=stack.back();
      stack.pop_back();
      break;
    case CONST:
      stack.push_back((++op_iter)->val);
      break;
    case NEG:
      stack[stack.size()-1]=(-stack[stack.size()-1]);
      break;
    case ADD:
      stack[stack.size()-2]=stack[stack.size()-2]+stack[stack.size()-1];
      stack.pop_back();
      break;
    case SUB:
      stack[stack.size()-2]=stack[stack.size()-2]-stack[stack.size()-1];
      stack.pop_back();
      break;
    case MUL:
      stack[stack.size()-2]=stack[stack.size()-2]*stack[stack.size()-1];
      stack.pop_back();
      break;
    case DIV:
      stack[stack.size()-2]=stack[stack.size()-2]/stack[stack.size()-1];
      stack.pop_back();
      break;
    }
  }
}


template<class R>
size_type
System::GeneralFunction<R>::smoothness() const
{
  //return std::numerical_limits<size_type>::max();
  return 1; 
}


template<class R>
dimension_type
System::GeneralFunction<R>::argument_dimension() const
{
  return 2;
}


template<class R>
dimension_type
System::GeneralFunction<R>::result_dimension() const
{
  return 1;
}



}

