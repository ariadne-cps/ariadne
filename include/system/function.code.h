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

namespace {

using namespace Ariadne;
using namespace Ariadne::System;

template<class R>
void evaluate(const std::vector<System::ByteCode>& ops, R** args)
{
  std::vector<R> stack;
  int v;
  int i;
  for(std::vector<System::ByteCode>::const_iterator op_iter=ops.begin();
      op_iter!=ops.end(); ++op_iter)
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
    case POS:
      stack[stack.size()-1]=(+stack[stack.size()-1]);
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
    case POW:
      // FIXME: Implement this!
      break;
    case MIN:
      stack[stack.size()-2]=Numeric::min(stack[stack.size()-2],stack[stack.size()-1]);
      stack.pop_back();
      break;
    case MAX:
      stack[stack.size()-2]=Numeric::max(stack[stack.size()-2],stack[stack.size()-1]);
      stack.pop_back();
      break;
    case ABS:
      stack[stack.size()-1]=Numeric::abs(stack[stack.size()-1]);
      break;
    case EXP:
      stack[stack.size()-1]=Numeric::exp(stack[stack.size()-1]);
      break;
    case LOG:
      stack[stack.size()-1]=Numeric::log(stack[stack.size()-1]);
      break;
    case SIN:
      stack[stack.size()-1]=Numeric::sin(stack[stack.size()-1]);
      break;
    case COS:
      stack[stack.size()-1]=Numeric::cos(stack[stack.size()-1]);
      break;
    case TAN:
      stack[stack.size()-1]=Numeric::tan(stack[stack.size()-1]);
      break;
    case ASIN:
      stack[stack.size()-1]=Numeric::asin(stack[stack.size()-1]);
      break;
    case ACOS:
      stack[stack.size()-1]=Numeric::acos(stack[stack.size()-1]);
      break;
    case ATAN:
      stack[stack.size()-1]=Numeric::atan(stack[stack.size()-1]);
      break;
    }
  }
}


template<>
void evaluate<Numeric::Rational>(const std::vector<System::ByteCode>& ops, Numeric::Rational** args)
{
  std::vector<Numeric::Rational> stack;
  int v;
  int i;
  for(std::vector<System::ByteCode>::const_iterator op_iter=ops.begin();
      op_iter!=ops.end(); ++op_iter)
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
    case POS:
      stack[stack.size()-1]=(+stack[stack.size()-1]);
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
    case POW:
      // FIXME: Implement this!
      break;
    case MIN:
      stack[stack.size()-2]=Numeric::min(stack[stack.size()-2],stack[stack.size()-1]);
      stack.pop_back();
      break;
    case MAX:
      stack[stack.size()-2]=Numeric::max(stack[stack.size()-2],stack[stack.size()-1]);
      stack.pop_back();
      break;
    case ABS:
      stack[stack.size()-1]=Numeric::abs(stack[stack.size()-1]);
      break;
    default:
      throw std::runtime_error("evaluate: operation not permitted for rational numbers");
    }
  }
}


}



namespace Ariadne {

namespace System {
std::ostream& 
operator<<(std::ostream& os, const Variable& var) 
{
  return os << var.name << "(" << var.num << ")[" << var.ind << "]";
}
std::ostream& 
operator<<(std::ostream& os, const FunctionParser::TokenType& type) 
{
  switch(type) {
  case FunctionParser::STR: os << "STR"; break;
  case FunctionParser::NUM: os << "NUM"; break;
  case FunctionParser::INT: os << "INT"; break;
  case FunctionParser::REAL: os << "REAL"; break;
  case FunctionParser::END: os << "END"; break;
  case FunctionParser::ASSIGN: os << "ASSIGN"; break;
  case FunctionParser::LP: os << "LP"; break;
  case FunctionParser::RP: os << "RP"; break;
  case FunctionParser::LSP: os << "LSP"; break;
  case FunctionParser::RSP: os << "RSP"; break;
  case FunctionParser::PLUS: os << "PLUS"; break;
  case FunctionParser::MINUS: os << "MINUS"; break;
  case FunctionParser::TIMES: os << "TIMES"; break;
  case FunctionParser::DIVIDES: os << "DIVIDES"; break;
  case FunctionParser::POWER: os << "POWER"; break;
  case FunctionParser::MIN: os << "MIN"; break;
  case FunctionParser::MAX: os << "MAX"; break;
  case FunctionParser::ABS: os << "ABS"; break;
  case FunctionParser::EXP: os << "EXP"; break;
  case FunctionParser::LOG: os << "LOG"; break;
  case FunctionParser::SIN: os << "SIN"; break;
  case FunctionParser::COS: os << "COS"; break;
  case FunctionParser::TAN: os << "TAN"; break;
  case FunctionParser::ASIN: os << "ASIN"; break;
  case FunctionParser::ACOS: os << "ACOS"; break;
  case FunctionParser::ATAN: os << "ATAN"; break;
  }
  return os << "(" << int(type) << ")";
}

std::ostream& 
operator<<(std::ostream& os, const FunctionParser::Token& tok) 
{
  os << "Token( type=" << tok.type;
  if(tok.type==FunctionParser::STR) { os << ", str=\"" << tok.str<<"\""; }
  if(tok.type==FunctionParser::INT) { os << ", int=" << tok.val; }
  if(tok.type==FunctionParser::REAL) { os << ", rl=" << tok.val; }
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
std::istream&
System::GeneralFunction<R>::read(std::istream& is)
{
  FunctionParser parser(is);
  this->_variable_names=parser.variable_names();
  this->_operations=parser.operations();

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
    case POS:
      os << "POS " << std::flush;
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
    case POW:
      os << "POW " << std::flush;
      break;
    case MIN:
      os << "MIN " << std::flush;
      break;
    case MAX:
      os << "MAX " << std::flush;
      break;
    case ABS:
      os << "ABS " << std::flush;
      break;
    case EXP:
      os << "EXP " << std::flush;
      break;
    case LOG:
      os << "LOG " << std::flush;
      break;
    case SIN:
      os << "SIN " << std::flush;
      break;
    case COS:
      os << "COS " << std::flush;
      break;
    case TAN:
      os << "TAN " << std::flush;
      break;
    case ASIN:
      os << "ASIN " << std::flush;
      break;
    case ACOS:
      os << "ACOS " << std::flush;
      break;
    case ATAN:
      os << "ATAN " << std::flush;
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
  ::evaluate(this->_operations, args);
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
  return 3;
}


template<class R>
dimension_type
System::GeneralFunction<R>::result_dimension() const
{
  return 2;
}





System::FunctionParser::FunctionParser(std::istream& is)
{
  this->read_function(is);
}



System::FunctionParser::FunctionParser(const std::string& str)
{
  std::stringstream ss(str);
  this->read_function(ss);
}






System::FunctionParser::Token
System::FunctionParser::peek_token(std::istream& is)
{
  Token result;
  char c;
  is>>c;
  is.putback(c);
  if(is.eof()) {
    c=';';
  }
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
    case '^': result=Token(POWER); break;
    default: ARIADNE_THROW(SyntaxError,"peek_token","Unknown token type '"<<c<<"'");
    }
  }
  std::cerr << "Peek" << result << std::endl;
  return result;
}




System::FunctionParser::Token
System::FunctionParser::get_token(std::istream& is)
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
    if(str=="min") {
      result=Token(MIN); 
    } else if(str=="max") {
      result=Token(MAX); 
    } else if(str=="abs") {
      result=Token(ABS);
    } else if(str=="exp") {
      result=Token(EXP);
    } else if(str=="log") {
      result=Token(LOG);
    } else if(str=="sin") {
      result=Token(SIN);
    } else if(str=="cos") {
      result=Token(COS);
    } else if(str=="tan") {
      result=Token(TAN);
    } else if(str=="asin") {
      result=Token(ASIN);
    } else if(str=="acos") {
      result=Token(ACOS);
    } else if(str=="atan") {
      result=Token(ATAN);
    } else {
      result=Token(str);
    }
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




std::string
System::FunctionParser::read_string(std::istream& is)
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


int
System::FunctionParser::read_int(std::istream& is)
{
  int n;
  is >> n;
  return n;
}



System::Variable
System::FunctionParser::get_variable(std::istream& is)
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





uint
System::FunctionParser::variable_number(const std::string& name) const
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



std::vector<std::string>
System::FunctionParser::variable_names()
{
  return this->_variable_names;
}


std::vector<System::ByteCode>
System::FunctionParser::operations()
{
  return this->_operations;
}


void
System::FunctionParser::read_function(std::istream& is)
{
  this->_variable_names.push_back("y");
  this->_variable_names.push_back("x");
  this->_variable_names.push_back("p");
  this->_variable_names.push_back("c");

  // Read statements 
  Token tok=this->peek_token(is);
  while(tok.type!=END) {
    this->read_statement(is);
    tok=this->peek_token(is);
  }
  // Read end-of-function token
  this->get_token(is);
}



void
System::FunctionParser::read_statement(std::istream& is)
{
  ByteCode code;
  Variable lhs=this->get_variable(is);
  
  this->read_assignment(is);
  
  this->read_expression(is);
  
  Token tok=this->get_token(is);
  if(tok.type!=END) {
    ARIADNE_THROW(SyntaxError,"read_statement","expected ';' at end of statement");
  }

  code.op=PULL;
  this->_operations.push_back(code);
  code.var=lhs.num;
  this->_operations.push_back(code);
  code.ind=lhs.ind;
  this->_operations.push_back(code);

}


void
System::FunctionParser::read_assignment(std::istream& is)
{
  char c;
  is >> c;
  if(c!='=') {
    ARIADNE_THROW(SyntaxError,"read_assignement","expected '=' assignment");
  }
}
  


void
System::FunctionParser::read_expression(std::istream& is)
{
  // expression: 
  //   term ( [ PLUS | MINUS ] term )*

  std::cerr << "read_expression"<<std::endl;
  ByteCode code;
  Token tok;
  read_term(is);
  tok=this->peek_token(is);
  while(tok.type == PLUS || tok.type == MINUS) {
    tok=this->get_token(is);
    code.op = (tok.type==PLUS) ? ADD : SUB;
    read_term(is);
    this->_operations.push_back(code);
    tok=this->peek_token(is);
  }
}


void
System::FunctionParser::read_term(std::istream& is)
{
  //   term: factor ( [ TIMES | DIVIDES ] factor ) *

  std::cerr << "read_term"<<std::endl;
  ByteCode code;
  Token tok;
  read_factor(is);
  tok=this->peek_token(is);
  while(tok.type == TIMES || tok.type == DIVIDES) {
    tok=this->get_token(is);
    code.op = (tok.type==TIMES) ? MUL : DIV;
    read_factor(is);
    this->_operations.push_back(code);
    tok=this->peek_token(is);
  }
}


void
System::FunctionParser::read_factor(std::istream& is)
{
  std::cerr << "read_factor"<<std::endl;
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
  case EXP: case LOG:
  case SIN: case COS: case TAN: 
  case ASIN: case ACOS: case ATAN:
    switch(tok.type) {
    case EXP: code.op=System::EXP; break;
    case LOG: code.op=System::LOG; break;
    case SIN: code.op=System::SIN; break;
    case COS: code.op=System::COS; break;
    case TAN: code.op=System::TAN; break;
    case ASIN: code.op=System::ASIN; break;
    case ACOS: code.op=System::ACOS; break;
    case ATAN: code.op=System::ATAN; break;
    default: { }
    }
    tok=this->get_token(is);
    if(tok.type!=LP) {
      ARIADNE_THROW(SyntaxError,"read_factor","expected '(' after function");
    }
    this->read_expression(is);
    tok=this->get_token(is);
    if(tok.type!=RP) {
      ARIADNE_THROW(SyntaxError,"read_factor","expected ')' to close function argument");
    }
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
  tok=this->peek_token(is);
  if(tok.type==POWER) {
    tok=this->get_token(is);
    code.op=POW;
    this->_operations.push_back(code);
    tok=this->get_token(is);
    if(tok.type==INT) {
      code.val=code.val;
      this->_operations.push_back(code);
    } else {
      ARIADNE_THROW(SyntaxError,"read_factor","expected integer after '^'");
    }
  }
}



void
System::FunctionParser::read_primary(std::istream& is)
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









}

