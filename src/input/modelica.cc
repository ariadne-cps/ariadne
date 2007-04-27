/***************************************************************************
 *            modelica.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include "input/modelica.h"

#include <iostream>
#include <sstream>
#include <cassert>
#include <cstring>

#include "base/stlio.h"
#include "base/exceptions.h"
#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "output/logging.h"


namespace Ariadne {

using namespace System;

std::ostream& 
Input::operator<<(std::ostream& os, const Token& token)
{
  switch(token.type) {
  case Parser::IDNT: os << "\"" << token.identifier << "\""; break;
  case Parser::INT: os << token.integer; break;
  case Parser::RL: os << token.constant; break;
  
  case Parser::FUNCTION: os << "FUNCTION"; break; 
  case Parser::PROTECTED: os << "PROTECTED"; break; 
  case Parser::ALGORITHM: os << "ALGORITHM"; break; 
  case Parser::EQUATIONS: os << "EQUATIONS"; break; 
  case Parser::END: os << "END"; break; 
  case Parser::INTEGER: os << "INTEGER"; break; 
  case Parser::REAL: os << "REAL"; break; 
  case Parser::OUTPUT: os << "OUTPUT"; break; 
  case Parser::INPUT: os << "INPUT"; break; 
  case Parser::PARAMETER: os << "PARAMETER"; break; 
  case Parser::CONSTANT: os << "CONSTANT"; break; 

  case Parser::MIN: os << "MIN"; break;
  case Parser::MAX: os << "MAX"; break;
  case Parser::ABS: os << "ABS"; break;
  case Parser::EXP: os << "EXP"; break;
  case Parser::LOG: os << "LOG"; break;
  case Parser::SIN: os << "SIN"; break;
  case Parser::COS: os << "COS"; break;
  case Parser::TAN: os << "TAN"; break;
  case Parser::ASIN: os << "ASIN"; break;
  case Parser::ACOS: os << "ACOS"; break;
  case Parser::ATAN: os << "ATAN"; break;
  
  case Parser::EoF: os << "EoF"; break;
  case Parser::ASGN: os << "':='"; break;
  default: os << '\'' << char(token.type) << '\'';
  }
  return os;
}


std::ostream& 
operator<(std::ostream& os, const Input::Token& token)
{
  using namespace Input;
  switch(token.type) {
  case Parser::IDNT: os << token.identifier; break;
  case Parser::INT: os << token.integer; break;
  case Parser::RL: os << token.constant; break;
  
  case Parser::FUNCTION: os << "function "; break; 
  case Parser::PROTECTED: os << "protected\n"; break; 
  case Parser::ALGORITHM: os << "algorithm\n"; break; 
  case Parser::EQUATIONS: os << "equations\n"; break; 
  case Parser::END: os << "end "; break; 
  case Parser::INTEGER: os << "Integer"; break; 
  case Parser::REAL: os << "Real"; break; 
  case Parser::OUTPUT: os << "output "; break; 
  case Parser::INPUT: os << "input "; break; 
  case Parser::PARAMETER: os << "parameter "; break; 
  case Parser::CONSTANT: os << "constant "; break; 

  case Parser::MIN: os << "min"; break;
  case Parser::MAX: os << "max"; break;
  case Parser::ABS: os << "abs"; break;
  case Parser::EXP: os << "exp"; break;
  case Parser::LOG: os << "log"; break;
  case Parser::SIN: os << "sin"; break;
  case Parser::COS: os << "cos"; break;
  case Parser::TAN: os << "tan"; break;
  case Parser::ASIN: os << "asin"; break;
  case Parser::ACOS: os << "acos"; break;
  case Parser::ATAN: os << "atan"; break;
  
  case Parser::EoF: os << "\n"; break;
  case Parser::ASGN: os << ":="; break;
  default: os << char(token.type);
  }
  return os;
}


const std::string& 
Input::ModelicaParser::function_name() const
{
  return this->_name;
}



std::ostream& 
Input::ModelicaParser::write(std::ostream& os) const
{
  bool newline=true;
  for(std::vector<Token>::const_iterator tok_iter = this->_tokens.begin();
      tok_iter!=this->_tokens.end(); ++tok_iter)
  {
    const Token& tok=*tok_iter;
    if(tok.type==Parser::FUNCTION) {
      ++tok_iter; os << "function " << tok_iter->identifier << "\n"; newline=true;
    } else if (tok.type==Parser::PROTECTED || tok.type==Parser::ALGORITHM || tok.type==Parser::EQUATIONS) { 
      os < tok; newline=true;
    } else if (tok.type==Parser::END) {
      os < tok; newline=false;
    } 
    else if (tok.type==Parser::SC) {
      os << ";\n"; newline=true;
    } else {
      if(newline) { os << "  "; newline=false; }
      os < tok;
    }
    if((tok.type==Parser::RSP || tok.type==Parser::REAL) && tok_iter[1].type==Parser::IDNT) {
      os << " ";
    }
  }
  return os;
}





void
Input::ModelicaParser::tokenize(std::istream& is)
{
  Token tok=this->get_token(is);
  if(tok.type!=Parser::FUNCTION) {
    ARIADNE_THROW(SyntaxError,"ModelicaParser::tokenize(istream)",": expected \"function\" at beginning of function input");
  }
  
  while(!is.eof() && tok.type!=Parser::EoF) {
    this->_tokens.push_back(tok);
    tok=this->get_token(is);
  }

  
  ARIADNE_LOG(3,this->_tokens<<"\n"<<std::endl;);
  if(verbosity>=3) { this->write(std::clog)<<std::endl; }
}



Input::Token
Input::ModelicaParser::get_token(std::istream& is)
{
  Token result;
  char c;
  is >> c;
  if(is.eof()) {
    result=Token(Parser::EoF);
  } else if(isalpha(c)) {
    std::string str;
    while(isalpha(c) || isdigit(c) || c=='_') {
      str+=c;
      c=is.get();
    }
    is.putback(c);
    if(str=="EoF") {
      result=Token(Parser::EoF);
    } else if(str=="function") {
      result=Token(Parser::FUNCTION); 
    } else if(str=="protected") {
      result=Token(Parser::PROTECTED); 
    } else if(str=="algorithm") {
      result=Token(Parser::ALGORITHM); 
    } else if(str=="equations") {
      result=Token(Parser::EQUATIONS); 
    } else if(str=="end") {
      result=Token(Parser::END); 
    } else if(str=="input") {
      result=Token(Parser::INPUT); 
    } else if(str=="output") {
      result=Token(Parser::OUTPUT); 
    } else if(str=="parameter") {
      result=Token(Parser::PARAMETER); 
    } else if(str=="constant") {
      result=Token(Parser::CONSTANT); 
    } else if(str=="Real") {
      result=Token(Parser::REAL); 
    } else if(str=="min") {
      result=Token(Parser::MIN); 
    } else if(str=="max") {
      result=Token(Parser::MAX); 
    } else if(str=="abs") {
      result=Token(Parser::ABS);
    } else if(str=="exp") {
      result=Token(Parser::EXP);
    } else if(str=="log") {
      result=Token(Parser::LOG);
    } else if(str=="sin") {
      result=Token(Parser::SIN);
    } else if(str=="cos") {
      result=Token(Parser::COS);
    } else if(str=="tan") {
      result=Token(Parser::TAN);
    } else if(str=="asin") {
      result=Token(Parser::ASIN);
    } else if(str=="acos") {
      result=Token(Parser::ACOS);
    } else if(str=="atan") {
      result=Token(Parser::ATAN);
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
    case '(': result=Token(Parser::LP); break;
    case ')': result=Token(Parser::RP); break;
    case '[': result=Token(Parser::LSP); break;
    case ']': result=Token(Parser::RSP); break;
    case '{': result=Token(Parser::LCP); break;
    case '}': result=Token(Parser::RCP); break;
    case '=': result=Token(Parser::EQ); break;
    case ';': result=Token(Parser::SC); break;
    case '+': result=Token(Parser::PL); break;
    case '-': result=Token(Parser::MN); break;
    case '*': result=Token(Parser::TM); break;
    case '/': result=Token(Parser::DV); break;
    case '^': result=Token(Parser::PW); break;
    case ':':  
      c=is.peek(); 
      if(c=='=') { c=is.get(); result=Token(Parser::ASGN); break; } 
      else { ARIADNE_THROW(SyntaxError,"","Unrecognised token '"<<c<<"'"); }
      break;
    default: ARIADNE_THROW(SyntaxError,"","Unrecognised token beginning with '"<<c<<"'";);  
    }
  }
  return result;
}


Input::Token
Input::ModelicaParser::next_token()
{
  ++this->cursor;
  return *this->cursor;
}

Input::Token
Input::ModelicaParser::get_token()
{
  if(this->cursor==this->_tokens.end()) {
    return Token(Parser::EoF);
  }
  return *(this->cursor++);

}

Input::Token
Input::ModelicaParser::peek_token()
{
  return *this->cursor;
}



Input::ModelicaParser::ModelicaParser(std::istream& is)
{
  this->tokenize(is);
  this->cursor=this->_tokens.begin();
  this->read_function();
}



Input::ModelicaParser::ModelicaParser(const std::string& str)
{
  std::stringstream ss(str);
  this->tokenize(ss);
  this->cursor=this->_tokens.begin();
  this->read_function();
}




System::FunctionVariable
Input::ModelicaParser::variable(const std::string& name) const
{
  for(std::vector<FunctionVariable>::const_iterator var_iter=this->_variables.begin();
      var_iter!=this->_variables.end(); ++var_iter)
  {
    if(var_iter->name==name) {
      return *var_iter;
    }
  }
  ARIADNE_THROW(ParseError,""," variable \""<<name<<"\" not found");
}


void
Input::ModelicaParser::new_variable(const std::string& name, FunctionVariable::Type type, uint size)
{
  size_type start=0;
  for(std::vector<FunctionVariable>::const_iterator var_iter=this->_variables.begin(); 
      var_iter!=this->_variables.end(); ++var_iter)
  {
    if(var_iter->name==name) {
      ARIADNE_THROW(std::runtime_error,"new_variable","identifier of name \""<<name<<"\" already exists");
    }
    if(var_iter->type==type) {
      start+=var_iter->size;
    }
  }

  FunctionVariable var;
  var.name=name;
  var.type=type;
  var.start=start;
  if(size==0) {
    var.array_flag=false;
    var.size=1; 
  } else {
    var.array_flag=true;
    var.size=size;
  }

  this->_variables.push_back(var);
  ARIADNE_LOG(4,"new variable: " << var);
}



VirtualMachine::Index
Input::ModelicaParser::get_variable_index()
{
  Token tok;
  FunctionVariable var;
  VirtualMachine::Index index;
  Numeric::Rational q;
  
  tok=this->get_token();
  assert(tok.type==Parser::IDNT);
  var=this->variable(tok.identifier);
  index[0]=var.type;
  index[1]=var.start;
  if(var.array_flag) {
    tok=this->get_token();
    assert(tok.type==Parser::LSP);
    tok=this->get_token();
    if(tok.integer<0 || uint(tok.integer)>=var.size) {
      ARIADNE_THROW(ParseError,"read_variable","x["<<tok.integer<<"] for variable of size "<<var.size);
    }
    index[1]+=tok.integer;
    tok=this->get_token();
    assert(tok.type==Parser::RSP);
  } else {
    tok=this->peek_token();
    assert(tok.type!=Parser::LSP);
  }
  return index;
}






const std::vector<System::FunctionVariable>&
Input::ModelicaParser::variables() const
{
  return this->_variables;
}


const std::vector<Numeric::Rational>&
Input::ModelicaParser::constants() const
{
  return this->_constants;
}


const std::vector<VirtualMachine::ByteCode>&
Input::ModelicaParser::operations() const
{
  return this->_operations;
}


void
Input::ModelicaParser::read_function()
{
  Token tok;

  tok=this->get_token();
  if(tok.type!=Parser::FUNCTION) {
    ARIADNE_THROW(SyntaxError,"ModelicaParser::read_function()",": expected \"function\" at beginning of input; received "<<tok);
  }
  
  tok=this->get_token();
  if(tok.type!=Parser::IDNT) {
    ARIADNE_THROW(SyntaxError,"ModelicaParser::read_function()",": expected identifier");
  }
  this->_name=tok.identifier;
  
  tok=this->peek_token();
  while(tok.type!=Parser::ALGORITHM) {
    if(tok.type==Parser::PROTECTED) {
      tok=this->get_token();
      tok=this->peek_token();
    } else {
      this->read_declaration();
    }
    tok=*cursor;
  }
  // Read statements 
  tok=this->get_token();
  tok=this->peek_token();
  while(tok.type!= Parser::END) {
    this->read_statement();
    tok=this->peek_token();
  }
  // Read end-of-function token
  tok=this->get_token();
  if(tok.type!=Parser::END) {
    ARIADNE_THROW(SyntaxError,"ModelicaParser::read_function()",": expected \"end\" at end-of-function; obtained "<<tok);
  }
  tok=this->get_token();
  if(tok.type!=Parser::IDNT || tok.identifier!=this->_name) {
    ARIADNE_THROW(SyntaxError,"ModelicaParser::read_function()",": expected \""<<this->_name<<"\" after \"end\" at end-of-function; obtained "<<tok);
  }
  tok=this->get_token();
  if(tok.type!=Parser::SC) {
    ARIADNE_THROW(SyntaxError,"ModelicaParser::read_function()",": expected ';' at end of function input; obtained "<<tok);
  }
  }



void
Input::ModelicaParser::read_declaration()
{
  std::string name;
  FunctionVariable::Type type;
  Numeric::Rational value;
  uint size;
  Token tok=this->get_token();
  switch(tok.type) {
  case Parser::INPUT: 
    type=FunctionVariable::INPUT; tok=this->get_token(); break;
  case Parser::OUTPUT:
    type=FunctionVariable::OUTPUT; tok=this->get_token(); break;
  case Parser::CONSTANT:
    type=FunctionVariable::CONSTANT; tok=this->get_token(); break;
  default:
    type=FunctionVariable::INTERMEDIATE; 
  }
  switch(tok.type) {
  case Parser::REAL: 
    break;
  default:
    ARIADNE_THROW(SyntaxError,"read_declaration","expected \"Real\"; obtained "<<tok);
  }
  tok=this->peek_token(); 
  switch(tok.type) {
  case Parser::LSP: 
    tok=this->get_token();
    tok=this->get_token();
    if(tok.type!=Parser::INT) {
      ARIADNE_THROW(SyntaxError,"read_declaration","expected integer after '['; obtained "<<tok);
    }
    size=tok.integer;
    if(size<=0) {
      ARIADNE_THROW(SyntaxError,"read_declaration","variable size must be strictly positive");
    }
    tok=this->get_token();
    if(tok.type!=Parser::RSP) {
      ARIADNE_THROW(SyntaxError,"read_declaration","expected ']' after '['"<<size<<"; obtained "<<tok);
    }
    break;
  default:
    size=0;
  }
  tok=this->get_token();
  switch(tok.type) {
  case Parser::IDNT:
    name=tok.identifier; break;
  default:
    ARIADNE_THROW(SyntaxError,"read_declaration",": expected identifier; obtained "<<tok);
  }
  if(type==FunctionVariable::CONSTANT) {
    if(size!=0) {
      ARIADNE_THROW(SyntaxError,"read_declaration",": constant "<<name<<" declared as array");
    }
    tok=this->get_token();
    if(tok.type!=Parser::EQ) {
      ARIADNE_THROW(SyntaxError,"read_declaration",": expected '=' after constant-declaration; obtained "<<tok);
    }
    tok=this->get_token();
    switch(tok.type) {
    case Parser::INT: case Parser::RL:
      this->_constants.push_back(tok.constant);
      break;
    default:
      ARIADNE_THROW(SyntaxError,"read_declaration","expected constant after \"constant Real "<<name<<"=\"; obtained "<<tok);
    }
    this->new_variable(name,type,size);
  } else {
    this->new_variable(name,type,size);
  } 
  tok=this->get_token();
  switch(tok.type) {
  case Parser::SC: break;
  default: ARIADNE_THROW(SyntaxError,"read_declaration","expected ';' after variable-declaration; obtained "<<tok);
  }
}






void
Input::ModelicaParser::read_statement()
{
  VirtualMachine::ByteCode code;
  VirtualMachine::Index index=this->get_variable_index();
  
  this->read_assignment();
  
  this->read_expression();
  
  Token tok=this->get_token();
  if(tok.type!=Parser::SC) {
    ARIADNE_THROW(SyntaxError,"read_statement","expected ';' at end of statement");
  }

  code.op=VirtualMachine::PULL;
  this->_operations.push_back(code);
  code.ind=index;
  this->_operations.push_back(code);
}


void
Input::ModelicaParser::read_assignment()
{
  Token tok=this->get_token();
  if(tok.type!=Parser::ASGN && tok.type!=Parser::EQ) {
    ARIADNE_THROW(SyntaxError,"read_assignement","expected '=' or ':=' as assignment");
  }
}
  


void
Input::ModelicaParser::read_expression()
{
  // expression: 
  //   { [ PLUS | MINUS ] } term ( [ PLUS | MINUS ] term )*

  ARIADNE_LOG(4,"read_expression");
  VirtualMachine::ByteCode code;
  Token tok;
  read_term();
  tok=this->peek_token();
  while(tok.type == Parser::PL || tok.type == Parser::MN) {
    tok=this->get_token();
    code.op = (tok.type==Parser::PL) ? VirtualMachine::ADD : VirtualMachine::SUB;
    read_term();
    this->_operations.push_back(code);
    tok=this->peek_token();
  }
}


void
Input::ModelicaParser::read_term()
{
  //   term: factor ( [ TIMES | DIVIDES ] factor ) *

  ARIADNE_LOG(4,"read_term");
  VirtualMachine::ByteCode code;
  Token tok;
  read_factor();
  tok=this->peek_token();
  while(tok.type == Parser::TM || tok.type == Parser::DV) {
    tok=this->get_token();
    code.op = (tok.type==Parser::TM) ? VirtualMachine::MUL : VirtualMachine::DIV;
    read_factor();
    this->_operations.push_back(code);
    tok=this->peek_token();
  }
}


void
Input::ModelicaParser::read_factor()
{
  // factor:  { PLUS | MINUS } primary | primary ^ integer

  ARIADNE_LOG(4,"read_factor");
  VirtualMachine::ByteCode pncode;
  VirtualMachine::ByteCode code;
  Token tok;
  Variable var;
  Numeric::Rational q;
  
  tok=this->peek_token();
  if(tok.type==Parser::PL || tok.type==Parser::MN) {
    this->get_token();
    pncode.op=(tok.type==Parser::PL) ? VirtualMachine::POS : VirtualMachine::NEG;
  } else {
    pncode.op=VirtualMachine::POS;
  }
  read_primary();
  tok=this->peek_token();
  if(tok.type==Parser::PW) {
    tok=this->get_token();
    code.op=VirtualMachine::POW;
    this->_operations.push_back(code);
    tok=this->get_token();
    if(tok.type==Parser::INT) {
      code.val=tok.integer;
      this->_operations.push_back(code);
    } else {
      ARIADNE_THROW(SyntaxError,"read_factor","expected integer after '^'");
    }
  }
  if(pncode.op==VirtualMachine::NEG) {
      this->_operations.push_back(pncode);
  }
    
}



void
Input::ModelicaParser::read_primary()
{
  // primary:  variable | constant | function(expression) | (expression)
  ARIADNE_LOG(4,"read_primary");
  VirtualMachine::ByteCode code;
  Token tok;
  Variable var;
  Numeric::Rational q;
  
  tok=this->peek_token();
  switch(tok.type) {
  case Parser::IDNT: 
    code.op=VirtualMachine::PUSH;
    this->_operations.push_back(code);
    read_variable();
    break;
  case Parser::INT: case Parser::RL: 
    read_constant();
    break;
  case Parser::EXP: case Parser::LOG:
  case Parser::SIN: case Parser::COS: case Parser::TAN: 
  case Parser::ASIN: case Parser::ACOS: case Parser::ATAN:
    get_token();
    switch(tok.type) {
    case Parser::EXP: code.op=VirtualMachine::EXP; break;
    case Parser::LOG: code.op=VirtualMachine::LOG; break;
    case Parser::SIN: code.op=VirtualMachine::SIN; break;
    case Parser::COS: code.op=VirtualMachine::COS; break;
    case Parser::TAN: code.op=VirtualMachine::TAN; break;
    case Parser::ASIN: code.op=VirtualMachine::ASIN; break;
    case Parser::ACOS: code.op=VirtualMachine::ACOS; break;
    case Parser::ATAN: code.op=VirtualMachine::ATAN; break;
    default: { }
    }
    tok=this->get_token();
    if(tok.type!=Parser::LP) {
      ARIADNE_THROW(SyntaxError,"read_factor","expected '(' after function");
    }
    this->read_expression();
    tok=this->get_token();
    if(tok.type!=Parser::RP) {
      ARIADNE_THROW(SyntaxError,"read_factor","expected ')' to close function argument");
    }
    this->_operations.push_back(code);
    break;
  case Parser::LP:
    tok=this->get_token();
    this->read_expression();
    tok=this->get_token();
    if(tok.type!=Parser::RP) {
      ARIADNE_THROW(SyntaxError,"read_expression","expected ')'");
    }
    break;
  default:
    ARIADNE_THROW(SyntaxError,"read_primary","tok=" << tok << "; expected variable | constant | function(expression) | (expression)");
  } 
}

void
Input::ModelicaParser::read_variable()
{
  ARIADNE_LOG(4,"read_variable");
  VirtualMachine::ByteCode code;
  code.ind=this->get_variable_index();
  this->_operations.push_back(code);
}

void
Input::ModelicaParser::read_constant()
{
  ARIADNE_LOG(4,"read_constant");
  VirtualMachine::ByteCode code;
  Token tok=get_token();
  assert(tok.type==Parser::INT || tok.type==Parser::RL);
  code.op=VirtualMachine::CONST;
  this->_operations.push_back(code);
  if(tok.type==Parser::RL) {
    ARIADNE_THROW(InvalidInput,"ModelicaParser::read_constant()",": Non-integer undeclared constants not currently supported");
  } else {
    code.val=tok.integer;
  }
  this->_operations.push_back(code);
}



}
