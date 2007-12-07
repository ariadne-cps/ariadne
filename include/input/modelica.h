/***************************************************************************
 *            modelica.h
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
 
/*! \file system/modelica.h
 *  \brief Modelica function parser.
 */
 
#ifndef ARIADNE_INPUT_MODELICA_H
#define ARIADNE_INPUT_MODELICA_H

#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

#include "../base/types.h"
#include "../numeric/traits.h"
#include "../numeric/rational.h"
#include "../function/interpreted_function.h"
#include "../function/virtual_machine.h"

namespace Ariadne {
  namespace Input {

    class SyntaxError : public std::runtime_error { public: SyntaxError(const std::string& what) : std::runtime_error(what) { } };
    class ParseError : public std::runtime_error { public: ParseError(const std::string& what) : std::runtime_error(what) { } };

    namespace Parser {
      enum TokenType { IDNT, INT, RL, 
                       SC=';', LP='(', RP=')', LSP='[', RSP=']', LCP='{', RCP='}', EQ='=', PL='+', MN='-', TM='*', DV='/', PW='^', ASGN, EoF,
                       FUNCTION, PROTECTED, ALGORITHM, EQUATIONS, END, 
                       INTEGER, REAL, 
                       OUTPUT, INPUT, PARAMETER, CONSTANT,
                       MIN, MAX, ABS, EXP, LOG, SIN, COS, TAN, ASIN, ACOS, ATAN };
    
    }

    /*! \brief A token in the Modelica language. */
    struct Token { 
      Parser::TokenType type; 
      std::string identifier;
      int integer; 
      Numeric::Rational constant; 
      
      Token() { }
      Token(Parser::TokenType t) : type(t) { }
      Token(std::string s) : type(Parser::IDNT), identifier(s) { }
      Token(int n) : type(Parser::INT), integer(n), constant(n) { }
      Token(double x) : type(Parser::RL), constant(x) { }
      
      bool operator==(Parser::TokenType t) { return this->type==t; }
    };
      

    /*! \brief A parser for a function. */
    class ModelicaParser {
     public:
      
      //const char * functions[3]={ "function", "end", "algorithm" };
      
      ModelicaParser(std::istream&);
      ModelicaParser(const std::string&);
      
      const std::string& function_name() const;
      const std::vector<Function::VirtualMachine::ByteCode>& operations() const;
      const std::vector<Function::FunctionVariable>& variables() const;
      const std::vector<Numeric::Rational>& constants() const;

     private:
      void tokenize(std::istream& is);
      Token get_token(std::istream& is);

      void new_variable(const std::string& name, Function::FunctionVariable::Type type, uint size);
      Function::FunctionVariable variable(const std::string& name) const;
      
      Token peek_token();
      Token next_token();
      Token get_token();
      Token token();

      Function::VirtualMachine::Index get_variable_index();

      void read_function();
      void read_declaration();
      void read_statement();
      void read_assignment();
      void read_expression();
      void read_term();
      void read_factor();
      void read_primary();
      void read_variable();
      void read_constant();

      /*! \brief Write the Modelica function to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      std::string _name;
      std::vector<Token> _tokens;
      std::vector<Function::FunctionVariable> _variables;
      std::vector<Numeric::Rational> _constants;
      std::vector<Function::VirtualMachine::ByteCode> _operations;
      
      std::vector<Token>::const_iterator cursor;
    };
  
    std::ostream& operator<<(std::ostream& os, const Token& tok);

  }
}

#endif /* ARIADNE_INPUT_MODELICA_H */
