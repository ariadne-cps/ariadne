/***************************************************************************
 *            ciftypes.h
 *
 *  Copyright  2009  Davide Bresolin, Mirko Massignani
 *
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

/*! \file ciftypes.h
 *  \brief Datatypes and functions for CIF parsing.
 */

#ifndef ARIADNE_CIF_TYPES_H
#define ARIADNE_CIF_TYPES_H

#include <string>
#include <iostream>
#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include "hybrid_io_automaton.h"
#include "expression.h"
#include "formula.h"

using namespace Ariadne;

class CifModel;

typedef boost::shared_ptr<CifModel> CifModelRef;
typedef boost::shared_ptr<HybridIOAutomaton> HybridIOAutomatonRef;

enum CifModelType {
    MT_AUT,     // Single automaton
    MT_INST,    // Instantiation of another model
    MT_COMP,    // Parallel composition of models
    MT_UNDEF    // Type not defined yet
};

class CifModel
{
  private:
    std::string _name;
    CifModelType _type;
    HybridIOAutomatonRef _automaton;
    std::string _model;
    std::list< CifModelRef > _components;
  
  public:
    CifModel();
    CifModel(CifModelType type);
    CifModel(const std::string& name);
    
    std::string name() const;
    CifModelType type() const;
    bool is_automaton() const;
    bool is_instantiation() const;
    bool is_composition() const;
    bool is_undefined() const;
    HybridIOAutomaton& automaton();
    const HybridIOAutomaton& automaton() const; 
    const std::string& model() const;
    std::list< CifModelRef >& components();
    const std::list< CifModelRef >& components() const;

    void set_name(const std::string& name);
    void set_automaton(HybridIOAutomaton& automaton);
    void set_automaton(HybridIOAutomatonRef automaton);
    void set_model(const std::string& model);
    void add_component(CifModel& model);
    void add_component(CifModelRef model);
    void add_component(CifModel* model);
    void add_component(const std::string& name);
    void add_components(const std::list<CifModel*>& modlist);
};

//! \brief Generate a new name for an unnamed component
std::string generate_name();

//! \brief Generate a new event for a unlabelled transition
std::string generate_event();

//! \brief Merge two CifModel into one.
CifModel merge(CifModel& model1, CifModel& model2);

//! \brief Output a CifModel to a stream.
std::ostream& operator<<(std::ostream& os, const CifModel& model);

//! \brief Output a CifModelType to a stream.
std::ostream& operator<<(std::ostream& os, const CifModelType& type);

enum CifExprType{
    ET_BOOL,    // boolean predicate
    ET_REAL,    // real expression
    ET_TUPLE,   // tuple of expressions
    ET_FUNCID,  // function identifier
    ET_ASGN,    // variable assignment
    ET_VAR      // single real variable (dotted or non dotted).
};

struct CifExpression
{
    CifExprType type;
    
    RealVariable* realvar;              // Real variable
    bool dotted;                        // set to true when the variable is dotted
    ContinuousPredicate* boolexpr;      
    RealExpression* realexpr;
    std::list<CifExpression*>* tuple;
    Operator op;
};

struct CifIdentifier 
{
    int pos;
    int len;
    int line;
    char *name;
};

enum DeclType {
    DE_INVAR,       // Input variable (unspecified type)
    DE_DISC_OUTVAR, // Output discrete variable
    DE_CONT_OUTVAR, // Output continuous variable
    DE_DISC_VAR,    // Internal discrete variable
    DE_CONT_VAR,    // Internal continuous variable
    DE_CLOCK,       // Internal clock
    DE_INCHAN,      // Input channel
    DE_OUTCHAN,     // Output channel    
    DE_CONNECT,     // Connect statement
    DE_ACTION       // Discrete action label
};

struct Declaration
{
    DeclType type;
    char* id;
    char* id2;
};

//! \brief Convert an expression of type ET_VAR to ET_REAL. Do nothing if expr is of different type
void to_realexpr(CifExpression* exptr);

//! \brief Make a RealExpression from a function operator and a list of CifExpression
RealExpression make_expression(Operator op, const std::list<CifExpression*>& tuple);

//! \brief Output a CifExpression to a stream.
std::ostream& operator<<(std::ostream& os, const CifExpression& model);

//! \brief Output a CifExprType to a stream.
std::ostream& operator<<(std::ostream& os, const CifExprType& type);


#endif // ARIADNE_CIF_TYPES_H

