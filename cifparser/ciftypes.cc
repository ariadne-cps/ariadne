/***************************************************************************
 *            ciftypes.c
 *
 *  Copyright  2009  Davide Bresolin
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

#include "ciftypes.h"

//
// Static variables holding the number of the last generated name and event.
//
static int lastname = 0;
static int lastevent = 0;

//
// Methods of CifModel class
//
CifModel::CifModel()
    : _name(),
      _type(MT_UNDEF),
      _components()
{
}

CifModel::CifModel(CifModelType type)
    : _name(),
      _type(type),
      _components()
{
}

CifModel::CifModel(const std::string& name)
    : _name(name),
      _type(MT_UNDEF),
      _components()
{
}

std::string 
CifModel::name() const
{
    return this->_name;
}

CifModelType
CifModel::type() const
{
    return this->_type;
}

bool 
CifModel::is_automaton() const
{
    return this->_type == MT_AUT;
}

bool 
CifModel::is_instantiation() const
{
    return this->_type == MT_INST;
}

bool 
CifModel::is_composition() const
{
    return this->_type == MT_COMP;
}

bool 
CifModel::is_undefined() const
{
    return this->_type == MT_UNDEF;
}

HybridIOAutomaton& 
CifModel::automaton()
{
    ARIADNE_ASSERT_MSG(this->is_automaton(),
        "CifModel " << this->_name << " is not a single automaton.");
        
    return *(this->_automaton);        
}

const HybridIOAutomaton& 
CifModel::automaton() const
{
    ARIADNE_ASSERT_MSG(this->is_automaton(),
        "CifModel " << this->_name << " is not a single automaton.");
        
    return *(this->_automaton);        
}

const std::string& 
CifModel::model() const
{
    ARIADNE_ASSERT_MSG(this->is_instantiation(),
        "CifModel " << this->_name << " is not an instantiation of a model.");
        
    return this->_model;        
}

std::list<CifModelRef>& 
CifModel::components()
{
    ARIADNE_ASSERT_MSG(this->is_composition(),
        "CifModel " << this->_name << " is not a composition of models.");
        
    return this->_components;        
}

const std::list<CifModelRef>& 
CifModel::components() const
{
    ARIADNE_ASSERT_MSG(this->is_composition(),
        "CifModel " << this->_name << " is not a composition of models.");
        
    return this->_components;        
}
void 
CifModel::set_name(const std::string& name)
{
    this->_name = name;
    if(this->is_automaton())
        this->_automaton->set_name(name);
}

void 
CifModel::set_automaton(HybridIOAutomaton& automaton)
{
    HybridIOAutomatonRef autref(&automaton);
    this->set_automaton(autref);
}

void 
CifModel::set_automaton(HybridIOAutomatonRef automaton)
{
    ARIADNE_ASSERT_MSG(this->is_undefined(),
        "CifModel " << this->_name << " has type different from undefined.");

    this->_type = MT_AUT;
    this->_automaton = automaton;
    // change the automaton's name to be the same of the model
    if(this->_name == "") {
        this->_name = automaton->name();
    } else {
        automaton->set_name(this->_name);
    }
}


void 
CifModel::set_model(const std::string& model)
{
    ARIADNE_ASSERT_MSG(this->is_undefined(),
        "CifModel " << this->_name << " has type different from undefined.");

    this->_type = MT_INST;
    this->_model = model;    
}


void 
CifModel::add_component(CifModel& model)
{
    CifModelRef modref(&model);
    this->add_component(modref);
}

void 
CifModel::add_component(CifModelRef model)
{
    ARIADNE_ASSERT_MSG(this->is_undefined() || this->is_composition(),
        "CifModel " << this->_name << " is not a composition of models, cannot add a new component.");

    this->_type = MT_COMP;      
    this->_components.push_back(model);
}

void 
CifModel::add_component(CifModel* modptr)
{
    ARIADNE_ASSERT_MSG(this->is_undefined() || this->is_composition(),
        "CifModel " << this->_name << " is not a composition of models, cannot add a new component.");
    ARIADNE_ASSERT_MSG(modptr != NULL,
        "CifModel pointer is NULL: cannot add the component.");       
    this->_type = MT_COMP;      
    CifModelRef modref(modptr);
    this->_components.push_back(modref);
}

void 
CifModel::add_component(const std::string& name)
{
    ARIADNE_ASSERT_MSG(this->is_undefined() || this->is_composition(),
        "CifModel " << this->_name << " is not a composition of models, cannot add a new component.");

    this->_type = MT_COMP;        
    CifModelRef model(new CifModel(name));     
    this->_components.push_back(model);
}

void 
CifModel::add_components(const std::list<CifModel*>& modlist)
{
    for(std::list<CifModel*>::const_iterator iter=modlist.begin();
        iter != modlist.end(); iter++)
    {
        this->add_component(*iter);    
    }
}


CifModel merge(CifModel& model1, CifModel& model2) {
    ARIADNE_ASSERT_MSG(!model1.is_undefined(),
        "Model " << model1 << " is undefined, cannot merge with " << model2);
    ARIADNE_ASSERT_MSG(!model2.is_undefined(),
        "Model " << model2 << " is undefined, cannot merge with " << model1);

    CifModel result("");
   
    switch(model1.type()) 
    {
        // model1 is a single component: add it to result
        case MT_AUT:
        case MT_INST:   
            result.add_component(model1);
            break;
        // model1 is a composition: add all components to result
        case MT_COMP:
            for(std::list< CifModelRef >::const_iterator iter=model1.components().begin();
                iter != model1.components().end(); iter++)
            {
                result.add_component(*iter);
            }
            break;
        default:
            ARIADNE_FAIL_MSG("Model " << model1 << " is of unknown type " << model1.type());
    }
    // Finally, add model2 as a component of result
    result.add_component(model2);
    
    return result;
}

std::string generate_name()
{
    lastname++;
    std::stringstream out;
    out << "comp" << lastname;
    return out.str();
}

std::string generate_event()
{
    lastevent++;
    std::stringstream out;
    out << "dummy" << lastevent;
    return out.str();
}

void 
to_realexpr(CifExpression* exptr)
{
    if(exptr == NULL) {
        ARIADNE_FAIL_MSG("NULL expression pointer in to_realexpr.");
    }
    if(exptr->type == ET_VAR && !exptr->dotted) {
        exptr->type = ET_REAL;
        exptr->realexpr = new RealExpression(*exptr->realvar);
        delete exptr->realvar;
    }
}


RealExpression 
make_expression(Operator op, const std::list<CifExpression*>& tuple)
{
    RealExpression result;
    CifExpression* arg1;
    
    switch(op)
    {
        case SQRT:
            if(tuple.size() != 1) {
                throw std::runtime_error("Arity of function SQRT must be one.");                
            }
            arg1 = tuple.front();
            to_realexpr(arg1);
            if(arg1->type != ET_REAL) {
                throw std::runtime_error("Only real expression allowed inside function SQRT.");                
            }            
            result = sqrt(*arg1->realexpr);
            break;
        case EXP:   // Natural exponential
            if(tuple.size() != 1) {
                throw std::runtime_error("Arity of function EXP must be one.");               
            }
            arg1 = tuple.front();
            to_realexpr(arg1);
            if(arg1->type != ET_REAL) {
                throw std::runtime_error("Only real expression allowed inside function EXP.");                
            }            
            result = exp(*arg1->realexpr);
            break;
        case LOG:   // Natural logarithm
            if(tuple.size() != 1) {
                throw std::runtime_error("Arity of function LOG must be one.");                
            }
            arg1 = tuple.front();
            to_realexpr(arg1);
            if(arg1->type != ET_REAL) {
                throw std::runtime_error("Only real expression allowed inside function LOG.");                
            }            
            result = log(*arg1->realexpr);
            break;
        case SIN:   // Sine
            if(tuple.size() != 1) {
                throw std::runtime_error("Arity of function SIN must be one.");                
            }
            arg1 = tuple.front();
            to_realexpr(arg1);
            if(arg1->type != ET_REAL) {
                throw std::runtime_error("Only real expression allowed inside function SIN.");                
            }            
            result = sin(*arg1->realexpr);
            break;
        case COS:   // Cosine
            if(tuple.size() != 1) {
                throw std::runtime_error("Arity of function COS must be one.");                
            }
            arg1 = tuple.front();
            to_realexpr(arg1);
            if(arg1->type != ET_REAL) {
                throw std::runtime_error("Only real expression allowed inside function COS.");                
            }            
            result = cos(*arg1->realexpr);
            break;
        case TAN:   // Tangent
            if(tuple.size() != 1) {
                throw std::runtime_error("Arity of function TAN must be one.");                
            }
            arg1 = tuple.front();
            to_realexpr(arg1);
            if(arg1->type != ET_REAL) {
                throw std::runtime_error("Only real expression allowed inside function TAN.");                
            }            
            result = tan(*arg1->realexpr);
            break;
        case ABS:   // Absolute value
            if(tuple.size() != 1) {
                throw std::runtime_error("Arity of function ABS must be one.");               
            }
            arg1 = tuple.front();
            to_realexpr(arg1);
            if(arg1->type != ET_REAL) {
                throw std::runtime_error("Only real expression allowed inside function ABS.");                
            }            
            result = abs(*arg1->realexpr);
            break;
        default:
            throw std::runtime_error("Unsupported function.");            
    }
    return result;
}


std::ostream&
operator<<(std::ostream& os, const CifModel& mod)
{
    os << "CifModel( name=" << mod.name() << ", type=" << mod.type();
    
    switch(mod.type())
    {
        case MT_AUT:
            os << ", automaton=" << mod.automaton() << ")";
            break;
        case MT_INST:
            os << ", model=" << mod.model() << ")";
            break;
        case MT_COMP:
            os << ", components=[";
            for(std::list<CifModelRef>::const_iterator iter= mod.components().begin();
                iter != mod.components().end(); iter++)
            {
                os << std::endl << "\t" << **iter;
            }
            os << "] )";
            break;
        default:
            os << ", UNDEFINED )";        
    } 

    return os;    
 }

std::ostream&
operator<<(std::ostream& os, const CifModelType& type)
{
    switch(type)
    {
        case MT_AUT:
            os << "MT_AUT";
            break;
        case MT_INST:
            os << "MT_INST";
            break;
        case MT_COMP:
            os << "MT_COMP";
            break;
        default:
            os << "MT_UNDEF";
    }
    
    return os;
}

std::ostream&
operator<<(std::ostream& os, const CifExpression& expr)
{
    os << "Expression( type=" << expr.type;
    
    switch(expr.type)
    {
        case ET_BOOL:
            os << ", boolexpr=" << *expr.boolexpr << ")";
            break;
        case ET_REAL:
            os << ", realexpr=" << *expr.realexpr << ")";
            break;
        case ET_TUPLE:
            os << ", tuple=[";
            for(std::list<CifExpression*>::const_iterator iter= expr.tuple->begin();
                iter != expr.tuple->end(); iter++)
            {
                os << std::endl << "\t" << **iter;
            }
            os << "] )";
            break;
        case ET_FUNCID:
            os << ", function=" << expr.op << ")";
            break;
        case ET_ASGN:
            os << ", " << (expr.dotted ? "dotted " : "" ) 
               << "variable=" << *expr.realvar << ", expression=" << *expr.realexpr << ")";
            break;         
        case ET_VAR:
            os << ", " << (expr.dotted ? "dotted " : "" ) 
               << "variable=" << *expr.realvar << ")";
            break;         
        default:
            os << ", UNDEFINED )";        
    } 
    return os;    
 }

std::ostream&
operator<<(std::ostream& os, const CifExprType& type)
{
    switch(type)
    {
        case ET_BOOL:
            os << "ET_BOOL";
            break;
        case ET_REAL:
            os << "ET_REAL";
            break;
        case ET_TUPLE:
            os << "ET_TUPLE";
            break;
        case ET_FUNCID:
            os << "ET_FUNCID";
            break;
        case ET_ASGN:
            os << "ET_ASGN";
            break;         
        case ET_VAR:
            os << "ET_VAR";
            break;         
        default:
            os << "UNDEFINED";        
    } 
    
    return os;
}


