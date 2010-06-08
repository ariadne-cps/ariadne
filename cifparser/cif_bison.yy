/*******
CIF GRAMMAR
*******/

// prologue

%{
#define YYDEBUG 1

#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <math.h>
#include <string.h>
#include <assert.h>

#include <map>
#include <list>
#include <stack>

#include "ciftypes.h"

using namespace std;

int	yylex (void); 
extern int	yylineno; //is the number of line that we now parse
extern ostream* msgStream; //stream use to send tool information to user

void yyerror(char const* str){

	if(strcmp(str, "WARNING") == 0){
		(*msgStream) << "WARNING at line " << yylineno << ": " << str << std::endl;
		exit(0);
	}
	else{
		(*msgStream) << "ERROR at line " << yylineno << ": " << str << std::endl;
	}
}


#ifdef DEBUG_ENABLED
	#define DEBUG_HERE  printf ("[DBG %i]",__LINE__);
#else
	#define DEBUG_HERE ;
#endif

#define YYERROR_VERBOSE 1

#define CLEAR(_dest_) 			{memset(& _dest_ , 0, sizeof(_dest_) );}

#define TRANSFER(_dest_, _source_) 	{assert(sizeof(_source_)==sizeof(_dest_)); memcpy(& _dest_, & _source_, sizeof(_dest_));}

// Create the main CifModel
CifModel main_model;

// Create the map of model definitions
std::map< std::string, CifModelRef > moddef;

%}    

/*
 * Definition of the union YYSTYPE, needed by Bison
 */
%union{
	char character;
	char* value;

	CifIdentifier Identifier_data;

    std::list< CifIdentifier >* Identifier_list;

	struct{
		int pos;
		int len;
		int line;
		double value;
	} Constant_data;

	struct{
		int pos;
		int len;
		int line;
		char *name;
		Operator op;
	} FuncId_data;
	
    CifModel* Model_data;               // A single CIF component
    std::list<CifModel*>* Model_list;   // A parallel composition of components 
    HybridIOAutomaton* Automaton_data;  // An atomic automaton
    
    // List of modes and transitions
    struct{
        std::list<DiscreteIOMode*>* modes;
        std::list<DiscreteIOTransition*>* trans;
    } Modes_data;
    
    struct{
        std::list<DiscreteIOTransition*>* trans;
        std::map<RealVariable, RealExpression>* dynamics;
        std::list<RealExpression>* invariants;
    } ModeDecl_data;
    
    DiscreteIOTransition* Transition_data;
    
    CifExpression* Expression_data;
    
    std::list< CifExpression* >* ExprList_data;
    
    RealExpression* RealExpr_data;
    
    std::list< RealVariable >* VarList_data;
    
    std::map< RealVariable, RealExpression >* Reset_data; 
    
    Declaration Declaration_data;
    std::list< Declaration >* Declaration_list;
    
}

%type < Identifier_data > Y_IDENTIFIER Y_STRINGCONST autInst action chanAddressable comLabel
%type < Identifier_list > identifiers
%type < Constant_data > Y_NUMDOTNUM Y_NUMBERCONST Y_REALCONST
%type < FuncId_data > Y_FUNCIDENTIFIER
%type < Model_data > model closedScope automaton cAutomaton2 oAutomaton2 openScope
%type < Model_list > oAutomaton cAutomaton
%type < Automaton_data > atomicAut
%type < Modes_data > modes mode
%type < ModeDecl_data > modeDecl dyns dyn edges
%type < Transition_data > edge
%type < RealExpr_data > guard
%type < ExprList_data > preds exprs
%type < Expression_data > pred expr expr13 expr12 expr11 expr10 expr9 expr8 expr7 expr6 expr5 
                          expr4 expr3 expr2 expr0 constantExpr
%type < VarList_data > addressables addressable addressable1 addressable0     
%type < Reset_data > update
%type < Declaration_list > cScopeDecls cScopeDecl decls decl inputVarDecls varDecls connectSets 
                           chanDecls varDecl chanDecl
                    


%initial-action 
{ 
    // nothing
}
; 



%start spec
%token
Y_ABS
Y_ACOS
Y_ACOSH
Y_ASIN
Y_ASINH
Y_ATAN
Y_ATANH
Y_B2S
Y_BERNOULLI_B
Y_BERNOULLI_N
Y_BETA
Y_BINOMIAL
Y_CBRT
Y_CEIL
Y_CONSTANT
Y_COS
Y_COSH
Y_DQUOTE
Y_DROP
Y_ERLANG
Y_EXP
Y_EXPONENTIAL
Y_FLOOR
Y_GAMMA
Y_GEOMETRIC
Y_HD
Y_HR
Y_I2N
Y_I2R
Y_I2S
Y_INSERT
Y_LEN
Y_LN
Y_LOG
Y_LOGNORMAL
Y_N2I
Y_N2R
Y_N2S
Y_NL
Y_NONE
Y_NORMAL
Y_POISSON
Y_R2S
Y_RANDOM
Y_RESET
Y_RMOD
Y_ROUND
Y_S2B
Y_S2I
Y_S2N
Y_S2R
Y_SETSEED
Y_SIN
Y_SINH
Y_SIZE
Y_SORT
Y_SQRT
Y_TAB
Y_TAKE
Y_TAN
Y_TANH
Y_TL
Y_TR
Y_TRIANGLE
Y_UNIFORM_I
Y_UNIFORM_N
Y_UNIFORM_R
Y_WEIBULL

%token
Y_ACT
Y_ALG
Y_AND
Y_AUTOMATON
Y_BOOL
Y_CHAN
Y_CLOCK
Y_CONNECT
Y_CONT
Y_DISC
Y_DIV
Y_DO
Y_DOT
Y_EXTERN
Y_FALSE
Y_FLOW
Y_FUNC
Y_GOTO
Y_IN
Y_INIT
Y_INPUT
Y_INT
Y_INTERN
Y_INV
Y_MAX
Y_MIN
Y_MOD
Y_MODE
Y_MODEL
Y_NAT
Y_NOT
Y_NOW
Y_OLD
Y_OR
Y_OUTPUT
Y_PICK
Y_REAL
Y_SAMPLE
Y_STRING
Y_SUB
Y_TAU
Y_TCP
Y_TIME
Y_TRUE
Y_VAR
Y_VOID
Y_WHEN

%token
Y_IDENTIFIER 
Y_NUMBERCONST
Y_NUMDOTNUM  
Y_REALCONST  
Y_STRINGCONST
Y_FUNCIDENTIFIER 


%nonassoc Y_NE Y_LT Y_LE Y_GT Y_GE
%left Y_PLUS Y_MINUS Y_AMPERSAND

/*t_Space */
%token 
Y_AMPERSAND
Y_ARROW
Y_ASGN
Y_ASTERISK
Y_ATSIGN
Y_CARET
Y_COLON
Y_COLONCOLON
Y_COMMA
Y_DOTTOKEN
Y_DOUBLEARROW
Y_EQ
Y_GE
Y_GT
Y_INTERSECT
Y_LBRACK
Y_LCSCOPE
Y_LCUR
Y_LE
Y_LOSCOPE
Y_LSQ
Y_LT
Y_MINUS
Y_MINUSMINUS
Y_NE
Y_NUMBERSIGN
Y_PARA
Y_PIPE
Y_PLUS
Y_PLUSPLUS
Y_RBRACK
Y_RCSCOPE
Y_RCUR
Y_RECV
Y_ROSCOPE
Y_RSQ
Y_SEMICOLON
Y_SEND
Y_SLASH
Y_UNION

%%

spec            : nonModelDecls model nonModelDecls
{
                  main_model = *$2;
                  
                  (*msgStream) << "CIF model before instantiation = " << main_model << endl << endl;
                  
                  (*msgStream) << "Final list of automata definitions = " << moddef << endl;
}
;

annotations      : Y_ATSIGN annos Y_ATSIGN
                 | Y_ATSIGN Y_ATSIGN
                 |
;

annos            : anno
                 | annos Y_COMMA anno
;

anno             : type annoIden Y_EQ expr
;

annoIden         : Y_IDENTIFIER Y_COLON Y_IDENTIFIER
                 | annoIden Y_COLON Y_IDENTIFIER
;

nonModelDecls   : nonModelDecls funcDecl
                {
                    yyerror("Function declaration unsupported.");
                    YYERROR;
                }
                | nonModelDecls autDef
                |
;

funcDecl        : Y_FUNC Y_IDENTIFIER Y_LBRACK Y_RBRACK Y_ARROW type Y_EQ declarationBody
                | Y_FUNC Y_IDENTIFIER Y_LBRACK paramDecls Y_RBRACK Y_ARROW type Y_EQ declarationBody
;

declarationBody : Y_LCUR filenames Y_RCUR Y_COLONCOLON Y_STRINGCONST
;

filenames       : filenames Y_COMMA Y_STRINGCONST
                | Y_STRINGCONST
;

autDef          : Y_AUTOMATON Y_IDENTIFIER Y_LBRACK paramDecls Y_RBRACK Y_EQ closedScope
                {
                    yyerror("Automata declaration with parameters unsupported.");
                    YYERROR;
                }
                | Y_AUTOMATON Y_IDENTIFIER Y_LBRACK Y_RBRACK Y_EQ closedScope
                {
                    if(moddef.find($2.name) != moddef.end()) {  
                        // Automaton defined twice, error
                        std::stringstream err;
                        err << "Automaton " << $2.name << " already defined.";
                        yyerror(err.str().c_str());
                        YYERROR;
                    }
                    // Insert closedScope into model definitions.
                    $6->set_name($2.name);
                    CifModelRef modref($6);
                    moddef[$2.name] = modref;
                }   
;

paramDecls      : paramDecl
                | paramDecls Y_COMMA paramDecl
;

paramDecl       : identifiers Y_COLON type
;

identifiers     : Y_IDENTIFIER
                {
                    $$ = new std::list< CifIdentifier >;
                    $$->push_back($1);
                }
                | identifiers Y_COMMA Y_IDENTIFIER
                {
                    $$ = $1;
                    $$->push_back($3);
                }
;

closedScope     : Y_LCSCOPE cScopeDecls Y_COLONCOLON automaton Y_RCSCOPE
                {
                    $$ = $4;    // The definition of the model is given by automaton
                }
                | Y_LCSCOPE automaton Y_RCSCOPE
                {
                    $$ = $2;    // The definition of the model is given by automaton                
                }
;

cScopeDecls     : cScopeDecl
                {
                    $$ = $1;
                }
                | cScopeDecls cScopeDecl
                {
                    $$ = $1;
                    $$->insert($$->end(), $2->begin(), $2->end());
                }
;

cScopeDecl      : Y_EXTERN decls
                {
                    for(std::list< Declaration >::const_iterator iter = $2->begin();
                        iter != $2->end() ; iter++)
                    {
                        if(iter->type != DE_INCHAN && iter->type != DE_OUTCHAN)
                        {
                            yyerror("EXTERN declarations can contain only input or ouput channels.");
                            YYERROR;
                        }
                    }
                    $$ = $2;
                }
                | Y_INPUT Y_VAR inputVarDecls
                {
                    $$  = $3;
                }
                | Y_OUTPUT Y_VAR varDecls
                {
                    for(std::list< Declaration >::iterator iter = $3->begin();
                        iter != $3->end() ; iter++)
                    {
                        switch(iter->type)
                        {
                            case DE_DISC_VAR:
                                iter->type = DE_DISC_OUTVAR;
                                break;
                            case DE_CONT_VAR:
                                iter->type = DE_CONT_OUTVAR;
                                break;
                            default:
                                yyerror("Unsupported output var declaration.");
                                YYERROR;
                        }
                    }
                    $$  = $3;                    
                }
                | Y_INTERN decls
                {
                    $$ = $2;
                }
                | Y_CONNECT connectSets
                {
                    $$ = $2;
                }
;

decls           : decl
                {
                    $$ = $1;
                }
                | decls decl
                {
                    $$ = $1;
                    $$->insert($$->end(), $2->begin(), $2->end());
                }
;

decl            : Y_VAR varDecls
                {
                    $$ = $2;
                }
                | Y_CLOCK identifiers
                {
                    $$ = new std::list< Declaration >;
                    for(std::list< CifIdentifier >::const_iterator iter = $2->begin();
                        iter != $2->end() ; iter++)
                    {
                        Declaration decl;
                        decl.type = DE_CLOCK;
                        decl.id = iter->name;
                        $$->push_back(decl);
                    }
                }
                | Y_CHAN chanDecls
                {
                    $$ = $2;
                }
                | Y_ACT identifiers
                {
                    $$ = new std::list< Declaration >;
                    for(std::list< CifIdentifier >::const_iterator iter = $2->begin();
                        iter != $2->end() ; iter++)
                    {
                        Declaration decl;
                        decl.type = DE_ACTION;
                        decl.id = iter->name;
                        $$->push_back(decl);
                    }
                }
;

varDecls        : varDecl
                {
                    $$ = $1;
                }
                | varDecls Y_SEMICOLON varDecl
                {
                    $$ = $1;
                    $$->insert($$->end(), $3->begin(), $3->end());
                }
;

varDecl         : identifiers annotations Y_COLON Y_DISC type Y_EQ expr
                {
                    $$ = new std::list< Declaration >;
                    for(std::list< CifIdentifier >::const_iterator iter = $1->begin();
                        iter != $1->end(); iter++)
                    {
                        Declaration decl;
                        decl.type = DE_DISC_VAR;
                        decl.id = iter->name;
                        $$->push_back(decl);
                    }                    
                }
                | identifiers annotations Y_COLON Y_CONT type Y_EQ expr
                {
                    $$ = new std::list< Declaration >;
                    for(std::list< CifIdentifier >::const_iterator iter = $1->begin();
                        iter != $1->end(); iter++)
                    {
                        Declaration decl;
                        decl.type = DE_CONT_VAR;
                        decl.id = iter->name;
                        $$->push_back(decl);
                    }                    
                }
                | identifiers annotations Y_COLON Y_ALG  type Y_EQ expr
                {
                    yyerror("Algebraic variables not supported.");
                    YYERROR;
                }
                | identifiers annotations Y_COLON Y_DISC type
                {
                    $$ = new std::list< Declaration >;
                    for(std::list< CifIdentifier >::const_iterator iter = $1->begin();
                        iter != $1->end(); iter++)
                    {
                        Declaration decl;
                        decl.type = DE_DISC_VAR;
                        decl.id = iter->name;
                        $$->push_back(decl);
                    }                    
                }
                | identifiers annotations Y_COLON Y_CONT type
                {
                    $$ = new std::list< Declaration >;
                    for(std::list< CifIdentifier >::const_iterator iter = $1->begin();
                        iter != $1->end(); iter++)
                    {
                        Declaration decl;
                        decl.type = DE_CONT_VAR;
                        decl.id = iter->name;
                        $$->push_back(decl);
                    }                    
                }
                | identifiers annotations Y_COLON Y_ALG  type
                {
                    yyerror("Algebraic variables not supported.");
                    YYERROR;
                }
;

chanDecls       : chanDecl
                {
                    $$ = $1;
                }
                | chanDecls Y_SEMICOLON chanDecl
                {
                    $$ = $1;
                    $$->insert($$->end(), $3->begin(), $3->end());
                }
;

chanDecl        : identifiers Y_SEND Y_COLON chanType
                {
                    $$ = new std::list< Declaration >;
                    for(std::list< CifIdentifier >::const_iterator iter = $1->begin();
                        iter != $1->end(); iter++)
                    {
                        Declaration decl;
                        decl.type = DE_OUTCHAN;
                        decl.id = iter->name;
                        $$->push_back(decl);
                    }                    
                }
                | identifiers Y_RECV Y_COLON chanType
                {
                    $$ = new std::list< Declaration >;
                    for(std::list< CifIdentifier >::const_iterator iter = $1->begin();
                        iter != $1->end(); iter++)
                    {
                        Declaration decl;
                        decl.type = DE_INCHAN;
                        decl.id = iter->name;
                        $$->push_back(decl);
                    }                    
                }
                | identifiers Y_COLON chanType
                {
                    yyerror("Mixed input/output channels not supported.");
                    YYERROR;
                }
;

inputVarDecls   : inputVarDecl
                | inputVarDecls Y_SEMICOLON inputVarDecl
;

inputVarDecl    : identifiers annotations Y_COLON type
;

connectSets     : Y_LCUR connectors Y_RCUR
                {
                    $$ = new std::list<Declaration>;
                }
                | connectSets Y_COMMA Y_LCUR connectors Y_RCUR
                {
                    $$ = new std::list<Declaration>;
                }
;

connectors      : connector
                | connectors Y_COMMA connector
;

connector       : Y_IDENTIFIER Y_DOTTOKEN Y_IDENTIFIER
                | Y_IDENTIFIER
;

automaton       : cAutomaton    
                {
                    if($1->size() == 1) {    // List of a single element
                        $$ = $1->back();
                    } else {
                        $$ = new CifModel(generate_name());
                        $$->add_components(*$1);
                    }
                }
                | oAutomaton    
                {
                    if($1->size() == 1) {    // List of a single element
                        $$ = $1->back();
                    } else {
                        $$ = new CifModel(generate_name());
                        $$->add_components(*$1);
                    }
                }
;

cAutomaton      : cAutomaton2
                {
                    $$ = new std::list<CifModel*>();
                    $$->push_back($1);
                }
                | cAutomaton Y_PARA cAutomaton2
                {
                    $$ = $1;
                    $$->push_back($3);
                }
;

cAutomaton2     : Y_IDENTIFIER Y_COLON closedScope
                {
                    $$ = $3;
                    $$->set_name($1.name);                   
                }
                | Y_IDENTIFIER Y_COLON autInst
                {
                    $$ = new CifModel($1.name);
                    $$->set_model($3.name);
                }
                | autInst
                {
                    $$ = new CifModel(generate_name());
                    $$->set_model($1.name);
                }
                | closedScope
                {
                    $$ = $1;
                }
;

autInst         : Y_IDENTIFIER Y_LBRACK Y_RBRACK
                {
                    $$ = $1;
                }
                | Y_IDENTIFIER Y_LBRACK exprs Y_RBRACK
                {
                    yyerror("Instantiation of automata with parameters unsupported.");
                    YYERROR;
                }
;

oAutomaton      : oAutomaton2
                {
                    $$ = new std::list<CifModel*>();
                    $$->push_back($1);
                }
                | oAutomaton Y_PARA oAutomaton2
                {
                    $$ = $1;
                    $$->push_back($3);
                }
;

oAutomaton2     : atomicAut
                {
                    $$ = new CifModel(generate_name());
                    $$->set_automaton(*($1));
                }
                | openScope
                {
                    $$ = $1;
                }
                | Y_IDENTIFIER Y_COLON openScope
                {
                    $$ = $3;
                    $$->set_name($1.name);
                }
                | dyns
                {
                     yyerror("UNSUPPORTED: dynamics only.");
                     YYERROR;                    
                }
;

atomicAut       : Y_LOSCOPE init Y_COMMA Y_MODE modes Y_COLONCOLON Y_IDENTIFIER Y_ROSCOPE annotations
                {
                    $$ = new HybridIOAutomaton();
                    for(std::list<DiscreteIOMode*>::const_iterator iter=$5.modes->begin();
                        iter != $5.modes->end(); iter++)
                    {
                        $$->new_mode(**iter);
                    }
                    for(std::list<DiscreteIOTransition*>::const_iterator iter=$5.trans->begin();
                        iter != $5.trans->end(); iter++)
                    {
                        $$->new_transition(**iter);
                    }
                }
                | Y_LOSCOPE Y_MODE modes Y_COLONCOLON Y_IDENTIFIER Y_ROSCOPE annotations
                {
                    $$ = new HybridIOAutomaton();
                    for(std::list<DiscreteIOMode*>::const_iterator iter=$3.modes->begin();
                        iter != $3.modes->end(); iter++)
                    {
                        $$->new_mode(**iter);
                    }
                    for(std::list<DiscreteIOTransition*>::const_iterator iter=$3.trans->begin();
                        iter != $3.trans->end(); iter++)
                    {
                        $$->new_transition(**iter);
                    }
                }
;

init            : Y_INIT preds
;

preds           : pred
                {
                    $$ = new std::list<CifExpression*>;
                    $$->push_back($1);
                }
                | preds Y_AMPERSAND pred
                {
                    $$ = $1;
                    $$->push_back($3);
                }               
;

modes           : mode
                {
                    $$ = $1;
                }
                | modes Y_COMMA mode
                {
                    $$ = $1;
                    $$.modes->merge(*$3.modes);
                    $$.trans->merge(*$3.trans);
                }
;

mode            : Y_IDENTIFIER Y_EQ annotations modeDecl 
                {
                    DiscreteState loc($1.name);
                    DiscreteIOMode* mode = new DiscreteIOMode(loc, *($4.dynamics), *($4.invariants));
                    $$.modes = new std::list<DiscreteIOMode*>;                    
                    $$.modes->push_back(mode);
                    $$.trans = $4.trans;
                    // the source mode of the transitions is the current mode
                    for(std::list<DiscreteIOTransition*>::const_iterator iter = $$.trans->begin();
                        iter != $$.trans->end(); iter++)
                    {
                        (*iter)->set_source(loc);
                    }
                }
;

modeDecl        : dyns edges
                {
                    $$=$1;
                    $$.trans = $2.trans;
                }
                | dyns 
                {
                    $$ = $1;
                }
                | edges 
                {
                    $$ = $1;
                }
;

dyns            : dyn
                {
                    $$ = $1;
                }
                | dyns dyn
                {
                    $$ = $1;
                    $$.dynamics->insert($2.dynamics->begin(),$2.dynamics->end());
                    $$.invariants->insert($$.invariants->end(),$2.invariants->begin(),$2.invariants->end());
                }
;

dyn             : dyntype preds
                {
                    $$.trans = new std::list<DiscreteIOTransition*>;
                    $$.dynamics = new std::map<RealVariable, RealExpression>;
                    $$.invariants = new std::list<RealExpression>;
                    for(std::list<CifExpression*>::const_iterator iter = $2->begin();
                        iter != $2->end(); iter++)
                    {
                        switch((*iter)->type)
                        {
                            case ET_BOOL:  // expression is an invariant
                                $$.invariants->push_back(indicator(*(*iter)->boolexpr, negative));
                                break; 
                            case ET_ASGN:  // expression is an assignment
                                if((*iter)->dotted == false) {
                                    yyerror("Algebraic expressions for the dynamics are not supported.");
                                    YYERROR;
                                }
                                if($$.dynamics->find(*(*iter)->realvar) != $$.dynamics->end()) {
                                    yyerror("Dynamics defined twice for a variable.");
                                    YYERROR;
                                }
                                $$.dynamics->insert(make_pair(*(*iter)->realvar,*(*iter)->realexpr));
                                break;
                            default:
                                yyerror("Expression is not an invariant nor a dynamics.");
                                YYERROR;                           
                        }
                    }
                }
;

dyntype         : Y_FLOW
                | Y_INV
                | Y_TCP
;

edges           : edge
                {
                    $$.trans = new std::list<DiscreteIOTransition*>;
                    $$.dynamics = new std::map<RealVariable, RealExpression>;
                    $$.invariants = new std::list<RealExpression>;               
                    $$.trans->push_back($1);
                }
                | edges edge
                {
                    $$ = $1;
                    $$.trans->push_back($2);
                }
;   

edge            : Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(generate_event(), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($2.name), false);
                }
                | guard update Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(generate_event(), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($4.name), 
                                *$2,
                                *$1,
                                false);
                }
                | guard action Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(DiscreteEvent($2.name), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($4.name), 
                                *$1,
                                false);
                }
                | guard action update Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(DiscreteEvent($2.name), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($5.name), 
                                *$3,
                                *$1,
                                false);
                }
                | guard Y_NOW Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(generate_event(), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($4.name), 
                                *$1,
                                true);
                }
                | guard Y_NOW update Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(generate_event(), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($5.name), 
                                *$3,
                                *$1,
                                true);
                }
                | guard Y_NOW action Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(DiscreteEvent($3.name), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($5.name), 
                                *$1,
                                true);
                }
                | guard Y_NOW action update Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(DiscreteEvent($3.name), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($6.name), 
                                *$4,
                                *$1,
                                true);
                }
                | update Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(generate_event(), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($3.name), *$1, false);
                }
                | action Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(DiscreteEvent($1.name), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($3.name), false);
                }
                | action update Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(DiscreteEvent($1.name), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($4.name), *$2, false);
                }
                | Y_NOW Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(generate_event(), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($3.name), true);
                }
                | Y_NOW update Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(generate_event(), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($4.name), *$2, true);
                }
                | Y_NOW action Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(DiscreteEvent($2.name), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($4.name), true);
                }
                | Y_NOW action update Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(DiscreteEvent($2.name), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($5.name), *$3, true);
                }
                | guard Y_GOTO Y_IDENTIFIER
                {
                    $$ = new DiscreteIOTransition(generate_event(), 
                                DiscreteState(""),      // source mode will be specified later
                                DiscreteState($3.name), 
                                *$1,
                                false);
                }
;

guard           : Y_WHEN pred
                {
                    if($2->type != ET_BOOL) {
                        yyerror("Only boolean predicates are allowed in guards.");
                        YYERROR;
                    }
                    $$ = new RealExpression(indicator(*($2->boolexpr), positive));
                }
;

action          : Y_ACT Y_IDENTIFIER
                {
                    $$ = $2;                    
                }
                | Y_ACT Y_TAU
                {
                    std::string event=generate_event();
                    $$.name = new char[event.size()+1];
                    strcpy($$.name,event.c_str());
                }
                | Y_ACT comLabel
                {
                    $$ = $2;
                }
;

comLabel        : chanAddressable Y_SEND exprs
                {
                    yyerror("Data exchange through channels not supported.");
                    YYERROR;
                }
                | chanAddressable Y_SEND
                {
                    $$ = $1;
                }
                | chanAddressable Y_RECV addressables
                {
                    yyerror("Data exchange through channels not supported.");
                    YYERROR;
                }
                | chanAddressable Y_RECV
                {
                    $$ = $1;
                }
                | chanAddressable Y_SEND Y_RECV addressables Y_ASGN exprs
                {
                    yyerror("Bidirectional communication through channels not supported.");
                    YYERROR;
                }
                | chanAddressable Y_RECV Y_SEND addressables Y_ASGN exprs
                {
                    yyerror("Bidirectional communication through channels not supported.");
                    YYERROR;
                }
                | chanAddressable Y_SEND Y_RECV
                {
                    yyerror("Bidirectional communication through channels not supported.");
                    YYERROR;
                }
                | chanAddressable Y_RECV Y_SEND
                {
                    yyerror("Bidirectional communication through channels not supported.");
                    YYERROR;
                }
;

update          : Y_DO addressables Y_ASGN exprs
                {
                    std::list< CifExpression* >* elptr;
                    if($4->size() == 1 && $4->front()->type == ET_TUPLE)
                    {
                        // exprs is a sigle tuple of expressions, converto to a list
                        elptr = $4->front()->tuple;
                    } else {
                        // exprs is a list of expressions, do nothing
                        elptr = $4;
                    }
                    
                    if($2->size() != elptr->size()) 
                    {
                        yyerror("Syntax error in reset.");
                        YYERROR;
                    }
                    
                    $$ = new std::map< RealVariable, RealExpression >;
                    std::list< CifExpression* >::const_iterator eiter = elptr->begin();
                    std::list< RealVariable >::const_iterator viter = $2->begin();
                    while(eiter != elptr->end())
                    {
                        to_realexpr(*eiter);
                        if((*eiter)->type != ET_REAL)
                        {
                            yyerror("Only real expression are allowed in resets.");
                            YYERROR;                           
                        }
                        if($$->find(*viter) != $$->end())
                        {                            
                            yyerror("Reset defined twice for a variable.");
                            YYERROR;                           
                        }
                        (*$$)[*viter] = *((*eiter)->realexpr);
                        viter++;
                        eiter++;
                    }        
                }
                | Y_DO jumpVarsSet Y_COLON preds
                {
                    yyerror("Reset predicates with Jump var sets not supported.");
                    YYERROR;
                }
;

addressables    : addressable
                {
                    $$ = $1;
                }
                | addressables Y_COMMA addressable
                {
                    if($3->size() != 1) {
                        yyerror("Syntax error in assignment.");
                        YYERROR;
                    }                  
                    $$ = $1;
                    $$->push_back($3->front());                   
                }
;

addressable     : addressable1
                {
                    $$ = $1;
                }
;

addressable1    : addressable1 Y_DOTTOKEN expr2
                {
                   yyerror("Operator . not allowed in assignments.");
                   YYERROR;                    
                }
                | addressable1 Y_DOTTOKEN Y_NUMDOTNUM
                {
                   yyerror("Operator . not allowed in assignments.");
                   YYERROR;                    
                }
                | addressable0
                {
                   $$ = $1;
                }
;

addressable0    : Y_IDENTIFIER
                {
                    $$ = new std::list<RealVariable>;
                    $$->push_back(RealVariable($1.name));
                }
                | Y_LT addressables Y_GT
                {
                    $$ = $2;
                }
                | Y_LBRACK addressables Y_RBRACK
                {
                    $$ = $2;
                }
;

jumpVarsSet     : Y_LCUR Y_RCUR
                | Y_LCUR jumpVars Y_RCUR
;

jumpVars        : addressable
                | jumpVars Y_COMMA addressable
;

projections     : Y_DOTTOKEN Y_NUMDOTNUM
                | Y_DOTTOKEN expr2
                | projections Y_DOTTOKEN Y_NUMDOTNUM
                | projections Y_DOTTOKEN expr2
;

chanAddressable : Y_IDENTIFIER
                {
                    $$ = $1;
                }
                | Y_IDENTIFIER projections
                {
                    yyerror("Projection operators not supported.");
                    YYERROR;
                }
;

openScope       : Y_LOSCOPE oScopeDecls Y_COLONCOLON oAutomaton Y_ROSCOPE
                {
                    if($4->size() == 1) {    // List of a single element
                        $$ = $4->back();
                    } else {
                        $$ = new CifModel(generate_name());
                        $$->add_components(*$4);
                    }
                }
                | Y_LOSCOPE oAutomaton Y_ROSCOPE
                {
                    if($2->size() == 1) {    // List of a single element
                        $$ = $2->back();
                    } else {
                        $$ = new CifModel(generate_name());
                        $$->add_components(*$2);
                    }
                }
;

oScopeDecls     : oScopeDecl
                | oScopeDecls oScopeDecl
;

oScopeDecl      : Y_INTERN decls
;

model           : Y_MODEL Y_IDENTIFIER Y_LBRACK Y_RBRACK Y_EQ annotations closedScope
                {                   
                    $$ = $7;    // The model is defined by the closedScope
                    string modname = $2.name;
                    $$->set_name(modname);  // Set the model name according to Y_IDENTIFIER
                    (*msgStream) << "Main model name " << modname << std::endl;
                }
;

types           : type
                | types Y_COMMA type
;

type            : basicType
                | containerType
                {
                    yyerror("Container types not supported.");
                    YYERROR;
                }                
                | functionType
                {
                    yyerror("Function types not supported.");
                    YYERROR;
                }                
                | distributionType
                {
                    yyerror("Distribution types not supported.");
                    YYERROR;
                }                
                | Y_LBRACK types Y_RBRACK
                {
                    yyerror("Multiple types not supported.");
                    YYERROR;
                }                
;

basicType       : Y_BOOL
                {
                    yyerror("Type BOOL not supported.");
                    YYERROR;
                }
                | Y_NAT
                | Y_INT
                | Y_REAL
                | Y_STRING
                {
                    yyerror("Type STRING not supported.");
                    YYERROR;
                }
;

containerType   : Y_LSQ type Y_RSQ
                | constantExpr Y_ASTERISK type
                | Y_LCUR type Y_RCUR
;

functionType    : Y_LBRACK Y_RBRACK Y_ARROW type
                | Y_LBRACK types Y_RBRACK Y_ARROW type
;

distributionType : Y_ARROW type
;

chanType        : Y_VOID
                | type
                {
                    yyerror("Channels with type different from VOID not supported.");
                    YYERROR;
                }
                | bundleSizes Y_VOID
                {
                    yyerror("Channels with type different from VOID not supported.");
                    YYERROR;
                }
                | bundleSizes type
                {
                    yyerror("Channels with type different from VOID not supported.");
                    YYERROR;
                }
;

bundleSizes     : Y_NUMBERCONST Y_NUMBERSIGN
                | bundleSizes Y_NUMBERCONST Y_NUMBERSIGN
;

pred            : expr
                {
                    $$ = $1;
                }
;


exprs           : expr
                {
                    $$ = new std::list<CifExpression*>;
                    $$->push_back($1);
                }
                | exprs Y_COMMA expr
                {
                    $$ = $1;
                    $$->push_back($3);
                }
;

expr            : expr Y_OR expr13
                {
                    if($1->type != ET_BOOL || $3->type != ET_BOOL)
                    {
                        yyerror("OR operator allowed only in boolean predicates.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(*$1->boolexpr || *$3->boolexpr);
                }
                | expr13
                {
                    $$ = $1;
                }
;

expr13          : expr13 Y_AND expr12
                {
                    if($1->type != ET_BOOL || $3->type != ET_BOOL)
                    {
                        yyerror("AND operator allowed only in boolean predicates.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(*$1->boolexpr && *$3->boolexpr);
                }                    
                | expr12
                {
                    $$ = $1;
                }                
;

expr12          : expr12 Y_DOUBLEARROW expr11
                {
                    if($1->type != ET_BOOL || $3->type != ET_BOOL)
                    {
                        yyerror("=> operator allowed only in boolean predicates.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(!(*$1->boolexpr) || *$3->boolexpr);
                }                    
                | expr11
;

expr11          : expr11 Y_LT expr10
                {
                    to_realexpr($1);
                    to_realexpr($3);
                    if($1->type != ET_REAL || $3->type != ET_REAL)
                    {
                        yyerror("Comparison operator allowed only between real predicates.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(*$1->realexpr < *$3->realexpr);                    
                }                
                | expr11 Y_LE expr10
                {
                    to_realexpr($1);
                    to_realexpr($3);
                    if($1->type != ET_REAL || $3->type != ET_REAL)
                    {
                        yyerror("Comparison operator allowed only between real predicates.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(*$1->realexpr <= *$3->realexpr);                    
                }                
                | expr11 Y_EQ expr10
                {
                    to_realexpr($3);
                    if($1->type != ET_VAR || $3->type != ET_REAL)
                    {
                        (*msgStream) << *$1 << " = " << *$3 << endl;
                        yyerror("Equality operator allowed only in assignments.");
                        YYERROR;
                    }
                    $$ = $3;
                    $$->type = ET_ASGN;
                    // extract the RealVariable from the expression:
                    // it is the first (and only) element of arguments()
                    $$->realvar = $1->realvar;
                    $$->dotted = $1->dotted;
                }                
                | expr11 Y_NE expr10
                {
                    yyerror("Inequality operator not allowed.");
                    YYERROR;
                }                
                | expr11 Y_GT expr10
                {
                    to_realexpr($1);
                    to_realexpr($3);
                    if($1->type != ET_REAL || $3->type != ET_REAL)
                    {
                        yyerror("Comparison operator allowed only between real predicates.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(*$1->realexpr > *$3->realexpr);                    
                }                              
                | expr11 Y_GE expr10
                {
                    to_realexpr($1);
                    to_realexpr($3);
                    if($1->type != ET_REAL || $3->type != ET_REAL)
                    {
                        yyerror("Comparison operator allowed only between real predicates.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(*$1->realexpr >= *$3->realexpr);                    
                }                
                | expr11 Y_IN expr10
                {   
                    yyerror("IN operator not supported.");
                    YYERROR;
                }
                | expr11 Y_SUB expr10
                {   
                    yyerror("SUB operator not supported.");
                    YYERROR;
                }
                | expr10
                {
                    $$ = $1;
                }
;

vectorExprs     : vectorExpr
                | vectorExprs Y_COMMA vectorExpr
;

vectorExpr      : vectorExpr Y_OR vectorExpr13
                | vectorExpr13
;

vectorExpr13    : vectorExpr13 Y_AND vectorExpr12
                | vectorExpr12
;

vectorExpr12    : vectorExpr12 Y_DOUBLEARROW vectorExpr11
                | vectorExpr11
;

vectorExpr11    : vectorExpr11 Y_LT expr10
                | vectorExpr11 Y_LE expr10
                | vectorExpr11 Y_EQ expr10
                | vectorExpr11 Y_NE expr10
                | vectorExpr11 Y_GE expr10
                | vectorExpr11 Y_IN expr10
                | vectorExpr11 Y_SUB expr10
                | expr10
;

expr10          : Y_NOT expr10
                {                     
                    if($2->type != ET_BOOL)
                    {
                        yyerror("NOT operator allowed only in boolean predicates.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(!*$2->boolexpr);                    
                }                
                | expr9
                {
                    $$ = $1;
                }
;

expr9           : expr9 Y_MIN expr8
                {   
                    yyerror("MIN operator not supported.");
                    YYERROR;
                }
                | expr9 Y_MAX expr8
                {   
                    yyerror("MAX operator not supported.");
                    YYERROR;
                }
                | expr8
                {   
                    $$ = $1;
                }
;

expr8           : expr8 Y_PLUS expr7
                {
                    to_realexpr($1);
                    to_realexpr($3);
                    if($1->type != ET_REAL)
                    {
                        yyerror("Operator + allowed only in real expressions.");
                        YYERROR;
                    }
                    if($3->type != ET_REAL)
                    {
                        yyerror("Operator + allowed only in real expressions.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_REAL;
                    $$->realexpr = new RealExpression(*$1->realexpr + *$3->realexpr);                    
                }                
                | expr8 Y_UNION expr7
                {   
                    yyerror("UNION operator not supported.");
                    YYERROR;
                }
                | expr8 Y_MINUS expr7
                {
                    to_realexpr($1);
                    to_realexpr($3);
                    if($1->type != ET_REAL)
                    {
                        yyerror("Operator - allowed only in real expressions.");
                        YYERROR;
                    }
                    if($3->type != ET_REAL)
                    {
                        yyerror("Operator - allowed only in real expressions.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_REAL;
                    $$->realexpr = new RealExpression(*$1->realexpr - *$3->realexpr);                    
                }
                | expr8 Y_PLUSPLUS expr7
                {   
                    yyerror("Operator ++ not supported.");
                    YYERROR;
                }
                | expr8 Y_MINUSMINUS expr7
                {   
                    yyerror("Operator -- not supported.");
                    YYERROR;
                }
                | expr7
                {   
                    $$ = $1;
                }
;

expr7           : expr7 Y_ASTERISK expr6
                {
                    to_realexpr($1);
                    to_realexpr($3);
                    if($1->type != ET_REAL)
                    {
                        yyerror("Operator * allowed only in real expressions.");
                        YYERROR;
                    }
                    if($3->type != ET_REAL)
                    {
                        yyerror("Operator * allowed only in real expressions.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_REAL;
                    $$->realexpr = new RealExpression(*$1->realexpr * *$3->realexpr);                    
                }                
                | expr7 Y_INTERSECT expr6
                {   
                    yyerror("INTERSECT operator not supported.");
                    YYERROR;
                }
                | expr7 Y_SLASH expr6
                {
                    to_realexpr($1);
                    to_realexpr($3);
                    if($1->type != ET_REAL)
                    {
                        yyerror("Operator / allowed only in real expressions.");
                        YYERROR;
                    }
                    if($3->type != ET_REAL)
                    {
                        yyerror("Operator / allowed only in real expressions.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_REAL;
                    $$->realexpr = new RealExpression(*$1->realexpr / *$3->realexpr);                    
                }                
                | expr7 Y_DIV expr6
                {   
                    yyerror("DIV operator not supported.");
                    YYERROR;
                }
                | expr7 Y_MOD expr6
                {   
                    yyerror("MOD operator not supported.");
                    YYERROR;
                }
                | expr6
                {   
                    $$ = $1;
                }
;

expr6           : expr5 Y_CARET expr6
                {   
                    yyerror("operator ^ not supported.");
                    YYERROR;
                }
                | expr5
                {   
                    $$ = $1;
                }
;

expr5           : Y_SAMPLE expr5
                {   
                    yyerror("SAMPLE operator not supported.");
                    YYERROR;
                }
                | Y_PICK expr5
                {   
                    yyerror("PICK operator not supported.");
                    YYERROR;
                }
                | Y_PLUS expr5
                {
                    to_realexpr($2);
                    if($2->type != ET_REAL)
                    {
                        yyerror("Operator + allowed only in real expressions.");
                        YYERROR;
                    }
                    $$ = $2;                    
                }
                | Y_MINUS expr5
                {
                    to_realexpr($2);
                    if($2->type != ET_REAL)
                    {
                        yyerror("Operator - allowed only in real expressions.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_REAL;
                    $$->realexpr = new RealExpression(- *$2->realexpr);                    
                }
                | expr4
                {   
                    $$ = $1;
                }
;

expr4           : Y_DOT expr4
                {
                    if($2->type != ET_VAR) {
                        yyerror("DOT operator can be used only in front of assignments.");
                        YYERROR;
                    }
                    if($2->dotted) {
                        yyerror("Multiple occurrences of DOT are not allowed.");
                        YYERROR;
                    }                    
                    $$ = $2;
                    $$->dotted = true;
                }
                | expr3
                {   
                    $$ = $1;
                }
;

expr3           : expr3 Y_DOTTOKEN expr2
                {   
                    yyerror("operator . not supported.");
                    YYERROR;
                }
                | expr3 Y_DOTTOKEN Y_NUMDOTNUM
                {   
                    yyerror("operator . not supported.");
                    YYERROR;
                }
                | Y_NUMDOTNUM
                {
                    $$ = new CifExpression;
                    $$->type = ET_REAL;
                    $$->realexpr = new RealExpression($1.value);
                }
                | expr2
                {   
                    $$ = $1;
                }
;

expr2           : expr2 Y_LBRACK Y_RBRACK
                {
                    if($1->type != ET_FUNCID)
                    {   
                        yyerror("Incorrect use of parenthesis.");
                        YYERROR;
                    }
                    std::list< CifExpression*> empty;
                    $$ = new CifExpression;
                    $$->type = ET_REAL;               
                    try {
                        $$->realexpr = new RealExpression(make_expression($1->op, empty));
                    } catch(const std::runtime_error& e) {
                        yyerror(e.what());
                        YYERROR;
                    }                    
                }
                | expr2 Y_LBRACK exprs Y_RBRACK
                {
                    if($1->type != ET_FUNCID)
                    {   
                        yyerror("Incorrect use of parenthesis.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_REAL;                    
                    try {
                        $$->realexpr = new RealExpression(make_expression($1->op, *$3));
                    } catch(const std::runtime_error& e) {
                        yyerror(e.what());
                        YYERROR;
                    }
                }
                | expr0
                {   
                    $$ = $1;
                }
;

expr0           : constantExpr
                {
                    $$ = $1;
                }
                | Y_LT vectorExprs Y_GT
                {   
                    yyerror("operator < > not supported.");
                    YYERROR;
                }
                | Y_LBRACK condAlts Y_RBRACK
                {
                    yyerror("Conditional alternatives operator | not supported.");
                    YYERROR;
                }
                | Y_IDENTIFIER
                {
                    $$ = new CifExpression;
                    $$->type = ET_VAR;
                    $$->realvar = new RealVariable($1.name);
                    $$->dotted = false;
                }
                | Y_TIME
                {   
                    yyerror("TIME operator not supported.");
                    YYERROR;
                }
                | Y_FUNCIDENTIFIER
                {
                    if($1.op == CNST)
                    {
                        yyerror("Function not supported.");
                        YYERROR;
                    }
                    $$ = new CifExpression;
                    $$->type = ET_FUNCID;
                    $$->op = $1.op;
                }                
                | Y_OLD Y_LBRACK expr Y_RBRACK
                {   
                    yyerror("operator OLD not supported.");
                    YYERROR;
                }
                | Y_LBRACK exprs Y_RBRACK
                {                    
                    if($2->size() != 1) 
                    {
                        $$ = new CifExpression;
                        $$->type = ET_TUPLE;
                        $$->tuple = $2;
                    } else {
                        $$ = $2->front();
                    }
                }                                                
                | Y_LSQ Y_RSQ
                {   
                    yyerror("operator [ ] not supported.");
                    YYERROR;
                }
                | Y_LSQ exprs Y_RSQ
                {   
                    yyerror("operator [ ] not supported.");
                    YYERROR;
                }
                | Y_LCUR Y_RCUR
                {   
                    yyerror("operator { } not supported.");
                    YYERROR;
                }
                | Y_LCUR exprs Y_RCUR
                {   
                    yyerror("operator { } not supported.");
                    YYERROR;
                }
;

constantExpr    : Y_NUMBERCONST
                {
                    $$ = new CifExpression;
                    $$->type = ET_REAL;
                    $$->realexpr = new RealExpression($1.value);
                }
                | Y_REALCONST
                {
                    $$ = new CifExpression;
                    $$->type = ET_REAL;
                    $$->realexpr = new RealExpression($1.value);
                }
                | Y_STRINGCONST
                {   
                    yyerror("String constants not supported.");
                    YYERROR;
                }
                | Y_TRUE
                {
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(true);
                }
                | Y_FALSE
                {
                    $$ = new CifExpression;
                    $$->type = ET_BOOL;
                    $$->boolexpr = new ContinuousPredicate(false);
                }
;

condAlts        : condAlt
                | condAlts Y_PIPE condAlt
;

condAlt         : expr Y_ARROW expr
;

%%
