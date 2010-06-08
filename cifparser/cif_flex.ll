%option   yylineno

%{
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <search.h>
#include "ciftypes.h"
#include "cif_bison.h"

using namespace std;

extern char *yysptr;
extern char *yysbuf;
extern int yylineno;
extern int yycolumno;
extern char yytchar;

#define YYLMAX 1024

#define MVL_LG_MC 15
#define MVL_NB_MC 47
#define MVL_FN_MC 67

typedef struct {
        char id[MVL_LG_MC];
        int kval;
        } el_mc;

typedef struct {
        char id[MVL_LG_MC];
        int kval;
        Operator op;
        } fn_id;

/**
  struct of reserved function of language, 
  WARNING: when add an entry here must be increment MVL_FN_MC value
 */
static fn_id tab_fn []=
  {
	{"abs"		    , Y_ABS         , ABS},
	{"acos"		    , Y_ACOS        , CNST},
	{"acosh"	    , Y_ACOSH       , CNST},
	{"asin"		    , Y_ASIN        , CNST},
	{"asinh"	    , Y_ASINH       , CNST},
	{"atan"		    , Y_ATAN        , CNST},
	{"atanh"	    , Y_ATANH       , CNST},
	{"b2s"		    , Y_B2S         , CNST},
	{"bernoulli_b"	, Y_BERNOULLI_B , CNST},
	{"bernoulli_n"	, Y_BERNOULLI_N , CNST},
	{"beta"		    , Y_BETA        , CNST},
	{"binomial"	    , Y_BINOMIAL    , CNST},
	{"cbrt"		    , Y_CBRT        , CNST},
	{"ceil"		    , Y_CEIL        , CNST},
	{"constant"	    , Y_CONSTANT    , CNST},
	{"cos"		    , Y_COS         , COS},
	{"cosh"		    , Y_COSH        , CNST},
	{"dquote"	    , Y_DQUOTE      , CNST},
	{"drop"		    , Y_DROP        , CNST},
	{"erlang"	    , Y_ERLANG      , CNST},
	{"exp"		    , Y_EXP         , EXP},
	{"exponential"	, Y_EXPONENTIAL , CNST},
	{"floor"	    , Y_FLOOR       , CNST},
	{"gamma"	    , Y_GAMMA       , CNST},
	{"geometric"	, Y_GEOMETRIC   , CNST},
	{"hd"		    , Y_HD          , CNST},
	{"hr"		    , Y_HR          , CNST},
	{"i2n"		    , Y_I2N         , CNST},
	{"i2r"		    , Y_I2R         , CNST},
	{"i2s"		    , Y_I2S         , CNST},
	{"insert"	    , Y_INSERT      , CNST},
	{"len"		    , Y_LEN         , CNST},
	{"ln"		    , Y_LN          , CNST},
	{"log"		    , Y_LOG         , LOG},
	{"lognormal"	, Y_LOGNORMAL   , CNST},
	{"n2i"		    , Y_N2I         , CNST},
	{"n2r"		    , Y_N2R         , CNST},
	{"n2s"		    , Y_N2S         , CNST},
	{"nl"		    , Y_NL          , CNST},
	{"none"		    , Y_NONE        , CNST},
	{"normal"	    , Y_NORMAL      , CNST},
	{"poisson"	    , Y_POISSON     , CNST},
	{"r2s"		    , Y_R2S         , CNST},
	{"random"	    , Y_RANDOM      , CNST},
	{"reset"	    , Y_RESET       , CNST},
	{"rmod"		    , Y_RMOD        , CNST},
	{"round"	    , Y_ROUND       , CNST},
	{"s2b"		    , Y_S2B         , CNST},
	{"s2i"		    , Y_S2I         , CNST},
	{"s2r"		    , Y_S2R         , CNST},
	{"setseed"	    , Y_SETSEED     , CNST},
	{"sin"		    , Y_SIN         , SIN},
	{"sinh"		    , Y_SINH        , CNST},
	{"size"		    , Y_SIZE        , CNST},
	{"sort"		    , Y_SORT        , CNST},
	{"sqrt"		    , Y_SQRT        , SQRT},
	{"tab"		    , Y_TAB         , CNST},
	{"take"		    , Y_TAKE        , CNST},
	{"tan"		    , Y_TAN         , TAN},
	{"tanh"		    , Y_TANH        , CNST},
	{"tl"		    , Y_TL          , CNST},
	{"tr"		    , Y_TR          , CNST},
	{"triangle"	    , Y_TRIANGLE    , CNST},
	{"uniform_i"	, Y_UNIFORM_I   , CNST},
	{"uniform_n"	, Y_UNIFORM_N   , CNST},
	{"uniform_r"	, Y_UNIFORM_R   , CNST},
	{"weibull"	    , Y_WEIBULL     , CNST}
  };

/**
  struct of reserved operators of language, 
  WARNING: when add an entry here must be increment MVL_NB_MC value
 */
static el_mc tab_mc []=
  {
	{"act"		, Y_ACT},
	{"alg"		, Y_ALG},
	{"and"		, Y_AND},
	{"automaton"	, Y_AUTOMATON},
	{"bool"		, Y_BOOL},
	{"chan"		, Y_CHAN},
	{"clock"	, Y_CLOCK},
	{"connect"	, Y_CONNECT},
	{"cont"		, Y_CONT},
	{"disc"		, Y_DISC},
	{"div"		, Y_DIV},
	{"do"		, Y_DO},
	{"dot"		, Y_DOT},
	{"extern"	, Y_EXTERN},
	{"false"	, Y_FALSE},
	{"flow"		, Y_FLOW},
	{"func"		, Y_FUNC},
	{"goto"		, Y_GOTO},
	{"in"		, Y_IN},
	{"init"		, Y_INIT},
	{"input"	, Y_INPUT},
	{"int"		, Y_INT},
	{"intern"	, Y_INTERN},
	{"inv"		, Y_INV},
	{"max"		, Y_MAX},
	{"min"		, Y_MIN},
	{"mod"		, Y_MOD},
	{"mode"		, Y_MODE},
	{"model"	, Y_MODEL},
	{"nat"		, Y_NAT},
	{"not"		, Y_NOT},
	{"now"		, Y_NOW},
	{"old"		, Y_OLD},
	{"or"		, Y_OR},
	{"output"	, Y_OUTPUT},
	{"pick"		, Y_PICK},
	{"real"		, Y_REAL},
	{"sample"	, Y_SAMPLE},
	{"string"	, Y_STRING},
	{"sub"		, Y_SUB},
	{"tau"		, Y_TAU},
	{"tcp"		, Y_TCP},
	{"time"		, Y_TIME},
	{"true"		, Y_TRUE},
	{"var"		, Y_VAR},
	{"void"		, Y_VOID},
	{"when"		, Y_WHEN}
  };


typedef int compare_function_type(const void *, const void *);

static int find_mc(char * s)
{
  char  loc[YYLMAX];
  int   l;
  el_mc *pt;
  l=strlen(s);
  strcpy(loc,s);
  while(l--) loc[l]=tolower(loc[l]);	// conversion to lower_case 

  pt= (el_mc *) bsearch(loc, (void *)tab_mc, MVL_NB_MC, sizeof(el_mc), (compare_function_type *)strcmp);
  if (pt==NULL) return(-1);
  else return(pt->kval);
}

static int find_fn(char * s, Operator& op)
{
  char  loc[YYLMAX];
  int   l;
  fn_id *fn;
  l=strlen(s);
  strcpy(loc,s);
  while(l--) loc[l]=tolower(loc[l]);	// conversion to lower_case 

  fn= (fn_id *) bsearch(loc, (void *)tab_fn, MVL_FN_MC, sizeof(fn_id), (compare_function_type *)strcmp);
  if (fn==NULL) return(-1);
  op = fn->op;
  return (fn->kval);
}
%}

space_character 	  		[ \t]
digit		  			[0-9]
upper_case_letter 	  		[A-Z]
lower_case_letter 	  		[a-z]
letter	   	  			({upper_case_letter}|{lower_case_letter})
letter_or_digit	   	  		({letter}|{digit})

Y_IDENTIFIER 	     		({letter})\w*     	
Y_NUMBERCONST        		{digit}+
Y_NUMDOTNUM                     {digit}+\.{digit}+
Y_REALCONST          		{digit}+(\.{digit})?[Ee]?[\-\+]?{digit}+
Y_STRINGCONST        		\"[^"]*\"

%%
{space_character}	{ yycolumno++;   }
"&"       { yycolumno++;  return(Y_AMPERSAND); }
"->"      { yycolumno+=2;  return(Y_ARROW); }
":="      { yycolumno+=2;  return(Y_ASGN); }
"*"       { yycolumno++;  return(Y_ASTERISK); }
"@"	  { yycolumno++;  return(Y_ATSIGN); }
"^"       { yycolumno++;  return(Y_CARET); }
":"       { yycolumno++;  return(Y_COLON); }
"::"      { yycolumno+=2;  return(Y_COLONCOLON); }
","       { yycolumno++;  return(Y_COMMA); }
"."       { yycolumno++;  return(Y_DOTTOKEN); }
"=>"      { yycolumno+=2;  return(Y_DOUBLEARROW); }
"="	  { yycolumno++;  return(Y_EQ); }
">="      { yycolumno+=2;  return(Y_GE); }
">"       { yycolumno++;  return(Y_GT); }
"/\\"     { yycolumno+=2;  return(Y_INTERSECT); }
"("	  { yycolumno++;  return(Y_LBRACK); }
"|["      { yycolumno+=2;  return(Y_LCSCOPE); }
"{"       { yycolumno++;  return(Y_LCUR); }
"<="      { yycolumno+=2;  return(Y_LE); }
"|("      { yycolumno+=2;  return(Y_LOSCOPE); }
"["       { yycolumno++;  return(Y_LSQ); }
"<"       { yycolumno++;  return(Y_LT); }
"-"	  { yycolumno++;  return(Y_MINUS); }
"--"      { yycolumno+=2;  return(Y_MINUSMINUS); }
"/="      { yycolumno+=2;  return(Y_NE); }
"#"       { yycolumno++;  return(Y_NUMBERSIGN); }
"||"      { yycolumno+=2;  return(Y_PARA); }
"|"       { yycolumno++;  return(Y_PIPE); }
"+"       { yycolumno++;  return(Y_PLUS); }
"++"      { yycolumno+=2;  return(Y_PLUSPLUS); }
")"       { yycolumno++;  return(Y_RBRACK); }
"]|"      { yycolumno+=2;  return(Y_RCSCOPE); }
"}"	  { yycolumno++;  return(Y_RCUR); }
"?"       { yycolumno++;  return(Y_RECV); }
")|"      { yycolumno+=2;  return(Y_ROSCOPE); }
"]"	  { yycolumno++;  return(Y_RSQ); }
";"       { yycolumno++;  return(Y_SEMICOLON); }
"!"	  { yycolumno++;  return(Y_SEND); }
"/"       { yycolumno++;  return(Y_SLASH); }
"\\/"     { yycolumno+=2;  return(Y_UNION); }


{letter}(_?{letter_or_digit})* {
	int itoken;                        
    int fntoken;	
    Operator op;
	itoken=find_mc(yytext); 
	fntoken=find_fn(yytext, op);
	
	if (itoken== -1 && fntoken== -1)                  
	{                                   
		yylval.Identifier_data.pos = yycolumno;
		yylval.Identifier_data.len = strlen(yytext);
		yylval.Identifier_data.line= yylineno;
		yylval.Identifier_data.name= strdup(yytext);

		yycolumno += strlen(yytext);
		return (Y_IDENTIFIER);

	} else if (itoken != -1 ){
		yycolumno += strlen(yytext);
	       	return itoken;
	}
	else                    
	{                                   
		yylval.FuncId_data.pos = yycolumno;
		yylval.FuncId_data.len = strlen(yytext);
		yylval.FuncId_data.line= yylineno;
		yylval.FuncId_data.name= strdup(yytext);
		yylval.FuncId_data.op= op;

		yycolumno += strlen(yytext);
		return (Y_FUNCIDENTIFIER);
	}
}

{Y_NUMDOTNUM}	{
	yylval.Constant_data.pos = yycolumno;
	yylval.Constant_data.len = strlen(yytext);
	yylval.Constant_data.line = yylineno;
	yylval.Constant_data.value = atof(yytext);
	yycolumno += strlen(yytext);

	return (Y_NUMDOTNUM);
}

{Y_NUMBERCONST} {	
	yylval.Constant_data.pos = yycolumno;
	yylval.Constant_data.len = strlen(yytext);
	yylval.Constant_data.line = yylineno;
	yylval.Constant_data.value = atoi(yytext);
	yycolumno += strlen(yytext);
	
	return (Y_NUMBERCONST);
}

{Y_REALCONST}	{	
	yylval.Constant_data.pos = yycolumno;
	yylval.Constant_data.len = strlen(yytext);
	yylval.Constant_data.line = yylineno;
	yylval.Constant_data.value = atof(yytext);
	yycolumno += strlen(yytext); 
	
	return (Y_REALCONST);
}

{Y_STRINGCONST}	{
	yylval.Identifier_data.pos = yycolumno;
	yylval.Identifier_data.len = strlen(yytext);
	yylval.Identifier_data.line = yylineno;
	yycolumno += strlen(yytext);
	
	yylval.Identifier_data.name= strdup(yytext+1);
	if (yylval.Identifier_data.name[strlen(yytext)-2] != '"')
		printf("WARNING Y_STRINGCONST\n");
	else
		yylval.Identifier_data.name[strlen(yytext)-2] = 0;
	return (Y_STRINGCONST);
		}

\n	{
//	yylineno++;
//	yycolumno = 1;
	
	}

%%
