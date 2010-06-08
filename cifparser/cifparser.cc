/**
 * This is the main code for cif2hif tool. Permitt translation of 
 * code from cif to hif/hbf format
 *
 * @version 1.0
 */

//library inclusion
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "ciftypes.h"
#include "cif_bison.h"

using namespace std;

//default create and setting of variables used by parser and lexer
extern FILE *yyin;
extern FILE *yyout;
extern int yy_flex_debug;
extern int yydebug;
int yycolumno = 1;

//void yyparse(SystemObject*);
int yyparse();

//stream use to send tool information to user
ostream* msgStream; 

int main(int argc, char** argv){
    msgStream = &cout;
        
	yydebug = 0;  //bison debug
	cout<<"\n\t\t ---- cif2ariadne 0.1 ----"<<endl;

    if(argc < 2) {
        cout << "USAGE: cifparser filename.cif [0|1]" << std::endl;
        return 0;
    }
    
    if(argc >= 3) { // second argument specify the vale for yydebug
        yydebug = atoi(argv[2]);
    }   

	yyin = fopen(argv[1], "r");
	if (!yyin) {
		cout << "ERROR: Could not open Cif file " << argv[1] << std::endl;
		return 0;
	}

    // disable output

	yy_flex_debug = 0;   // flex debug
	if(yyparse( ) == 0) {
	    cout<<"Operation successfully completed." <<endl;	
	} else{
	    cout<<"ERROR in parsing Cif file." <<endl;
    }	
	fclose(yyin);
}
