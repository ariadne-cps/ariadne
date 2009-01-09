#include <iostream>
#include <iomanip>
#include <fenv.h>

#if defined C_ROUNDING
typedef int rounding_mode_t;
#else
typedef unsigned short rounding_mode_t;
#endif

#define MCW_EM          0x003f          /* interrupt exception masks   */
#define EM_INVALID      0x0001          /*   invalid                   */
#define EM_DENORMAL     0x0002          /*   denormal                  */
#define EM_ZERODIVIDE   0x0004          /*   zero divide               */
#define EM_OVERFLOW     0x0008          /*   overflow                  */
#define EM_UNDERFLOW    0x0010          /*   underflow                 */
#define EM_INEXACT      0x0020          /*   inexact (precision)       */

#define MCW_IC          0x1000          /* infinity control            */
#define IC_AFFINE       0x1000          /*   affine                    */
#define IC_PROJECTIVE   0x0000          /*   projective                */

#define MCW_RC          0x0c00          /*  rounding control           */
#define RC_CHOP         0x0c00          /*    chop                     */
#define RC_UP           0x0800          /*    up                       */
#define RC_DOWN         0x0400          /*    down                     */
#define RC_NEAR         0x0000          /*    near                     */

#define MCW_PC          0x0300          /*  precision control          */
#define PC_24           0x0000          /*    24 bits                  */
#define PC_53           0x0200          /*    53 bits                  */
#define PC_64           0x0300          /*    64 bits                  */

/**** math coprocessor default control word and rounding modes (80x87) */

#define CW_DEFAULT\
        (IC_AFFINE      | RC_NEAR       | PC_64         |\
         EM_DENORMAL    | EM_OVERFLOW   | EM_UNDERFLOW  | EM_INEXACT)

#define CW_ROUND_CHOP   ((CW_DEFAULT & ~MCW_RC) | RC_CHOP)
#define CW_ROUND_UP     ((CW_DEFAULT & ~MCW_RC) | RC_UP)
#define CW_ROUND_DOWN   ((CW_DEFAULT & ~MCW_RC) | RC_DOWN)
#define CW_ROUND_NEAR   ((CW_DEFAULT & ~MCW_RC) | RC_NEAR)

//unsigned short ROUND_UP      = CW_ROUND_UP;
//unsigned short ROUND_DOWN    = CW_ROUND_DOWN;
//unsigned short ROUND_NEAREST = CW_ROUND_NEAR;
unsigned short ROUND_UP      = 2943;
unsigned short ROUND_DOWN    = 1919;
unsigned short ROUND_NEAREST = 895;

static volatile rounding_mode_t grnd;
inline rounding_mode_t asm_get_control_word() { asm volatile ("fstcw grnd"); return grnd; }
inline rounding_mode_t c_get_rounding_mode() { return fegetround(); }
inline void c_set_rounding_mode(rounding_mode_t rnd) { fesetround(rnd); }

inline void c_set_round_nearest() { fesetround(FE_TONEAREST); }
double div(volatile double x, volatile double y);

#if defined C_ROUNDING

inline void set_round_up() { fesetround(FE_UPWARD);  }
inline void set_round_down() { fesetround(FE_DOWNWARD);  }
inline void set_round_nearest() { fesetround(FE_TONEAREST);  }

inline void set_rounding_mode(rounding_mode_t rnd) { fesetround(rnd); }
inline rounding_mode_t get_rounding_mode() { return fegetround(); }

#elif defined ASM_ROUNDING

unsigned int rnd_up = 0x1a3f;
void set_round_up() { asm ("fldcw rnd_up"); }

//inline void set_round_up() { asm volatile ("fwait\n\t" "fldcw ROUND_UP\n\t" "fwait\n\t"); }
inline void set_round_down() { asm volatile ("fldcw ROUND_DOWN"); }
inline void set_round_nearest() { asm volatile ("fldcw ROUND_NEAREST"); }

//inline void set_rounding_mode(rounding_mode_t rnd) { asm volatile ("fldcw (%0)" : : "r" (rnd) ); }
//inline rounding_mode_t get_rounding_mode() { rounding_mode_t rnd asm volatile ("fldcw (%0)" : : "r" (rnd) ); return rnd; }

inline void set_rounding_mode(rounding_mode_t rnd) { grnd=rnd; asm volatile ("fldcw grnd"); } 
inline rounding_mode_t get_rounding_mode() { asm volatile ("fstcw grnd"); return grnd; }

#elif defined ASM_MACRO_ROUNDING

#define set_round_up()           asm volatile ("fldcw ROUND_UP")   
#define set_round_down()         asm volatile ("fldcw ROUND_DOWN")   
#define set_round_nearest()      asm volatile ("fldcw ROUND_NEAREST")   

inline void set_rounding_mode(rounding_mode_t rnd) { fesetround(rnd); }
inline rounding_mode_t get_rounding_mode() { return fegetround(); }

/*
#elif defined ASM_LIB_ROUNDING

void set_round_up();
void set_round_down();
void set_round_nearest();

void set_rounding_mode(rounding_mode_t rnd);
rounding_mode_t get_rounding_mode();
*/
#endif

void test_rounding(volatile double p, volatile double q) 
{
    std::cout<<"Testing correct rounding\n";

    rounding_mode_t fcw=asm_get_control_word();
    std::cout<<"Initial control word="<<fcw<<"\n";
    rounding_mode_t rnd=get_rounding_mode();
    std::cout<<"Initial rounding mode="<<rnd<<"\n";
    
    set_round_down();
    std::cout<<"Down rounding mode="<<get_rounding_mode()<<"\n";
    std::cout<<"Down C rounding mode="<<c_get_rounding_mode()<<"\n";
    volatile double xl=p/q;
    std::cout<<"  Computed "<<p<<"/"<<q<<"="<<xl<<std::endl;
    
    set_round_up();
    std::cout<<"Up rounding mode="<<get_rounding_mode()<<"\n";
    std::cout<<"Up C rounding mode="<<c_get_rounding_mode()<<"\n";
    volatile double xu=p/q;
    std::cout<<"  Computed "<<p<<"/"<<q<<"="<<xu<<std::endl;
    
    set_round_nearest();
    std::cout<<"Nearest rounding mode="<<get_rounding_mode()<<"\n";
    volatile double xn=p/q;
    std::cout<<"  Computed "<<p<<"/"<<q<<"="<<xn<<std::endl;
    
    std::cout<<p<<"/"<<q<<" ~ "<<xn<<std::endl;
    std::cout<<xl<<" < "<<p<<"/"<<q<<" < "<<xu<<std::endl;
    set_rounding_mode(rnd);
    std::cout<<"Restored rounding mode="<<get_rounding_mode()<<"\n";
    std::cout<<std::endl<<std::endl;
}


int main(int argc, const char* argv[]) {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    std::cout<<"sizeof(rounding_mode_t)="<<sizeof(rounding_mode_t)<<"\n";
    std::cout<<"sizeof(fenv_t)="<<sizeof(fenv_t)<<"\n";

    double x=atoi(argv[1]);
    double y=atoi(argv[2]);

    test_rounding(x,y);
}
