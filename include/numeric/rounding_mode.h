/***************************************************************************
 *            numeric/rounding_mode.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file numeric/set_rounding_mode.h
 *  \brief Type definitions and conversion operators for 64-bit fixed precision floating point numbers.
 */

#ifndef ARIADNE_NUMERIC_ROUNDING_MODE_H
#define ARIADNE_NUMERIC_ROUNDING_MODE_H

#include "numeric/rounding.h"


/**** math coprocessor control word (80x87) ****************************/

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

#define ARIADNE_CW_DEFAULT\
        (IC_AFFINE      | RC_NEAR       | PC_64         |\
         EM_DENORMAL    | EM_OVERFLOW   | EM_UNDERFLOW  | EM_INEXACT)

#define ARIADNE_CW_ROUND_CHOP   ((ARIADNE_CW_DEFAULT & ~MCW_RC) | RC_CHOP)
#define ARIADNE_CW_ROUND_UP     ((ARIADNE_CW_DEFAULT & ~MCW_RC) | RC_UP)
#define ARIADNE_CW_ROUND_DOWN   ((ARIADNE_CW_DEFAULT & ~MCW_RC) | RC_DOWN)
#define ARIADNE_CW_ROUND_NEAR   ((ARIADNE_CW_DEFAULT & ~MCW_RC) | RC_NEAR)

namespace Ariadne {
namespace Numeric {

typedef unsigned short rounding_mode_type;

const rounding_mode_type CwRoundUp      = ARIADNE_CW_ROUND_UP;
const rounding_mode_type CwRoundDown    = ARIADNE_CW_ROUND_DOWN;
const rounding_mode_type CwRoundNear    = ARIADNE_CW_ROUND_NEAR;
const rounding_mode_type CwRoundChop    = ARIADNE_CW_ROUND_CHOP;

template<class Rnd> void set_rounding_mode();

}
}

#endif
