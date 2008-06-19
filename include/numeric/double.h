/***************************************************************************
 *            numeric/double.h
 *
 *  Copyright  2004-7  Pieter Collins
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
 
/*! \file numeric/double.h
 *  \brief Rounded arithmetic for double precision
 */

#ifndef ARIADNE_NUMERIC_DOUBLE_H
#define ARIADNE_NUMERIC_DOUBLE_H

namespace Ariadne {


bool initialise();

double med(double x, double y);
double rad(double x, double y);
double pow(double x, unsigned int n);
double hypot(double x, double y);

}

#include "double.inline.h"

#endif /* ARIADNE_NUMERIC_DOUBLE_H */
