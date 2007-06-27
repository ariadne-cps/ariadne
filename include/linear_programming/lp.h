/***************************************************************************
 *            lp.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
 ****************************************************************************/
/*
 * Based on the linear programming algorithms in PPL-0.8
 *   Copyright (C) 2001-2006 Roberto Bagnara <bagnara@cs.unipr.it>
 */

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
 
/*! \file lp.h
 *  \brief Linear programming include file.
 */

#ifndef ARIADNE_LP_H
#define ARIADNE_LP_H

#include <iosfwd>
#include <cassert>
#include <map>

#include "../base/tribool.h"
#include "../linear_algebra/matrix.h"
#include "../linear_algebra/permutation.h"
#include "../linear_algebra/exceptions.h"
#include "../output/logging.h"
#include "../numeric/float64.h"

#include "exceptions.h"
#include "lputil.h"
#include "lpstp.h"
#include "lptst.h"
#include "lpfsp.h"
#include "lpslv.h"

#endif /* ARIADNE_LP_H */
