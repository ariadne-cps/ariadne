/***************************************************************************
 *            numeric/logical.decl.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/logical.decl.h
 *  \brief
 */



#ifndef ARIADNE_LOGICAL_DECL_H
#define ARIADNE_LOGICAL_DECL_H

#include "numeric/paradigm.h"

namespace Ariadne {

typedef Void Void;

enum class LogicalValue : char;
template<class P=Void> class Logical;

//typedef Logical<Exact> Boolean;
//typedef Logical<Validated> Kleenean;
//typedef Logical<Upper> Sierpinski;
//typedef Logical<Approximate> Fuzzy;

class Boolean;
class Kleenean;
class Sierpinski;
class NegSierpinski;
class Fuzzy;

template<class P> using LogicalType = Logical<P>;
using ExactLogic = Logical<Exact>;
using EffectiveLogic = Logical<Effective>;
using ValidatedLogic = Logical<Validated>;
using UpperLogic = Logical<Upper>;
using LowerLogic = Logical<Lower>;
using ApproximateLogic = Logical<Approximate>;

using Decidable = Logical<Exact>;
using Quasidecidable = Logical<Effective>;
using Verifyable = Logical<EffectiveUpper>;
using Falsifyable = Logical<EffectiveLower>;

}

#endif
