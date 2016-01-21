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

//typedef Logical<ExactTag> Boolean;
//typedef Logical<ValidatedTag> ValidatedKleenean;
//typedef Logical<UpperTag> ValidatedSierpinskian;
//typedef Logical<ApproximateTag> Fuzzy;

class Boolean;
class ValidatedKleenean;
class ValidatedSierpinskian;
class ValidatedNegatedSierpinskian;
class Fuzzy;

template<class P> using LogicalType = Logical<P>;
using ExactLogic = Logical<ExactTag>;
using EffectiveLogic = Logical<EffectiveTag>;
using EffectiveUpperLogic = Logical<EffectiveUpperTag>;
using EffectiveLowerLogic = Logical<EffectiveLowerTag>;
using ValidatedLogic = Logical<ValidatedTag>;
using ValidatedUpperLogic = Logical<ValidatedUpperTag>;
using ValidatedLowerLogic = Logical<ValidatedLowerTag>;
using ApproximateLogic = Logical<ApproximateTag>;

using Decidable = Logical<ExactTag>;
using Quasidecidable = Logical<EffectiveTag>;
using Verifyable = Logical<EffectiveUpperTag>;
using Falsifyable = Logical<EffectiveLowerTag>;

}

#endif
