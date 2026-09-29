/***************************************************************************
 *            numeric/logical_sequence.hpp
 *
 *  Copyright  2026  Ariadne contributors
 *
 ****************************************************************************/
#ifndef ARIADNE_NUMERIC_LOGICAL_SEQUENCE_HPP
#define ARIADNE_NUMERIC_LOGICAL_SEQUENCE_HPP
#include "foundation/logical.hpp"
namespace Ariadne {
template<class X> class Sequence;
LowerKleenean disjunction(Sequence<LowerKleenean> const&);
UpperKleenean conjunction(Sequence<UpperKleenean> const&);
}
#endif
