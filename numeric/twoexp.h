/***************************************************************************
 *            numeric/twoexp.h
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

/*! \file numeric/twoexp.h
 *  \brief 
 */



#ifndef ARIADNE_TWOEXP_H
#define ARIADNE_TWOEXP_H

namespace Ariadne {

class TwoExp {
    Int _exp;
  public:
    TwoExp(Int exp) : _exp(n);
    Int exponent() const { return this->_exp; }
};

} // namespace Ariadne

#endif
