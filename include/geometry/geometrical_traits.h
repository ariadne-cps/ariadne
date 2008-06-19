/***************************************************************************
 *            geometrical_traits.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file geometrical_traits.h
 *  \brief Traits classes to define properties of gemetric types.
 */

#ifndef ARIADNE_GEOMETRICAL_TRAITS_H
#define ARIADNE_GEOMETRICAL_TRAITS_H


namespace Ariadne {
  
    
    class basic_set_tag { };
    class denotable_set_tag { };
    class list_set_tag : public denotable_set_tag { };
    class abstract_set_tag { };

    template<class S> class geometrical_traits;

    template<class S>
    class geometrical_traits
    {
     public:
      typedef typename S::set_category set_category;
    };

} // namespace Ariadne
  

#endif /* ARIADNE_GEOMETRICAL_TRAITS_H */
