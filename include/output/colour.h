/****************************************************************************
 *            colour.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#ifndef ARIADNE_COLOUR_H
#define ARIADNE_COLOUR_H

/*! \file colour.h
 *  \brief Colours
 */
 
#include <string>

namespace Ariadne {
  
    
     /*!\brief A class representing a colour, with a \a name and \a rgb values. 
      * The default constructor creates a fully transparant colour. */
     class Colour {
     public:
      Colour()
        : _name("transparant"), _transparant(true) { }
      Colour(const char* name);
      Colour(const char* name, unsigned char red, unsigned char green, unsigned char blue) 
        : _name(name), _red(red), _green(green), _blue(blue), _transparant(false) { }
      std::string name() const { return this->_name; }
      unsigned short int red() const { return this->_red; }
      unsigned short int green() const { return this->_green; }
      unsigned short int blue() const { return this->_blue; }
      bool transparant() const { return this->_transparant; }
     private:
      std::string _name;
      unsigned char _red, _green, _blue;
      bool _transparant;
    };

    std::ostream& operator<<(std::ostream& os, const Colour& c);

    extern const Colour transparant;

    extern const Colour white;
    extern const Colour black;
    extern const Colour red;
    extern const Colour green;
    extern const Colour blue;
    extern const Colour yellow;
    extern const Colour cyan;
    extern const Colour magenta;

  
} // namespace Ariadne


#endif /* ARIADNE_COLOUR_H */
