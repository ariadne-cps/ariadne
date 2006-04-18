/***************************************************************************
 *            generator_system.tpl
 *
 *  Thu Feb  16 08:54:28 2005
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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

#include "../numeric/arithmetic.h"

#include "generator_system.h"

namespace Ariadne {
  namespace LinearAlgebra {
  
    template <typename R>
    GeneratorSystem<R>::GeneratorSystem(const Parma_Polyhedra_Library::Generator_System& ppl_gen)
    {
      Ariadne::LinearAlgebra::GeneratorSystem<R>& gen(*this);
      typedef R real_type;
            
      Parma_Polyhedra_Library::Generator_System::const_iterator j_gen, begin, end;
            
      begin=ppl_gen.begin();
      end=ppl_gen.end();
                
      size_t space_dim=ppl_gen.space_dimension();
      size_t number_of_points=0,ray_nb=0, i,j_ray=0, j_point=0;
                
      for (j_gen=begin; j_gen!=end; j_gen++) {
        if ((j_gen->is_line())||(j_gen->is_ray())) {
          ray_nb++;
        } else {
          number_of_points++;  
        }
      }
      
      gen.reset_dimensions(number_of_points,ray_nb,space_dim);
    
      for (j_gen=begin; j_gen!=end; j_gen++) {
        if (j_gen->is_line()) {
          for (i=0; i< space_dim; ++i) {
            gen._rays(i,j_ray)=j_gen->coefficient(Parma_Polyhedra_Library::Variable(i));
          }
          gen.set_ray_type(j_ray,LINE);
          j_ray++;
        } 
        else {
          if (j_gen->is_ray()) {
            for (i=0; i< space_dim; ++i) {
              gen._rays(i,j_ray)=j_gen->coefficient(Parma_Polyhedra_Library::Variable(i));
            }
            gen.set_ray_type(j_ray,RAY);
            j_ray++;
          } 
          else {
            if (j_gen->is_point()) {
              real_type den=j_gen->divisor();
              for (i=0; i< space_dim; ++i) {
                gen._points(i,j_point)=j_gen->coefficient(Parma_Polyhedra_Library::Variable(i))/den;
              }
              gen.set_point_type(j_point,POINT);
              j_point++;
            } 
            else { /*j_gen is a closure _points */
              real_type den=j_gen->divisor();
              for (i=0; i< space_dim; ++i) {
                gen._points(i,j_point)=j_gen->coefficient(Parma_Polyhedra_Library::Variable(i))/den;
              }
              gen.set_point_type(j_point,CLOSURE_POINT);
              j_point++;
            }
          }
        }
      }
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif  
    }
    
    
    
    
    
    
    template <typename R>
    Parma_Polyhedra_Library::Generator_System
    GeneratorSystem<R>::ppl_generator_system() const
    {
      Parma_Polyhedra_Library::Generator_System ppl_gen;
      Parma_Polyhedra_Library::Linear_Expression ppl_lin_expr;
      const GeneratorSystem<R>& gen(*this);
      
      size_t dim=dimension();
      
      #ifdef DEBUG
        std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      #endif      
    
      R den;
            
      for (size_t i=0; i< gen.number_of_points(); ++i) {
        den=denominator(gen._points(0,i));
        for (size_t j=1; j< dim; ++j) {
          den=lcm(numerator(den),denominator(gen._points(j,i)));
        }
        ppl_lin_expr=Rational(den*gen._points(0,i))*Parma_Polyhedra_Library::Variable(0);
        for (size_t j=1; j< dim; ++j) {
          ppl_lin_expr+=Rational(den*gen._points(j,i))*Parma_Polyhedra_Library::Variable(j);
        }
        if (gen.is_point(i)) {    
          ppl_gen.insert(Parma_Polyhedra_Library::Generator::point(ppl_lin_expr,Rational(den)));
        } 
        else {
          ppl_gen.insert(Parma_Polyhedra_Library::Generator::closure_point(ppl_lin_expr,Rational(den)));
        }
      }
      
      for (size_t i=0; i< gen.number_of_rays(); ++i) {
        ppl_lin_expr=Rational(gen._rays(0,i))*Parma_Polyhedra_Library::Variable(0);
        for (size_t j=1; j<dim; ++j) {
          ppl_lin_expr+=Rational(gen._rays(j,i))*Parma_Polyhedra_Library::Variable(j);
        }
        if (gen.is_ray(i)) {    
          ppl_gen.insert(Parma_Polyhedra_Library::Generator::ray(ppl_lin_expr));
        } 
        else {
          ppl_gen.insert(Parma_Polyhedra_Library::Generator::line(ppl_lin_expr));
        }
      }
      
      return ppl_gen;
    }
    
    
  }
}
