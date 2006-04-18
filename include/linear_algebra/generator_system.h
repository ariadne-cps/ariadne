/***************************************************************************
 *            generator.h
 *
 *  Thu Feb  16 08:54:28 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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

/*! \file
 *  \brief Points and rays generating a polyhedron.
 */
 
#ifndef _ARIADNE_GENERATOR_SYSTEM_H
#define _ARIADNE_GENERATOR_SYSTEM_H

#include <ppl.hh>

#include "../declarations.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"

namespace Ariadne { namespace Geometry {
    template<typename R>
    Polyhedron<R> 
    minkowski_sum(const Polyhedron<R>& A, 
                  const Polyhedron<R>& B);
}}

namespace Ariadne {
  namespace LinearAlgebra {
    
    
    /*! \brief A generating structure for a polyhedron in terms of points and rays. */
    template <typename R> 
    class GeneratorSystem {
     private:
      friend Geometry::Polyhedron<R> 
      Geometry::minkowski_sum<>(const Geometry::Polyhedron<R>& A, 
                                const Geometry::Polyhedron<R>& B);
     public:
      enum GeneratorType {
        POINT,
        LINE,
        RAY,
        CLOSURE_POINT
      };
      
      typedef Matrix<R> matrix_type;
      typedef Vector<R> vector_type;
      typedef Geometry::Point<R> state_type;
    
      typedef std::vector<GeneratorType> generator_type_list_type;
      
      GeneratorSystem() { }
      
      
      GeneratorSystem(const matrix_type& p_Matrix, 
                      const generator_type_list_type& ptype_vector,
                      const matrix_type& r_Matrix,
                      const generator_type_list_type& rtype_vector)
        : _points(p_Matrix), _point_types(ptype_vector),
          _rays(r_Matrix), _ray_types(rtype_vector) 
      { }
      
      GeneratorSystem(const Parma_Polyhedra_Library::Generator_System& ppl_gen);
      Parma_Polyhedra_Library::Generator_System ppl_generator_system() const;
      operator Parma_Polyhedra_Library::Generator_System() const {
        return ppl_generator_system(); }
        
      bool empty() const {
        return ( (this->number_of_points()==0) && (this->number_of_rays()==0) );
      }
      
      size_type dimension() const {
        return this->_points.size1();
      }
      
      size_type number_of_points() const {
        return this->_points.size2();    
      }
      
      size_type number_of_rays() const {
        return this->_rays.size2();
      }
      
      void set_point_type(const size_type& n, const GeneratorType& gt) {
        this->_point_types[n]=gt;
      }
      
      void set_ray_type(const size_type& n, const GeneratorType& gt) {
        this->_ray_types[n]=gt;
      }
     
      bool is_point(const size_type& n) const {
        return (this->_point_types[n]==POINT);
      }
      
      bool is_closure_point(const size_type& n) const {
        return (this->_point_types[n]==CLOSURE_POINT);
      }
      
      bool is_ray(const size_type& n) const {
        return (this->_ray_types[n]==RAY);
      }
      
      bool is_line(const size_type& n) const {
        return (this->_ray_types[n]==LINE);
      }
      
      const GeneratorSystem<R>& sum_vector_to_all_points(const vector_type &v) {
        #ifdef DEBUG
          std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif
        for (size_type i=0; i<this->dimension(); ++i) {
          const R& c=v(i);
          for (size_type j=0; j<this->number_of_points(); ++j) {
            this->_points(i,j) += c;
          }
        }
        #ifdef DEBUG
            std::cout << __FILE__ << ":" << __LINE__ << std::endl;
        #endif
        return *this;
      }
      
      GeneratorSystem<R>& operator=(const GeneratorSystem<R>& gen) {
        if(this!=&gen) {
          this->_points =  gen._points;
          this->_point_types =  gen._point_types;
          this->_rays =  gen._rays;
          this->_ray_types =  gen._ray_types;
        }
        return *this;
      }
      
      void reset_dimensions(const size_type& number_of_points, 
                                   const size_type& ray_nb, 
                                   const size_type& dim) 
      {
        matrix_type new_point_Matrix(dim,number_of_points);
        generator_type_list_type new_point_type_vector;
        
        new_point_type_vector.resize(number_of_points);
        
        matrix_type new_ray_Matrix(dim,ray_nb);
        generator_type_list_type new_ray_type_vector;
        
        new_ray_type_vector.resize(ray_nb);
        
        this->_points=new_point_Matrix;
        this->_point_types=new_point_type_vector;
        
        this->_rays=new_ray_Matrix;
        this->_ray_types=new_ray_type_vector;
      }
     private:  
      friend class Ariadne::Geometry::Polyhedron<R>;
    /*  friend Ariadne::Geometry::Polyhedron<R> 
          Ariadne::Geometry::minkowski_sum<>(
              const Ariadne::Geometry::Polyhedron<R>& A,
              const Ariadne::Geometry::Polyhedron<R>& B);*/
     private:  
      matrix_type _points;
      generator_type_list_type _point_types;
      
      matrix_type _rays;
      generator_type_list_type _ray_types;
    };

  }
}

#endif /* _ARIADNE_GENERATOR_SYSTEM_H */
