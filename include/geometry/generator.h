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
 
#ifndef _ARIADNE_GENERATOR_H
#define _ARIADNE_GENERATOR_H

  
namespace Ariadne {

namespace Geometry {
  template<typename R> class Point;
  template<typename R> class Polyhedron;
  template<typename R> Polyhedron<R> minkowski_sum(const Polyhedron<R>& A,
                                                   const Polyhedron<R>& B);
}
  
namespace LinearAlgebra {

/*! \brief A generating structure for a polyhedron in terms of points and rays. */
template <typename R> 
class GeneratorSystem {
 private:
 public:
  enum GeneratorType {
    POINT,
    LINE,
    RAY,
    CLOSURE_POINT
  };
  
  typedef boost::numeric::ublas::matrix<R> matrix_type;
  typedef boost::numeric::ublas::vector<R> vector_type;
  typedef Ariadne::Geometry::Point<R> state_type;

  typedef std::vector<GeneratorType> generator_type_list_type;
  
  GeneratorSystem() { }
  
  
  GeneratorSystem(const matrix_type& p_matrix, 
                  const generator_type_list_type& ptype_vector,
                  const matrix_type& r_matrix,
                  const generator_type_list_type& rtype_vector)
    : _points(p_matrix), _point_types(ptype_vector),
      _rays(r_matrix), _ray_types(rtype_vector) 
  { }
  
  GeneratorSystem(const Parma_Polyhedra_Library::Generator_System& ppl_gen);
  Parma_Polyhedra_Library::Generator_System ppl_generator_system() const;
  operator Parma_Polyhedra_Library::Generator_System() const {
    return ppl_generator_system(); }
    
  bool empty() const {
    return ( (this->number_of_points()==0) && (this->number_of_rays()==0) );
  }
  
  size_type space_dimension() const {
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
    matrix_type new_point_matrix(dim,number_of_points);
    generator_type_list_type new_point_type_vector;
    
    new_point_type_vector.resize(number_of_points);
    
    matrix_type new_ray_matrix(dim,ray_nb);
    generator_type_list_type new_ray_type_vector;
    
    new_ray_type_vector.resize(ray_nb);
    
    this->_points=new_point_matrix;
    this->_point_types=new_point_type_vector;
    
    this->_rays=new_ray_matrix;
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
  
  size_t dim=space_dimension();
  
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

#endif /* _ARIADNE_GENERATOR_H */
