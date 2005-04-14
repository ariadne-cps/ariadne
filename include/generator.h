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
 
#ifndef _GENERATOR_H
#define _GENERATOR_H

namespace Ariadne {	
namespace LinearAlgebra {
	
template <typename R> 
class AriadneGeneratorSystem {
	
	private:

	
	public:
		
		enum GeneratorType{
			POINT,
			LINE,
			RAY,
			CLOSURE_POINT
		};
		
		typedef boost::numeric::ublas::matrix<R> Matrix;
		typedef boost::numeric::ublas::vector<R> Vector;
		typedef std::vector<GeneratorType> TypeVector;
		
		AriadneGeneratorSystem(const Matrix &pmatrix, 
					const TypeVector &ptype_vector,const Matrix &r_matrix,
					const TypeVector &rtype_vector):
				point(pmatrix), point_type(ptype_vector),
				ray(r_matrix),ray_type(rtype_vector){}
						
		AriadneGeneratorSystem() {}
		
		inline bool empty() const{
			return ((this->point_nb()==0)&&(this->ray_nb()==0));
		}
		
		inline size_t space_dim() const{
				
			return this->point.size1();
		}
		
		inline size_t point_nb() const{
				
			return this->point.size2();
			
		}
		
		inline size_t ray_nb() const{
				
			return this->ray.size2();
			
		}
		
		inline void set_point(const size_t &j){
				
			this->point_type[j]=POINT;
			
		}
		
		inline void set_closure_point(const size_t &j){
				
			this->point_type[j]=CLOSURE_POINT;
			
		}
		
		inline void set_ray(const size_t &j){
				
			this->ray_type[j]=RAY;
			
		}
		
		inline void set_line(const size_t &j){
				
			this->ray_type[j]=LINE;
			
		}
		
		inline bool is_point(const size_t &j) const {
				
			return (this->point_type[j]==POINT);
			
		}
		
		inline bool is_closure_point(const size_t &j) const {
				
			return (this->point_type[j]==CLOSURE_POINT);
			
		}
		
		inline bool is_ray(const size_t &j) const {
				
			return (this->ray_type[j]==RAY);
			
		}
		
		inline bool is_line(const size_t &j) const {
				
			return (this->ray_type[j]==LINE);
			
		}
		
		inline const AriadneGeneratorSystem<R> 
				&sum_vector_to_all_points(const Vector &v) {
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
					
			size_t i,j;

			for (j=0; j< this->space_dim(); j++) {
				
				const R &c=v(j);
				
				for (i=0; i< this->point_nb(); i++) {
					this->point(j,i)+=c;
				}
			}						
					
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif	
			
			return *this;
		}
		
		inline const AriadneGeneratorSystem<R> &operator=(
					const AriadneGeneratorSystem<R> &gen){
				
				(this->point) = (cs.point);
				(this->point_type) = (cs.point_type);
				(this->ray) = (cs.ray);
				(this->ray_type) = (cs.ray_type);
						
				return *this;
		}
		
		inline void reset_dimensions(const size_t &point_nb, const size_t &ray_nb, const size_t &dim) {
		
			Matrix new_point_matrix(dim,point_nb);
			TypeVector new_point_type_vector;
			
			new_point_type_vector.resize(point_nb);
			
			Matrix new_ray_matrix(dim,ray_nb);
			TypeVector new_ray_type_vector;
		
			new_ray_type_vector.resize(ray_nb);
			
			this->point=new_point_matrix;
			this->point_type=new_point_type_vector;
			
			this->ray=new_ray_matrix;
			this->ray_type=new_ray_type_vector;
		}
		
		Matrix point;
		TypeVector point_type;

		Matrix ray;
		TypeVector ray_type;	
		
};
		
}
}

#endif /* _CONSTRAINT_H */
