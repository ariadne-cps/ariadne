/***************************************************************************
 *            constraint.h
 *
 *  Thu Feb  3 09:31:28 2005
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
 
#ifndef _CONSTRAINT_H
#define _CONSTRAINT_H

namespace Ariadne {	
namespace LinearAlgebra {
	
template <typename R> 
class ConstrainSystem {
	
	public:
		
		enum RelType{
			EQUALITY,
			INEQUALITY,
			STRICT_INEQUALITY
		};
		
		typedef boost::numeric::ublas::matrix<R> Matrix;
		typedef boost::numeric::ublas::vector<R> Vector;
		typedef std::vector<RelType> Rel_Vector;
		
				
		Matrix C;
		Vector d;	
		
		Rel_Vector r_type;
		
		ConstrainSystem(Matrix new_C, Vector new_d , Rel_Vector new_r_type):
					C(new_C),d(new_d),r_type(new_r_type){}
						
		ConstrainSystem() {}
		
		inline size_t nb_equality() const{
			
			size_t equality_num=0;
			
			for (size_t i=0; i< this->nb_constraints(); i++) {
				
				if (this->r_type[i]==EQUALITY) equality_num++;
			}
		
			return equality_num;
		}
		
		inline bool already_at_precision(const R &delta) const{			
			
			size_t i,j;
			R min;
			R too_big=delta*delta*delta;
			
			for (j=0; j< this->nb_constraints(); j++) {
				
				min=0.0;
				
				if (abs(this->d(j))!=0)
					min=abs(this->d(j));
				
				for (i=0; i< this->space_dim(); i++) {
					if  ((abs(this->C(j,i))!=0)&&((min==0)||(abs(this->C(j,i))<min))) {
						min=abs(this->C(j,i));
					}
				}
				
				if (min>too_big) {
					return false;	
				}
			}
			
			return true;
		}
		
		inline void reduce_precision_to_expanding(const R &delta){			
			
			size_t i,j;
			
			R too_big=delta*delta*delta;
			
			R min;
			
			for (j=0; j< this->nb_constraints(); j++) {
				
				min=0.0;
				
				if (abs(this->d(j))!=0)
					min=abs(this->d(j));
				
				for (i=0; i< this->space_dim(); i++) {
					if  ((abs(this->C(j,i))!=0)&&((min==0)||(abs(this->C(j,i))<min))) {
						min=abs(this->C(j,i));
					}
				}
				
				if (min>too_big) {
					
					while (min>too_big) {
						min=numerator(min)/numerator(delta);
						
						for (i=0; i< this->space_dim(); i++) {
							this->C(j,i)=numerator(this->C(j,i))/numerator(delta);
						}
						
						this->d(j)=numerator(this->d(j))/numerator(delta);
					
					}
					
					this->d(j)--;

				}
				
			}
		}
		
		inline void reduce_precision_to_shrinking(const R &delta){			
			
			size_t i,j;
			
			R min;
			
			for (j=0; j< this->nb_constraints(); j++) {
				
				min=0.0;
				
				if (abs(this->d(j))!=0)
					min=abs(this->d(j));
				
				for (i=0; i< this->space_dim(); i++) {
					if  ((abs(this->C(j,i))!=0)&&((min==0)||(abs(this->C(j,i))<min))) {
						min=abs(this->C(j,i));
					}
				}
				
				if (min>delta) {
					
					while (min>delta) {
						min=numerator(min)/numerator(delta);
						
						for (i=0; i< this->space_dim(); i++) {
							this->C(j,i)=numerator(this->C(j,i))/numerator(delta);
						}
						
						this->d(j)=numerator(this->d(j))/numerator(delta);
					
					}
					
					this->d(j)--;

				}
				
			}
		}
		
		inline void divide_all_by(const R &delta){
			
			size_t i,j;
			
			for (j=0; j< this->nb_constraints(); j++) {
				
				for (i=0; i< this->space_dim(); i++) {
					this->C(j,i)=numerator(this->C(j,i))/numerator(delta);
				}
				
				this->d(j)=numerator(this->d(j))/numerator(delta);
				
			}
		}
		
		inline const ConstrainSystem<R> &expand_by(const R &delta) {
			
			if (delta==0) return *this;
			
			const size_t &dim=this->space_dim();
			
			if (!this->any_equality()) {
				
				for (size_t j=0; j< this->nb_constraints(); j++) {
					this->d(j)+=delta;
				}
				
				return *this;
			}
			
			if (delta<0) {
				Matrix new_C(1, dim);
				Vector new_d(1);
				
				this->r_type.resize(1);
				
				for (size_t i=0; i< dim; i++) {
					new_C(1,i)=0;
				}
				
				new_d(1)=-1;
				
				this->C=new_C;
				this->d=new_d;
				
				return *this;				
			} else {
				
				return this->expand_equality_by(delta);
			}
			
		}
		
		
		inline bool any_equality() const{
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			for (size_t i=0; i< this->nb_constraints(); i++) {
				
				if (this->r_type[i]==EQUALITY) {
					
					#ifdef DEBUG
						std::cout << __FILE__ << ":" << __LINE__ << std::endl;
					#endif
		
					return true;
				}
			}
			
			#ifdef DEBUG
				std::cout << __FILE__ << ":" << __LINE__ << std::endl;
			#endif
			
			return false;
		}
		
		inline bool is_equality(const size_t &i) const { 
			return(this->r_type[i]==EQUALITY); 
		}
		
		inline void set_equality(const size_t &i) {
			this->r_type[i]=EQUALITY;
		}
		
		inline bool empty() const{
			return (this->C.size1()==0);
		}
		
		inline size_t space_dim() const{
				
			return this->C.size2();
		}
		
		inline size_t nb_constraints() const{
				
			return this->C.size1();
			
		}
		
		inline const ConstrainSystem<R> &operator=(
					const ConstrainSystem<R> &cs){
				
			(this->C) = (cs.C);
			(this->d) = (cs.d);
			(this->r_type) = (cs.r_type);
					
			return *this;
		}
		
		inline void reset_dimensions(const size_t &dim, const size_t &nb_constraints) {
		
			Matrix new_C(nb_constraints,dim);
			Vector new_d(nb_constraints);
		
			Rel_Vector new_r_type(nb_constraints);

			C=new_C;
			d=new_d;

			r_type=	new_r_type;
		}
		
		inline const ConstrainSystem<R> &expand_equality_by(const R &delta) {
			
			const size_t &dim=this->space_dim();
			size_t new_constr=this->nb_equality()+this->nb_constraints();
			size_t new_j=0, i,j;
			
			R den=denumerator(delta), num=numerator(delta);
			
			
			Matrix new_C(new_constr, dim);
			Vector new_d(new_constr);
			
			/* TODO: REIMPLEMENT */
			for (j=0; j< this->nb_constraints(); j++) {
				
				if (!this->is_equality(j)) {
					
					for (i=0; i< dim; i++) {
						new_C(new_j,i)=this->C(j,i);
					}
					
					new_d(new_j)=this->d(j);
					
					new_j++;
					
				} else {
				
					for (i=0; i< dim; i++) {
						new_C(new_j,i)=-this->C(j,i)*den;
						new_C(new_j+1,i)=this->C(j,i)*den;
					}
					
					new_d(new_j)=-this->d(j)*den-num;
					new_d(new_j+1)=this->d(j)*den-num;
					
					new_j+=2;
				}
				
			}
			
			this->C=new_C;
			this->d=new_d;
			
			this->r_type.resize(new_constr);
			
			for (i=0; i< dim; i++) {
				this->r_type[i]=INEQUALITY;
			}
			
			return *this;
			
		}
	
	private:
		
};
		
}
}

#endif /* _CONSTRAINT_H */
