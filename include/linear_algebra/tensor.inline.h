/***************************************************************************
 *            tensor.inline.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it Pieter.Collins@cwi.nl
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
 

namespace Ariadne {
      
template<class R> inline
LinearAlgebra::Tensor<R>::Tensor(const array<size_type>& sz) 
  : _sizes(sz), _strides(compute_strides(sz)), _elements() 
{
  _elements.resize(this->number_of_elements()); 
}

      
template<class R> inline
LinearAlgebra::Tensor<R>::Tensor(const array<size_type>& sz, const R* ptr) 
  : _sizes(sz), _strides(compute_strides(sz)), _elements(ptr,ptr+this->number_of_elements()) 
{
}

      
template<class R> template<class Rl> inline
LinearAlgebra::Tensor<R>::Tensor(const Tensor<Rl>& t)
  : _sizes(t.sizes()), _strides(compute_strides(t.sizes())), _elements(t.begin(),t.begin()+t.number_of_elements()) 
{
}


template<class R> inline
LinearAlgebra::Tensor<R>::Tensor(const Tensor<R>& t) : _sizes(t._sizes), _strides(t._strides), _elements(t._elements) 
{
}


template<class R> inline
LinearAlgebra::Tensor<R>& 
LinearAlgebra::Tensor<R>::operator=(const Tensor<R>& t) 
{
  if(this!=&t) {
    this->_sizes=t._sizes;
    this->_strides=t._strides;
    this->_elements=t._elements;
  }
  return *this;
}        

template<class R> inline
size_type
LinearAlgebra::Tensor<R>::rank() const 
{
  return this->_sizes.size(); 
}

      
template<class R> inline
size_type
LinearAlgebra::Tensor<R>::number_of_elements() const 
{
  return this->_strides[0]*this->_sizes[0]; 
}


template<class R> inline
const array<size_type>&
LinearAlgebra::Tensor<R>::Tensor<R>::sizes() const 
{
  return this->_sizes; 
}


template<class R>
size_type
LinearAlgebra::Tensor<R>::size(const size_type& i) const 
{
  return this->_sizes[i]; 
}


template<class R>
const R&
LinearAlgebra::Tensor<R>::operator() (const index_array_type ind) const 
{
  return this->_elements[this->position(ind)]; 
}


      
template<class R>
R&
LinearAlgebra::Tensor<R>::operator() (const index_array_type ind) 
{
  return this->_elements[this->position(ind)]; 
}


    
template<class R>
R*
LinearAlgebra::Tensor<R>::begin() 
{
  return _elements.begin(); 
}


      
template<class R>
const R*
LinearAlgebra::Tensor<R>::begin() const 
{
  return _elements.begin(); 
}


      
template<class R> inline
LinearAlgebra::Tensor<R>::Tensor(size_type m, size_type n1, size_type n2)
  : _sizes(3), _strides(3), _elements(m*n1*n2,R(0)) 
{
  _sizes[0]=m; _sizes[1]=n1; _sizes[2]=n2; 
  _strides=compute_strides(_sizes); 
}


template<class R> inline
const R& 
LinearAlgebra::Tensor<R>::operator() (const size_type& i) const 
{
  if(this->rank()!=1) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return this->_elements[i];
}


template<class R>
R&
LinearAlgebra::Tensor<R>::operator() (const size_type& i) 
{
  if(this->rank()!=1) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return this->_elements[i];
}

      
template<class R>
const R&
LinearAlgebra::Tensor<R>::operator() (const size_type& i, const size_type& j) const 
{
  if(this->rank()!=2) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return this->_elements[i*this->_sizes[1]+j];
}
      

template<class R>
R&
LinearAlgebra::Tensor<R>::operator() (const size_type& i, const size_type& j) 
{
  if(this->rank()!=2) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return this->_elements[i*this->_sizes[1]+j];
}
      

template<class R>
const R&
LinearAlgebra::Tensor<R>::operator() (const size_type& i, const size_type& j, const size_type& k) const 
{
  if(this->rank()!=3) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return this->_elements[(i*_sizes[1]+j)*_sizes[2]+k];
}
      

template<class R>
R&
LinearAlgebra::Tensor<R>::operator() (const size_type& i, const size_type& j, const size_type& k) 
{
  if(this->rank()!=3) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return this->_elements[(i*_sizes[1]+j)*_sizes[2]+k];
}


template<class R>
size_type
LinearAlgebra::Tensor<R>::position(const array<size_type>& ind) const 
{
  size_type result=ind[0];
  for(size_type d=1; d!=this->rank(); ++d) { 
    result=this->_sizes[d-1]+ind[d]; 
  }
  return result;
}


template<class R>
array<size_type>
LinearAlgebra::Tensor<R>::compute_strides(const array<size_type>& sz) 
{
  size_type rnk=sz.size(); array<size_type> result(rnk);
  result[rnk-1]=1;
  for(size_type d=rnk; d!=0; --d) { 
    result[d-1]=sz[d]*result[d]; 
  }
  return result;
}


 
    
    
template<class R> inline
LinearAlgebra::SymmetricTensor<R>::SymmetricTensor(size_type n, size_type deg)
  : _arg_size(n), _degree(deg), _elements(Numeric::bin(n+deg-1,deg))
{ 
  //std::cerr << "SymmetricTensor(size_type m, size_type n, size_type deg)" << std::endl;
}
      

template<class R> inline
LinearAlgebra::SymmetricTensor<R>::SymmetricTensor(size_type n, size_type deg, const R* ptr) 
  : _arg_size(n), _degree(deg), _elements(Numeric::bin(n+deg-1,deg))
{
  _elements.fill(ptr); 
}
   
   
template<class R> template<class E> inline
LinearAlgebra::SymmetricTensor<R>::SymmetricTensor(const SymmetricTensorExpression<E>& t)
  : _arg_size(t().space_size()), _degree(t().degree()), 
    _elements(t().number_of_elements()) 
{ 
  MultiIndexIterator curr(this->_arg_size,this->_degree);
  MultiIndexIterator end(this->_arg_size,this->_degree+1);
  for( ; curr!=end; ++curr) {
    this->_elements[curr->position()]=t(*curr);
  }
}


template<class R> template<class E> inline
LinearAlgebra::SymmetricTensor<R>& 
LinearAlgebra::SymmetricTensor<R>::operator=(const SymmetricTensorExpression<E>& t)
{
  if(this!=&t) {
    const E& e=t();
    this->_arg_size=e.argument_size();
    this->_degree=e.degree();
    this->_elements().resize(this->number_of_elements());
    MultiIndexIterator curr(this->_arg_size,this->_degree);
    MultiIndexIterator end(this->_arg_size,this->_degree+1);
    for( ; curr!=end; ++curr) {
      this->_elements[curr->position()]=t(*curr);
    }
  }     
  return *this;
}        
      
template<class R> inline
size_type
LinearAlgebra::SymmetricTensor<R>::rank() const
{ 
  return this->degree();
}
      

template<class R> inline
size_type
LinearAlgebra::SymmetricTensor<R>::degree() const
{ 
  return this->_degree;
}


template<class R> inline
size_type
LinearAlgebra::SymmetricTensor<R>::number_of_elements() const
{ 
  return this->_arg_size^this->_degree;
}
        

template<class R> inline
size_type
LinearAlgebra::SymmetricTensor<R>::number_of_independent_elements() const
{ 
  return number_of_independent_elements(this->_arg_size,this->_degree);
}
        

template<class R> inline
size_type
LinearAlgebra::SymmetricTensor<R>::number_of_independent_elements(const size_type& n, const size_type& d) 
{
  return Numeric::bin(n+d-1,d); 
}
 
        

template<class R> inline
const size_type&
LinearAlgebra::SymmetricTensor<R>::argument_size() const
{ 
  return this->_arg_size;
}


template<class R> inline
size_type
LinearAlgebra::SymmetricTensor<R>::size(const size_type& i) const
{ 
  return this->_arg_size;
}
      

template<class R> inline
const R&
LinearAlgebra::SymmetricTensor<R>::operator() (const multi_index_type j) const
{ 
  return this->_elements[this->position(j)];
}
      

template<class R> inline
R&
LinearAlgebra::SymmetricTensor<R>::operator() (const multi_index_type j)
{ 
  return this->_elements[this->position(j)];
}
    

template<class R> inline
const R&
LinearAlgebra::SymmetricTensor<R>::operator() (const index_array_type i) const
{ 
  return this->_elements[this->position(i)];
}
      

template<class R> inline
R&
LinearAlgebra::SymmetricTensor<R>::operator() (const index_array_type i)
{ 
  return this->_elements[this->position(i)];
}
      

template<class R> inline
R*
LinearAlgebra::SymmetricTensor<R>::begin()
{ 
  return _elements.begin();
}
      

template<class R> inline
const R*
LinearAlgebra::SymmetricTensor<R>::begin() const
{ 
  return _elements.begin();
}
  

template<class R> inline
size_type 
LinearAlgebra::SymmetricTensor<R>::position(const multi_index_type& j) const
{
  if(j.number_of_variables()!=this->_arg_size) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  if(j.degree()!=this->_degree) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return j.position();
}
    
template<class R> inline
size_type 
LinearAlgebra::SymmetricTensor<R>::position(const index_array_type& i) const
{
  if(i.size()!=this->_degree) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return multi_index_type(this->_arg_size,i).position();
}
    




template<class R> inline
LinearAlgebra::DerivativeTensor<R>::DerivativeTensor(size_type m, size_type n, size_type deg)
  : _res_size(m), _arg_size(n), _degree(deg), _elements(number_of_independent_elements(m,n,deg))
{ 
  //std::cerr << "DerivativeTensor(size_type m, size_type n, size_type deg)" << std::endl;
}
      
template<class R> inline
LinearAlgebra::DerivativeTensor<R>::DerivativeTensor(size_type m, size_type n, size_type deg, const R* ptr) 
  : _res_size(m), _arg_size(n), _degree(deg), _elements(number_of_independent_elements(m,n,deg))
{ 
  _elements.fill(ptr); 
}
      
      
template<class R> template<class E> inline
LinearAlgebra::DerivativeTensor<R>::DerivativeTensor(const DerivativeTensorExpression<E>& t)
  : _res_size(t().result_size()), _arg_size(t().argument_size()), _degree(t().degree()), 
    _elements(t().number_of_independent_elements()) 
{ 
  this->assign(t());
}
  
    
template<class R> template<class E> inline
LinearAlgebra::DerivativeTensor<R>& 
LinearAlgebra::DerivativeTensor<R>::operator=(const DerivativeTensorExpression<E>& t) 
{
  if(this!=&t) {
    const E& e=t();
    this->_res_size=e.result_size();
    this->_arg_size=e.argument_size();
    this->_degree=e.degree();
    this->_elements().resize(this->number_of_independent_elements());
    this->assign(e);
  }     
  return *this;
}        
      

template<class R> inline
size_type
LinearAlgebra::DerivativeTensor<R>::rank() const
{ 
  return this->degree()+1;
}
      

template<class R> inline
size_type
LinearAlgebra::DerivativeTensor<R>::degree() const
{ 
  return this->_degree;
}
      

template<class R> inline
size_type
LinearAlgebra::DerivativeTensor<R>::number_of_elements() const
{ 
  return this->_res_size*this->_arg_size^this->_degree;
}
        

template<class R> inline
size_type
LinearAlgebra::DerivativeTensor<R>::number_of_independent_elements() const
{ 
  return number_of_independent_elements(this->_res_size,this->_arg_size,this->_degree);
}
        

template<class R> inline
const size_type&
LinearAlgebra::DerivativeTensor<R>::result_size() const
{ 
  return this->_res_size;
}


template<class R> inline
const size_type&
LinearAlgebra::DerivativeTensor<R>::argument_size() const
{ 
  return this->_arg_size;
}


template<class R> inline
size_type
LinearAlgebra::DerivativeTensor<R>::size(const size_type& i) const
{ 
  return i==0 ? this->_res_size : this->_arg_size;
}
      

template<class R> inline
const R&
LinearAlgebra::DerivativeTensor<R>::operator() (const size_type& i, const multi_index_type j) const
{
  if(j.degree()!=this->_degree) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  if(i>=this->_res_size) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return this->_elements[this->position(i,j)]; 
}
      

template<class R> inline
R&
LinearAlgebra::DerivativeTensor<R>::operator() (const size_type& i, const multi_index_type j)
{
  if(j.degree()!=this->degree()) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return this->_elements[this->position(i,j)]; 
}
    

template<class R> inline
R*
LinearAlgebra::DerivativeTensor<R>::begin()
{ 
  return _elements.begin();
}
      

template<class R> inline
const R*
LinearAlgebra::DerivativeTensor<R>::begin() const
{ 
  return _elements.begin();
}
      

template<class R> inline
size_type
LinearAlgebra::DerivativeTensor<R>::number_of_independent_elements(const size_type& nr, const size_type& na, const size_type& d)
{ 
  return nr*Numeric::bin(na+d-1,d);
}



 
template<class R> inline
size_type 
LinearAlgebra::DerivativeTensor<R>::position(const size_type& i, const multi_index_type& j) const
{
  if(i>=this->_res_size) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  if(j.number_of_variables()!=this->_arg_size) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  if(j.degree()!=this->degree()) { throw InvalidIndex(__PRETTY_FUNCTION__); }
  return i*bin(this->_arg_size+this->_degree-1,this->_degree)+j.position();
}


template<class R> template<class E> inline
void 
LinearAlgebra::DerivativeTensor<R>::assign(const E& e) 
{
  MultiIndexIterator k(this->argument_size(),this->degree()+1);
  for(size_type i=0; i!=this->result_size(); ++i) {
    for(MultiIndexIterator j(this->argument_size(),this->degree()); j!=k; ++j) {
      (*this)(i,j)=e(i,j);
    }
  }
}
    

template<class R> inline
std::ostream&
LinearAlgebra::operator<< (std::ostream& os, const Tensor<R>& T) 
{
  return T.write(os);
} 


template<class R> inline
std::ostream&
LinearAlgebra::operator<< (std::ostream& os, const SymmetricTensor<R>& T) 
{
  return T.write(os);
} 
    

template<class R> inline
std::ostream&
LinearAlgebra::operator<< (std::ostream& os, const DerivativeTensor<R>& T) 
{
  return T.write(os);
} 
    

} // namespace Ariadne

