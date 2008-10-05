#include <iostream>
#include <vector>
#include <algorithm>

#include "zonotope.h"

#include "array.h"
#include "vector.h"
#include "matrix.h"
#include "point.h"
#include "box.h"
#include "list_set.h"

//#include "linear_program.h"
//#include "affine_model.h"


namespace Ariadne {

template<class X> class LinearProgram {
 public:
  LinearProgram(const Matrix<X>& A) { assert(false); }
  bool is_feasible() { assert(false); }
};


template<class BS>
ListSet<BS> 
subdivide(const BS& bs, const Float& r)
{
  ListSet<BS> result;
  ListSet<BS> working(bs);
  while(!working.size()==0) {
    BS set=working.pop();
    if(set.radius()<r) {
      result.adjoin(set);
    } else {
      working.adjoin(split(set));
    }
  }
  return result;
}



void 
accumulate(Float& value, Float& error, uint n, const Float* aptr, const Float* bptr) 
{
  Interval v=value;
  for(uint i=0; i!=n; ++i) {
    v+=aptr[i]*bptr[i];
  }
  value=v.midpoint();
  error=add_up(error,v.radius());
}

Vector<Float>
row_norms(const Matrix<Interval>& A)
{
  uint const& m=A.row_size();
  uint const& n=A.column_size();
  Vector<Float> e(m);
  for(uint i=0; i!=m; ++i) {
    for(uint j=0; j!=n; ++j) {
      e[i]=add_up(e[i],mag(A(i,j)));
    }
  }
  return e;
}

Vector<Float>
row_errors(const Matrix<Interval>& A)
{
  uint const& m=A.row_size();
  uint const& n=A.column_size();
  Vector<Float> e(m);
  for(uint i=0; i!=m; ++i) {
    for(uint j=0; j!=n; ++j) {
      e[i]=add_up(e[i],A(i,j).radius());
    }
  }
  return e;
}

Vector<Float>
errors(const Vector<Interval>& pt)
{  
  Vector<Float> result(pt.size());
  for(uint i=0; i!=pt.size(); ++i) {
    result[i]=pt[i].radius();
  }
  return result;
}


Vector<Float>
row_errors(const Vector<Interval>& pt, const Matrix<Interval>& A)
{
  assert(pt.size()==A.row_size());
  Vector<Float> result(pt.size());
  for(uint i=0; i!=A.row_size(); ++i) {
    result[i]=pt[i].radius();
    for(uint j=0; j!=A.column_size(); ++j) {
      result[i]=add_up(result[i],A(i,j).radius());
    }
  }
  return result;
}
  
Vector<Float>
add_up(const Vector<Float>& v1, const Vector<Float>& v2) 
{
  Vector<Float> result;
  for(uint i=0; i!=v1.size(); ++i) {
    result[i]=add_up(v1[i],v2[i]);
  }
  return result;
}

tribool 
norm_grtr(const Vector<Float>& v1, const Vector<Float>& v2) 
{
  return norm(v1)>norm(v2);
}






 
Zonotope::Zonotope()
  : _centre(), _generators(), _error()
{
}

 
Zonotope::Zonotope(uint d)
  : _centre(d), _generators(d,0), _error(d)
{
}

 
Zonotope::Zonotope(uint d, uint m)
  : _centre(d), _generators(d,m), _error(d)
{
}

 
Zonotope::Zonotope(const Vector<Float>& c, const Matrix<Float>& G, const Vector<Float>& e)
  : _centre(c), _generators(G), _error(e)
{
  assert(c.size()==G.row_size());
  assert(c.size()==e.size());
}

 
Zonotope::Zonotope(const Vector<Float>& c, const Matrix<Float>& G)
  : _centre(c), _generators(G), _error(c.size())
{
  assert(c.size()==G.row_size());
}

 
Zonotope::Zonotope(const Vector<Interval>& c, const Matrix<Float>& G)
  : _centre(midpoint(c)), _generators(G), _error(errors(c))
{
  assert(c.size()==G.row_size());
}

 
Zonotope::Zonotope(const Vector<Float>& c, const Matrix<Interval>& G)
  : _centre(c), _generators(midpoint(G)), _error(row_errors(G))
{
  assert(c.size()==G.row_size());
}

 
Zonotope::Zonotope(const Vector<Interval>& c, const Matrix<Interval>& G)
  : _centre(midpoint(c)), _generators(midpoint(G)), _error(row_errors(c,G))
{
  assert(c.size()==G.row_size());
}


       
Zonotope::Zonotope(const Zonotope& z) 
  : _centre(z._centre), _generators(z._generators), _error(z._error)
{
}

       
Zonotope&
Zonotope::operator=(const Zonotope& z) 
{ 
  if(this!=&z) {
    this->_centre=z._centre;
    this->_generators=z._generators;
    this->_error=z._error;
  }
  return *this;
}

       
uint
Zonotope::dimension() const
{
  return this->_centre.size();
}

       
uint
Zonotope::number_of_generators() const
{
  return this->_generators.column_size();
}

       
const Vector<Float>&
Zonotope::centre() const
{
  return this->_centre;
}

       
const Matrix<Float>&
Zonotope::generators() const
{
  return this->_generators;
}

       
const Vector<Float>&
Zonotope::error() const
{
  return this->_error;
}

       
Vector<Interval>
Zonotope::domain() const
{
  return Vector<Interval>(this->number_of_generators(),Interval(-1,1));
}

       
Vector<Interval>
Zonotope::bounding_box() const
{
  const Zonotope& z=*this;
  Vector<Interval> b=z.centre()+Matrix<Interval>(z.generators())*z.domain()+z.error()*Interval(-1,1);
  return b;
}

       
Float
Zonotope::radius() const
{
  return Ariadne::radius(this->bounding_box());
}

       
tribool
Zonotope::contains(const Point& pt) const
{
  return Ariadne::contains(*this,pt);
}




       
tribool
empty(const Zonotope& z) 
{
  return false;
}

       
tribool
bounded(const Zonotope& z) 
{
  return true;
}

       

       
Float 
radius(const Zonotope& z) 
{
  return Ariadne::radius(z.centre()+z.generators()*z.domain()+z.error()*Interval(-1,1));
}


       







Box
bounding_box(const Zonotope& z)
{
  return z.bounding_box();
}



ListSet< Zonotope >
split(const Zonotope& z)
{
  // FIXME: Not quite guarenteed to give an over-approximation
  typedef Interval I;
  
  
  ListSet< Zonotope  > result;
  
  uint d=z.dimension();
  uint m=z.number_of_generators();
  Vector<Float> const& c=z.centre();
  Matrix<Float> const& G=z.generators();
  Vector<Float> const& e=z.error();
  
  array<Float> norms(m,0);
  for(uint j=0; j!=m; ++j) {
    norms[j]=norm(Vector<Float>(column(G,j)));
  }

  Float max_norm=0;
  uint longest_generator=0;
  for(uint j=0; j<m; ++j) {
    if(norms[j]>max_norm) {
      max_norm=norms[j];
      longest_generator=j;
    }
  }
  for(uint k=0; k<d; ++k) {
    if(e[k]>max_norm) {
      max_norm=e[k];
      longest_generator=m+k;
    }
  }
  
  if(longest_generator<m) {
    Matrix<Float> new_generators=z.generators();
    uint j=longest_generator;
    for(uint i=0; i!=d; ++i) {
      new_generators(i,j)=div_up(new_generators(i,j),2);
    }
    
    Vector<Float> v=column(new_generators,j);
    Vector<Float> new_centre=sub_approx(c,v);
    result.adjoin(Zonotope(new_centre,new_generators,e));
    new_centre=add_approx(c,v);
    result.adjoin(Zonotope(new_centre,new_generators,e));
 } else {
    uint k=longest_generator-m;
    Vector<Float> new_centre = z.centre();
    const Matrix<Float>& new_generators = z.generators();
    Vector<Float> new_error=e;
    new_error[k]=div_up(new_error[k],2);
    new_centre[k]=add_approx(z.centre()[k],new_error[k]);
    result.adjoin(Zonotope(new_centre,new_generators,new_error));
    new_centre[k]=sub_approx(z.centre()[k],new_error[k]);
    result.adjoin(Zonotope(new_centre,new_generators,new_error));
  }
  return result;
} 




       
Zonotope::Zonotope(const Box& r) 
  : _centre(r.size()), _generators(r.size(),r.size()), _error(r.size())
{
  uint d=r.size();
  Vector<Float>& c=this->_centre;
  Matrix<Float>& G=this->_generators;
  Vector<Float>& e=this->_error;
  for(uint i=0; i!=d; ++i) {
    c[i]=med_approx(r[i].lower(),r[i].upper());
    for(uint j=0; j!=d; ++j) {
      G(i,j)=0;
    }
    G(i,i)=rad_up(r[i].lower(),r[i].upper());
    e(i)=0;
  }
}



/*
Zonotope
apply(const AffineModel& am,
      const Zonotope& z)
{
  
  typedef Interval I;
  
  assert(z.centre()==am.centre());
  if(!subset(z,am.domain())) {
    std::cerr<<"z="<<z<<"\nz.bounding_box()="<<z.bounding_box()<<"\nam.domain()="<<am.domain()<<std::endl;
  }
  assert(possibly(subset(z,am.domain())));
  
  uint d=z.size();
  uint m=z.number_of_generators();
  uint nd=am.result_size();

  //Vector<Float> const& c=z.centre();
  Matrix<Float> const& G=z.generators();
  Vector<Float> const& e=z.error();

  Vector<Interval> const& nic=am.value();
  Matrix<Interval> const& iDf=am.jacobian();
  Matrix<Interval> niG=iDf*G;
  
  Vector<Float> nc=midpoint(nic);
  Matrix<Float> nG=midpoint(niG);
  Vector<Float> ne(nd);

  for(uint i=0; i!=nd; ++i) {
    R& err=ne[i];
    err=add_up(err,nic[i].radius());
    for(uint j=0; j!=m; ++j) {
      err=add_up(err,niG(i,j).radius());
    }
    for(uint k=0; k!=d; ++k) {
      err=add_up(err,mul_up(mag(iDf(i,k)),e[k]));
    }
  }
  
  return Zonotope(nc,nG,ne);

}

*/


 
Zonotope 
approximation(const Zonotope& z)
{
  return Zonotope(z.centre(),z.generators());
}

 
Zonotope
over_approximation(const Zonotope& z)
{
  return z;
}

 
Zonotope
error_free_over_approximation(const Zonotope& z)
{
  uint d=z.dimension();
  uint m=z.number_of_generators();
  Matrix<Float> nG(d,m+d);
  project(nG,range(0,d),range(0,m))=z.generators();
  for(uint i=0; i!=d; ++i) {
    nG(i,i+m)=z.error()[i];
  }
  return Zonotope(z.centre(),nG);
}

 
/*
Zonotope
nonsingular_over_approximation(const Zonotope& z)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}  
*/


 
 /*
Zonotope
orthogonal_over_approximation(const Zonotope& z)
{
  //assert(iz.size()==iz.number_of_generators());
  typedef Interval I;
  Zonotope ez=error_free_over_approximation(z);

  const Vector<Float>& c=ez.centre();
  const Matrix<Float>& G=ez.generators();
  
  Matrix<Float> aQ,aR;
  make_lpair(aQ,aR)=qr_approx(G);

  Matrix<Interval> aQinv=inverse(aQ);
  Matrix<Interval> iR=aQinv*G;
  DiagonalMatrix<Float> aD(::row_norms(iR));

  Matrix<Interval> niG=aQ*aD;

  return Zonotope(c,niG);
}  
 */


Zonotope
cascade_over_approximation(const Zonotope& z, uint cs)
{
  using namespace std;
  
  if(z.number_of_generators()<=z.dimension()*cs) { return z; }  

  assert(z.number_of_generators()%z.dimension()==0);

  uint d=z.dimension();
  uint nb=z.number_of_generators()/z.dimension(); // number of generator blocks
   

  const Matrix<Float>& G=z.generators();
  array<Float> norms(nb);
  for(uint i=0; i!=nb; ++i) {
    norms[i]=Ariadne::norm(Matrix<Float>(project(G,range(0,d),range(i*d,(i+1)*d))));
  }
  
  // Compute the new number of blocks
  uint nnb=cs;
  Float sum=0;
  for(uint i=nb-1; i!=0; --i) {
    sum=add_approx(sum,norms[i]);
    if(sum>norms[i-1]) {
      nnb=i;
    }
  }
  nnb=min(nnb,cs);
  // Reduce generators
  Matrix<Float> rG(d,d*nnb);
  project(rG,range(0,d),range(0,d*(nnb-1)))=project(G,range(0,d),range(0,d*(nnb-1)));
  for(uint i=0; i!=d; ++i) {
    Float& err=rG(i,d*(nnb-1)+i);
    for(uint j=d*(nnb-1); j!=G.column_size(); ++j) {
      err=add_up(err,std::abs(G(i,j)));
    }
  }
  return Zonotope(z.centre(),rG);
}



/*
 
Zonotope<Interval,R> 
orthogonal_over_approximation(const Zonotope<R,R>& z)
{
  // FIXME: Subdivide in zero order as well!
  static bool warn=true;
  if(warn) {
    std::cerr << std::endl << "WARNING: orthogonal_over_approximation(Zonotope<I,R>) does not over-approximate roundoff errors." << std::endl;
    warn=false;
  }
  Zonotope<R,R> oaz=over_approximation(z);
  
  QRMatrix< Interval > QR(oaz.generators());
  Point< Interval > c(oaz.centre());
  Matrix<Float> G(z.size(),z.number_of_generators());

  Matrix< Interval > q=QR.Q();
  Matrix< Interval > r=QR.R();
  for(uint i=0; i!=z.size();++i) {
    Interval a=0;
    for(uint j=i; j!=z.number_of_generators(); ++j) {
      a+=r(i,j);
    }
    for(uint k=0; k!=z.size(); ++k) {
      Interval b=q(k,i)*a;
      G(k,i)=b.midpoint();
      c[k]+=(b-b.midpoint());
    }
  }
  return Zonotope<R,R>(midpoint(c),G);
}

 
Zonotope<Interval,R> 
orthogonal_over_approximation(const Zonotope<Interval,R>& z)
{
  Zonotope<R,R> oaz=over_approximation(z);
  
  QRMatrix< Interval > QR(oaz.generators());
  Point< Interval > c(oaz.centre());
  Matrix<Float> G(z.size(),z.number_of_generators());

  Matrix< Interval > q=QR.Q();
  Matrix< Interval > r=QR.R();
  for(uint i=0; i!=z.size();++i) {
    Interval a=0;
    for(uint j=i; j!=z.number_of_generators(); ++j) {
      a+=r(i,j);
    }
    for(uint k=0; k!=z.size(); ++k) {
      Interval b=q(k,i)*a;
      G(k,i)=b.midpoint();
      c[k]+=(b-b.midpoint());
    }
  }
  return Zonotope<Interval,R>(c,G);
}

 
Zonotope< Interval > 
orthogonal_over_approximation(const Zonotope< Interval >& z)
{
  Zonotope<R,R> oaz=over_approximation(z);
  
  QRMatrix< Interval > QR(oaz.generators());
  Point< Interval > c(oaz.centre());
  Matrix<Float> G(z.size(),z.number_of_generators());

  Matrix< Interval > q=QR.Q();
  Matrix< Interval > r=QR.R();
  for(uint i=0; i!=z.size();++i) {
    Interval a=0;
    for(uint j=i; j!=z.number_of_generators(); ++j) {
      a+=r(i,j);
    }
    for(uint k=0; k!=z.size(); ++k) {
      Interval b=q(k,i)*a;
      G(k,i)=b.midpoint();
      c[k]+=(b-b.midpoint());
    }
  }
  return Zonotope< Interval >(c,G);
}
*/






std::ostream&
operator<<(std::ostream& os, const Zonotope& z) 
{
  os << "["<<z.centre();
  for(uint j=0; j!=z.number_of_generators(); ++j) {
    for(uint i=0; i!=z.dimension(); ++i) {
      os << (i==0 ? ';' : ',') << z.generators()[i][j];
    }
  }
  os << "]";
  return os;
}




std::istream& 
operator>>(std::istream& is, Zonotope& z)
{
  Vector<Float> centre;
  Matrix<Float> generators;
  char c0,c1,c2;
  is >> c0 >> centre >> c1 >> generators >> c2;
  z = Zonotope(centre,generators);
  return is;
}





// Test vertices individually. Highly inefficient!! 



/* Set up linear program to solve 
 *   \f[x=c+Ge;\ l\leq x\leq u;\ -1\leq e\leq1\f].
 *
 * Change variables to normalize \f$x\f$ and \f$e\f$
 *   \f[x'=x-l,\ e'=e+1;   x'-Ge' = c-G1-l;  0\leq x\leq u-l; \ 0\leq e\leq 2.\f] 
 * 
 * Introduce slack variables sx and se, and artificial variables ax. Problem in standard form
 *   \f[ \begin{matrix}I&0\\0&I\\\pm I&\mp G\end{matrix} \begin{matrix}x'\\e'\end{matrix}
 *        + \begin{matrix}I&&\\&I&\\&&I\end{matrix}\begin{matrix}sx\\se\\ax\end{matrix}
 *             = \begin{matrix}u-l\\2\\\pm(c-G1-l)\end{matrix} \f]
 * 
 */
tribool
disjoint(const Zonotope& z, const Box& r)
{
  assert(z.dimension()==r.size());
  uint d=z.dimension();
  uint m=z.number_of_generators();
  
  // Construct tableau for testing intersection of zonotope and rectangle
  // Box  l<=x<=u
  // Zonotope  x==c+Ge,  -1<=e<=1
  // 
  // Translate x'=x-l,  e'=e+1
  //   0<=x'<=u-l      ->  x' +     + sx'               == u-l
  //   0<=e'<=2        ->     +  e' +     + se'         == 2
  //   x'+l==c+G(e'-1) ->  x' + Ge'             +/- ax' == c-l-G1
  //  
  // Change sign of RHS of first equality if necessary
  // Introduce slack variables for last two inequalities
  Matrix<Float> T(2*d+m+1,d+m+1);
  
  const Vector<Float> l=Ariadne::lower(r)-z.error();
  const Vector<Float> u=Ariadne::upper(r)+z.error();
  const Vector<Float>& c=z.centre();
  const Matrix<Float>& G=z.generators();
  //const Vector<Float>& e=z.error();
  
  const Vector<Float> qo(m,1.0);
  const Vector<Float> ql=l;
  const Vector<Float> qu=u;
  const Vector<Float> qd=qu-ql;
  const Vector<Float> qc=c;
  const Matrix<Float> qG=G;
  const Vector<Float> qrhs=qc-ql-qG*qo;
  
  { std::clog << "ql=" << ql << ", qd=" << qd <<", qc=" << qc << ", qrhs=" << qrhs << std::endl; }
  
  // Set up constraints x+sx=u-l
  for(uint i=0; i!=d; ++i) {
    T(i,i)=1;
    T(i,d+m)=qu(i)-ql(i);
  }
  
  // Set up constraints e+se=2
  for(uint j=0; j!=m; ++j) {
    T(d+j,d+j)=1;
    T(d+j,d+m)=2;
  }
  
  // Set up constraints x-Ge \pm ax=c-l-G1 
  for(uint i=0; i!=d; ++i) {
    if(qrhs(i)>=0.0) {
      T(i+d+m,i)=1;
      for(uint j=0; j!=m; ++j) {
        T(i+d+m,d+j)=-qG(i,j);
      }
      T(i+d+m,d+m)=qrhs(i);
    }
    else {
      T(i+d+m,i)=-1;
      for(uint j=0; j!=m; ++j) {
        T(i+d+m,d+j)=qG(i,j);
      }
      T(i+d+m,d+m)=-qrhs(i);
    }
  } 
  
  // Set up cost function ax^T1
  for(uint i=0; i!=d; ++i) {
    T(2*d+m,i) -= T(i+d+m,i);
    for(uint j=0; j!=m; ++j) {
      T(2*d+m,d+j) -= T(i+d+m,d+j);
    }
    T(2*d+m,d+m) -= T(i+d+m,d+m);
  }
  
  LinearProgram<Float> lp(T);
  tribool result=!lp.is_feasible();
  
  return result;
}



tribool
disjoint(const Zonotope& z1, const Zonotope& z2)
{
  assert(z1.dimension()==z2.dimension());
  
  uint d=z1.dimension();
  const Float zero=0;
  const Float one=1;
  uint m1=z1.number_of_generators();
  uint m2=z2.number_of_generators();
  
  Matrix<Float> T(m1+m2+d+1,m1+m2+1);
  
  const Vector<Float> qo1(m1,one);
  const Vector<Float> qo2(m2,one);
  
  const Vector<Float>& qc1=z1.centre();
  const Matrix<Float>& qG1=z1.generators();
  const Vector<Float>& qc2=z2.centre();
  const Matrix<Float>& qG2=z2.generators();
  Vector<Float> qrhs = qG1*qo1 - qG2*qo2 + (qc2 - qc1);
  
  // Set up constraints e1 + se1 = 2
  for(uint j1=0; j1!=m1; ++j1) {
    T(j1,j1)=1;
    T(j1,m1+m2)=2;
  }
  
  // Set up constraints e2 + se2 = 2
  for(uint j2=0; j2!=m2; ++j2) {
    T(m1+j2,m1+j2)=1;
    T(m1+j2,m1+m2)=2;
  }
  
  // Set up constraints G1*e1 - G2*e2 = (c2 - G2*1) - (c1 - G1*1)
  for(uint i=0; i!=d; ++i) {
    if(qrhs(i)>=zero) {
      for(uint j1=0; j1!=m1; ++j1) {
        T(m1+m2+i,j1)=qG1(i,j1);
      }
      for(uint j2=0; j2!=m2; ++j2) {
        T(m1+m2+i,m1+j2)=qG2(i,j2);
      }
      T(m1+m2+i,m1+m2)=qrhs(i);
    }
    else {
      for(uint j1=0; j1!=m1; ++j1) {
        T(m1+m2+i,j1)=-qG1(i,j1);
      }
      for(uint j2=0; j2!=m2; ++j2) {
        T(m1+m2+i,m1+j2)=-qG2(i,j2);
      }
      T(m1+m2+i,m1+m2)=-qrhs(i);
    }
  } 
  
  // Set up cost function ax^T1
  for(uint i=0; i!=d; ++i) {
    for(uint j1=0; j1!=m1; ++j1) {
      T(m1+m2+d,j1) -= T(m1+m2+i,j1);
    }
    for(uint j2=0; j2!=m2; ++j2) {
      T(m1+m2+d,m1+j2) -= T(m1+m2+i,m1+j2);
    }
    T(m1+m2+d,m1+m2) -= T(m1+m2+i,m1+m2);
  }
  
  LinearProgram<Float> lp(T);
  
  tribool result=!lp.is_feasible();
  
  //std::clog << "disjoint(" << z1 << "," << z2 << ")=" << result << std::endl;
  return result;
}



/* Set up LP problem to solve \f$c+Ge=p\f$; \f$-1<=e<=1\f$.
 * Change variables so that the problem becomes \f$Ge=p-c-G1;\ 0\leq e\leq2\f$.
 * Change sign of \f$ Ge=p-c-G1\f$ to make right-hand side positive.
 */
tribool 
contains(const Zonotope& z, const Point& pt)
{ 
  //std::clog << "Zonotope::contains(const Vector<Float>& )" << std::endl;
  assert(z.dimension()==pt.dimension());
  uint d=z.dimension();
  uint m=z.number_of_generators();
  
  Float zero=0;
  Float one=1;
  Float two=2;
  
  const Vector<Float>& qc=z.centre();
  const Vector<Float>& qp=pt;
  const Matrix<Float>& qG=z.generators();
  const Vector<Float> qo(m,one);
  const Vector<Float> zv(m,zero);
  const Vector<Float> tv(m,two);
  
  Vector<Float> qrhs=qp-qc+qG*qo;
  
  Matrix<Float> T(d+m+1,m+1);
  
  // Set up constraints e+se=2
  for(uint j=0; j!=m; ++j) {
    T(j,j)=1;
    T(j,m)=2;
  }
  
  // Set up constraints Ge \pm ax = p-c+G1
  for(uint i=0; i!=d; ++i) {
    if(qrhs(i)>=zero) {
      for(uint j=0; j!=m; ++j) {
        T(m+i,j)=qG(i,j); 
      }
      T(m+i,m)=qrhs(i);
    } else {
      for(uint j=0; j!=m; ++j) {
        T(m+i,j)=-qG(i,j); 
      }
      T(m+i,m)=-qrhs(i);
    }
  }
  
  // Set up cost function ax^T 1
  for(uint i=0; i!=d; ++i) {
    for(uint j=0; j!=m; ++j) {
      T(m+d,j)-=T(m+i,j);
    }
    T(m+d,m)-=T(m+i,m);
  }
  
  LinearProgram<Float> lp(T);
  //std::clog << lp.tableau() << std::endl;
  tribool result=lp.is_feasible();
  //std::clog << lp.tableau() << std::endl;
  return result;
}

} // namespace Ariadne

