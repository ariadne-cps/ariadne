/***************************************************************************
 *            linear_programming.cc
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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
 
#include "config.h"

#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "linear_programming.h"


namespace Ariadne {

enum VariableType { LOWER, UPPER, BASIS, INFEASIBLE };

std::ostream& operator<<(std::ostream& os, VariableType t) {
    return os << (t==INFEASIBLE ? 'I' : t==BASIS ? 'B' : t==LOWER ? 'L' : 'U');  
}

template<class X>
tribool feasible(const Matrix<X>& A, const Vector<X>& b, const Vector<X>& l, const Vector<X>& u)
{
    const size_t m=A.row_size();
    const size_t n=A.column_size();
    ARIADNE_ASSERT(b.size()==m);
    ARIADNE_ASSERT(l.size()==n);
    ARIADNE_ASSERT(u.size()==n);
    Vector<X> c(n);

    array<int> p(m); // Array of basis indices
    array<int> q(n-m); // Array of non-basis indices
    Matrix<X> B(m,m);
    Vector<X> x(n);
    Vector<X> y(m);
    Vector<X> z(n-m);
    
    array<VariableType> t(n);

    // TODO: Find basis properly

    std::cerr<<"\nfeasible(A,b,l,u) with\nA="<<A<<" b="<<b<<" l="<<l<<" u="<<u<<"\n"<<std::endl;

    // Assume that last m columns of A form identity matrix
    for(size_t i=0; i!=m; ++i) {
        assert(A[i][i+(n-m)]==1);
        for(size_t j=0; j!=m; ++j) {
            if(i!=j) {
                assert(A[i][j+(n-m)]==0);
            }
        }
    }
    
    // Set initial basis and inverse matrix
    for(size_t i=0; i!=m; ++i) { p[i]=i+(n-m);  }
    for(size_t i=0; i!=n-m; ++i) { q[i]=i;  }
    B=Matrix<X>::identity(m);
    
    std::cerr<<"  p="<<p<<" q="<<q<<"\n  B="<<B<<std::endl;
    for(uint i=0; i!=n; ++i) { 
        x[i]=l[i];
    }

    for(uint i=0; i!=n-m; ++i) { t[i]=LOWER; }
    for(uint i=0; i!=m; ++i) { t[i+(n-m)]=INFEASIBLE; }

     {
        std::cerr << "Setting up LP problem\n";
        // Compute the vector w = b - A_N x_N
        Vector<X> w(m);
        for(uint i=0; i!=m; ++i) {
            w[i]=b[i];
            for(uint j=0; j!=n-m; ++j) {
                uint k=q[j];
                w[i]-= A[i][k] * x[k];
            }
        }
        
        // Compute x_B = B w
        for(uint i=0; i!=m; ++i) {
            uint k=p[i];
            x[k]=0;
            for(uint j=0; j!=m; ++j) {
                x[k]+=B[i][j]*w[j];
            }
        }
        
        Matrix<X> A_B(m,m);
        for(size_t i=0; i!=m; ++i) {
            for(size_t j=0; j!=m; ++j) {
                A_B[i][j]=A[i][p[j]];
            }
        }
        assert(norm(Matrix<X>(prod(A_B,B)-Matrix<X>::identity(m)))<0.03125);
     }
     
     for(uint v=0; v!=m; ++v) {
         std::cerr << "  t="<<t<<" p="<<p<<" q="<<q<<" B="<<B<<"\n  x="<<x<<" Ax="<<prod(A,x)<<"\n";
         size_t k=p[v];
         std::cerr << "Trying to put variable "<<k<<" into constraint set, solving constraint"<<v<<"\n";
         Vector<X> c(n);
         if(x[k]<l[k]) { c[k]=-1; } else if(x[k]>u[k]) { c[k]=+1; } else { c[k]=0; }
         std::cerr << "  A="<<A<<" b="<<b<<" c="<<c<<" t="<<t<<" x="<<x<<std::endl;
         while(x[k]<l[k] || x[k]>u[k]) {
             std::cerr << "Starting LP step\n";
             // Try to decrease c[k];
             // Compute dual variables y=B^Tc_B
             for(uint i=0; i!=m; ++i) {
                 y[i]=0;
                 for(uint j=0; j!=m; ++j) {
                     y[i]+=B[j][i]*c[p[j]];
                 }
             }
             Vector<X> c_B(m); for(uint i=0; i!=m; ++i) { c_B[i]=c[p[i]]; }
             std::cout << "      y="<<y<<"="<<prod(c_B,B)<<"\n";
             
             // Compute reduced costs z_N=c_N-A_N^Ty
             for(uint i=0; i!=n-m; ++i) {
                 z[i]=c[q[i]];
                 for(uint j=0; j!=m; ++j) {
                     z[i]+=A[j][q[i]]*y[j];
                 }
             }
             std::cerr << "    x="<<x<<" y="<<y<<" z="<<z<<std::endl;
             
             size_t pi=n-m;
             for(size_t i=0; i!=n-m; ++i) {
                 if(z[i]<0) { pi=i; }
             }
             
             if(pi==(n-m)) {
                 std::cerr<<"  can't add constraint "<<k<<": "<<l[k]<<"<=x["<<k<<"]<="<<u[k]<<"; currently x["<<k<<"]="<<x[k]<<std::endl;
                 // Compute certificate of infeasibility
                 if(x[k]<l[k]) { t[k]=LOWER; } else if(x[k]>u[k]) { t[k]=UPPER; } else { t[k]=BASIS; }
                 Vector<X> x_N(n);
                 for(uint i=0; i!=n; ++i) {
                     if(t[i]==LOWER) { x_N[i]=l[i]; }  else if(t[i]==UPPER) { x_N[i]=u[i]; } else { x_N[i]=0; }
                 }
                 Vector<X> w(m);
                 w=b-prod(A,x_N);
                 Vector<X> x_B=prod(B,w);
                 Vector<X> c_B(m); for(uint i=0; i!=m; ++i) { c_B[i]=c[p[i]]; }
                 Vector<X> y=prod(c_B,B);
                 Vector<X> yA=prod(y,A);
                 X d=dot(y,w);
                 std::cerr<<"  x="<<x<<" w="<<w<<" y="<<y<<" y^Tw="<<d<<" t="<<t<<" yA="<<yA<<"\n";
                 assert(d>0);
                 for(uint i=0; i!=n; ++i) {
                     if(t[i]==LOWER) { if(yA[i]>=0) { std::cerr<<"  t["<<i<<"]=L, yA["<<i<<"]="<<yA[i]<<std::endl; } assert(yA[i]<0); }
                     if(t[i]==UPPER) { if(yA[i]<=0) { std::cerr<<"  t["<<i<<"]=U, yA["<<i<<"]="<<yA[i]<<std::endl; } assert(yA[i]>0); }
                 }
                 std::cerr << "  t="<<t<<"  y="<<y<<std::endl;
                 return false; 
             }
         }
         t[k]=BASIS;
     }
         
     // Check solution
     for(size_t i=0; i!=n; ++i) {
         assert(x[i]>=l[i]);
         assert(x[i]<=u[i]);
     }
     Vector<X> Ax=prod(A,x);
     for(size_t i=0; i!=m; ++i) {
         assert(Ax[i]==b[i]);
     }

     return true;
}



template tribool feasible<Rational>(const Matrix<Rational>&, const Vector<Rational>&, const Vector<Rational>&, const Vector<Rational>&);


} // namespace Ariadne

