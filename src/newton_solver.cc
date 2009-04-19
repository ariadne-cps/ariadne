/***************************************************************************
 *            newton_solver.cc
 *
 *  Copyright  2006-8  Pieter Collins
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
 
#include "newton_solver.h"

#include "logging.h"
#include "vector.h"
#include "matrix.h"
#include "function_interface.h"
#include "point.h"
#include "box.h"

namespace Ariadne {

Vector<Interval>
IntervalNewtonSolver::solve(const FunctionInterface& f, 
                            const Vector<Interval>& ix)
{
  const double& e=this->maximum_error();
  uint n=this->maximum_number_of_steps();
  ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
  Vector<Interval> x(ix);
  Vector<Interval> r(x);
  while(n>0) {
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<radius(x)<<"  x="<<x<<"\n");
    Vector<Float> m=midpoint(x);
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<Interval> im(m);
    Vector<Interval> w=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<w<<"\n");
    Matrix<Interval> A=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<A<<"\n");
    Matrix<Interval> Ainv=inverse(A);
    ARIADNE_LOG(5,"  inverse(Df(r))="<<Ainv<<"\n");
    Vector<Interval> dx=prod(Ainv , w);
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<Interval> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<Interval> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");

    if(subset(nr,r) && radius(nx) < e) {
      return nr;
    }
    if(disjoint(nr,r)) {
      throw EvaluationException("No result found -- disjoint");
    }
    r=intersection(nr,r);
    x=r;
    n=n-1;
  }
  throw EvaluationException("No result found -- maximum number of steps reached");
}

Vector<Interval>
IntervalNewtonSolver::fixed_point(const FunctionInterface& f, 
                                  const Vector<Interval>& ix)
{
  const double& e=this->maximum_error();
  uint n=this->maximum_number_of_steps();
  ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
  Vector<Interval> x(ix);
  Vector<Interval> r(x);
  while(n>0) {
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<radius(x)<<"  x="<<x<<"\n");
    Vector<Float> m=midpoint(x);
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<Interval> im(m);
    Vector<Interval> w=f.evaluate(im)-im;
    ARIADNE_LOG(5,"  f(m)="<<w<<"\n");
    Matrix<Interval> A=f.jacobian(x)-Matrix<Float>::identity(n);
    ARIADNE_LOG(5,"  Df(r)="<<A<<"\n");
    Matrix<Interval> Ainv=inverse(A);
    ARIADNE_LOG(5,"  inverse(Df(r))="<<Ainv<<"\n");
    Vector<Interval> dx=prod(Ainv , w);
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<Interval> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<Interval> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");

    if(subset(nr,r) && radius(nx) < e) {
      return nr;
    }
    if(disjoint(nr,r)) {
      throw EvaluationException("No result found -- disjoint");
    }
    r=intersection(nr,r);
    x=r;
    n=n-1;
  }
  throw EvaluationException("No result found -- maximum number of steps reached");
}

Vector<Interval>
KrawczykSolver::solve(const FunctionInterface& f,
                      const Vector<Interval>& ix)
{
  const double& e=this->maximum_error();
  uint n=this->maximum_number_of_steps();
  ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
  Vector<Interval> x(ix);
  Vector<Interval> r(x);
  Matrix<Float> I=Matrix<Float>::identity(ix.size());
  while(n>0) {
    ARIADNE_LOG(4,"Testing for root in "<<x<<"\n");
    ARIADNE_LOG(5,"  e="<<radius(x)<<"  x="<<x<<"\n");
    Vector<Float> m=midpoint(x);
    ARIADNE_LOG(5,"  m="<<m<<"\n");
    Vector<Interval> im(m);
    Vector<Interval> fm=f.evaluate(im);
    ARIADNE_LOG(5,"  f(m)="<<fm<<"\n");
    Matrix<Interval> J=f.jacobian(x);
    ARIADNE_LOG(5,"  Df(r)="<<J<<"\n");
    Matrix<Interval> M=inverse(midpoint(J));
    ARIADNE_LOG(5,"  inverse(Df(m))="<<M<<"\n");
    Vector<Interval> dx=prod(M,fm)-prod(Matrix<Interval>(I-prod(M,J)),Vector<Interval>(x-m));
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<Interval> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    Vector<Interval> nr(nx);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");

    if(subset(nr,r) && radius(nx) < e) {
      return nr;
    }
    if(disjoint(nr,r)) {
      throw EvaluationException("No result found -- disjoint");
    }
    r=intersection(nr,r);
    x=r;
    n=n-1;
  }
  throw EvaluationException("No result found -- maximum number of steps reached");
}

Vector<Interval>
KrawczykSolver::fixed_point(const FunctionInterface& f,
                            const Vector<Interval>& ix)
{
    ARIADNE_NOT_IMPLEMENTED;
}
} // namespace Ariadne
