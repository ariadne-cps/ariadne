/***************************************************************************
 *            solver.cc
 *
 *  Copyright  2006-9  Pieter Collins
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
 
#include "solver.h"

#include "logging.h"
#include "vector.h"
#include "matrix.h"
#include "function_interface.h"

namespace Ariadne {

namespace {

void
solve_all(Set< Vector<Interval> >& r,
          const SolverInterface& s,
          const FunctionInterface& f,
          const Vector<Interval>& ix)
{
    // Test for no solution
    if(disjoint(f.evaluate(ix),ix)) {
        //std::cerr<<"No solution in"<<ix<<"\n";
        return;
    }

    // If radius is too small, assume solution is not verified
    if(radius(ix)<s.maximum_error()) {
        std::cerr<<"Warning: Cannot verify solution in "<<ix<<"\n";
        return;
    }

    bool invertible_jacobian=true;
    try {
        Matrix<Interval> Jinv=inverse(f.jacobian(ix));
    }
    catch(const SingularMatrixException& e) {
        invertible_jacobian=false;
    }

    if(invertible_jacobian) {
        //std::cerr<<"Nonsingular matrix -- applying contractor\n";
        try {
            r.insert(s.solve(f,ix));
        }
        catch(const EvaluationException& e) {
            std::cerr<<"Evaluation exception -- No solution in"<<ix<<"\n";
            // No solution
        }
    } else {
        //std::cerr<<"Singular matrix over "<<ix<<" -- splitting\n";
        std::pair< Vector<Interval>, Vector<Interval> > splt=split(ix);
        solve_all(r,s,f,splt.first);
        solve_all(r,s,f,splt.second);
    }

}

} // namespace


SolverBase::SolverBase(double max_error, uint max_steps)
  : _max_error(max_error), _max_steps(max_steps) 
{
}


Set< Vector<Interval> >
SolverBase::solve_all(const FunctionInterface& f,
                      const Vector<Interval>& ix) const
{
    Set< Vector<Interval> > r;
    Ariadne::solve_all(r,*this,f,ix);
    return r;
}



Vector<Interval>
SolverBase::solve(const FunctionInterface& f,
                  const Vector<Interval>& ix) const
{
  const double& e=this->maximum_error();
  uint n=this->maximum_number_of_steps();
  ARIADNE_LOG(1,"verbosity="<<verbosity<<"\n");
  Vector<Interval> r(ix);
  bool has_solution=false;
  while(n>0) {
    Vector<Interval> nr=this->step(f,r);
    ARIADNE_LOG(5,"  nr="<<nr<<"\n");

    if(!has_solution && subset(nr,r)) {
        has_solution=true;
    }

    if(has_solution && radius(nr) < e) {
      return nr;
    }

    if(disjoint(nr,r)) {
      throw EvaluationException("No result found -- disjoint");
    }
    r=intersection(nr,r);
    n=n-1;
  }
  throw EvaluationException("No result found -- maximum number of steps reached");
}



Vector<Interval>
SolverBase::fixed_point(const FunctionInterface& f, const Vector<Interval>& pt) const
{
  return Vector<Interval>(this->solve(DifferenceFunction(f),pt)); 
}




Vector<Interval>
IntervalNewtonSolver::step(const FunctionInterface& f,
                           const Vector<Interval>& x) const
{
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
    Vector<Interval> dx=prod(Ainv, w);
    ARIADNE_LOG(5,"  dx="<<dx<<"\n");
    Vector<Interval> nx= m - dx;
    ARIADNE_LOG(5,"  nx="<<nx<<"\n");
    return nx;
}

Vector<Interval>
KrawczykSolver::step(const FunctionInterface& f,
                     const Vector<Interval>& x) const
{
    Matrix<Interval> I=Matrix<Interval>::identity(x.size());
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
    return nr;
}



} // namespace Ariadne
