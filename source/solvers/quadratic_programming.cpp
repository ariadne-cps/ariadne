/***************************************************************************
 *            quadratic_programming.cpp
 *
 *          Copyright 2018  Nicola Dessì
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without e ven the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "../config.hpp"
#include "../function/functional.hpp"

#include "../algebra/diagonal_matrix.hpp"
#include "../algebra/matrix.hpp"
#include "../algebra/vector.hpp"
#include "../function/affine.hpp"
#include "../numeric/numeric.hpp"
#include "../solvers/quadratic_programming.hpp"
#include "../utility/tuple.hpp"

#include "../output/logging.hpp"
#include "../utility/macros.hpp"

<<<<<<< HEAD
namespace Ariadne {
//----------------------------------------------------------------------------//
//            DEBUGGER PRINT UTILITY
//----------------------------------------------------------------------------//
void printVector(const Vector<FloatDP> &vc) {
  const unsigned n = vc.size();
  for (unsigned i = 0; i < n; ++i)
    std::cout << vc[i] << " ";
  std::cout << "\n";
}
void printMatrix(const Matrix<FloatDP> &mt) {
  const unsigned n = mt.row_size();
  const unsigned m = mt.column_size();
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < m; ++j)
      std::cout << mt[i][j] << " ";
    std::cout << "\n";
  }
  std::cout << "\n";
}
//----------------------------------------------------------------------------//

struct ASMQPSolver::StepData {
  friend std::ostream &operator<<(std::ostream &out,
                                  const struct StepData &stepData) {
    return out << "W: " << stepData.W << "\nS: " << stepData.S
               << "\ns: " << stepData.s << "\ny: " << stepData.y
               << "\nx: " << stepData.x << "\np: " << stepData.p
               << "\ng: " << stepData.g << "\npinvS: " << stepData.pinvS
               << "\n\nH: " << stepData.H << "\nd: " << stepData.d
               << "\nA_ori: " << stepData.A_ori << "\na_ori: " << stepData.A_ori
               << "\nB_ori: " << stepData.B_ori << "\nb_ori: " << stepData.b_ori
               << "\n\nu_H: " << stepData.u_H << "\nv_H: " << stepData.v_H
               << "\nR: " << stepData.R << "\n\nm_act: " << stepData.m_act
               << "\nm_tot: " << stepData.m_tot << "\nm_in: " << stepData.m_in
               << "\nn: " << stepData.n << "\n\nrtol: " << stepData.rtol
               << "\nalpha: " << stepData.alpha << "\niter: " << stepData.iter
               << "\n";
  }

  void print() { std::cout << *this << "\n"; }
  //! Structures used in algorithm
  std::vector<unsigned> W; //<! Working set
  Matrix<FloatDP> S;       //<! Active set (constraint on left)
  Vector<FloatDP> s;       //<! Active set (constraint on right)
  Vector<FloatDP> y;       //<! Lagrange multipliers
  Vector<FloatDP> y_b;     //<! Lagrange multipliers only for inequalities
  Vector<FloatDP> x;       //<! Current value of variables
  Vector<FloatDP> p;       //<! Direction
  Vector<FloatDP> g;       //<! q+Hx
  Matrix<FloatDP> pinvS;   //<! pinvS used

  //! Original problem variables
  Matrix<FloatDP> H;     //<! H matrix of the problem
  Matrix<FloatDP> Hinv;  //<! inverse of H, the matrix of the problem
  Vector<FloatDP> d;     //<! d vector of the problem
  Matrix<FloatDP> A_ori; //<! A orginal
  Vector<FloatDP> a_ori; //<! a original
  Matrix<FloatDP> B_ori; //<! B original
  Vector<FloatDP> b_ori; //<! b original

  //! Important characteristics of original problem
  Vector<FloatDP> u_H; //<! Minimum eigenvector of H
  FloatDP v_H;         //<! Minimum eigenvalue of H
  bool is_positive_H;  //<! TODO: delete this, not necessary
  Matrix<FloatDP>
      R; //<! Matrix used in inverse and other (cholesky decomposition)

  //! Counters
  unsigned m_act; //<! Active constraint (inequalities)
  unsigned m_tot; //<! Total constraint (both: inequalities and equalities)
  unsigned m_in;  //<! Inequalities constraint
  unsigned m_eq;  //<! Equalities contraint
  unsigned n;     //<! Variables

  //! Variables of algorithm
  FloatDP rtol; //<! Tolerance

  FloatDP alpha; //<! alpha value used in linesearch

  unsigned iter = 0; //<! Iteration
};

void ASMQPSolver::initialize_step_data(
    const Matrix<FloatDP> &H, const Vector<FloatDP> &d,
    const Matrix<FloatDP> &A, const Vector<FloatDP> &a,
    const Matrix<FloatDP> &B, const Vector<FloatDP> &b,
    const Vector<FloatDP> &x, struct StepData &v) const {
  v.H = H;
  v.d = d;
  v.A_ori = A;
  v.a_ori = a;
  v.B_ori = B;
  v.b_ori = b;
  v.x = x;

  v.rtol = static_cast<FloatDP>(sqrt(std::numeric_limits<double>::epsilon()));
  v.n = H.column_size();
  v.m_eq = A.row_size();
  v.m_in = B.row_size();
  v.m_act = v.m_eq;
  v.m_tot = v.m_eq + v.m_in;

  Tuple<Vector<FloatDP>, FloatDP> eigH = eigen_eigs(H);

  v.u_H = std::get<0>(eigH);
  v.v_H = std::get<1>(eigH);

  if (v.v_H > 0)
    v.Hinv = inverse(v.H);

  v.y_b.resize(v.m_in);
  v.y = Vector<FloatDP>::one(v.m_in + v.m_eq);

  v.S = A;
  v.s = a;

  if (!feasible(A, a, B, b, v.x, v.rtol)) {
    /// Phase I: Feasibility phase
    feasible_hotstart(v.x, A, a, B, b, v.rtol);
  }

  if (v.m_in > 0) {
    Vector<FloatDP> res = B * v.x - b;
    for (unsigned i = 0; i < v.m_in; ++i) {
      res[i] /= (1.0 + abs(b[i]));

      if (res[i] < v.rtol) {
        v.m_act++;
        v.S.resize(v.m_act, v.n);
        v.s.resize(v.m_act);
        v.W.resize(v.m_act);
        for (unsigned j = 0; j < v.n; ++j) {
          v.S[v.m_act - 1][j] = B[i][j];
        }
        v.s[v.m_act - 1] = b[i];
        v.W[v.m_act - 1] = i;
=======
namespace Ariadne
{
  //----------------------------------------------------------------------------//
  //            DEBUGGER PRINT UTILITY
  //----------------------------------------------------------------------------//
  void printVector(const Vector<FloatDP>& vc)
  {
    const unsigned n = vc.size();
    for(unsigned i = 0;i<n;++i)
      std::cout<<vc[i]<<" ";
    std::cout<<"\n";
  }
  void printMatrix(const Matrix<FloatDP>& mt)
  {
    const unsigned n = mt.row_size();
    const unsigned m = mt.column_size();
    for(unsigned i = 0;i<n;++i)
    {
      for(unsigned j=0;j<m;++j)
        std::cout<<mt[i][j]<<" ";
      std::cout<<"\n";
    }
    std::cout<<"\n";
  }
  //----------------------------------------------------------------------------//

struct ASMQPSolver::StepData
{
  friend std::ostream& operator<<(std::ostream& out, const struct StepData& stepData)
  {
   return out << "W: "<<stepData.W<<"\nS: "<<stepData.S<<"\ns: "<<stepData.s
              <<"\ny: "<<stepData.y<<"\nx: "<<stepData.x<<"\np: "<<stepData.p
              <<"\ng: "<<stepData.g<<"\npinvS: "<<stepData.pinvS<<"\n\nH: "
              <<stepData.H<<"\nd: "<<stepData.d<<"\nA_ori: "<<stepData.A_ori
              <<"\na_ori: "<<stepData.A_ori<<"\nB_ori: "<<stepData.B_ori
              <<"\nb_ori: "<<stepData.b_ori<<"\n\nu_H: "<<stepData.u_H
              <<"\nv_H: "<<stepData.v_H<<"\nR: "<<stepData.R<<"\n\nm_act: "
              <<stepData.m_act<<"\nm_tot: "<<stepData.m_tot<<"\nm_in: "
              <<stepData.m_in<<"\nn: "<<stepData.n<<"\n\nrtol: "<<stepData.rtol
              <<"\nalpha: "<<stepData.alpha<<"\niter: "<<stepData.iter<<"\n";
  }

  void print()
  {
    std::cout<<*this<<"\n";
  }
  //! Structures used in algorithm
  std::vector<unsigned>   W;      //<! Working set
  Matrix<FloatDP>         S;      //<! Active set (constraint on left)
  Vector<FloatDP>         s;      //<! Active set (constraint on right)
  Vector<FloatDP>         y;      //<! Lagrange multipliers
  Vector<FloatDP>         y_b;    //<! Lagrange multipliers only for inequalities
  Vector<FloatDP>         x;      //<! Current value of variables
  Vector<FloatDP>         p;      //<! Direction
  Vector<FloatDP>         g;      //<! q+Hx
  Matrix<FloatDP>         pinvS;  //<! pinvS used

  //! Original problem variables
  Matrix<FloatDP>         H;      //<! H matrix of the problem
  Matrix<FloatDP>         Hinv;   //<! inverse of H, the matrix of the problem
  Vector<FloatDP>         d;      //<! d vector of the problem
  Matrix<FloatDP>         A_ori;  //<! A orginal
  Vector<FloatDP>         a_ori;  //<! a original
  Matrix<FloatDP>         B_ori;  //<! B original
  Vector<FloatDP>         b_ori;  //<! b original

  //! Important characteristics of original problem
  Vector<FloatDP>         u_H;    //<! Minimum eigenvector of H
  FloatDP                 v_H;    //<! Minimum eigenvalue of H
  bool                    is_positive_H; //<! TODO: delete this, not necessary
  Matrix<FloatDP>         R;      //<! Matrix used in inverse and other (cholesky decomposition)

  //! Counters
  unsigned                m_act;  //<! Active constraint (inequalities)
  unsigned                m_tot;  //<! Total constraint (both: inequalities and equalities)
  unsigned                m_in;   //<! Inequalities constraint
  unsigned                m_eq;   //<! Equalities contraint
  unsigned                n;      //<! Variables

  //! Variables of algorithm
  FloatDP                 rtol;   //<! Tolerance

  FloatDP                 alpha;  //<! alpha value used in linesearch

  unsigned                iter = 0;   //<! Iteration
};

void
ASMQPSolver::initialize_step_data(const Matrix<FloatDP> &H,
                                  const Vector<FloatDP> &d,
                                  const Matrix<FloatDP> &A,
                                  const Vector<FloatDP> &a,
                                  const Matrix<FloatDP> &B,
                                  const Vector<FloatDP> &b,
                                  const Vector<FloatDP> &x,
                                  struct StepData &v) const
{
  v.H           = H;
  v.d           = d;
  v.A_ori       = A;
  v.a_ori       = a;
  v.B_ori       = B;
  v.b_ori       = b;
  v.x           = x;

  v.rtol        = static_cast<FloatDP>(sqrt(std::numeric_limits<double>::epsilon()));
  v.n           = H.column_size();
  v.m_eq        = A.row_size();
  v.m_in        = B.row_size();
  v.m_act       = v.m_eq;
  v.m_tot       = v.m_eq+v.m_in;

  Tuple<Vector<FloatDP>, FloatDP>
    eigH        = eigen_eigs(H);

  v.u_H         = std::get<0>(eigH);
  v.v_H         = std::get<1>(eigH);

  if(v.v_H>0)
    v.Hinv = inverse(v.H);



  v.y_b.resize(v.m_in);
  v.y = Vector<FloatDP>::one(v.m_in+v.m_eq);

  v.S           = A;
  v.s           = a;


  if(!feasible(A,a,B,b,v.x,v.rtol))
  {
    /// Phase I: Feasibility phase
    feasible_hotstart(v.x,A,a,B,b,v.rtol);
  }

  if(v.m_in>0)
  {
    Vector<FloatDP> res = B*v.x-b;
    for(unsigned i=0;i<v.m_in;++i)
    {
      res[i]/=(1.0+abs(b[i]));

      if(res[i]<v.rtol)
      {
        v.m_act++;
        v.S.resize(v.m_act,v.n);
        v.s.resize(v.m_act);
        v.W.resize(v.m_act);
        for(unsigned j = 0;j<v.n;++j)
        {
          v.S[v.m_act-1][j]=B[i][j];
        }
        v.s[v.m_act-1]=b[i];
        v.W[v.m_act-1]=i;
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
      }
    }
  }
}

Tuple<FloatDP, Vector<FloatDP>, Vector<FloatDP>>
ASMQPSolver::minimise(const Matrix<FloatDP> &Q, const Vector<FloatDP> &d,
<<<<<<< HEAD
                      const Vector<FloatDP> &xl, const Vector<FloatDP> &xu,
                      const Matrix<FloatDP> &C, const Vector<FloatDP> &l,
                      const Vector<FloatDP> &u, int &caller_status,
                      const Vector<FloatDP> x0) const {
  const unsigned n = Q.column_size();

  // const unsigned        m = C.row_size();
  // assert(Q.row_size()==Q.column_size());
  // assert(d.size()==n);
  // assert(xl.size()==n || xl.size()==0);
  // assert(xu.size()==n || xu.size()==0);
  // assert(C.column_size()==n);
  // assert(l.size()==m);
  // assert(u.size()==m);

  RawFloatVector x(n, 0);
  if (x0.size() > 0)
    x = x0;

  RawFloatMatrix A, B;
  RawFloatVector a, b;
  make_ltuple(A, a, B, b) = normalize_problem(C, l, u, x, xl, xu);

  ARIADNE_LOG(4, "\nNormalized problem:\n\tA=" << A << ", a=" << a << "\n\tB="
                                               << B << ", b=" << b << "\n");

  struct StepData v;
  initialize_step_data(Q, d, A, a, B, b, x, v);

  ARIADNE_LOG(4, "Starting from: " << v.x << "\n");

  return _minimise(v, caller_status);
}

void ASMQPSolver::null_space(struct StepData &v) const {
  // -----------------------------------------------------------------------------
  // ########################################
  //            NULL SPACE METHOD
  // ########################################
  auto nullvS = eigen_null(v.S);
  Matrix<FloatDP> Z = std::get<0>(nullvS);
  unsigned zSize = v.n - std::get<1>(nullvS);
  Matrix<FloatDP> rH = transpose(Z) * v.H * Z;
  v.pinvS = eigen_pinv(v.S);
  bool positive = true;
  if (zSize > 0) {
    // Z^THZ must be positive definite
    try {
      Matrix<FloatDP> tR = chol(rH);
      // v.R = tR;
    } catch (IndefiniteMatrixException &ime) {
      positive = false;
    }
  }
  if (positive) {
    if (zSize > 0) {
      Matrix<FloatDP> rHinv = inverse(rH);
      Vector<FloatDP> pz = -rHinv * transpose(Z) * v.g;

      // Global step
      v.p = Z * pz;
    } else {
      v.p = Vector<FloatDP>(v.n);
    }
  } else {
    // Compute the most negative curvature
    auto eigs_rH = eigen_eigs(rH);
    Vector<FloatDP> ru_H = std::get<0>(eigs_rH);
    // FloatDP rv_H = std::get<1>(eigs_rH);
    v.p = Z * ru_H;
    if (transpose(v.p) * v.g > v.rtol)
      v.p = -v.p;
  }
=======
  const Vector<FloatDP> &xl, const Vector<FloatDP> &xu,
  const Matrix<FloatDP> &C, const Vector<FloatDP> &l,
  const Vector<FloatDP> &u, int &caller_status,
  const Vector<FloatDP> x0) const
  {
    const unsigned        n = Q.column_size();

    // const unsigned        m = C.row_size();
    // assert(Q.row_size()==Q.column_size());
    // assert(d.size()==n);
    // assert(xl.size()==n || xl.size()==0);
    // assert(xu.size()==n || xu.size()==0);
    // assert(C.column_size()==n);
    // assert(l.size()==m);
    // assert(u.size()==m);

    RawFloatVector        x(n,0);
    if(x0.size()>0)
      x=x0;

    RawFloatMatrix      A,B;
    RawFloatVector      a,b;
    make_ltuple(A,a,B,b)=normalize_problem(C,l,u,x,xl,xu);

    ARIADNE_LOG(4, "\nNormalized problem:\n\tA="<<A<<", a="<<a<<"\n\tB="<<B<<", b="<<b<<"\n");

    struct StepData v;
    initialize_step_data(Q,d,A,a,B,b,x,v);

    ARIADNE_LOG(4, "Starting from: "<<v.x<<"\n");

    return _minimise(v, caller_status);
  }

void
ASMQPSolver::null_space(struct StepData &v) const
{
  // -----------------------------------------------------------------------------
      // ########################################
      //            NULL SPACE METHOD
      // ########################################
      auto              nullvS = eigen_null(v.S);
      Matrix<FloatDP>   Z     = std::get<0>(nullvS);
      unsigned          zSize = v.n-std::get<1>(nullvS);
      Matrix<FloatDP>   rH    = transpose(Z)*v.H*Z;
      v.pinvS                 = eigen_pinv(v.S);
      bool              positive = true;
      if(zSize>0)
      {
        // Z^THZ must be positive definite
        try
        {
          Matrix<FloatDP> tR  = chol(rH);
          // v.R = tR;
        }
        catch(IndefiniteMatrixException ime)
        {
          positive = false;
        }
      }
      if(positive)
      {
        if(zSize>0)
        {
          Matrix<FloatDP> rHinv = inverse(rH);
          Vector<FloatDP> pz = -rHinv * transpose(Z)*v.g;

          // Global step
          v.p = Z*pz;
        }
        else
        {
          v.p = Vector<FloatDP>(v.n);
        }
      }
      else
      {
        // Compute the most negative curvature
        auto eigs_rH = eigen_eigs(rH);
        Vector<FloatDP> ru_H = std::get<0>(eigs_rH);
        FloatDP rv_H = std::get<1>(eigs_rH);
        v.p = Z*ru_H;
        if(transpose(v.p)*v.g>v.rtol)
          v.p=-v.p;
      }
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
  // -----------------------------------------------------------------------------
}

Tuple<Matrix<FloatDP>, Vector<FloatDP>, Matrix<FloatDP>, Vector<FloatDP>>
ASMQPSolver::normalize_problem(const Matrix<FloatDP> &A_in,
                               const Vector<FloatDP> &A_lb,
                               const Vector<FloatDP> &A_ub,
                               const Vector<FloatDP> &x,
                               const Vector<FloatDP> &x_lb,
<<<<<<< HEAD
                               const Vector<FloatDP> &x_ub) const {
  const Nat n = A_in.column_size();
=======
                               const Vector<FloatDP> &x_ub) const
{
  const Nat     n = A_in.column_size();
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85

  // OPTIMIZEME: implement a better solution to normalize the problem

  std::vector<unsigned> eq;
  std::vector<unsigned> in;
<<<<<<< HEAD
  RawFloatMatrix A, B;
  RawFloatVector a, b;

  if (x_lb.size() > 0) {
    for (unsigned i = 0; i < n; ++i) {
      if (x_lb[i] == x_ub[i]) {
=======
  RawFloatMatrix  A,B;
  RawFloatVector  a,b;

  if(x_lb.size()>0)
  {
    for(unsigned i = 0; i < n;++i)
    {
      if(x_lb[i]==x_ub[i])
      {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
        eq.push_back(i);
        continue;
      }
      in.push_back(i);
    }

<<<<<<< HEAD
    // Equalities on matrix A
    for (auto e : eq) {
      bool add = true;

      if (add) {
        const unsigned aSize = A.row_size();
        A.resize(aSize + 1, n), a.resize(aSize + 1);
        A[aSize][e] = 1;
        a[aSize] = x_lb[e];
      }
    }

    // Inequalities on matrix B
    for (auto e : in) {
      bool add = true;

      if (add) {
        unsigned bSize = B.row_size();
        if (!is_inf(x_lb[e])) {
          B.resize(bSize + 1, n), b.resize(bSize + 1);
          B[bSize][e] = 1;
          b[bSize] = x_lb[e];
          bSize++;
        }
        if (!is_inf(x_ub[e])) {
          B.resize(bSize + 1, n), b.resize(bSize + 1);
          for (unsigned i = 0; i < n; ++i) {
            B[bSize][e] = -1;
          }
          b[bSize] = -x_ub[e];
=======
    //Equalities on matrix A
    for(auto e : eq)
    {
      bool add = true;

      if(add)
      {
        const unsigned aSize = A.row_size();
        A.resize(aSize+1,n),a.resize(aSize+1);
        A[aSize][e]=1;
        a[aSize]=x_lb[e];
      }
    }

    //Inequalities on matrix B
    for(auto e : in)
    {
      bool add = true;

      if(add)
      {
        unsigned bSize = B.row_size();
        if(!is_inf(x_lb[e]))
        {
          B.resize(bSize+1,n),b.resize(bSize+1);
          B[bSize][e]=1;
          b[bSize]=x_lb[e];
          bSize++;
        }
        if(!is_inf(x_ub[e]))
        {
          B.resize(bSize+1,n),b.resize(bSize+1);
          for(unsigned i = 0; i<n;++i)
          {
            B[bSize][e]=-1;
          }
          b[bSize]=-x_ub[e];
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
        }
      }
    }
    std::vector<unsigned> empty_eq;
    std::vector<unsigned> empty_in;
<<<<<<< HEAD
    std::swap(empty_eq, eq);
    std::swap(empty_in, in);
  }

  for (unsigned i = 0; i < A_in.row_size(); ++i) {
    if (!(is_inf(-A_lb[i]) || is_inf(A_ub[i])) && A_lb[i] == A_ub[i]) {
=======
    std::swap(empty_eq,eq);
    std::swap(empty_in,in);
  }


  for(unsigned i = 0; i < A_in.row_size();++i)
  {
    if(!(is_inf(-A_lb[i]) || is_inf(A_ub[i]))&&A_lb[i]==A_ub[i])
    {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
      eq.push_back(i);
      continue;
    }
    in.push_back(i);
  }

<<<<<<< HEAD
  // A matrix
  for (auto e : eq) {
    bool add = true;

    if (add) {
      const unsigned aSize = A.row_size();
      A.resize(aSize + 1, n), a.resize(aSize + 1);
      for (unsigned i = 0; i < n; ++i)
        A[aSize][i] = A_in[e][i];
      a[aSize] = A_lb[e];
    }
  }
  // B matrix
  for (auto e : in) {
    bool add = true;

    if (add) {
      unsigned bSize = B.row_size();
      if (!is_inf(A_lb[e])) {
        B.resize(bSize + 1, n), b.resize(bSize + 1);
        for (unsigned i = 0; i < n; ++i) {
          B[bSize][i] = A_in[e][i];
        }
        b[bSize] = A_lb[e];
        bSize++;
      }
      if (!is_inf(A_ub[e])) {
        B.resize(bSize + 1, n), b.resize(bSize + 1);
        for (unsigned i = 0; i < n; ++i) {
          B[bSize][i] = -A_in[e][i];
        }
        b[bSize] = -A_ub[e];
      }
    }
  }
  return make_tuple(A, a, B, b);
}

void ASMQPSolver::feasible_hotstart(Vector<FloatDP> &x,
                                    const Matrix<FloatDP> &A,
                                    const Vector<FloatDP> &a,
                                    const Matrix<FloatDP> &B,
                                    const Vector<FloatDP> &b,
                                    const FloatDP &rtol) const {
  const unsigned n = x.size();
  const unsigned m_eq = A.row_size();
  const unsigned m_in = B.row_size();
  const unsigned m_tot = m_in + n - m_eq;
  Vector<FloatDP> x_bar = x;
  Vector<FloatDP> e(m_tot);
  unsigned zSize = 0;
  Matrix<FloatDP> I(m_in, m_in);
  Vector<FloatDP> x_lb(m_tot, -inf);
  Vector<FloatDP> x_ub(m_tot, inf);
  Matrix<FloatDP> IN;
  Vector<FloatDP> in;
  // if m_eq > 0
  Matrix<FloatDP> Z;

  bool eq_feasible = true, in_feasible = true;
  if (m_eq > 0) {
    Vector<FloatDP> Ax = A * x_bar;

    for (unsigned i = 0; i < m_eq; ++i) {
      if (Ax[i] != a[i]) {
=======
  //A matrix
  for(auto e : eq)
  {
    bool add = true;

    if(add)
    {
      const unsigned aSize = A.row_size();
      A.resize(aSize+1,n),a.resize(aSize+1);
      for(unsigned i=0;i<n;++i)
        A[aSize][i]=A_in[e][i];
      a[aSize]=A_lb[e];
    }
  }
  //B matrix
  for(auto e : in)
  {
    bool add = true;

    if(add)
    {
      unsigned bSize = B.row_size();
      if(!is_inf(A_lb[e]))
      {
        B.resize(bSize+1,n),b.resize(bSize+1);
        for(unsigned i = 0; i<n;++i)
        {
          B[bSize][i]=A_in[e][i];
        }
        b[bSize]=A_lb[e];
        bSize++;
      }
      if(!is_inf(A_ub[e]))
      {
        B.resize(bSize+1,n),b.resize(bSize+1);
        for(unsigned i = 0; i<n;++i)
        {
          B[bSize][i]=-A_in[e][i];
        }
        b[bSize]=-A_ub[e];
      }
    }
  }
  return make_tuple(A,a,B,b);
}

void
ASMQPSolver::feasible_hotstart(Vector<FloatDP> &x, const Matrix<FloatDP> &A,
                       const Vector<FloatDP> &a, const Matrix<FloatDP> &B,
                       const Vector<FloatDP> &b, const FloatDP &rtol) const
{
  const unsigned    n = x.size();
  const unsigned    m_eq = A.row_size();
  const unsigned    m_in = B.row_size();
  const unsigned    m_tot = m_in+n-m_eq;
  Vector<FloatDP>   x_bar = x;
  Vector<FloatDP>   e(m_tot);
  unsigned          zSize = 0;
  Matrix<FloatDP>   I(m_in,m_in);
  Vector<FloatDP>   x_lb(m_tot,-inf);
  Vector<FloatDP>   x_ub(m_tot,inf);
  Matrix<FloatDP>   IN;
  Vector<FloatDP>   in;
  // if m_eq > 0
  Matrix<FloatDP>   Z;

  bool eq_feasible = true,in_feasible=true;
  if(m_eq>0)
  {
    Vector<FloatDP> Ax = A*x_bar;

    for(unsigned i = 0; i<m_eq;++i)
    {
      if(Ax[i]!=a[i])
      {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
        eq_feasible = false;
        break;
      }
    }
  }
<<<<<<< HEAD
<<<<<<< HEAD
  if (!eq_feasible)
    x_bar = eigen_pinv(A) * a;

  if (m_in > 0) {
    if (m_eq > 0) {
      auto nullA = eigen_null(A);
      Z = std::get<0>(nullA);
      zSize = std::get<1>(nullA);
      if (zSize == 0)
        throw InfeasibleQuadraticProgram(
            "A is square and full rank, but x_bar is not feasible");
    }
    Vector<FloatDP> Bx = B * x_bar - b;
    for (unsigned i = 0; i < m_in; ++i) {
      if (Bx[i] < -rtol * (1 + abs(b[i]))) {
=======
=======
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
  if(!eq_feasible)
    x_bar = eigen_pinv(A)*a;

  if(m_in>0)
  {
    if(m_eq>0)
    {
      auto nullA = eigen_null(A);
      Z = std::get<0>(nullA);
      zSize = std::get<1>(nullA);
      if(zSize==0)
      throw InfeasibleQuadraticProgram("A is square and full rank, but x_bar is not feasible");
    }
    Vector<FloatDP> Bx = B*x_bar - b;
    for(unsigned i = 0; i<m_in;++i)
    {
      if(Bx[i]<-rtol*(1+abs(b[i])))
      {
<<<<<<< HEAD
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.
=======
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
        in_feasible = false;
        break;
      }
    }
  }

<<<<<<< HEAD
  if (in_feasible) {
    x = x_bar;
=======
  if(in_feasible)
  {
    x=x_bar;
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
    ARIADNE_LOG(5, "X is feasible, starting the algorithm");
    return;
  }

  // Try to solve the problem of feasible region with Simplex method
  //    position on a vertex of polytope
<<<<<<< HEAD
  for (unsigned i = n - m_eq; i < m_tot; ++i)
    e[i] = 1.0;
  for (unsigned i = 0; i < m_in; ++i)
    I[i][i] = 1.0;

  if (m_eq > 0) {
    IN.resize(m_in, m_in + Z.column_size());
    in.resize(m_in);
    Matrix<FloatDP> BZ = B * Z;
    for (unsigned i = 0; i < m_in; ++i) {
      for (unsigned j = 0; j < Z.column_size(); ++j)
        IN[i][j] = BZ[i][j];
      for (unsigned j = 0; j < m_in; ++j)
        IN[i][j + Z.column_size()] = I[i][j];
    }
    in = -B * x_bar - b;
  } else {
    IN.resize(m_in, n + m_in);
    in.resize(m_in);
    for (unsigned i = 0; i < m_in; ++i) {
      for (unsigned j = 0; j < n; ++j)
        IN[i][j] = B[i][j];
      for (unsigned j = 0; j < m_in; ++j)
        IN[i][j + n] = I[i][j];
=======
  for(unsigned i=n-m_eq ;i<m_tot;++i)
    e[i]=1.0;
  for(unsigned i=0;i<m_in;++i)
    I[i][i]=1.0;

  if(m_eq>0)
  {
    IN.resize(m_in,m_in+Z.column_size());
    in.resize(m_in);
    Matrix<FloatDP> BZ = B*Z;
    for(unsigned i=0;i<m_in;++i)
    {
      for(unsigned j=0;j<Z.column_size();++j)
        IN[i][j]=BZ[i][j];
      for(unsigned j=0;j<m_in;++j)
        IN[i][j+Z.column_size()]=I[i][j];
    }
    in = -B*x_bar - b;
  }
  else
  {
    IN.resize(m_in,n+m_in);
    in.resize(m_in);
    for(unsigned i=0;i<m_in;++i)
    {
      for(unsigned j=0;j<n;++j)
        IN[i][j]=B[i][j];
      for(unsigned j=0;j<m_in;++j)
        IN[i][j+n]=I[i][j];
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
      in[i] = b[i];
    }
  }

<<<<<<< HEAD
  for (unsigned i = n - m_eq; i < m_tot; ++i)
    x_lb[i] = 0.0;
  // printVector(x_lb);
  int errnum = 0;
<<<<<<< HEAD
  Vector<FloatDP> p = lp_min(e, IN, in, x_lb, errnum);
  // bool p_l_rtol = true;
  // FloatDP n2 = norm2(in);
  // FloatDP tolerance = rtol*(1+n2);
  // for(unsigned i=n-m_eq+1;i<p.size();++i)
  // {
  //   if(p[i]>=tolerance)
  //   {
  //     p_l_rtol=false;
  //     break;
  //   }
  // }
  if (errnum != 0) // || !p_l_rtol)
  {
=======
=======
  for(unsigned i=n-m_eq;i<m_tot;++i)
    x_lb[i]=0.0;
  // printVector(x_lb);
  int errnum = 0;
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
  Vector<FloatDP> p = lp_min(e,IN,in,x_lb,errnum);
  bool p_l_rtol = true;
  FloatDP n2 = norm2(in);
  FloatDP tolerance = rtol*(1+n2);
  for(unsigned i=n-m_eq+1;i<p.size();++i)
  {
    if(p[i]>=tolerance)
    {
      p_l_rtol=false;
      break;
    }
  }
  if(errnum!=0 || !p_l_rtol)
  {
<<<<<<< HEAD
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.
    // ARIADNE_WARN("QP subproblem is infeasible!");
    throw InfeasibleQuadraticProgram(std::to_string(errnum));
  }
  if (m_eq > 0) {
    Vector<FloatDP> p_n(n - m_eq);
    for (unsigned i = 0; i < n - m_eq; ++i)
      p_n[i] = p[i];
    x = x_bar + (Z * p_n);
  } else {
    for (unsigned i = 0; i < n; ++i)
      x[i] = p[i];
  }
  ARIADNE_LOG(4, "Found x feasible: " << x << "\n");
}

<<<<<<< HEAD
bool ASMQPSolver::feasible(const RawFloatMatrix &A, const RawFloatVector &a,
                           const RawFloatMatrix &B, const RawFloatVector &b,
                           const RawFloatVector &x, const FloatDP &rtol) const {
  bool ax_eq_a = true, bx_geq_b = true;
  if (A.row_size() > 0) {
    Vector<FloatDP> Ax = A * x;
    for (unsigned i = 0; i < a.size(); ++i) {
      if (Ax[i] != a[i]) {
=======
=======
    // ARIADNE_WARN("QP subproblem is infeasible!");
    throw InfeasibleQuadraticProgram(std::to_string(errnum));
  }
  if(m_eq>0)
  {
    Vector<FloatDP> p_n(n-m_eq);
    for(unsigned i=0;i<n-m_eq;++i)
      p_n[i]=p[i];
    x=x_bar+(Z*p_n);
  }
  else
  {
    for(unsigned i=0;i<n;++i)
      x[i]=p[i];
  }
  ARIADNE_LOG(4,"Found x feasible: "<<x<<"\n");
}

>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
bool
ASMQPSolver::feasible(const RawFloatMatrix &A, const RawFloatVector &a,
              const RawFloatMatrix &B, const RawFloatVector &b,
              const RawFloatVector &x,const FloatDP &rtol) const
{
  bool ax_eq_a = true,bx_geq_b=true;
  if(A.row_size()>0)
  {
    Vector<FloatDP> Ax = A*x;
    for(unsigned i = 0; i<a.size();++i)
    {
      if(Ax[i]!=a[i])
      {
<<<<<<< HEAD
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.
=======
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
        ax_eq_a = false;
        break;
      }
    }
  }
<<<<<<< HEAD
  if (B.row_size() > 0) {
    Vector<FloatDP> Bx = B * x;
    for (unsigned i = 0; i < b.size(); ++i) {
      if (Bx[i] < b[i]) {
=======
  if(B.row_size()>0)
  {
    Vector<FloatDP> Bx = B*x;
    for(unsigned i = 0; i<b.size();++i)
    {
      if(Bx[i]<b[i])
      {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
        bx_geq_b = false;
        break;
      }
    }
  }
<<<<<<< HEAD
  if (!ax_eq_a || !bx_geq_b) {
=======
  if(!ax_eq_a || !bx_geq_b)
  {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
    return false;
  }
  return true;
}

Tuple<FloatDP, Vector<FloatDP>, Vector<FloatDP>>
<<<<<<< HEAD
ASMQPSolver::_minimise(struct StepData &v, int &caller_staus) const {

  v.iter = 0;

  ARIADNE_LOG(3, "Starting from " << v.x << "\n");
=======
ASMQPSolver::_minimise(struct StepData &v, int &caller_staus) const
{

  v.iter = 0;

  ARIADNE_LOG(3, "Starting from "<<v.x<<"\n");
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85

  const unsigned max_k = 200;
  /// Phase II: Optimality phase
  QuadraticProgramStatus status =
      QuadraticProgramStatus::INDETERMINATE_FEASIBILITY;
  while (status == QuadraticProgramStatus::INDETERMINATE_FEASIBILITY &&
         v.iter < max_k) {

<<<<<<< HEAD
    status = _minimisation_step(v);
    v.iter++;

    ARIADNE_LOG(2, "x=" << v.x << ", f(x)="
                        << ((0.5 * (transpose(v.x) * v.H * v.x)) +
                            (transpose(v.x) * v.d))
                        << ", y=" << v.y_b << ", W=" << v.W << "\n");
  }

  if (v.iter == max_k)
    caller_staus = -2;

  v.y_b = (v.m_act > 0) ? transpose(v.pinvS) * (v.g + v.H * v.p) : EMPTY_VEC;
  Vector<FloatDP> lambda(v.m_tot);

  for (unsigned i = 0; i < v.m_act; ++i) {
    unsigned wi = v.W[i];
    lambda[wi] = v.y_b[i];
=======
      status = _minimisation_step(v);
      v.iter++;

      ARIADNE_LOG(2,
        "x=" << v.x << ", f(x)="
        << ((0.5 * (transpose(v.x) * v.H * v.x)) + (transpose(v.x) * v.d))
        << ", y=" << v.y_b <<", W="<<v.W<< "\n");
  }

  if (v.iter == max_k)
      caller_staus = -2;

  v.y_b = (v.m_act>0)?transpose(v.pinvS)*(v.g + v.H * v.p):EMPTY_VEC;
  Vector<FloatDP> lambda(v.m_tot);

  for(unsigned i = 0; i<v.m_act;++i)
  {
    unsigned wi = v.W[i];
    lambda[wi]=v.y_b[i];
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
  }

  FloatDP f_obj = ((0.5 * (transpose(v.x) * v.H * v.x)) + transpose(v.x) * v.d);

  return make_tuple(f_obj, v.x, lambda);
}

<<<<<<< HEAD
QuadraticProgramStatus
ASMQPSolver::_minimisation_step(struct StepData &v) const {

  v.g = v.d + (v.H * v.x);

  if (v.m_act == 0) {
    if (v.v_H > 0) //<! H is positive definite
    {
      // v.R=eigen_chol(v.H);
      // Matrix<FloatDP> Hinv = inverse(v.H);
      v.p = -(v.Hinv * v.g);
    } else {
      v.p = v.u_H;
      if (transpose(v.p) * v.g > v.rtol)
        v.p = -v.p;
    }
    v.y_b = Vector<FloatDP>(v.m_in);
  } else {
=======

QuadraticProgramStatus ASMQPSolver::_minimisation_step(struct StepData &v) const
{

  v.g = v.d+(v.H*v.x);

  if(v.m_act == 0)
  {
    if(v.v_H>0) //<! H is positive definite
    {
      // v.R=eigen_chol(v.H);
      // Matrix<FloatDP> Hinv = inverse(v.H);
      v.p = -(v.Hinv*v.g);
    }
    else
    {
      v.p = v.u_H;
      if(transpose(v.p)*v.g>v.rtol)
        v.p=-v.p;
    }
    v.y_b = Vector<FloatDP>(v.m_in);
  }
  else
  {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
    null_space(v);
  }

  bool p_eq_0 = true;
<<<<<<< HEAD
  for (unsigned i = 0; i < v.n; ++i) {
    if (abs(v.p[i]) > v.rtol) {
=======
  for(unsigned i = 0; i<v.n;++i)
  {
    if(abs(v.p[i])>v.rtol)
    {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
      p_eq_0 = false;
      break;
    }
  }

<<<<<<< HEAD
  if (p_eq_0) {
    if (v.m_act - v.m_eq == 0) {
      return QuadraticProgramStatus::PRIMAL_FEASIBLE;
    }
    // Computing multipliers only for active inequalities
    Matrix<FloatDP> pinvS_in(v.m_act - v.m_eq, v.n);
    for (unsigned i = 0; i < v.m_act - v.m_eq; ++i) {
      for (unsigned j = 0; j < v.n; ++j)
        pinvS_in[i][j] = v.pinvS[j][i + v.m_eq];
    }

    Vector<FloatDP> y_tmp = pinvS_in * (v.g + (v.H * v.p));
    for (unsigned i = 0; i < v.m_act; ++i) {
      v.y_b[v.W[i]] = y_tmp[i];
    }
    unsigned index_min = 0;
    FloatDP min_y_b = y_tmp[index_min];
    for (unsigned i = 1; i < v.m_act - v.m_eq; ++i) {
      if (min_y_b > y_tmp[i]) {
=======
  if(p_eq_0)
  {
    if(v.m_act-v.m_eq == 0)
    {
      return QuadraticProgramStatus::PRIMAL_FEASIBLE;
    }
    // Computing multipliers only for active inequalities
    Matrix<FloatDP> pinvS_in(v.m_act-v.m_eq, v.n);
    for(unsigned i = 0; i<v.m_act-v.m_eq;++i)
    {
      for(unsigned j = 0;j<v.n;++j)
        pinvS_in[i][j]=v.pinvS[j][i+v.m_eq];
    }

    Vector<FloatDP> y_tmp = pinvS_in*(v.g+(v.H*v.p));
    for(unsigned i = 0; i<v.m_act;++i)
    {
      v.y_b[v.W[i]]=y_tmp[i];
    }
    unsigned index_min = 0;
    FloatDP min_y_b = y_tmp[index_min];
    for(unsigned i = 1; i<v.m_act-v.m_eq;++i)
    {
      if(min_y_b>y_tmp[i])
      {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
        min_y_b = y_tmp[i];
        index_min = i;
      }
    }
<<<<<<< HEAD
    if (min_y_b >= 0) {
=======
    if(min_y_b>=0)
    {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
      return QuadraticProgramStatus::PRIMAL_FEASIBLE;
    }

    // One multiplier is negative, remove it from the Set
<<<<<<< HEAD
    //  starting from that point, slit back all others
    v.m_act--;
    for (unsigned i = index_min; i < v.m_act - v.m_eq; ++i) {
      v.W[i] = v.W[i + 1];
      for (unsigned j = 0; j < v.n; ++j)
        v.S[i + v.m_eq][j] = v.S[i + v.m_eq + 1][j];
      v.s[i + v.m_eq] = v.s[i + v.m_eq + 1];
    }
    v.W.resize(v.m_act - v.m_eq);
    v.s.resize(v.m_act);
    v.S.resize(v.m_act, v.n);
  } else {
    // if all cons are active, take a full step
    if (v.m_act - v.m_eq == v.m_in) {
      v.x += v.p;
    } else {
      v.alpha = linesearch(v);
      v.x += v.alpha * v.p;
      ARIADNE_LOG(4, "Done a step of size: " << v.alpha << "*" << v.p << "\n");
=======
    //  starting from that point, shit back all others
    v.m_act--;
    for(unsigned i = index_min;i<v.m_act-v.m_eq;++i)
    {
      v.W[i]=v.W[i+1];
      for(unsigned j= 0; j<v.n;++j)
        v.S[i+v.m_eq][j]=v.S[i+v.m_eq+1][j];
      v.s[i+v.m_eq]=v.s[i+v.m_eq+1];
    }
    v.W.resize(v.m_act-v.m_eq);
    v.s.resize(v.m_act);
    v.S.resize(v.m_act,v.n);
  }
  else
  {
    // if all cons are active, take a full step
    if(v.m_act -v.m_eq == v.m_in)
    {
      v.x+=v.p;
    }
    else
    {
      v.alpha = linesearch(v);
      v.x+=v.alpha*v.p;
      ARIADNE_LOG(4,"Done a step of size: "<<v.alpha<<"*"<<v.p<<"\n");
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
    }
  }

  return QuadraticProgramStatus::INDETERMINATE_FEASIBILITY;
}

<<<<<<< HEAD
FloatDP ASMQPSolver::linesearch(struct StepData &v) const {
  FloatDP alpha = 1.0;
  bool block = false;
  unsigned uBlock = 0;

  for (unsigned i = 0; i < v.m_in; ++i) {
=======
FloatDP
ASMQPSolver::linesearch(struct StepData &v) const
{
  FloatDP   alpha = 1.0;
  bool      block = false;
  unsigned  uBlock = 0;

  for(unsigned i = 0; i<v.m_in;++i)
  {
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
    // OPTIMIZEME: this can be optimized
    if (std::count(v.W.begin(), v.W.end(), i))
      continue;

<<<<<<< HEAD
    Vector<FloatDP> b_i(v.n);
    for (unsigned k = 0; k < v.n; ++k)
      b_i[k] = v.B_ori[i][k];
    FloatDP b_i_p = transpose(b_i) * v.p;
    FloatDP b_i_x = transpose(b_i) * v.x;

    if (b_i_p < 0.0) {
      FloatDP alpha_bar = (v.b_ori[i] - b_i_x) / b_i_p;
      if (alpha_bar < alpha) {
        alpha = alpha_bar;
        uBlock = i;
        block = true;
      }
    }
  }
  if (block) // found a blocking cons, add it!
  {
    v.m_act++;
    v.W.resize(v.m_act - v.m_eq, uBlock);
    v.s.resize(v.m_act);
    v.s[v.m_act - 1] = uBlock;
    v.S.resize(v.m_act, v.n);
    for (unsigned i = 0; i < v.n; ++i) {
      v.S[v.m_act - 1][i] = v.B_ori[uBlock][i];
=======

      Vector<FloatDP> b_i(v.n);
      for(unsigned k = 0;k<v.n;++k)
        b_i[k]=v.B_ori[i][k];
      FloatDP b_i_p = transpose(b_i)*v.p;
      FloatDP b_i_x = transpose(b_i)*v.x;

      if(b_i_p<0.0)
      {
        FloatDP alpha_bar = (v.b_ori[i]-b_i_x)/b_i_p;
        if(alpha_bar < alpha)
        {
          alpha = alpha_bar;
          uBlock = i;
          block=true;
        }
      }

  }
  if(block) // found a blocking cons, add it!
  {
    v.m_act++;
    v.W.resize(v.m_act-v.m_eq,uBlock);
    v.s.resize(v.m_act);
    v.s[v.m_act-1]=uBlock;
    v.S.resize(v.m_act,v.n);
    for(unsigned i = 0; i<v.n;++i)
    {
      v.S[v.m_act-1][i]=v.B_ori[uBlock][i];
>>>>>>> 681346c6af58fdfff85dfaa109ead700efe84d85
    }
  }
  return alpha;
}
} // namespace Ariadne
