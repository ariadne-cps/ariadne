#pragma once

namespace Ariadne {
#if defined HAVE_GLPK_H
template<class X>
Vector<X>
lp_min(const Vector<X>& C,
       const Matrix<X>& A,
       const Vector<X>& b,
       const Vector<X>& lb,
       int& errnum)
{
  // constants
  const int n = C.size();
  const int m = b.size();
  //  Fill the problem
  int* ia = new int[static_cast<unsigned>(n * m + 1)];
  int* ja = new int[static_cast<unsigned>(n * m + 1)];
  double* ar = new double[static_cast<unsigned>(n * m + 1)];

  glp_prob* lp = glp_create_prob();
  glp_smcp smcp;

  // Result
  Vector<X> xmin(static_cast<unsigned>(n));

  int index = 1;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      ia[static_cast<unsigned>(index)] = i + 1,
      ja[static_cast<unsigned>(index)] = j + 1,
      ar[static_cast<unsigned>(index)] = static_cast<double>(
        A[static_cast<unsigned>(i)][static_cast<unsigned>(j)].dbl);
      index++;
    }
  }

  glp_set_obj_dir(lp, GLP_MIN);

  glp_add_cols(lp, n);
  for (int i = 0; i < n; ++i) {
    if (is_inf(lb[static_cast<unsigned>(i)])) {
      glp_set_col_bnds(lp,
                       i + 1,
                       GLP_FR,
                       -std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity());
    } else {
      glp_set_col_bnds(lp,
                       i + 1,
                       GLP_LO,
                       static_cast<double>(lb[static_cast<unsigned>(i)].dbl),
                       std::numeric_limits<double>::infinity());
    }
    glp_set_obj_coef(
      lp, i + 1, static_cast<double>(C[static_cast<unsigned>(i)].dbl));
  }

  glp_add_rows(lp, m);
  for (int i = 0; i < m; ++i) {
    glp_set_row_bnds(lp,
                     i + 1,
                     GLP_LO,
                     static_cast<double>(b[static_cast<unsigned>(i)].dbl),
                     static_cast<double>(b[static_cast<unsigned>(i)].dbl));
  }

  glp_load_matrix(lp, index - 1, ia, ja, ar);

  glp_init_smcp(&smcp);
  smcp.msg_lev = 1;
  smcp.meth = 1;
  smcp.pricing = 34;
  smcp.r_test = 34;
  smcp.tol_bnd = 1e-7;
  smcp.tol_dj = 1e-7;
  smcp.tol_piv = 1e-10;
  smcp.obj_ll = std::numeric_limits<double>::max();
  smcp.obj_ul = std::numeric_limits<double>::max();
  smcp.it_lim = std::numeric_limits<int>::max();
  smcp.tm_lim = std::numeric_limits<int>::max();
  smcp.out_frq = 200;
  smcp.out_dly = 0;
  smcp.presolve = 1;
  errnum = glp_simplex(lp, &smcp);

  for (int i = 0; i < n; ++i)
    xmin[static_cast<unsigned>(i)] =
      static_cast<X>(glp_get_col_prim(lp, i + 1));

  glp_delete_prob(lp);
  glp_free_env();
  delete[] ia;
  delete[] ja;
  delete[] ar;

  return xmin;
}

#endif

} // namespace Ariadne

#include "glpk_interface.impl.hpp"