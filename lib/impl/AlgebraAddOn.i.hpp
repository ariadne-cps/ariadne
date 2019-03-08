namespace Ariadne {
#define ci(a) (static_cast<int>(a))
#define cd(a) (static_cast<double>(a))
#define cx(a) (static_cast<X>(a))
#define cu(a) (static_cast<unsigned>(a))

template <class X>
Vector<X> lp_min(const Vector<X> &C, const Matrix<X> &A, const Vector<X> &b,
                 const Vector<X> &lb, int &errnum) {
  // constants
  const int n = C.size();
  const int m = b.size();
  //  Fill the problem
  int *ia = new int[cu(n * m + 1)];
  int *ja = new int[cu(n * m + 1)];
  double *ar = new double[cu(n * m + 1)];

  glp_prob *lp = glp_create_prob();
  glp_smcp smcp;

  // Result
  Vector<X> xmin(cu(n));

  int index = 1;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      ia[cu(index)] = i + 1, ja[cu(index)] = j + 1,
      ar[cu(index)] = cd(A[cu(i)][cu(j)].dbl);
      index++;
    }
  }

  glp_set_obj_dir(lp, GLP_MIN);

  glp_add_cols(lp, n);
  for (int i = 0; i < n; ++i) {
    if (is_inf(lb[cu(i)])) {
      glp_set_col_bnds(lp, i + 1, GLP_FR,
                       -std::numeric_limits<double>::infinity(),
                       std::numeric_limits<double>::infinity());
    } else {
      glp_set_col_bnds(lp, i + 1, GLP_LO, cd(lb[cu(i)].dbl),
                       std::numeric_limits<double>::infinity());
    }
    glp_set_obj_coef(lp, i + 1, cd(C[cu(i)].dbl));
  }

  glp_add_rows(lp, m);
  for (int i = 0; i < m; ++i) {
    glp_set_row_bnds(lp, i + 1, GLP_LO, cd(b[cu(i)].dbl), cd(b[cu(i)].dbl));
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
    xmin[cu(i)] = cx(glp_get_col_prim(lp, i + 1));

  glp_delete_prob(lp);
  glp_free_env();
  delete[] ia;
  delete[] ja;
  delete[] ar;

  return xmin;
}

template <class X>
Tuple<Matrix<X>, Matrix<X>, PivotMatrix>
orthogonal_decomposition(const Matrix<X> &A, Bool allow_pivoting) {
  X zero = A.zero_element();
  X one = zero + 1;
  SizeType m = A.row_size();
  SizeType n = A.column_size();
  Matrix<X> Q(m, m);
  Matrix<X> R(A);
  PivotMatrix P(n);
  Array<X> p(n);
  Vector<X> u(m);
  for (SizeType i = 0; i != m; ++i) {
    for (SizeType j = 0; j != m; ++j) {
      Q[i][j] = zero;
    }
    Q[i][i] = one;
  }
  for (SizeType k = 0; k != std::min(m, n); ++k) {
    // std::cerr<<"k="<<k<<" Q="<<Q<<" R="<<R<<std::flush;
    if (allow_pivoting) {
      // Find a pivot column
      SizeType pivot_column = k;
      X max_column_norm = zero;
      for (SizeType j = k; j != n; ++j) {
        X column_norm = zero;
        for (SizeType i = k; i != m; ++i) {
          column_norm += abs(R[i][j]);
        }
        if (decide(column_norm > max_column_norm)) {
          pivot_column = j;
          max_column_norm = column_norm;
        }
      }
      SizeType l = pivot_column;
      // Swap working column and pivot column
      for (SizeType i = 0; i != m; ++i) {
        std::swap(R[i][l], R[i][k]);
      }
      // Set pivot column in result
      P[k] = l;
    }
    // Compute |a| where a is the working column
    X nrmas = zero;
    for (SizeType i = k; i != m; ++i) {
      nrmas += R[i][k] * R[i][k];
    }
    X nrma = sqrt(nrmas);
    // Compute u=a +/- |a|e
    for (SizeType i = 0; i != k; ++i) {
      u[i] = zero;
    }
    for (SizeType i = k; i != m; ++i) {
      u[i] = R[i][k];
    }
    if (decide(u[k] >= 0)) {
      u[k] += nrma;
    } else {
      u[k] -= nrma;
    }
    // Compute -2/u.u
    X nrmus = zero;
    for (SizeType i = k; i != m; ++i) {
      nrmus += sqr(u[i]);
    }
    X mtdnu = (-2) / nrmus;
    // Compute H=(1-2uu'/u'u)
    // Matrix<X> H(n,n); for(SizeType i=0; i!=n; ++i) {
    // H[i][i]=1.0; for(SizeType j=0; j!=n; ++j) { H[i][j]+=u[i]*u[j]*mtdnu; } }
    // For each column b of R, compute b-=2u(u.b)/u.u
    for (SizeType j = k; j != n; ++j) {
      X udtb = zero;
      for (SizeType i = k; i != m; ++i) {
        udtb += u[i] * R[i][j];
      }
      X scl = udtb * mtdnu;
      for (SizeType i = k; i != m; ++i) {
        R[i][j] += scl * u[i];
      }
    }
    // For the kth column, set R[k][k]=-/+ |a|
    // and R[i][k]=0 for i>k
    for (SizeType i = k + 1; i != m; ++i) {
      R[i][k] = zero;
    }
    if (decide(u[k] >= 0)) {
      R[k][k] = -nrma;
    } else {
      R[k][k] = nrma;
    }
    // Update Q'=QH = Q(I-2uu'/u'u)
    // For each row q, compute q-=2u(u.q)/(u.u)
    for (SizeType i = 0; i != m; ++i) {
      X qdtu = zero;
      for (SizeType j = k; j != m; ++j) {
        qdtu += Q[i][j] * u[j];
      }
      X scl = qdtu * mtdnu;
      for (SizeType j = k; j != m; ++j) {
        Q[i][j] += scl * u[j];
      }
    }
  }
  return std::make_tuple(Q, R, P);
}

template <class X> Matrix<X> pinv(const Matrix<X> &G) {
  unsigned int m = G.row_size();
  unsigned int n = G.column_size();
  bool rect = false;
  Matrix<X> A;
  if (m < n) {
    rect = true;
    A = G * transpose(G);
    n = m;
  } else {
    A = transpose(G) * G;
  }
  X minDiag = A[0][0];
  for (unsigned int i = 0; i < A.row_size(); ++i) {
    if (minDiag > A[i][i])
      minDiag = A[i][i];
  }
  X tol = minDiag * 1e-9;
  Matrix<X> L = Matrix<X>::zero(A.row_size(), A.column_size());
  unsigned int r = 0;
  for (unsigned int k = 0; k < n; ++k) {
    r++;
    for (unsigned int i = k; i < n; ++i) {
      L[i][r - 1] = A[i][k];
      X l_wt_l = 0;
      for (unsigned int j = 0; j < r - 1; ++j)
        l_wt_l += L[i][j] * L[k][j];
      L[i][r - 1] -= l_wt_l;
    }
    if (L[k][r - 1] > tol) {
      L[k][r - 1] = sqrt(L[k][r - 1]);
      if (k < n) {
        for (unsigned int j = k + 1; j < n; j++) {
          L[j][r - 1] = L[j][r - 1] / L[k][r - 1];
        }
      }
    } else {
      r--;
    }
  }
  Matrix<X> I(L.row_size(), r);
  for (unsigned int i = 0; i < I.row_size(); ++i) {
    for (unsigned int j = 0; j < r; ++j) {
      I[i][j] = L[i][j];
    }
  }
  Matrix<X> M = inverse(transpose(I) * I);
  if (rect) {
    return transpose(G) * I * M * M * transpose(I);
  } else {
    return I * M * M * transpose(I) * transpose(G);
  }
}

template <class X> Tuple<Matrix<X>, unsigned> null(const Matrix<X> &A) {
  unsigned n = A.column_size();
  unsigned m = A.row_size();
  Matrix<X> Q;
  Matrix<X> R;
  PivotMatrix P;
  // Compute N(A) as the left null space of A' that is the null space of A
  Matrix<X> A_T = transpose(A);
  make_ltuple(Q, R, P) = orthogonal_decomposition(A_T, false);
  // Assumption A has full row rank
  unsigned rank_A = m;

  if (rank_A >= n) {
    return make_tuple(Matrix<X>::zero(n, m), 0);
  }
  // Last n-rank(A) cols are the left null space of A' = N(A)
  Matrix<X> Z(n, n - rank_A);
  unsigned zSize = Z.column_size();
  for (unsigned i = 0; i < Z.row_size(); ++i) {
    for (unsigned j = 0; j < zSize; ++j)
      Z[i][j] = Q[i][j + rank_A];
  }
  return make_tuple(Z, zSize);
}

template <class X> Matrix<X> chol(const Matrix<X> &A) {
  unsigned n = A.row_size();
  Matrix<X> chol(n, n, 0);
  unsigned i, j, k;

  for (i = 0; i < n; i++) {
    chol[i][i] = A[i][i];
    for (k = 0; k < i; k++)
      chol[i][i] -= chol[k][i] * chol[k][i];
    if (chol[i][i] <= 0) {
      throw IndefiniteMatrixException("The matrix is not positive definite");
      return chol;
    }
    chol[i][i] = sqrt(chol[i][i]);

    for (j = i + 1; j < n; j++) {
      chol[i][j] = A[i][j];
      for (k = 0; k < i; k++)
        chol[i][j] -= chol[k][i] * chol[k][j];
      chol[i][j] /= chol[i][i];
    }
  }

  return chol;
}

template <class X>
Tuple<Matrix<X>, Vector<X>> jacobi_eigs(const Matrix<X> &matrix) {
  Matrix<X> A = matrix;
  unsigned n = A.row_size();
  unsigned it_max = 100;
  unsigned it_num;
  unsigned rot_num;

  Vector<X> bw(n);
  Vector<X> zw(n);
  X thresh;
  X gapq;
  X termp;
  X h;
  X c;
  X g;
  X s;
  X t;
  X tau;
  X term;
  X termq;
  X theta;
  X w;

  Matrix<X> eigs = Matrix<X>::identity(n, n);
  Vector<X> d(n); //=diag of A

  unsigned i, j, k, l, m, p, q;
  for (unsigned index = 0; index < n; ++index)
    d[index] = A[index][index];
  // r8mat_identity ( n, v );
  // r8mat_diag_get_vector ( n, a, d );

  for (i = 0; i < n; i++) {
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  it_num = 0;
  rot_num = 0;

  while (it_num < it_max) {
    it_num = it_num + 1;
    //
    //  The convergence threshold is based on the size of the elements in
    //  the strict upper triangle of the matrix.
    //
    thresh = 0.0;
    for (j = 0; j < n; j++) {
      for (i = 0; i < j; i++) {
        thresh = thresh + A[i][j] * A[i][j];
      }
    }

    thresh = sqrt(thresh) / static_cast<X>(4 * n);

    if (thresh == 0.0) {
      break;
    }

    for (p = 0; p < n; p++) {
      for (q = p + 1; q < n; q++) {
        gapq = 10.0 * abs(A[p][q]);
        termp = gapq + abs(d[p]);
        termq = gapq + abs(d[q]);
        //
        //  Annihilate tiny offdiagonal elements.
        //
        if (4 < it_num && termp == abs(d[p]) && termq == abs(d[q])) {
          A[p][q] = 0.0;
        }
        //
        //  Otherwise, apply a rotation.
        //
        else if (thresh <= abs(A[p][q])) {
          h = d[q] - d[p];
          term = abs(h) + gapq;

          if (term == abs(h)) {
            t = A[p][q] / h;
          } else {
            theta = 0.5 * h / A[p][q];
            t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0) {
              t = -t;
            }
          }
          c = 1.0 / sqrt(1.0 + t * t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * A[p][q];
          //
          //  Accumulate corrections to diagonal elements.
          //
          zw[p] = zw[p] - h;
          zw[q] = zw[q] + h;
          d[p] = d[p] - h;
          d[q] = d[q] + h;

          A[p][q] = 0.0;
          //
          //  Rotate, using information from the upper triangle of A only.
          //
          for (j = 0; j < p; j++) {
            g = A[j][p];
            h = A[j][q];
            A[j][p] = g - s * (h + g * tau);
            A[j][q] = h + s * (g - h * tau);
          }

          for (j = p + 1; j < q; j++) {
            g = A[p][j];
            h = A[j][q];
            A[p][j] = g - s * (h + g * tau);
            A[j][q] = h + s * (g - h * tau);
          }

          for (j = q + 1; j < n; j++) {
            g = A[p][j];
            h = A[q][j];
            A[p][j] = g - s * (h + g * tau);
            A[q][j] = h + s * (g - h * tau);
          }
          //
          //  Accumulate information in the eigenvector matrix.
          //
          for (j = 0; j < n; j++) {
            g = eigs[j][p];
            h = eigs[j][q];
            eigs[j][p] = g - s * (h + g * tau);
            eigs[j][q] = h + s * (g - h * tau);
          }
          rot_num = rot_num + 1;
        }
      }
    }

    for (i = 0; i < n; i++) {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
  //
  //  Restore upper triangle of input matrix.
  //
  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      A[i][j] = A[j][i];
    }
  }
  //
  //  Ascending sort the eigenvalues and eigenvectors.
  //
  for (k = 0; k < n - 1; k++) {
    m = k;
    for (l = k + 1; l < n; l++) {
      if (d[l] < d[m]) {
        m = l;
      }
    }

    if (m != k) {
      t = d[m];
      d[m] = d[k];
      d[k] = t;
      for (i = 0; i < n; i++) {
        w = eigs[i][m];
        eigs[i][m] = eigs[i][k];
        eigs[i][k] = w;
      }
    }
  }

  return make_tuple(eigs, d);
}

template <class X>
Tuple<Vector<X>, X> power_eigs(const Matrix<X> &matrix, const unsigned itmax) {
  const unsigned n = matrix.row_size();
  Matrix<X> A = matrix;
  Vector<X> lambda = Vector<X>::one(n);
  Vector<X> orthogonal = A * lambda;
  bool orthogonal_eq_0 = true;
  for (unsigned i = 0; i < n; ++i) {
    if (orthogonal[i] != 0) {
      orthogonal_eq_0 = false;
      break;
    }
  }
  if (orthogonal_eq_0)
    throw std::runtime_error("Initial guess of power_eigs is orthogonal!");

  for (unsigned i = 0; i < itmax; ++i) {
    lambda = (A * lambda);
    lambda /= norm(lambda);
  }
  return make_tuple(lambda, norm2(A * lambda));
}

template <class X> X norm2(const Vector<X> &V) {
  const unsigned n = V.size();
  X sum = 0;
  for (unsigned i = 0; i < n; ++i)
    sum += V[i] * V[i];
  return sqrt(sum);
}

template <class X> X norm1(const Vector<X> &V) {
  const unsigned n = V.size();
  X sum = 0;
  for (unsigned i = 0; i < n; ++i)
    sum += abs(V[i]);
  return sum;
}

template <class X>
Tuple<Vector<X>, X> inverse_power_eigs(const Matrix<X> &matrix, const X &mu,
                                       const unsigned itmax) {
  const X rtol = static_cast<X>(std::numeric_limits<double>::epsilon());
  const unsigned n = matrix.row_size();
  Matrix<X> A;
  Matrix<X> I = Matrix<X>::zero(n, n);
  Vector<X> lambda = Vector<X>::one(n);
  Vector<X> orthogonal;
  bool orthogonal_eq_0 = true;
  X ck;
  Vector<X> Y;
  unsigned it = 0;

  bool is_null = true;
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) {
      if (matrix[i][j] != 0) {
        is_null = false;
        break;
      }
    }
  }
  if (is_null) //<! Matrix is null, then e_1 is the smallest and 0 is the
               //eigenvalue
  {
    Vector<X> null_ret = Vector<X>::zero(n);
    null_ret[0] = static_cast<X>(1.0);
    return make_tuple(null_ret, static_cast<X>(0));
  }

  for (unsigned i = 0; i < n; ++i)
    I[i][i] = 1.0;

  A = inverse(matrix - (mu * I));
  orthogonal = A * lambda;

  for (unsigned i = 0; i < n; ++i) {
    if (orthogonal[i] != 0) {
      orthogonal_eq_0 = false;
      break;
    }
  }
  if (orthogonal_eq_0)
    throw std::runtime_error("Initial guess of power_eigs is orthogonal!");

  Y = A * lambda;
  ck = (transpose(Y) * lambda) / (transpose(lambda) * lambda);
  X ck_1;
  do {
    it++;
    lambda = Y / norm2(Y);
    Y = A * lambda;
    ck_1 = ck;
    ck = (transpose(Y) * lambda) / (transpose(lambda) * lambda);
  } while (abs(ck - ck_1) > rtol && it < itmax);
  return make_tuple(lambda, 1 / ck + mu);
}

// Eigen computation of eigenvalues
template <class X> Tuple<Vector<X>, X> eigen_eigs(const Matrix<X> &A) {
  const unsigned n = A.row_size();
  const unsigned m = A.column_size();
  Eigen::MatrixXd eigenA(n, m);
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < m; ++j)
      eigenA(i, j) = static_cast<double>(A[i][j].dbl);
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(
      eigenA, Eigen::ComputeEigenvectors);
  auto eigenValues = solver.eigenvalues();
  auto eigenVectors = solver.eigenvectors();
  unsigned index = 0;
  for (unsigned i = 0; i < eigenValues.size(); ++i) {
    if (eigenValues(i) < eigenValues(index))
      index = i;
  }
  Vector<X> uRet(n);
  X vRet = static_cast<X>(eigenValues(index));
  for (unsigned i = 0; i < n; ++i)
    uRet[i] = static_cast<X>(eigenVectors(i, index));

  return make_tuple(uRet, vRet);
}

template <class X> Tuple<Matrix<X>, unsigned> eigen_null(const Matrix<X> &G) {
  const unsigned n = G.row_size();
  const unsigned m = G.column_size();

  Eigen::MatrixXd A(n, m);

  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < m; ++j)
      A(i, j) = static_cast<double>(G[i][j].dbl);
  }
  Eigen::FullPivLU<Eigen::MatrixXd> lu(A);

  Eigen::MatrixXd A_null_space = lu.kernel();
<<<<<<< HEAD
  const unsigned n_null = A_null_space.rows();
  const unsigned m_null = A_null_space.cols();
  Matrix<X> G_null_space(n_null, m_null);
  bool infeasible = true;
  for (unsigned i = 0; i < n_null; ++i) {
    for (unsigned j = 0; j < m_null; ++j)
      G_null_space[i][j] = static_cast<X>(A_null_space(i, j));
    if (infeasible && G_null_space[i][0] != 0)
      infeasible = false;
  }

  // lu.rank -> FIX to A(i,:)=[0...0] bug!
  return make_tuple(G_null_space,
                    (infeasible) ? 0u : static_cast<unsigned>(lu.rank()));
=======
  const unsigned  n_null = A_null_space.rows();
  const unsigned  m_null = A_null_space.cols();
  Matrix<X>       G_null_space(n_null,m_null);
  bool            infeasible=true;
  for(unsigned i = 0; i<n_null;++i)
  {
    for(unsigned j=0;j<m_null;++j)
      G_null_space[i][j]=static_cast<X>(A_null_space(i,j));
    if(infeasible && G_null_space[i][0]!=0)
      infeasible=false;
  }

  // lu.rank -> FIX to A(i,:)=[0...0] bug!
  return make_tuple(G_null_space, (infeasible)?0u:static_cast<unsigned>(lu.rank()));
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.
}

template <class X> Matrix<X> eigen_pinv(const Matrix<X> &G) {
  const unsigned n = G.row_size();
  const unsigned m = G.column_size();
  Eigen::MatrixXd A(n, m);
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < m; ++j)
      A(i, j) = static_cast<double>(G[i][j].dbl);
  }
  Eigen::MatrixXd pinv = A.completeOrthogonalDecomposition().pseudoInverse();

  const unsigned nRows = pinv.rows();
  const unsigned nCols = pinv.cols();
  Matrix<X> ret(nRows, nCols);
  for (unsigned i = 0; i < nRows; ++i) {
    for (unsigned j = 0; j < nCols; ++j)
      ret[i][j] = static_cast<X>(pinv(i, j));
  }

  return ret;
}

// compute the Cholesky factorization
template <class X> Matrix<X> eigen_chol(const Matrix<X> &A) {
  const unsigned nRow = A.row_size();
  const unsigned nCol = A.column_size();

  Eigen::MatrixXd eG(nRow, nCol);
  Eigen::MatrixXd eL;

  Matrix<X> L(nRow, nCol);

  for (unsigned i = 0; i < nRow; ++i) {
    for (unsigned j = 0; j < nCol; ++j)
      eG(i, j) = static_cast<double>(A[i][j].dbl);
  }

  Eigen::LLT<Eigen::MatrixXd> eLLTA(eG);
  eL = eLLTA.matrixL();

  for (unsigned i = 0; i < nRow; ++i) {
    for (unsigned j = 0; j < nCol; ++j)
      L[i][j] = static_cast<X>(eL(j, i));
  }

  return L;
}
<<<<<<< HEAD
} // namespace Ariadne
=======
} //namesapce Ariadne
>>>>>>> Small fixes. Implemented temporary __feasible__ function to test barrier method.
