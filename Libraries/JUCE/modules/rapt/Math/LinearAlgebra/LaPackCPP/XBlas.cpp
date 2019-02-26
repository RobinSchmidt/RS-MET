
namespace LaPackCPP {

void blas_dgbmv_x(enum blas_order_type order, enum blas_trans_type trans,
  int m, int n, int kl, int ku, double alpha,
  const double *a, int lda, const double *x, int incx,
  double beta, double *y, int incy, enum blas_prec_type prec)
{
  // LOCAL VARIABLES 
  //
  // As an example, these variables are described on the mxn, column *  major, banded matrix 
  // described in section 2.2.3 of the specification  
  // astart      indexes first element in A where computation begins
  // incai1      indexes first element in row where row is less than lbound
  // incai2      indexes first element in row where row exceeds lbound
  // lbound      denotes the number of rows before  first element shifts 
  // rbound      denotes the columns where there is blank space
  // ra          index of the rightmost element for a given row
  // la          index of leftmost  elements for a given row
  // ra - la     width of a row
  //                         rbound 
  //             la   ra    ____|_____ 
  //              |    |   |          |
  //          |  a00  a01   *    *   *
  //  lbound -|  a10  a11  a12   *   *
  //          |  a20  a21  a22  a23  *
  //              *   a31  a32  a33 a34
  //              *    *   a42  a43 a44
  // Varations on order and transpose have been implemented by modifying these local variables. 

  static const char routine_name[] = "BLAS_dgbmv_x";

  switch(prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
  {
    int ky, iy, kx, jx, j, i, rbound, lbound, ra, la, lenx, leny;
    int incaij, aij, incai1, incai2, astart, ai;
    double *y_i = y;
    const double *a_i = a;
    const double *x_i = x;
    double alpha_i = alpha;
    double beta_i = beta;
    double tmp1;
    double tmp2;
    double result;
    double sum;
    double prod;
    double a_elem;
    double x_elem;
    double y_elem;

    if(order != blas_colmajor && order != blas_rowmajor)
      BLAS_error(routine_name, -1, order, NULL);
    if(trans != blas_no_trans &&
      trans != blas_trans && trans != blas_conj_trans) {
      BLAS_error(routine_name, -2, trans, NULL);
    }
    if(m < 0)
      BLAS_error(routine_name, -3, m, NULL);
    if(n < 0)
      BLAS_error(routine_name, -4, n, NULL);
    if(kl < 0 || kl >= m)
      BLAS_error(routine_name, -5, kl, NULL);
    if(ku < 0 || ku >= n)
      BLAS_error(routine_name, -6, ku, NULL);
    if(lda < kl + ku + 1)
      BLAS_error(routine_name, -9, lda, NULL);
    if(incx == 0)
      BLAS_error(routine_name, -11, incx, NULL);
    if(incy == 0)
      BLAS_error(routine_name, -14, incy, NULL);

    if((m == 0) || (n == 0) || (((alpha_i == 0.0) && (beta_i == 1.0))))
      return;

    if(trans == blas_no_trans) {
      lenx = n;
      leny = m;
    }
    else {			/* change back */
      lenx = m;
      leny = n;
    }

    if(incx < 0) {
      kx = -(lenx - 1) * incx;
    }
    else {
      kx = 0;
    }

    if(incy < 0) {
      ky = -(leny - 1) * incy;
    }
    else {
      ky = 0;
    }

    /* if alpha = 0, return y = y*beta */
    if((order == blas_colmajor) && (trans == blas_no_trans)) {
      astart = ku;
      incai1 = 1;
      incai2 = lda;
      incaij = lda - 1;
      lbound = kl;
      rbound = n - ku - 1;
      ra = ku;
    }
    else if((order == blas_colmajor) && (trans != blas_no_trans)) {
      astart = ku;
      incai1 = lda - 1;
      incai2 = lda;
      incaij = 1;
      lbound = ku;
      rbound = m - kl - 1;
      ra = kl;
    }
    else if((order == blas_rowmajor) && (trans == blas_no_trans)) {
      astart = kl;
      incai1 = lda - 1;
      incai2 = lda;
      incaij = 1;
      lbound = kl;
      rbound = n - ku - 1;
      ra = ku;
    }
    else {			/* rowmajor and blas_trans */
      astart = kl;
      incai1 = 1;
      incai2 = lda;
      incaij = lda - 1;
      lbound = ku;
      rbound = m - kl - 1;
      ra = kl;
    }

    la = 0;
    ai = astart;
    iy = ky;
    for(i = 0; i < leny; i++) {
      sum = 0.0;
      aij = ai;
      jx = kx;

      for(j = ra - la; j >= 0; j--) {
        x_elem = x_i[jx];
        a_elem = a_i[aij];
        prod = x_elem * a_elem;
        sum = sum + prod;
        aij += incaij;
        jx += incx;
      }


      tmp1 = sum * alpha_i;
      y_elem = y_i[iy];
      tmp2 = beta_i * y_elem;
      result = tmp1 + tmp2;
      y_i[iy] = result;
      iy += incy;
      if(i >= lbound) {
        kx += incx;
        ai += incai2;
        la++;
      }
      else {
        ai += incai1;
      }
      if(i < rbound) {
        ra++;
      }
    }


  }
  break;
  case blas_prec_extra:
  {
    int ky, iy, kx, jx, j, i, rbound, lbound, ra, la, lenx, leny;
    int incaij, aij, incai1, incai2, astart, ai;
    double *y_i = y;
    const double *a_i = a;
    const double *x_i = x;
    double alpha_i = alpha;
    double beta_i = beta;
    double head_tmp1, tail_tmp1;
    double head_tmp2, tail_tmp2;
    double result;
    double head_sum, tail_sum;
    double head_prod, tail_prod;
    double a_elem;
    double x_elem;
    double y_elem;
    FPU_FIX_DECL;

    if(order != blas_colmajor && order != blas_rowmajor)
      BLAS_error(routine_name, -1, order, NULL);
    if(trans != blas_no_trans &&
      trans != blas_trans && trans != blas_conj_trans) {
      BLAS_error(routine_name, -2, trans, NULL);
    }
    if(m < 0)
      BLAS_error(routine_name, -3, m, NULL);
    if(n < 0)
      BLAS_error(routine_name, -4, n, NULL);
    if(kl < 0 || kl >= m)
      BLAS_error(routine_name, -5, kl, NULL);
    if(ku < 0 || ku >= n)
      BLAS_error(routine_name, -6, ku, NULL);
    if(lda < kl + ku + 1)
      BLAS_error(routine_name, -9, lda, NULL);
    if(incx == 0)
      BLAS_error(routine_name, -11, incx, NULL);
    if(incy == 0)
      BLAS_error(routine_name, -14, incy, NULL);

    if((m == 0) || (n == 0) || (((alpha_i == 0.0) && (beta_i == 1.0))))
      return;

    if(trans == blas_no_trans) {
      lenx = n;
      leny = m;
    }
    else {			/* change back */
      lenx = m;
      leny = n;
    }

    if(incx < 0) {
      kx = -(lenx - 1) * incx;
    }
    else {
      kx = 0;
    }

    if(incy < 0) {
      ky = -(leny - 1) * incy;
    }
    else {
      ky = 0;
    }

    FPU_FIX_START;

    /* if alpha = 0, return y = y*beta */
    if((order == blas_colmajor) && (trans == blas_no_trans)) {
      astart = ku;
      incai1 = 1;
      incai2 = lda;
      incaij = lda - 1;
      lbound = kl;
      rbound = n - ku - 1;
      ra = ku;
    }
    else if((order == blas_colmajor) && (trans != blas_no_trans)) {
      astart = ku;
      incai1 = lda - 1;
      incai2 = lda;
      incaij = 1;
      lbound = ku;
      rbound = m - kl - 1;
      ra = kl;
    }
    else if((order == blas_rowmajor) && (trans == blas_no_trans)) {
      astart = kl;
      incai1 = lda - 1;
      incai2 = lda;
      incaij = 1;
      lbound = kl;
      rbound = n - ku - 1;
      ra = ku;
    }
    else {			/* rowmajor and blas_trans */
      astart = kl;
      incai1 = 1;
      incai2 = lda;
      incaij = lda - 1;
      lbound = ku;
      rbound = m - kl - 1;
      ra = kl;
    }

    la = 0;
    ai = astart;
    iy = ky;
    for(i = 0; i < leny; i++) {
      head_sum = tail_sum = 0.0;
      aij = ai;
      jx = kx;

      for(j = ra - la; j >= 0; j--) {
        x_elem = x_i[jx];
        a_elem = a_i[aij];
        {
          /* Compute double_double = double * double. */
          double a1, a2, b1, b2, con;

          con = x_elem * split;
          a1 = con - x_elem;
          a1 = con - a1;
          a2 = x_elem - a1;
          con = a_elem * split;
          b1 = con - a_elem;
          b1 = con - b1;
          b2 = a_elem - b1;

          head_prod = x_elem * a_elem;
          tail_prod =
            (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
        }
        {
          /* Compute double-double = double-double + double-double. */
          double bv;
          double s1, s2, t1, t2;

          /* Add two hi words. */
          s1 = head_sum + head_prod;
          bv = s1 - head_sum;
          s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

          /* Add two lo words. */
          t1 = tail_sum + tail_prod;
          bv = t1 - tail_sum;
          t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

          s2 += t1;

          /* Renormalize (s1, s2)  to  (t1, s2) */
          t1 = s1 + s2;
          s2 = s2 - (t1 - s1);

          t2 += s2;

          /* Renormalize (t1, t2)  */
          head_sum = t1 + t2;
          tail_sum = t2 - (head_sum - t1);
        }
        aij += incaij;
        jx += incx;
      }


      {
        /* Compute double-double = double-double * double. */
        double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

        con = head_sum * split;
        a11 = con - head_sum;
        a11 = con - a11;
        a21 = head_sum - a11;
        con = alpha_i * split;
        b1 = con - alpha_i;
        b1 = con - b1;
        b2 = alpha_i - b1;

        c11 = head_sum * alpha_i;
        c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

        c2 = tail_sum * alpha_i;
        t1 = c11 + c2;
        t2 = (c2 - (t1 - c11)) + c21;

        head_tmp1 = t1 + t2;
        tail_tmp1 = t2 - (head_tmp1 - t1);
      }
      y_elem = y_i[iy];
      {
        /* Compute double_double = double * double. */
        double a1, a2, b1, b2, con;

        con = beta_i * split;
        a1 = con - beta_i;
        a1 = con - a1;
        a2 = beta_i - a1;
        con = y_elem * split;
        b1 = con - y_elem;
        b1 = con - b1;
        b2 = y_elem - b1;

        head_tmp2 = beta_i * y_elem;
        tail_tmp2 = (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
      }
      {
        /* Compute double-double = double-double + double-double. */
        double bv;
        double s1, s2, t1, t2;

        /* Add two hi words. */
        s1 = head_tmp1 + head_tmp2;
        bv = s1 - head_tmp1;
        s2 = ((head_tmp2 - bv) + (head_tmp1 - (s1 - bv)));

        /* Add two lo words. */
        t1 = tail_tmp1 + tail_tmp2;
        bv = t1 - tail_tmp1;
        t2 = ((tail_tmp2 - bv) + (tail_tmp1 - (t1 - bv)));

        s2 += t1;

        /* Renormalize (s1, s2)  to  (t1, s2) */
        t1 = s1 + s2;
        s2 = s2 - (t1 - s1);

        t2 += s2;

        /* Renormalize (t1, t2)  */
        result = t1 + t2;
      }
      y_i[iy] = result;
      iy += incy;
      if(i >= lbound) {
        kx += incx;
        ai += incai2;
        la++;
      }
      else {
        ai += incai1;
      }
      if(i < rbound) {
        ra++;
      }
    }

    FPU_FIX_STOP;
  }
  break;
  }
}	/* end GEMV_NAME(d, d, d, _x) */

//-------------------------------------------------------------------------------------------------

void blas_dgbmv2_x(enum blas_order_type order, enum blas_trans_type trans, int m, int n, int kl,
  int ku, double alpha, const double *a, int lda, const double *head_x, const double *tail_x,
  int incx, double beta, double *y, int incy, enum blas_prec_type prec)
{
  static const char routine_name[] = "BLAS_dgbmv2_x";

  switch(prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
  {
    int iy0, iy, ix0, jx, j, i, rbound, lbound, ra, la, lenx, leny;
    int incaij, aij, incai1, incai2, astart, ai;
    double *y_i = y;
    const double *a_i = a;
    const double *head_x_i = head_x;
    const double *tail_x_i = tail_x;
    double alpha_i = alpha;
    double beta_i = beta;
    double tmp1;
    double tmp2;
    double tmp3;
    double tmp4;
    double result;
    double sum1;
    double sum2;
    double prod;
    double a_elem;
    double x_elem;
    double y_elem;


    if(order != blas_colmajor && order != blas_rowmajor)
      BLAS_error(routine_name, -1, order, NULL);
    if(trans != blas_no_trans &&
      trans != blas_trans && trans != blas_conj_trans) {
      BLAS_error(routine_name, -2, trans, NULL);
    }
    if(m < 0)
      BLAS_error(routine_name, -3, m, NULL);
    if(n < 0)
      BLAS_error(routine_name, -4, n, NULL);
    if(kl < 0 || kl >= m)
      BLAS_error(routine_name, -5, kl, NULL);
    if(ku < 0 || ku >= n)
      BLAS_error(routine_name, -6, ku, NULL);
    if(lda < kl + ku + 1)
      BLAS_error(routine_name, -9, lda, NULL);
    if(incx == 0)
      BLAS_error(routine_name, -12, incx, NULL);
    if(incy == 0)
      BLAS_error(routine_name, -15, incy, NULL);

    if(m == 0 || n == 0)
      return;
    if((alpha_i == 0.0) && (beta_i == 1.0))
      return;

    if(trans == blas_no_trans) {
      lenx = n;
      leny = m;
    }
    else {
      lenx = m;
      leny = n;
    }

    ix0 = (incx > 0) ? 0 : -(lenx - 1) * incx;
    iy0 = (incy > 0) ? 0 : -(leny - 1) * incy;



    /* if alpha = 0, return y = y*beta */
    if((order == blas_colmajor) && (trans == blas_no_trans)) {
      astart = ku;
      incai1 = 1;
      incai2 = lda;
      incaij = lda - 1;
      lbound = kl;
      rbound = n - ku - 1;
      ra = ku;
    }
    else if((order == blas_colmajor) && (trans != blas_no_trans)) {
      astart = ku;
      incai1 = lda - 1;
      incai2 = lda;
      incaij = 1;
      lbound = ku;
      rbound = m - kl - 1;
      ra = kl;
    }
    else if((order == blas_rowmajor) && (trans == blas_no_trans)) {
      astart = kl;
      incai1 = lda - 1;
      incai2 = lda;
      incaij = 1;
      lbound = kl;
      rbound = n - ku - 1;
      ra = ku;
    }
    else {			/* rowmajor and blas_trans */
      astart = kl;
      incai1 = 1;
      incai2 = lda;
      incaij = lda - 1;
      lbound = ku;
      rbound = m - kl - 1;
      ra = kl;
    }









    la = 0;
    ai = astart;
    iy = iy0;
    for(i = 0; i < leny; i++) {
      sum1 = 0.0;
      sum2 = 0.0;
      aij = ai;
      jx = ix0;

      for(j = ra - la; j >= 0; j--) {
        x_elem = head_x_i[jx];
        a_elem = a_i[aij];
        prod = x_elem * a_elem;
        sum1 = sum1 + prod;
        x_elem = tail_x_i[jx];
        prod = x_elem * a_elem;
        sum2 = sum2 + prod;
        aij += incaij;
        jx += incx;
      }


      tmp1 = sum1 * alpha_i;
      tmp2 = sum2 * alpha_i;
      tmp3 = tmp1 + tmp2;
      y_elem = y_i[iy];
      tmp4 = beta_i * y_elem;
      result = tmp4 + tmp3;
      y_i[iy] = result;

      iy += incy;
      if(i >= lbound) {
        ix0 += incx;
        ai += incai2;
        la++;
      }
      else {
        ai += incai1;
      }
      if(i < rbound) {
        ra++;
      }
    }


  }
  break;
  case blas_prec_extra:
  {
    int iy0, iy, ix0, jx, j, i, rbound, lbound, ra, la, lenx, leny;
    int incaij, aij, incai1, incai2, astart, ai;
    double *y_i = y;
    const double *a_i = a;
    const double *head_x_i = head_x;
    const double *tail_x_i = tail_x;
    double alpha_i = alpha;
    double beta_i = beta;
    double head_tmp1, tail_tmp1;
    double head_tmp2, tail_tmp2;
    double head_tmp3, tail_tmp3;
    double head_tmp4, tail_tmp4;
    double result;
    double head_sum1, tail_sum1;
    double head_sum2, tail_sum2;
    double head_prod, tail_prod;
    double a_elem;
    double x_elem;
    double y_elem;
    FPU_FIX_DECL;

    if(order != blas_colmajor && order != blas_rowmajor)
      BLAS_error(routine_name, -1, order, NULL);
    if(trans != blas_no_trans &&
      trans != blas_trans && trans != blas_conj_trans) {
      BLAS_error(routine_name, -2, trans, NULL);
    }
    if(m < 0)
      BLAS_error(routine_name, -3, m, NULL);
    if(n < 0)
      BLAS_error(routine_name, -4, n, NULL);
    if(kl < 0 || kl >= m)
      BLAS_error(routine_name, -5, kl, NULL);
    if(ku < 0 || ku >= n)
      BLAS_error(routine_name, -6, ku, NULL);
    if(lda < kl + ku + 1)
      BLAS_error(routine_name, -9, lda, NULL);
    if(incx == 0)
      BLAS_error(routine_name, -12, incx, NULL);
    if(incy == 0)
      BLAS_error(routine_name, -15, incy, NULL);

    if(m == 0 || n == 0)
      return;
    if((alpha_i == 0.0) && (beta_i == 1.0))
      return;

    if(trans == blas_no_trans) {
      lenx = n;
      leny = m;
    }
    else {
      lenx = m;
      leny = n;
    }

    ix0 = (incx > 0) ? 0 : -(lenx - 1) * incx;
    iy0 = (incy > 0) ? 0 : -(leny - 1) * incy;

    FPU_FIX_START;

    /* if alpha = 0, return y = y*beta */
    if((order == blas_colmajor) && (trans == blas_no_trans)) {
      astart = ku;
      incai1 = 1;
      incai2 = lda;
      incaij = lda - 1;
      lbound = kl;
      rbound = n - ku - 1;
      ra = ku;
    }
    else if((order == blas_colmajor) && (trans != blas_no_trans)) {
      astart = ku;
      incai1 = lda - 1;
      incai2 = lda;
      incaij = 1;
      lbound = ku;
      rbound = m - kl - 1;
      ra = kl;
    }
    else if((order == blas_rowmajor) && (trans == blas_no_trans)) {
      astart = kl;
      incai1 = lda - 1;
      incai2 = lda;
      incaij = 1;
      lbound = kl;
      rbound = n - ku - 1;
      ra = ku;
    }
    else {			/* rowmajor and blas_trans */
      astart = kl;
      incai1 = 1;
      incai2 = lda;
      incaij = lda - 1;
      lbound = ku;
      rbound = m - kl - 1;
      ra = kl;
    }









    la = 0;
    ai = astart;
    iy = iy0;
    for(i = 0; i < leny; i++) {
      head_sum1 = tail_sum1 = 0.0;
      head_sum2 = tail_sum2 = 0.0;
      aij = ai;
      jx = ix0;

      for(j = ra - la; j >= 0; j--) {
        x_elem = head_x_i[jx];
        a_elem = a_i[aij];
        {
          /* Compute double_double = double * double. */
          double a1, a2, b1, b2, con;

          con = x_elem * split;
          a1 = con - x_elem;
          a1 = con - a1;
          a2 = x_elem - a1;
          con = a_elem * split;
          b1 = con - a_elem;
          b1 = con - b1;
          b2 = a_elem - b1;

          head_prod = x_elem * a_elem;
          tail_prod =
            (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
        }
        {
          /* Compute double-double = double-double + double-double. */
          double bv;
          double s1, s2, t1, t2;

          /* Add two hi words. */
          s1 = head_sum1 + head_prod;
          bv = s1 - head_sum1;
          s2 = ((head_prod - bv) + (head_sum1 - (s1 - bv)));

          /* Add two lo words. */
          t1 = tail_sum1 + tail_prod;
          bv = t1 - tail_sum1;
          t2 = ((tail_prod - bv) + (tail_sum1 - (t1 - bv)));

          s2 += t1;

          /* Renormalize (s1, s2)  to  (t1, s2) */
          t1 = s1 + s2;
          s2 = s2 - (t1 - s1);

          t2 += s2;

          /* Renormalize (t1, t2)  */
          head_sum1 = t1 + t2;
          tail_sum1 = t2 - (head_sum1 - t1);
        }
        x_elem = tail_x_i[jx];
        {
          /* Compute double_double = double * double. */
          double a1, a2, b1, b2, con;

          con = x_elem * split;
          a1 = con - x_elem;
          a1 = con - a1;
          a2 = x_elem - a1;
          con = a_elem * split;
          b1 = con - a_elem;
          b1 = con - b1;
          b2 = a_elem - b1;

          head_prod = x_elem * a_elem;
          tail_prod =
            (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
        }
        {
          /* Compute double-double = double-double + double-double. */
          double bv;
          double s1, s2, t1, t2;

          /* Add two hi words. */
          s1 = head_sum2 + head_prod;
          bv = s1 - head_sum2;
          s2 = ((head_prod - bv) + (head_sum2 - (s1 - bv)));

          /* Add two lo words. */
          t1 = tail_sum2 + tail_prod;
          bv = t1 - tail_sum2;
          t2 = ((tail_prod - bv) + (tail_sum2 - (t1 - bv)));

          s2 += t1;

          /* Renormalize (s1, s2)  to  (t1, s2) */
          t1 = s1 + s2;
          s2 = s2 - (t1 - s1);

          t2 += s2;

          /* Renormalize (t1, t2)  */
          head_sum2 = t1 + t2;
          tail_sum2 = t2 - (head_sum2 - t1);
        }
        aij += incaij;
        jx += incx;
      }


      {
        /* Compute double-double = double-double * double. */
        double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

        con = head_sum1 * split;
        a11 = con - head_sum1;
        a11 = con - a11;
        a21 = head_sum1 - a11;
        con = alpha_i * split;
        b1 = con - alpha_i;
        b1 = con - b1;
        b2 = alpha_i - b1;

        c11 = head_sum1 * alpha_i;
        c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

        c2 = tail_sum1 * alpha_i;
        t1 = c11 + c2;
        t2 = (c2 - (t1 - c11)) + c21;

        head_tmp1 = t1 + t2;
        tail_tmp1 = t2 - (head_tmp1 - t1);
      }
      {
        /* Compute double-double = double-double * double. */
        double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

        con = head_sum2 * split;
        a11 = con - head_sum2;
        a11 = con - a11;
        a21 = head_sum2 - a11;
        con = alpha_i * split;
        b1 = con - alpha_i;
        b1 = con - b1;
        b2 = alpha_i - b1;

        c11 = head_sum2 * alpha_i;
        c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

        c2 = tail_sum2 * alpha_i;
        t1 = c11 + c2;
        t2 = (c2 - (t1 - c11)) + c21;

        head_tmp2 = t1 + t2;
        tail_tmp2 = t2 - (head_tmp2 - t1);
      }
      {
        /* Compute double-double = double-double + double-double. */
        double bv;
        double s1, s2, t1, t2;

        /* Add two hi words. */
        s1 = head_tmp1 + head_tmp2;
        bv = s1 - head_tmp1;
        s2 = ((head_tmp2 - bv) + (head_tmp1 - (s1 - bv)));

        /* Add two lo words. */
        t1 = tail_tmp1 + tail_tmp2;
        bv = t1 - tail_tmp1;
        t2 = ((tail_tmp2 - bv) + (tail_tmp1 - (t1 - bv)));

        s2 += t1;

        /* Renormalize (s1, s2)  to  (t1, s2) */
        t1 = s1 + s2;
        s2 = s2 - (t1 - s1);

        t2 += s2;

        /* Renormalize (t1, t2)  */
        head_tmp3 = t1 + t2;
        tail_tmp3 = t2 - (head_tmp3 - t1);
      }
      y_elem = y_i[iy];
      {
        /* Compute double_double = double * double. */
        double a1, a2, b1, b2, con;

        con = beta_i * split;
        a1 = con - beta_i;
        a1 = con - a1;
        a2 = beta_i - a1;
        con = y_elem * split;
        b1 = con - y_elem;
        b1 = con - b1;
        b2 = y_elem - b1;

        head_tmp4 = beta_i * y_elem;
        tail_tmp4 = (((a1 * b1 - head_tmp4) + a1 * b2) + a2 * b1) + a2 * b2;
      }
      {
        /* Compute double-double = double-double + double-double. */
        double bv;
        double s1, s2, t1, t2;

        /* Add two hi words. */
        s1 = head_tmp4 + head_tmp3;
        bv = s1 - head_tmp4;
        s2 = ((head_tmp3 - bv) + (head_tmp4 - (s1 - bv)));

        /* Add two lo words. */
        t1 = tail_tmp4 + tail_tmp3;
        bv = t1 - tail_tmp4;
        t2 = ((tail_tmp3 - bv) + (tail_tmp4 - (t1 - bv)));

        s2 += t1;

        /* Renormalize (s1, s2)  to  (t1, s2) */
        t1 = s1 + s2;
        s2 = s2 - (t1 - s1);

        t2 += s2;

        /* Renormalize (t1, t2)  */
        result = t1 + t2;
      }
      y_i[iy] = result;

      iy += incy;
      if(i >= lbound) {
        ix0 += incx;
        ai += incai2;
        la++;
      }
      else {
        ai += incai1;
      }
      if(i < rbound) {
        ra++;
      }
    }

    FPU_FIX_STOP;
  }
  break;
  }
}				/* end BLAS_dgbmv2_x */

//-------------------------------------------------------------------------------------------------

void BLAS_error(const char *rname, int iflag, int ival, char *form, ...)
{
#if !defined(CONFIG_USE_XERBLA)
{
  va_list argptr;
  va_start(argptr, form);
  fprintf(stderr, "Error #%d from routine %s:\n", iflag, rname);
  if(form)
    vfprintf(stderr, form, argptr);
  else if(iflag < 0)
    fprintf(stderr,
      "  Parameter number %d to routine %s had the illegal value %d\n",
      -iflag, rname, ival);
  else
    fprintf(stderr, "  Unknown error code %d from routine %s\n",
      iflag, rname);
  exit(iflag);
}
#else
{
  int ln, argno;
  ln = strlen(rname);
  argno = -iflag;
  xerbla_array(rname, &ln, &argno);
}
#endif
}


}