#include "copyright.h"
#include "matrix_ops.h"

namespace stormm {
namespace math {

//-------------------------------------------------------------------------------------------------
void checkArraySizeToMatrixRank(const size_t actual_rank, const size_t rank, const size_t s_a,
                                const size_t s_v, const ExceptionResponse policy) {
  if (actual_rank == 0LLU) {
    rtErr("A matrix of zero rank is invalid.", "jacobiEigensolver");
  }
  if (actual_rank * actual_rank != s_a || actual_rank * actual_rank != s_v) {
    if (s_a < actual_rank * actual_rank) {
      if (rank == 0LLU) {
        rtErr("Inferred rank " + std::to_string(actual_rank) + " suggests that the input matrix "
              "does not contain enough data (" + std::to_string(s_a) + " elements in the "
              "object).", "jacobiEigensolver");
      }
      else {
        rtErr("The input matrix does not contain enough data (" + std::to_string(s_a) +
              " elements) to be of rank " + std::to_string(rank) + ".", "jacobiEigensolver");
      }
    }
    else if (s_v < actual_rank * actual_rank) {
      rtErr("Eigenvector matrix contains too little space to store eigenvectors of a rank " +
            std::to_string(actual_rank) + " matrix.", "jacobiEigensolver");
    }
    else {
      if (rank == 0LLU) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("Input matrix or eigenvector space are too large for a rank " +
                std::to_string(actual_rank) + " matrix.", "jacobiEigensolver");
        case ExceptionResponse::WARN:
          rtWarn("Input matrix or eigenvector space are too large for a rank " +
                 std::to_string(actual_rank) + " matrix.", "jacobiEigensolver");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      else {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("Input matrix or eigenvector space are too large for a matrix of inferred rank " +
                std::to_string(actual_rank) + ".", "jacobiEigensolver");
        case ExceptionResponse::WARN:
          rtWarn("Input matrix or eigenvector space are too large for a matrix of inferred rank " +
                 std::to_string(actual_rank) + ".", "jacobiEigensolver");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
    }
  }
}

// CHECK
//-------------------------------------------------------------------------------------------------
void TRED2(double** A, int n, double* d, double* e)
{
  int i, j, k, l;
  double scale, hh, h, g, f;
  double* tmp;
  double* tm2p;

  for (i = n-1; i >= 1; i--) {
    l = i - 1;
    h = 0.0;
    scale = 0.0;
    tmp = A[i];
    if (l > 0) {
      for (k = 0; k <= l; k++) {
        scale += fabs(tmp[k]);
      }
      if (scale == 0.0) {
        e[i] = tmp[l];
      }
      else {
        for (k = 0; k <= l; k++) {
          tmp[k] /= scale;
          h += tmp[k]*tmp[k];
        }
        f = tmp[l];
        g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale*g;
        h -= f*g;
        tmp[l] = f - g;
        f = 0.0;
        for (j = 0; j <= l; j++) {
          tm2p = A[j];
          tm2p[i] = tmp[j]/h;
          g = 0.0;
          for (k = 0; k <= j; k++) {
            g += tm2p[k]*tmp[k];
          }
          for (k = j+1; k <=l; k++) {
            g += A[k][j]*tmp[k];
          }
          e[j] = g/h;
          f += e[j]*tmp[j];
        }
        hh = f/(h + h);
        for (j = 0; j <= l; j++) {
          f = tmp[j];
          e[j] = g = e[j] - hh*f;
          tm2p = A[j];
          for (k = 0; k <= j; k++) {
            tm2p[k] -= (f*e[k] + g*tmp[k]);
          }
        }
      }
    }
    else {
      e[i] = tmp[l];
    }
    d[i] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;
  for (i = 0; i < n; i++) {
    tmp = A[i];
    l = i - 1;
    if (d[i]) {
      for (j = 0; j <= l; j++) {
        g = 0.0;
        for (k = 0; k <= l; k++) {
          g += tmp[k]*A[k][j];
        }
        for (k = 0; k <= l; k++) {
          A[k][j] -= g*A[k][i];
        }
      }
    }
    d[i] = tmp[i];
    tmp[i] = 1.0;
    for (j = 0; j <= l; j++) {
      A[j][i] = tmp[j] = 0.0;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void TQLI(double* d, double* e, int n, double** z)
{
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  double* tmp;

  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;
  for (l = 0; l < n; l++) {
    iter = 0;
    m = l - 1;
    while (m != l) {
      for (m = l; m < n - 1; m++) {
        dd = fabs(d[m]) + fabs(d[m+1]);
        if (fabs(e[m]+dd) == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) {
          printf("TQLI >> Error: Too many iterations in tqli\n");
        }
        g = (d[l+1]-d[l])/(2.0*e[l]);
        r = hypotenuse<double>(g, 1.0);
        const double sign_result = (g >= 0.0) ? std::abs(r) : std::abs(-r);
        g = d[m]-d[l]+e[l]/(g + sign_result);
        c = 1.0;
        s = 1.0;
        p = 0.0;
        for (i = m-1; i >= l; i--) {
          f = s*e[i];
          b = c*e[i];
          e[i+1] = (r = hypotenuse<double>(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f/r;
          c = g/r;
          g = d[i+1]-p;
          r = (d[i] - g)*s + 2.0*c*b;
          d[i+1] = g + (p = s*r);
          g = c*r - b;
          for (k = 0; k < n; k++) {
            tmp = z[k];
            f = tmp[i+1];
            tmp[i+1] = s*tmp[i] + c*f;
            tmp[i] = c*tmp[i] - s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    }
  }
}
// END CHECK
  
} // namespace math
} // namespace stormm
