#ifndef ABSTRACT_MAT_C_H
#define ABSTRACT_MAT_C_H

#include "error.hpp"


template <class T>
class Vector;  // Forward declaration of class Vector

template <class T>
class AbsMtx    // Base Matrix, an Abstract Class
{
protected:
  int nrows;                // Number of rows in the matrix

public:
  virtual Vector<T> operator*(const Vector<T>&) const = 0;  // matrix vector multiply

  int CG(Vector<T>& x, const Vector<T>& b, double& eps, int& iter, int pn=0);
      // preconditioned Conjugate gradient method for Ax=b. A: sym pos def
      // x:  on entry:  initial guess; on return: approximate solution
      // b:  right side vector
      // eps: on entry:  stopping criterion, epsilon
      //      on return: absolute residual in two-norm for approximate solution
      // iter: on entry:  max number of iterations allowed;
      //       on return: actual number of iterations taken.
      // pn: =0 if no preconditioner, =1 if diag precr, =2 if SSOR precr
      // it returns 0 for sucessful return and 1 for breakdowns
};

// *** definitions of iterative solvers for AbsMtx

template<class T>
int AbsMtx<T>::CG(Vector<T>& x, const Vector<T>& b, double& eps, int& iter, int prec) {
  // Conjugate gradient method.
  // x: on entry, initial guess; on retrun: approximate solution
  // b: right hand side vector as in A x = b;
  // eps:  on entry, tolerance; on retrun: absolute residual in Euclid norm
  // iter: on entry, max number of iterations allowed;
  //       on return, actual number of iterations used
  // prec= 0 if no preconditioning, 1 if diag prec, 2 if SSOR prec

  if (nrows != b.size()) error("matrix and vector sizes do not match");
  const int maxiter = iter;
  Vector<T> r = b - (*this)*x;                   // initial residual
  //Vector<T> z = preconding(r,prec);              // solve the precond system
  Vector<T> z = r;
  Vector<T> p = z;                               // p: search direction
  T zr = dot(z,r);                            // inner prod of z and r
  const double stp = eps*b.twonorm();         // stopping criterion

  if (r.maxnorm() == 0.0) {                   // if intial guess is true soln,
    eps = 0.0;                                // return. Otherwise division by
    iter = 0;                                 // zero would occur.
    return 0;
  }

  for (iter = 0; iter < maxiter; iter++) {     // main loop of CG method
    Vector<T> mp = (*this)*p;                     // one matrix-vector multiply
    T pap = dot(mp,p);                         // one of two inner products
    T alpha = zr/pap;                          // pap is 0 only when r is 0
    x += alpha*p;                              // update the iterative soln
    r -= alpha*mp;                             // update residual
    if (r.twonorm() <= stp) break;             // stop if convergence achieved
    //z = preconding(r,prec);
	  z = r;                   // preconditioning
    T zrold = zr;
    zr = dot(z,r);                             // another of two inner products
    T beta= zr/zrold;                          // zrold = 0 only when r is 0
    p = z + beta*p;                            // update search direction
  }

  eps = r.twonorm();
  if (iter == maxiter) return 1;
  else return 0;
}
#endif
