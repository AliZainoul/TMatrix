#ifndef ABSTRACT_MAT_C_H
#define ABSTRACT_MAT_C_H
#include "error.hpp"
template <class T> class Vector;   // Forward declaration of class Vector
template <class T> class AbsMtx    // Base Matrix, an Abstract Class
{
protected:
  int nrows;    // Number of rows in the matrix

public:
  // Matrix-Vector Product
  virtual Vector<T> operator* (const Vector<T>&) const = 0;

  // Conjuguate Gradient Method
  int CG(Vector<T>& x, const Vector<T>& b, T& eps, int& iter, int pn=0);
  // Preconditioned Conjugate Gradient Method for Ax=b.
  // A: symmetric positive definite
  // x:  On Entry:  Initial guess; On return: approximate solution
  // b:  Right Side Vector
  // eps: On Entry: (stopping criterion, epsilon)
  //  On return: absolute residual in two-norm for approximate solution
  // iter: On Entry:  max number of iterations allowed;
  //  On return: actual number of iterations taken.
  // pn: =0 if no preconditioner, =1 if diag precr, =2 if SSOR precr
  // it returns 0 for sucessful return and 1 for breakdowns
};

// ************************
//  Implementation
// ************************
template<class T>
int AbsMtx<T>::CG(Vector<T>& x, const Vector<T>& b, T& eps, int& iter, int prec)
{
  if (nrows != b.size()) error("Matrix and Vector sizes do not match.");
  const int maxiter = iter;
  Vector<T> r = b - (*this)*x;                 // Initial residual
  //Vector<T> z = preconding(r,prec);          // Solve the precond system
  Vector<T> z = r;
  Vector<T> p = z;                             // p: search direction
  T zr = dot(z,r);                             // inner prod of z and r
  const double stp = eps*b.twonorm();          // stopping criterion

  if (r.maxnorm() == 0.0) {                  // If intial guess is true solution
    eps = 0.0;                               // return. Otherwise division by
    iter = 0;                                // Division by Zero would occur.
    return 0;
  }

  for (iter = 0; iter < maxiter; iter++)      // Main loop of CG method
  {
    Vector<T> mp = (*this)*p;                 // one matrix-vector multiply
    T pap = dot(mp,p);                        // One of two inner products
    T alpha = zr/pap;                         // pap is 0 only when r is 0
    x += alpha*p;                             // Update the iterative soln
    r -= alpha*mp;                            // Update residual
    if (r.twonorm() <= stp) break;            // Stop if convergence achieved
    // z = preconding(r,prec);
	  z = r;                                    // Preconditioning
    T zrold = zr;
    zr = dot(z,r);                            // Another of two inner products
    T beta= zr/zrold;                         // zrold = 0 only when r is 0
    p = z + beta*p;                           // Update search direction
  }
  eps = r.twonorm();
  if (iter == maxiter) return 1;
  else return 0;
}
// ************************

#endif
