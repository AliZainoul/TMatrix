#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
#include <cassert>
#include <cmath>
using namespace std;
#include "error.hpp"

template <class T>
class Vector
{
 private:
  T* M_data;
  int M_size;

 public:
  // Constructors
  Vector();                           // Constructor 0: Default Constructor
  Vector(int i);                      // Constructor 1: Size Constructor
  Vector(int i, T v);                 // Constructor 2: Size-Value Constructor
  Vector(initializer_list<T> lst);    // Constructor 3: List Constructor
  Vector(Vector const &v);            // Constructor 4: Copy Constructor
  Vector(Vector &&v);                 // Constructor 5: Move constructor

  // Destructors
  ~Vector();

  // void reset(double x);
  void affiche();
  void reset();
  int size() const;    // Return Size of the Vector

  // Norms
  T maxnorm() const;                // L-infinity norm
  T twonorm() const;                // L-2 norm

  // Operators
  Vector<T>& operator=(Vector<T> const &v);    // Assignment Operator
  Vector<T>& operator=(Vector<T>&& v);         // Move Assignment Operator
  Vector<T>& operator+=(const Vector<T>& vec); // Binary Operator, V += W
  Vector<T>& operator-=(const Vector<T>& vec); // Binary Operator, V -= W
  T& operator[](int i){return M_data[i];} // Subscript Operator left value
  T operator[](int i) const{return M_data[i];} // Subscript Operator const

  // Friends Functions
  template <class S>
  friend S dot(const Vector<S>& v1, const Vector<S>& v2);  // Dot product
  template <class S>
  friend Vector<S> operator*(S, const Vector<S>&);   // Scalar-Vector Product
  template <class S>
  friend Vector<S> operator*(const Vector<S>&, S);   // Vector-Scalar Product
  template <class S>
  friend Vector<S> operator+(const Vector<S>&);                 // Unary Operator, V = +W
  template <class S>
  friend Vector<S> operator-(const Vector<S>&);                // Unary Operator, V = -W
  template <class S>
  friend Vector<S> operator+(const Vector<S>&, const Vector<S>&); // Binary Operator V = V1 + V2
  template <class S>
  friend Vector<S> operator-(const Vector<S>&, const Vector<S>&); // Binary Operator V = V1 - V2
  template <class S>
  friend std::ostream& operator<<(std::ostream&, const Vector<S>&); // '<<' Overload
};

// IMPLEMENTATION

// Constructor 0 - Default Constructor
// Usage: Vector V();
template<class T>
Vector<T>::Vector(){
	M_data = nullptr;
	M_size = 0;
}

// Constructor 1 - Size Constructor
// Usage: Vector V(10) => V = [0.0 ; 0.0 ; ... ; 0.0], size=10;
template<class T>
Vector<T>::Vector(int i){
	M_data = new T[M_size=i];
}

// Constructor 2 - Size and value Constructor
// Usage: Vector V(2,3) => V = [3;3]
template<class T>
Vector<T>::Vector(int i, T v){
	M_data = new T[M_size=i];
  for(int k=0;k<M_size;++k) M_data[k] = v;
}

// Constructor 3 - List Constructor
// Usage: Vector V(list) => V = [list.begin; ... ; list.end];
template<class T>
Vector<T>::Vector(initializer_list<T> lst){
	M_data = new T[lst.size()];
	M_size = static_cast<int>(lst.size());
 	copy(lst.begin(), lst.end(), M_data); // std::copy
}

// Constructor 4 - Copy Constructor
// Usage: Vector V(W) => V = W;
template<class T>
Vector<T>::Vector(Vector const &v){
	v.M_size==0? M_data = nullptr : M_data = new T[v.M_size];
	M_size = v.M_size;
	for (int k=0;k<M_size;++k)M_data[k] = v.M_data[k];
}

// Constructor 5 - Move Constructor
// Usage: Vector V(W) => V = W;
template<class T>
Vector<T>::Vector(Vector && v){
	M_data=v.M_data;
	M_size=v.M_size;
	v.M_size = 0;
	v.M_data = nullptr;
}

// Destructor
template<class T>
Vector<T>::~Vector(){delete[] M_data;}

// Show Function
template<class T>
void Vector<T>::affiche(){
  // cout << "Taille du vecteur: "<< M_size<< "\n";
  for(int i=0; i<M_size; i++) cout<<i<<" "<<M_data[i]<<endl;
}

// Show Function
template<class T>
void Vector<T>::reset(){
  for (unsigned int k= 0; k< M_size ; k++) M_data[k]= 0.0;
}



template<class T>
int Vector<T>::size() const { return M_size; }    // Return Size of the Vector

/*Vector Vector::reset(double x) {
for(int i=0; i<M_size; i++) M_data[i]=x;
};*/



// Assignment Operator
template<class T>
Vector<T>& Vector<T>::operator=(Vector<T> const &v){
  if (this != &v){
	  delete[] M_data;
	  M_data = new double[M_size=v.M_size];
	  for(int k=0;k<M_size;++k) M_data[k]=v.M_data[k];
  }
  return *this;
}

// Move Assignment Operator
template<class T>
Vector<T>& Vector<T>::operator=(Vector<T> &&v){
	delete[] M_data;
	M_data = v.M_data;
	M_size = v.M_size;
	v.M_data = nullptr;
	v.M_size = 0;
	return *this;
}

// Surcharge Of Operator 'Outstream'
	template<class T>
	std::ostream& operator<<(std::ostream& s, const Vector<T> &v){
  //s << "Surcharge opérateur: taille du vecteur: " << v.M_size << endl;
  for (int i = 0; i < v.M_size; i++)
	s<< " | "<< "Index : "<< i+1<< " | "<< " "<< v.M_data[i]<<" | "<< endl;
  return s;
}

// L-infinity Norm for Vectors
template<class T>
T Vector<T>::maxnorm() const{
  T max_norm=abs(M_data[0]);
  for (int i = 1; i < M_size; i++){
		T vi = abs(M_data[i]);
		if (max_norm < vi) max_norm = vi;
	}
  return max_norm;
}

// L-2 Norm for Vectors
template<class T>
T Vector<T>::twonorm() const{
  T two_norm=0.0;
  for (int i = 0; i < M_size; i++) two_norm += pow(M_data[i],2);
  return sqrt(two_norm);
}

// Dot Operator
template<class T>
T dot(const Vector<T>& v1, const Vector<T>& v2) {
  if (v1.M_size != v2.M_size ) error("Bad vector sizes in dot product.");
  T tm = v1[0]*v2[0];
  for (int i = 1; i < v1.M_size; i++) tm += v1[i]*v2[i];
  return tm;
}

// Scalar-Vector Product
template<class T>
Vector<T> operator*(T scalar, const Vector<T>& vec) {
  Vector<T> tm(vec.M_size);
  for (int i = 0; i < vec.M_size; i++) tm[i] = scalar*vec[i];
  return tm;
}

// Vector-Scalar Product
template<class T>
Vector<T> operator*(const Vector<T> &vec, T scalar) {
  Vector<T> tm(vec.M_size);
  for (int i = 0; i < vec.M_size; i++) tm[i] = scalar*vec[i];
  return tm;
}

// Unary Operator "+" overload - Usage: Vector V1 = +V2;
template<class T>
Vector<T> operator+(const Vector<T>& vec) {return vec;}

// Unary Operator "-" overload - Usage: Vector V1 = -V2;
template<class T>
Vector<T> operator-(const Vector<T> &vec) {
  Vector<T> minus = vec;
  for (int i = 0; i < vec.M_size; i++) minus[i] = (-1) * vec[i];
  return minus;
}
// Binary Operator "+" overload - Usage: Vector V = V1+V2;
template<class T>
Vector<T> operator+(const Vector<T> &v1, const Vector<T> &v2) {
  if (v1.M_size != v2.M_size ) error("Bad vector sizes in Vector addition.");
  Vector<T> sum = v1;
  sum += v2;
  return sum;
}

// Binary Operator "-" overload - Usage: Vector V = V1-V2;
template<class T>
Vector<T> operator-(const Vector<T> &v1, const Vector<T> &v2) {
  if (v1.M_size != v2.M_size ) error("Bad vector sizes in Vector substraction.");
  Vector<T> sum = v1; // It would cause problem without copy constructor
  sum -= v2;
  return sum;
}
// Unary Operator "+=" overload - Usage: V += W
template<class T>
Vector<T>& Vector<T>::operator+=(const Vector<T>& vec) {
  if (M_size != vec.M_size) error("Bad vector sizes in Vector::operator+=()");
  for (int i = 0; i < M_size; i++) M_data[i] += vec.M_data[i];
  return *this;
}

// Unary Operator "-=" overload - Usage: V -= W
template<class T>
Vector<T>& Vector<T>::operator-=(const Vector<T>& vec) {
  if (M_size != vec.M_size) error("Bad vector sizes in Vector::operator-=()");
  for (int i = 0; i < M_size; i++) M_data[i] -= vec.M_data[i];
  return *this;
}

#endif
