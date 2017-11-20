#ifndef RS_POLYNOMIAL_H
#define RS_POLYNOMIAL_H

namespace RSLib
{

  /**

  This is a class for representing polynomials

  */

  template<class T>
  class rsPolynomial
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. You have to pass the array of coefficients in coeffs which has to be of
    length order+1 */
    rsPolynomial(T *coeffs, int order);

    /** Copy constructor. */
    rsPolynomial(const rsPolynomial& other);

    /** Destructor. */
    ~rsPolynomial();


    /** \name Operators */

    rsPolynomial& operator=(const rsPolynomial& other);

    rsPolynomial operator-() const;

    bool operator==(const rsPolynomial& other) const;
    bool operator!=(const rsPolynomial& other) const;

    rsPolynomial operator+(const rsPolynomial& other);
    rsPolynomial operator-(const rsPolynomial& other);
    rsPolynomial operator*(const rsPolynomial& other);
    rsPolynomial operator/(const rsPolynomial& other);

    rsPolynomial& operator+=(const rsPolynomial& other);
    rsPolynomial& operator-=(const rsPolynomial& other);
    rsPolynomial& operator*=(const rsPolynomial& other);
    rsPolynomial& operator/=(const rsPolynomial& other);

  protected:

    //void setCoefficients(T *coeffs, int order);
    //void allocate();
    //void free();




    /** \name Data */

    int N;  // order
    T  *a;  // array of coefficients of length N+1

  };

  /** Explicit template instantiation, to be used by the rsPow template-function. */
  template<class T>
  rsPolynomial<T> rsUnityValue(rsPolynomial<T> value);

}

#endif
