//  lopt.cpp -- Optimum 'L' Filter algorithm.
//  (C) 2004, C. Bond.
//
//  Based on discussion in Kuo, "Network Analysis and Synthesis",
//  pp. 379-383. Original method due to A.Papoulis."On Monotonic
//  Response Filters", Proc. IRE, 47, Feb. 1959.
//
#include <math.h>

//  This routine calculates the coefficients of the Legendre polynomial
//  of the 1st kind. It uses a recursion relation. The first few polynomials
//  are hard coded and the rest are found by recursion.
//
//  (n+1)Pn+1 = (2n+1)xPn - nPn-1 	Recursion relation.
//
void legendre(double *p,int n)
{
    double *a,*b;
    int i,j;

    if (n == 0) {
        p[0] = 1.0;
        return;
    }
    if (n == 1) {
        p[0] = 0.0;
        p[1] = 1.0;
        return;
    }
    p[0] = -0.5;
    p[1] = 0.0;
    p[2] = 1.5;

    if (n == 2) return;

    a = new double [n+1];
    b = new double [n+1];

    for (i=0;i<=n;i++) {
        a[i] = b[i] = 0.0;
    }
    b[1] = 1.0;

    for (i=3;i<=n;i++) {
        for (j=0;j<=i;j++) {
            a[j] = b[j];
            b[j] = p[j];
            p[j] = 0.0;
        }
        for (j=i-2;j>=0;j-=2) {
            p[j] -= (i-1)*a[j]/i;
        }
        for (j=i-1;j>=0;j-=2) {
            p[j+1] += (2*i-1)*b[j]/i;
        }

    }

    delete [] b;
    delete [] a;
}
//
//
//  In the following routine n = 2k + 1 for odd 'n' and n = 2k + 2 for
//  even 'n'.
//
//
//      n   k
//      -----
//      1   0
//      2   0
//      3   1
//      4   1
//      5   2
//      6   2
//

void lopt(double *w,int n)
{
    double *a,*p,*s,*v,c0,c1;
    int i,j,k;

    a = new double [n+1];
    p = new double [2*n+1];
    s = new double [2*n+1];
    v = new double [2*n+4];

    k = (n-1)/2;
//
//  form vector of 'a' constants
//
    if (n & 1) {                // odd
        for (i=0;i<=k;i++) {
            a[i] = (2.0*i+1.0)/(M_SQRT2*(k+1.0));
        }
    }                           // even
    else {
        for (i=0;i<k+1;i++) {
            a[i] = 0.0;
        }
        if (k & 1) {
            for (i=1;i<=k;i+=2) {
                a[i] = (2*i+1)/sqrt((k+1)*(k+2));
            }
        }
        else {
            for (i=0;i<=k;i+=2) {
                a[i] = (2*i+1)/sqrt((k+1)*(k+2));
            }
        }
    }
    for (i=0;i<=n;i++){
        s[i] = 0.0;
        w[i] = 0.0;
    }
//
// form s[] = sum of a[i]*P[i]
//
    s[0] = a[0];
    s[1] = a[1];
    for (i=2;i<=k;i++) {
        legendre(p,i);
        for (j=0;j<=i;j++) {
            s[j] += a[i]*p[j];
        }
    }
//
//  form v[] = square of s[]
//
    for (i=0;i<=2*k+2;i++) {
        v[i] = 0.0;
    }
    for (i=0;i<=k;i++) {
        for (j=0;j<=k;j++) {
            v[i+j] += s[i]*s[j];    
        }
    }
//
//  modify integrand for even 'n'
//
    v[2*k+1] = 0.0;
    if ((n & 1) == 0) {
        for (i=n;i>=0;i--) {
            v[i+1] += v[i];
        }
    }
//
//  form integral of v[]
//
    for (i=n+1;i>=0;i--) {
        v[i+1] = v[i]/(double)(i+1.0);
    }
    v[0] = 0.0;
//
// clear s[] for use in computing definite integral
//
    for (i=0;i<(n+2);i++){ 
        s[i] = 0.0;
    }
    s[0] = -1.0;
    s[1] = 2.0;
//
//  calculate definite integral
//
    for (i=1;i<=n;i++) {
        if (i > 1) {
            c0 = -s[0];
            for (j=1;j<i+1;j++) {
                c1 = -s[j] + 2.0*s[j-1];
                s[j-1] = c0;
                c0 = c1;
            }
            c1 = 2.0*s[i];
            s[i] = c0;
            s[i+1] = c1;
        }
        for (j=i;j>0;j--) {
            w[j] += (v[i]*s[j]);
        }
    }
    if ((n & 1) == 0) w[1] = 0.0;
    delete [] v;
    delete [] p;
    delete [] s;
    delete [] a;
}

#include <iostream.h>

int main()
{
    double w[20];
    int i,n;

    cout << "Order of Optimal (L) Filter (n < 20): ";	// This limit is arbitrary!
    cin >> n;

    if (n > 19) return 1;
    lopt(w,n);
    for (i=1;i<=n;i++) {
        if (w[i])
        cout << w[i] << "w^" << 2*i << endl;
    }
    return 0;
}

