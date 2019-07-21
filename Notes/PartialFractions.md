We want to derive a practical algorithm for finding the partial fraction expansion of a strictly
proper rational function: R(x) = B(x)/A(x) where A(x) may have repeated poles. The algorithm should
not rely on solving a linear system of equations - instead, we want to generalize the Heaviside 
cover-up method (https://en.wikipedia.org/wiki/Heaviside_cover-up_method) suitably, so as to allow 
also for repeated poles. The residue method 
(https://en.wikipedia.org/wiki/Partial_fraction_decomposition#Residue_method) doesn't seem to be
attractive for an algorithm either because it involves taking derivatives of R(x) whose
computation would require invoking the quotient rule, which would make the degrees of involved
polynomials blow up exponentially with the number of derivatives that have to be taken, which
equals the multiplicity of the given pole (although an algorithm may (or may not) possible that
doesn't involve creating the coefficient arrays explicitly (maybe research further, if this is in
fact possible). We start with an example and then generalize by recognizing the patterns.

Consider the rational function:

        B(x)     2       3         1       1    3x^3-8x^2+5x-1     3x^3-8x^2+5x-1 
R(x) = ------ = --- + ------- + ------- + --- = -------------- = ------------------
        A(x)    x-1   (x-1)^2   (x-1)^3   x-2   (x-1)^3*(x-2)    x^4-5x^3+9x^2-7x+2

with a triple pole at x=1 and a simple pole at x=2. We express R(x) as a partial fraction expansion
with as-of-yet unknown residues:

3x^3-8x^2+5x-1   r11     r12       r13     r21
-------------- = --- + ------- + ------- + ---             (Eq.1)
(x-1)^3*(x-2)    x-1   (x-1)^2   (x-1)^3   x-2

The residue r21 for the pole at x=2 can be easily found by the regular cover-up method. We are 
interested in the resiudes for the triple pole. The residue r13 for the term with the highest
power of x-1 can also still be found by the regular cover-up method - namely by multiplying the
equation by (x-1)^3. This gives:

3x^3-8x^2+5x-1                                    r21*(x-1)^3
-------------- =  r11*(x-1)^2 + r12*(x-1) + r13 + -----------
    (x-2)                                             x-2

This equation must hold for all values of x, so we strategically plug in the pole x=1 to make all 
terms except r13 vanish on the right hand side (rhs):

3-8+5-1   -1
------- = --- = 1 = r13
  1-2     -1

We have found r13. As said, we could find r21 in the same way, by multiplying the equation by
(x-2) - this would cancel all terms except r21 on the rhs. This is all just the regular cover-up 
method, so far. But now let's find r12 - we do this by multiplying the (Eq.1) by (x-1)^2:

3x^3-8x^2+5x-1                      r13    r21*(x-1)^2
-------------- = r11*(x-1) + r12 + ----- + -----------
 (x-1)*(x-2)                       (x-1)       x-2
 
 Let's solve for r12:
 
       3x^3-8x^2+5x-1    r13                r21*(x-1)^2
 r12 = -------------- - ----- - r11*(x-1) - -----------
        (x-1)*(x-2)     (x-1)                   x-2
        
If we now just plug in x=1, the last two terms would vanish as before but the first two terms would 
evaluate to infinity - so we have a problem here. We are dealing with an infinity-minus-infinity 
situation here. Let's transform it into a zero-divided-by-zero by bringing the first and second 
term over the common denominator (x-1)*(x-2):

       3x^3-8x^2+5x-1    r13*(x-2)        3x^3-8x^2+5x-1 - r13*x + r13*2
 r12 = -------------- - ----------- - Z = ------------------------------ - Z
        (x-1)*(x-2)     (x-1)*(x-2)                 (x-1)*(x-2)
        
where the Z has been used for the second two terms which are of no interest because they will 
soon vanish (the letter Z indicates, that this part will eventually become zero, although it isn't 
zero just yet). Note that r13 is already known from our previous step: r13=1 - so we just plug 
that value in and obtain:

       3x^3-8x^2+4x+1       3x^3-8x^2+4x+1
 r12 = -------------- - Z = -------------- - Z             (Eq.2)
        (x-1)*(x-2)            x^2-3x+2

If we now plug in x=1, we arrive at a 0/0 expression and Z would vanish. That screams for taking 
the limit as x approaches 1 and invoking l'Hospital's rule. So let's take derivatives of numerator
and denominator separately and *then* plug in x=1:

      9x^2-16x+4|      9-16+4   -3
r12 = ----------|    = ------ = --- = 3
         2x-3   |x=1    2-3     -1

and there we have it: r12=3, which is the correct result! Note that - as opposed to the resiude 
method - we didn't have to take the derivative of the full rational function R(x) (times some 
factor). We only had to take derivatives of polynomials - which is easy - no quotient rule, no 
blowing up of polynomial degrees. However, instead of taking the limit and invoking l'Hospital, we 
could also factor out (x-1) from 3x^3-8x^2+4x+1 in (Eq.2) and cancel it with (x-1) in the 
denominator.

      3x^3-8x^2+4x+1       (x-1)*(3x^2-5x-1)       3x^2-5x-1
r12 = -------------- - Z = ----------------- - Z = --------- - Z
       (x-1)*(x-2)            (x-1)*(x-2)             x-2
       
After canceling the (x-1) factors and plugging in x=1, the "Z" term dutyfully vanishes and we are 
left with:
       
      3x^2-5x-1|      3-5-1   -3
r12 = ---------|    = ----- = --- = 3
         x-2   |x=1    -1     -1

which is again the correct result. Invoking l'Hospital is always justified but for the factoring 
out the (x-1)-factor method, we have assumed that 3x^3-8x^2+4x+1 actually *has* such a factor. 
This turned out to be the case here, but can we assume that to be the case in general? ...i think 
so - it would otherwise be a weird coincidence....but i have no proof for that yet....i think, it 
must be the case because we really need such a pole-zero cancellation - otherwise r12 could not be 
a finite value - todo: work out the details

Now for the last one: r11. We get it by multiplying (Eq.1) by x-1:

3x^3-8x^2+5x-1         r12     r13     r21*(x-1)
-------------- = r11 + --- + ------- + ---------
(x-1)^2*(x-2)          x-1   (x-1)^2      x-2

and solve for r11:

      3x^3-8x^2+5x-1   r12    r13      r21*(x-1)
r11 = -------------- - --- - ------- - ---------
      (x-1)^2*(x-2)    x-1   (x-1)^2      x-2
      
now, we call only the last term Z and bring the r12 and r13 terms both over the common denominator
with the first term:


      3x^3-8x^2+5x-1   r12*(x-1)*(x-2)     r13*(x-2)   
r11 = -------------- - --------------- - ------------- - Z
      (x-1)^2*(x-2)     (x-1)^2*(x-2)    (x-1)^2*(x-2) 
      
and now plug in the already known values r13=1 and r12=3:

      3x^3-8x^2+5x-1 - 3*(x-1)*(x-2) - 1*(x-2)
r11 = ---------------------------------------- - Z
                  (x-1)^2*(x-2)    

      3x^3-11x^2+13x-5
r11 = ---------------- - Z
       x^3-4x^2+5x-2

Now, plugging in x=1 will again give 0/0, using l'Hopital once will again lead to 0/0 and using it 
twice and then plugging in x=1 will give:

      18x-22|      -4
r11 = ------|    = --- = 2
       6x-8 |x=1   -2
       
Again, instead of using l'Hospital, we could have factored out the common factor (x-1)^2 from 
numerator and denominator to get:

      (x-1)^2*(3x-5)       3x-5
r11 = -------------- - Z = ---- - Z
       (x-1)^2*(x-2)       x-2
       
and evaluating it at x=1:

      3x-5|      -2
r11 = ----|    = --- = 2
      x-2 |x=1   -1

gives the same result, as it should. Done. We have found all 3 coeffs for the triple pole at x=1 
without invoking the residue rule and without solving a linear system of equations. All we needed 
to do is to add polynomials and find derivatives of polynomials or do divisions of a polynomial by 
a monomial. It seems that the route via factoring-out method is more attractive algorithmically 
than the derivative route because the involved polynomials tended to have smaller coeffs ...but 
does that actually matter - maybe implement both and do precision- and performance comparisons..

---------------------------------------------------------------------------------------------------

So much for the example - now, we want to generalize what we have done and derive an algorithm. For
the computation for r13, r12, r11 we have used the following formulas:

     B(x)       B(x) - r13*(x-2)       B(x) - r13*(x-2) - r12*(x-1)*(x-2)
r13: ----, r12: ----------------, r11: ----------------------------------
     x-2          (x-2)*(x-1)                   (x-2)*(x-1)^2
     
i'm not using an equals sign here because that wouldn't make sense - the rhs are all still 
functions of x. To compute the rij values, we had to plug the pole at x=1 into these expressions. 
The first one works as is. The second one requires invoking l'Hospital once or cancel a common 
factor of (x-1) in numerator and denominator once. The third requires to take these steps twice. 
We can see a pattern here. To establish a general formula, we introduce the notation:

       B(x)   b0 + b1*x + ... + bM*x^M   b0 + b1*x + ... + bM*x^M                 r_i,j
R(x) = ---- = ------------------------ = ------------------------ = sum_i sum_j ---------
       A(x)   a0 + a1*x + ... +  1*x^N     product_i (x-p_i)^m_i                (x-p_i)^j
       
the index i runs from 1,...,I where I is the number of distinct poles and j in the inner sum runs 
from 1,...,m_i where m_i is the multiplicity of pole i. So, i is the pole index and j is the 
exponent. In the first step above, we wanted to find r13 and we did this by multiplying (Eq.1) by
(x-1)^3. In this case, the 1 comes from the poles index i (i=1) and the 3 is the multiplicity of 
the i-th (i.e. 1st) pole (m_i=3). So, we have multiplied by (x-p_i)^mi. We define the (N-j)th 
degree polynomial A_i,j as:

              A(x)
A_i,j (x) = ---------
            (x-p_i)^j

i claim that the general formula for the residue r_i,j is given by:

        B(x) - r_i,j+1 * A_i,j+1 (x) - r_i,j+2 * A_i,j+2 (x) - ... - r_i,mi * A_i,mi (x) |
r_i,j = ---------------------------------------------------------------------------------|
                                   A_i,j (x)                                             |x=p_i

The formula is recursive in the sense that the computation of r_i,j requires knowledge of all 
higher r_i,j+k coeffs, where j+k runs up to mi (the i-th pole's mutliplicity) - and these coeffs 
are those which were computed in previous steps of the algorithm. The recursion starts with 
computing r_i,mi for which no higher coeffs are needed and then always computes the next lower 
residue. For j=mi, the evaluation of the formula at x=p_i poses no problem, but for j < mi, the 
evaluation at x=p_i will result in a 0/0 expression, so it requires either taking the limit as x 
approaches p_i and then invoking l'Hospital's rule mi-j times or cancelling a common factor of 
(x-p_i)^(mi-j) from numerator and denominator. In the case of cancelling, note that in the 
denominator, we'll always end up at A_i,mi - so we may compute that denominator once. So, we may 
also write:

        [B(x) - r_i,j+1 * A_i,j+1(x) - ... - r_i,mi * A_i,mi (x)] / (x-p_i)^(mi-j) |
r_i,j = ---------------------------------------------------------------------------|
                                   A_i,mi (x)                                      |x=p_i

An unoptimized prototype implementation in C++ (using RAPT::rsPolynomial) could look like this:

std::vector<double> partialFractions(
  const RAPT::rsPolynomial<double>& B, 
  const RAPT::rsPolynomial<double>& A, 
  const std::vector<double>& p, 
  const std::vector<int>& m)  
{
  typedef RAPT::rsPolynomial<double> Poly;
  int N = A.getDegree();                     // degree of denominator = number of residues
  int k, j0 = 0;                             // flat array index into r and base-index
  std::vector<double> r(N);                  // array of residues
  Poly Li, Bij, Aij, Lij, Cij;               // the involved polynomials
  double num, den;                           // numerator and denominator, evaluated at the pole
  for(int i = 0; i < (int)p.size(); i++) {   // loop over distinct poles
    Li  = std::vector<double>{ -p[i], 1.0 }; // linear factor (x - p[i])
    Bij = B;                                 // needed for computation of numerator
    Aij = A / (Li^m[i]);                     // parentheses needed, ^ has lower precedence than /
    den = Aij(p[i]);                         // we need to evaluate this only once
    for(int j = m[i]; j >= 1; j--) {         // loop over the exponents for pole i
      k    = j0+j-1;                         // index into r-array
      Lij  = Li^(m[i]-j);                    // Shlemiel the painter strikes again
      Cij  = Bij / Lij;                      // cancel common factor with denominator
      num  = Cij(p[i]);                      // evaluate numerator
      r[k] = num / den;                      // denominator stays the same inside the j-loop
      Bij  = Bij - Aij * r[k];               // establish B-polynomial for next iteration
      Aij  = Aij * Li;                       // establish A-polynomial for next iteration
    }
    j0 += m[i];                              // increment base-index
  }
  return r;
}

The implementaton above uses the factoring-out approach. A variant using l'Hospital looks like:

std::vector<double> partialFractions2( 
  const RAPT::rsPolynomial<double>& B, 
  const RAPT::rsPolynomial<double>& A, 
  const std::vector<double>& p, 
  const std::vector<int>& m)  
{
  typedef RAPT::rsPolynomial<double> Poly;
  int N = A.getDegree();
  int k, j0 = 0;
  std::vector<double> r(N);
  Poly Li, Bij, Aij, Lij;                        // no Cij needed anymore
  double num, den;
  for(int i = 0; i < (int)p.size(); i++) { 
    Li  = std::vector<double>{ -p[i], 1.0 };
    Bij = B; 
    Aij = A / (Li^m[i]);
    for(int j = m[i]; j >= 1; j--) {
      k    = j0+j-1;
      Lij  = Li^(m[i]-j);
      num  = Bij.derivativeAt(p[i], m[i]-j);
      den  = Aij.derivativeAt(p[i], m[i]-j);     // den must now be evaluated for each j
      r[k] = num / den;
      Bij  = Bij - Aij * r[k];
      Aij  = Aij * Li;
    }
    j0 += m[i]; 
  }
  return r;
} 
