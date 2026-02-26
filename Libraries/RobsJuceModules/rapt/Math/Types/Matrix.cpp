


/*

ToDo:

- Generalize matrix multiplication for matrices of shapes that do not fit togther natrurally by
  zero padding the rows and columns of the arguments as needed. The shape of the product C = A*B 
  should be: #rows(C) = max(#rows(A), #cols(B)), #cols(C) = max(#cols(A), #rows(B)), I think. We 
  first need to zero-pad A,B to the shape #rows(C) x #cols(C) and then we can do the 
  multiplication of two square matrices as usual. Then verify if (A*B)^T = B^T * A^T still holds
  with this extended definition. What about associativity and distributivity over addition? Maybe
  we need to extend the definition of matrix addition, too. Just pad the matrices with zeros. In 
  the case of addition, we need: #rows(C) = max(#rows(A), #rows(B)), 
  #cols(C) = max(#cols(A), #cols(B)), I think.

- Implement a class to represent triangular matrices which may be used with minor modifictions for
  symmetric matrices, too: just instead of returning 0 when j > i, return A(j,i). Maybe 
  antisymmetric matrices can also be represented by a minor modificiation: we need a way to remove 
  the diagonal. Maybe when the user requests a size MxN matrix, we internally use a triangular 
  array of size (M-1)x(N-1) or something. Of course, the formula to compute the flat index also 
  needs modification. See:
  https://en.wikipedia.org/wiki/Triangular_matrix
  https://en.wikipedia.org/wiki/Triangular_array 
  https://en.wikipedia.org/wiki/Category:Triangles_of_numbers






Resources:

- This video:
  https://www.youtube.com/watch?v=PdxWhTX8aKo  Eigenvectors before Eigenvalues | solve(x)
  shows a method for computing the eigenvectors of a 2x2 matrix without computing the eigenvalues 
  first. May that could be useful when we are only interested in the eigenvectors but not in the 
  eigenvalues.

- This video:
  https://www.youtube.com/watch?v=GHctcSBd6Z4
  Matrix Multiplication Deep Dive || Cache Blocking, SIMD & Parallelization - Aliaksei Sala - CppCon
  explains various ways of optimizing matrix multiplication.

*/
