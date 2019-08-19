// from here https://www.kvraudio.com/forum/viewtopic.php?p=7494722#p7494722
// you can compute filter coefficient easily like that:
// g = tanpi_blt(freq*sample_rate_inverse)
// Error is around 3.6 ULPs


// Approximates tan(x*pi), x on (1e-35, 0.4999999). Useful for BLT prewrap.
// SSE2 instructions, 4-way parallel
inline __m128 tanpi_blt(__m128 x) 
{
  // R(x^2)*x is relative error minmax rational fit to tan(x*pi):
  // R(x) =
  // (-3.13832354545593x+3.31913113594055)/(x^2+-4.47475242614746x+1.05651223659515)
  // ~4 ULP error
  const __m128 num0 = _mm_set1_ps(3.319131136e+00);
  const __m128 num1 = _mm_set1_ps(-3.138323545e+00);
  const __m128 den0 = _mm_set1_ps(1.056512237e+00);
  const __m128 den1 = _mm_set1_ps(-4.474752426e+00);

  const __m128 half = _mm_set1_ps(0.5f);
  const __m128 quater = _mm_set1_ps(0.25f);
  const __m128 signmask = _mm_set1_ps(-0.0f);

  // We use this useful identity to compute on finite range of values:
  // tan(x - pi/2)  == 1 / -tan(x)
  // Since we compute tan(x*pi), identity is tanpi(x-0.5) == 1 / -tanpi(x)

  // Force the value into the finite range (-0.25,0.25)

  // Compute the mask, whether we compute using the identity (x>0.25) or
  // directly (x<0.25)
  __m128 inverse_mask = _mm_cmpgt_ps(x, quater);
  x = _mm_sub_ps(x, _mm_and_ps(inverse_mask, half));

  __m128 f2 = _mm_mul_ps(x, x);
  // Numerator, Horner's scheme
  __m128 num = _mm_add_ps(_mm_mul_ps(num1, f2), num0);
  // Andnot with sign mask to produce only positive results
  num = _mm_mul_ps(_mm_andnot_ps(signmask, x), num);

  // Denominator, Horner's scheme
  __m128 denom = _mm_add_ps(f2, den1);
  denom = _mm_add_ps(den0, _mm_mul_ps(denom, f2));

  // Since we already compute a rational function, we just need to swap
  // numerator and denominator to produce an inverse

  // If (mask) then swap(a,b)
  // Conditional swap with one simple trick...
  __m128 swap_xor = _mm_and_ps(inverse_mask, _mm_xor_ps(num, denom));
  num = _mm_xor_ps(num, swap_xor);
  denom = _mm_xor_ps(denom, swap_xor);

  return _mm_div_ps(num, denom);
}

// Approximates tan(x*pi), x on (1e-35, 0.4999999). Useful for BLT prewrap.
// AVX2+FMA instructions, 8-way parallel
inline __m256 tanpi_blt_fma(__m256 x) {
  const __m256 num0 = _mm256_set1_ps(3.319131136e+00);
  const __m256 num1 = _mm256_set1_ps(-3.138323545e+00);
  const __m256 den0 = _mm256_set1_ps(1.056512237e+00);
  const __m256 den1 = _mm256_set1_ps(-4.474752426e+00);
  const __m256 half = _mm256_set1_ps(0.5f);
  const __m256 quater = _mm256_set1_ps(0.25f);
  const __m256 signmask = _mm256_set1_ps(-0.0f);

  __m256 inverse_mask = _mm256_cmp_ps(x, quater, _CMP_GT_OS);
  x = _mm256_sub_ps(x, _mm256_and_ps(inverse_mask, half));

  __m256 f2 = _mm256_mul_ps(x, x);

  __m256 num = _mm256_fmadd_ps(f2, num1, num0);
  num = _mm256_mul_ps(_mm256_andnot_ps(signmask, x), num);

  __m256 denom = _mm256_add_ps(f2, den1);
  denom = _mm256_fmadd_ps(denom, f2, den0);

  __m256 denom_select = _mm256_blendv_ps(denom, num, inverse_mask);
  __m256 num_select = _mm256_blendv_ps(num, denom, inverse_mask);

  return _mm256_div_ps(num_select, denom_select);
}