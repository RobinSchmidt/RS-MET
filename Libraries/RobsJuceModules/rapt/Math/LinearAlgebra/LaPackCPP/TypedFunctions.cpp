#pragma once // this should avoid multiple definition errors...hopefully


#include <stdarg.h>
#include <stdio.h>
#include <algorithm> // min/max

#include "LaPack.hpp"

/** For some of the functions in LibF2C, Blas and LaPack, it didn't make sense to templatize them.
These are all put into this file because it may be problematic in certain contexts to have mixed 
templated and typed code in a single implementation file - if you include the same implementation
file in two different unity-build compilation units, you will get "multiple definition" linker 
errors. You may have to include the same implementation file into different compilation unit, if
these compilation units instantiate different templates or the same template for different 
datatypes. I know - it's a mess :-( */

namespace LaPackCPP {

template<class T> inline T min(T x, T y) {  return std::min(x, y); }
template<class T> inline T max(T x, T y) {  return std::max(x, y); }

//=================================================================================================
// LibF2C

integer i_len(char *s, ftnlen n)
{
  return(n);
}

integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
  unsigned char *a, *aend, *b, *bend;
  a = (unsigned char *)a0;
  b = (unsigned char *)b0;
  aend = a + la;
  bend = b + lb;

  if(la <= lb)
  {
    while(a < aend)
      if(*a != *b)
        return( *a - *b );
      else
      { ++a; ++b; }

    while(b < bend)
      if(*b != ' ')
        return( ' ' - *b );
      else	++b;
  }

  else
  {
    while(b < bend)
      if(*a == *b)
      { ++a; ++b; }
      else
        return( *a - *b );
    while(a < aend)
      if(*a != ' ')
        return(*a - ' ');
      else	++a;
  }
  return(0);
}

void s_copy(char *a, char *b, ftnlen la, ftnlen lb)
{
  char *aend, *bend;

  aend = a + la;

  if(la <= lb)
#ifndef NO_OVERWRITE
    if(a <= b || a >= b + la)
#endif
      while(a < aend)
        *a++ = *b++;
#ifndef NO_OVERWRITE
    else
      for(b += la; a < aend; )
        *--aend = *--b;
#endif

  else {
    bend = b + lb;
#ifndef NO_OVERWRITE
    if(a <= b || a >= bend)
#endif
      while(b < bend)
        *a++ = *b++;
#ifndef NO_OVERWRITE
    else {
      a += lb;
      while(b < bend)
        *--a = *--bend;
      a += lb;
    }
#endif
    while(a < aend)
      *a++ = ' ';
  }
}


//=================================================================================================
// Blas

// translated from lsame, Reference BLAS level1 routine (version 3.1)
logical lsame(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
  // System generated locals
  logical ret_val;

  // Local variables 
  static integer inta, intb, zcode;

  // Test if the characters are equal
  ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
  if(ret_val) {
    return ret_val;
  }

  // Now test for equivalence if both characters are alphabetic.
  zcode = 'Z';

  // Use 'Z' rather than 'A' so that ASCII can be detected on Prime 
  // machines, on which ICHAR returns a value with bit 8 set. 
  // ICHAR('A') on Prime machines returns 193 which is the same as 
  // ICHAR('A') on an EBCDIC machine.

  inta = *(unsigned char *)ca;
  intb = *(unsigned char *)cb;

  if(zcode == 90 || zcode == 122) {

    // ASCII is assumed - ZCODE is the ASCII code of either lower or upper case 'Z'. 
    if(inta >= 97 && inta <= 122) {
      inta += -32;
    }
    if(intb >= 97 && intb <= 122) {
      intb += -32;
    }

  }
  else if(zcode == 233 || zcode == 169) {

    // EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or 
    // upper case 'Z'.
    if(inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta
      >= 162 && inta <= 169) {
      inta += 64;
    }
    if(intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb
      >= 162 && intb <= 169) {
      intb += 64;
    }

  }
  else if(zcode == 218 || zcode == 250) {

    // ASCII is assumed, on Prime machines - ZCODE is the ASCII code
    // plus 128 of either lower or upper case 'Z'. 

    if(inta >= 225 && inta <= 250) {
      inta += -32;
    }
    if(intb >= 225 && intb <= 250) {
      intb += -32;
    }
  }
  ret_val = inta == intb;

  // End of LSAME

  return ret_val;
}


int xerbla(char *srname, integer *info, ftnlen srname_len)
{
  // Code of the function has been commented out by Robin Schmidt - at the moment, xerbla is just
  // a dummy function - i should probably set a debug-breakpoint here and re-implement it 
  // completely - it deals with some weird I/O functions from the f2c library which seems to be 
  // useless clutter...

  /*
  // Table of constant values
  static integer c__1 = 1;

  // Format strings 
  static char fmt_9999[] = "(\002 ** On entry to \002,a,\002 parameter num"
    "ber \002,i2,\002 had \002,\002an illegal value\002)";

  // Builtin functions 
  integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
  int s_stop(char *, ftnlen);

  // Local variables
  extern integer len_trim__(char *, ftnlen);
  // Note by Robin Schmidt:
  // this was declared as an intrinsic function in the original xerbla.f file, but it stumped the 
  // f2c translator, giving an error about an unknown intrinsic function, so i changed the 
  // "INTRINSIC" keyword to "EXTERNAL" - that allowed the translation, but i'm not sure, if it 
  // gives the correct behavior 

  // Fortran I/O blocks
  static cilist io___1 = { 0, 6, 0, fmt_9999, 0 };

  s_wsfe(&io___1);
  do_fio(&c__1, srname, len_trim__(srname, srname_len));
  do_fio(&c__1, (char *)&(*info), (ftnlen)sizeof(integer));
  e_wsfe();

  s_stop("", (ftnlen)0);

  // End of XERBLA
  */

  return 0;
} // xerbla

int xerbla(const char* srname, integer* info, ftnlen srname_len)
{
  return 0;
}


//=================================================================================================
// XBlas

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

//=================================================================================================
// LaPack

// LAPACK computational routine (version 3.7.0)
void chla_transtype(char *ret_val, ftnlen ret_val_len, integer *trans)
{
  if (*trans == 111) {
    *(unsigned char *)ret_val = 'N';
  } else if (*trans == 112) {
    *(unsigned char *)ret_val = 'T';
  } else if (*trans == 113) {
    *(unsigned char *)ret_val = 'C';
  } else {
    *(unsigned char *)ret_val = 'X';
  }
  return ;
} 

// LAPACK auxiliary routine (version 3.7.0) 
integer ieeeck(integer *ispec, f2c_real *zero, f2c_real *one)
{
  // System generated locals
  integer ret_val;

  // Local variables
  static f2c_real nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, negzro, 
    newzro;

  ret_val = 1;

  posinf = *one / *zero;
  if (posinf <= *one) {
    ret_val = 0;
    return ret_val;
  }
  neginf = -(*one) / *zero;
  if (neginf >= *zero) {
    ret_val = 0;
    return ret_val;
  }
  negzro = *one / (neginf + *one);
  if (negzro != *zero) {
    ret_val = 0;
    return ret_val;
  }
  neginf = *one / negzro;
  if (neginf >= *zero) {
    ret_val = 0;
    return ret_val;
  }
  newzro = negzro + *zero;
  if (newzro != *zero) {
    ret_val = 0;
    return ret_val;
  }
  posinf = *one / newzro;
  if (posinf <= *one) {
    ret_val = 0;
    return ret_val;
  }
  neginf *= posinf;
  if (neginf >= *zero) {
    ret_val = 0;
    return ret_val;
  }
  posinf *= posinf;
  if (posinf <= *one) {
    ret_val = 0;
    return ret_val;
  }

  // Return if we were only asked to check infinity arithmetic
  if (*ispec == 0) {
    return ret_val;
  }

  nan1 = posinf + neginf;
  nan2 = posinf / neginf;
  nan3 = posinf / posinf;
  nan4 = posinf * *zero;
  nan5 = neginf * negzro;
  nan6 = nan5 * *zero;

  if (nan1 == nan1) {
    ret_val = 0;
    return ret_val;
  }
  if (nan2 == nan2) {
    ret_val = 0;
    return ret_val;
  }
  if (nan3 == nan3) {
    ret_val = 0;
    return ret_val;
  }
  if (nan4 == nan4) {
    ret_val = 0;
    return ret_val;
  }
  if (nan5 == nan5) {
    ret_val = 0;
    return ret_val;
  }
  if (nan6 == nan6) {
    ret_val = 0;
    return ret_val;
  }

  return ret_val;
} // ieeeck

//-------------------------------------------------------------------------------------------------

// LAPACK auxiliary routine (version 3.8.0) 
integer ilaenv(integer *ispec, char *name__, char *opts, integer *n1, 
  integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len)
{
  // Table of constant values
  static integer c__1 = 1;
  static f2c_real c_b173 = 0.f;
  static f2c_real c_b174 = 1.f;
  static integer c__0 = 0;

  // System generated locals
  integer ret_val;

  // Local variables 
  static logical twostage;
  static integer i__;
  static char c1[1], c2[2], c3[3], c4[2];
  static integer ic, nb, iz, nx;
  static logical cname;
  static integer nbmin;
  static logical sname;
  extern integer ieeeck(integer *, f2c_real *, f2c_real *);
  static char subnam[16];
  extern integer iparmq(integer *, char *, char *, integer *, integer *, 
    integer *, integer *, ftnlen, ftnlen);

  switch (*ispec) {
  case 1:  goto L10;
  case 2:  goto L10;
  case 3:  goto L10;
  case 4:  goto L80;
  case 5:  goto L90;
  case 6:  goto L100;
  case 7:  goto L110;
  case 8:  goto L120;
  case 9:  goto L130;
  case 10:  goto L140;
  case 11:  goto L150;
  case 12:  goto L160;
  case 13:  goto L160;
  case 14:  goto L160;
  case 15:  goto L160;
  case 16:  goto L160;
  }

  // Invalid value for ISPEC

  ret_val = -1;
  return ret_val;

L10:

  // Convert NAME to upper case if the first character is lower case
  ret_val = 1;
  s_copy(subnam, name__, (ftnlen)16, name_len);
  ic = *(unsigned char *)subnam;
  iz = 'Z';
  if (iz == 90 || iz == 122) {

    // ASCII character set
    if (ic >= 97 && ic <= 122) {
      *(unsigned char *)subnam = (char) (ic - 32);
      for (i__ = 2; i__ <= 6; ++i__) {
        ic = *(unsigned char *)&subnam[i__ - 1];
        if (ic >= 97 && ic <= 122) {
          *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
        }
        // L20: 
      }
    }

  } else if (iz == 233 || iz == 169) {

    // EBCDIC character set 

    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
      ic <= 169) {
      *(unsigned char *)subnam = (char) (ic + 64);
      for (i__ = 2; i__ <= 6; ++i__) {
        ic = *(unsigned char *)&subnam[i__ - 1];
        if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
          162 && ic <= 169) {
          *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
        }
        // L30:
      }
    }

  } else if (iz == 218 || iz == 250) {

    // Prime machines:  ASCII+128

    if (ic >= 225 && ic <= 250) {
      *(unsigned char *)subnam = (char) (ic - 32);
      for (i__ = 2; i__ <= 6; ++i__) {
        ic = *(unsigned char *)&subnam[i__ - 1];
        if (ic >= 225 && ic <= 250) {
          *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
        }
        // L40:
      }
    }
  }

  *(unsigned char *)c1 = *(unsigned char *)subnam;
  sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
  cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
  if (! (cname || sname)) {
    return ret_val;
  }
  s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
  s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
  s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);
  twostage = i_len(subnam, (ftnlen)16) >= 11 && *(unsigned char *)&subnam[
    10] == '2';

  switch (*ispec) {
  case 1:  goto L50;
  case 2:  goto L60;
  case 3:  goto L70;
  }

L50:

  // ISPEC = 1:  block size

  // In these examples, separate code is provided for setting NB for 
  // real and complex.  We assume that NB will take the same value in 
  // single or double precision. 

  nb = 1;

  if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    } else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
      "RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
        3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
      == 0) {
      if (sname) {
        nb = 32;
      } else {
        nb = 32;
      }
    } else if (s_cmp(c3, "QR ", (ftnlen)3, (ftnlen)3) == 0) {
      if (*n3 == 1) {
        if (sname) {
          //  M*N 
          if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
            nb = *n1;
          } else {
            nb = 32768 / *n2;
          }
        } else {
          if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
            nb = *n1;
          } else {
            nb = 32768 / *n2;
          }
        }
      } else {
        if (sname) {
          nb = 1;
        } else {
          nb = 1;
        }
      }
    } else if (s_cmp(c3, "LQ ", (ftnlen)3, (ftnlen)3) == 0) {
      if (*n3 == 2) {
        if (sname) {
          // M*N 
          if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
            nb = *n1;
          } else {
            nb = 32768 / *n2;
          }
        } else {
          if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
            nb = *n1;
          } else {
            nb = 32768 / *n2;
          }
        }
      } else {
        if (sname) {
          nb = 1;
        } else {
          nb = 1;
        }
      }
    } else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 32;
      } else {
        nb = 32;
      }
    } else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 32;
      } else {
        nb = 32;
      }
    } else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    }
  } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    }
  } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        if (twostage) {
          nb = 192;
        } else {
          nb = 64;
        }
      } else {
        if (twostage) {
          nb = 192;
        } else {
          nb = 64;
        }
      }
    } else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 32;
    } else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 64;
    }
  } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (twostage) {
        nb = 192;
      } else {
        nb = 64;
      }
    } else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 32;
    } else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 64;
    }
  } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
      }
    } else if (*(unsigned char *)c3 == 'M') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
      }
    }
  } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
      }
    } else if (*(unsigned char *)c3 == 'M') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
      }
    }
  } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        if (*n4 <= 64) {
          nb = 1;
        } else {
          nb = 32;
        }
      } else {
        if (*n4 <= 64) {
          nb = 1;
        } else {
          nb = 32;
        }
      }
    }
  } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        if (*n2 <= 64) {
          nb = 1;
        } else {
          nb = 32;
        }
      } else {
        if (*n2 <= 64) {
          nb = 1;
        } else {
          nb = 32;
        }
      }
    }
  } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    } else if (s_cmp(c3, "EVC", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    }
  } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    }
  } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 1;
    }
  } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
    nb = 32;
    if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 32;
      } else {
        nb = 32;
      }
    }
  }
  ret_val = nb;
  return ret_val;

L60:

  // ISPEC = 2:  minimum block size 

  nbmin = 2;
  if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
      ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
        ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
    {
      if (sname) {
        nbmin = 2;
      } else {
        nbmin = 2;
      }
    } else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nbmin = 2;
      } else {
        nbmin = 2;
      }
    } else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nbmin = 2;
      } else {
        nbmin = 2;
      }
    } else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nbmin = 2;
      } else {
        nbmin = 2;
      }
    }
  } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nbmin = 8;
      } else {
        nbmin = 8;
      }
    } else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nbmin = 2;
    }
  } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nbmin = 2;
    }
  } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
      }
    } else if (*(unsigned char *)c3 == 'M') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
      }
    }
  } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
      }
    } else if (*(unsigned char *)c3 == 'M') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
      }
    }
  } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
    nbmin = 2;
    if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
      nbmin = 2;
    }
  }
  ret_val = nbmin;
  return ret_val;

L70:

  // ISPEC = 3:  crossover point 

  nx = 0;
  if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
      ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
        ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
    {
      if (sname) {
        nx = 128;
      } else {
        nx = 128;
      }
    } else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nx = 128;
      } else {
        nx = 128;
      }
    } else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nx = 128;
      } else {
        nx = 128;
      }
    }
  } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
    if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nx = 32;
    }
  } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nx = 32;
    }
  } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nx = 128;
      }
    }
  } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nx = 128;
      }
    }
  } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
    nx = 128;
    if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
      nx = 128;
    }
  }
  ret_val = nx;
  return ret_val;

L80:

  // ISPEC = 4:  number of shifts (used by xHSEQR)

  ret_val = 6;
  return ret_val;

L90:

  // ISPEC = 5:  minimum column dimension (not used)

  ret_val = 2;
  return ret_val;

L100:

  // ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)

  ret_val = (integer) ((f2c_real) min(*n1,*n2) * 1.6f);
  return ret_val;

L110:

  // ISPEC = 7:  number of processors (not used) 

  ret_val = 1;
  return ret_val;

L120:

  // ISPEC = 8:  crossover point for multishift (used by xHSEQR)

  ret_val = 50;
  return ret_val;

L130:

  // ISPEC = 9:  maximum size of the subproblems at the bottom of the 
  //             computation tree in the divide-and-conquer algorithm 
  //             (used by xGELSD and xGESDD) 

  ret_val = 25;
  return ret_val;

L140:

  // ISPEC = 10: ieee NaN arithmetic can be trusted not to trap 

  // ILAENV = 0 
  ret_val = 1;
  if (ret_val == 1) {
    ret_val = ieeeck(&c__1, &c_b173, &c_b174);
  }
  return ret_val;

L150:

  // ISPEC = 11: infinity arithmetic can be trusted not to trap 

  // ILAENV = 0 
  ret_val = 1;
  if (ret_val == 1) {
    ret_val = ieeeck(&c__0, &c_b173, &c_b174);
  }
  return ret_val;

L160:

  // 12 <= ISPEC <= 16: xHSEQR or related subroutines.

  ret_val = iparmq(ispec, name__, opts, n1, n2, n3, n4, name_len, opts_len)
    ;
  return ret_val;

  // End of ILAENV 

} // ilaenv

//-------------------------------------------------------------------------------------------------

// LAPACK computational routine (version 3.7.0) 
integer ilaprec(char *prec, ftnlen prec_len)
{
  // System generated locals 
  integer ret_val;


  if (lsame(prec, "S", (ftnlen)1, (ftnlen)1)) {
    ret_val = 211;
  } else if (lsame(prec, "D", (ftnlen)1, (ftnlen)1)) {
    ret_val = 212;
  } else if (lsame(prec, "I", (ftnlen)1, (ftnlen)1)) {
    ret_val = 213;
  } else if (lsame(prec, "X", (ftnlen)1, (ftnlen)1) || lsame(prec, "E", (
    ftnlen)1, (ftnlen)1)) {
    ret_val = 214;
  } else {
    ret_val = -1;
  }
  return ret_val;

  // End of ILAPREC 

} // ilaprec

//-------------------------------------------------------------------------------------------------

//  LAPACK computational routine (version 3.7.0) 
integer ilatrans(char *trans, ftnlen trans_len)
{
  integer ret_val;
  if (lsame(trans, "N", (ftnlen)1, (ftnlen)1)) {
    ret_val = 111;
  } else if (lsame(trans, "T", (ftnlen)1, (ftnlen)1)) {
    ret_val = 112;
  } else if (lsame(trans, "C", (ftnlen)1, (ftnlen)1)) {
    ret_val = 113;
  } else {
    ret_val = -1;
  }
  return ret_val;
} 

//-------------------------------------------------------------------------------------------------

// LAPACK auxiliary routine (version 3.7.1)
integer iparmq(integer *ispec, char *name__, char *opts, integer *n, integer 
  *ilo, integer *ihi, integer *lwork, ftnlen name_len, ftnlen opts_len)
{
  // System generated locals 
  integer ret_val, i__1, i__2;
  f2c_real r__1;

  // Local variables 
  static integer i__, ic, nh, ns, iz;
  static char subnam[6];

  if (*ispec == 15 || *ispec == 13 || *ispec == 16) {

    // ==== Set the number simultaneous shifts ==== 

    nh = *ihi - *ilo + 1;
    ns = 2;
    if (nh >= 30) {
      ns = 4;
    }
    if (nh >= 60) {
      ns = 10;
    }
    if (nh >= 150) {
      // Computing MAX
      r__1 = (f2c_real) (log((f2c_real) nh) / log(2.f)); // outer cast to f2c_real added by Robin Schmidt
      i__1 = 10, i__2 = nh / i_nint(&r__1);
      ns = max(i__1,i__2);
    }
    if (nh >= 590) {
      ns = 64;
    }
    if (nh >= 3000) {
      ns = 128;
    }
    if (nh >= 6000) {
      ns = 256;
    }
    // Computing MAX 
    i__1 = 2, i__2 = ns - ns % 2;
    ns = max(i__1,i__2);
  }

  if (*ispec == 12) {


    // ===== Matrices of order smaller than NMIN get sent 
    // .     to xLAHQR, the classic double shift algorithm.
    // .     This must be at least 11. ==== 

    ret_val = 75;

  } else if (*ispec == 14) {

    // ==== INIBL: skip a multi-shift qr iteration and 
    // .    whenever aggressive early deflation finds 
    // .    at least (NIBBLE*(window size)/100) deflations. ==== 

    ret_val = 14;

  } else if (*ispec == 15) {

    // ==== NSHFTS: The number of simultaneous shifts =====

    ret_val = ns;

  } else if (*ispec == 13) {

    // ==== NW: deflation window size.  ==== 

    if (nh <= 500) {
      ret_val = ns;
    } else {
      ret_val = ns * 3 / 2;
    }

  } else if (*ispec == 16) {

    // ==== IACC22: Whether to accumulate reflections 
    // .     before updating the far-from-diagonal elements 
    // .     and whether to use 2-by-2 block structure while 
    // .     doing it.  A small amount of work could be saved 
    // .     by making this choice dependent also upon the 
    // .     NH=IHI-ILO+1. 

    // Convert NAME to upper case if the first character is lower case.
    ret_val = 0;
    s_copy(subnam, name__, (ftnlen)6, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

      // ASCII character set 

      if (ic >= 97 && ic <= 122) {
        *(unsigned char *)subnam = (char) (ic - 32);
        for (i__ = 2; i__ <= 6; ++i__) {
          ic = *(unsigned char *)&subnam[i__ - 1];
          if (ic >= 97 && ic <= 122) {
            *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
          }
        }
      }

    } else if (iz == 233 || iz == 169) {

      // EBCDIC character set 
      if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 
        && ic <= 169) {
        *(unsigned char *)subnam = (char) (ic + 64);
        for (i__ = 2; i__ <= 6; ++i__) {
          ic = *(unsigned char *)&subnam[i__ - 1];
          if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || 
            ic >= 162 && ic <= 169) {
            *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
          }
        }
      }

    } else if (iz == 218 || iz == 250) {

      // Prime machines:  ASCII+128
      if (ic >= 225 && ic <= 250) {
        *(unsigned char *)subnam = (char) (ic - 32);
        for (i__ = 2; i__ <= 6; ++i__) {
          ic = *(unsigned char *)&subnam[i__ - 1];
          if (ic >= 225 && ic <= 250) {
            *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
          }
        }
      }
    }

    if (s_cmp(subnam + 1, "GGHRD", (ftnlen)5, (ftnlen)5) == 0 || s_cmp(
      subnam + 1, "GGHD3", (ftnlen)5, (ftnlen)5) == 0) {
      ret_val = 1;
      if (nh >= 14) {
        ret_val = 2;
      }
    } else if (s_cmp(subnam + 3, "EXC", (ftnlen)3, (ftnlen)3) == 0) {
      if (nh >= 14) {
        ret_val = 1;
      }
      if (nh >= 14) {
        ret_val = 2;
      }
    } else if (s_cmp(subnam + 1, "HSEQR", (ftnlen)5, (ftnlen)5) == 0 || 
      s_cmp(subnam + 1, "LAQR", (ftnlen)4, (ftnlen)4) == 0) {
      if (ns >= 14) {
        ret_val = 1;
      }
      if (ns >= 14) {
        ret_val = 2;
      }
    }

  } else {
    // ===== invalid value of ispec ===== 
    ret_val = -1;

  }

  // ==== End of IPARMQ ==== 

  return ret_val;
} // iparmq

//-------------------------------------------------------------------------------------------------

// todo: in some cases, i guessed what function should map to a particular member of 
// numeric_limits - verify everything! ..also in smach.c the functions also return double precision
// values - check that - maybe we should return double here, too instead of T

template<class T>
integer minexponent_(T *dummy)
{
  return std::numeric_limits<T>::min_exponent; // there's also a min_exponent10, but they probably mean that
}
template<class T>
integer maxexponent_(T *dummy)
{
  return std::numeric_limits<T>::max_exponent; // dito
}
template<class T>
T huge_(T *dummy)
{
  return std::numeric_limits<T>::max();  // guessed!!
}
template<class T>
T tiny_(T *dummy)
{
  return std::numeric_limits<T>::min();  // guessed!!
}
template<class T>
T radix_(T *dummy)
{
  return std::numeric_limits<T>::radix;  // this is probably right
}
template<class T>
T digits_(T *dummy)
{
  return std::numeric_limits<T>::digits;  // this probably too
}
template<class T>
T epsilon_(T *dummy)
{
  return std::numeric_limits<T>::epsilon(); // this is very probably right
}

// DLAMC3 is intended to force A and B to be stored prior to doing the addition of A and B, for use
// in situations where optimizers might hold one of these in a register.
// from dlamch.f - LAPACK auxiliary routine (version 3.7.0)
// templatize!
doublereal lamc3(doublereal *a, doublereal *b)
{
  doublereal ret_val;
  ret_val = *a + *b;
  return ret_val;
}

// from dlamch - LAPACK auxiliary routine (version 3.7.0) --
doublereal lamch(char *cmach, ftnlen cmach_len)
{
  static double c_b2 = 0.;  // hmm..it seems, this is the dummy in the original code

  // System generated locals
  doublereal ret_val;

  // Local variables 
  static double rnd, eps;
  static double rmach;
  static double small, sfmin;

  // Assume rounding, not chopping. Always.
  rnd = 1.;

  if (1. == rnd) {
    eps = epsilon_(&c_b2) * .5f;
  } else {
    eps = epsilon_(&c_b2);
  }

  if (lsame(cmach, "E", (ftnlen)1, (ftnlen)1)) {
    rmach = eps;
  } else if (lsame(cmach, "S", (ftnlen)1, (ftnlen)1)) {
    sfmin = tiny_(&c_b2);
    small = 1. / huge_(&c_b2);
    if (small >= sfmin) {

      // Use SMALL plus a bit, to avoid the possibility of rounding 
      // causing overflow when computing  1/sfmin. 

      sfmin = small * (eps + 1.);
    }
    rmach = sfmin;
  } else if (lsame(cmach, "B", (ftnlen)1, (ftnlen)1)) {
    rmach = radix_(&c_b2);
  } else if (lsame(cmach, "P", (ftnlen)1, (ftnlen)1)) {
    rmach = eps * radix_(&c_b2);
  } else if (lsame(cmach, "N", (ftnlen)1, (ftnlen)1)) {
    rmach = digits_(&c_b2);
  } else if (lsame(cmach, "R", (ftnlen)1, (ftnlen)1)) {
    rmach = rnd;
  } else if (lsame(cmach, "M", (ftnlen)1, (ftnlen)1)) {
    rmach = (doublereal) minexponent_(&c_b2);
  } else if (lsame(cmach, "U", (ftnlen)1, (ftnlen)1)) {
    rmach = tiny_(&c_b2);
  } else if (lsame(cmach, "L", (ftnlen)1, (ftnlen)1)) {
    rmach = (doublereal) maxexponent_(&c_b2);
  } else if (lsame(cmach, "O", (ftnlen)1, (ftnlen)1)) {
    rmach = huge_(&c_b2);
  } else {
    rmach = 0.;
  }

  ret_val = rmach;
  return ret_val;

} // dlamch



}