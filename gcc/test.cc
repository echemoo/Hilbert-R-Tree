
#include "hrtree_headers.h"
/* See LICENSE below for information on rights to use, modify and distribute
   this code. */

/*
 * hilbert.c - Computes Hilbert space-filling curve coordinates, without
 * recursion, from integer index, and vice versa, and other Hilbert-related
 * calculations.  Also known as Pi-order or Peano scan.
 *
 * Author:      Doug Moore
 *              Dept. of Computational and Applied Math
 *              Rice University
 *              http://www.caam.rice.edu/~dougm
 * Date:        Sun Feb 20 2000
 * Copyright (c) 1998-2000, Rice University
 *
 * Acknowledgement:
 * This implementation is based on the work of A. R. Butz ("Alternative
 * Algorithm for Hilbert's Space-Filling Curve", IEEE Trans. Comp., April,
 * 1971, pp 424-426) and its interpretation by Spencer W. Thomas, University
 * of Michigan (http://www-personal.umich.edu/~spencer/Home.html) in his widely
 * available C software.  While the implementation here differs considerably
 * from his, the first two interfaces and the style of some comments are very
 * much derived from his work. */

/* implementation of the hilbert functions */

#define adjust_rotation(rotation,nDims,bits)                            \
  do {                                                                    \
    /* rotation = (rotation + 1 + ffs(bits)) % nDims; */              \
    bits &= -bits & nd1Ones;                                          \
    while (bits)                                                      \
    bits >>= 1, ++rotation;                                         \
    if ( ++rotation >= nDims )                                        \
    rotation -= nDims;                                              \
  } while (0)

// 生成指定数量的1. 
#define ones(T,k) ((((T)2) << (k-1)) - 1)

#define rdbit(w,k) (((w) >> (k)) & 1)

#define rotateRight(arg, nRots, nDims)                                  \
  ((((arg) >> (nRots)) | ((arg) << ((nDims)-(nRots)))) & ones(bitmask_t,nDims))

#define rotateLeft(arg, nRots, nDims)                                   \
  ((((arg) << (nRots)) | ((arg) >> ((nDims)-(nRots)))) & ones(bitmask_t,nDims))

#define DLOGB_BIT_TRANSPOSE
  static bitmask_t
bitTranspose(unsigned nDims, unsigned nBits, bitmask_t inCoords)
#if defined(DLOGB_BIT_TRANSPOSE)
{
  unsigned const nDims1 = nDims-1;
  unsigned inB = nBits;
  unsigned utB;
  bitmask_t inFieldEnds = 1;
  bitmask_t inMask = ones(bitmask_t,inB);
  bitmask_t coords = 0;

  while ((utB = inB / 2))
  {
    unsigned const shiftAmt = nDims1 * utB;
    bitmask_t const utFieldEnds =
      inFieldEnds | (inFieldEnds << (shiftAmt+utB));
    bitmask_t const utMask =
      (utFieldEnds << utB) - utFieldEnds;
    bitmask_t utCoords = 0;
    unsigned d;
    if (inB & 1)
    {
      bitmask_t const inFieldStarts = inFieldEnds << (inB-1);
      unsigned oddShift = 2*shiftAmt;
      for (d = 0; d < nDims; ++d)
      {
        bitmask_t in = inCoords & inMask;
        inCoords >>= inB;
        coords |= (in & inFieldStarts) << oddShift++;
        in &= ~inFieldStarts;
        in = (in | (in << shiftAmt)) & utMask;
        utCoords |= in << (d*utB);
      }
    }
    else
    {
      for (d = 0; d < nDims; ++d)
      {
        bitmask_t in = inCoords & inMask;
        inCoords >>= inB;
        in = (in | (in << shiftAmt)) & utMask;
        utCoords |= in << (d*utB);
      }
    }
    inCoords = utCoords;
    inB = utB;
    inFieldEnds = utFieldEnds;
    inMask = utMask;
  }
  coords |= inCoords;
  return coords;
}
#else
{
  bitmask_t coords = 0;
  unsigned d;
  for (d = 0; d < nDims; ++d)
  {
    unsigned b;
    bitmask_t in = inCoords & ones(bitmask_t,nBits);
    bitmask_t out = 0;
    inCoords >>= nBits;
    for (b = nBits; b--;)
    {
      out <<= nDims;
      out |= rdbit(in, b);
    }
    coords |= out << d;
  }
  return coords;
}
#endif

/*****************************************************************
 * hilbert_i2c
 *
 * Convert an index into a Hilbert curve to a set of coordinates.
 * Inputs:
 *  nDims:      Number of coordinate axes.
 *  nBits:      Number of bits per axis.
 *  index:      The index, contains nDims*nBits bits
 *              (so nDims*nBits must be <= 8*sizeof(bitmask_t)).
 * Outputs:
 *  coord:      The list of nDims coordinates, each with nBits bits.
 * Assumptions:
 *      nDims*nBits <= (sizeof index) * (bits_per_byte)
 */
  void
hilbert_i2c(unsigned nDims, unsigned nBits, bitmask_t index, bitmask_t coord[])
{
  if (nDims > 1)
  {
    bitmask_t coords;
    halfmask_t const nbOnes = ones(halfmask_t,nBits);
    unsigned d;

    if (nBits > 1)
    {
      unsigned const nDimsBits = nDims*nBits;
      halfmask_t const ndOnes = ones(halfmask_t,nDims);
      halfmask_t const nd1Ones= ndOnes >> 1; /* for adjust_rotation */
      unsigned b = nDimsBits;
      unsigned rotation = 0;
      halfmask_t flipBit = 0;
      bitmask_t const nthbits = ones(bitmask_t,nDimsBits) / ndOnes;
      index ^= (index ^ nthbits) >> 1;
      coords = 0;
      do
      {
        halfmask_t bits = (index >> (b-=nDims)) & ndOnes;
        coords <<= nDims;
        coords |= rotateLeft(bits, rotation, nDims) ^ flipBit;
        flipBit = (halfmask_t)1 << rotation;
        adjust_rotation(rotation,nDims,bits);
      } while (b);
      for (b = nDims; b < nDimsBits; b *= 2)
        coords ^= coords >> b;
      coords = bitTranspose(nBits, nDims, coords);
    }
    else
      coords = index ^ (index >> 1);

    for (d = 0; d < nDims; ++d)
    {
      coord[d] = coords & nbOnes;
      coords >>= nBits;
    }
  }
  else
    coord[0] = index;
}

/*****************************************************************
 * hilbert_c2i
 *
 * Convert coordinates of a point on a Hilbert curve to its index.
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBits:      Number of bits/coordinate.
 *  coord:      Array of n nBits-bit coordinates.
 * Outputs:
 *  index:      Output index value.  nDims*nBits bits.
 * Assumptions:
 *      nDims*nBits <= (sizeof bitmask_t) * (bits_per_byte)
 */
  bitmask_t //unsigned long long
hilbert_c2i(unsigned nDims, unsigned nBits, bitmask_t const coord[])//convert coordinate to hilbert index
{
  if (nDims > 1)
  {
    unsigned const nDimsBits = nDims*nBits;
    bitmask_t index;
    unsigned d;
    bitmask_t coords = 0;
    for (d = nDims; d--; )//copy coord[] into coords. format:"z,y,x"
    {
      coords <<= nBits;//high dimension bits move to left
      coords |= coord[d];
    }

    if (nBits > 1)
    {
      halfmask_t const ndOnes = ones(halfmask_t,nDims);
      halfmask_t const nd1Ones= ndOnes >> 1; /* for adjust_rotation */
      unsigned b = nDimsBits;
      unsigned rotation = 0;
      halfmask_t flipBit = 0;
      bitmask_t const nthbits = ones(bitmask_t,nDimsBits) / ndOnes;
      coords = bitTranspose(nDims, nBits, coords);
      coords ^= coords >> nDims;
      index = 0;
      do
      {
        halfmask_t bits = (coords >> (b-=nDims)) & ndOnes;
        bits = rotateRight(flipBit ^ bits, rotation, nDims);
        index <<= nDims;
        index |= bits;
        flipBit = (halfmask_t)1 << rotation;
        adjust_rotation(rotation,nDims,bits);
      } while (b);
      index ^= nthbits >> 1;
    }
    else
      index = coords;
    for (d = 1; d < nDimsBits; d *= 2)
      index ^= index >> d;
    return index;
  }
  else
    return coord[0];
}

/*****************************************************************
 * Readers and writers of bits
 */

typedef bitmask_t (*BitReader) (unsigned nDims, unsigned nBytes,
    char const* c, unsigned y);
typedef void (*BitWriter) (unsigned d, unsigned nBytes,
    char* c, unsigned y, int fold);


#if defined(sparc)
#define __BIG_ENDIAN__
#endif

#if defined(__BIG_ENDIAN__)
#define whichByte(nBytes,y) (nBytes-1-y/8)
#define setBytes(dst,pos,nBytes,val) \
  memset(&dst[pos+1],val,nBytes-pos-1)
#else
#define whichByte(nBytes,y) (y/8)
#define setBytes(dst,pos,nBytes,val) \
  memset(&dst[0],val,pos)
#endif

  static bitmask_t
getIntBits(unsigned nDims, unsigned nBytes, char const* c, unsigned y)
{
  unsigned const bit = y%8;
  unsigned const offs = whichByte(nBytes,y);
  unsigned d;
  bitmask_t bits = 0;
  c += offs;
  for (d = 0; d < nDims; ++d)
  {
    bits |= rdbit(*c, bit) << d;
    c += nBytes;
  }
  return bits;
}

#include <string.h>
  static void
propogateIntBits(unsigned d, unsigned nBytes,
    char* c, unsigned y, int fold)
{
  unsigned const byteId = whichByte(nBytes,y);
  unsigned const b = y%8;
  char const bthbit = 1 << b;
  char* const target = &c[d*nBytes];
  target[byteId] ^= bthbit;
  if (!fold)
  {
    char notbit = ((target[byteId] >> b) & 1) - 1;
    if (notbit)
      target[byteId] |= bthbit-1;
    else
      target[byteId] &=  -bthbit;
    setBytes(target,byteId,nBytes,notbit);
  }
}

/* An IEEE double is treated as a 2100 bit number.  In particular, 0 is treated
   as a 1 followed by 2099 zeroes, and negative 0 as a 0 followed by 2099 ones.
   Only 53 bits differ between a number and a zero of the same sign, with the
   position of the 53 determined by the exponent, and the values of the 53 by
   the significand (with implicit leading 1 bit).  Although IEEE 754 uses the
   maximum exponent for NaN's and infinities, this implementation ignores that
   decision, so that infinities and NaN's are treated as very large numbers.
   Note that we do not explicitly construct a 2100 bit bitmask in the IEEE
   routines below. */

enum { IEEEexpBits = 11 };
enum { IEEEsigBits = 52 };
enum { IEEErepBits = (1 << IEEEexpBits) + IEEEsigBits };

typedef union ieee754_double
{
  double d;

  /* This is the IEEE 754 double-precision format.  */
  struct
  {
#if defined(__BIG_ENDIAN__)
    unsigned int negative:1;
    unsigned int exponent:11;
    /* Together these comprise the mantissa.  */
    unsigned int mantissa0:20;
    unsigned int mantissa1:32;
#else       /* Big endian.  */
    /* Together these comprise the mantissa.  */
    unsigned int mantissa1:32;
    unsigned int mantissa0:20;
    unsigned int exponent:11;
    unsigned int negative:1;
#endif        /* Little endian.  */
  } ieee;
} ieee754_double;

  static bitmask_t
getIEEESignBits(unsigned nDims, double const* c)
{
  unsigned d;
  ieee754_double x;
  bitmask_t bits = 0;
  for (d = 0; d < nDims; ++d)
  {
    x.d = c[d];
    bits |= x.ieee.negative << d;
  }
  return bits;
}

  static bitmask_t
getIEEEBits(unsigned nDims,
    unsigned ignoreMe, /* ignored */
    char const* cP,
    unsigned y)
  /* retrieve bits y of elements of double array c, where an expanded IEEE
     double has 2100 bits. */
{
  unsigned d;
  double const* c = (double const*) cP;
  ieee754_double x;
  bitmask_t bits = 0;
  for (x.d = c[d=0]; d < nDims; x.d = c[++d])
  {
    bitmask_t bit = x.ieee.negative;
    unsigned normalized = (x.ieee.exponent != 0);
    unsigned diff = y - (x.ieee.exponent - normalized);
    if (diff <= 52)
      bit ^= 1 & ((diff <  32)? x.ieee.mantissa1 >> diff:
          (diff <  52)? x.ieee.mantissa0 >> (diff - 32):
          /* else */    normalized);
    else
      bit ^= (y == IEEErepBits-1);

    bits |= bit << d;
  }
  return bits;
}

  static void
propogateIEEEBits(unsigned d, unsigned nBytes,
    char* cP, unsigned y, int fold)
{
  ieee754_double* x = d + (ieee754_double*) cP;
  unsigned normalized = (x->ieee.exponent != 0);
  unsigned diff = y - (x->ieee.exponent - normalized);
  if (diff < 32)
  {
    unsigned b = 1 << diff;
    unsigned bit = x->ieee.mantissa1 & b;
    x->ieee.mantissa1 &= ~(b-1);
    x->ieee.mantissa1 |= b;
    if (bit)
      --x->ieee.mantissa1;
  }
  else if (diff < 52)
  {
    unsigned b = 1 << (diff - 32);
    unsigned bit = x->ieee.mantissa0 & b;
    x->ieee.mantissa0 &= ~(b-1);
    x->ieee.mantissa0 |= b;
    if (bit)
      --x->ieee.mantissa0;
    x->ieee.mantissa1 = bit?-1: 0;
  }
  else if (diff == 52) /* "flip" the implicit 1 bit */
  {
    if (normalized)
      --x->ieee.exponent;
    else
      x->ieee.exponent = 1;
    x->ieee.mantissa0 = -normalized;
    x->ieee.mantissa1 = -normalized;
  }
  else if (diff < IEEErepBits)
  {
    if (y == IEEErepBits-1)
    {
      x->ieee.negative ^= 1;
      x->ieee.exponent = 0;
    }
    else
      x->ieee.exponent = y - 51;
    x->ieee.mantissa0 = 0;
    x->ieee.mantissa1 = 0;
  }
}

  static unsigned
getIEEEexptMax(unsigned nDims, double const* c)
{
  unsigned max = 0;
  unsigned d;
  for (d = 0; d < nDims; ++d)
  {
    ieee754_double x;
    x.d = c[d];
    if (max < x.ieee.exponent)
      max = x.ieee.exponent;
  }
  if (max) --max;
  return max;
}

  static void
getIEEEinitValues(double const* c1,
    unsigned y,
    unsigned nDims,
    unsigned* rotation,
    bitmask_t* bits,
    bitmask_t* index)
{
  bitmask_t const one = 1;
  unsigned d;
  bitmask_t signBits = getIEEESignBits(nDims, c1);
  unsigned signParity, leastZeroBit, strayBit;

  /* compute the odd/evenness of the number of sign bits */
  {
    bitmask_t signPar = signBits;
    for (d = 1; d < nDims; d *= 2)
      signPar ^= signPar >> d;
    signParity = signPar & 1;
  }

  /* find the position of the least-order 0 bit in among signBits and adjust it
     if necessary */
  for (leastZeroBit = 0; leastZeroBit < nDims; ++leastZeroBit)
    if (rdbit(signBits, leastZeroBit) == 0)
      break;
  strayBit = 0;
  if (leastZeroBit == nDims-2)
    strayBit = 1;
  else if (leastZeroBit == nDims)
    leastZeroBit = nDims-1;

  if (y % 2 == 1)
  {
    *rotation = (IEEErepBits - y + 1 + leastZeroBit) % nDims;
    if (y < IEEErepBits-1)
    {
      *bits = signBits ^ (one << ((*rotation + strayBit) % nDims));
      *index = signParity;
    }
    else /* y == IEEErepBits-1 */
    {
      *bits = signBits ^ (ones(bitmask_t,nDims) &~ 1);
      *index =  signParity ^ (nDims&1);
    }
  }
  else /* y % 2 == 0 */
    if (y < IEEErepBits)
    {
      unsigned shift_amt = (IEEErepBits - y + leastZeroBit) % nDims;
      *rotation = (shift_amt + 2 + strayBit) % nDims;
      *bits = signBits ^ (one << shift_amt);
      *index = signParity ^ 1;
    }
    else /* y == IEEErepBits */
    {
      *rotation = 0;
      *bits = one << (nDims-1);
      *index = 1;
    }
}

/*****************************************************************
 * hilbert_cmp, hilbert_ieee_cmp
 *
 * Determine which of two points lies further along the Hilbert curve
 * Inputs:
 *  nDims:      Number of coordinates.
 *  nBytes:     Number of bytes of storage/coordinate (hilbert_cmp only)
 *  nBits:      Number of bits/coordinate. (hilbert_cmp only)
 *  coord1:     Array of nDims nBytes-byte coordinates (or doubles for ieee_cmp).
 *  coord2:     Array of nDims nBytes-byte coordinates (or doubles for ieee_cmp).
 * Return value:
 *      -1, 0, or 1 according to whether
 coord1<coord2, coord1==coord2, coord1>coord2
 * Assumptions:
 *      nBits <= (sizeof bitmask_t) * (bits_per_byte)
 */

