/* [ file: ffnrm.c ] */
# define DO_FFNRM "ffnrm(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C Program fftrn(*)                                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Normalization routine for fast fourier fransformation frogram FFTRN.C      *
*                                                                              *
*   This function normalizes the distributions and spectra in such a way       *
*   that the integral over the absolute values equals C <> ZERO.               *
*   The value of C is transferred from the calling program by                  *
*   writing it into double fpt->nor.                                           *
*   For  C=nor=ZERO , the integrals over the real and imaginary                *
*   parts are normalized to ZERO by substracting one complex                   *
*   constant from all sample values.                                           *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define _POSIX_SOURCE 0 /* set to 1: if POSIX.1 standard headers will be used*/
/*----------------------------------------------------------------------------*/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <time.h>
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma  OPTIMIZE ON
   # pragma  OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# include "./maths.h"  /* 'my' computation environment headers */
# include "./consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header */
# include "../CONFIG.H"
/*----------------------------------------------------------------------------*/
# include "./math/FFT.CONF"
# include "./math/fftyp.h"
/*----------------------------------------------------------------------------*/
# ifndef READ_REVERSE_
   # define READ_REVERSE_ 1 /* 1: reversed spectral domain convention:        */
# endif                     /*    [ 0. ,..., pos. bound , neg. bound,...,0. ] */
/*============================================================================*/

FFT *
ffnrm( FFT *fpt )
{
/* allusions: */

/* declarations: */

   static double 
      norm = ZERO,
      rnrm = ZERO,
      inrm = ZERO,
      xx = ZERO,
      yy = ZERO,
      dx = ZERO;

   static long 
      ii = null,
      jj = null,
      nn = null;

   static short 
      hh = null;

/* function prototypes: */

   double sqrt( double r );
/*----------------------------------------------------------------------------*/
   printf( "\n normalization subroutine started " );
/*----------------------------------------------------------------------------*/
   for ( hh=(fpt->p); hh<=(fpt->q) ; hh++ )
   {
      if (( *(fpt->opt) != 'b' )&&( *(fpt->opt) != 'B' )) 
      {                            /* distribution, pulse ... [ time domain ] */
         nn = (fpt->ttlg[hh]); 
         dx = (fpt->dt[hh]);
      }
      else if (( *(fpt->opt) == 'b')||( *(fpt->opt) == 'B' ))
      {         /* characteristic function, spectrum ... [ frequency domain ] */

         nn = (fpt->stlg[hh]);
         dx = (fpt->ds[hh]);
      };
/*----------------------------------------------------------------------------*/
# if READ_REVERSE_ == 1
      jj = ( long ) ( nn / TWO );
# else
      jj = null;
# endif
/*----------------------------------------------------------------------------*/
      xx = (fpt->r[hh][jj]);
      yy = (fpt->i[hh][jj]);

      if ( (fpt->nor) == null )
      {
         rnrm = xx;
         inrm = yy;
      }
      else
      {
         norm = sqrt( xx*xx + yy*yy )/2.;
      }; 

      ii=ONE; 
      do
      {
         xx = (fpt->r[hh][jj]);
         yy = (fpt->i[hh][jj]);

         if ( (fpt->nor) == ZERO )
         {
            rnrm += xx;
            inrm += yy;
         }
         else
         {
            norm += sqrt( xx*xx + yy*yy );
         };

         jj++ ;
         if ( nn <= jj )
         {
             jj -= nn;
         };
         ii++ ;
      }  while ( ii < nn );

      xx = (fpt->r[hh][jj]);
      yy = (fpt->i[hh][jj]);

      if ( fpt->nor == ZERO )
      {
         rnrm += xx;
         inrm += yy;

         rnrm /= nn;
         inrm /= nn; 

         for ( ii=null ; ii<nn ; ii++ )
         {
            (fpt->r[hh][ii]) -= rnrm;
            (fpt->i[hh][ii]) -= inrm;
         };
      }
      else
      {
         norm += sqrt( xx*xx + yy*yy )/2.;
         norm *= ( dx/(fpt->nor) );

         for ( ii=null ; ii<nn ; ii++ )
         {
            (fpt->r[hh][ii]) /= norm;
            (fpt->i[hh][ii]) /= norm;
         };
      };
   };/* next hh */
/*----------------------------------------------------------------------------*/
/* fende: */

   printf( "\n normalization subroutine terminated \n" );
   (fpt->rtn) = null;
   return fpt;
}
/*============================================================================*/
# undef READ_REVERSE_
# pragma  OPTIMIZE OFF
/***************** end of normalization function ffnrm(*) *********************/
