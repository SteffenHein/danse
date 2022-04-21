/* [ file: fftrf.c ] */
# define DO_FFTRF "fftrf(*)" 
/*******************************************************************************
*                                                                              *
*   ANSI C function fftrf(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   Fast Fourier Transformation                                                *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <time.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma  OPT_LEVEL 2
   # pragma  OPTIMIZE ON
# endif
# ifdef NO_SIDE_EFFECTS
   # pragma  NO_SIDE_EFFECTS fftrf.c
# endif
/*----------------------------------------------------------------------------*/
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H"
/*----------------------------------------------------------------------------*/
/* Edit and customize this header for POSTER configuration: */
# include "../poster/POSTER.CONF"
/*----------------------------------------------------------------------------*/
/*
# include "../math/FFT.CONF"
*/
/*----------------------------------------------------------------------------*/
/* the following two constants may yet been defined [ in main program FFTRN.C,*/
/* for instance ] */

/* FTR_NMBR = maximum number of distinct distributions to be simultaneously */
/* transformed [ e.g. for convolutions ]: */

# ifndef FTR_NMBR
   # define FTR_NMBR 10 
# endif

/* FTR_SIZE = maximum number of sample points */
/* [ maximum array length for complex values, preferably FTR_SIZE = 2^N */
/* e.g. FTR_SIZE = 2^16 = 65536 ]: */

# ifndef FTR_SIZE         
   # define FTR_SIZE 65536 
# endif

/* FTR_NORM = 1 yields normalized output ( fpt->r ), ( fpt->i ) */
# ifndef FTR_NORM
   # define FTR_NORM 0
# endif

/* FTR_DISP = 1 displays additional comments */
# ifndef FTR_DISP
   # define FTR_DISP 1
# endif
/*----------------------------------------------------------------------------*/
typedef struct
{                     
   signed char
      rtn;

   char
      opt[SHS_SIZE];

   short 
      p, q,
      mult[FTR_NMBR+ONE];

   long
      ttlg[FTR_NMBR+ONE], stlg[FTR_NMBR+ONE];

/* distributions, real and imaginary parts: */

   double
      r[FTR_NMBR+ONE][FTR_SIZE+ONE],
      i[FTR_NMBR+ONE][FTR_SIZE+ONE],

      t[FTR_NMBR+ONE], tt[FTR_NMBR+ONE], dt[FTR_NMBR+ONE],
      s[FTR_NMBR+ONE], ss[FTR_NMBR+ONE], ds[FTR_NMBR+ONE];
/* normalization constant: */

   double
      nor;

} FFT;
/*============================================================================*/

FFT * \
fftrf( FFT *fpt )
{
/* declarations: */

   static double
      wr = ZERO,
      wi = ZERO,
      wpr = ZERO,
      wpi = ZERO,
      real = ZERO,
      imag = ZERO,
      theta = ZERO,
      dx = ZERO,
      nyquist = ZERO;

# if FTR_NORM == 1  
   static double
      rnorm = ZERO;
# endif

   static long
      ii = null,
      jj = null,
      mm = null,
      nn = null,
      mmax = null, 
      istep = null,
      samp = null;

   static short
      hh = null, 
      kk = null;

   static signed char isign = null;

/* function prototypes: */

   double cos( double x );
   double sin( double x );
   double sqrt( double x );

/*----------------------------------------------------------------------------*/
/* parameter transfer errors: */

   if (( fpt->p ) < ONE )
   {
      printf( "\n\n Error message from function %s:", __func__ );
      printf( "\n Illegal transfered index ( fpt->p ) = %ld !!!",
         ( long ) ( fpt->p ));
      printf( "\n [ ( fpt->p ) must be a positive integer.]\n " );

      ( fpt->rtn ) = ONE;
      return fpt;
   };

   if (( fpt->q ) < ONE )
   {
      printf( "\n\n Error message from function %s:", __func__ );
      printf( "\n Illegal transfered index ( fpt->q ) = %ld !!!",
         ( long ) ( fpt->q ));
      printf( "\n [( fpt->q ) must be a positive integer.]\n " );

      ( fpt->rtn ) = ONE;
      return fpt;
   };

   if (( fpt->q ) < ( fpt->p ))
   {       
      printf( "\n\n Error message from function %s:", __func__ );
      printf( "\n Illegal transfered index ( fpt->q ) = %ld < "
         "( fpt->p ) = %ld !!!", ( long )( fpt->q ), ( long )( fpt->p ));
      printf( "\n [ required: ( fpt->p ) <= ( fpt->q ).]\n " );

      ( fpt->rtn ) = ONE;
      return fpt;
   };
/*----------------------------------------------------------------------------*/
/* option: forward transform */ 

   if (( *( fpt->opt ) != 'b' )
     &&( *( fpt->opt ) != 'B' ))
   {
      isign = -ONE;

      samp = ( fpt->ttlg[( fpt->p )] ); /* samp = number of values */
                                        /* = number of intervals + 1 */
      nn = TWO ; 
      while ( nn < samp )
      {
       /* nn <<= 1; *//* equivalent to ( nn *= 2 ) [ i.e. ( nn = 2*nn ) ] */
          nn *= TWO;
      };

      if ( FTR_SIZE < nn )
      {
         printf( "\n\n Message from function %s:", __func__ );
         printf( "\n Too many values on input file !!!\n " );

         ( fpt->rtn ) = ONE;
         return fpt;
      };

      dx = ( fpt->dt[( fpt->p )] );
      nyquist = .5/dx; 

# if FTR_NORM == 1  
      rnorm = dx;
# endif

      ( fpt->stlg[null] ) = nn; /* nn = intervals + 1 = no. of sample points */
      ( fpt->ds[null] ) = 2.*nyquist/nn;
      ( fpt->s[null] ) = - nyquist;
      ( fpt->ss[null] ) = nyquist;

      hh = ( fpt->p ); do
      {
         if ((( fpt->ttlg[hh] ) < TWO )
           ||(( fpt->dt[hh] ) <= ZERO )
           ||(( fpt->tt[hh] ) <= ( fpt->t[hh] )))
         {
            printf( "\n Message from function %s:", __func__ );
            printf( "\n Pathological distribution no. %d !!!", hh );
            printf( "\n [ please verify ]. " );

            ( fpt->rtn ) = ONE;
            return fpt;
         };

         if ((( fpt->ttlg[hh] ) != samp )
           ||(( fpt->dt[hh] ) != dx ))
         {
            printf( "\n Message from function %s:", __func__ );
            printf( "\n Incompatible distributions %ld and %ld !!!",
               ( long ) ( fpt->p ), ( long ) hh );
            printf( "\n [ Check transferred distribution parameters.]\n " );

            ( fpt->rtn ) = ONE;
            return fpt;
         };

         ( fpt->stlg[hh] ) = ( fpt->stlg[null] );
         ( fpt->ds[hh] ) = ( fpt->ds[null] );
         ( fpt->ss[hh] ) = ( fpt->ss[null] );
         ( fpt->s[hh] ) = ( fpt->s[null] );
      } while (( ++hh ) <= ( fpt->q ));
   }
/*----------------------------------------------------------------------------*/
/* option: back transform */

   else if (( *( fpt->opt ) == 'b' )||( *( fpt->opt ) == 'B' ))
   {
      isign = ONE;

      samp = ( fpt->stlg[( fpt->p )] ); /* samp = number of values */
                                        /* = number of intervals + 1 */
      nn = TWO ;
      while ( nn < samp )
      {
         /* nn <<= 1; */ /* equivalent to ( nn *= 2 ) [ i.e. ( nn = 2*nn ) ] */
         nn *= TWO;
      };

      if ( FTR_SIZE < nn )
      {
         printf( "\n\n Message from function %s:", __func__ );
         printf( "\n Too many values on input file !!!\n " );

         ( fpt->rtn ) = ONE;
         return fpt;
      };

      dx = ( fpt->ds[( fpt->p )] );
      nyquist = .5/dx;

# if FTR_NORM == 1  
      rnorm = dx;
# endif

      ( fpt->ttlg[null] ) = nn; 
      ( fpt->dt[null] ) = 2.*nyquist/nn;
      ( fpt->tt[null] ) = nyquist;
      ( fpt->t[null] ) = -nyquist;

      hh = ( fpt->p ); do
      {
         if ((( fpt->stlg[hh] ) < TWO )
           ||(( fpt->ds[hh] ) <= ZERO )
           ||(( fpt->ss[hh] ) <= ( fpt->s[hh] )))
         {
            printf( "\n Message from function %s:", __func__ );
            printf( "\n Pathological spectrum no. %ld !!!", ( long ) hh );
            printf( "\n [ Check transferred spectrum parameters.]\n " );

            ( fpt->rtn ) = ONE;
            return fpt;
         };

         if ((( fpt->stlg[hh] ) != samp )
           ||(( fpt->ds[hh] ) != dx ))
         {
            printf( "\n Message from function '%s':", __func__ );
            printf( "\n Incompatible spectra %ld and %ld !!!",
               ( long ) ( fpt->p ), ( long ) hh );
            printf( "\n [ Check transferred spectrum parameters.]\n " );

            ( fpt->rtn ) = ONE;
            return fpt;
         };

         ( fpt->ttlg[hh] ) = ( fpt->ttlg[null] );
         ( fpt->dt[hh] ) = ( fpt->dt[null] );
         ( fpt->tt[hh] ) = ( fpt->tt[null] );
         ( fpt->t[hh] ) = ( fpt->t[null] );
      } while (( ++hh ) <= ( fpt->q ));
   };
/*............................................................................*/
/* common FFtransform section: */

   hh = ( fpt->p ); do
   {
      ii = samp ;
      while ( ii < nn )
      {
         ( fpt->r[hh][ii] ) = ZERO;
         ( fpt->i[hh][ii] ) = ZERO;
         ii++ ;
      };
/*............................................................................*/
# if FTR_DISP == 1
      if (( *( fpt->opt ) != 'b' )&&( *( fpt->opt ) != 'B' ))
      {
         printf( "\n --------- %2d. Fourier transformation started; "
            "please don't disturb ! ---------", hh );
      }
      else if (( *( fpt->opt ) == 'b' )||( *( fpt->opt ) == 'B' ))
      {
         printf( "\n ---- %2d. Fourier ( back - ) transformation started; "
            "please don't disturb ! ---", hh );
      };
# endif
/*............................................................................*/
/*  Fast Fourier Transform algorithm [ Cooley - Tukey ]: */

      jj = null;
/*............................................................................*/
/* bit reversal section: */

      ii = null; do
      {
         if ( ii < jj ) /* Exchange the two complex numbers: */
         {
            real = ( fpt->r[hh][jj] ); 
            imag = ( fpt->i[hh][jj] );

            ( fpt->r[hh][jj] ) = ( fpt->r[hh][ii] );
            ( fpt->i[hh][jj] ) = ( fpt->i[hh][ii] ); 
            ( fpt->r[hh][ii] ) = real;
            ( fpt->i[hh][ii] ) = imag;
         }; 

         mm = nn/TWO;
            
         while (( ONE <= mm )
              &&( mm <= jj ))
         {
            jj -= mm;
           /* mm >>= 1;*/ /* equivalent to ( mm /= 2 ) [ i.e. ( mm = mm/2 ) ] */
            mm /= TWO;
         };
         jj += mm;
      } while (( ++ii ) < nn );
/*............................................................................*/
/* Danielson - Lanczos section: */

      mmax = ONE; do
      {
         istep = TWO*mmax;
         theta = PI/( isign*mmax );
           wpr = cos( theta );
           wpi = sin( theta );

         wr = 1.; 
         wi = ZERO;

         mm = null; do
         {
            ii = mm; do 
            {
               jj = ii + mmax;

               real = wr*( fpt->r[hh][jj] ) - wi*( fpt->i[hh][jj] );
               imag = wr*( fpt->i[hh][jj] ) + wi*( fpt->r[hh][jj] ); 

               ( fpt->r[hh][jj] ) = ( fpt->r[hh][ii] ) - real;
               ( fpt->i[hh][jj] ) = ( fpt->i[hh][ii] ) - imag; 
               ( fpt->r[hh][ii] ) += real;
               ( fpt->i[hh][ii] ) += imag;

               ii += istep; 
            } while ( ii < nn );

            real = wr;                      
              wr =  wr*wpr - wi*wpi;
              wi =  wi*wpr + real*wpi;
         } while (( ++mm ) < mmax );
         mmax = istep;
      } while ( mmax < nn );

# if FTR_NORM == 1  
      ii = null; do 
      {
         ( fpt->r[hh][ii] ) *= rnorm;
         ( fpt->i[hh][ii] ) *= rnorm;
      } while (( ++ii ) < nn );
# endif
/*............................................................................*/
# if FTR_DISP == 1
      if (( *( fpt->opt ) != 'b')&&( *( fpt->opt ) != 'B' ))
      {
         printf( "\r ------------------- %2d. Fourier transformation "
            "terminated --------------------", hh );
      } 
      else
      { 
/*
         printf( "\r -------------- %2d. Fourier ( back - ) "
            "transformation terminated --------------", hh );
*/
         printf( "\r ---------------- %2d. reverse Fourier "
            "transformation terminated ---------------", hh );
      };
# endif
   } while (( ++hh ) <= ( fpt->q ));
/*............................................................................*/
/* convolution: */

   ii = null; do       
   {
      ( fpt->r[null][ii] ) = 1.;
      ( fpt->i[null][ii] ) = ZERO;

      hh=( fpt->p ); do
      {
	 for ( kk=ONE; kk<=( fpt->mult[hh] ); kk++ )
	 {
	    real = ( fpt->r[null][ii] );
	    imag = ( fpt->i[null][ii] );

	    ( fpt->r[null][ii] ) = \
               real*( fpt->r[hh][ii] ) - imag*( fpt->i[hh][ii] );
	    ( fpt->i[null][ii] ) = \
               real*( fpt->i[hh][ii] ) + imag*( fpt->r[hh][ii] );
         }; /* next kk */ 
      } while (( ++hh ) <= ( fpt->q ));
   } while (( ++ii ) < nn );
/*............................................................................*/

   ( fpt->rtn ) = null;
   return fpt;
}
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma  OPTIMIZE OFF
# endif
/************ end of fast Fourier transformation function fftrf(*) ************/

