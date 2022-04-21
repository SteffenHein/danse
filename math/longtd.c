/* [ file: longtd.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function longtd(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Computes the component of any given real vector lgt[] in the direction     *
*   of any other given real vector dir[]. Thereby lgt[] is overwritten by      *
*   that component, i.e.                                                       *
*          lgt[] := dir[] * < lgt[] | dir[] >/< dir[] | dir[] >.               *
*   Also, the scalar product of the originally transferred vector lgt[]        *
*   with the directional vector dir[],                                         *
*                     scpr = < lgt[] | dir[] >,                                *
*   is returned.                                                               *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define _POSIX_SOURCE 0 /* set to 1: if POSIX.1 standard headers will be used*/
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifndef SMALL_VAL
   # define SMALL_VAL ( double )( 1.00e-277 )
# endif
# ifndef PRECISION
   # define PRECISION ( double )( 1.000e-15 )
# endif
# ifndef RNDOFFBND
   # define RNDOFFBND ( 5.*PRECISION )
# endif
/*============================================================================*/

double longtd( double *lgt, double *dir, int nn )
{
   static short
      kk = null;

   static double
      nrm2 = ZERO,
      scpr = ZERO;
/*----------------------------------------------------------------------------*/
   if ( nn <= null )
      return ZERO;
   else
   {
      nrm2 = ZERO;
      scpr = ZERO;
      kk = null; do
      {
         nrm2 += ( dir[kk]*dir[kk] );
         scpr += ( lgt[kk]*dir[kk] );
      } while(( ++kk ) < nn );
      if ( SMALL_VAL < nrm2 )
      {
         nrm2 = scpr/nrm2;
         while( null < ( kk-- ))
            lgt[kk] = dir[kk]*nrm2;
      }
      else
      {
         while( null < ( kk-- ))
            lgt[kk] = ZERO;
      };
   };

   return scpr;
}
/*============================================================================*/
/****************************** end of file longtd.c **************************/
