/* [ file: trnsvs.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function trnsvs(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Computes the projection of any given real vector trv[] onto the plane      *
*   transverse to the direction of another given real vector dir[].            *
*   Thereby trv[] is overwritten by that projection, i.e.                      *
*          trv[] := trv[] - dir[] * < trv[] | dir[] >/< dir[] | dir[] >.       *
*   Also, the scalar product of the originally transferred vector trv[]        *
*   with the directional vector dir[],                                         *
*                     scpr = < trv[] | dir[] >,                                *
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

double trnsvs( double *trv, double *dir, int nn )
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
         scpr += ( trv[kk]*dir[kk] );
      } while(( ++kk ) < nn );
      if ( SMALL_VAL < nrm2 )
      {
         nrm2 = scpr/nrm2;
         while( null < ( kk-- ))
            trv[kk] -= ( dir[kk]*nrm2 );
      }
      else
      {
         while( null < ( kk-- ))
            trv[kk] = ZERO;
      };
   };

   return scpr;
}
/*============================================================================*/
/***************************** end of file trnsvs.c ***************************/
