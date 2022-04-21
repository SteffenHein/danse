/* [ file: rndoff.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function rndoff(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   This function returns ( double ) xx round off to ( short ) nn digits       *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifndef SMALL_VAL
   # define SMALL_VAL ( 1.e-277 )
# endif
/*----------------------------------------------------------------------------*/
# define USE_ROUND 0
/*============================================================================*/

double rndoff( double xx, short nn )
{
   static double
      ss = ZERO;

   static short
      ii = null;

   double pow( double xx, double yy );

# if USE_ROUND == 1
   double round( double xx );
# endif
/*----------------------------------------------------------------------------*/
   ss = fabs( xx );

   if ( ss < SMALL_VAL )
      return ZERO;

   ii = log10( ss );

   ss *= ( pow( 10., -ii ));
   ss *= ( pow( 10., nn ));

# if USE_ROUND == 1
   ss = round( ss );
# else
   ss += (( double ) HALF );
   ss = ( double )( long )( ss );
# endif

   ss *= ( pow( 10., -nn ));
   ss *= ( pow( 10., ii ));

   if ( xx < ZERO )
      ss = -ss;

   return ss;
}
/*============================================================================*/
/************************ end of function rndoff(*) ***************************/

