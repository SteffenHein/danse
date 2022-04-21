/* [ file: zexpr.h ] */
# define DO_ZEXPN "zexpr(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function zexpr(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   This function computes z^rr for any real rr and complex                    *
*                                                                              *
*                   zz = ( cpt->r ) + j*( cpt->i ).                            *
*                                                                              *
*   The result is written into a structure of type COMPLEX,                    *
*   a pointer to that is returned to the calling program.                      *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
typedef struct
{
   double r, i;

   double arg, nrm;
}  COMPLEX;
# include "../math/consts.h"
*/
# include "../math/argc.h"
/*============================================================================*/

COMPLEX *
zexpr( COMPLEX *cpt, double rr ) 
{
   static COMPLEX z = {null},
                 *zpt = &z;

   static double 
     uu = ZERO,
     vv = ZERO,
     xx = ZERO,
     yy = ZERO;

   double 
      exp( double x ),
      sin( double x ), 
      cos( double x ); 

   COMPLEX *
      argc( COMPLEX *zpt, short n );

/*----------------------------------------------------------------------------*/

   zpt->r = cpt->r;
   zpt->i = cpt->i;  

/*............................................................................*/
   zpt = argc( zpt, null );   /* determines modulus and phase of z            */
/*..........................*/

   uu = zpt->nrm;

   if ( 1.e-277 < uu )
   {
      vv = zpt->arg;
      uu = exp( rr * log(uu ));
      zpt->r = uu*cos( rr*vv );
      zpt->i = uu*sin( rr*vv );

/*............................................................................*/
      zpt = argc( zpt, null );   /* determines modulus and phase of z^rr      */
/*.............................*/
   }
   else
   {
      zpt->r = ZERO;
      zpt->i = ZERO; 
   };

   return zpt;
}
/*============================================================================*/
/************************* end of function zexpr(*) ***************************/
