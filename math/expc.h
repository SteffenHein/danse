/* [ file: expc.h ] */
# define DO_EXPC "expc(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function expc(*)                                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   This function computes the exponential function exp(z), for any complex    *
*   number                                                                     *
*                               z = uu + j*vv                                  *
*                                                                              *
*   The result is written into a structure of type COMPLEX,                    *
*   and a pointer to that structure is returned to the calling program.        *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "argc.h"
/*============================================================================*/

COMPLEX *
expc( COMPLEX *ipt )
{
   static COMPLEX cc = {null},
                 *rpt = &cc;

   static double 
     xx = ZERO,
     yy = ZERO,
     rr = ZERO;

   double exp( double x );
   double sin( double x ); 
   double cos( double x ); 
   double fabs( double x );

   COMPLEX *argc( COMPLEX *ipt, short n );
/*----------------------------------------------------------------------------*/

   xx = cos( ipt->i );
   yy = sin( ipt->i );
   rr = exp( ipt->r );

   (rpt->r) = rr*xx;
   (rpt->i) = rr*yy;

   xx = ( argc( rpt, null ) -> arg );
   (rpt->arg) = xx;
   rr = fabs( rr );
   (rpt->nrm) = rr;

   return rpt;
}
/*============================================================================*/
/************************* end of function expc(*) ****************************/
