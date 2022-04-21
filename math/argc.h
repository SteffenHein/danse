/* [ file: argc.h ] */
# define DO_ARGC "argc(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function argc(*)                                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   This function returns arg(z) ( the argument or 'phase' of z ),             *
*   for any complex  z = uu + i*vv.                                            *
*                                                                              *
*   The value  of  arg(z)  is returned in the interval                         *
*   [ - PI ,  PI [ + 2PI*int(n/2), if nn is even, or in                        *
*   [    0 , 2PI [ + 2PI*int(n/2), if nn is odd and nn>0,                      *
*   [ -2PI ,   0 [ + 2PI*int(n/2), if nn is odd and nn<0.                      *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifndef EPSILON
   # define EPSILON (( double )( 1.000e-277 ))
# endif
/*----------------------------------------------------------------------------*/
# ifndef PRECISION 
   # define PRECISION (( double )( 1.000e-15 ))
# endif
/*============================================================================*/

COMPLEX *\
argc( COMPLEX *ipt, short nn )
{
   static COMPLEX cc = {null},
                *rpt = &cc;

   static double 
   uu = ZERO,
   vv = ZERO;

   double 
      fabs( double s ),
      atan( double s ),
      sqrt( double s );
/*----------------------------------------------------------------------------*/
   uu = ( ipt->r );
   vv = ( ipt->i );

   ( rpt->r ) = uu;
   ( rpt->i ) = vv;

   ( rpt->nrm ) = sqrt( uu*uu + vv*vv );
     
   if (( rpt->nrm ) < EPSILON )
   {
      ( rpt->arg ) = ZERO;
      goto value; 
   };

   uu /= ( rpt->nrm );
   vv /= ( rpt->nrm );
 
   if ( EPSILON < uu ) 
   {
      ( rpt->arg ) = atan( vv/uu );
      goto value;
   } 
   else if ( uu < -EPSILON )
   {
      if ( ZERO < vv ) 
         ( rpt->arg ) = atan( vv/uu ) + PI;
      if ( vv <= ZERO ) 
         ( rpt->arg ) = atan( vv/uu ) - PI;
      goto value;
   }
   else
   {
      if ( EPSILON < vv )
      {
         ( rpt->arg ) = PI/2.;
         goto value;
      }
      else if ( vv < -EPSILON )
      {
         ( rpt->arg ) = -PI/2.;
         goto value;
      }
      else
      {
         ( rpt->arg ) = ZERO;
         goto value;
      };
   };
/*............................................................................*/
  value:
    
   if ( nn != TWO*( short )( nn/TWO ))         /* case: nn odd */
   {
      if (( null < nn )&&(( rpt->arg ) <  ZERO )) 
         ( rpt->arg ) += 2.*PI;
      if (( nn < null )&&( ZERO <= ( rpt->arg ))) 
         ( rpt->arg ) -= 2.*PI;
   };
   ( rpt->arg ) += ( 2.*PI*( short )( nn/TWO ) );

   return rpt;
}
/*============================================================================*/
# undef EPSILON
# undef PRECISION
/************************* end of function argc(*) ****************************/
