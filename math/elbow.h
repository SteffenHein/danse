/* [ file: elbow.h ] */
# define DO_ELBOW "elbow(*)"
/*******************************************************************************
*                                                                              *
*   ANSI-C function elbow(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given a line by two points P, Q, and an ellipse by its center C, the two   *
*   half axes a, b, and the direction arg ( argument ) of the principle axis a,*
*   this function returns, for any transferred point X = ( xx, yy ) :          *
*                                                                              *
*   (i)  the intersection point S1 of the ray ( C, X ) with the ellipse        *
*        determined by ( a, b, C, arg ) ,                                      *
*                                                                              *
*   (ii) the intersection point S2 of said ray with the line ( P, Q ),         *
*                                                                              *
*   as well as                                                                 *
*                                                                              *
*   (iii) the point R such that                                                *
*                                                                              *
*                     CX * CS1  =  CR * CS2   .                                *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
# include "../math/argz.h"
*/
typedef struct
{
   char
      rtn;

   double 
      a, b, arg, cx, cy,
      xx, yy,
      px, py,
      qx, qy,
      rx, ry,
      sx1, sy1,
      sx2, sy2,
      phi,
      nm1,nm2;
} ELBOW;
static ELBOW elb = {null};
/*============================================================================*/

ELBOW *\
elbow ( ELBOW *ept )
{
/* allusions: */
/*
   extern ELBOW elb;
*/
   static ELBOW *rpt = &elb;

   static double 
      xx1 = ZERO,
      yy1 = ZERO,

      sx1 = ZERO,
      sy1 = ZERO,
      sx2 = ZERO,
      sy2 = ZERO,

      px = ZERO,
      py = ZERO,
      qx = ZERO,
      qy = ZERO,
      rx = ZERO,
      ry = ZERO,

      sina = ZERO,
      cosa = ZERO,
      sinb = ZERO,
      cosb = ZERO,
      sinc = ZERO,
      cosc = ZERO,

      uu = ZERO,
      vv = ZERO,

      epsilon = ZERO;

   double sin ( double x );
   double cos ( double x );
   double fabs( double x );

/*----------------------------------------------------------------------------*/

   if ( ept == null ) /* initialize: */
   { 
      rpt->a = ZERO;
      rpt->b = ZERO;
      rpt->arg = ZERO;
      rpt->cx = ZERO;
      rpt->cy = ZERO;

      rpt->xx = ZERO;
      rpt->yy = ZERO;

      rpt->px = ZERO;
      rpt->py = ZERO;
      rpt->qx = ZERO;
      rpt->qy = ZERO;

      rpt->rx = ZERO;
      rpt->ry = ZERO;
      rpt->sx1 = ZERO;
      rpt->sy1 = ZERO;
      rpt->sx2 = ZERO;
      rpt->sy2 = ZERO;

      rpt->phi = ZERO;
      rpt->nm1 = ZERO;
      rpt->nm2 = ZERO;

      return rpt;
   };
/*............................................................................*/

   rpt = ept;

   epsilon = 1.e-277;

   if (( epsilon < ept->a )&&( epsilon < ept->b ))
   {
      uu = ept->xx - ept->cx;
      vv = ept->yy - ept->cy;

      sina = sin( ept->arg );
      cosa = cos( ept->arg );

      xx1 =   uu * cosa + vv * sina;
      yy1 = - uu * sina + vv * cosa;

      ( rpt->phi ) = argz( xx1, yy1, null );

      sinb = sin( rpt->phi );
      cosb = cos( rpt->phi );

      uu = ( ept->a ) * sinb;
      vv = ( ept->b ) * cosb;
   
      ( rpt->nm1 ) = (( ept->a ) * ( ept->b )) / sqrt( uu*uu + vv*vv );

      sx1 = rpt->nm1 * cosb;
      sy1 = rpt->nm1 * sinb;
      
      ( rpt->sx1 ) = sx1 * cosa - sy1 * sina + ept->cx;
      ( rpt->sy1 ) = sx1 * sina + sy1 * cosa + ept->cy;
   }
   else
   {
      printf( "\n\n Error message from function '%s': ", DO_ELBOW );
      printf( "\n Unspecified or illegal half axes !!! " );
      printf( "\n [ specify half axes before function call, " );
      printf( "\n   in struct ELBOW *ept: ept->a := a, ept->b := b .] \n\n " );
      ( rpt->rtn ) = ONE;
      return rpt;
   };
   
   uu = rpt->phi + PI/2.;

   sinc = sin( uu );
   cosc = cos( uu );

   uu = ept->px - ept->cx;
   vv = ept->py - ept->cy;

   px =   uu * cosa + vv * sina;
   py = - uu * sina + vv * cosa;

   uu = ept->qx - ept->cx;
   vv = ept->qy - ept->cy;

   qx =   uu * cosa + vv * sina;
   qy = - uu * sina + vv * cosa;

   uu = ( qx - px )*cosc + ( qy - py )*sinc;

   if( epsilon < fabs( uu ))
   {
      vv = qx * cosc + qy * sinc;
      vv /= uu;

      sx2 = qx + vv * ( px - qx );
      sy2 = qy + vv * ( py - qy );

      rpt->nm2 = sqrt( sx2 * sx2 + sy2 * sy2 );

      rpt->sx2 = sx2 * cosa - sy2 * sina + ept->cx;
      rpt->sy2 = sx2 * sina + sy2 * cosa + ept->cy;

      rx = xx1 * rpt->nm1 / rpt->nm2;
      ry = yy1 * rpt->nm1 / rpt->nm2;
      
      rpt->rx = rx * cosa - ry * sina + ept->cx;
      rpt->ry = rx * sina + ry * cosa + ept->cy;
   };

   ( rpt->rtn ) = null;

   return rpt;
}
/*============================================================================*/
/************************* end of function elbow(*) ***************************/
