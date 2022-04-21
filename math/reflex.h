/* [ file: reflex.h ] */
# define DO_REFLEX
/*******************************************************************************
*                                                                              *
*   ANSI C function reflex(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given a point p = ( px, py ), another point v = ( vx, vy ), and an angle   *
*   alfa [ in radians ] this routine determines the point q = ( qx, qy ) in    *
*   the plane that lies symmetrical to p with respect to the straight line     *
*   through v in direction alfa [ measured from the positive x axis ].         *
*   [ All variables in struct REFLEX reflx.]                                   *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*----------------------------------------------------------------------------*/
typedef struct
{
   char 
      rtn; 

   double
      px, py, qx, qy, vx, vy, alfa;

} REFLEX;
static REFLEX 
   reflx = {null};
/*============================================================================*/

REFLEX *reflex( REFLEX *rxp )
{
   static REFLEX
     *rpt = &reflx;

   static double 
      alfa = ZERO,
      xx1 = ZERO,
      yy1 = ZERO,
      xx2 = ZERO,
      yy2 = ZERO,
      ss = ZERO,
      tt = ZERO;

   static const double
      bound = 1.e-17;
/*----------------------------------------------------------------------------*/
/* initialize if rxp is null pointer: */

   if( rxp == null )
   {
      ( rpt->px ) = ZERO;
      ( rpt->py ) = ZERO;
      ( rpt->qx ) = ZERO;
      ( rpt->qy ) = ZERO;
      ( rpt->vx ) = ZERO;
      ( rpt->vy ) = ZERO;
      ( rpt->alfa ) = ZERO;

      ( rpt->rtn ) = null;
      return rpt;
   };
/*............................................................................*/
   alfa = rxp->alfa;

   while( alfa < ZERO ) /* shift alfa into interval [0, 2*PI) */
   {
      alfa += 2.*PI;
   }
   while( 2.*PI <= alfa )
   {
      alfa -= 2.*PI;
   };

   xx1 = ( rxp->px ) - ( rxp->vx );
   yy1 = ( rxp->py ) - ( rxp->vy );
   ( rpt->px ) = rxp->px;
   ( rpt->py ) = rxp->py;
      
   if(( fabs( alfa ) < bound )||
      ( fabs( alfa - PI ) < bound ))
   {
      xx2 = xx1;
      yy2 = -yy1;
      alfa = ZERO;
   }
   else if(( fabs( alfa - PI/2. ) < bound )||
           ( fabs( alfa - 3.*PI/2. ) < bound ))
   {
      xx2 = -xx1;
      yy2 = yy1;
      alfa = PI/2.;
   }
   else
   {
      ss = tan( alfa );
      tt = -1./ss;  /* implies: ss - tt != ZERO */

      xx2 = ( 2.*yy1 - ( ss + tt )*xx1 )/( ss - tt );
      yy2 = ss*( xx2 + xx1 ) - yy1;
   };

   ( rpt->alfa ) = alfa;

   ( rpt->px ) = rxp->px;
   ( rpt->py ) = rxp->py;
   ( rpt->qx ) = xx2 + rxp->vx;
   ( rpt->qy ) = yy2 + rxp->vy;
   
   return rpt;
}
/*============================================================================*/
# undef DO_REFLEX
/************************* end of function reflex(*) **************************/
