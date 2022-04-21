/* [ file: subarc.h ] */
# define DO_ARGZ "argz(*)"
/*******************************************************************************
*                                                                              *
*   ANSI-C function argz(*)                                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   This function returns arg(z) ( the argument or 'phase' of z ),             *
*   for any complex  z = x + i*y.                                              *
*                                                                              *
*   The value  of  arg(z)  is returned in the interval                         *
*   [ - PI ,  PI [ + 2PI*int(n/2), if n is even, or in                         *
*   [    0 , 2PI [ + 2PI*int(n/2), if n is odd and n>0,                        *
*   [ -2PI ,   0 [ + 2PI*int(n/2), if n is odd and n<0.                        *
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
   # define PRECISION (( double)( 1.000e-15 ))
# endif
/*============================================================================*/

double argz( double x, double y, short n ) 
{
   static double 
      u = ZERO,
      v = ZERO,
      norm = ZERO,
      gamma = ZERO;

   double fabs( double s );
   double atan( double s );
   double sqrt( double s );
/*----------------------------------------------------------------------------*/

   norm = sqrt( x*x + y*y );
     
   if ( norm < EPSILON )
   {
      gamma = ZERO;
      goto value; 
   };

   u = x/norm;
   v = y/norm;
 
   if ( EPSILON < u ) 
   {
      gamma = atan( v/u );
      goto value;
   } 
   else if ( u < -EPSILON )
   {
      if ( ZERO < v ) 
         gamma = atan( v/u ) + PI;
      if ( v <= ZERO ) 
         gamma = atan( v/u ) - PI;
      goto value;
   }
   else
   {
      if ( EPSILON < v )
      {
         gamma = PI/2.;
         goto value;
      }
      else if ( v < -EPSILON )
      {
         gamma = - PI/2.;
         goto value;
      }
      else
      {
         gamma = ZERO;
         goto value;
      };
   };
/*............................................................................*/
  value:
    
   if ( n != TWO*( int )( n/TWO ))          /* case: n odd */
   {
      if (( null < n )&&( gamma <  ZERO )) 
         gamma += 2.*PI;
      if (( n < null )&&( ZERO <= gamma )) 
         gamma -= 2.*PI;
   };
   gamma += 2.*PI*( int )( n/TWO );

   return gamma;
}
/*============================================================================*/
# undef EPSILON
# undef PRECISION
/************************** end of function arg(z) ****************************/





# define DO_SUBARC "subarc(*)"
/*******************************************************************************
*                                                                              *
*   ANSI-C  function subarc(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given two points P=(px,py) , Q=(qx,qy), an angle alfa in the  o p e n      *
*   interval ]-2PI,2PI[ and given any real ratio r in the interval [0,1],      *
*   this  subroutine  returns  the point  R = (arc.rx,arc.ry)  on the arc      *
*   from P to Q of total curvature alfa, which devides PQ into two subarcs     *
*   PR and RQ , the lengths ratio PR to RQ being equal to r / (1-r) .          *
*   - Further geometric objects are returned in structure arc.                 *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/* header files of called functions: */
/*
# include "../math/argz.h" */ /* function argz(*) */
/*----------------------------------------------------------------------------*/
# ifndef EPSILON
   # define EPSILON (( double )( 1.000e-277 ))
# endif
/*----------------------------------------------------------------------------*/
# ifndef PRECISION
   # define PRECISION (( double)( 1.000e-15 ))
# endif
/*----------------------------------------------------------------------------*/
struct subarc   
{
   double r, s, u, rx, ry, mx, my, nx, ny, xx, yy, arg;
}; 
static struct subarc arc = {ZERO};
/*============================================================================*/

double subarc( double px, double py, double qx, double qy, 
               double alfa, double rr )
{
/* allusions: */
/*
   extern struct 
      subarc arc;
*/
/* declarations: */

   static double 
      ss = ZERO, 
      xx = ZERO, 
      yy = ZERO,
      beta = ZERO,
      gamma = ZERO;

   static const double
      bound0 = 7.77*PRECISION,
      bound1 = 77.7*PRECISION,
      bound2 = 77.7*PRECISION;

   static signed char 
      sgn = null;

   double fabs( double x );
   double sqrt( double x );
   double  cos( double x );
   double  sin( double x ); 
   double  tan( double x );
   double atan( double x );
   double argz( double x , double y , short n );
/*----------------------------------------------------------------------------*/

   if (( alfa < - 2*PI + EPSILON )||( 2*PI - EPSILON < alfa )) 
   {
      fprintf( stderr, "\n\n Error message from function %s:", DO_SUBARC );
      fprintf( stderr, "\n Illegal transferred angle, |alfa| >= 2*PI !!!" );  
      fprintf( stderr, "\n [ required:  -2*PI < alfa < 2*PI.]\n " ); 

      exit( EXIT_FAILURE );
   };

/* euclidean distance |Q-P|: arc.s */

   xx = qx - px;
   yy = qy - py;

   arc.s = sqrt( xx*xx + yy*yy );
                                            
/*............................................................................*/
   arc.arg = argz ( xx, yy, ONE );    /* argument ( direction ) of vector Q-P */
/*..................................*/

/* unit vector orthogonal to Q-P [ in positive orientation ] */

   gamma = arc.arg + PI/2.;

   arc.nx = cos( gamma );
   arc.ny = sin( gamma );

   if ( fabs( alfa ) < bound1 )  /* straight line (P,Q) */
   {
      arc.r  = ( double ) HUGE_VAL;
      arc.u  = ( double ) HUGE_VAL;
      arc.mx = ( double ) HUGE_VAL;
      arc.my = ( double ) HUGE_VAL;

      arc.rx = px + rr*( qx - px );
      arc.ry = py + rr*( qy - py );
         
      return rr;
   }
   else /* if bound1 <= fabs( alfa ) */ /* true arc PQ */
   {

/* arc.u: distance, center ( arc.mx, arc.my ) of circle through P,Q to segmt.PQ
   arc.r: radius of that circle */

      if (( fabs( alfa - PI ) < bound2 )||( fabs( alfa + PI ) < bound2 ))
      {
         arc.u = ZERO;
         arc.r = arc.s/2.;
      }
      else
      {
         ss = arc.s/2.;
         arc.u = ss / tan( alfa/2. );
         arc.r = sqrt( arc.u*arc.u + ss*ss );
      };

/* circle midpoint: ( arc.mx, arc.my ) */

      arc.mx = ( px + qx )/2. + arc.u*arc.nx ;
      arc.my = ( py + qy )/2. + arc.u*arc.ny ;

      sgn = alfa/fabs( alfa );
      beta = arc.arg - sgn*PI/2. + ( rr -.5 )*alfa;

      if ( arc.s < EPSILON )
      {
         arc.rx = px;
         arc.ry = py;

         return rr;
      }
      else if ( fabs( rr ) < bound0 )
      {
         arc.rx = px;
         arc.ry = py;
      }
      else if ( fabs( rr - 1. ) < bound0 )
      {
         arc.rx = qx;
         arc.ry = qy;
      }
      else
      {
         arc.rx = arc.mx + arc.r*cos( beta );
         arc.ry = arc.my + arc.r*sin( beta );
      };

/* z_ratio: */
/* projection of vector R-P on Q-P: */

      ss = arc.rx - px;
      yy = arc.ry - py;
      xx = ss + ( ss*arc.nx + yy*arc.ny )*arc.nx;
      yy += ( ( ss*arc.nx + yy*arc.ny )*arc.ny );

      ss = sqrt( xx*xx + yy*yy )/arc.s ;      /* length ratio  |R'-P| / |Q-P| */

      return ss;
   };
}
/*============================================================================*/
# undef EPSILON
# undef PRECISION
/************************* end of funtion subarc(*) ***************************/
