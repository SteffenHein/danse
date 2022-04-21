/* [ file: arcsct.h ] */
# define DO_ARGZ "argz(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function argz(*)                                                    *
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
# include "./maths.h"  /* 'my' computation environment headers */
# include "./consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifndef EPSILON
   # define EPSILON (( double )( 1.000e-277 ))
# endif
/*----------------------------------------------------------------------------*/
# ifndef PRECISION
   # define PRECISION (( double )( 1.000e-15 ))
# endif
/*============================================================================*/

double argz( double x , double y , short n ) 
{
   static double 
      u = ZERO,
      v = ZERO,
      norm = ZERO,
      gamma = ZERO;

   double fabs( double s );
   double atan( double s );
   double sqrt( double s );
        
/*============================================================================*/

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
/************************ end of function 'arg(z)' ****************************/





# define DO_SUBARC "subarc(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function subarc(*)                                                  *
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
# include "argz.h" */ /* function argz(*) */
/*----------------------------------------------------------------------------*/
# ifndef EPSILON
   # define EPSILON (( double )( 1.000e-277 ))
# endif
/*----------------------------------------------------------------------------*/
# ifndef PRECISION
   # define PRECISION (( double )( 1.000e-15 ))
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
/* declaration: */

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

/*============================================================================*/

   if (( alfa < - 2*PI + EPSILON )||( 2*PI - EPSILON < alfa )) 
   {
      fprintf( stderr, "\n\n Error message from function %s :", DO_SUBARC );
      fprintf( stderr, "\n Illegal transferred angle, |alfa| >= 2*PI !!!" );  
      fprintf( stderr, "\n [ required:  -2*PI < alfa < 2*PI . ]\n" ); 

      exit( EXIT_FAILURE );
   };

/* euclidean distance |Q-P|: arc.s */

   xx = qx - px;
   yy = qy - py;

   arc.s = sqrt( xx*xx + yy*yy );
                                            
/*............................................................................*/
   arc.arg = argz ( xx , yy , ONE );  /* argument ( direction ) of vector Q-P */
/*..................................*/

/* unit vector orthogonal to Q-P [ in positive orientation ] */

   gamma = arc.arg + PI/2.;

   arc.nx = cos( gamma );
   arc.ny = sin( gamma );

   if ( fabs( alfa ) < bound1 )  /* straight line (P,Q) */
   {
      arc.r  = HUGE_VALF;
      arc.u  = HUGE_VALF;
      arc.mx = HUGE_VALF;
      arc.my = HUGE_VALF;

      arc.rx = px + rr*( qx - px );
      arc.ry = py + rr*( qy - py );
         
      return rr;
   }
   else /* if fabs( alfa ) >= bound1 */  /* true arc PQ */
   {

/* arc.u: distance, center ( arc.mx,arc.my ) of circle through P,Q to segmt.PQ
   arc.r: radius of that circle                                               */

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
      beta = arc.arg - sgn*PI/2. + ( rr - .5 )*alfa;

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
/*************************** end of funtion subarc ****************************/










# define DO_ARCSCT "arcsct(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function arcsct(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given two arcs  AB  and  CD  with end points  A=(ax,ay), B=(bx,by) and     *
*   C=(cx,cy), D=(dx,dy), respectively, and with total curvatures alfa and     *
*   beta, this function computes the point of intersection S of AB and CD,     *
*   if this point exists and is unique.                                        *
*   Also in this case, the length ratio  r = |S'-A| / |B-A|  is computed,      *
*   where S' denotes the orthogonal projection of point S onto line (A,B).     *
*   The coordinates of S are written into  ( arc.xx, arc.yy ) of structure     *
*   subarc  [ cf. function subarc(*) ]  and the length ratio r is returned     *
*   to the calling program.                                                    *
*                                                                              *
*   Otherwise, if no point of intersection exists or if it is  non-unique,     *
*   -ONE is returned to the calling program [ and a message printed on the     *
*   standard output ].                                                         *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define MAXITER   1000
/*----------------------------------------------------------------------------*/
# ifndef EPSILON
   # define EPSILON (( double )( 1.000e-277 ))
# endif
/*----------------------------------------------------------------------------*/
# ifndef PRECISION
   # define PRECISION (( double )( 1.000e-15 ))
# endif
/*============================================================================*/

double arcsct( double ax, double ay, double bx, double by, double alfa,
               double cx, double cy, double dx, double dy, double beta )
{
/* allusions: */
/*
   extern struct 
      subarc arc;
*/
/* declarations: */

   static double
      rr = ZERO,
      ss = ZERO,

      xx = ZERO,
      yy = ZERO,

      x0 = ZERO,
      x1 = ZERO,
      x2 = ZERO,
      x3 = ZERO,

      y0 = ZERO,
      y1 = ZERO, 
      y2 = ZERO,
      y3 = ZERO,

      d0 = ZERO,
      d1 = ZERO,
      d2 = ZERO,
      d3 = ZERO,

      aa00 = ZERO,
      aa01 = ZERO,
      aa02 = ZERO, 
      aa10 = ZERO,
      aa11 = ZERO,
      aa12 = ZERO,

      s0 = ZERO,
      s1 = ZERO,

      mx = ZERO,
      my = ZERO,

      gamma0 = ZERO,
      gamma1 = ZERO,

      scale  = ZERO;

   static const double
      bound0 = 2.*PI*PRECISION,        /* absolute angular bound 'zero'       */
      bound1 = 7.*PRECISION,           /* point coincidence bound             */
      bound2 = 49.*PRECISION*PRECISION,/* square of ''                        */
      bound3 = 77.*PRECISION;          /* bound for Newton [secant] iteration */

   static short
       nn = null; 

   static signed char
      sgn = null;

/* prototypes: */

   double cos( double x );
   double sin( double x );
   double sqrt( double x );
   double fabs( double x );
   double argz( double x, double y, short n );

/*----------------------------------------------------------------------------*/
/* distances: */

   scale = ZERO;

   xx = bx - ax;
   yy = by - ay;

   rr = xx*xx + yy*yy;

   if ( scale < rr )
      scale = rr;

   xx = dx - cx;
   yy = dy - cy;

   rr = xx*xx + yy*yy;

   if ( scale < rr )
      scale = rr;

/*............................................................................*/

   xx = cx - ax;
   yy = cy - ay;

   d0 = xx*xx + yy*yy; /* |C-A|^2 */

   if ( scale < d0 )
      scale = d0;  

   xx = dx - ax;
   yy = dy - ay;

   d1 = xx*xx + yy*yy; /* |D-A|^2 */

   if ( scale < d1 )
      scale = d1;

   xx = cx - bx;
   yy = cy - by;

   d2 = xx*xx + yy*yy; /* |C-B|^2 */

   if ( scale < d2 )
      scale = d2;

   xx = dx - bx;
   yy = dy - by;

   d3 = xx*xx + yy*yy; /* |D-B|^2 */

   if ( scale < d3 )
      scale = d3;

/* check point coincidences: */

   if ( scale < EPSILON )
   {
      fprintf( stdout, "\n\n Warning message from function %s:", DO_ARCSCT ); 
      fprintf( stdout, "\n Coincident points A, B, C, D  !!!" );  
      fprintf( stdout, "\n [ Function returns point A.]\n " );
      arc.xx = ax;
      arc.yy = ay;

      return ZERO;
   };

   if ( d0/scale < bound2 )
   {
      arc.xx = ax;
      arc.yy = ay;

      return ZERO;
   };

   if ( d1/scale < bound2 )
   {
      arc.xx = ax;
      arc.yy = ay;

      return ZERO;
   };

   if ( d2/scale < bound2 )
   {
      arc.xx = bx;
      arc.yy = by;

      return 1.;
   };

   if ( d3/scale < bound2 )
   {
      arc.xx = bx;
      arc.yy = by;

      return 1.;
   };

   scale = sqrt( scale );
/*............................................................................*/
/* copy points: */

   x0 = ax;
   y0 = ay;
   x1 = bx;
   y1 = by;
   x2 = cx;
   y2 = cy;
   x3 = dx;
   y3 = dy;

   gamma0 = alfa;
   gamma1 = beta;

/*............................................................................*/
/* straight lines (A,B), (C,D) : */

   if (( fabs( gamma0 ) < bound0 )&&( fabs( gamma1 ) < bound0 ))   
   {
/*............................................................................*/
      ss = argz( x1-x0, y1-y0, ONE );        /* argument ( phase )            */
/*.........................................*//* of vector B-A                 */

      ss += ( PI/2. );

      xx = cos( ss );
      yy = sin( ss );

      aa00 = xx;
      aa01 = yy;
      aa02 = x0*xx + y0*yy;         /* distance, line (A,B) from origin (0,0) */

/*............................................................................*/
      ss = argz( x3-x2, y3-y2, ONE );        /* argument ( direction )        */
/*.........................................*//* of vector D-C                 */

      ss += ( PI/2. );

      xx = cos( ss );
      yy = sin( ss );

      aa10 = xx;
      aa11 = yy;
      aa12 = x2*xx + y2*yy;         /* distance, line (C,D) from origin (0,0) */

      ss = aa00 * aa11 - aa01 * aa10; /* determinant ( aa0, aa1 )             */
         
      if ( fabs( ss ) < EPSILON ) /* parallel lines (A,B) and (C,D)           */
      {
         fprintf( stdout,
            "\n\n Message from function %s:", DO_ARCSCT );
         fprintf( stdout,
            "\n Intersection points not existing or non-unique." );

         return -ONE;
      };
/* 
*//*  Hesses's normal form of line (A,B)  -->  aa00*x + aa01*y = aa02
*//*  "        "      "    "  "    (C,D)  -->  aa10*x + aa11*y = aa12
*//*  Intersection point, i.e. solution of these linear equations:
*/
      arc.xx = ( aa11 * aa02 - aa01 * aa12 )/ss;
      arc.yy = ( aa00 * aa12 - aa10 * aa02 )/ss;

      goto z_ratio; 
   };
/*............................................................................*/

   sgn = ONE;

   if (( fabs( gamma1 ) < fabs( gamma0 ))) /* interchange roles of AB and CD  */
   {
      sgn *= ( -ONE );                         

      ss = x0;
      x0 = x2;
      x2 = ss;

      ss = x1;
      x1 = x3;
      x3 = ss;

      ss = y0;
      y0 = y2;
      y2 = ss;

      ss = y1;
      y1 = y3;
      y3 = ss;

      ss = gamma0;
      gamma0 = gamma1;
      gamma1 = ss;
   }; 

/*............................................................................*/
   subarc( x2, y2, x3, y3, gamma1, ZERO );   /*                               */
/*.........................................*//* returns center (mx,my) etc.   */

   rr = arc.r;       /* radius of circle through arc CD                       */

   mx = arc.mx;      /* center M of this circle                               */
   my = arc.my;

   xx = x0 - mx;
   yy = y0 - my;

   d0 = sqrt( xx*xx + yy*yy )/rr - 1.; /* |A-M|/r - ONE */

   if( fabs( d0 ) < bound1 )  
   {
      arc.xx = x0;
      arc.yy = y0;

      goto z_ratio;
   };
      
   xx = x1 - mx;
   yy = y1 - my;

   d1 = sqrt( xx*xx + yy*yy )/rr - 1.; /* |B-M|/r - ONE */

   if ( fabs( d1 ) < bound1 )
   {
      arc.xx = x1;
      arc.yy = y1;

      goto z_ratio;
   };

   if ( EPSILON < ( d0*d1 ) )
   {
      fprintf( stdout, "\n\n Message from function %s:", DO_ARCSCT );
      fprintf( stdout, "\n Illegal constellation !" );

      if ( sgn == ONE )
      {
         fprintf( stdout, "\n [ Points A, B are not"
            " separated by circle CD.]\n " );
      }
      else if ( sgn == - ONE )
      {
         fprintf( stdout, "\n [ Points C, D are not"
            " separated by circle AB.]\n " );
      };
      return -ONE;
   };
/*............................................................................*/
/* check illegal constellations: */

   if ( fabs( gamma0 ) < bound0 )             /* straight line AB             */
   {
/*............................................................................*/
      ss = argz( x1-x0, y1-y0, ONE );        /* argument ( direction )        */
/*.........................................*//* of vector AB                  */

      ss += ( PI/2. );

      xx  = cos( ss );
      yy  = sin( ss );

      d2  =  ( x2 - x0 )*xx + ( y2 - y0 )*yy;
      d2 /= scale;

      if( fabs( d2 ) < bound1 )
      {
         arc.xx = x2;
         arc.yy = y2;

         goto z_ratio;
      };

      d3  = ( x3 - x0 )*xx + ( y3 - y0 )*yy;
      d3 /= scale;

      if( fabs( d3 ) < bound1 )
      {
         arc.xx = x3;
         arc.yy = y3;

         goto z_ratio;
      };

      if ( EPSILON < ( d2*d3 ) )
      {
         fprintf( stdout, "\n\n Message from function %s:", DO_ARCSCT );
         fprintf( stdout, "\n Illegal constellation !" );

         if ( sgn == ONE )
         {
            fprintf( stdout, "\n [ Points C, D are not "
               "separated by line AB.]\n " );
         }
         else if ( sgn == -ONE )
         {
            fprintf( stdout, "\n [ Points A, B are not "
               "separated by line CD.]\n " );
         };
         return -ONE;
      };
   }
   else /* if bound0 <= fabs( gamma0 ) */        /* true arc AB               */
   {
/*............................................................................*/
      subarc( x0, y0, x1, y1, gamma0, 0. );     /*                            */
/*............................................*/

      ss = arc.r*arc.r;

      xx = x2 - arc.mx;
      yy = y2 - arc.my;

      d2 = sqrt( xx*xx + yy*yy )/ss - 1.;

      if( fabs( d2 ) < bound1 )
      {
         arc.xx = x2;
         arc.yy = y2;

         goto z_ratio;
      };

      xx = x3 - arc.mx;
      yy = y3 - arc.my;
      d3 = sqrt( xx*xx + yy*yy )/ss - 1.;

      if( fabs( d3 ) < bound1 )
      {
         arc.xx = x3;
         arc.yy = y3;

         goto z_ratio;
      };

      if ( EPSILON < ( d2*d3 ))
      {
         fprintf( stdout, "\n\n Message from function %s:", DO_ARCSCT );
         fprintf( stdout, "\n Illegal constellation !" );

         if ( sgn == ONE )
         {
            fprintf( stdout, "\n [ Points A, B are not "
               "sepatated by circle CD.]\n " );
         }
         else if ( sgn == - ONE )
         {
            fprintf( stdout, "\n [ Points C, D are not "
               "separated by circle AB.]\n " );
         };
         return -ONE;
      };
   };

/*............................................................................*/
/* compute intersection point by Newton [ secant ] iteration:                 */

   arc.rx = x0;
   arc.ry = y0;

   s0 = 0.;                   /* d0 = sqrt((ax - mx)^2 + (ay - my)^2)/rr - 1. */
   s1 = 1.;                   /* d1 = sqrt((bx - mx)^2 + (by - my)^2)/rr - 1. */

   nn = null;
   while ( nn < MAXITER )
   {
      if ( EPSILON <= fabs( d1 - d0 ))
      {
         ss = ( s0*d1 - s1*d0 )/( d1 - d0 );  
         s0 = s1;
         s1 = ss;
/*............................................................................*/
         subarc( x0, y0, x1, y1, gamma0, ss );   /*                           */
/*.............................................*/  
         xx = arc.rx - mx;
         yy = arc.ry - my;

         d0 = d1;
         d1 = sqrt( xx*xx + yy*yy )/rr - 1.;

         if ( fabs( d1 ) < bound3 )
         {
/*
            fprintf( stdout, "\n Newton approximation stopped "
               "on iteration no. %ld.", ( nn + ONE ));
*/
            arc.xx = arc.rx;
            arc.yy = arc.ry;

            goto z_ratio;
         };
      }
      else if(( fabs( d0 - d1 ) < EPSILON )&&( fabs( d0 + d1 ) < 2.*bound3 ))
      {
         arc.xx = arc.rx;
         arc.yy = arc.ry;

         goto z_ratio;
      }
      else
      {
         fprintf( stdout, "\n\n Message from function %s:", DO_ARCSCT );
         fprintf( stdout,
            "\n Can't find intersection point of arcs AB and CD !!!" ); 
         fprintf( stdout, "\n [ Unknown reason.]\n " );
         return -ONE;
      };
      nn++;
   };

   fprintf( stderr, "\n\n Message from function %s :", DO_ARCSCT );
   fprintf( stderr, "\n Too many iterations !!!" );
   fprintf( stderr, "\n [ The number exceeds %d = macro MAXITER "
      "in %s.]\n", MAXITER, DO_ARCSCT ); 
   exit( EXIT_FAILURE );

  z_ratio:

   xx = bx - ax; 
   yy = by - ay;

   rr = xx*xx + yy*yy;    /* |B-A|^2 */

   if ( rr < EPSILON )
   {
      return ZERO;
   }
   else /* if EPSILON <= rr */
   {
/*............................................................................*/
      ss = argz( xx, yy, ONE );                    /* argument (angle ) of    */
/*...............................................*//* vector AB               */
      ss += ( PI/2. );
      mx = cos( ss );       /* unit vector orthogonal to B-A [ pos. orient. ] */
      my = sin( ss );

/* projection of vector X-A on B-A: */

      ss = arc.xx - ax;
      yy = arc.yy - ay;
      xx = ss + ( ss*mx + yy*my )*mx;  
      yy += ( ( ss*mx + yy*my )*my );

      ss = sqrt( ( xx*xx + yy*yy )/rr );       /* length ratio |S'-A| / |B-A| */

      return ss;
   };
}
/*============================================================================*/
# undef MAXITER
/*********************** end of function arcsct(*) ****************************/
