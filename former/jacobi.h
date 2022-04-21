/* [ file: jacobi.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function jacobi(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given any complex selfadjoint matrix H = jac.h.r[][] + j*jac.h.i[][]       *
*   of order n <=JAC_MAXRNK, this subroutine returns the eigenvalues and       *
*   eigenvectors of H.                                                         *
*   The eigenvalues are thereby written into the diagonal jac.h.r[k][k]        *
*   of H ( k = 0,..., n-1 )  and the associated eigenvectors  into  the        *
*   columns of the matrix  E = jac.e.r[][] + j*jac.e.i[][].                    *
*   Also, the Hilbert ( spectral ) norm of H is returned in jac.hn .           *
*                                                                              *
*   [ The Jacobi algorithm is treated, e.g., in: R. Zurmuehl; 'Matrizen',      *
*     4th. ed.; Springer Verlag, Berlin 1964 ]                                 *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*----------------------------------------------------------------------------*/
/*
# include <stdio.h>
*/
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifndef JAC_STRUCT
   # define JAC_STRUCT 1
# endif

# if JAC_STRUCT == 0
   # include "jacbtp.h"
# else
/*----------------------------------------------------------------------------*/
/* Fix these bounds: */

# define JAC_MAXITR 1000  /* maximum number of iterations                     */
# define JAC_BOUND1 1.0e-277 /*                                               */
# define JAC_BOUND2 1.3e-27  /* changed from 5.3e-28:  09-07-1999             */
                          /* [ which proved to be too restrictive             */
                          /*                      for DSC model 'mod_rj99d' ] */

# define JAC_RNDOFF 1.0e-16  /* roundoff limit for eigenvalues jac.h.r,i      */
/*----------------------------------------------------------------------------*/
                                                  /*---- structure types ---> */
/*----------------------------------------------------------------------------*/
# define JAC_MAXRNK 10    /* maximum order of matrix H = jac.h.r + j*jac.h.i  */
/*----------------------------------------------------------------------------*/
struct cmatrix  
{
   double 
      r[JAC_MAXRNK][JAC_MAXRNK], i[JAC_MAXRNK][JAC_MAXRNK];
};
/*----------------------------------------------------------------------------*/
typedef struct
{
   signed char
      rtn;

   double
      max, bnd1, bnd2, mod, phi, hn, c, s,
      sn[JAC_MAXRNK];

   short
      n, i1, i2;

   unsigned char skew;

   struct cmatrix 
      h, e;

}  JACOBI_EV;
static JACOBI_EV jac = {null};
# endif /* JAC_STRUCT != 0 */
/*----------------------------------------------------------------------------*/
                                       /*---- system function prototypes ---> */
void abort( void );
int scanf( const char *format,...);

# define SLFADJ_CHECK   1   /* SLFADJ_CHECK = 1:  check selfadjointness       */
                            /* and print error message for non-selfadjoint    */
                            /* input matrix ( jac.h.r + i*jac.h.i )           */
/*============================================================================*/

JACOBI_EV *\
jacobi( short n )
{
/* allusions: */
/*
   extern JACOBI_EV jac;
*/
/* declarations: */

   static JACOBI_EV *jev = &jac;

   static short  ii = null,
                 jj = null,
                 mm = null;
   
   static signed char ind = null;

   static double  u = ZERO,
                  v = ZERO;

/* math. function prototypes: */

   double  fabs( double x );
   double  sqrt( double x );
   double   cos( double x );
   double   sin( double x );
   double  atan( double x );
   double atan2( double x, double y );
   int abs( int h );

/* user defined function prototypes: */

   short selidx( void );
   short realze( void );
   short reduce( void );
/*............................................................................*/
   jac.n     = n;
   jac.skew  = null;
/*............................................................................*/
/* diag: */
     
   for ( ii=null ; ii<jac.n ; ii++ )
   {
      for ( jj=null ; jj<jac.n ; jj++ )
      {
         jac.e.r[ii][jj] = ( double) ( ii == jj ); 
         jac.e.i[ii][jj] = ZERO;
      };
   };
/*............................................................................*/
/* norm: */

   jac.max = ZERO;

   for ( ii=null ; ii<jac.n ; ii++ )  
   {
      jac.sn[ii] = ZERO;
      for ( jj=null ; jj<jac.n ; jj++ )
      {
         if ( ii != jj )
         {
            u = jac.h.r[ii][jj];
            v = jac.h.i[ii][jj];
            jac.sn[ii] += ( u*u + v*v );
         };
      };
      u  = jac.h.r[ii][ii];
      v  = jac.sn[ii]  +  u*u;
      if ( jac.max < v )  
         jac.max = v;  
   };
      
   if ( JAC_BOUND1 <= jac.max ) 
   {
/*----------------------------------------------------------------------------*/
/* bounds: */

      jac.bnd1 = jac.max*JAC_BOUND2;
      jac.bnd2 = sqrt(jac.bnd1);

/* iterate: */ 

      mm = null; do
      {
         ind = selidx( );
       
         if ( ind == null ) 
         {
            printf( "\n\n Error on function 'jacobi(*)'"
		    "[ unsatisfied 'selidx(*)' ] " ); 
            return null;
         };

         if ( jac.bnd1 <= jac.max )
         {
            ind = realze( );

            if ( ind == null )
            {
               printf( "\n\n Error on function 'jacobi(*)'"
	          "[ unsatisfied 'realze(*)' ]" );
               return null;
            };

            ind = reduce( );

            if ( ind == null )
            {
               printf( "\n\n Error on function 'jacobi(*)'"
		  "[ unsatisfied 'reduce(*)' ] " );
               return null;
            };
         } /* end if ( jac.bnd1 <= jac.max ) */
         else
            goto hilbert_norm;

         mm++;

      } while ( mm < JAC_MAXITR ); 
/*............................................................................*/
/* excess: */

      printf( " \n\n Message from subroutine  %s :", __func__ );
      printf( " \n Too many iterations !!!" );
      printf( " \n [ Jacobi algorithm interrupted: "
         "returning to calling program.] " );   
      return null;
   }; /* end if ( JAC_BOUND1 <= jac.max ) */
/*............................................................................*/
/* compute Hilbert ( spectral ) norm: */

  hilbert_norm:

   jac.hn = ZERO;

   for ( ii=null ; ii<jac.n ; ii++ )    
   {
      u = jac.h.r[ii][ii];
      u *= u;
      if ( jac.hn < u ) 
         jac.hn = u; 
   };
   jac.hn = sqrt( jac.hn );

# ifdef JAC_RNDOFF

   u = JAC_RNDOFF*jac.hn;

   ii=null;
   do
   {
      if ( fabs( jac.h.r[ii][ii] ) < u )
         jac.h.r[ii][ii] = ZERO;
      ii++;
   }  while ( ii<jac.n );

# endif

   return jev;
}
/*============================================================================*/

short selidx( void )
{ 
/* allusions: */
/*
   extern JACOBI_EV jac;
*/
/* declarations: */

   static double u = ZERO,
                 v = ZERO;
     
   static short ii = null;
                       
   jac.i2 = null;
   jac.max  = jac.sn[null];

   for ( ii=ONE ; ii<jac.n ; ii++ )
   {
      if ( jac.max < jac.sn[ii] )
      {
         jac.i2 = ii;
         jac.max  = jac.sn[ii];
      };
   };

   jac.max = ZERO;
   jac.i1  = null;

   for ( ii=null ; ii<jac.n ; ii++ )
   {
      if ( ii != jac.i2 ) 
      {
         u = jac.h.r[ii][jac.i2];
         v = jac.h.i[ii][jac.i2];
         u = u*u + v*v ;
         if ( jac.max < u )
         {
            jac.i1   = ii;
            jac.max  = u;
         };
      }; 
   };

   return ONE;
}
/*============================================================================*/

short realze( void )
{ 
/* allusions: */
/*
   extern JACOBI_EV jac;
*/
/* declarations: */

   static double re = ZERO,
                 im = ZERO,
                 x1 = ZERO,
                 x2 = ZERO,
                 y1 = ZERO,
                 y2 = ZERO; 

   static short ii = null;

   static signed char sgr = null,
                      sgi = null;
/*............................................................................*/
/* phase: */ 

   re = jac.h.r[jac.i1][jac.i2];
   im = jac.h.i[jac.i1][jac.i2];
   jac.mod = sqrt( re*re + im*im );

   if ( jac.mod < JAC_BOUND1 )
   {
      jac.phi = ZERO;

      return ONE;
   };
           
   if      ( re > ZERO ) { sgr =  1; }
   else if ( re < ZERO ) { sgr = -1; };

   if      ( im > ZERO ) { sgi =  1; }
   else if ( im < ZERO ) { sgi = -1; };
     
   if ( fabs(re) < JAC_BOUND1 )
   {
      jac.phi = .5*sgi*PI;
      goto mltply;
   };

   if ( fabs(im) < JAC_BOUND1 )
   {
      jac.phi = .5*( 1.- sgr )*PI;
   }
   else
   {
      jac.phi = atan2( im, re ) + .5*( 1.- sgr )*sgi*PI;
   };

   if (( fabs( jac.phi ) < JAC_BOUND1 )||( fabs( jac.phi - PI ) < JAC_BOUND1 ))
                                         /* arg(H) = phi; H = mod*exp(j*phi) */
   {
      return ONE;
   };

  mltply:

   jac.c  = cos(.5*jac.phi);
   jac.s  = sin(.5*jac.phi);

   for ( ii=null; ii<jac.n; ii++)
   {
      x1 = jac.h.r[ii][jac.i1];
      y1 = jac.h.i[ii][jac.i1];
      x2 = jac.h.r[ii][jac.i2];
      y2 = jac.h.i[ii][jac.i2];

      jac.h.r[ii][jac.i1] =  x2*jac.c + y2*jac.s;
      jac.h.i[ii][jac.i1] = -x2*jac.s + y2*jac.c; 
      jac.h.r[ii][jac.i2] =  x1*jac.c - y1*jac.s;
      jac.h.i[ii][jac.i2] =  x1*jac.s + y1*jac.c; 

      x1 = jac.e.r[ii][jac.i1];
      y1 = jac.e.i[ii][jac.i1];
      x2 = jac.e.r[ii][jac.i2];
      y2 = jac.e.i[ii][jac.i2];

      jac.e.r[ii][jac.i1] =  x2*jac.c + y2*jac.s;
      jac.e.i[ii][jac.i1] = -x2*jac.s + y2*jac.c;
      jac.e.r[ii][jac.i2] =  x1*jac.c + y1*jac.s;
      jac.e.i[ii][jac.i2] =  x1*jac.s + y1*jac.c;
   }; 

   for ( ii=null; ii<jac.n; ii++)
   {
      x1 = jac.h.r[jac.i1][ii];
      y1 = jac.h.i[jac.i1][ii];
      x2 = jac.h.r[jac.i2][ii];
      y2 = jac.h.i[jac.i2][ii];

      jac.h.r[jac.i1][ii] = x2*jac.c - y2*jac.s;
      jac.h.i[jac.i1][ii] = y2*jac.c + x2*jac.s;
      jac.h.r[jac.i2][ii] = x1*jac.c + y1*jac.s;
      jac.h.i[jac.i2][ii] = y1*jac.c - x1*jac.s;
   };
   jac.h.i[jac.i1][jac.i2] = ZERO;  /*  already ZERO by algorithm  */
   jac.h.i[jac.i2][jac.i1] = ZERO;  /*  --------- dito ----------- */

   return ONE; 
}
/*============================================================================*/

short reduce( void )
{
/* allusions: */
/*
   extern JACOBI_EV jac;
*/
/* declarations: */

   static double u = ZERO,
                 v = ZERO,
                 t = ZERO, 
                x1 = ZERO,
                x2 = ZERO,
                y1 = ZERO,
                y2 = ZERO, 
               x11 = ZERO,
               x12 = ZERO,
               x21 = ZERO,
               x22 = ZERO,
               y11 = ZERO,
               y22 = ZERO;

   static short ii = null;
          
   static char ptr[STS_SIZE] = {null};

   static signed char sgn = null;

/* allocate: */
/*
   ptr = ( char *) calloc( STS_SIZE, ONE );
*/
/*............................................................................*/
/* reduce: */
       
   x11 = jac.h.r[jac.i1][jac.i1];
   x12 = jac.h.r[jac.i1][jac.i2];
   x21 = jac.h.r[jac.i2][jac.i1];
   x22 = jac.h.r[jac.i2][jac.i2];
   y11 = jac.h.i[jac.i1][jac.i1];
   y22 = jac.h.i[jac.i2][jac.i2];

# if SLFADJ_CHECK == 1

   if ( jac.bnd2 < fabs(x12-x21) )  
      goto skew;
   if ( jac.bnd2 < fabs(y11) )     
      goto skew;
   if ( jac.bnd2 < fabs(y22) )     
      goto skew;

# endif

   if ( fabs(x11-x22) < JAC_BOUND1 )
   {
      jac.c = 1. / sqrt(2.);
      if  ( x12 == ZERO ) 
      {
         jac.s = jac.c ;
      }
      else  
      {
         jac.s = ( x12 / fabs(x12) ) * jac.c;
      };
   }
   else 
   {
      t = 2.*x12 / (x11 - x22);
      if ( JAC_BOUND1 < fabs(t) ) 
      {
         if ( ZERO < t )  sgn = 1;
         if ( t < ZERO ) sgn = -1; 
         u = 1. / sqrt( 1. + t*t );
         jac.c = sqrt( .5*(1.+ u));
         jac.s = sgn * sqrt(.5*(1.- u));
      }
      else
      {
         jac.s = .5 * t;
         jac.c = 1. - jac.s*jac.s/2; 
      };
   };

   for ( ii=null ; ii<jac.n ; ii++ )
   {
      x1 = jac.h.r[ii][jac.i1];
      y1 = jac.h.i[ii][jac.i1];
      x2 = jac.h.r[ii][jac.i2];
      y2 = jac.h.i[ii][jac.i2];

      jac.h.r[ii][jac.i1] = x1*jac.c + x2*jac.s;
      jac.h.i[ii][jac.i1] = y1*jac.c + y2*jac.s;
      jac.h.r[ii][jac.i2] = x1*jac.s - x2*jac.c;
      jac.h.i[ii][jac.i2] = y1*jac.s - y2*jac.c;

      x1 = jac.e.r[ii][jac.i1];
      y1 = jac.e.i[ii][jac.i1];
      x2 = jac.e.r[ii][jac.i2];
      y2 = jac.e.i[ii][jac.i2];

      jac.e.r[ii][jac.i1] = x1*jac.c + x2*jac.s;
      jac.e.i[ii][jac.i1] = y1*jac.c + y2*jac.s;
      jac.e.r[ii][jac.i2] = x1*jac.s - x2*jac.c;
      jac.e.i[ii][jac.i2] = y1*jac.s - y2*jac.c;
   };

   for ( ii=null ; ii<jac.n ; ii++ )
   {
      x1 = jac.h.r[jac.i1][ii];
      y1 = jac.h.i[jac.i1][ii];
      x2 = jac.h.r[jac.i2][ii];
      y2 = jac.h.i[jac.i2][ii];

      jac.h.r[jac.i1][ii] = x1*jac.c + x2*jac.s;
      jac.h.i[jac.i1][ii] = y1*jac.c + y2*jac.s;
      jac.h.r[jac.i2][ii] = x1*jac.s - x2*jac.c;
      jac.h.i[jac.i2][ii] = y1*jac.s - y2*jac.c;
   };

   jac.h.r[jac.i1][jac.i2] = ZERO; /* ZERO by algorithm [ and here imposed ]  */
   jac.h.r[jac.i2][jac.i1] = ZERO; /* ----------------- "  -------------------*/

/* stephans_norm: */               /* redefine Stephan's norm                 */

   jac.sn[jac.i1] = ZERO;  
   jac.sn[jac.i2] = ZERO;

   for ( ii=null; ii<jac.n; ii++ )
   {
      if ( ii != jac.i1 ) 
      {
         u = jac.h.r[ii][jac.i1];
         v = jac.h.i[ii][jac.i1];
         jac.sn[jac.i1] += ( u*u + v*v );
      };

      if ( ii != jac.i2 ) 
      {
         u = jac.h.r[ii][jac.i2];
         v = jac.h.i[ii][jac.i2];
         jac.sn[jac.i2] += ( u*u + v*v );
      };
   };
     
   return ONE;

/*............................................................................*/
# if SLFADJ_CHECK == 1
  skew: 
       
   if ( jac.skew == ONE ) 
      return ONE;

   printf( "\n\n WARNING from subroutine %s :", __func__ );
   printf( "\n Transferred matrix H = jac.h.r,i[][] non-selfadjoint !!!" );
   printf( "\n Continue Jacobi algorithm ?   [ y/n ] ....................: " );
   scanf( "%s", ptr );

   if ( *ptr == 'y' )
   {
      jac.skew = ONE;
      printf( "\n [ Warning ignored: cont'd Jacobi algorithm. ] \n" );
   }
   else
   {
      printf( "\n Jacobi algorithm interrupted"
	      " [ returning to calling program]. \n" );            
      return null; 
   };
  
   return ONE;
# endif

} 
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/*----------------------------------------------------------------------------*/
# undef JAC_MAXRNK
# undef JAC_MAXITR
# undef JAC_BOUND1
# undef JAC_BOUND2
# undef JAC_RNDOFF
/********************* end of Jacobi subroutine jacobi(*) *********************/
