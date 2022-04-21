/* [ file: jacobi.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function jacobi(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given any complex selfadjoint matrix H = jac.h.r[][] + j*jac.h.i[][]       *
*   of order n <= JAC_MAXRNK, this subroutine returns the eigenvalues and      *
*   eigenvectors of H.                                                         *
*   Thereby, the eigenvalues are written into the diagonal jac.h.r[k][k]       *
*   of H ( k = 0,..., n-1 )  and the associated eigenvectors  into  the        *
*   columns of the matrix  E = jac.e.r[][] + j*jac.e.i[][].                    *
*   Also, the Hilbert ( spectral ) norm of H is returned in jac.hn .           *
*                                                                              *
*   [ The Jacobi algorithm is treated, e.g., in: R. Zurmuehl; 'Matrizen',      *
*     4th. ed.; Springer Verlag, Berlin 1964 ]                                 *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <string.h>
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H"
/*----------------------------------------------------------------------------*/
# ifndef JAC_STRUCT
   # define JAC_STRUCT 0
# endif

# if JAC_STRUCT == 0
   # include "jacbtp.h"
# else
/*----------------------------------------------------------------------------*/
/* Fix these bounds: */

# define JAC_MAXITR 1000  /* maximum number of iterations                     */
# define JAC_BOUND1 1.0e-277 /*                                               */
# define JAC_BOUND2 5.3e-27       /* changed from 5.3e-28:  09-07-1999        */
                                  /* [ which proved to be too restrictive     */
                                  /*              for TLM model 'mod_rj99d' ] */
# define JAC_RNDOFF 1.0e-15  /* roundoff limit for eigenvalues jac.h.r,i      */
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
      rank, i1, i2;

   unsigned char
      skew;

   struct cmatrix 
      h, e;

}  JACOBI_EV;
# endif /* JAC_STRUCT != 0 */
/*----------------------------------------------------------------------------*/
static JACOBI_EV jac = {null};
/*----------------------------------------------------------------------------*/
                                       /*---- system function prototypes ---> */
void abort( void );
int scanf( const char *format,...);

# define SLFADJ_CHECK   1   /* SLFADJ_CHECK = 1:  check selfadjointness       */
                            /* and print error message for non-selfadjoint    */
                            /* input matrix ( jac.h.r + i*jac.h.i )           */
/*============================================================================*/

JACOBI_EV *\
jacobi( JACOBI_EV *jcp )
{
/* allusions: */
/*
   extern JACOBI_EV jac;
*/
/* declarations: */

   static JACOBI_EV
     *rpt = &jac;

   static short
      ii = null,
      jj = null,
      mm = null,
      nn = null;
   
   static double
      u = ZERO,
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

   JACOBI_EV 
     *selidx( JACOBI_EV *jcp ),
     *realze( JACOBI_EV *jcp ),
     *reduce( JACOBI_EV *jcp );
/*----------------------------------------------------------------------------*/
/* initialize: */ 

   if ( jcp == NULL )
   {
      rpt = &jac;
      ( rpt->rtn ) = ONE;
      return rpt;
   };
/*............................................................................*/
   rpt = jcp;
   ( rpt->rtn ) = ONE;
   ( rpt->skew ) = null;
   nn = jcp->rank;
/*............................................................................*/
/* diag: */
     
   for ( ii=null; ii<nn; ii++ )
   {
      for ( jj=null; jj<nn; jj++ )
      {
         rpt->e.r[ii][jj] = ( double) ( ii == jj );
         rpt->e.i[ii][jj] = ZERO;
      };
   };
/*............................................................................*/
/* norm: */

   rpt->max = ZERO;

   for ( ii=null; ii<nn; ii++ )  
   {
      rpt->sn[ii] = ZERO;
      for ( jj=null; jj<nn; jj++ )
      {
         if ( ii != jj )
         {
            u = rpt->h.r[ii][jj];
            v = rpt->h.i[ii][jj];
            rpt->sn[ii] += ( u*u + v*v );
         };
      };
      u  = rpt->h.r[ii][ii];
      v  = rpt->sn[ii] + u*u;

      if ( rpt->max < v )  
         rpt->max = v;  
   };

   if ( JAC_BOUND1 <= ( rpt->max ))
   {
/*............................................................................*/
/* bounds: */

      ( rpt->bnd1 ) = ( rpt->max )*JAC_BOUND2;
      ( rpt->bnd2 ) = sqrt( rpt->bnd1 );

/*............................................................................*/
/* iterate: */ 

      mm = null; do
      {
         rpt = selidx( rpt );
       
         if ( rpt == NULL ) 
         {
            printf( "\n\n Error message from function jacobi(*):"
               "\n unsatisfied function selidx(*)\n " );
            return NULL;
         };

         if (( rpt->bnd1 ) <= ( rpt->max ))
         {
/*............................................................................*/
            rpt = realze( rpt );      /*                                      */
/*..................................*/
            if ( rpt == NULL )
            {
               printf( "\n\n Error message from function jacobi(*):"
                  "\n unsatisfied function realze(*)\n " );
               return NULL;
            };
/*............................................................................*/
            rpt = reduce( rpt );      /*                                      */
/*..................................*/
            if ( rpt == NULL )
            {
               printf( "\n\n Error message from function jacobi(*):"
                  "\n unsatisfied function reduce(*)\n " );
               return NULL;
            };
         } /* end if (( rpt->bnd1 ) <= ( rpt->max )) */
         else
            goto hilbert_norm;

         mm++;

      } while ( mm < JAC_MAXITR );
/*............................................................................*/
/* excess: */

      printf( " \n\n Message from subroutine %s :", __func__ );
      printf( " \n Too many iterations !!!" );
      printf( " \n [ Jacobi algorithm interrupted:"
         " returning to calling program.] " );   
      return rpt;

   }; /* end if ( JAC_BOUND1 <= ( rpt->max )) */ 
/*............................................................................*/
/* compute Hilbert [ spectral ] norm: */

  hilbert_norm:
   
   rpt->hn = ZERO;

   for ( ii=null; ii<nn; ii++ )    
   {
      u = rpt->h.r[ii][ii];
      u *= u;
      if (( rpt->hn ) < u ) 
         ( rpt->hn ) = u; 
   };
   ( rpt->hn ) = sqrt( rpt->hn );
/*............................................................................*/
# ifdef JAC_RNDOFF

   u = JAC_RNDOFF*( rpt->hn );

   ii = null;
   while ( ii < nn )
   {
      if ( fabs( rpt->h.r[ii][ii] ) < u )
         rpt->h.r[ii][ii] = ZERO;
      ii++;
   };

# endif
/*............................................................................*/

   return rpt;
}
/*============================================================================*/

JACOBI_EV *\
selidx( JACOBI_EV *jcp )
{ 
/* allusions: */
/*
   extern JACOBI_EV jac;
*/
/* declarations: */

   static JACOBI_EV
     *rpt;

   static double
      u = ZERO,
      v = ZERO,
      max = ZERO;
     
   static short
     ii = null,
     i1 = null,
     i2 = null,
     nn = null;
/*----------------------------------------------------------------------------*/
   rpt = jcp;
   nn = jcp->rank;
   i1 = jcp->i1;
   i2 = jcp->i2;
   ( rpt->rtn ) = ONE;
/*............................................................................*/
   i2 = null;
   max = rpt->sn[null];

   for ( ii=ONE; ii<nn; ii++ )
   {
      if ( max < ( rpt->sn[ii] ))
      {
         i2 = ii;
         max = rpt->sn[ii];
      };
   };

   max = ZERO;
   i1 = null;

   for ( ii=null; ii<nn; ii++ )
   {
      if ( ii != i2 )
      {
         u = rpt->h.r[ii][i2];
         v = rpt->h.i[ii][i2];
         u = ( u*u + v*v ) ;

         if ( max < u )
         {
            i1 = ii;
            max = u;
         };
      };
   };

   ( rpt->max ) = max;
   ( rpt->i1 ) = i1;
   ( rpt->i2 ) = i2;

   return rpt;
}
/*============================================================================*/

JACOBI_EV *\
realze( JACOBI_EV *jcp )
{ 
/* allusions: */
/*
   extern JACOBI_EV jac;
*/
/* declarations: */

   static JACOBI_EV 
     *rpt;

   static double 
      re = ZERO,
      im = ZERO,
      x1 = ZERO,
      x2 = ZERO,
      y1 = ZERO,
      y2 = ZERO; 

   static short
      ii = null,
      i1 = null,
      i2 = null,
      nn = null;

   static signed char
      sgr = null,
      sgi = null;
/*...........................................................................*/
   rpt = jcp;
   nn = jcp->rank;
   i1 = jcp->i1;
   i2 = jcp->i2;
   ( rpt->rtn ) = ONE;
/*............................................................................*/
/* phase: */ 

   re = rpt->h.r[i1][i2];
   im = rpt->h.i[i1][i2];

   ( rpt->mod ) = sqrt( re*re + im*im );

   if (( rpt->mod ) < JAC_BOUND1 )
   {
      ( rpt->phi ) = ZERO;

      return rpt;
   };
/*............................................................................*/
   if ( ZERO < re )
      sgr = 1;
   else if ( re < ZERO )
      sgr = -1;

   if ( ZERO < im )
      sgi = 1;
   else if ( im < ZERO )
      sgi = -1;
     
   if (( fabs( re )) < JAC_BOUND1 )
      ( rpt->phi ) = .5*sgi*PI;
   else /* if ( JAC_BOUND1 <= ( fabs( re ))) */
   {
      if (( fabs( im )) < JAC_BOUND1 )
         ( rpt->phi ) = .5*( 1.- sgr )*PI;
      else
         ( rpt->phi ) = atan2( im, re ) + .5*( 1.- sgr )*sgi*PI;

      if ((( fabs( rpt->phi )) < JAC_BOUND1 )||
          (( fabs(( rpt->phi )) - PI ) < JAC_BOUND1 ))
      return rpt;

   }; /* end if ( JAC_BOUND1 <= fabs( re )) */
/*............................................................................*/
/* [ phi = arg(H); H = mod*exp(j*phi) ] */

   ( rpt->c ) = cos(.5*( rpt->phi ));
   ( rpt->s ) = sin(.5*( rpt->phi ));

   for ( ii=null; ii<nn; ii++)
   {
      x1 = rpt->h.r[ii][i1];
      y1 = rpt->h.i[ii][i1];
      x2 = rpt->h.r[ii][i2];
      y2 = rpt->h.i[ii][i2];

      rpt->h.r[ii][i1] =  x2*( rpt->c ) + y2*( rpt->s );
      rpt->h.i[ii][i1] = -x2*( rpt->s ) + y2*( rpt->c ); 
      rpt->h.r[ii][i2] =  x1*( rpt->c ) - y1*( rpt->s );
      rpt->h.i[ii][i2] =  x1*( rpt->s ) + y1*( rpt->c ); 

      x1 = rpt->e.r[ii][i1];
      y1 = rpt->e.i[ii][i1];
      x2 = rpt->e.r[ii][i2];
      y2 = rpt->e.i[ii][i2];

      rpt->e.r[ii][i1] =  x2*( rpt->c ) + y2*( rpt->s );
      rpt->e.i[ii][i1] = -x2*( rpt->s ) + y2*( rpt->c );
      rpt->e.r[ii][i2] =  x1*( rpt->c ) + y1*( rpt->s );
      rpt->e.i[ii][i2] =  x1*( rpt->s ) + y1*( rpt->c );
   }; 

   for ( ii=null; ii<nn; ii++)
   {
      x1 = rpt->h.r[i1][ii];
      y1 = rpt->h.i[i1][ii];
      x2 = rpt->h.r[i2][ii];
      y2 = rpt->h.i[i2][ii];

      rpt->h.r[i1][ii] = x2*( rpt->c ) - y2*( rpt->s );
      rpt->h.i[i1][ii] = y2*( rpt->c ) + x2*( rpt->s );
      rpt->h.r[i2][ii] = x1*( rpt->c ) + y1*( rpt->s );
      rpt->h.i[i2][ii] = y1*( rpt->c ) - x1*( rpt->s );
   };
/*............................................................................*/
/* the following should yet be ZERO by algorithm [and is strengthened here]: */

   ( rpt->h.i[i1][i2] ) = ZERO;
   ( rpt->h.i[i2][i1] ) = ZERO;

   return rpt; 
}
/*============================================================================*/

JACOBI_EV * \
reduce( JACOBI_EV *jcp )
{
/* allusions: */
/*
   extern JACOBI_EV jac;
*/
/* declarations: */

   static JACOBI_EV
     *rpt = &jac;

   static double
      u = ZERO,
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

   static short
      ii = null,
      i1 = null,
      i2 = null,
      nn = null;
          
   static char
      ptr[STS_SIZE] = {null};

   static signed char
      sgn = null;
/*----------------------------------------------------------------------------*/
/* allocate: */
/*
   ptr = ( char *) calloc( STS_SIZE, ONE );
*/
/*............................................................................*/
   rpt = jcp;
   nn = jcp->rank;
   i1 = jcp->i1;
   i2 = jcp->i2;
   ( rpt->rtn ) = ONE;
/*............................................................................*/
/* reduce: */
       
   x11 = rpt->h.r[i1][i1];
   x12 = rpt->h.r[i1][i2];
   x21 = rpt->h.r[i2][i1];
   x22 = rpt->h.r[i2][i2];
   y11 = rpt->h.i[i1][i1];
   y22 = rpt->h.i[i2][i2];

# if SLFADJ_CHECK == 1

   if ( rpt->bnd2 < fabs(x12-x21) )  
      goto skew;
   if ( rpt->bnd2 < fabs(y11) )     
      goto skew;
   if ( rpt->bnd2 < fabs(y22) )     
      goto skew;

# endif

   if (( fabs( x11-x22 )) < JAC_BOUND1 )
   {
      ( rpt->c ) = 1./ sqrt(2.);

      if  ( x12 == ZERO ) 
         ( rpt->s ) = rpt->c ;
      else  
         ( rpt->s ) = ( x12 / fabs( x12 )) * ( rpt->c );
   }
   else 
   {
      t = 2.*x12 / ( x11 - x22 );
      if ( JAC_BOUND1 < fabs( t ))
      {
         if ( ZERO < t ) sgn = 1;
         if ( t < ZERO ) sgn = -1;

         u = 1./ sqrt( 1.+ t*t );
         ( rpt->c ) = sqrt(.5*( 1.+ u ));
         ( rpt->s ) = sgn * sqrt(.5*( 1.- u ));
      }
      else
      {
         ( rpt->s ) = .5 * t;
         ( rpt->c ) = 1. - .5*( rpt->s )*( rpt->s );
      };
   };

   for ( ii=null; ii<nn; ii++ )
   {
      x1 = rpt->h.r[ii][i1];
      y1 = rpt->h.i[ii][i1];
      x2 = rpt->h.r[ii][i2];
      y2 = rpt->h.i[ii][i2];

      rpt->h.r[ii][i1] = x1*( rpt->c ) + x2*( rpt->s );
      rpt->h.i[ii][i1] = y1*( rpt->c ) + y2*( rpt->s );
      rpt->h.r[ii][i2] = x1*( rpt->s ) - x2*( rpt->c );
      rpt->h.i[ii][i2] = y1*( rpt->s ) - y2*( rpt->c );

      x1 = rpt->e.r[ii][i1];
      y1 = rpt->e.i[ii][i1];
      x2 = rpt->e.r[ii][i2];
      y2 = rpt->e.i[ii][i2];

      rpt->e.r[ii][i1] = x1*( rpt->c ) + x2*( rpt->s );
      rpt->e.i[ii][i1] = y1*( rpt->c ) + y2*( rpt->s );
      rpt->e.r[ii][i2] = x1*( rpt->s ) - x2*( rpt->c );
      rpt->e.i[ii][i2] = y1*( rpt->s ) - y2*( rpt->c );
   };

   for ( ii=null; ii<nn; ii++ )
   {
      x1 = rpt->h.r[i1][ii];
      y1 = rpt->h.i[i1][ii];
      x2 = rpt->h.r[i2][ii];
      y2 = rpt->h.i[i2][ii];

      rpt->h.r[i1][ii] = x1*( rpt->c ) + x2*( rpt->s );
      rpt->h.i[i1][ii] = y1*( rpt->c ) + y2*( rpt->s );
      rpt->h.r[i2][ii] = x1*( rpt->s ) - x2*( rpt->c );
      rpt->h.i[i2][ii] = y1*( rpt->s ) - y2*( rpt->c );
   };
/*............................................................................*/
/* the following should yet be ZERO by algorithm [ and is strenghened here ]: */

   rpt->h.r[i1][i2] = ZERO;
   rpt->h.r[i2][i1] = ZERO;
/*............................................................................*/
/* Stephans_norm [ here redefined ]: */  

   rpt->sn[i1] = ZERO;  
   rpt->sn[i2] = ZERO;

   for ( ii=null; ii<nn; ii++ )
   {
      if ( ii != i1 )
      {
         u = rpt->h.r[ii][i1];
         v = rpt->h.i[ii][i1];
         ( rpt->sn[i1] ) += ( u*u + v*v );
      };

      if ( ii != i2 )
      {
         u = rpt->h.r[ii][i2];
         v = rpt->h.i[ii][i2];
         ( rpt->sn[i2] ) += ( u*u + v*v );
      };
   };
     
   return rpt;
/*............................................................................*/
# if SLFADJ_CHECK == 1
  skew: 
       
   if ( rpt->skew == ONE )
      return rpt;

   printf( "\n\n WARNING from subroutine %s :", __func__ );
   printf( "\n Transferred matrix H = jac.r,i[][] non-selfadjoint !!!" );
   printf( "\n Continue Jacobi algorithm ?   [ y/n ] ....................: " );
   scanf( "%s", ptr );

   if ( *ptr == 'y' )
   {
      rpt->skew = ONE;
      printf( "\n [ Warning ignored: cont'd Jacobi algorithm. ] \n" );
   }
   else
   {
      printf( "\n Jacobi algorithm interrupted"
	      " [ returning to calling program]. \n" );            
      return rpt; 
   };
  
   return rpt;

# endif
/*............................................................................*/

   return rpt;
} 
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif

# undef JAC_MAXRNK
# undef JAC_MAXITR
# undef JAC_BOUND1
# undef JAC_BOUND2
# undef JAC_RNDOFF
/********************* end of Jacobi subroutine jacobi(*) *********************/
