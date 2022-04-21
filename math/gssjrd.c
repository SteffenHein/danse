/* [ file: gssjrd.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function gssjrd(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Gauss-Jordan elimination                                                   *
*                                                                              *
*   Given any non-singular complex square matrix M of rank n > 0 and           *
*   an array of n complex numbers a[i] ( i=0,.., n-1 ), this function          *
*   returns the solution z[k] of the linear equation                           *
*                                                                              *
*                sum( k=0,...,n-1)  M[i][k] * z[k]  =  a[i] .                  *
*                                                                              *
*   To this end, M and a[i] ( i=0,...,n-1 ) are transferred to the function    *
*   as the columns k=0,...,n-1 and k=n, respectively, of the structure         *
*   matrix                                                                     *
*                                                                              *
*          ( gjp->mr[i][k] + j*( gjp->mi[i][k] ))  ( i=0,...,n-1 ).            *
*                                                                              *
*   The function the writes the solution z[k] into                             *
*                                                                              *
*          ( gjp->zr[k][0] + j*( gjp->zi[k][0] ))  ( k=0,...,n-1 )             *
*                                                                              *
*   ( option opt = 'e' ).                                                      *
*                                                                              *
*   Also, in option opt = 'i', the inverse matrix of M,  M^(-1),               *
*   is computed and written into                                               *
*                                                                              *
*          ( gjp->zr[i][k] + j*( gjp->zi[i][k] ))  ( i,k=0,...,n-1 )           *
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
/* operational constants:                                                     */

# ifndef GSS_STRUCT
   # define GSS_STRUCT 0   /* 1: define struct GAUSS_JRD [may be done elsw.] */
# endif

# if GSS_STRUCT == 0
   # include "gssjtp.h"
# else
/*----------------------------------------------------------------------------*/
/* The structure type definition of Gauss-Jordan algorithm functions          */
/* gssjrd(*) [line pivoting] and gssjpv(*) [full pivoting ]                   */
/*----------------------------------------------------------------------------*/
# ifndef GSS_DISP
   # define GSS_DISP  1     /* 1: display if coeff. mtx.gss.m[][] is singular */
# endif
# ifndef GSS_IMPLCT
   # define GSS_IMPLCT 1    /* 1: implicit pivoting  [ normalizing equations ]*/
# endif
# ifndef GSS_STRUCT
   # define GSS_STRUCT 1    /* 1: define struct gaussc gss [may be done elsw.]*/
# endif
# ifndef GSS_LONG_DBL
   # define GSS_LONG_DBL 0  /* 1: long double prexcision */
# endif
# ifndef GSS_MAXRNK
   # define GSS_MAXRNK 10   /* maximum rank of matrix M */
# endif
# ifndef GSS_SINGLR 
   # define GSS_SINGLR (1.e-277) /* bound for singularity check */
# endif
/*----------------------------------------------------------------------------*/
typedef struct
{
   signed char
      rtn;

   char 
      opt;

   short
      rank, /* rank of matrix M = ( mr[][]+j*mi[][] ) */
      neqs; /* number of linear equations to be solved [ simultaneously ] */
            /* - only used in function gssjpv(*) */

   # if GSS_LONG_DBL  == 1
      long double
   # else
      double
   # endif

         dtr,                         /* determinant  ( real  part )          */
         dti,                         /* "            ( imag. part )          */
         mr[GSS_MAXRNK][2*GSS_MAXRNK],/* r = real part of matrix              */
         mi[GSS_MAXRNK][2*GSS_MAXRNK],/* i = imag."    "  "                   */
         zr[GSS_MAXRNK][GSS_MAXRNK],
         zi[GSS_MAXRNK][GSS_MAXRNK];

}  GAUSS_JRD;
# endif /* GSS_STRUCT != 0 */
/*----------------------------------------------------------------------------*/
static GAUSS_JRD gss = {null};
/*============================================================================*/

GAUSS_JRD *\
gssjrd( GAUSS_JRD *gjp )
{
/* allusions: */
/*
   extern GAUSS_JRD gss;
*/
/* declarations: */

   static GAUSS_JRD
     *rpt = &gss;

# if GSS_LONG_DBL  == 1
   static long double
# else
   static double
# endif
      uu = ZERO,
      vv = ZERO,
      rr = ZERO,
      ss = ZERO, 
      u1 = ZERO,
      v1 = ZERO,
      u2 = ZERO,
      v2 = ZERO;

   static short 
      hh = null,
      ii = null,
      jj = null,
      kk = null,
      ll = null,
      mm = null,
      nn = null,
      h1 = null,
      n1 = null;
/*----------------------------------------------------------------------------*/
/* options:                                                                   */
/*                     compute
   opt = 'i'[inversion]   -->   matrix inversion
   opt = 'e'[equation]    -->   linear equation  
   opt = 'd'[determinant] -->   determinant
*/
/*...........................................................................*/
/* initialize: */

   if ( gjp == NULL )
   {
      rpt = &gss;
      ( rpt->rtn ) = ONE;
      return rpt;
   };
/*...........................................................................*/
   rpt = gjp;
   nn = ( gjp->rank );
/*...........................................................................*/

   if ( GSS_MAXRNK < nn )
   {
      printf( "\n\n Error message from funcion %s:", __func__ );
      printf( "\n\n Rank of matrix 'rpt->m[][]' too high !!!" );
      printf( "\n [ The maximum rank is %d = macro GSS_MAXRNK", GSS_MAXRNK );
      printf( "\n   - change macro in compliance "
         "with memory resources. ]\n " );
      return NULL;
   };
/*............................................................................*/
/* option: matrix inversion: */

   ( rpt->rtn ) = ONE;
   ( rpt->dtr ) = ONE;
   ( rpt->dti ) = ZERO;

   if ((( gjp->opt ) == 'i' )
     ||(( gjp->opt ) == 'I' ))
   {
      n1 = TWO*( nn );
      for ( hh=null; hh<nn; hh++ ) /* write unit matrix */
      {
         for ( ii=nn; ii<n1; ii++ )
         {   
            ll = ii-nn;

# if GSS_LONG_DBL == 1
            rpt->mr[hh][ii] = ( long double ) ( hh == ll );
# else
            rpt->mr[hh][ii] = ( double ) ( hh == ll );
# endif
            rpt->mi[hh][ii] = ZERO; 
         };
      };
   }
/*............................................................................*/
/* option: linear equation  */

   else if ((( gjp->opt ) == 'e' )
          ||(( gjp->opt ) == 'E' ))
   {
      n1 = nn + ONE;
   }
/*............................................................................*/
/* option: determinant */

   else if ((( gjp->opt ) == 'd' )
          ||(( gjp->opt ) == 'D' ))
   {
      n1 = nn;
   }
/*............................................................................*/
/* error message for unknown or unspecified option: */

   else
   {
      printf( "\n\n Error message from function %s:", __func__ );
      printf( "\n\n Unknown or unspecified option !!!" );
      printf( "\n The legal options are [ case insensitive ]:"
         "\n opt = 'e' - linear equation"
         "\n opt = 'i' - matrix inversion"
         "\n opt = 'd' - determinant."
         "\n [ Please specify option on function call.]\n " );
      return NULL;
   };
/*............................................................................*/
/* pivote: */

   for ( hh=null; hh<nn; hh++ )  
   {
     h1 = hh + ONE;
     ll = nn - ONE;

      if ( hh == ll ) 
         goto check_sgl;

      uu = rpt->mr[hh][hh];
      vv = rpt->mi[hh][hh];
      uu = uu*uu + vv*vv;
/*............................................................................*/
/* compare magnitude of elements in lines hh and jj: */

      jj = hh;
      for ( ii=h1; ii<nn; ii++ )
      {
         rr = ( rpt->mr[ii][hh] );
         ss = ( rpt->mi[ii][hh] );
         rr = rr*rr + ss*ss; 
         if ( uu < rr )
         {
            uu = rr;
            jj = ii;
         };
      };

      if ( jj == hh ) 
         goto check_sgl;

      ( rpt->dtr ) *= ( -ONE ); /* interchange lines  hh <-> jj */
      ( rpt->dti ) *= ( -ONE ); /* [ changes determinant sign ] */

      for ( mm=hh; mm<n1; mm++ )
      {
         rr = ( rpt->mr[hh][mm] );
         ss = ( rpt->mi[hh][mm] );

         ( rpt->mr[hh][mm] ) = ( rpt->mr[jj][mm] );
         ( rpt->mi[hh][mm] ) = ( rpt->mi[jj][mm] );
         ( rpt->mr[jj][mm] ) = rr;
         ( rpt->mi[jj][mm] ) = ss;
      };
/*............................................................................*/
/* Gauss - Jordan iteration: */

     check_sgl:

      uu = ( rpt->mr[hh][hh] );
      vv = ( rpt->mi[hh][hh] );
      uu = uu*uu + vv*vv;
     
      if ( uu < GSS_SINGLR ) 
         goto singular;

      for ( ii=h1; ii<nn; ii++ )
      {
         rr = (( rpt->mr[ii][hh] )*( rpt->mr[hh][hh] ) +\
               ( rpt->mi[ii][hh] )*( rpt->mi[hh][hh] )) / uu;
         ss = (( rpt->mi[ii][hh] )*( rpt->mr[hh][hh] ) -\
               ( rpt->mr[ii][hh] )*( rpt->mi[hh][hh] )) / uu; 

         ( rpt->mr[ii][hh] ) = ZERO;
         ( rpt->mi[ii][hh] ) = ZERO;

         for ( mm=h1; mm<n1; mm++ )
         {
            u1 = ( rpt->mr[ii][mm] );
            v1 = ( rpt->mi[ii][mm] );
            u2 = ( rpt->mr[hh][mm] );
            v2 = ( rpt->mi[hh][mm] );

            ( rpt->mr[ii][mm] ) = u1 - rr*u2 + ss*v2;
            ( rpt->mi[ii][mm] ) = v1 - ss*u2 - rr*v2;
         };
      };  
   }; /* next hh */
/*............................................................................*/
/* compute solution zr + j*zi */

   for ( hh=nn-ONE; null<=hh; hh-- ) 
   {
      uu = ( rpt->dtr );
      vv = ( rpt->dti );

      ( rpt->dtr ) = uu*( rpt->mr[hh][hh] ) - vv*( rpt->mi[hh][hh] );
      ( rpt->dti ) = vv*( rpt->mr[hh][hh] ) + uu*( rpt->mi[hh][hh] );

      for ( ii=nn; ii<n1; ii++ )
      {
         jj = ii - nn;

         ( rpt->zr[hh][jj] ) = ( rpt->mr[hh][ii] );
         ( rpt->zi[hh][jj] ) = ( rpt->mi[hh][ii] );

         for ( kk=hh+ONE; kk<nn; kk++ )
         {
            uu = ( rpt->zr[hh][jj] );
            vv = ( rpt->zi[hh][jj] );
            u1 = ( rpt->zr[kk][jj] );
            v1 = ( rpt->zi[kk][jj] );

            ( rpt->zr[hh][jj] ) = uu - u1*( rpt->mr[hh][kk] ) +\
                                       v1*( rpt->mi[hh][kk] );
            ( rpt->zi[hh][jj] ) = vv - v1*( rpt->mr[hh][kk] ) -\
                                       u1*( rpt->mi[hh][kk] );
         };
             
         rr = ( rpt->mr[hh][hh] );
         ss = ( rpt->mi[hh][hh] );
         rr = rr*rr + ss*ss ;

         uu = ( rpt->zr[hh][jj] );
         vv = ( rpt->zi[hh][jj] );

         ( rpt->zr[hh][jj] ) =\
             ( uu*( rpt->mr[hh][hh] ) + vv*( rpt->mi[hh][hh] )) / rr ; 
         ( rpt->zi[hh][jj] ) =\
             ( vv*( rpt->mr[hh][hh] ) - uu*( rpt->mi[hh][hh] )) / rr ; 
      };
   }; /* next hh */

   return rpt;
/*............................................................................*/
/* display singular matrix: */

  singular: 

   ( rpt->dtr ) = ZERO;
   ( rpt->dti ) = ZERO;

   if ((( gjp->opt ) == 'd' )||
       (( gjp->opt ) == 'D' ))
      return rpt;
/*............................................................................*/
# if GSS_DISP == 1
   printf( "\n\n Message from function %s:\n", __func__ );
   ii = null; do
   {
      printf( "\n Singular coefficient matrix "
         "M = ( rpt->mr[][] + i*rpt->mi[][] ) !!!" ); 
      ii++ ;
   } while ( ii < THREE );
   printf( "\n\n [ Abnormal end of function %s.]\n ", __func__ );
# endif

   ( rpt->rtn ) = null;
   return rpt;
} 
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/*============================================================================*/
# undef GSS_DISP
# undef GSS_STRUCT
# undef GSS_LONG_DBL
# undef GSS_MAXRNK
# undef GSS_SINGLR
/*********************** end of function gssjrd(*) ****************************/
