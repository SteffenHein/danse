/* [ file: gaussc.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C  function gaussc(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given any non-singular complex matrix  M  of degree  n > 0  and an         *
*   array of  n  complex numbers  a[i] ( i=0,..,n-1 ), this subroutine         *
*   returns the solution  z  of the linear equation                            *
*                                                                              *
*                sum( k=0,...,n-1)  M[i][k] * z[k]  =  a[i] .                  *
*                                                                              *
*   To this end M and a[i] ( i=0,...,n-1 ) are transferred to the subroutine   *
*   as  columns  k=0,...,n-1 and k=n, respectively, of the structure           *
*   matrix                                                                     *
*                                                                              *
*              ( gss.mr[i][k] + j*gss.mi[i][k] )  ( i=0,...,n-1 ).             *
*                                                                              *
*   The program writes the solution  z[k]  into                                *
*                                                                              *
*              ( gss.zr[k][0] + j*gss.zi[k][0] )  ( k=0,...,n-1 )              *
*                                                                              *
*   ( option *opt = 'e' ) .                                                    *
*                                                                              *
*   Also, in option *opt = 'i', the inverse matrix of M,  M^(-1),              *
*   is computed and written into                                               *
*                                                                              *
*              ( gss.zr[i][k] + j*gss.zi[i][k] )  ( i,k=0,...,n-1 )            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
# include <stdio.h>
*/
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# include "./maths.h"  /* 'my' computation environment headers */
# include "./consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H"
/*----------------------------------------------------------------------------*/
/* operational constants: */

# ifndef GSS_STRUCT
   # define GSS_STRUCT 0   /* 1: define struct gaussc gss [may be done elsw.] */
# endif

# if GSS_STRUCT == 0
   # include "gssjtp.h"
# else
/*----------------------------------------------------------------------------*/
/* The structure type definition of Gauss-Jordan algorithm functions          */
/* gaussc(*) [line pivoting] and gsscpv(*) [full pivoting ]                   */
/* Update: 11 August 2001                                                     */
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
gaussc( short rank, char *opt )
{
/* allusions: */
/*
   extern GAUSS_JRD gss;
*/
/* declarations: */

   static GAUSS_JRD *stp = &gss; 

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
      h1 = null,
      n1 = null;

/*----------------------------------------------------------------------------*/
/* options:                                                                   */
/*                     compute
   opt = 'i'[inversion]   -->   matrix inversion
   opt = 'e'[equation]    -->   linear equation  
   opt = 'd'[determinant] -->   determinant
*/
/*----------------------------------------------------------------------------*/
   gss.dtr = ONE;
   gss.dti = ZERO;

   if ( GSS_MAXRNK < rank )
   {
      printf( "\n\n Error message from funcion %s:", __func__ );
      printf( "\n\n Rank of matrix 'gss.m[][]' too high !!!" );
      printf( "\n [ The maximum rank is %d = macro GSS_MAXRNK", GSS_MAXRNK );
      printf( "\n   - change macro in compliance "
         "with memory resources. ]\n " );
      return null;
   };

   if (( *opt == 'i' )||( *opt == 'I' ))         /* option: matrix inversion */
   {
      n1 = TWO*rank;
      for ( hh=null ; hh<rank ; hh++ )            /* [ write unit matrix ]    */
      {
         for ( ii=rank ; ii<n1 ; ii++ )
         {   
            ll = ii-rank;

# if GSS_LONG_DBL == 1
            gss.mr[hh][ii] = ( long double ) ( hh == ll );
# else
            gss.mr[hh][ii] = ( double ) ( hh == ll );
# endif
            gss.mi[hh][ii] = ZERO; 
         };
      };
   }
   else if (( *opt == 'e' )||( *opt == 'E' ))   /* option: linear equation  */
   {
      n1 = rank + ONE;
   }
   else if (( *opt == 'd' )||( *opt == 'D' ))   /* option: determinant      */
   {
      n1 = rank;
   }
   else
   {
      printf( "\n\n Error message from function %s:", __func__ );
      printf( "\n\n Unknown or unspecified option on calling function !!!" );
      printf( "\n Legal options are [ case insensitive ]:"
         "\n opt = 'e' - linear equation"
         "\n opt = 'i' - matrix inversion"
         "\n opt = 'd' - determinant."
         "\n [ Please specify option on function call.]\n " );
      return null; 
   };
/*----------------------------------------------------------------------------*/
/* pivote: */

   for ( hh=null; hh<rank; hh++ )  
   {
     h1 = hh + ONE;
     ll = rank - ONE;

      if ( hh == ll ) 
         goto check_sgl;

      uu = gss.mr[hh][hh];
      vv = gss.mi[hh][hh];
      uu = uu*uu + vv*vv;

/* compare elements of line hh and jj: */

      jj = hh;
      for ( ii=h1; ii<rank; ii++ )
      {
         rr = gss.mr[ii][hh];
         ss = gss.mi[ii][hh];
         rr = rr*rr + ss*ss; 
         if ( uu < rr )
         {
            uu = rr;
            jj = ii;
         };
      };

      if ( jj == hh ) 
         goto check_sgl;

      gss.dtr *= ( -ONE ); /* interchange line  hh <--> jj */
      gss.dti *= ( -ONE ); /* change sign of determinant   */

      for ( mm=hh; mm<n1; mm++ )
      {
         rr = gss.mr[hh][mm];
         ss = gss.mi[hh][mm];

         gss.mr[hh][mm] = gss.mr[jj][mm];
         gss.mi[hh][mm] = gss.mi[jj][mm];
         gss.mr[jj][mm] = rr;
         gss.mi[jj][mm] = ss;
      };

/* Gauss - Jordan iteration: */

     check_sgl:

      uu = gss.mr[hh][hh];
      vv = gss.mi[hh][hh];
      uu = uu*uu + vv*vv;
     
      if ( uu < GSS_SINGLR ) 
         goto singular;

      for ( ii=h1; ii<rank; ii++ )
      {
         rr = ( gss.mr[ii][hh]*gss.mr[hh][hh] + gss.mi[ii][hh]*gss.mi[hh][hh] )/
             uu ;
         ss = ( gss.mi[ii][hh]*gss.mr[hh][hh] - gss.mr[ii][hh]*gss.mi[hh][hh] )/
             uu ;

         gss.mr[ii][hh] = ZERO;
         gss.mi[ii][hh] = ZERO;

         for ( mm=h1; mm<n1; mm++ )
         {
            u1 = gss.mr[ii][mm];
            v1 = gss.mi[ii][mm];
            u2 = gss.mr[hh][mm];
            v2 = gss.mi[hh][mm];

            gss.mr[ii][mm] = u1 - rr*u2 + ss*v2 ;
            gss.mi[ii][mm] = v1 - ss*u2 - rr*v2 ;
         };
      };  
   };/* next hh */

/* compute solution zr + j*zi: */

   for ( hh=rank-ONE; null<=hh; hh-- ) 
   {
      uu = gss.dtr;
      vv = gss.dti;

      gss.dtr = uu*gss.mr[hh][hh] - vv*gss.mi[hh][hh];
      gss.dti = vv*gss.mr[hh][hh] + uu*gss.mi[hh][hh];

      for ( ii=rank; ii<n1; ii++ )
      {
         jj = ii - rank;

         gss.zr[hh][jj] = gss.mr[hh][ii];
         gss.zi[hh][jj] = gss.mi[hh][ii];

         for ( kk=hh+ONE ; kk<rank ; kk++ )
         {
            uu = gss.zr[hh][jj];
            vv = gss.zi[hh][jj];
            u1 = gss.zr[kk][jj];
            v1 = gss.zi[kk][jj];

            gss.zr[hh][jj] = uu - u1*gss.mr[hh][kk] + v1*gss.mi[hh][kk];
            gss.zi[hh][jj] = vv - v1*gss.mr[hh][kk] - u1*gss.mi[hh][kk];
         };
             
         rr = gss.mr[hh][hh];
         ss = gss.mi[hh][hh];
         rr = rr*rr + ss*ss ;

         uu = gss.zr[hh][jj];
         vv = gss.zi[hh][jj];

         gss.zr[hh][jj] = ( uu*gss.mr[hh][hh] + vv*gss.mi[hh][hh] ) / rr ; 
         gss.zi[hh][jj] = ( vv*gss.mr[hh][hh] - uu*gss.mi[hh][hh] ) / rr ; 
      };
   };/* next hh  */

   return stp; 

/* singular matrix: */

  singular: 

   gss.dtr = ZERO;
   gss.dti = ZERO;

   if (( *opt == 'd' )||( *opt == 'D' ))
      return stp;

# if GSS_DISP == 1
   printf( "\n\n Message from function %s:\n", __func__ );
   ii = null; do
   {
      printf( "\n Singular coefficient matrix M"
         " = ( gss.mr[][] + i*gss.mi[][] ) !!!" ); 
      ii++ ;
   } while ( ii < THREE ); 
   printf( "\n\n [ Abnormal end of function %s.]\n ", __func__ );    
# endif

   return null;
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
/*********************** end of function gaussc(*) ****************************/
