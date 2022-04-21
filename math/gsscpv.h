/* [ file: gsscpv.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI-C Gauss-Jordan elimination with full pivoting                         *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given any non-singular complex matrix m[][] of degree n > 0 and an         *
*   array of  m complex vectors a[][j] ( j=0,..,m-1 ), this subroutine         *
*   returns the m solutions z[][j] of the linear equations                     *
*                                                                              *
*        sum( k=0,...,n-1)  m[i][k] * z[k][j] = a[i][j] .                      *
*                                                                              *
*   To this end, m[][] and a[][][] are transferred to the function as          *
*   columns  k=0,...,n-1 and k=n,...,n+m-1, respectively, of the stru-         *
*   ture matrix                                                                *
*                                                                              *
*              ( gss.mr[i][k] + j*gss.mi[i][k] )  ( i=0,...,n-1 ).             *
*                                                                              *
*   The program writes the solutions z[][j] into                               *
*                                                                              *
*               ( gss.zr[][j] + j*gss.zi[][j] )   ( j=0,...,m-1 )              *
*                                                                              *
*   ( option  *opt = 'e' ) .                                                   *
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
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
# include "../math/CPRODUCT.H"
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
static GAUSS_JRD gss = { ZERO };
/*============================================================================*/

GAUSS_JRD *\
gsscpv( short rank, short neqs, char *opt )
{
/* allusions: */
/*
   extern GAUSS_JRD gss;
*/
/* declarations: */

   static GAUSS_JRD *stp = &gss;

# if GSS_LONG_DBL == 1
   static long double
# else
   static double  
# endif
      uu = ZERO,
      vv = ZERO,
      rr = ZERO,
      ss = ZERO,
      sq = ZERO,
      mx = ZERO,
      pr = ZERO,
      pi = ZERO;

   static short 
      ii = null,
      jj = null,
      kk = null,
      nn = null,
      lin = null,
      col = null;

   static short
      ipv[GSS_MAXRNK],
      ixl[GSS_MAXRNK],
      ixc[GSS_MAXRNK];

/*----------------------------------------------------------------------------*/
/* options:                                                                   */
/*                     compute
   opt = 'i'[inversion]   -->   matrix inversion
   opt = 'e'[equation]    -->   linear equation  
   opt = 'd'[determinant] -->   determinant
*/
/*----------------------------------------------------------------------------*/
   if ( GSS_MAXRNK < rank )
   {
      printf( "\n\n Error message from funcion %s:", __func__ );
      printf( "\n\n Rank of matrix 'gss.m[][]' too high !!!" );
      printf( "\n [ The maximum rank is %d = macro GSS_MAXRNK", GSS_MAXRNK );
      printf( "\n   - change macro in compliance "
         "with memory resources. ]\n " );
      return null;
   };

   switch( *opt )
   {
     case 'i': /* option: matrix inversion */
     case 'I':
                                     
      neqs = rank;

      for ( ii=null ; ii<rank ; ii++ )
      {
         for ( jj=null ; jj<rank ; jj++ ) /* [ write unit matrix ] */
         {   
            kk = rank + jj;

# if GSS_LONG_DBL == 1
            gss.mr[ii][kk] = ( long double ) ( ii == jj );
# else
            gss.mr[ii][kk] = ( double ) ( ii == jj );
# endif
            gss.mi[ii][kk] = ZERO;
         };
      };
      break;

     case 'e': /* option: linear equation */
     case 'E':

      if ( neqs <= null )
         neqs = ONE;

      if ( GSS_MAXRNK < neqs )
      {
         printf( "\n\n Error message from funcion %s:", __func__ );
         printf( "\n\n Too many equations !!!" );
         printf( "\n [ The maximum number is %d = macro GSS_MAXRNK ", 
            GSS_MAXRNK );
         printf( "\n   - change macro in compliance "
            "with memory resources. ]\n " );
         return null;
      };
      break;                                    

     case 'd': /* option: determinant */
     case 'D': 

      neqs = null;

      break; 

     default:
      printf( "\n\n Error message from function %s:", __func__ );
      printf( "\n\n Unknown or unspecified option on calling function !!!" );
      printf( "\n Legal are: *option ="
         "\n 'e' - linear equation"
         "\n 'i' - matrix inversion"
         "\n 'd' - determinant."
         "\n [ Please specify option in function call.]\n " );
      return null; 
   };

   nn = rank + neqs;

   gss.dtr = 1.;
   gss.dti = ZERO;

/*............................................................................*/
/* store coefficient matrix [ for recovery ]                                  */

   ii = null;
   while( ii < rank )
   {
      jj = null;
      while( jj < nn )
      {
         gss.zr[ii][jj] = gss.mr[ii][jj];
         gss.zi[ii][jj] = gss.mi[ii][jj];
         jj++ ;
      };
      ii++ ;
   };

/*............................................................................*/
# if GSS_IMPLCT == 1
/* normalize equations:                                                       */

   while( ii < rank )
   {
      col = - ONE;
      mx = ZERO;
      jj = null;
      while( jj < rank )
      {
         uu = gss.mr[ii][jj];
         vv = gss.mi[ii][jj];
         sq = uu*uu + vv*vv; 
         if( mx < sq )
         { 
            col = jj;
            mx = sq;
         };
         jj++ ;
      };
      if( mx < GSS_SINGLR )
      {
         goto restore;
      }
      else
      {
         uu = gss.mr[ii][col];
         vv = gss.mi[ii][col];

         CPRODUCT( uu, vv, gss.dtr, gss.dti, rr, ss );
         gss.dtr = rr;
         gss.dti = ss;

         gss.mr[ii][col] = 1.;
         gss.mi[ii][col] = 0.; 

         jj = null;
         while( jj < nn )
         {
            if( jj != col )
            {
               CQUOTIENT( gss.mr[ii][jj], gss.mi[ii][jj], uu, vv, rr, ss );
               gss.mr[ii][jj] = rr;
               gss.mi[ii][jj] = ss;
            };
            jj++ ;
         };
      };/* end if */
      ii++ ;
   };
# endif 
/*............................................................................*/
/* pivote: */

   jj = null;
   while( jj < rank )
   {
      ipv[jj] = - ONE;
      jj++ ;
   };

   ii = null;
   while( ii < rank )
   {
      mx = ZERO;
      jj = null;
      while( jj < rank )
      {
         if ( ipv[jj] != null ) 
         {
            kk = null;
            while( kk < rank )
            {
               if( ipv[kk] == - ONE )
               {
                  uu = gss.mr[jj][kk];
                  vv = gss.mi[jj][kk];

                  sq = uu*uu + vv*vv; 
                  if ( mx < sq )
                  {
                     mx = sq;
                     lin = jj;
                     col = kk; 
                  };
               }
               else if ( null < ipv[kk] ) /* singular matrix */
               {
                  mx = ZERO;
                  goto restore;
               }; 

               kk++ ;
            }; /* end while( kk < rank ) */  
         }; /* end if ipv[jj] != - ONE */    
         jj++ ;
      }; /* end while( jj < rank ) */

      ipv[col]++ ; 

/* Having got the pivot element, lines may now be interchanged, if necessary, */
/* to put the pivot element on the diagonal. The columns are not physically   */
/* interchanged, only relabeled: ixc[I], the column of the Ith pivot element, */
/* is the Ith column that is reduced, while ixl[I] is the line in which that  */
/* pivot element was originally located. If ixl[I] != ixc[I] than there is an */
/* implied column interchange. With this form of bookkeeping, the solution    */
/* ( gss.zr[][] , gss.zi[][] ) will end up in correct order, and the inverse  */
/* matrix will be scrambled by culumns.                                       */

      if ( lin != col )
      {
         gss.dtr = - gss.dtr;
         gss.dti = - gss.dti;

         jj = null;  
         while( jj < nn )
         {
            uu = gss.mr[lin][jj];
            vv = gss.mi[lin][jj];
            gss.mr[lin][jj] = gss.mr[col][jj];
            gss.mi[lin][jj] = gss.mi[col][jj];
            gss.mr[col][jj] = uu;
            gss.mi[col][jj] = vv;
            jj++ ;
         }; /* end while ... */
      };/* end if lin != col */

      ixl[ii] = lin; /* Now ready to divide the pivot element, */
      ixc[ii] = col; /* located at [lin][col]                  */

      uu = gss.mr[col][col];
      vv = gss.mi[col][col];
      sq = uu*uu + vv*vv; 

      if ( sq < GSS_SINGLR ) /* singular matrix */ 
      {
         mx = ZERO;
         goto restore;
      };

      CPRODUCT( uu, vv, gss.dtr, gss.dti, rr, ss );

      gss.dtr = rr;
      gss.dti = ss;

      pr =   uu / sq;
      pi = - vv / sq; 

      gss.mr[col][col] = 1.;
      gss.mi[col][col] = ZERO;

      jj = null; 
      while( jj < nn )
      {
         CPRODUCT( gss.mr[col][jj], gss.mi[col][jj], pr, pi, uu, vv );
         gss.mr[col][jj] = uu;
         gss.mi[col][jj] = vv;
         jj++ ; 
      }; /* end while ... */

      jj = null;
      while( jj < rank )
      {
         if( jj != col )
         {
            uu = gss.mr[jj][col];
            vv = gss.mi[jj][col];

            kk = null;  
            while( kk < nn )
            {
               CPRODUCT( uu, vv, gss.mr[col][kk], gss.mi[col][kk], rr, ss );
               gss.mr[jj][kk] -= rr;
               gss.mi[jj][kk] -= ss;
               kk++ ;
            };
         };/* end if jj != col */
         jj++ ;
      }; /* end while( jj < rank ) */
      ii++ ;
   };/* end while( ii < rank ) [ end of main loop ] */

   ii = rank - ONE;
   while( null <= ii )
   {
      if( ixl[ii] != ixc[ii] )
      {
         jj = null;
         while( jj < rank )
         {
            uu = gss.mr[jj][ixl[ii]];
            vv = gss.mi[jj][ixl[ii]];
            gss.mr[jj][ixl[ii]] = gss.mr[jj][ixc[ii]]; 
            gss.mi[jj][ixl[ii]] = gss.mi[jj][ixc[ii]];
            gss.mr[jj][ixc[ii]] = uu;
            gss.mi[jj][ixc[ii]] = vv;
            jj++ ;
         };
      };/* end if ... */
      ii-- ;
   };  

/*............................................................................*/
/* restore coeficcient matrix m[i][j]  ( i,j = 0,...,rank+neqs ) */
/* and write m solutions into  z[][k]  ( k = 0, ..., m-1 ) */

  restore:

   ii = null;
   while( ii < rank )
   {
      jj = null;
      while( jj < rank )
      {
         kk = jj + rank;
         gss.mr[ii][jj] = gss.zr[ii][jj];
         gss.mi[ii][jj] = gss.zi[ii][jj];
         gss.zr[ii][jj] = gss.mr[ii][kk];
         gss.zi[ii][jj] = gss.mi[ii][kk];
         gss.mr[ii][kk] = gss.zr[ii][kk];
         gss.mi[ii][kk] = gss.zi[ii][kk];
         jj++ ;
      };
      ii++ ;
   };

   if ( GSS_SINGLR <= mx )
      return stp;

/*............................................................................*/
/* singular coeficient matrix: */

   gss.dtr = ZERO;
   gss.dti = ZERO;

   if (( *opt == 'd' )||( *opt == 'D' ))
      return stp;

# if GSS_DISP == 1
   printf( "\n\n Message from function %s:\n", __func__ );
   ii = null; do
   {
       printf( "\n Singular coefficient matrix "
          "M = ( gss.mr[][] + i*gss.mi[][] ) !!!" );
      ii++ ;
   } while ( ii < THREE );
   printf( "\n\n [ Abnormal end of function %s.]\n ",
      __func__ );
# endif

   return null;
} 
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/*============================================================================*/
# undef CPRODUCT
# undef CQUOTIENT
# undef GSS_LONG_DBL
# undef GSS_MAXRNK
# undef GSS_IMPLCT
# undef GSS_STRUCT
# undef GSS_SINGLR
# undef GSS_DISP
/************************ end of function gsscpv(*) ***************************/
