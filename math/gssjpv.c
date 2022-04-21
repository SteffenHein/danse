/* [ file: gssjpv.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C Gauss-Jordan elimination with full pivoting                         *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
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
*         ( gjp->zr[i][k] + j*( gjp->zi[i][k] ))  ( i,k=0,...,n-1 )            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*----------------------------------------------------------------------------*/
# include <stdio.h>
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
   # define GSS_SINGLR ( 1.e-277 ) /* bound for singularity check */
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
gssjpv( GAUSS_JRD *gjp )
{
/* allusions: - */

/* declarations: */

   static GAUSS_JRD 
     *rpt = NULL;

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
      pi = ZERO,
      dtr = ZERO,
      dti = ZERO;

   static short 
      ii = null,
      jj = null,
      kk = null,
      mm = null,
      nn = null,
      lin = null,
      col = null,
      rank = null;

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
   rank = ( gjp->rank );
/*...........................................................................*/
   if ( GSS_MAXRNK < rank )
   {
      printf( "\n\n Error message from funcion %s:", __func__ );
      printf( "\n\n Rank of matrix ( gjp->m[][] ) too high !!!" );
      printf( "\n [ The maximum rank is %d = macro GSS_MAXRNK", GSS_MAXRNK );
      printf( "\n   - change macro in compliance "
         "with memory resources. ]\n " );

      ( rpt->rtn ) = null;
      return rpt;
   };
/*............................................................................*/
/* the options: */

   switch( gjp->opt )
   {
     case 'i': /* option: matrix inversion */
     case 'I':
                                     
      mm = ( gjp->rank );

      for ( ii=null; ii<rank; ii++ )
      {
         for ( jj=null; jj<rank ; jj++ ) /* write unit matrix */
         {   
            kk = rank + jj;

# if GSS_LONG_DBL == 1
            ( rpt->mr[ii][kk] ) = ( long double ) ( ii == jj );
# else
            ( rpt->mr[ii][kk] ) = ( double ) ( ii == jj );
# endif
            ( rpt->mi[ii][kk] ) = ZERO;
         };
      };
      break;

     case 'e': /* option: linear equation */
     case 'E':

      mm = ( gjp->neqs );

      if ( mm <= null )
         mm = ONE;

      if ( GSS_MAXRNK < mm )
      {
         printf( "\n\n Error message from funcion %s:", __func__ );
         printf( "\n\n Too many equations !!!" );
         printf( "\n [ The maximum number is %d = macro GSS_MAXRNK", 
            GSS_MAXRNK );
         printf( "\n   - change macro in compliance "
            "with memory resources. ]\n " );
         return NULL;
      };
      break;                                    

     case 'd':                                     /* option: determinant     */
     case 'D': 

      mm = null;

      break; 

     default:
      printf( "\n\n Error message from function %s:", __func__ );
      printf( "\n\n Unknown or unspecified option on calling function !" );
      printf( "\n The legal options are [ case insensitive ]:"
         "\n opt = 'e' - linear equation"
         "\n opt = 'i' - matrix inversion"
         "\n opt = 'd' - determinant."
         "\n [ Please specify option on function call.]\n " );

      return NULL; 
   };
/*............................................................................*/
   nn = rank + mm;

   dtr = 1.;
   dti = ZERO;
/*............................................................................*/
/* store coefficient matrix [for recovery in the case of multiple equations]: */

   ii = null;
   while( ii < rank )
   {
      jj = null;
      while( jj < nn )
      {
         ( rpt->zr[ii][jj] ) = ( rpt->mr[ii][jj] );
         ( rpt->zi[ii][jj] ) = ( rpt->mi[ii][jj] );
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
         uu = ( rpt->mr[ii][jj] );
         vv = ( rpt->mi[ii][jj] );
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
         uu = ( rpt->mr[ii][col] );
         vv = ( rpt->mi[ii][col] );

         CPRODUCT( uu, vv, dtr, dti, rr, ss );
         dtr = rr;
         dti = ss;

         ( rpt->mr[ii][col] ) = 1.;
         ( rpt->mi[ii][col] ) = 0.; 

         jj = null;
         while( jj < nn )
         {
            if( jj != col )
            {
               CQUOTIENT(( rpt->mr[ii][jj] ), ( rpt->mi[ii][jj] ), \
                  uu, vv, rr, ss );
               ( rpt->mr[ii][jj] ) = rr;
               ( rpt->mi[ii][jj] ) = ss;
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
                  uu = ( rpt->mr[jj][kk] );
                  vv = ( rpt->mi[jj][kk] );

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
/*............................................................................*/
/* Having got the pivot element, lines may now be interchanged, if necessary, */
/* to put the pivot element on the diagonal. The columns are not physically   */
/* interchanged, only relabeled: ixc[I], the column of the Ith pivot element, */
/* is the Ith column that is reduced, while ixl[I] is the line in which that  */
/* pivot element was originally located. If ixl[I] != ixc[I] than there is an */
/* implied column interchange. With this form of bookkeeping, the solution    */
/* ( rpt->zr[][], rpt->zi[][] ) will end up in correct order, and the inverse */
/* matrix will be scrambled by culumns.                                       */

      if ( lin != col )
      {
         dtr = - dtr;
         dti = - dti;

         jj = null;  
         while( jj < nn )
         {
            uu = ( rpt->mr[lin][jj] );
            vv = ( rpt->mi[lin][jj] );
            ( rpt->mr[lin][jj] ) = ( rpt->mr[col][jj] );
            ( rpt->mi[lin][jj] ) = ( rpt->mi[col][jj] );
            ( rpt->mr[col][jj] ) = uu;
            ( rpt->mi[col][jj] ) = vv;
            jj++ ;
         }; /* end while ... */
      };/* end if lin != col */

      ixl[ii] = lin; /* now ready to divide the pivot element, */
      ixc[ii] = col; /* located at [lin][col]                  */

      uu = ( rpt->mr[col][col] );
      vv = ( rpt->mi[col][col] );
      sq = uu*uu + vv*vv; 

      if ( sq < GSS_SINGLR ) /* singular matrix */ 
      {
         mx = ZERO;
         goto restore;
      };

      CPRODUCT( uu, vv, dtr, dti, rr, ss );

      dtr = rr;
      dti = ss;

      pr =   uu / sq;
      pi = - vv / sq; 

      ( rpt->mr[col][col] ) = 1.;
      ( rpt->mi[col][col] ) = ZERO;

      jj = null; 
      while( jj < nn )
      {
         CPRODUCT(( rpt->mr[col][jj] ), ( rpt->mi[col][jj] ), \
            pr, pi, uu, vv );
         ( rpt->mr[col][jj] ) = uu;
         ( rpt->mi[col][jj] ) = vv;
         jj++ ; 
      }; /* end while ... */

      jj = null;
      while( jj < rank )
      {
         if( jj != col )
         {
            uu = ( rpt->mr[jj][col] );
            vv = ( rpt->mi[jj][col] );

            kk = null;  
            while( kk < nn )
            {
               CPRODUCT( uu, vv, ( rpt->mr[col][kk] ), ( rpt->mi[col][kk] ), \
                  rr, ss );
               ( rpt->mr[jj][kk] ) -= rr;
               ( rpt->mi[jj][kk] ) -= ss;
               kk++ ;
            };
         }; /* end if jj != col */
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
            uu = ( rpt->mr[jj][ixl[ii]] );
            vv = ( rpt->mi[jj][ixl[ii]] );
            ( rpt->mr[jj][ixl[ii]] ) = ( rpt->mr[jj][ixc[ii]] ); 
            ( rpt->mi[jj][ixl[ii]] ) = ( rpt->mi[jj][ixc[ii]] );
            ( rpt->mr[jj][ixc[ii]] ) = uu;
            ( rpt->mi[jj][ixc[ii]] ) = vv;
            jj++ ;
         };
      };/* end if ... */
      ii-- ;
   };  
/*............................................................................*/
/* restore coeficcient matrix m[i][j]  ( i,j = 0,..., rank+mm )               */
/* and write mm solutions into  z[][k]  ( k = 0, ..., mm-1 )                  */

  restore:

   ii = null;
   while( ii < rank )
   {
      jj = null;
      while( jj < rank )
      {
         kk = jj + gjp->rank;
         ( rpt->mr[ii][jj] ) = ( rpt->zr[ii][jj] );
         ( rpt->mi[ii][jj] ) = ( rpt->zi[ii][jj] );
         ( rpt->zr[ii][jj] ) = ( rpt->mr[ii][kk] );
         ( rpt->zi[ii][jj] ) = ( rpt->mi[ii][kk] );
         ( rpt->mr[ii][kk] ) = ( rpt->zr[ii][kk] );
         ( rpt->mi[ii][kk] ) = ( rpt->zi[ii][kk] );
         jj++ ;
      };
      ii++ ;
   };
/*............................................................................*/
   if ( GSS_SINGLR <= mx )
   {
      ( rpt->dtr ) = dtr;
      ( rpt->dti ) = dti;
      ( rpt->rtn ) = ONE;

      return rpt;
   };
/*............................................................................*/
/* singular coeficient matrix:                                                */

   ( rpt->dtr ) = ZERO;
   ( rpt->dti ) = ZERO;

   if ((( gjp->opt ) == 'd' )||(( gjp->opt ) == 'D' ))
   {
      ( rpt->rtn ) = ONE;
      return rpt;
   };

# if GSS_DISP == 1
   printf( "\n\n Message from function %s:\n", __func__ );
   ii = null; do
   {
       printf( "\n Singular coefficient matrix "
          "M = ( gjp->mr[][] + i*gjp->mi[][] ) !!!" );
      ii++ ;
   } while ( ii < THREE );
   printf( "\n\n [ Abnormal end of function %s.]\n ",
      __func__ );
# endif

   ( rpt->rtn ) = null;
   return rpt;
} 
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/*----------------------------------------------------------------------------*/
# undef CPRODUCT
# undef CQUOTIENT
# undef GSS_LONG_DBL
# undef GSS_MAXRNK
# undef GSS_IMPLCT
# undef GSS_STRUCT
# undef GSS_SINGLR
# undef GSS_DISP
/************************ end of function gssjpv(*) ***************************/
