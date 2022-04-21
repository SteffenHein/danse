/* [ file: expmtx.h ] */
# define DO_EXPNT "expmtx(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function expmtx(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given any complex matrix  U = exu.ur[i][j] + j* exu.ui[i][j] of rank       *
*   n > 0  this subroutine returns the exponential matrix                      *
*                                                                              *
*             exp(U)  =  exu.er[i][j] + j* exu.ei[i][j].                       *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
# include <stdio.h>
*/
/*----------------------------------------------------------------------------*/
# include "./maths.h"  /* 'my' computation environment headers */
# include "./consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# define EXP_RANK     10          /* maximum order of matrix U                */
# define EXP_BOUND    ( 1.e-277 ) /* accuracy                                 */
# define EXP_MAXITER  5000        /* maximum number of iterations             */
/*----------------------------------------------------------------------------*/
struct expmtx 
{
   struct expmtx *pt;

   double
    norm, 
    ur[EXP_RANK][EXP_RANK],                 /* r = real part of matrix U      */
    ui[EXP_RANK][EXP_RANK],                 /* i = imag."    "  "             */
    er[EXP_RANK][EXP_RANK],                 /* r = real part of matrix E      */
    ei[EXP_RANK][EXP_RANK];                 /*  ...                           */
};
struct expmtx exu = {null};
/*============================================================================*/

int expmtx( short rank )
{
/* allusions: */
/*
   extern struct expmtx exu;
*/
/* declarations: */

   static double 
      ur0[EXP_RANK][EXP_RANK] = {{null}},
      ui0[EXP_RANK][EXP_RANK] = {{null}},
      ur1[EXP_RANK][EXP_RANK] = {{null}},
      ui1[EXP_RANK][EXP_RANK] = {{null}};

   static double 
      der = ZERO,
      dei = ZERO,
      norm = ZERO,
      bound = ZERO;

   static short 
      i = null,
      j = null,
      k = null;

   static long 
      n = null; 

   double 
      fabs( double x ),
      sqrt( double x ); 
 
/*----------------------------------------------------------------------------*/

   if ( EXP_RANK < rank )
   {
      printf( "\n\n Message from function 'expmtx(*)' : " );
      printf( "\n\n Rank of matrix U = ( exu.ur[][] + j*exu.ui[][] )"
              " exceeds maximum rank ! " );
      printf( "\n [ Change macro EXP_RANK to actual rank: %d ,"
              " in compliance with memory resources. ]", rank );  
      return ONE;
   };

/*............................................................................*/

   exu.norm = ZERO;
    
   for ( i=null ; i<rank ; i++ )
   {
      for ( j=null ; j<rank ; j++ )
      {   
         exu.er[i][j] = ( double ) ( i == j );     /*  write unit matrix   */
         exu.ei[i][j] = ZERO;

         ur0[i][j] = ( double ) ( i == j );
         ui0[i][j] = ZERO;

         norm  =  exu.ur[i][j]*exu.ur[i][j]; 
         norm +=  exu.ui[i][j]*exu.ui[i][j];

         if ( exu.norm < norm )  exu.norm = norm; 
      };
   };

   bound = EXP_BOUND*exu.norm;
   exu.norm = sqrt( exu.norm );

/*............................................................................*/
/* iterate: */

   n = ONE ; do
   {
      norm = ZERO;

      for ( i=null; i<rank; i++ )
      {
         for ( j=null; j<rank; j++ )
         {    
            der = ZERO;  
            dei = ZERO;

            for ( k=null; k<rank ; k++ )
            {
               der += ur0[i][k]*exu.ur[k][j] - ui0[i][k]*exu.ui[k][j];
               dei += ur0[i][k]*exu.ui[k][j] + ui0[i][k]*exu.ur[k][j];
            };

            der /= n;
            dei /= n;
                
            exu.er[i][j] += der;
            exu.ei[i][j] += dei;

            ur1[i][j] = der;
            ui1[i][j] = dei;

            norm += ( der*der + dei*dei ); 
         };/* next j */
            
         for ( j=null; j<rank; j++ )
         {
            ur0[i][j] = ur1[i][j];
            ui0[i][j] = ui1[i][j];
         };/* next j */
      };/* next i */

      n++;

      if ( EXP_MAXITER < n ) 
      {
         printf( "\n\n Message from function 'expmtx(*)': " );
         printf( "\n\n Too many iterations: n > %d"
                 " [ = macro EXP_MAXITER ] !!! ", EXP_MAXITER );
         return ONE;
      }; 

   }  while ( bound < norm );

   return null;
}
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/*============================================================================*/
# undef EXP_MAXITER
/*********************** end of subroutine expmtx.h ***************************/
