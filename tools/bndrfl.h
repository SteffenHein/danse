/* [ file: bndrfl.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function bndrfl(*); DANSE release 1.0.                              *
*   [ in DANSE.C program package ]                                             *
*                                                                              *
*   Given any cell label cll of a defined mesh cell  [ i.e. with  yet          *
*   determined vertex points  cor.c[top.cm[cll][]][] ] -  and given a          *
*   cell face label fce [ 0 <= fce < 6 ] and any surface conductivity          *
*   sigma  [ unit:  Siemens = 1./Ohms ],  this  function  returns the          *
*   pertinent face reflection coefficients bnd.r00[i],...,bnd.r11[i],         *
*   where i=bnd.n is the actual boundary face label.                           *
*   The cell and face labels , cll and fce, are transferred to struct          *
*   boundaries 'bnd' [ cf. function formdrv(*) ] as                            *
*                  bnd.m[i] = cll  and  bnd.f[i] = fce ,                       *
*   respectively , for  i = bnd.n ;                                            *
*   Past these operations , label bnd.n is raised by ONE before being          *
*   returned to the calling program .                                          *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifndef EPS_VAC
   # define EPS_VAC ( 8.8541878170e-12 ) /* vacuum permittivity [A*sec/(V*m)] */
# endif

# ifndef MY_VAC_
   # define MY_VAC_ ( 1.2566370614e-06 ) /* vacuum permeability [V*sec/(A*m)] */
# endif
/*----------------------------------------------------------------------------*/
# define BNDRFL_RNDOFF ( 1.e-14 ) /* if def: roundoff bound for refl.coeffs.*/
/*----------------------------------------------------------------------------*/
# include "../math/gssjtp.h"
/*----------------------------------------------------------------------------*/
static GAUSS_JRD gss = {null};
/*============================================================================*/

long bndrfl( long cll, char fce, double sigma )
{
/* allusions: */
/*   
   extern FORMSTATE *spt; 
   extern GAUSS_JRD gss;
*/
/* declarations: */

   static GAUSS_JRD
     *gsj = &gss;

   static double 
      adm = ZERO,  /* char. vacuum field admittance */ 
     scpr = ZERO,  /* scalar product < p1 | p2 >    */ 
     qscp = ZERO,  /* scpr^2                        */
    qnrm1 = ZERO,  /* < p1 | p1 >                   */
    qnrm2 = ZERO,  /* < p2 | p2 >                   */
    alfa1 = ZERO,
    alfa2 = ZERO,
    beta1 = ZERO,
    beta2 = ZERO,
       xx = ZERO,
       yy = ZERO;

   static long     
    bdr = null,
    cp0 = null,
    cp1 = null,
    cp2 = null,
    cp3 = null,
    cp4 = null,
    cp5 = null,
    cp6 = null,
    cp7 = null;

   static int
    ii = null,
    jj = null;

/* functions [ prototypes ] : */ 
   
   double sqrt( double x );

   GAUSS_JRD
     *gssjrd( GAUSS_JRD *gjp );
/*----------------------------------------------------------------------------*/
   if ( BNDAP <= ( spt->bpt->n ))
   {
      fprintf( stderr, "\n\n Message from function %s:", __func__ ); 
      fprintf( stderr, "\n Too many boundary faces !!!" );
      fprintf( stderr, "\n [ number exceeds macro BNDAP-1 = %ld "
         "in function %s.", ( long ) BNDAP, "FORMER.CONF" );
      fprintf( stderr, "\n - Change macro in compliance with memory resources.]\n" );

      exit( EXIT_FAILURE );
   };
/*............................................................................*/
/* vacuum field admittance: */

   adm = sqrt( EPS_VAC / MY_VAC_ );
/*............................................................................*/
/* vertex point indices on boundary face fce:                                 */

   switch( fce )
   {
     case 0: 
      cp0 = ( spt->tpt->cm[cll][2] );
      cp1 = ( spt->tpt->cm[cll][6] );
      cp2 = ( spt->tpt->cm[cll][0] );
      cp3 = ( spt->tpt->cm[cll][4] );
      cp4 = ( spt->tpt->cm[cll][4] );
      cp5 = ( spt->tpt->cm[cll][6] );
      cp6 = ( spt->tpt->cm[cll][0] );
      cp7 = ( spt->tpt->cm[cll][2] );
      break;

     case 1:    
      cp0 = ( spt->tpt->cm[cll][5] );
      cp1 = ( spt->tpt->cm[cll][7] );
      cp2 = ( spt->tpt->cm[cll][1] );
      cp3 = ( spt->tpt->cm[cll][3] );
      cp4 = ( spt->tpt->cm[cll][3] );
      cp5 = ( spt->tpt->cm[cll][7] );
      cp6 = ( spt->tpt->cm[cll][1] );
      cp7 = ( spt->tpt->cm[cll][5] );
      break;

     case 2:
      cp0 = ( spt->tpt->cm[cll][4] );
      cp1 = ( spt->tpt->cm[cll][5] );
      cp2 = ( spt->tpt->cm[cll][0] );
      cp3 = ( spt->tpt->cm[cll][1] );
      cp4 = ( spt->tpt->cm[cll][1] );
      cp5 = ( spt->tpt->cm[cll][5] );
      cp6 = ( spt->tpt->cm[cll][0] );
      cp7 = ( spt->tpt->cm[cll][4] );
      break;

     case 3:
      cp0 = ( spt->tpt->cm[cll][3] );
      cp1 = ( spt->tpt->cm[cll][7] );
      cp2 = ( spt->tpt->cm[cll][2] );
      cp3 = ( spt->tpt->cm[cll][6] );
      cp4 = ( spt->tpt->cm[cll][6] );
      cp5 = ( spt->tpt->cm[cll][7] );
      cp6 = ( spt->tpt->cm[cll][3] );
      cp7 = ( spt->tpt->cm[cll][2] );
      break;

     case 4:
      cp0 = ( spt->tpt->cm[cll][1] );
      cp1 = ( spt->tpt->cm[cll][3] );
      cp2 = ( spt->tpt->cm[cll][0] );
      cp3 = ( spt->tpt->cm[cll][2] );
      cp4 = ( spt->tpt->cm[cll][2] );
      cp5 = ( spt->tpt->cm[cll][3] );
      cp6 = ( spt->tpt->cm[cll][0] );
      cp7 = ( spt->tpt->cm[cll][1] );
      break;

     case 5:
      cp0 = ( spt->tpt->cm[cll][6] );
      cp1 = ( spt->tpt->cm[cll][7] );
      cp2 = ( spt->tpt->cm[cll][4] );
      cp3 = ( spt->tpt->cm[cll][5] );
      cp4 = ( spt->tpt->cm[cll][5] );
      cp5 = ( spt->tpt->cm[cll][7] );
      cp6 = ( spt->tpt->cm[cll][4] );
      cp7 = ( spt->tpt->cm[cll][6] );
      break;
   };
/*............................................................................*/
/* compute scpr = < p1 | p2 >  and  qnrm[i] = < p[i] | p[i] > , i=1,2         */
/* port vectors p1, p2 in 'parcel twines' scheme */
/* ( p1, p2, internal face normal ) in positive orientation */

   scpr  = ZERO;    
   qnrm1 = ZERO;
   qnrm2 = ZERO;

   for ( ii=null; ii<DIMNS; ii++ )
   {
      xx =  ( spt->ppt->cpt->c[cp0][ii]+spt->ppt->cpt->c[cp1][ii] );
      xx -= ( spt->ppt->cpt->c[cp2][ii]+spt->ppt->cpt->c[cp3][ii] );
      xx /= 2.;
      qnrm1 += ( xx * xx );
     
      yy =  ( spt->ppt->cpt->c[cp4][ii]+spt->ppt->cpt->c[cp5][ii] );
      yy -= ( spt->ppt->cpt->c[cp6][ii]+spt->ppt->cpt->c[cp7][ii] );
      yy /= 2.;
      qnrm2 += ( yy * yy );
      scpr += ( xx * yy );
   };

   qscp = scpr * scpr;      /* qscp = < p1 | p2 >^2 */

   xx = 1./ sqrt ( qnrm1*qnrm2 - qscp );

   alfa2 = xx * scpr;
   alfa1 = - alfa2;
   beta1 = xx * qnrm1;
   beta2 = - xx * qnrm2;
/* ...........................................................................*/
/* set up the equations for boundary refl. factors x1=r11,..., x4=r22         */
/* coefficient matrix:                                                        */

   for ( ii=null; ii<FOUR ; ii++ )                            /* clear matrix */
   {
      for ( jj=null; jj<FIVE; jj++ )
      {
         gss.mr[ii][jj] = ZERO;
         gss.mi[ii][jj] = ZERO;
      };
   };

   gss.mr[0][0] = sigma*alfa1;
   gss.mr[0][2] = sigma*beta1 + adm; 
   gss.mr[0][4] = - sigma*alfa1;
    
   gss.mr[1][1] = sigma*alfa1;
   gss.mr[1][3] = sigma*beta1 + adm;
   gss.mr[1][4] = - sigma*beta1 + adm;

   gss.mr[2][0] = sigma*beta2 - adm;
   gss.mr[2][2] = sigma*alfa2;
   gss.mr[2][4] = - sigma*beta2 - adm;

   gss.mr[3][1] = sigma*beta2 - adm;
   gss.mr[3][3] = sigma*alfa2;
   gss.mr[3][4] = - sigma*alfa2;

/*............................................................................*/
/* solve these linear equations: */

   ( gsj->rank ) = FOUR;
   ( gsj->opt ) = 'e'; /* option = 'e'quation */
/*............................................................................*/
   gsj = gssjrd( gsj );    /* Gauss-Jordan elimination                        */
/*.......................*/
/* transfer solution: */

# ifdef BNDRFL_RNDOFF

   ii = null;
   do
   {
      if ( fabs( gss.zr[ii][null] ) < BNDRFL_RNDOFF )
         gss.zr[ii][null] = ZERO;
   } while (( ++ii ) < FOUR );

# endif

   bdr = ( spt->bpt->n );

   ii = fce - 2*( int ) ( fce/2 );          /* ii = 0 [1] if fce even [ odd ] */
   
   switch( ii )
   {
     case 0:
      ( spt->bpt->r00[bdr] ) = gss.zr[0][null];
      ( spt->bpt->r01[bdr] ) = gss.zr[1][null]; 
      ( spt->bpt->r10[bdr] ) = gss.zr[2][null];
      ( spt->bpt->r11[bdr] ) = gss.zr[3][null];
      break;

     case 1:
      ( spt->bpt->r11[bdr] ) = gss.zr[0][null];
      ( spt->bpt->r10[bdr] ) = gss.zr[1][null];
      ( spt->bpt->r01[bdr] ) = gss.zr[2][null];
      ( spt->bpt->r00[bdr] ) = gss.zr[3][null];
      break;
   };

   ( spt->bpt->m[bdr] ) = cll;
   ( spt->bpt->f[bdr] ) = fce;

   ( spt->bpt->n ) = ( ++bdr );

   return bdr;
}
# undef BNDRFL_RNDOFF
/************************* end of function bndrfl(*) **************************/
