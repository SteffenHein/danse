/* [ file: excfld.h ] */
# define DO_EXCFLD "excfld(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function excfld(*), DANSE release 1.0.                              *
*   [ in DANSE program package ]                                               *
*                                                                              *
*   Given any cell label cll of a defined mesh cell  [ i.e. with  yet          *
*   determined vertex points cor.c[top.cm[cll][*]][*] ] - and given a          *
*   cell face label fce [ 0 <= fce < 6 ] , an E or H-field identifier          *
*   ftype and a field strength vector ( xx, yy, zz ) , this function           *
*   determines the associated excited port identifiers pe[*], ph[*],           *
*   and excitation amplitudes e[*], h[*], which all are written into           *
*   the structure of type EXCITATION pointed to by ept in FORMSTATE            *
*   *spt [ cf. formdrv(*) ].                                                   *
*                                                                              *
*   Through function execution, the pertinent excitation port counter,         *
*   [ ept->ne or ept->nh ] are updated [ adding excited ports up ].            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*----------------------------------------------------------------------------*/
/* macros:                                                                    */

# define E_EXCESS( ) \
{ \
   if (( EXCEP <= spt->ept->ne )&&( ftype == 'e' )) \
   { \
      printf( "\n\n Message from function %s:", DO_EXCFLD ); \
      printf( "\n\n Too many excited E-ports !!!" ); \
      printf( "\n [ number exceeds macro EXCEP-1 = %ld " \
         "in '%s' . ", (long) EXCEP, "FORMER.CONF" ); \
      printf( "\n - Change macro in compliance with memory resources.]\n " ); \
      exit( EXIT_FAILURE ); \
   }; \
}
/*............................................................................*/
# define H_EXCESS( ) \
{ \
   if (( EXCHP <= spt->ept->nh )&&( ftype != 'e' )) \
   { \
      printf( "\n\n Message from function %s:", DO_EXCFLD ); \
      printf( "\n\n Too many excited H-ports !!!" ); \
      printf( "\n [ number exceeds macro EXCHP-1 = %ld " \
         "in '%s' . ", (long) EXCHP, "FORMER.CONF" ); \
      printf( "\n - Change macro in compliance with memory resources. ]\n "); \
      exit( EXIT_FAILURE ); \
   }; \
}
/*............................................................................*/
# define CPRODUCT( a, b, c, d, z, w )                                          \
{                                                                              \
   z = (a)*(c) - (b)*(d);                                                      \
   w = (b)*(c) + (a)*(d);                                                      \
}
/*============================================================================*/
int excfld( long cll, char fce, char ftype,
            double xx, double yy, double zz, double phase )
{
/* allusions: */
/*   
   extern FORMSTATE *spt; 
*/
/* declarations: */

   static double 
    scpr1 = ZERO,  /* scalar product < p1 | field >    */ 
    scpr2 = ZERO,  /* "      "       < p2 | field >    */
    qnrm1 = ZERO,  /* < p1 | p1 >                      */
    qnrm2 = ZERO,  /* < p2 | p2 >                      */
       uu = ZERO,
       vv = ZERO,
       rr = ZERO, 
       ss = ZERO,
       tt = ZERO;

   static double
      vct[THREE] = {ZERO};

   static const double
      bound1 = 1.e-277,
      bound2 = 7.e-14;

   static long     
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
    p1 = null,
    p2 = null;

/* functions [ prototypes ] : */ 
   
   double sqrt( double x );
   double fabs( double x );
   double sin( double x );
   double cos( double x );
/*----------------------------------------------------------------------------*/
   uu = cos( phase );
   vv = sin( phase );

   if (( ftype != 'h' )&&( ftype != 'H' ))
      ftype = 'e';

   E_EXCESS( );

   H_EXCESS( );
/*............................................................................*/
/* vertex point indices on cell face fce:                                     */

   vct[null] = xx;
   vct[ONE]  = yy;
   vct[TWO]  = zz;

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
      p1 = 8;
      p2 = 11;
      break;

     case 1:    
      cp0 = ( spt->tpt->cm[cll][3] );
      cp1 = ( spt->tpt->cm[cll][7] );
      cp2 = ( spt->tpt->cm[cll][1] );
      cp3 = ( spt->tpt->cm[cll][5] );
      cp4 = ( spt->tpt->cm[cll][5] );
      cp5 = ( spt->tpt->cm[cll][7] );
      cp6 = ( spt->tpt->cm[cll][1] );
      cp7 = ( spt->tpt->cm[cll][3] );
      p1 = 6;
      p2 = 9;
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
      p1 = 12;
      p2 =  3;
      break;

     case 3:
      cp0 = ( spt->tpt->cm[cll][6] );
      cp1 = ( spt->tpt->cm[cll][7] );
      cp2 = ( spt->tpt->cm[cll][2] );
      cp3 = ( spt->tpt->cm[cll][3] );
      cp4 = ( spt->tpt->cm[cll][3] );
      cp5 = ( spt->tpt->cm[cll][7] );
      cp6 = ( spt->tpt->cm[cll][2] );
      cp7 = ( spt->tpt->cm[cll][6] );
      p1 = 10;
      p2 =  1;
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
      p1 =  4;
      p2 =  7;
      break;

     case 5:
      cp0 = ( spt->tpt->cm[cll][5] );
      cp1 = ( spt->tpt->cm[cll][7] );
      cp2 = ( spt->tpt->cm[cll][4] );
      cp3 = ( spt->tpt->cm[cll][6] );
      cp4 = ( spt->tpt->cm[cll][6] );
      cp5 = ( spt->tpt->cm[cll][7] );
      cp6 = ( spt->tpt->cm[cll][4] );
      cp7 = ( spt->tpt->cm[cll][5] );
      p1 =  2;
      p2 =  5;
      break;
   };
/*............................................................................*/
/* compute scpr[i] = < p[i], field >  and  qnrm[i] = < p[i] | p[i] > , i=1,2  */
                                 /* port vectors p1,p2: 'parcel twines' scheme*/
                                 /* ( p1, p2 , int. face normal ) pos. orient.*/
   scpr1 = ZERO;    
   qnrm1 = ZERO;
   scpr2 = ZERO;
   qnrm2 = ZERO;

   ii = null; do
   {
      tt  = (( spt->ppt->cpt->c[cp0][ii] ) + ( spt->ppt->cpt->c[cp1][ii] ));
      tt -= (( spt->ppt->cpt->c[cp2][ii] ) + ( spt->ppt->cpt->c[cp3][ii] ));
      tt /= 2.;
      qnrm1 += ( tt * tt );
      scpr1 += ( tt * vct[ii] );

      tt  = (( spt->ppt->cpt->c[cp4][ii] ) + ( spt->ppt->cpt->c[cp5][ii] ));
      tt -= (( spt->ppt->cpt->c[cp6][ii] ) + ( spt->ppt->cpt->c[cp7][ii] ));
      tt /= 2.;
      qnrm2 += ( tt * tt );
      scpr2 += ( tt * vct[ii] );
      ii++;
   } while ( ii < DIMNS );

   qnrm1  = sqrt( qnrm1 );
   qnrm2  = sqrt( qnrm2 );

   if ( bound1 < fabs( scpr1 ))
      if ( fabs( scpr2 / scpr1 ) < bound2 )
         scpr2 = ZERO;

   if ( bound1 < fabs( scpr2 ))
      if ( fabs( scpr1 / scpr2 ) < bound2 )
         scpr1 = ZERO;

   if ( ftype == 'e' ) 
   {
      tt = fabs( scpr1 );
      if ( bound1 < tt )
      {
         ( spt->ept->me[( spt->ept->ne )] ) = cll;
         ( spt->ept->pe[( spt->ept->ne )] ) = p1;

         CPRODUCT( uu, vv, scpr1, ZERO, rr, ss );

         ( spt->ept->er[( spt->ept->ne )] ) = rr;
         ( spt->ept->ei[( spt->ept->ne )] ) = ss;
         ( spt->ept->ne )++;
      };

      tt = fabs( scpr2 );
      if ( bound1 < tt )
      {
         E_EXCESS( );

         ( spt->ept->me[( spt->ept->ne )] ) = cll;
         ( spt->ept->pe[( spt->ept->ne )] ) = p2;

         CPRODUCT( uu, vv, scpr2, ZERO, rr, ss );

         ( spt->ept->er[( spt->ept->ne )] ) = rr;
         ( spt->ept->ei[( spt->ept->ne )] ) = ss;
         ( spt->ept->ne )++;
      };
   }
   else if ( ftype != 'e' )
   {
      tt = fabs( scpr1 );
      if ( bound1 < tt )
      {
         ( spt->ept->mh[( spt->ept->nh )] ) = cll;
         ( spt->ept->ph[( spt->ept->nh )] ) = p1;

         CPRODUCT( uu, vv, scpr1, ZERO, rr, ss );

         ( spt->ept->hr[( spt->ept->nh )] ) = rr;
         ( spt->ept->hi[( spt->ept->nh )] ) = ss;
         ( spt->ept->nh )++;
      };

      tt = fabs( scpr2 );
      if ( bound1 < tt )
      {
         H_EXCESS( );

         ( spt->ept->mh[( spt->ept->nh )] ) = cll;
         ( spt->ept->ph[( spt->ept->nh )] ) = p2;

         CPRODUCT( uu, vv, scpr2, ZERO, rr, ss );

         ( spt->ept->hr[( spt->ept->nh )] ) = rr;
         ( spt->ept->hi[( spt->ept->nh )] ) = ss;
         ( spt->ept->nh )++;
      };
   };
/* ...........................................................................*/
   return ONE;
}
/*============================================================================*/
# undef CPRODUCT
# undef E_EXCESS
# undef H_EXCESS
/************************* end of function excfld(*) **************************/
