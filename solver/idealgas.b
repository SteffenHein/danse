/* [ file: idealgas.b ] */
/*******************************************************************************
*                                                                              *
*   Function body idealgas(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   Returns nodal pressure as a function of temperature and density            *
*   [ Boyle-Mariotte law of ideal gas ]                                        *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# if (( CMPRSSBL != 0 )\
    &&( BSQAPRX == 0 ))
/*----------------------------------------------------------------------------*/
/* defaults: ... */
/*============================================================================*/

DSC_HCRRTS
*idealgas( struct solverstat *state )
{
/*----------------------------------------------------------------------------*/
/* allusions: */
/*     
   extern struct topology top;
   extern struct hcrsmx hcs;
   extern struct boundary bnd;
*/
/* declarations: */
/*............................................................................*/
/* structure pointers: */
/*
*//* presently not used [ but maybe in a future implementation ... ]
*//*
*//*   static struct topology
*//*     *tpt = NULL;
*//*
*//*   static struct tlmsmx
*//*     *spt = NULL;
*//*
*//*   static struct boundary
*//*     *bpt = NULL;
*/
   static struct hcrsmx
     *hsp = NULL;

   static DSC_HCRRTS
     *hre = NULL;
/*............................................................................*/
/*
*//* presently not used [ but maybe in a future implementation ... ]
*//*
*//*   static CLUSTER
*//*     *tmp = NULL,
*//*     *prs = NULL,
*//*     *dns = NULL,
*//*     *flw[THREE] = {NULL},
*//*     *vel[THREE] = {NULL};
*/
/*............................................................................*/
/* register type */

   register long
      hm = null,
      mm = null;

   register signed char 
      jj = null,
      kk = null,
      fc = null;
/*............................................................................*/
/* ... and others */
/*............................................................................*/
# if CNN_DISPL != 0
   static char
      ptr[STS_SIZE] = {null};
# endif
/*............................................................................*/
   static long
      lwfcl = null,
      upfcl = null;
/*............................................................................*/
   static short
      cc = null; /* fluid [ connected ] component label */

   static double
      dr = ZERO,
      ds = ZERO,
      pnupwd = 1.000e+00, /* 0 <  pnupwd <= 1.0 */
      pndnwd; = 0.000e+00, /* will be set to ( 1. - pnupwd ) */
/*............................................................................*/
/*
   double
      fmod( double x, double y );

   char
     *lotos( long, char );
*/
/*
*//* presently not used [ but maybe in a future implementation ... ]
*//*
*//*   CLUSTER 
*//*      *gradnt( CLUSTER *cls, long hh ),
*//*      *gradfc( CLUSTER *cls, long hh );
*/
/*----------------------------------------------------------------------------*/
/* overtake solverstate: */
/*
*//* presently not used [ but maybe in a future implementation ... ]
*//*
*//*   tpt = ( state->tpt );
*//*   spt = ( state->spt );
*//*   bpt = ( state->bpt );
*/
   hsp = ( state->hsp );
   kk = ( state->hclbl );
   hre = ( state->hre[( int )kk] );
/*............................................................................*/
   pndnwd = 1. - pnupwd;
/*............................................................................*/
   tmp = ( state->tmp );
   prs = ( state->prs );
   dns = ( state->dns );
/*
*//*   kk = null; do
*//*   {
*//*      flw[kk] = ( state->flw[kk] );
*//*   } while(( ++kk ) < THREE );
*/
/*............................................................................*/
/* here starts the job: */
/*............................................................................*/
   dr = ZERO;
   ds = ZERO;
/*............................................................................*/
   lwfcl = ( state->lwfcl );
   upfcl = ( state->upfcl );

   cc = null; /* cc: index, fluid connected component [ 0 < cc <= NFCNN ] */
   while(( cc++ ) < ( hsp->cnn[null] ))
   {
      mm = lwfcl;
      while( mm <= upfcl )
      {
         if ( cc == ( hsp->cnn[mm] )) /* fluid type cell within */
         {
            hm = ( hsp->hh[mm] ); /* hm: thermal s-parameter index */

         }; /* end if ( cc == ( hsp->cnn[mm] )) */
         mm++;
      }; /* end while( mm <= upfcl ) */
   }; /end while(( cc++ ) < ( hsp->cnn[null] )) */

   return hre;
 }
/*============================================================================*/
# endif /* (( CMPRSSBL != 0 )\
          &&( BSQAPRX == 0 )) */
/************************* end of function idelgas(*) *************************/
/*************************** end of file idealgas.b ***************************/
