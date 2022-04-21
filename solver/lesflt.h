/* [ file: lesflt.b ] */
/*******************************************************************************
*                                                                              *
*   Function body lesfilter(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   A simple Large Eddy Simulation filter                                      *
*   [ based on cellular coarse graining ]                                      *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/* defaults */
/*----------------------------------------------------------------------------*/
# ifndef LES_EPTLTOT
/* [ don't ] compute total pot. energy in fluid components: [0] 1 */
/* default: 0  */
   # define LES_EPTLTOT 1
# endif /* not defined LES_EPTLTOT */
/*----------------------------------------------------------------------------*/
# ifndef LES_PRSSHFT
/* [ don't ] enable pressure shiftr: [0] 1 */
/* in general recommended: 1 [ enabled ] */
   # define LES_PRSSHFT 1
# endif /* not defined LES_PRSSHFT */
/*----------------------------------------------------------------------------*/
# if CMPRSSBL != 0
/*............................................................................*/
# ifndef LES_CRSNDNS
/* [ don't ] enable density coarsening: [0] 1 */
/* highly recommended: 1 [ enabled ] */
   # define LES_CRSNDNS 1
# endif /* not defined LES_CRSNDNS */
/*----------------------------------------------------------------------------*/
# ifndef LES_MASSCRR
/* [ don't ] enable mass correction [0] 1 */
/* in general recommended: 1 [ enabled ] */
   # define LES_MASSCRR 1
# endif /* not defined LES_MASSCRR */
/*............................................................................*/
# elif CMPRSSBL == 0
   # undef LES_CRSNDNS
   # define LES_CRSNDNS 0 
   # undef LES_MASSCRR
   # define LES_MASSCRR 0
# endif /* CMPRSSBL == 0 */
/*============================================================================*/

DSC_HCRRTS
*lesfilter( struct solverstat *state )
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

   static signed char
      crsmrk = 124; /* ASCII character '|'; at coarsening times replaced */
                    /* by 'C' and written into dsc.log */
/*............................................................................*/
   static short
      cc = null; /* fluid [ connected ] component label */
/*
   static double
      uf[FACES][THREE] = {{ZERO}};
*/
   static double
      dr = ZERO,
      unupwd = 1.000e+00, /* 0 <  unupwd <= 1.0 */
      undnwd = 0.000e+00, /* will be set to ( 1. - unupwd ) */
      crstme[NFCNN+ONE] = {ZERO};

# if CNN_LESFLTR == 4 
   static double
      ds = ZERO;
# endif

/*............................................................................*/
# if LES_PRSSHFT != 0
   static double
      prss[NFCNN+ONE] = {ZERO};
# endif /* LES_PRSSHFT != 0 */
/*............................................................................*/
# if CMPRSSBL != 0
   static double
      mass[NFCNN+ONE] = {ZERO};
/*............................................................................*/
# if LES_CRSNDNS == 1
   static double
      rnupwd = 1.000e+00, /* 0 <  rnupwd <= 1.0 */
      rndnwd = 0.000e+00; /* will be set to ( 1. - rnupwd ) */
# endif /* LES_CRSNDNS == 1 */
/*............................................................................*/
# endif /* CMPRSSBL != 0 */
/*............................................................................*/
/* prototypes: */

   double
      fmod( double x, double y );
/*
   char
     *lotos( long, char );
*/
/*----------------------------------------------------------------------------*/
   kk = ( state->hclbl );
   hre = ( state->hre[( int )kk] );
   hsp = ( state->hsp );
/*............................................................................*/
   undnwd = 1. - unupwd;
   dr = ZERO;
/*............................................................................*/
# if (( CMPRSSBL != 0 )\
    &&( LES_CRSNDNS == 1 ))
   rndnwd = 1. - rnupwd;
# endif /* (( CMPRSSBL != 0 )&& */
/*............................................................................*/
/* here starts the job: */
/*............................................................................*/
   lwfcl = ( state->lwfcl );
   upfcl = ( state->upfcl );

   if (( state->fldnn ) == ONE )
   {
      crstme[null] = ZERO;
      cc = null;
      while(( cc++ ) < ( hsp->cnn[null] ))
      {
         crstme[cc] = ZERO;
      };
   } /* end if (( state->fldnn ) == ONE ) INITIALIZATION TERMINATED */
/*............................................................................*/
/* COARSENING [ cellular coarse graining - a simple large eddy filter ] */
/*............................................................................*/
   else /* if ( ONE < ( state->fldnn )) */
   {
      crsmrk = 124; /* ASCII character '|'; at coarsening times replaced */

      cc = null; /* cc: index, fluid connected component [ 0 < cc <= NFCNN ] */
      while(( cc++ ) < ( hsp->cnn[null] ))
      {
         if ( crstme[cc] < ( state->hcdt ))
         {
/*............................................................................*/
# if LES_PRSSHFT != 0
            prss[cc] = ZERO;
# endif /* LES_PRSSHFT != 0 */
/*............................................................................*/
# if LES_EPTLTOT != 0
            ( hsp->eptlt[cc] ) = ZERO;
# endif /* LES_EPTLTOT != 0 */
/*............................................................................*/
# if CMPRSSBL != 0
            mass[cc] = ZERO;
# endif
/*............................................................................*/
            mm = lwfcl;
            while( mm <= upfcl )
            {
               if ( cc == ( hsp->cnn[mm] )) /* fluid type cell within */
               {
                  hm = ( hsp->hh[mm] ); /* hm: thermal s-parameter index */

                  jj = null; do
                  {
/*............................................................................*/
# if CNN_LESFLTR == 1 /* arithmetic mean - simple [ and fast ], yet quite */
                      /* artificial [ and boundary type insensitive ] */ 

                     dr = ZERO;
                     fc = null; do
                     {
                        dr += ( hre->uf[mm][fc][jj] );
                     } while(( ++fc ) < FACES );

                     ( hre->un[mm][jj] ) *= undnwd;
                     ( hre->un[mm][jj] ) += ( unupwd*dr/FACES );

# elif CNN_LESFLTR == 2 /* face area weighted mean - rather natural */
                        /* [ yet still boundary type insensitive ] */

                     dr = ZERO;
                     fc = null; do
                     {
                        dr += (( hre->uf[mm][fc][jj] )*( hsp->fm[hm][fc] ));
                     } while(( ++fc ) < FACES );

                     ( hre->un[mm][jj] ) *= undnwd;
                     ( hre->un[mm][jj] ) += ( unupwd*dr/( hsp->sf[hm] ));

# elif CNN_LESFLTR == 3 /* arithmetic mean - artificial */ 
                        /* [ yet boundary type sensitive ] */ 

                     dr = ZERO;
                     kk = null;
                     fc = null; do
                     {
                        if (( hsp->fctype[mm][fc] ) == THREE ) /* !=free slip */
                           continue;

                        kk++;
                        dr += ( hre->uf[mm][fc][jj] );
                     } while(( ++fc ) < FACES );

                     if ( kk == null )
                        continue;

                     ( hre->un[mm][jj] ) *= undnwd;
                     ( hre->un[mm][jj] ) += ( unupwd*dr/kk );

# elif CNN_LESFLTR == 4 /* face area weighted mean - natural */
                      /* [ and boundary type sensitive ] */ 

                     dr = ZERO;
                     ds = ZERO;
                     fc = null; do
                     {
                        if (( hsp->fctype[mm][fc] ) == THREE ) /* !=free slip */
                           continue;

                        dr += (( hre->uf[mm][fc][jj] )*( hsp->fm[hm][fc] ));
                        ds += ( hsp->fm[hm][fc] );
                     } while(( ++fc ) < FACES );

                     if ( ds <= CNN_TRIVIAL )
                        continue;

                     ( hre->un[mm][jj] ) *= undnwd;
                     ( hre->un[mm][jj] ) += ( unupwd*dr/ds );
# endif /* CNN_LESFLTR == ... */
/*............................................................................*/
                  } while(( ++jj ) < THREE );
/*............................................................................*/
# if (( LES_EPTLTOT != 0 )\
    ||( LES_PRSSHFT != 0 ))

                  dr = ( hre->pn[mm] )*( hsp->vol[hm] );

# endif /* (( LES_EPTLTOT != 0 )\
          ||( LES_PRSSHFT != 0 )) */
/*............................................................................*/
# if LES_EPTLTOT != 0
                  ( hsp->eptlt[cc] ) += fabs( dr );
# endif /* LES_EPTLTOT != 0 */
/*............................................................................*/
# if LES_PRSSHFT != 0
                  prss[cc] += dr;
# endif /* LES_PRSSHFT != 0 */
/*............................................................................*/
# if CMPRSSBL != 0
/*............................................................................*/
# if LES_CRSNDNS == 1
/*............................................................................*/
# if CNN_LESFLTR == 1 /* arithmetic mean - simple [ and fast ], yet quite */
                      /* artificial [ and boundary type insensitive ] */ 

                  dr = ZERO;
                  fc = null; do
                  {
                     dr += ( hre->rf[mm][fc] );
                  } while(( ++fc ) < FACES );

                  ( hre->rn[mm] ) *= rndnwd;
                  ( hre->rn[mm] ) += ( rnupwd*dr/FACES );

# elif CNN_LESFLTR == 2 /* face area weighted mean - rather natural */
                        /* [ yet still boundary type insensitive ] */

                  dr = ZERO;
                  fc = null; do
                  {
                     dr += (( hre->rf[mm][fc] )*( hsp->fm[hm][fc] ));
                  } while(( ++fc ) < FACES );

                  ( hre->rn[mm] ) *= rndnwd;
                  ( hre->rn[mm] ) += ( rnupwd*dr/( hsp->sf[hm] ));

# elif CNN_LESFLTR == 3 /* arithmetic mean - artificial */ 
                        /* [ yet boundary type sensitive ] */ 

                  dr = ZERO;
                  kk = null;
                  fc = null; do
                  {
                     if (( hsp->fctype[mm][fc] ) == THREE ) /* !=free slip */
                        continue;

                     kk++;
                     dr += ( hre->rf[mm][fc] );
                  } while(( ++fc ) < FACES );

                  if ( kk == null )
                     continue;

                  ( hre->rn[mm] ) *= rndnwd;
                  ( hre->rn[mm] ) += ( rnupwd*dr/kk );

# elif CNN_LESFLTR == 4 /* face area weighted mean - natural */
                        /* [ and boundary type sensitive ] */ 

                  dr = ZERO;
                  ds = ZERO;
                  fc = null; do
                  {
                     if (( hsp->fctype[mm][fc] ) == THREE ) /* !=free slip */
                        continue;

                     dr += (( hre->rf[mm][fc] )*( hsp->fm[hm][fc] ));
                     ds += ( hsp->fm[hm][fc] );
                  } while(( ++fc ) < FACES );

                  if ( ds <= CNN_TRIVIAL )
                     continue;

                  ( hre->rn[mm] ) *= rndnwd;
                  ( hre->rn[mm] ) += ( rnupwd*dr/ds );

# endif /* CNN_LESFLTR == ... */
/*............................................................................*/
# endif /* LES_CRSNDNS == 1 */
/*............................................................................*/
                  mass[cc] += ( hre->rn[mm] )*( hsp->vol[hm] );

# endif /* CMPRSSBL != 0 */
/*............................................................................*/
               }; /* end if ( cc == ( hsp->cnn[mm] )) */
               mm++;
            }; /* end while( mm <= upfcl ) */
/*............................................................................*/
# if CMPRSSBL != 0

            if ( ZERO < ( hsp->mass[cc] ))
	       ( hsp->msdr[cc] ) = mass[cc]/( hsp->mass[cc] );

# endif /* CMPRSSBL != 0 */
/*............................................................................*/
# if LES_PRSSHFT != 0
            prss[cc] /= ( hsp->volume[cc] );
# endif /* LES_PRSSHFT != 0 */
/*............................................................................*/
# if (( LES_PRSSHFT != 0 )\
    ||( LES_MASSCRR != 0 ))
/* [ pressure shift and/or mass corrections: ] */

            mm = lwfcl;
            while( mm <= upfcl )
            {
               if ( cc == ( hsp->cnn[mm] ))
               {
/*............................................................................*/
# if LES_PRSSHFT != 0
                  ( hre->pn[mm] ) -=  prss[cc];
# endif /* LES_PRSSHFT != 0 */
/*............................................................................*/
# if LES_MASSCRR != 0
                  ( hre->rn[mm] ) /= ( hsp->msdr[cc] );
# endif /* LES_MASSCRR != 0 */
/*............................................................................*/
                  fc = null; do 
                  {
/*............................................................................*/
# if LES_PRSSHFT != 0
                     ( hre->pf[mm][fc] ) -=  prss[cc]; 
# endif /* LES_PRSSHFT != 0 */
/*............................................................................*/
# if LES_MASSCRR != 0
                     ( hre->rf[mm][fc] ) /= ( hsp->msdr[cc] );
# endif /* LES_MASSCRR != 0 */
/*............................................................................*/
                  } while(( ++fc ) < FACES );
               }; /* end if ( cc == ( hsp->cnn[mm] )) */

               mm++;
            }; /* end while( mm <= upfcl ) */
# endif /* (( LES_PRSSHFT != 0 )\
          ||( LES_MASSCRR != 0 )) */
/*............................................................................*/
            crsmrk = 67; /* ASCII character 'C', marks coarsening in dsc.log */
         }; /* end if ( crstme[cc] < ( state->hcdt )) */

         crstme[cc] += ( state->hcdt );

/* reset crstme[ ]: */

         if (( hsp->crsdt[cc] ) <= crstme[cc] )
            crstme[cc] = fmod( crstme[cc], ( state->hcdt ));

      }; /* end while(( ++cc ) < ( hsp->cnn[null] )) */
   }; /* end if ( ONE < ( state->fldnn )) */

   ( state->lesfltr ) = crsmrk;
   return hre;
 }
/*============================================================================*/
/************************* end of function lesfilter(*) ***********************/
/***************************** end of file lesflt.b ***************************/
