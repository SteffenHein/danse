/* [ file: conflow-1.0r3.a ] */
/*******************************************************************************
*                                                                              *
*   Function body conflow(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC fluid flow connection map                                              *
*   [ diffusive heat currents are separately connected in function             *
*     cncthcr(*) ]                                                             *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: March 18, 2022 ]                      <contact@steffen-hein.org> *
*                                                                              *
*******************************************************************************/
/* Macro defaults: */
/*............................................................................*/
/* Different modes for pressure SOR iterations [0: ( 0 < mn ), 1: ( mm < mn ) */
/* 2: forward-backward cycles; forward: ( mm < mn ) - backward: ( mn < mm )] */
/* [ default: 2 ] */
# ifndef CNN_SORMODE
   # define CNN_SORMODE 2
# endif
/*............................................................................*/
/* Secure SOR re-iteration period: */ 
# ifndef CNN_SORREPT
   # define CNN_SORREPT 1
# endif
/*............................................................................*/
/* Maximum number of pressure SOR iterations: */
# ifndef CNN_MAXITRS
   # define CNN_MAXITRS 10
# endif
/*............................................................................*/
/* Minimum number of pressure SOR iterations: */
# ifndef CNN_MINITRS
   # define CNN_MINITRS 3
# endif
/*............................................................................*/
/* Gauss-Seidel: */
# ifndef CNN_RLXPRSS
   # define CNN_RLXPRSS ( 1.00e+00 )
# endif
/*............................................................................*/
/* SOR initial dynamic bound: */
# ifndef CNN_SORDBND
   # define CNN_SORDBND ( 1.00e-02 )
# endif
/*............................................................................*/
/* Lower bound for pressure corrections: */
# ifndef CNN_PRCRBND
   # define CNN_PRCRBND ( 1.00e-77 )
# endif
/*............................................................................*/
/* Dynamic SOR bound updating mode: */
# ifndef CNN_SORUPDT
   # define CNN_SORUPDT 3
# endif
/*............................................................................*/
/* Increasing and decreasing rates [ for dynamic SOR bound updating ] */
# if CNN_SORUPDT != 0
/* Increasing rate: */
# ifndef CNN_SORINCR
   # define CNN_SORINCR ( 1.0300e+00 )
# endif 
/* Decreasing rate: */
# ifndef CNN_SORDECR
/* # define CNN_SORDECR ( 9.0000e-01 ) */
   # define CNN_SORDECR ( 1./CNN_SORINCR )
# endif
/* Strengthening factor */
# ifndef CNN_SORRDCE
   # define CNN_SORRDCE ( 9.9997e-01 )
# endif
# endif /* CNN_SORUPDT != 0 */
/*............................................................................*/
/* von Neumann boundary conitions: */
# ifndef CNN_BDCPRSS
   # define CNN_BDCPRSS 1
# endif
/*............................................................................*/
/* LES filter enabled */
# ifndef CNN_LESFLTR
   # define CNN_LESFLTR 1
   # if CNN_LESFLTR != 0
/* Pressure reset disabled */
      # define CNN_RSTPRSS 0
   # endif
# endif
/*............................................................................*/
/* Flow clipping at domain boundaries disabled: */
# ifndef CNN_BNDCLIP
   # define CNN_BNDCLIP 0
# endif
/*............................................................................*/
/* No-slip mode: Set uf[fc] := ZERO */
# ifndef CNN_NSLPMDE
   # define CNN_NSLPMDE 0
# endif
/*............................................................................*/
/* Free slip mode: Cut normal component */
# ifndef CNN_FRSLMDE
   # define CNN_FRSLMDE 0
# endif
/*............................................................................*/
# ifndef CNN_TRIVIAL
   # define CNN_TRIVIAL ( DBL_MIN )
# endif
/*============================================================================*/

DSC_HCRRTS *\
conflow( struct solverstat *state )
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

   static struct topology
     *tpt = NULL;

   static struct hcrsmx
     *hsp = NULL;

   static struct boundary
     *bpt = NULL;

   static DSC_HCRRTS
     *hci = NULL,
     *hre = NULL;
/*............................................................................*/
   static CLUSTER
     *prs = NULL;
/*............................................................................*/
/* register type */

   register long
      hm = null,
      hn = null,
      mm = null,
      mn = null,
      nn = null;

   register signed char 
      jj = null,
      kk = null,
      ll = null,
      fc = null,
      fn = null;

   register double
      ss0 = ZERO,
      ss1 = ZERO,
      ss2 = ZERO,
      ss3 = ZERO,
      ss4 = ZERO,
      ss5 = ZERO,
      ss6 = ZERO,
      ss7 = ZERO,
      ss8 = ZERO;
/*............................................................................*/
/* ... and others */
/*
   static char
      ptr[STS_SIZE] = {null};
*/
   static long 
      lwfcl = null,
      upfcl = null;

   static signed char
      sgn = null,
      lesmrk = 124; /* ASCII character '|'; at coarsening times replaced */
                    /* by 'C' and written into dsc.log */
/*............................................................................*/
/*
# if (( CMPRSSBL == 0 )\
    ||( BSQAPRX != 0 ))
*/
   static signed char
      sorrept[NFCNN] = {null};

   static short
      cc = null, /* fluid [ connected ] component label */
      itr = null,
      maxitrs = CNN_MAXITRS, /* maximum number of relaxation cyles */
      minitrs = CNN_MINITRS; /* number of initial relaxation cycles */

/* # endif */
/* (( CMPRSSBL == 0 )\
  ||( BSQAPRX != 0 )) */
/*............................................................................*/
   static long
      max = null,
      cdfmx = null,
      cminp = null,
      cmaxp = null;

   static double
      uf[FACES][THREE] = {{ZERO}}; /* face fluid flows */

   static double
      dt = ZERO, /* thermal/fluid time step [ hcdt ] */
      df = ZERO,
      ds = ZERO,
      dv = ZERO,  
      dw = ZERO,
      dru = ZERO,    /* div( ro*u ) */
      dfm = ZERO,    /* maximum df displayed in dsc.log, e.g.*/
      pcrbnd = ZERO, /* tolerated ds */
      actbnd = ZERO, /* df absolute bound [ actual ] */
      dynbnd = ZERO, /* df absolute bound [ dynamic ] */
      minpr = ZERO,  /* pressure minimum [ displayed in prs.log file ] */
      maxpr = ZERO,  /* pressure maximum [ displayed in prs.log file ] */
      dwdw = ZERO;   /* sqrt( SUM_over_all_cells ( dp - dru/dt )^2 ) */

   static double
      rlxprs = CNN_RLXPRSS; /* SOR relaxation factor */
/*............................................................................*/
# if DSC_ADJCLL == 1
   static signed char
      pp = null;
/* 
   faces pertinent to ports 0,1,...,11:
   static const char 
      fce[PORTS] = { 3, 5, 2, 4, 5, 1, 4, 0, 1, 3, 0, 2 };
*/
# endif
/*............................................................................*/
# if CNN_FRSLMDE == 1
   static double
      vf[FACES][THREE] = {{ZERO}};
# endif
/*............................................................................*/
# if CNN_SORUPDT != 0
   static char
      updmrk = null;

   static short
      maxlimit = CNN_SORUPDT,
      minlimit = ( short )( -2*CNN_SORUPDT );

   static double
      sorincr = CNN_SORINCR,
      sordecr = CNN_SORDECR,
      sorrdce = CNN_SORRDCE;
# endif
/*............................................................................*/
/* prototypes: */

   double
      pow( double x, double y );

   double
      fmod( double x, double y );

   double
      trnsvs( double *trv, double *dir, int nn ),
      longtd( double *lgt, double *dir, int nn );
/*
   char
     *lotos( long, char );
*/
   CLUSTER
      *gradnt( CLUSTER *cls, long hh ),
      *gradfc( CLUSTER *cls, long hh );
/*............................................................................*/
# if CNN_LESFLTR != 0
   DSC_HCRRTS
      *lesfilter( struct solverstat *state );
# endif
/*----------------------------------------------------------------------------*/
/* overtake solverstate: */

   tpt = ( state->tpt );
   hsp = ( state->hsp );
   bpt = ( state->bpt );

   jj = ( state->hclbl );

   hci = ( state->hre[( int )jj] );
   hre = ( state->hci[( int )jj] );
/*............................................................................*/
   sgn = null; /* use this at least once [ to satisfy the compiler ] */
/*............................................................................*/
   prs = ( state->prs );
/*............................................................................*/
/* here starts the job: */
/*............................................................................*/
/* INITIALIZE [ parameter reset ]: */

   dt = ( state->hcdt );

   lwfcl = ( state->lwfcl );
   upfcl = ( state->upfcl );

   if (( state->fldnn ) == ONE )
   {
      dynbnd = CNN_SORDBND;  /* dynamic |df| bound */
      dynbnd /= dt;
      actbnd = dynbnd;
      pcrbnd = CNN_PRCRBND;  /* bound for pressure correction */
      pcrbnd /= dt;

      sorrept[null] = null;

      cc = null;
      while(( cc++ ) < ( hsp->cnn[null] ))
      {
         sorrept[cc] = CNN_SORREPT;
/*............................................................................*/
# if CNN_SORUPDT != 0
         updmrk = null;
# endif
/*............................................................................*/
/* initialize fluid density and */
/* reset boundary type identifiers: */

         mm = lwfcl;
         while( mm <= upfcl )
         {
            if ( cc == ( hsp->cnn[mm] ))
            {
               hm = ( hsp->hh[mm] ); /* hm: therm/fluid s-parameter index */
               fc = null; do
               {
/*............................................................................*/
/* default: no-slip [ type TWO ]: */
                 ( hsp->fctype[mm][fc] ) = TWO;
/*............................................................................*/
/* the adjacent mesh cell index */
# if DSC_ADJCLL == 0 /* [ neighbouring cells labelled with faces ] */
                  mn = ( tpt->m[mm][fc] );
# else /* if DSC_ADJCLL == 1: [ neighbouring cells labelled with ports ] */
                  pp = prt1[fc];
                  mn = ( tpt->m[mm][pp] );
# endif
/*............................................................................*/
                  if( null < mn )                /* adjacent cell in the same */
                     if( cc == ( hsp->cnn[mn] )) /* fluid connected component */
                        ( hsp->fctype[mm][fc] ) = null;
               } while(( ++fc ) < FACES );
            }; /* end if ( cc == ( hsp->cnn[mm] )) */
            mm++;
         }; /* end while( mm <= upfcl ) */
      }; /* end while(( cc++ ) < ( hsp->cnn[null] )) */
/*............................................................................*/
/* marking inflow boundary faces [ type ONE ]: */

      nn = null;
      while( nn < ( bpt->nif ))
      {
         mm = ( bpt->mif[nn] );
         fc = ( bpt->fif[nn] );

         ( hsp->fctype[mm][fc] ) = ONE;

         jj = null; do
         {
            ss0 = ( bpt->uf[nn][jj] ); /* initial inflow */

            ( hre->uf[mm][fc][jj] ) = ss0;
            ( hre->un[mm][jj] ) = ss0;

         } while(( ++jj ) < THREE );
         ss0 = ( bpt->ti[nn] ); /* inflow temperature */
         ( hre->tf[mm][fc] ) = ss0;
         ( hre->tn[mm] ) = ss0;
         nn++ ;
      }; /* end while( nn < ( bpt->nif )) */
/*............................................................................*/
/* marking outflow boundary faces [ type -ONE ]: */

      nn = null;
      while( nn < ( bpt->nof ))
      {
         mm = ( bpt->mof[nn] );
         fc = ( bpt->fof[nn] );
         ( hsp->fctype[mm][fc] ) = -ONE;

         jj = null; do
         {
            ss0 = ( bpt->vf[nn][jj] ); /* initial outflow */
            ( hre->uf[mm][fc][jj] ) = ss0;
            ( hre->un[mm][jj] ) = ss0;
         } while(( ++jj ) < THREE );
         nn++ ;
      }; /* end while( nn < ( bpt->nof )) */
/*............................................................................*/
/* marking no-slip boundary faces [ type TWO ] */

      nn = null;
      while( nn < ( bpt->nns ))
      {
         mm = ( bpt->mns[nn] );
         fc = ( bpt->fns[nn] );

         ( hsp->fctype[mm][fc] ) = TWO;
/*............................................................................*/
# if NUSSELT != 0
/* Boundary layer Nusselt number: */
         ( hre->nus[mm] ) = ( bpt->nus[nn] );
# endif
/*............................................................................*/
         jj = null; do
         {
            ( hre->uf[mm][fc][jj] ) = ZERO;
            ( hre->un[mm][jj] ) = ZERO;
         } while(( ++jj ) < THREE );
         nn++ ;
      }; /* end while( nn < ( bpt->nns )) */
/*............................................................................*/
/* marking free slip boundary faces [ type THREE ]: */

      nn = null;
      while( nn < ( bpt->nsl ))
      {
         mm = ( bpt->msl[nn] );
         fc = ( bpt->fsl[nn] );

         ( hsp->fctype[mm][fc] ) = THREE;

         jj = null; do
         {
            ( hre->uf[mm][fc][jj] ) = ZERO;
            ( hre->un[mm][jj] ) = ZERO;
         } while(( ++jj ) < THREE );
         nn++ ;
      }; /* end while( nn < ( bpt->nsl )) */
/*............................................................................*/
/* display initial bounds: */
/*............................................................................*/
/*    format --
*//*  -- SOR log file:
*//*
*//*  fprintf(( state->dscstr ),
*//*     "| iteration " );
*//*  fprintf(( state->dscstr ),
*//*     "| time_[sec] " ); 
*//*  fprintf(( state->dscstr ),
*//*     "| SOR_steps " );
*//*  fprintf(( state->dscstr ),
*//*     "| sum_dw*dw  " );
*//*  fprintf(( state->dscstr ),
*//*     "| max_|df|   " );
*//*  fprintf(( state->dscstr ),
*//*     "| cell_no  |\n" );
*/
/*............................................................................*/
/* store initial bounds */

      fprintf(( state->dscstr ),
         "|%10ld | %.4e | %-4s%5ld | %.4e | %-21s |\n", ( state->fldnn ),
            ( state->hctme ), "MAX:", ( long ) maxitrs , dynbnd,
               "<---INITIAL_BOUND--<<" );

      fflush( state->dscstr );
/*............................................................................*/
/* pressure log file */
/*
*//*  format:
*//*
*//*  fprintf(( state->prsstr ),
*//*     "| iteration " );
*//*  fprintf(( state->prsstr ),
*//*     "| time_[s]   " ); 
*//*  fprintf(( state->prsstr ),
*//*     "| minimum [Pa] " );
*//*  fprintf(( state->prsstr ),
*//*     "|  cell no " );
*//*  fprintf(( state->prsstr ),
*//*     "| maximum [Pa] " );
*//*  fprintf(( state->prsstr ),
*//*     "| cell no  |\n" );
*//*..........................................................................*/
/* store values: */

      fprintf(( state->prsstr ),
         "|%10ld | %.4e | %+.5e |%9ld | %+.5e |%9ld |\n", ( state->fldnn ),
            ( state->hctme ), minpr, cminp, maxpr, cmaxp );

      fflush( state->prsstr );
   } /* end if (( state->fldnn ) == ONE ), INITIALIZATION TERMINATED */
   else /* if ( ONE < ( state->fldnn )) */
   {
/*............................................................................*/
/* CONNECTION CYCLE */
/*............................................................................*/
      dwdw = ZERO;
      dfm = ZERO;
   /* lesmrk = 124; *//* ASCII character '|' */

      cc = null; /* cc: index, fluid connected component [ 0 < cc <= NFCNN ] */
      while(( cc++ ) < ( hsp->cnn[null] ))
      {
/*............................................................................*/
/* CONNECTION CYCLE, FLUID COMPONENT ( cc ) */
/*............................................................................*/
         mm = lwfcl;
         while( mm <= upfcl )
         {
            if ( cc == ( hsp->cnn[mm] ))
            {
               hm = ( hsp->hh[mm] ); /* hm: therm/fluid s-parameter index */
               fc = null; do
               {
/*............................................................................*/
/* the adjacent mesh cell index */
# if DSC_ADJCLL == 0 /* [ neighbouring cells labelled with faces ] */
                  mn = ( tpt->m[mm][fc] );
# else /* if DSC_ADJCLL == 1: [ neighbouring cells labelled with ports ] */
                  pp = prt1[fc];
                  mn = ( tpt->m[mm][pp] );
# endif
/*............................................................................*/
                  if ( mm < mn ) /* [ fc is an interface to neighbour mn ] */
                  {
/* check if neighbour of fluid type in same connected component: */
                     if ( cc == ( hsp->cnn[mn] ))
                     {
/* pertinent s-parameter label */
                        hn = ( hsp->hh[mn] );
/*............................................................................*/
/* adjacent face index: */
# if DSC_FCELBL == 0
                        fn = ( tpt->f[mm][fc] );
# elif DSC_FCELBL == 1
                        fn = ( tpt->f[mm][fc] ) - ONE;
# elif DSC_FCELBL == 2
                        fn = ( tpt->f[mm][fc] );

                        if ( fn < null )
                        {
                           fn = - fn;
                           sgn = - ONE;
                        }
                        else
                           sgn = ONE;

                        fn -= ONE;

# endif /* DSC_FCELBL == 2 */
/*............................................................................*/
                        kk = ( short ) ( fc/TWO );
                        ll = ( short ) ( fn/TWO );
/*............................................................................*/
/* [ s[fc][] = f[fc]*( adj(B)^-1 ) ] */
/*............................................................................*/
/* the mesh cell form vector at face fc: */
                        ss1 = ( hsp->s[hm][fc][kk] );
/* ... and at the adjacent face fn: */
                        ss2 = ( hsp->s[hn][fn][ll] );
/*............................................................................*/
# if CMPRSSBL != 0
                        ss0 = ss1 + ss2;

                        if ( CNN_TRIVIAL < fabs( ss0 ))
                        {
/*............................................................................*/
# if BSQAPRX == 0
/*............................................................................*/
/* CONTINUITY OF PRESSURE [ GRADIENT ] at cell interfaces: */

/* [ rt[fc] = .5*f[fc]*( adj(B)^-1 )*( dr/db ); the trunct'd density gradnt ] */
/*............................................................................*/
/* the truncated gradient at face fc: */
                           ss3 = ( hci->pt[mm][fc] );
/* ... and at the adjacent face fn: */
                           ss4 = ( hci->pt[mn][fn] ); 
/* ss5 = pf[fc] = pf[fn] */
                           ss5 = ( ss3+ss4 )/ss0;

                           ( hre->pf[mm][fc] ) = ss5;
                           ( hre->pf[mn][fn] ) = ss5;
# endif /* BSQAPRX == 0 */
/*............................................................................*/
/* CONTINUITY OF DENSITY [ GRADIENT ] at cell interfaces: */

/*............................................................................*/
/* [ rt[fc] = .5*f[fc]*( adj(B)^-1 )*( dr/db ); the trunct'd density gradnt ] */
/*............................................................................*/
/* the truncated gradient at face fc: */
                           ss3 = ( hci->rt[mm][fc] );
/* ... and at the adjacent face fn: */
                           ss4 = ( hci->rt[mn][fn] );
/* ss5 = rf[fc] = rf[fn] */
                           ss5 = ( ss3+ss4 )/ss0;

                           ( hre->rf[mm][fc] ) = ss5;
                           ( hre->rf[mn][fn] ) = ss5;
                        } /* end if ( CNN_TRIVIAL < fabs ( ss0 )) */
                        else /* trivial ss0 */
                        {
/*............................................................................*/
# if BSQAPRX == 0
                           ( hre->pf[mm][fc] ) = ZERO;
# endif /* BSQAPRX == 0 */
/*............................................................................*/
                           ( hre->rf[mm][fc] ) = ZERO;
                        }; /* trivial ss0 */
# endif /* CMPRSSBL != 0 */
/*............................................................................*/
/* CONTINUITY OF FLOW [ GRADIENTS ] at cell interfaces: */

/*............................................................................*/
# if TURBMOD != 0
/* myt = turbulent dynamic viscosity [ from any model ] */
/* [ Prandtl's e.g.] */
                        ss3 = ss1*( hci->myt[mm] );
                        ss4 = ss2*( hci->myt[mn] );
# else
                        ss3 = ss1*( hsp->ny[hm] );
                        ss4 = ss2*( hsp->ny[hn] );
# endif
/*............................................................................*/
                        ss0 = ss3+ss4;

                        if ( CNN_TRIVIAL < fabs( ss0 ))
                        {
                           jj = null; do
                           {
/*............................................................................*/
/* [ ut[fc][] = .5*f[fc]*( adj(B)^-1 )*( du[]~/db ); the truncated gradient ] */
/*............................................................................*/
/* the truncated gradient at face fc: */
                              ss5 = ( hci->ut[mm][fc][jj] );
/* ... and at the adjacent face fn: */
                              ss6 = ( hci->ut[mn][fn][jj] );
/* ss7 = uf[fc][jj] = uf[fn][jj] */
                              ss7 = ( ss5+ss6 )/ss0;
/* ss8 = < fc | grad( uf[jj] ) > = - < fn | grad( uf[jj] ) > */
                              ss8 = ss6 - ss5 + ss7*( ss3-ss4 );

                              ( hre->uf[mm][fc][jj] ) = ss7;
                              ( hre->uf[mn][fn][jj] ) = ss7;
                              ( hre->iu[mm][fc][jj] ) = ss8;
                              ( hre->iu[mn][fn][jj] ) = - ss8;
                           } while (( ++jj ) < THREE );
                        } /* end if ( CNN_TRIVIAL < fabs ( ss0 )) */
                        else /* trivial ss0 */
                        {
/*............................................................................*/
/* DEFAULT BOUNDARY CONDITIONS [ may be overwritten in section BOUNDARY */
/* CONDITIONS, below ]: */
/*............................................................................*/
# if CNN_BDCFLOW == 1
/*............................................................................*/
# if CNN_NSLPMDE == 0 /* no-slip mode 0 */

                           jj = null; do
                           { 
                              ( hre->uf[mm][fc][jj] ) = ZERO;
                              ( hre->iu[mm][fc][jj] ) = ZERO;
/* adjacent face in neighbouring cell: */
                              ( hre->uf[mn][fn][jj] ) = ZERO;
                              ( hre->iu[mn][fn][jj] ) = ZERO;
                           } while (( ++jj ) < THREE );

# elif (( CNN_NSLPMDE == 1 )\
      ||( CNN_NSLPMDE == 2 ))

                           jj = null; do
                           { 
                              ( hre->uf[mm][fc][jj] ) = - ( hci->un[mm][jj] );
                              ( hre->iu[mm][fc][jj] ) = ZERO;
/* adjacent face in neighbouring cell: */			
                              ( hre->uf[mn][fn][jj] ) = - ( hci->un[mn][jj] );
                              ( hre->iu[mn][fn][jj] ) = ZERO;
                           } while (( ++jj ) < THREE );
# if CNN_NSLPMDE == 2
                           trnsvs(( hre->uf[mm][fc] ),
                              ( hsp->f[hm][fc] ), THREE );
                           trnsvs(( hre->uf[mn][fn] ),
                              ( hsp->f[hn][fn] ), THREE );
# endif /* CNN_NSLPMDE == 2 */
# endif /* CNN_NSLPMDE != 0 */
/*............................................................................*/
# endif /* CNN_BDCFLOW == 1 */
/*............................................................................*/
                        };
                     } /* end if cc == ( hsp->cnn[mn] ) */
                     else /* non-fluid neighbour or other connected component */
                     {
/*............................................................................*/
/* DEFAULT BOUNDARY CONDITIONS [ may be overwritten in section BOUNDARY */
/* CONDITIONS, below ]: */
/*............................................................................*/
# if CNN_BDCFLOW == 1
/*............................................................................*/
# if CNN_NSLPMDE == 0 /* no-slip mode 0 */

                        jj = null; do
                        { 
                           ( hre->uf[mm][fc][jj] ) = ZERO;
                           ( hre->iu[mm][fc][jj] ) = ZERO;
                        } while (( ++jj ) < THREE );

# elif (( CNN_NSLPMDE == 1 )\
      ||( CNN_NSLPMDE == 2 ))

                        jj = null; do
                        { 
                           ( hre->uf[mm][fc][jj] ) = - ( hci->un[mm][jj] );
                           ( hre->iu[mm][fc][jj] ) = ZERO;
                        } while (( ++jj ) < THREE );
# if CNN_NSLPMDE == 2
                        trnsvs(( hre->uf[mm][fc] ),
                           ( hsp->f[hm][fc] ), THREE );
# endif /* CNN_NSLPMDE == 2 */
# endif /* CNN_NSLPMDE != 0 */
/*............................................................................*/
# endif /* CNN_BDCFLOW == 1 */
/*............................................................................*/
                    ;} /* end if ( cc != ... ) */
                  } /* end if ( mm < mn ) */
                  else if ( mn <= null ) /* [ fc is a mesh boundary face ] */
                  {
/*............................................................................*/
/* DEFAULT BOUNDARY CONDITIONS [ may be overwritten in section BOUNDARY */
/* CONDITIONS, below ]: */
/*............................................................................*/
# if CMPRSSBL != 0
                     ( hre->rf[mm][fc] ) = ( hci->rn[mm] );
# endif /* CMPRSSBL != 0 */
/*............................................................................*/
# if CNN_BDCFLOW == 1
/*............................................................................*/
# if CNN_NSLPMDE == 0 /* no-slip mode 0 */

                     jj = null; do
                     { 
                        ( hre->uf[mm][fc][jj] ) = ZERO;
                        ( hre->iu[mm][fc][jj] ) = ZERO;
                     } while (( ++jj ) < THREE );

# elif (( CNN_NSLPMDE == 1 )\
      ||( CNN_NSLPMDE == 2 ))

                     jj = null; do
                     { 
                        ( hre->uf[mm][fc][jj] ) = - ( hci->un[mm][jj] );
                        ( hre->iu[mm][fc][jj] ) = ZERO;
                     } while (( ++jj ) < THREE );
# if CNN_NSLPMDE == 2
                     trnsvs(( hre->uf[mm][fc] ),
                        ( hsp->f[hm][fc] ), THREE );
# endif /* CNN_NSLPMDE == 2 */
# endif /* CNN_NSLPMDE != 0 */
/*............................................................................*/
# endif /* CNN_BDCFLOW == 1 */
/*............................................................................*/
                 ;}; /* end if ( mn <= null ) [ fc is mesh boundary face ] */
               } while(( ++fc ) < FACES );
            }; /* end if ( cc == ( hsp->cnn[mm] )) */
            mm++;
         }; /* end while( mm <= upfcl ) */
/*............................................................................*/
# if (( CMPRSSBL == 0 )\
    ||( BSQAPRX == 1 ))
/*............................................................................*/
/* SUCCESIVE OVERRELAXATION [ solving Poisson's equation for the pressure ] */
/* reflection part of SOR cycle */
/*............................................................................*/
         if ( CNN_SORREPT < ( ++( sorrept[cc] )))
         {
            sorrept[null] = ONE;
            sorrept[cc] = ONE;

            actbnd = dynbnd;
            itr = null; do /* while (( ++itr ) <= maxitrs )... */
            {
               df = ZERO;
               dv = ZERO;
               minpr = DBL_MAX;
               maxpr = - minpr;

               mm = lwfcl;
               while( mm <= upfcl )
               {
                  if ( cc == ( hsp->cnn[mm] ))
                  {
                     hm = ( hsp->hh[mm] );
/*............................................................................*/
/* divergence of flow: */

                     dru = ZERO;
                     kk = null; do
                     {
                        fc = null; do
                        {
                           ds = ( hre->uf[mm][fc][kk] );
                           uf[fc][kk] = ds;
/*............................................................................*/
# if CMPRSSBL != 0
                           ds *= ( hre->rf[mm][fc] );
# endif /* if CMPRSSBL != 0 */
/*............................................................................*/
                           dru += ( ds*( hsp->f[hm][fc][kk] ));
                        } while(( ++fc ) < FACES );
                     } while(( ++kk ) < THREE );
/*............................................................................*/
# if CMPRSSBL != 0
                     dru /= dt;
# else /* if CMPRSSBL == 0 */
                     dru /= ( hsp->ft[hm] );
# endif /* if CMPRSSBL ... */
/*............................................................................*/
/* PRESSURE GRADIENT: */
/*............................................................................*/
/* actual pressure in the cell: */

                     if ( itr == null )
                     {
                        ( hre->pn[mm] ) = ( hci->pn[mm] );

                        fc = null; do
                        {
                           ( hre->pf[mm][fc] ) = ( hci->pf[mm][fc] );
                        } while(( ++fc ) < FACES );
                     };
/*............................................................................*/
/* copy actual pressure */
/* ... in the node: */
                     ( prs->node ) = ( hre->pn[mm] );

/* ... and on the faces: */
                     fc = null; do
                     {
                        ( prs->face[fc] ) = ( hre->pf[mm][fc] );
                     } while(( ++fc ) < FACES );
/*............................................................................*/
/* actual nodal pressure gradient, */
/* in cell coordinates [ with respect to basis b[*] ] */
/* db[k] = < b[k] | grad(p) >, */
/* and in canonical coordinates */
/* dp[k] = < e[k] | grad(p) > = ( adj(B)^-1 )*db */

                     prs = gradnt( prs, hm );

/* actual pressure gradients on the faces: */

                     prs = gradfc( prs, hm );
/*............................................................................*/
/* | div(p) - div(ro*u)/dt | per cell volume: */

                     dw = ( prs->dvgr ) - dru;
                     ds = fabs( dw/( hsp->vol[hm] ));

                     if ( df < ds )
                     {
                        df = ds;
                        max = mm;
                     };

                     if ( pcrbnd < ds )
                     {
                        dv += ( dw*dw );
/*............................................................................*/
/* compute nodal pressure that compensates dru(u) */
/* [ still with divergence dru as computed above ] */

                        dw = -dru;
                        ds = ZERO;
                        sgn = -ONE;

                        fc = null; do
                        {
                           jj = fc/TWO;
                           kk = null; do
                           {
                              if( kk != jj )
                                 dw += ( sgn*( hsp->s[hm][fc][kk] )*\
                                     ( prs->db[kk] ));
                           } while (( ++kk ) < THREE );
                           dw += ( 2.*( prs->face[fc] )*\
                              ( hsp->s[hm][fc][jj] ));
                           ds += ( hsp->s[hm][fc][jj] );
                           sgn *= ( -ONE );
                        } while (( ++fc ) < FACES );
                        dw /= ( 2.*ds );
/*............................................................................*/
/* update pressure [ with given relaxation factor rlxprs ]*/
# ifdef CNN_RLXPRSS
                        ( prs->nupd ) = dw*rlxprs + ( 1.-rlxprs )*( prs->node );
# else
                        ( prs->nupd ) = dw;
# endif /* CNN_RLXPRSS ... */
/*............................................................................*/
/* regularize [ delimit ] pressure, if enabled: */
# ifdef CNN_LMTPRSS
                        dw = fabs(( prs->nupd )*3.4/CNN_LMTPRSS );
                        ( prs->nupd ) *= ( CNN_LMTPRSS*( 1.-exp( -dw ))/3.4 );
# endif
/*............................................................................*/
/* copy updated nodal pressure: */

                        ( hre->pn[mm] ) = ( prs->nupd );
                        ( prs->node ) = ( prs->nupd );
/*............................................................................*/
/* retain absolute minimum and maximum of pressure: */

                        if (( prs->node ) < minpr )
                        {
                           cminp = mm;
                           minpr = ( prs->node );
                        };

                        if ( maxpr < ( prs->node ))
                        {
                           cmaxp = mm;
                           maxpr = ( prs->node );
                        };
/*............................................................................*/
/* transfer new [ truncated ] face pressure gradients: */

                        sgn = ONE;
                        fc = null; do
                        {
                           kk = ( short )( fc/TWO );
                           dw = ( hsp->s[hm][fc][kk] )*( prs->node );

                           ++kk; kk%=THREE; /* kk := ( kk+ONE ) modulo THREE */
                           dw += ( sgn*( hsp->s[hm][fc][kk] )*\
                                       ( prs->db[kk] )/2.);
                           ++kk; kk%=THREE;
                           dw += ( sgn*( hsp->s[hm][fc][kk] )*\
                                       ( prs->db[kk] )/2.);

                           ( hre->pt[mm][fc] ) = dw;

                           sgn *= ( -ONE );
                        } while(( ++fc ) < FACES );
                     }; /* end if ( pcrbnd < ds ) */
                  }; /* end, fluid neighbour in connected compt cc */
                  mm++;
               }; /* end while( mm <= upfcl ) */
/*............................................................................*/
/* connection part of SOR cycle */
/* ensure continuity of pressure gradient at cell interfaces */
/*............................................................................*/
               mm = lwfcl;
               while( mm <= upfcl )
               {
                  if ( cc == ( hsp->cnn[mm] ))
                  {
                     hm = ( hsp->hh[mm] );
                     ( prs->node ) = ( hre->pn[mm] );

/* pressure on cell faces: */

                     fc = null; do
                     {
                        ( prs->face[fc] ) = ( hre->pf[mm][fc] );
                     } while(( ++fc ) < FACES );
/*............................................................................*/
/* actual nodal pressure gradient, */
/* in cell coordinates [ with respect to basis b[*] ] */
/* db[k] = < b[k] | grad(p) >, */
/* and in canonical coordinates */
/* dp[k] = < e[k] | grad(p) > = ( adj(B)^-1 )*db */

                     prs = gradnt( prs, hm );

                     fc = null; do
                     {
/*............................................................................*/
/* the adjacent mesh cell index */
# if DSC_ADJCLL == 0 /* [ neighbouring cells labelled with faces ] */
                        mn = ( tpt->m[mm][fc] );
# else /* if DSC_ADJCLL == 1: [ neighbouring cells labelled with ports ] */
                        pp = prt1[fc];
                        mn = ( tpt->m[mm][pp] );
# endif
/*............................................................................*/
                        if ( mm < mn )
                        { 
                           if ( cc == ( hsp->cnn[mn] ))
                           {
                              hn = ( hsp->hh[mn] );
/*............................................................................*/
/* adjacent face index: */
# if DSC_FCELBL == 0
                              fn = ( tpt->f[mm][fc] );
# elif DSC_FCELBL == 1
                              fn = ( tpt->f[mm][fc] ) - ONE;
# elif DSC_FCELBL == 2
                              fn = ( tpt->f[mm][fc] );

		              if ( fn < null )
		              {
		                 fn = - fn;
		                 sgn = - ONE;
                              }
                              else
                                 sgn = ONE;

                              fn -= ONE;
# endif /* DSC_FCELBL == 2 */
/*............................................................................*/
/* at the cell interface make sure: */
/* CONTINUITY OF PRESSURE [ GRADIENT ] AT CELL INTERFACES */
/*............................................................................*/
                              kk = ( short ) ( fc/TWO );
                              ll = ( short ) ( fn/TWO );
/*............................................................................*/
/* [ s[fc][] = f[fc]*( adj(B)^-1 ) ] */
/*............................................................................*/
/* the mesh cell form vector at face fc: */
                              ss1 = ( hsp->s[hm][fc][kk] );
/* ... and at the adjacent face fn: */
                              ss2 = ( hsp->s[hn][fn][ll] );
                              ss0 = ss1 + ss2;

                              if ( CNN_TRIVIAL < fabs( ss0 ))
                              {
/*............................................................................*/
/* [ pt[fc][] = .5*f[fc]*( adj(B)^-1 )*( dp~/db ); the truncated gradient ] */
/*............................................................................*/
/* the truncated gradient at face fc: */
                                 ss3 = ( hre->pt[mm][fc] );
/* ... and at the adjacent face fn: */
                                 ss4 = ( hre->pt[mn][fn] );
                                 ss5 = ( ss3 + ss4 )/ss0;

                                 ( hre->pf[mm][fc] ) = ss5;
                                 ( hre->pf[mn][fn] ) = ss5;
                              }
                              else /* vanishing pressure [ Dirichlet BDC ] */
                              {    
/*............................................................................*/
# if CNN_BDCPRSS == 0 /* vanishing pressure [ Dirichlet BDC ] */

                                 ( hre->pf[mm][fc] ) = ZERO;
                                 ( hre->pf[mn][fn] ) = ZERO;

# elif CNN_BDCPRSS == 1 /* vanishing pressure gradient [ von Neumann BDC ] */

                                 kk = ( short ) ( fc/TWO );
                                 ss1 = ( hsp->s[hm][fc][kk] );

                                 if ( CNN_TRIVIAL < fabs( ss1 ))
                                 { 
                                    sgn = ONE - TWO*( fc%TWO ); /* = (-1)^fc */

                                    ++kk; kk %= THREE;
                                    ss2 = ( hsp->s[hm][fc][kk] )*\
                                          ( prs->db[kk] );
                                    ++kk; kk %= THREE;
                                    ss2 += (( hsp->s[hm][fc][kk] )*\
                                            ( prs->db[kk] ));

                                    ( hre->pf[mm][fc] ) = ( hre->pn[mm] ) + \
                                       .5*sgn*ss2/ss1;
                                    ( hre->pf[mn][fn] ) = ( hre->pf[mm][fc] );
                                 }
                                 else
                                 {                              
                                    ( hre->pf[mm][fc] ) = ( hre->pn[mm] );
                                    ( hre->pf[mn][fn] ) = ( hre->pn[mn] );
                                 };
# endif /* CNN_BDCPRSS == 1 */
/*............................................................................*/
                              }; /* if ( fabs( ss0 ) <= CNN_TRIVIAL ) */
                           } /* end if ( cc == ( hsp->cnn[mn] )) */
                           else /* mn non-fluid node or in other cnnectd cmpt */
                           { 
/*............................................................................*/
# if CNN_BDCPRSS == 0 /* vanishing pressure [ Dirichlet BDC ] */

                              ( hre->pf[mm][fc] ) = ZERO;

# elif CNN_BDCPRSS == 1 /* vanishing pressure gradient [ von Neumann BDC ] */

                              kk = ( short ) ( fc/TWO );
                              ss1 = ( hsp->s[hm][fc][kk] );

                              if ( CNN_TRIVIAL < fabs( ss1 ))
                              { 
                                 sgn = ONE - TWO*( fc%TWO ); /* = (-1)^fc */

                                 ++kk; kk %= THREE;
                                 ss2 = ( hsp->s[hm][fc][kk] )*\
                                       ( prs->db[kk] );
                                 ++kk; kk %= THREE;
                                 ss2 += (( hsp->s[hm][fc][kk] )*\
                                         ( prs->db[kk] ));

                                 ( hre->pf[mm][fc] ) = ( hre->pn[mm] ) +\
                                    .5*sgn*ss2/ss1;
                              }
                              else
                                 ( hre->pf[mm][fc] ) = ( hre->pn[mm] );
# endif /* CNN_BDCPRSS == 1 */
/*............................................................................*/
                           }; /* end non-fluid neighbour or other connct comp */
                        } /* end if ( mm < mn ) */
                        else if ( mn <= null )
                        {                 
/*............................................................................*/
# if CNN_BDCPRSS == 0 /* vanishing pressure [ Dirichlet BDC ] */

                           ( hre->pf[mm][fc] ) = ZERO;

# elif CNN_BDCPRSS == 1 /* vanishing pressure gradient [ von Neumann BDC ] */

                           kk = ( short ) ( fc/TWO );
                           ss1 = ( hsp->s[hm][fc][kk] );

                           if ( CNN_TRIVIAL < fabs( ss1 ))
                           { 
                              sgn = ONE - TWO*( fc%TWO ); /* = (-1)^fc */

                              ++kk; kk %= THREE;
                              ss2 = ( hsp->s[hm][fc][kk] )*\
                                    ( prs->db[kk] );
                              ++kk; kk %= THREE;
                              ss2 += (( hsp->s[hm][fc][kk] )*\
                                      ( prs->db[kk] ));

                              ( hre->pf[mm][fc] ) = ( hre->pn[mm] )+\
                                 .5*sgn*ss2/ss1;
                           }
                           else
                              ( hre->pf[mm][fc] ) = ( hre->pn[mm] );
# endif /* CNN_BDCPRSS == 1 */
/*............................................................................*/
                        }; /* endif ( mn <= null ) */
                     } while(( ++fc ) < FACES );
                  }; /* end if ( cc == ( hsp->cnn[mm] )) */
                  mm++;
               }; /* end while( mm <= upfcl ) */
	       itr++;
            } while ((( itr < maxitrs )\
                    &&( actbnd < dv ))\
                   ||( itr < minitrs ));

            if ( dwdw < dv )
               dwdw = dv;

            if ( dfm < df )
            {
               dfm = df;
               cdfmx = max;
            };
/*............................................................................*/
# if CNN_SORUPDT != 0
/* less than minitrs iterations -> dynbnd is slightly reduced, */
/* more than maxitrs iterations -> dynbnd slightly increased: */

            if ( itr == minitrs )
            {
               if (( --updmrk ) == minlimit )
               {
                  dynbnd *= sordecr;
                  updmrk = null;
               };
            }
            else if ( itr == maxitrs )
            {
               if (( ++updmrk ) == maxlimit )
               {
                  dynbnd *= sorincr;
                  updmrk = null;
               };
            }
            else
            {
/* continuously strengthen breaking conditions: */

               dynbnd *= sorrdce;
               updmrk = null;
            };
# endif /* CNN_SORUPDT != 0 */
/*............................................................................*/
         }; /* end if ( CNN_SORREPT < ( ++( sorrept[cc] ))) */
/*............................................................................*/
# endif /* (( CMPRSSBL == 0 )\
          ||( BSQAPRX != 0 )) */
/*............................................................................*/
# if CNN_BNDCLIP != 0
/* [ FLOW CLIPPING at domain boundaries ] */
/*............................................................................*/

         mm = lwfcl;
         while( mm <= upfcl )
         {
            if ( cc == ( hsp->cnn[mm] ))
            {
               hm = ( hsp->hh[mm] );

               fc = null; do
               {
                  switch( hsp->fctype[mm][fc] )
                  {
                    default: /* null: adjacent cell in the same */
                    break;         /* fluid connected component */
/*............................................................................*/
# if SLV_OUTMODE == 1
                    case -ONE: /* outflow: reduce to face-normal component */
                     longtd(( hre->un[mm] ), ( hsp->f[hm][fc] ), THREE );
                     longtd(( hre->uf[mm][fc] ),
                        ( hsp->f[hm][fc] ), THREE );
                    break;
# endif /* SLV_OUTMODE == 1 */
/*............................................................................*/
# if CNN_BNDCLIP == 1
                    case TWO: /* no-slip: cut face-normal component */
                     trnsvs(( hre->un[mm] ), ( hsp->f[hm][fc] ), THREE );
                    break;

                    case THREE: /* free slip: cut face-normal component */
                     trnsvs(( hre->un[mm] ), ( hsp->f[hm][fc] ), THREE );
                    break;
# endif /* CNN_BNDCLIP == 1 */
/*............................................................................*/
                  }; /* end switch */
               } while(( ++fc ) < FACES );
            }; /* end if ( cc == ( hsp->cnn[mm] )) */
            mm++;
         }; /* end while( mm <= upfcl ) */
# endif /* CNN_BNDCLIP != 0 [ FLOW CLIPPING at domain boundaries ] */
/*............................................................................*/
      }; /* end while(( ++cc ) < ( hsp->cnn[null] )) */
/*............................................................................*/
/* LARGE EDDY FILTERING [ if enabled ]                                        */
/*............................................................................*/
# if CNN_LESFLTR != 0  /* Call LES filter */

/*............................................................................*/
      hre = lesfilter( state ); /* Large Eddy Filter */
      lesmrk = ( state->lesfltr );
/*............................................................................*/
# endif /* CNN_LESFLTR != 0 */
/*............................................................................*/
/* DIVERGENCE AND PRESSURE [ log files ]                                      */
/*............................................................................*/
      if ( null < sorrept[null] )
      {
         sorrept[null] = null;
/*............................................................................*/
/* divergence log file */
/*
*//*     format
*//*
*//*     fprintf(( state->dscstr ),
*//*        "| iteration " );
*//*     fprintf(( state->dscstr ),
*//*        "| time_[s]   " ); 
*//*     fprintf(( state->dscstr ),
*//*        "| SOR_steps " );
*//*     fprintf(( state->dscstr ),
*//*        "| sum_dw*dw  " );
*//*     fprintf(( state->dscstr ),
*//*        "| max_|df|   " );
*//*     fprintf(( state->dscstr ),
*//*        "| cell_no  |\n" );
*//*..........................................................................*/
/* store values: */

         fprintf(( state->dscstr ),
            "|%10ld | %.4e | %9ld | %.4e | %.4e |%9ld %c\n", ( state->fldnn ),
            ( state->hctme ), ( long ) itr, dwdw, dfm, cdfmx, lesmrk );

         fflush( state->dscstr );
/*............................................................................*/
/* pressure log file */
/*
*//*     format:
*//*
*//*     fprintf(( state->prsstr ),
*//*        "| iteration " );
*//*     fprintf(( state->prsstr ),
*//*        "| time_[s]   " ); 
*//*     fprintf(( state->prsstr ),
*//*        "| minimum_[Pa] " );
*//*     fprintf(( state->prsstr ),
*//*        "| cell_no  " );
*//*     fprintf(( state->prsstr ),
*//*        "| maximum_[Pa] " );
*//*     fprintf(( state->prsstr ),
*//*        "| cell_no  |\n" );
*//*..........................................................................*/
/* store values: */

         fprintf(( state->prsstr ),
            "|%10ld | %.4e | %+.5e |%9ld | %+.5e |%9ld %c\n", ( state->fldnn ),
            ( state->hctme ), minpr, cminp, maxpr, cmaxp, lesmrk );

         fflush( state->prsstr );
      }; /* end if ( null < sorrept[null] ) */
   }; /* end if ( null < ( state->fldnn )) [ end of SOR ] */
/*............................................................................*/
/* BOUNDARY CONDITIONS */
/*............................................................................*/
/* inflow: */

   nn = null;
   while( nn < ( bpt->nif ))
   {
      mm = ( bpt->mif[nn] );
      fc = ( bpt->fif[nn] );
/*............................................................................*/
# if CMPRSSBL != 0
      hm = ( hsp->hh[mm] );
      ( hre->rn[mm] ) = ( hsp->rm[hm] );
      ( hre->rf[mm][fc] ) = ( hsp->rm[hm] );
# endif
/*............................................................................*/
      jj = null; do
      {
         ss0 = ( bpt->uf[nn][jj] );
/*
         if (( state->fldnn ) < 10000 )
	    ss0 *= (( double )( state->fldnn )/10000. );
*/
         ( hre->un[mm][jj] ) = ss0;
         ( hre->uf[mm][fc][jj] ) = ss0;
         ( hre->iu[mm][fc][jj] ) = ZERO;
      } while(( ++jj ) < THREE );

      ss0 = ( bpt->ti[nn] ); /* inflow temperature */
      ( hre->tn[mm] ) = ss0;
      ( hre->tf[mm][fc] ) = ss0;
      nn++ ;
   }; /* end while( nn < ( bpt->nif )) */
/*............................................................................*/
/* outflow: */

   nn = null;
   while( nn < ( bpt->nof ))
   {
      mm = ( bpt->mof[nn] );
      fc = ( bpt->fof[nn] );

      jj = null; do
      {
/*............................................................................*/
# if SLV_OUTMODE == 0
         ss0 = ( hre->un[mm][jj] );
         ( hre->uf[mm][fc][jj] ) = ss0;
# elif SLV_OUTMODE == 1 
         ss0 = ( hre->uf[mm][fc][jj] );
         ( hre->un[mm][jj] ) = ss0;
# elif SLV_OUTMODE == 2 
         ss0 = ( bpt->vf[nn][jj] ); /* imposed velocity */
         ( hre->uf[mm][fc][jj] ) = ss0;
         ( hre->un[mm][jj] ) = ss0;
# endif /* SLV_OUTMODE == 2 */
/*............................................................................*/
         ( hre->iu[mm][fc][jj] ) = ZERO;
      } while(( ++jj ) < THREE );
      nn++ ;
   }; /* end while( nn < ( bpt->nof )) */
/*............................................................................*/
/* no-slip: */

   nn = null;
   while( nn < ( bpt->nns ))
   {
      mm = ( bpt->mns[nn] );
      fc = ( bpt->fns[nn] );
/*............................................................................*/
# if CNN_NSLPMDE == 0

      jj = null; do
      {
         ( hre->uf[mm][fc][jj] ) = ZERO;
         ( hre->iu[mm][fc][jj] ) = ZERO;
      } while(( ++jj ) < THREE );

# elif (( CNN_NSLPMDE == 1 )\
      ||( CNN_NSLPMDE == 2 ))

      hm = ( hsp->hh[mm] );
      
      jj = null; do
      {
         ( hre->uf[mm][fc][jj] ) = - ( hci->un[mm][jj] );
         ( hre->iu[mm][fc][jj] ) = ZERO;
      } while(( ++jj ) < THREE );
/*............................................................................*/
# if CNN_NSLPMDE == 2
      trnsvs(( hre->uf[mm][fc] ), ( hsp->f[hm][fc] ), THREE );
# endif
/*............................................................................*/
# endif /* CNN_NSLPMDE != 0 */
/*............................................................................*/
      nn++ ;
   }; /* end while( nn < ( bpt->nns )) */
/*............................................................................*/
/* free slip: */

   nn = null;
   while( nn < ( bpt->nsl ))
   {
      mm = ( bpt->msl[nn] );
      fc = ( bpt->fsl[nn] );
      hm = ( hsp->hh[mm] );

      jj = null; do
      {
         uf[fc][jj] = ( hre->un[mm][jj] );
/*............................................................................*/
# if CNN_FRSLMDE == 1
         vf[fc][jj] = uf[fc][jj];
# endif
/*............................................................................*/
      } while(( ++jj ) < THREE );
/*............................................................................*/
/* the face-tangential [(->f)-transversal] component of uf[]: */

      trnsvs( uf[fc], ( hsp->f[hm][fc] ), THREE );
/*............................................................................*/
# if CNN_FRSLMDE == 1 /* reverse normal fluid component */
/* the face-orthogonal [(->f)-longitudinal] component of vf[]: */

      longtd( vf[fc], ( hsp->f[hm][fc] ), THREE );
# endif
/*............................................................................*/
      jj = null; do
      {
/*............................................................................*/
# if CNN_FRSLMDE == 1
         ( hre->uf[mm][fc][jj] ) = uf[fc][jj] - vf[fc][jj];
# else
         ( hre->uf[mm][fc][jj] ) = uf[fc][jj];
# endif
/*............................................................................*/
         ( hre->iu[mm][fc][jj] ) = ZERO;
      } while(( ++jj ) < THREE );
      nn++ ;
   }; /* end while( nn < ( bpt->nsl )) */
/* end: fluid flow boundary conditions */
/*............................................................................*/

   return hre;
}
/*============================================================================*/
/*********************** end of function body conflow(*) **********************/
/************************* end of file conflow-1.0r3.a ************************/
