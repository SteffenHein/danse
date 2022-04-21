/* [ file: cncthcr-1.0r3.a ] */
/*******************************************************************************
*                                                                              *
*   Function body cncthcr(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   DSC diffusive heat currents connection map                                 *
*   [ convective heat currents and fluid flows are separately connected        *
*     in function conflow(*) ]                                                 *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

DSC_HCRRTS *\
cncthcr( struct solverstat *state )
{
/*----------------------------------------------------------------------------*/
/* allusions: */
/*     
   extern struct topology top;
   extern struct hcrsmx hcs;
   extern struct boundary bnd;

# if DSC_DOMAIN == 1
   extern DSC_FIELDS fld;
# else
   extern DSC_FIELDS fld[];
# endif
*/
/* declarations: */
/*............................................................................*/
/* structure pointers: */

   static struct topology
     *tpt = NULL;

   static struct tlmsmx
     *spt = NULL;

   static struct hcrsmx
     *hsp = NULL;

   static struct boundary
     *bpt = NULL;

# if CNN_SMTHSRC == 1
   static struct evaluation
     *vpt = NULL;
# endif

   static DSC_HCRRTS
     *hci = NULL,
     *hre = NULL;

   static DSC_FIELDS
     *inc = NULL,
     *out = NULL;
/*............................................................................*/
/* register type: */

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
      fn = null,
      fc = null,
      pp = null,
      pn = null;

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

/* and others: */

   static double
      crr[TWO] = {ZERO},
      cri[TWO] = {ZERO},
      hsk[TWO] = {ZERO};

   static double
      khm = ZERO,
      khn = ZERO,
      scdnwd = ZERO,
      scupwd = ZERO;

   static signed char
      sgn = null;
/*
   static char
      ptr[STS_SIZE] = {null};
*/
   static const char /* 1st and 2nd ports pertinent to faces 0,...,5: */
      prt1[FACES] = { 7, 5, 11, 9, 3, 1 }, /* 1st port */
      prt2[FACES] = { 10, 8, 2, 0, 6, 4 }; /* 2nd port */
/*
   faces pertinent to ports 0,1,...,11:
   static const char 
      fce[PORTS] = { 3, 5, 2, 4, 5, 1, 4, 0, 1, 3, 0, 2 };
*/
/*............................................................................*/
   static CLUSTER
     *tmp = NULL;

/* prototypes: */

   double pow( double x, double y );

   CLUSTER 
      *gradnt( CLUSTER *cls, long hh );
/*----------------------------------------------------------------------------*/
/* overtake solverstate: */

   tpt = ( state->tpt );
   spt = ( state->spt );
   hsp = ( state->hsp );
   bpt = ( state->bpt );

# if CNN_SMTHSRC == 1
   vpt = ( state->vpt ); 
# endif

   inc = ( state->out );
   out = ( state->inc );

   jj = ( state->hclbl );

   hci = ( state->hre[( int )jj] );
   hre = ( state->hci[( int )jj] );

   tmp = ( state->tmp );
/*............................................................................*/
/* here starts the job: */
/*............................................................................*/
   sgn = null; /* at least use it once [ to satisfy the compiler ] */
/*............................................................................*/
/* smoothed heat source updating coefficients */
/* [ for skin effect lossy boundary ]: */

# if CNN_SMTHSRC == 0
   scdnwd = 0.;
# elif CNN_SMTHSRC == 1
   scdnwd = 1. - ( 1./(( double )(( state->vpt->rc ) + ONE )));
# else 
   scdnwd = 1. - ( 1./(( double ) CNN_SMTHSRC ));
# endif
   scupwd = 1. - scdnwd;
/*............................................................................*/
   if(( state->dmn ) == 't' ) /* time domain process */
   {
/*............................................................................*/
# if DSC_DOMAIN == 2
      printf( "\n Error message from function %s :", __func__ );
      printf( "\n\n Solver is compiled in frequency "
         "domain mode 'DSC_DOMAIN = 2'" );
      printf( "\n - but run in time domain mode 1 !!!" );
      printf( "\n\n [ Change time/frequency domain macro DSC_DOMAIN "
         "in SOLVER.CONF from 2 to" );
      printf( "\n   1: time domain, or"
         "\n   0: time & frequency domain." );
      printf( "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[solver_compiled_"
         "in_mode_DSC_DOMAIN=2] ", LGS_SIZE );
      return NULL;
/*............................................................................*/
# else /* DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
/* CONNECTION CYCLE */
/*............................................................................*/
      mm = null;
      while(( mm++ ) < ( tpt->n ))
      {
         hm = ( hsp->hh[mm] ); /* s-parameter index pertinent to cell mm */
         if (( hsp->ttyp[hm] ) != 't' ) /* non-trivial thermal node */
         {
            khm = ( hsp->kh[hm] );
/*............................................................................*/
# if DSC_FLDMDE != 0 
# if NUSSELT != 0
            if ((( hsp->ttyp[hm] ) == 'f' ) /* 'f'luid type cell */
              &&( CNN_TRIVIAL < ( hci->nus[mm] )))
            {
               khm *= ( hci->nus[mm] );
            };
# endif /* NUSSELT != 0 */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
/* copy actual face temperatures: */

            fc = null; do
            {
               ( tmp->face[fc] ) = ( hci->tf[mm][fc] );
            } while(( ++fc ) < FACES );
/*............................................................................*/
/* actual node temperature gradient, */
/* in cell coordinates [ with respect to basis b[] ] */
/* dtdb[k] = < b[k] | grad(T) >, */
/* and in canonical coordinates */
/* dt[k] = < e[k] | grad(T) > = adj(B^-1 )*dtdb */

            tmp = gradnt( tmp, hm );

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
               if ( mm < mn ) /* [ in particular, fc is an interface ] */ 
               {           
                  hn = ( hsp->hh[mn] ); /* the pertinent s-parameter index */
                  if (( hsp->ttyp[hn] ) != 't' ) /* non-trivial thermal */
                  {
                     khn = ( hsp->kh[hn] );
/*............................................................................*/
# if DSC_FLDMDE != 0 
# if NUSSELT != 0
                     if ((( hsp->ttyp[hn] ) == 'f' ) /* 'f'luid type cell */
                       &&( CNN_TRIVIAL < ( hci->nus[mn] )))
                     {
                        khn *= ( hci->nus[mn] );
                     };
# endif /* NUSSELT != 0 */
# endif /* DSC_FLDMDE != 0 */
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
/*............................................................................*/
# if CNN_CHKFCE == 1
                     if ( sgn != ( TWO*(( fc+fn ) % TWO ) - ONE ))
                     {
                        printf( "\n Error message from function %s :",
                           __func__ );
                        printf( "\n\n Face orientation inconsistent with "
                           "sense of ( tpt->f[%ld][%d] )", mm, fc );
                        printf( "\n [ as defined in elsy linker(*) "
                           "function ]." );
                        printf( "\n - Please check DSC mesh topology.\n " );

                        exit( EXIT_FAILURE );
                     };
# endif /* CNN_CHKFCE == 1 */
/*............................................................................*/
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
/* CONTINUITY OF HEAT CURRENT at cell interfaces: */

                     ss3 = ss1*khm;
                     ss4 = ss2*khn;
                     ss0 = ss3 + ss4;

                     if ( CNN_TRIVIAL < fabs( ss0 ))
                     {
/*............................................................................*/
/* [ rc[fc][] = .5*f[fc]*( adj(B)^-1 )*( dT~/db ); the truncated gradient ] */
/*............................................................................*/
/* the truncated temperature gradient at face fc: */
                        ss5 = ( hci->tt[mm][fc] ); 
/* ... and at the adjacent face fn: */
                        ss6 = ( hci->tt[mn][fn] );
/* ss7 = tf[fc] = tf[fn] */
                        ss7 = ( ss5 + ss6 )/ss0;
/* ss8 = < fc | grad( hc(fc)*tf ) > = - < fn | grad( hc(fn)*tf ) > */
                        ss8 = ss6 - ss5 + ss7*( ss3-ss4 );

                        ( hre->tf[mm][fc] ) = ss7;
                        ( hre->tf[mn][fn] ) = ss7;
                        ( hre->ic[mm][fc] ) = ss8;
                        ( hre->ic[mn][fn] ) = - ss8;
                     }
                     else /* von Neumann boundary conditions [ default ] */
                     {    
                        kk = ( short ) ( fc/TWO );
                        ss1 = ( hsp->s[hm][fc][kk] );

                        if ( CNN_TRIVIAL < fabs( ss1 ))
                        { 
                           sgn = ONE - TWO*( fc%TWO ); /* = ( -1 )^fc */

                           ++kk; kk %= THREE;
                           ss2 = ( hsp->s[hm][fc][kk] )*( tmp->db[kk] );
                           ++kk; kk %= THREE;
                           ss2 += (( hsp->s[hm][fc][kk] )*( tmp->db[kk] ));

                           ( hre->tf[mm][fc] ) = ( hre->tn[mm] )+\
                              .5*sgn*ss2/ss1;
                           ( hre->tf[mn][fn] ) = ( hre->tf[mm][fc] );
                        }
                        else
                        {                              
                           ( hre->tf[mm][fc] ) = ( hre->tn[mm] );
                           ( hre->tf[mn][fn] ) = ( hre->tn[mn] );
                        };

                        ( hre->ic[mm][fc] ) = ZERO;
                        ( hre->ic[mn][fn] ) = ZERO;
                     };
	          } /* end if (( hsp->ttyp[hn] ) != 't' ), thermal non-trivial*/
                  else /* non-thermal neighbouring cell -- no heat exchange: */
		  {    /* von Neumann boundary conditions [ default ] */
                     kk = ( short ) ( fc/TWO );
                     ss1 = ( hsp->s[hm][fc][kk] );

                     if ( CNN_TRIVIAL < fabs( ss1 ))
                     { 
                        sgn = ONE - TWO*( fc%TWO ); /* = ( -1 )^fc */

                        ++kk; kk %= THREE;
                        ss2 = ( hsp->s[hm][fc][kk] )*( tmp->db[kk] );
                        ++kk; kk %= THREE;
                        ss2 += (( hsp->s[hm][fc][kk] )*( tmp->db[kk] ));

                        ( hre->tf[mm][fc] ) = ( hre->tn[mm] )+\
                           .5*sgn*ss2/ss1;
                     }
                     else
                        ( hre->tf[mm][fc] ) = ( hre->tn[mm] );

                     ( hre->ic[mm][fc] ) = ZERO;
                  }; /* end, non-thermal neighbouring cell */
               } /* end if ( mm < mn ) */
               else if ( mn <= null ) /* no neighbour -- no heat exchange */
               {    /* von Neumann boundary conditions [ default ] */
                  kk = ( short ) ( fc/TWO );
                  ss1 = ( hsp->s[hm][fc][kk] );

                  if ( CNN_TRIVIAL < fabs( ss1 ))
                  { 
                     sgn = ONE - TWO*( fc%TWO ); /* = ( -1 )^fc */

                     ++kk; kk %= THREE;
                     ss2 = ( hsp->s[hm][fc][kk] )*( tmp->db[kk] );
                     ++kk; kk %= THREE;
                     ss2 += (( hsp->s[hm][fc][kk] )*( tmp->db[kk] ));

                     ( hre->tf[mm][fc] ) = ( hre->tn[mm] ) + \
                        .5*sgn*ss2/ss1;
                  }
                  else
                     ( hre->tf[mm][fc] ) = ( hre->tn[mm] );

                  ( hre->ic[mm][fc] ) = ZERO;
               }; /* end if ( mn <= null ) [ no neighbour ] */
/*............................................................................*/
            } while(( ++fc ) < FACES );
         }; /* end if (( hsp->ttyp[hm] ) != 't'rivial ) */
      }; /* end while(( mm++ ) < ( tpt->n )) */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
   } /* end if (( state->dmn ) == 't'ime domain ) */
   else /* if (( state->dmn ) != 't' ): frequency domain */
   {
/*............................................................................*/
# if DSC_DOMAIN == 1 
      printf( "\n Error message from function %s :", __func__ );
      printf( "\n\n Solver is compiled in time "
	 " domain mode 'DSC_DOMAIN = 1'" );
      printf( "\n - but run in frequency domain mode 2 !!!" );
      printf( "\n\n [ Change time/frequency domain macro DSC_DOMAIN "
	 "in SOLVER.CONF from 1 to" );
      printf( "\n   0: time & frequency domain, or" );
      printf( "\n   2: frequency domain." );
      printf( "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[solver_compiled_"
         "in_mode_DSC_DOMAIN=1] ", LGS_SIZE );
      return NULL;
/*............................................................................*/
# else /* the remaining cases: DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
      inc = ( state->inc );
      out = ( state->out );
/*............................................................................*/
/* CONNECTION CYCLE */
/*............................................................................*/
      mm = null;
      while(( mm++ ) < ( tpt->n ))
      {
         hm = ( hsp->hh[mm] ); /* s-parameter index pertinent to cell mm */
         if (( hsp->ttyp[hm] ) != 't' ) /* non-trivial thermal node */
         {
            khm = ( hsp->kh[hm] );
/*............................................................................*/
# if DSC_FLDMDE != 0 
# if NUSSELT != 0
            if ((( hsp->ttyp[hm] ) == 'f' ) /* 'f'luid type cell */
              &&( CNN_TRIVIAL < ( hci->nus[mm] )))
            {
               khm *= ( hci->nus[mm] );
            };
# endif /* NUSSELT != 0 */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
/* copy actual face temperatures: */

            fc = null; do
            {
               ( tmp->face[fc] ) = ( hci->tf[mm][fc] );
            } while(( ++fc ) < FACES );
/*............................................................................*/
/* actual node temperature gradient, */
/* in cell coordinates [ with respect to basis b[] ] */
/* dtdb[k] = < b[k] | grad(T) >, */
/* and in canonical coordinates */
/* dt[k] = < e[k] | grad(T) > = adj(B^-1 )*dtdb */

            tmp = gradnt( tmp, hm );

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
               if ( mm < mn ) /* [ in particular, fc is an interface ] */ 
               {           
                  hn = ( hsp->hh[mn] ); /* the pertinent s-parameter index */
                  if (( hsp->ttyp[hn] ) != 't' ) /* non-trivial thermal */
                  {
                     khn = ( hsp->kh[hn] );
/*............................................................................*/
# if DSC_FLDMDE != 0 
# if NUSSELT != 0
                     if ((( hsp->ttyp[hn] ) == 'f' ) /* 'f'luid type cell */
                       &&( CNN_TRIVIAL < ( hci->nus[mn] )))
                     {
                        khn *= ( hci->nus[mn] );
                     };
# endif /* NUSSELT != 0 */
# endif /* DSC_FLDMDE != 0 */
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
/*............................................................................*/
# if CNN_CHKFCE == 1
                     if ( sgn != ( TWO*(( fc+fn ) % TWO ) - ONE ))
                     {
                        printf( "\n Error message from function %s :",
                           __func__ );
                        printf( "\n\n Face orientation inconsistent with "
                           "sense of ( tpt->f[%ld][%d] )", mm, fc );
                        printf( "\n [ as defined in elsy linker(*) "
                           "function ]." );
                        printf( "\n - Please check DSC mesh topology.\n " );

                        exit( EXIT_FAILURE );
                     };
# endif /* CNN_CHKFCE == 1 */
/*............................................................................*/
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
/* CONTINUITY OF HEAT CURRENT at cell interfaces: */

                     ss3 = ss1*khm;
                     ss4 = ss2*khn;
                     ss0 = ss3 + ss4;

                     if ( CNN_TRIVIAL < fabs( ss0 ))
                     {
/*............................................................................*/
/* [ rc[fc][] = .5*f[fc]*( adj(B)^-1 )*( dT~/db ); the truncated gradient ] */
/*............................................................................*/
/* the truncated temperature gradient at face fc: */
                        ss5 = ( hci->tt[mm][fc] ); 
/* ... and at the adjacent face fn: */
                        ss6 = ( hci->tt[mn][fn] );
/* ss7 = tf[fc] = tf[fn] */
                        ss7 = ( ss5 + ss6 )/ss0;
/* ss8 = < fc | grad( hc(fc)*tf ) > = - < fn | grad( hc(fn)*tf ) > */
                        ss8 = ss6 - ss5 + ss7*( ss3-ss4 );

                        ( hre->tf[mm][fc] ) = ss7;
                        ( hre->tf[mn][fn] ) = ss7;
                        ( hre->ic[mm][fc] ) = ss8;
                        ( hre->ic[mn][fn] ) = - ss8;
                     }
                     else /* von Neumann boundary conditions [ default ] */
                     {    
                        kk = ( short ) ( fc/TWO );
                        ss1 = ( hsp->s[hm][fc][kk] );

                        if ( CNN_TRIVIAL < fabs( ss1 ))
                        { 
                           sgn = ONE - TWO*( fc%TWO ); /* = ( -1 )^fc */

                           ++kk; kk %= THREE;
                           ss2 = ( hsp->s[hm][fc][kk] )*( tmp->db[kk] );
                           ++kk; kk %= THREE;
                           ss2 += (( hsp->s[hm][fc][kk] )*( tmp->db[kk] ));

                           ( hre->tf[mm][fc] ) = ( hre->tn[mm] )+\
                              .5*sgn*ss2/ss1;
                           ( hre->tf[mn][fn] ) = ( hre->tf[mm][fc] );
                        }
                        else
                        {                              
                           ( hre->tf[mm][fc] ) = ( hre->tn[mm] );
                           ( hre->tf[mn][fn] ) = ( hre->tn[mn] );
                        };

                        ( hre->ic[mm][fc] ) = ZERO;
                        ( hre->ic[mn][fn] ) = ZERO;
                     };
	          } /* end if (( hsp->ttyp[hn] ) != 't' ), thermal non-trivial*/
                  else /* non-thermal neighbouring cell -- no heat exchange: */
		  {    /* von Neumann boundary conditions [ default ] */
                     kk = ( short ) ( fc/TWO );
                     ss1 = ( hsp->s[hm][fc][kk] );

                     if ( CNN_TRIVIAL < fabs( ss1 ))
                     { 
                        sgn = ONE - TWO*( fc%TWO ); /* = ( -1 )^fc */

                        ++kk; kk %= THREE;
                        ss2 = ( hsp->s[hm][fc][kk] )*( tmp->db[kk] );
                        ++kk; kk %= THREE;
                        ss2 += (( hsp->s[hm][fc][kk] )*( tmp->db[kk] ));

                        ( hre->tf[mm][fc] ) = ( hre->tn[mm] )+\
                           .5*sgn*ss2/ss1;
                     }
                     else
                        ( hre->tf[mm][fc] ) = ( hre->tn[mm] );

                     ( hre->ic[mm][fc] ) = ZERO;
                  }; /* end, non-thermal neighbouring cell */
               } /* end if ( mm < mn ) */
               else if ( mn <= null ) /* no neighbour -- no heat exchange */
               {    /* von Neumann boundary conditions [ default ] */
                  kk = ( short ) ( fc/TWO );
                  ss1 = ( hsp->s[hm][fc][kk] );

                  if ( CNN_TRIVIAL < fabs( ss1 ))
                  { 
                     sgn = ONE - TWO*( fc%TWO ); /* = ( -1 )^fc */

                     ++kk; kk %= THREE;
                     ss2 = ( hsp->s[hm][fc][kk] )*( tmp->db[kk] );
                     ++kk; kk %= THREE;
                     ss2 += (( hsp->s[hm][fc][kk] )*( tmp->db[kk] ));

                     ( hre->tf[mm][fc] ) = ( hre->tn[mm] ) + \
                        .5*sgn*ss2/ss1;
                  }
                  else
                     ( hre->tf[mm][fc] ) = ( hre->tn[mm] );

                  ( hre->ic[mm][fc] ) = ZERO;
               }; /* end if ( mn <= null ) [ no neighbour ] */
            } while(( ++fc ) < FACES );
         }; /* end if (( hsp->ttyp[hm] ) != 't'rivial ) */
      }; /* end while(( mm++ ) < ( tpt->n )) */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
   }; /* end if (( state->dmn ) != 't' ) */
/*............................................................................*/
/* BOUNDARY CONDITIONS */
/*............................................................................*/
/* [ skin effect ] heat source boundary faces: */

   if (( state->sws ) == ONE ) /* electric sources switched ON */
   {
      nn = null;
      while( nn < ( bpt->nsk ))
      {
         mm = ( bpt->msk[nn] );
         hm = ( hsp->hh[mm] );

         if (( hsp->ttyp[hm] ) != 't' ) /* cell mm: non trivial thermal */
         {
            fc = ( bpt->fsk[nn] );
            hn = ( spt->hh[mm] );

            if ((( spt->etyp[hn] ) != 't' ) /* cell mm: non-trivial electric */
              &&(( spt->mtyp[hn] ) != 't' ))
            {
               mn = mm;
               fn = fc;
            }
            else /* cell mm: trivial electric: evaluate neighbouring cell */
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
               if ( null < mn )
               {
/*............................................................................*/
/* adjacent face and/or ports: */
# if DSC_FCELBL == 0
                  fn = ( tpt->f[mm][fc] );
# elif DSC_FCELBL == 1
                  fn = ( tpt->f[mm][fc] ) - ONE;
# elif DSC_FCELBL == 2
                  fn = ( tpt->f[mm][fc] );

                  if ( fn < null )
                     fn = - fn;

                  fn -= ONE;
# endif
/*............................................................................*/
               }; /* end if ( null < mn ) */
            }; /* end if trivial electric cell mm */

            pp = prt1[fn];
            pn = prt2[fn];

            if (( state->dmn ) == 't' ) /* time domain */
            {
               crr[null] = ( inc->r[mn][pp] ) - ( out->i[mn][pp] );
               crr[ONE] = ( inc->r[mn][pn] ) - ( out->i[mn][pn] );

               ss0 = ZERO;
               kk = null; do
               {
                  hsk[kk] = ZERO;
                  ll = null; do
                  {
                     hsk[kk] += (( bpt->sk[nn][kk][ll] )*crr[ll] );
                  } while(( ++ll ) < TWO );
                  hsk[kk] *= hsk[kk];
                  ss0 += hsk[kk];
               } while(( ++kk ) < TWO );
            }
            else /* frequency domain */
            {
               crr[null] = ( inc->r[mn][pp] ) - ( out->r[mn][pp] );
               crr[ONE] = ( inc->r[mn][pn] ) - ( out->r[mn][pn] );
               cri[null] = ( inc->i[mn][pp] ) - ( out->i[mn][pp] );
               cri[ONE] = ( inc->i[mn][pn] ) - ( out->i[mn][pn] );

               ss0 = ZERO;
               kk = null; do
               {
                  hsk[kk] = ZERO;
                  ll = null; do
                  {
                     hsk[kk] += (( bpt->sk[nn][kk][ll] )*crr[ll] );
                  } while(( ++ll ) < TWO );
                  hsk[kk] *= hsk[kk];
                  ss0 += hsk[kk];
               } while(( ++kk ) < TWO );

               kk = null; do
               {
                  hsk[kk] = ZERO;
                  ll = null; do
                  {
                     hsk[kk] += (( bpt->sk[nn][kk][ll] )*cri[ll] );
                  } while(( ++ll ) < TWO );
                  hsk[kk] *= hsk[kk];
                  ss0 += hsk[kk];
               } while(( ++kk ) < TWO );
               ss0 /= 2.;
            }; /* end frequency domain */

            ( bpt->hs[nn] ) *= scdnwd;
            ( bpt->hs[nn] ) += ( ss0*scupwd );
            ( hre->ic[mm][fc] ) += ( bpt->hs[nn] ); 

         }; /* end if (( hsp->ttyp[hm] ) != 't' ): non trivial thermal cell*/

         nn++ ;
      }; /* end while( nn < ( bpt->nsk )) */
   }; /* end if (( state->sws ) == ONE ) [ electric sources switched ON ] */
/*............................................................................*/
/* [ imposed ] heat current boundary conditions: */

   nn = null;
   while( nn < ( bpt->nhc ))
   {
      mm = ( bpt->mhc[nn] );
      fc = ( bpt->fhc[nn] );

      if ( CNN_TRIVIAL < fabs ( bpt->hc[nn] ))
      {
         ( hre->ic[mm][fc] ) += ( bpt->hc[nn] );
      };

      if ( CNN_TRIVIAL < fabs ( bpt->rd[nn] ))
      {
         ss0 = pow((( hci->tn[mm] ) + CELSIUS_TO_KELVIN ), 4. );
         ( hre->ic[mm][fc] ) += ( ss0*( bpt->rd[nn] ));
      };

      nn++ ;
   }; /* end while( nn < ( bpt->nhc )) */
/*............................................................................*/
/* surface heat conductance for a reference temperature */
/* [ of environment, e.g.] */

   nn = null;
   while( nn < ( bpt->nsc ))
   {
      mm = ( bpt->msc[nn] );
      fc = ( bpt->fsc[nn] );

      if ( CNN_TRIVIAL < fabs ( bpt->sc[nn] ))
      {
         ss0 = ( bpt->tr[nn] ) - ( hci->tf[mm][fc] );
         ( hre->ic[mm][fc] ) += ( ss0*( bpt->sc[nn] ));
      };
      nn++ ;
   }; /* end while( nn < ( bpt->nsc )) */
/*............................................................................*/
/* fixed face temperature boundary conditions: */

   nn = null;
   while( nn < ( bpt->ntf ))
   {
      mm = ( bpt->mtf[nn] );
      fc = ( bpt->ftf[nn] );
      hm = ( hsp->hh[mm] );

      ( hre->tf[mm][fc] ) = ( bpt->tf[nn] );

      kk = ( short ) ( fc/TWO );
      sgn = ( TWO*( fc%TWO ) - ONE );

      khm = ( hsp->kh[hm] );
/*............................................................................*/
# if DSC_FLDMDE != 0
# if NUSSELT != 0
      if ((( hsp->ttyp[hm] ) == 'f' ) /* 'f'luid type cell */
        &&( CNN_TRIVIAL < ( hci->nus[mm] )))
      {
         khm *= ( hci->nus[mm] );
      };
# endif /* NUSSELT != 0 */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      ss0 = 2.*khm*( hsp->s[hm][fc][kk] )*\
   	       (( bpt->tf[nn] ) - ( hci->tn[mm] ));

      ll = null; do
      {
         ++kk;
         kk %= THREE;
         ss0 += ( sgn*khm*( hsp->s[hm][fc][kk] )*\
                   (( hre->tf[mm][2*kk+ONE] ) - ( hre->tf[mm][2*kk] )));
      } while(( ++ll ) < TWO );

      ( hre->ic[mm][fc] ) += ss0;

      nn++ ;
   }; /* end while( nn < ( bpt->ntf )) */
/*............................................................................*/
/* fixed node temperature [ internal conditions ]: */

   nn = null;
   while( nn < ( bpt->ntn ))
   {
      mm = ( bpt->mtn[nn] );
      ( hre->tn[mm] ) = ( bpt->tn[nn] );

      nn++ ;
   }; /* end while( nn < ( bpt->ntn )) */
/*............................................................................*/

   return hre;
}
/*============================================================================*/
/*********************** end of file cncthcr-1.0r3.a **************************/
