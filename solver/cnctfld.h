/* [ file: cnctfld.h ] */
/*******************************************************************************
*                                                                              *
*   Function body cnctfld(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC Maxwell field connection map                                           *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
/* macros: */
/*----------------------------------------------------------------------------*/
# define PERDC_TD( ) \
{ \
   uu =   ( bpt->r ) * ( inc->r[nn1][pn] ) - ( bpt->i ) * ( inc->r[nn3][pn] ); \
   vv =   ( bpt->i ) * ( inc->r[nn1][pn] ) + ( bpt->r ) * ( inc->r[nn3][pn] ); \
   ( out->i[nn2][pp] ) = sgn * uu; \
   ( out->i[nn4][pp] ) = sgn * vv; \
   uu =   ( bpt->r ) * ( inc->r[nn2][pp] ) + ( bpt->i ) * ( inc->r[nn4][pp] ); \
   vv = - ( bpt->i ) * ( inc->r[nn2][pp] ) + ( bpt->r ) * ( inc->r[nn4][pp] ); \
   ( out->i[nn3][pn] ) = sgn * uu; \
   ( out->i[nn1][pn] ) = sgn * vv; \
}
/*............................................................................*/
# define PERDC_FD( ) \
{ \
   uu =   ( bpt->r ) * ( inc->r[nn1][pn] ) - ( bpt->i ) * ( inc->i[nn1][pn] ); \
   vv =   ( bpt->i ) * ( inc->r[nn1][pn] ) + ( bpt->r ) * ( inc->i[nn1][pn] ); \
   ( out->r[nn2][pp] ) = sgn * uu; \
   ( out->i[nn2][pp] ) = cjg * sgn * vv; \
   uu =   ( bpt->r ) * ( inc->r[nn2][pp] ) + ( bpt->i ) * ( inc->i[nn2][pp] ); \
   vv = - ( bpt->i ) * ( inc->r[nn2][pp] ) + ( bpt->r ) * ( inc->i[nn2][pp] ); \
   ( out->r[nn1][pp] ) = sgn * uu; \
   ( out->i[nn1][pp] ) = cjg * sgn * vv; \
}
/*............................................................................*/
# define CNN_UNKNWN(AA,BB) \
{ \
   printf( "\n\n Error message from function %s : ", __func__ ); \
   printf( "\n Unknown/illegal macro option %s=%d ", (AA), (BB)); \
   printf( "\n [ Recompile program with a legal option ]\n\n " ); \
   return NULL; \
}
/*============================================================================*/

DSC_FIELDS * \
cnctfld( struct solverstat *state )
{
/* allusions: */
/*     
   extern struct topology top;
   extern struct tlmsmx smx;
   extern struct boundary bnd;

# if DSC_DOMAIN == 1
   extern DSC_FIELDS fld;
# else
   extern DSC_FIELDS fld[];
# endif
*/
/* declarations: */
/*............................................................................*/
/* register type: */

   register double
      uu = ZERO,
      vv = ZERO;

   register long
      hm = null,
      hn = null,
      mm = null,
      mn = null,
      nn = null,
      nn1 = null,
      nn2 = null;

   register signed char 
      fc = null,
      pn = null,
      pp = null,
      sgn = null;
/*
   static char
      ptr[STS_SIZE] = {null};
*/
/*............................................................................*/
# if DSC_DOMAIN != 1
   register signed char
      cjg = null;
# endif
/*............................................................................*/
# if DSC_DOMAIN != 2
   register long
     np = null,
     nn3 = null,
     nn4 = null;

   register signed char 
    cnt = null, 
    prd = null;
# endif
/*............................................................................*/
/* other types: */

   static const char /* 1st and 2nd ports pertinent to faces 0,...,5: */
      prt1[FACES] = { 7, 5, 11, 9, 3, 1 }, /* 1st port */
      prt2[FACES] = { 10, 8, 2, 0, 6, 4 }; /* 2nd port */
/*............................................................................*/
/* structure pointers: */

   static struct topology
     *tpt = NULL;

   static struct tlmsmx
     *spt = NULL;

   static struct boundary
     *bpt = NULL;

   static DSC_FIELDS
     *inc = NULL,
     *out = NULL;
/*----------------------------------------------------------------------------*/
/* overtake solverstate: */

   tpt = ( state->tpt );
   spt = ( state->spt );
   bpt = ( state->bpt );

   inc = ( state->out );
   out = ( state->inc );
/*............................................................................*/
/* here starts the job: */
/*............................................................................*/
   if (( state->dmn ) == 't' ) /* time domain process */
   {
/*............................................................................*/
# if DSC_DOMAIN == 2
      printf( "\n Error message from function %s :", __func__ );
      printf( "\n\n Program 'elfe' is compiled in frequency "
         "domain mode 'DSC_DOMAIN = 2'" );
      printf( "\n - but run in time domain mode 1 !!!" );
      printf( "\n\n [ Change time/frequency domain macro DSC_DOMAIN "
         "in program 'elfe' from 2 to" );
      printf( "\n   1: time domain, or"
         "\n   0: time & frequency domain." );
      printf( "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
         "in_mode_DSC_DOMAIN=2] ", LGS_SIZE );
      return NULL;
/*............................................................................*/
# else /* DSC_DOMAIN == 0 or 1 */

      prd = ( null < ( bpt->p ));

      np = null;
      cnt = null; do /* while cnt <= prd */
      {
         mm = null;
         while(( mm++ ) < ( tpt->n ))
         {
            hm = ( spt->hh[mm] );
            if ((( spt->etyp[hm] ) != 't' )
              &&(( spt->mtyp[hm] ) != 't' ))
            {
               nn1 = mm + np;
/*............................................................................*/
# if DSC_ADJCLL == 0

               fc = null; do
               {
                  mn = ( tpt->m[mm][fc] ); /* the adjacent mesh cell index */
                  if ( null < mn )
                  {
                     hn = ( spt->hh[mn] );

                     if ((( spt->etyp[hn] ) != 't' )
                       &&(( spt->mtyp[hn] ) != 't' ))
                     {
                        nn2 = mn + np;
/*............................................................................*/
/* 1st adjacent port: */
                        pp = prt1[fc];
                        pn = ( tpt->p[mm][pp] );
                        if ( pn < null )
                        {
                           pn = - pn;
                           sgn = - ONE;
                        }
                        else
                           sgn = ONE;

                        pn -= ONE;

                        ( out->i[nn1][pp] ) = sgn*( inc->r[nn2][pn] );
/*............................................................................*/
# if VICE_VERSA == 1
                        ( out->i[nn2][pn] ) = sgn*( inc->r[nn1][pp] );
# endif
/*............................................................................*/
/* 2nd adjacent port index: */
                        pp = prt2[fc];
                        pn = ( tpt->p[mm][pp] );

                        if ( pn < null )
                        {
                           pn = - pn;
                           sgn = - ONE;
                        }
                        else
                           sgn = ONE;

                        pn -= ONE;

                        ( out->i[nn1][pp] ) = sgn * ( inc->r[nn2][pn] );
/*............................................................................*/
# if VICE_VERSA == 1
                        ( out->i[nn2][pn] ) = sgn * ( inc->r[nn1][pp] );
# endif
/*............................................................................*/
                     }; /* end if (( spt->etyp[hn] ) != 't' ) */
                        /* [ non trivial cell ] */
                  } 
		  else /* if ( mn <= null ) */
		  { 
                     pp = prt1[fc];
                     pn = prt2[fc];

                     switch ( mn )
                     {
                       case ELECTRIC_WALL:
                        ( out->i[nn1][pp] ) = CNN_ERFL*\
                           ( inc->r[nn1][pp] );
                        ( out->i[nn1][pn] ) = CNN_ERFL*\
                           ( inc->r[nn1][pn] );
                        break;
                       
                       case MAGNETIC_WALL:
                        ( out->i[nn1][pp] ) = CNN_HRFL*\
                           ( inc->r[nn1][pp] );
                        ( out->i[nn1][pn] ) = CNN_HRFL*\
                           ( inc->r[nn1][pn] );
                        break;

                       default: /* case mn=null: */
/*............................................................................*/
# if DSC_HCRMDE == 0
/* [ in the heat propagation algorithm, that part conflicts with */
/*   skin effect lossy boundaries - hence made inactive, then.] */

                        ( out->i[nn1][pp] ) = ZERO;
                        ( out->i[nn1][pn] ) = ZERO;
# endif
/*............................................................................*/
                        break;
                     };
                  }; /* end if ( mn <= null ) */
               } while(( ++fc ) < FACES );
/*............................................................................*/
# else /* if DSC_ADJCLL == 1 */

               pp = null; do
	       {
		  mn = ( tpt->m[mm][pp] ); /* the adjacent mesh cell index */
		  if ( null < mn )
		  { 
                     hn = ( spt->hh[mn] );

                     if ((( spt->etyp[hn] ) != 't' )
                       &&(( spt->mtyp[hn] ) != 't' ))
                     {
                        nn2 = mn + np;
/*............................................................................*/
/* adjacent port index: */
		        pn = ( tpt->p[mm][pp] );
                        if ( pn < null )
                        {
                           pn = - pn;
                           sgn = - ONE;
                        }
                        else
                           sgn = ONE;

                        pn -= ONE;

                        ( out->i[nn1][pp] ) = sgn * ( inc->r[nn2][pn] );
/*............................................................................*/
# if VICE_VERSA == 1
                        ( out->i[nn2][pn] ) = sgn * ( inc->r[nn1][pp] );
# endif
/*............................................................................*/
                     }; /* end if (( spt->etyp[hn] ) != 't' ) */
                  } 
                  else /* if ( mn <= null ) */
                  { 
                     switch ( mn )
                     {
                       case ELECTRIC_WALL:
                        ( out->i[nn1][pp] ) = CNN_ERFL*\
                           ( inc->r[nn1][pp] );
                        break;
			  
                       case MAGNETIC_WALL:
			( out->i[nn1][pp] ) = CNN_HRFL*\
                           ( inc->r[nn1][pp] );
			break;

                       default: /* case null: */
/*............................................................................*/
# if DSC_HCRMDE == 0
/* [ in the heat propagation algorithm, that part conflicts with */
/*   skin effect lossy boundaries - hence made inactive, then.] */

			( out->i[nn1][pp] ) = ZERO;
# endif
/*............................................................................*/
                        break;
                     };
                  }; /* end if ( mn <= null ) */
               } while(( ++pp ) < PORTS );
/*............................................................................*/
# endif /* DSC_ADJCLL == 1 */
            }; /* end if ( ...typ != 't'rivial ) */
         }; /* end while(( mm++ ) < ( tpt->n )) */
/*............................................................................*/
/* [ time domain - real ] aperiodic boundary conditions: */

	 nn = null;
	 while( nn < ( bpt->n ))
	 {
	    mm = ( bpt->m[nn] );
	    fc = ( bpt->f[nn] );
	    pp = prt1[fc];
	    pn = prt2[fc];

	    nn1 = mm + np;

	    uu = ( bpt->r00[nn] ) * ( inc->r[nn1][pp] ) + \
		 ( bpt->r01[nn] ) * ( inc->r[nn1][pn] );

	    vv = ( bpt->r10[nn] ) * ( inc->r[nn1][pp] ) + \
		 ( bpt->r11[nn] ) * ( inc->r[nn1][pn] );

	    ( out->i[nn1][pp] ) = uu;
	    ( out->i[nn1][pn] ) = vv;

	    nn++;
	 };

	 np = ( tpt->n );
      } while(( ++cnt ) <= prd );
/*............................................................................*/
/* [ time domain ] periodic boundary conditions: */

      nn = null;
      while ( nn < ( bpt->p )) 
      {
	 nn1 = ( bpt->cb[nn] ); /* boundary cell */
	 nn2 = ( bpt->cp[nn] ); /* periodic cell */

	 nn3 = nn1 + ( tpt->n ); /* nn3: boundary field strength label */
	 nn4 = nn2 + ( tpt->n ); /* nn4: periodic field strength label */

	 pn = ( bpt->p0[nn] ) - ONE; /* 1st boundary cell port */
	 pp = ( bpt->pp0[nn] );      /* 1st periodic cell port */

	 if ( pp < null )
	 {
	    pp = - pp;
	    sgn = - ONE;
	 }
	 else
	    sgn = ONE;

	 pp -= ONE;

	 PERDC_TD( );

	 pn = ( bpt->p1[nn] ) - ONE; /* 2nd boundary cell port */
	 pp = ( bpt->pp1[nn] );      /* 2nd periodic cell port */

	 if ( pp < null )
	 {
	    pp = - pp;
	    sgn = - ONE;
	 }
	 else
	    sgn = ONE;

	 pp -= ONE;

	 PERDC_TD( );

	 nn++;
      };

      return out;
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
   } /* end if (( state->dmn ) == 't'ime domain ) */
   else /* if (( state->dmn ) != 't' ): frequency domain */
   {
/*............................................................................*/
# if DSC_DOMAIN == 1 
      printf( "\n Error message from function %s :", __func__ );
      printf( "\n\n Program 'elfe' is compiled in time "
         "domain mode 'DSC_DOMAIN = 1'" );
      printf( "\n - but run in frequency domain mode 2 !!!" );
      printf( "\n\n [ Change time/frequency domain macro DSC_DOMAIN "
         "in program elfe from 1 to" );
      printf( "\n   0: time & frequency domain, or" );
      printf( "\n   2: frequency domain." );
      printf( "\n   Then re-compile and restart program. ]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
      "in_mode_DSC_DOMAIN=1] ", LGS_SIZE );
      return NULL;
/*............................................................................*/
# else /* the remaining cases: DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
      cjg = ONE - 2*CNN_CONJGT;

      mm = null;
      while(( ++mm ) <= ( tpt->n ))
      {
         hm = ( spt->hh[mm] );
         if ((( spt->etyp[hm] ) != 't' )
           &&(( spt->mtyp[hm] ) != 't' ))
         {
/*............................................................................*/
# if DSC_ADJCLL == 0

            fc = null; do
	    {
	       mn = ( tpt->m[mm][fc] ); /* the adjacent mesh cell index */
	       if ( null < mn )
	       {
                  hn = ( spt->hh[mn] );

                  if ((( spt->etyp[hn] ) != 't' )
                    &&(( spt->mtyp[hn] ) != 't' ))
                  {
/*............................................................................*/
/* 1st adjacent port index: */
	             pp = prt1[fc];
	             pn = ( tpt->p[mm][pp] );

	             if ( pn < null )
	             {
		        pn = - pn;
		        sgn = - ONE;
	             }
	             else
		        sgn = ONE;

	             pn -= ONE;

	             ( out->r[mm][pp] ) = CNN_ONE_*sgn*( inc->r[mn][pn] );
	             ( out->i[mm][pp] ) = CNN_ONE_*cjg*sgn*( inc->i[mn][pn] );
/*............................................................................*/
# if VICE_VERSA == 1
	             ( out->r[mn][pn] ) = CNN_ONE_*sgn*( inc->r[mm][pp] );
	             ( out->i[mn][pn] ) = CNN_ONE_*cjg*sgn*( inc->i[mm][pp] );
# endif
/*............................................................................*/
/* 2nd adjacent port index */

		     pp = prt2[fc];
		     pn = ( tpt->p[mm][pp] );

		     if ( pn < null )
		     {
		        pn = - pn;
		        sgn = - ONE;
		     }
		     else
		        sgn = ONE;

		     pn -= ONE; 

		     ( out->r[mm][pp] ) = CNN_ONE_*sgn*( inc->r[mn][pn] );
		     ( out->i[mm][pp] ) = CNN_ONE_*cjg*sgn*( inc->i[mn][pn] );
/*............................................................................*/
# if VICE_VERSA == 1 
		     ( out->r[mn][pn] ) = CNN_ONE_*sgn*( inc->r[mm][pp] );
		     ( out->i[mn][pn] ) = CNN_ONE_*cjg*sgn*( inc->i[mm][pp] );
# endif 
/*............................................................................*/
                  }; /* end if ((( spt->etyp[hn] ) != 't' ) */
	       } 
	       else /* if ( mn <= null ) */
	       { 
		  pp = prt1[fc];
		  pn = prt2[fc];

		  switch ( mn )
		  {
		    case ELECTRIC_WALL:
		     ( out->r[mm][pp] ) = CNN_ERFL*\
			( inc->r[mm][pp] );
		     ( out->i[mm][pp] ) = CNN_ERFL*\
			cjg*( inc->i[mm][pp] );
		     ( out->r[mm][pn] ) = CNN_ERFL*\
			( inc->r[mm][pn] );
		     ( out->i[mm][pn] ) = CNN_ERFL*\
			cjg*( inc->i[mm][pn] );
		     break;
		       
		    case MAGNETIC_WALL:
		     ( out->r[mm][pp] ) = CNN_HRFL*\
                        ( inc->r[mm][pp] );
		     ( out->i[mm][pp] ) = CNN_HRFL*\
                        cjg*( inc->i[mm][pp] );
		     ( out->r[mm][pn] ) = CNN_HRFL*\
                        ( inc->r[mm][pn] );
		     ( out->i[mm][pn] ) = CNN_HRFL*\
                        cjg*( inc->i[mm][pn] );
		     break;

		    default: /* case null: */
/*............................................................................*/
# if DSC_HCRMDE == 0
/* [ in the heat propagation algorithm, that part conflicts with */
/*   skin effect lossy boundaries - hence made inactive, then.] */

		     ( out->r[mm][pp] ) = ZERO;
		     ( out->i[mm][pp] ) = ZERO;
		     ( out->r[mm][pn] ) = ZERO;
		     ( out->i[mm][pn] ) = ZERO;
# endif
/*............................................................................*/
		     break;
		  };
	       }; /* end if ( mn <= null ) */
	    } while(( ++fc ) < FACES );
/*............................................................................*/
# else /* if DSC_ADJCLL == 1 */

	    pp = null; do
	    {
	       mn = ( tpt->m[mm][pp] ); /* the adjacent mesh cell index */
	       if ( null < mn )
	       {  
                  hn = ( spt->hh[mn] );

                  if ((( spt->etyp[hn] ) != 't' )
                    &&(( spt->mtyp[hn] ) != 't' ))
                  {
/*............................................................................*/
/* adjacent port index: */
		     pn = ( tpt->p[mm][pp] );
		     if ( pn < null )
		     {
		        pn = - pn;
		        sgn = - ONE;
		     }
		     else
		        sgn = ONE;

		     pn -= ONE; 

		     ( out->r[mm][pp] ) = CNN_ONE_*sgn*( inc->r[mn][pn] );
		     ( out->i[mm][pp] ) = CNN_ONE_*cjg*sgn*( inc->i[mn][pn] );
/*............................................................................*/
# if VICE_VERSA == 1 
		     ( out->r[mn][pn] ) = CNN_ONE_*sgn*( inc->r[mm][pp] );
		     ( out->i[mn][pn] ) = CNN_ONE_*cjg*sgn*( inc->i[mm][pp] );  
# endif
/*............................................................................*/
                  }; /* end if ((( spt->etyp[hn] ) != 't' ) */
	       } 
	       else /* if ( mn <= null ) */
	       { 
		  switch ( mn )
		  {
		    case ELECTRIC_WALL:
		     ( out->r[mm][pp] ) = CNN_ERFL*\
                        ( inc->r[mm][pp] );
		     ( out->i[mm][pp] ) = CNN_ERFL*\
                        cjg*( inc->i[mm][pp] );
		     break;
		       
		    case MAGNETIC_WALL:
		     ( out->r[mm][pp] ) = CNN_HRFL*\
                        ( inc->r[mm][pp] );
		     ( out->i[mm][pp] ) = CNN_HRFL*\
                        cjg*( inc->i[mm][pp] );
		     break;

		    default: /* case null: */
/*............................................................................*/
# if DSC_HCRMDE == 0
/* [ in the heat propagation algorithm, that part conflicts with */
/*   skin effect lossy boundaries - hence made inactive, then.] */

		     ( out->r[mm][pp] ) = ZERO;
		     ( out->i[mm][pp] ) = ZERO;
# endif
/*............................................................................*/
		     break;
		  };
	       }; /* end if ( mn <= null ) */
	    } while(( ++pp ) < PORTS );
/*............................................................................*/
# endif /* DSC_ADJCLL == 1 */
	 }; /* end if ( ...typ != 't'rivial ) */
      }; /* end while(( mm++ ) < ( tpt->n )) */
/*............................................................................*/
/* [ frequency domain - complex ] aperiodic boundary conditions: */

      nn = null;
      while( nn < ( bpt->n ))
      {
	 mm = ( bpt->m[nn] );
	 fc = ( bpt->f[nn] );
	 pn = prt1[fc];
	 pp = prt2[fc];

	 uu = ( bpt->r00[nn] ) * ( inc->r[mm][pn] ) - \
	      ( bpt->i00[nn] ) * ( inc->i[mm][pn] ) + \
	      ( bpt->r01[nn] ) * ( inc->r[mm][pp] ) - \
	      ( bpt->i01[nn] ) * ( inc->i[mm][pp] );

	 vv = ( bpt->r00[nn] ) * ( inc->i[mm][pn] ) + \
	      ( bpt->i00[nn] ) * ( inc->r[mm][pn] ) + \
	      ( bpt->r01[nn] ) * ( inc->i[mm][pp] ) + \
	      ( bpt->i01[nn] ) * ( inc->r[mm][pp] );

	 ( out->r[mm][pn] ) = uu;
	 ( out->i[mm][pn] ) = cjg*vv;

	 uu = ( bpt->r10[nn] ) * ( inc->r[mm][pn] ) - \
	      ( bpt->i10[nn] ) * ( inc->i[mm][pn] ) + \
	      ( bpt->r11[nn] ) * ( inc->r[mm][pp] ) - \
	      ( bpt->i11[nn] ) * ( inc->i[mm][pp] );

	 vv = ( bpt->r10[nn] ) * ( inc->i[mm][pn] ) + \
	      ( bpt->i10[nn] ) * ( inc->r[mm][pn] ) + \
	      ( bpt->r11[nn] ) * ( inc->i[mm][pp] ) + \
	      ( bpt->i11[nn] ) * ( inc->r[mm][pp] );

	 ( out->r[mm][pp] ) = uu;
	 ( out->i[mm][pp] ) = cjg*vv;

	 nn++;
      };
/*............................................................................*/
/* [ frequency domain - complex ] periodic boundary conditions: */

      nn = null;
      while ( nn < ( bpt->p )) 
      {
	 nn1 = ( bpt->cb[nn] ); /* boundary cell */
	 nn2 = ( bpt->cp[nn] ); /* periodic cell */

	 pn = ( bpt->p0[nn] ) - ONE; /* 1st boundary port */
	 pp = ( bpt->pp0[nn] );      /* 1st periodic port */

	 if ( null < pp )
	 {
	    pp = - pp;
	    sgn = - ONE;
	 }
	 else
	    sgn = ONE;

	 pp -= ONE;

	 PERDC_FD( );

	 pn = ( bpt->p1[nn] ) - ONE; /* 2nd boundary port */
	 pp = ( bpt->pp1[nn] );      /* 2nd periodic port */

	 if ( null < pp )
	 {
	    pp = - pp;
	    sgn = - ONE;
	 }
	 else
	    sgn = ONE;

	 pp -= ONE;

	 PERDC_FD( );

	 nn++;
      }; 

      return out;
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
   }; /* end if (( state->dmn ) != 't' ) */

   return out;
}
/*============================================================================*/
# undef PERDC_TD
# undef PERDC_FD
# undef CNN_UNKNWN
# undef CNN_CONJGT
/*********************** end of file cnctfld-1.0r3.a **************************/
