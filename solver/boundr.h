/* [ file: boundr.h ] */
# define DO_BOUNDR "boundr(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function boundr(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   DSC system boundary input function                                         *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# define BND_DISP 1
/*----------------------------------------------------------------------------*/
/* the following macro should be defined in "../CONFIG.H"                     */
# ifndef DSC_ADJCLL      /* assign neighbouring cell index top.mn[][k] thus:  */
   # define DSC_ADJCLL 0 /* 0: to faces [ k is a face index; 0 <= k < 6 ]     */
# endif                  /* 1: to ports [ k is a port index; 0 <= k < 12 ]    */
/*----------------------------------------------------------------------------*/
static short bndlbl[DSC_JOBS] = {null};
/*============================================================================*/

short boundr( const short jj )
{
   static struct solverstat
      *state = &solver;

/* allusions: */
/*
   extern short bndlbl[];
   extern struct solverstat solver;
   extern struct topology top;
   extern struct boundary bnd;
*/
/* declarations: */

   static FILE
     *bndryfle = NULL;

   static const short
      lbnd =  0,
      cbnd =  8, /* <~ cbnd = clns - 2 */
      clns = 10; /* columns */

   static const char
      slsh  = 47,
      bslsh = 92;

# if DSC_ADJCLL == 0
   static const char /* faces pertinent to ports 0,1,...,11: */
      fce[PORTS] = { 3, 5, 2, 4, 5, 1, 4, 0, 1, 3, 0, 2 };
# endif

# if DSC_HCRMDE == 0
   static const char /* 1st and 2nd ports pertinent to faces 0,...,5: */
      prt1[FACES] = { 7, 5, 11, 9, 3, 1 }, /* 1st port */
      prt2[FACES] = { 10, 8, 2, 0, 6, 4 }; /* 2nd port */
# endif

   static signed char
      pp = null;

   static short
      ii = null,
      kk = null,
      qq = null,
      ind = null,
      llns = null;

   static long
      m = null,
      n = null;

   static double
      phi   = ZERO;

   static char 
      ptr[STS_SIZE] = {null},
     *bndptr = BOUNDARY_FILE, /* cf. main program SOLVER.C */
    **endp = NULL; 

   static const char 
     *scformat = "%80s";

/* prototypes: */

   double sin  ( double x ); 
   double cos  ( double x );
   double fabs ( double x );

   char 
     *lotos ( long lngint, char length );
/*----------------------------------------------------------------------------*/
/* memory allocations: - [void] */
/*----------------------------------------------------------------------------*/
# if DSC_LNGNAMES == 1
   strcpy( solver.bnd, solver.prfx );
   strcat( solver.bnd, bndptr );
# else
   strncpy( solver.bnd, solver.prfx, VSS_SIZE );
   strncat( solver.bnd, bndptr, ( SHS_SIZE - VSS_SIZE - THREE ));
# endif

   strcat( solver.bnd, lotos( bndlbl[jj], null ));

   bndryfle = fopen( solver.bnd, "r" );

   if ( bndryfle == null )
   {
      fprintf( display, "\n\n Error on opening boundary file %s \n ",
                                                               solver.bnd );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":unable_to_open_file_", 
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr );
   strncpy( bnd.name, ptr, ( SHS_SIZE - ONE ));
/*............................................................................*/
/* file_check: */    

   if ( null != strncmp( bnd.name, top.name, THREE ))
   {
      fclose( bndryfle );
      fprintf( display, "\n File error: Inconsistent system identifier"
	                "on file %s !!! ", ( state->bnd ));
      fprintf( display, "\n --- Please verify ! --- \n" );

      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":inconsistent_DSC_model_identifier_",
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( bndryfle, scformat, ptr );
   strncpy( bnd.text, ptr, STS_SIZE - ONE );
   fscanf( bndryfle, scformat, ptr ); /* string "________________..." */
/*............................................................................*/
   fscanf( bndryfle, scformat, ptr ); /* string "electric_operation" */
   fscanf( bndryfle, scformat, ptr ); /* string "TIME_DOMAIN", e.g. */ 

   if( null == strncmp( ptr, "TIME_DOMAIN", FOUR )) 
      pp = ONE;
   else if( null == strncmp( ptr, "FREQUENCY_DOMAIN", FOUR ))
      pp = TWO;
   else
   {
      fprintf( display, "\n\n Error message from function %s:", __func__ );
      fprintf( display, "\n\n Can't read time/frequency domain identifier "
         "%s !!!", ptr );
      fprintf( display, "\n [ Legal is 'TIME_DOMAIN' "
         "or 'FREQUENCY_DOMAIN' ]." );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":illegal_domain_",
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "Maxw.fld_boundaries", e.g. */

   fscanf( bndryfle, scformat, ptr ); /* string "aperiod_bndry_faces:" */
   fscanf( bndryfle, scformat, ptr ); /* string number <bnd.n> */
   bnd.n = strtol( ptr, endp, DEC );

   if ( BNDAP < bnd.n )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many aperiodic boundary faces defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDAP "
         "in SOLVER.CONF or %s.", ( long ) BNDAP, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_aperiodic_bounday_faces_", 
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "periodc_bndry_faces:" */
   fscanf( bndryfle, scformat, ptr ); /* string number <bnd.p> */
   bnd.p = strtol( ptr, endp, DEC );

   if ( BNDPR < bnd.p )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many periodic boundary faces defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDPR "
         "in SOLVER.CONF or %s.", ( long ) BNDPR, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_periodic_bounday_faces_", 
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };
/*............................................................................*/
# if DSC_HCRMDE != 0
/*............................................................................*/
# if DSC_FLDMDE == 0
   fscanf( bndryfle, scformat, ptr ); /* string "thermal_operation" */
# else
   fscanf( bndryfle, scformat, ptr ); /* string "thermal_&_fluid_operation" */
# endif
/*............................................................................*/
   fscanf( bndryfle, scformat, ptr ); /* string "TIME_DOMAIN", e.g. */ 
   fscanf( bndryfle, scformat, ptr ); /* string "thermal_boundaries */
   fscanf( bndryfle, scformat, ptr ); /* string "number_[nodes_or_faces]" */

   fscanf( bndryfle, scformat, ptr ); /* string "skin-eff_heat_sourc:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.nsk = strtol( ptr, endp, DEC );

   if ( BNDSK < bnd.nsk )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many [ skin effect ] heat source faces "
         "defined in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDSK "
         "in SOLVER.CONF or %s.", ( long ) BNDSK, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_heatsource_faces_", 
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "fixed_heat_current_:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.nhc = strtol( ptr, endp, DEC );

   if ( BNDHC < bnd.nhc )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many heatcurrent boundary faces defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDHC "
         "in SOLVER.CONF or %s.", ( long ) BNDHC, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_heatcurr_bounday_faces_", 
      ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "envrnmt-surfce_hcnd:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.nsc = strtol( ptr, endp, DEC );

   if ( BNDSC < bnd.nsc )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many envrmnt-surface heat cond faces "
         "defined in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDSC "
         "in SOLVER.CONF or %s.", ( long ) BNDSC, DO_BOUNDR );
      fprintf( display, "\n - Change macro in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_heatcurr_bounday_faces_", 
      ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "fixed_temperat_fces:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.ntf = strtol( ptr, endp, DEC );

   if ( BNDTF < bnd.ntf )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many boundary temperature faces defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDTF "
         "in SOLVER.CONF or %s.", ( long ) BNDTF, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_boundary_temp_faces_", 
      ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "fixed_temperat_ndes:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.ntn = strtol( ptr, endp, DEC );

   if ( BNDTN < bnd.ntn )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many boundary temperature nodes defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDTN "
         "in SOLVER.CONF or %s.", ( long ) BNDTN, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_boundary_temp_nodes_", 
      ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };
/*............................................................................*/
# if DSC_FLDMDE != 0

   fscanf( bndryfle, scformat, ptr ); /* string "fluid_flow_boundries" */
   fscanf( bndryfle, scformat, ptr ); /* string "number_of_faces" */


   fscanf( bndryfle, scformat, ptr ); /* string "no-slip____________:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.nns = strtol( ptr, endp, DEC );

   if ( BNDNS < bnd.nns )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many no-slip boundary faces defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDNS "
         "in SOLVER.CONF or %s.", ( long ) BNDNS, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_no_slip_fcs_",
      ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "free_slip__________:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.nsl = strtol( ptr, endp, DEC );

   if ( BNDSL < bnd.nsl )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many free slip boundary faces defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDSL "
         "in SOLVER.CONF or %s.", ( long ) BNDSL, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_slip_faces_",
      ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "inflow_____________:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.nif = strtol( ptr, endp, DEC );

   if ( BNDIF < bnd.nif )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many inflow boundary faces defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDIF "
         "in SOLVER.CONF or %s.", ( long ) BNDIF, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_inflw_faces_",
      ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };

   fscanf( bndryfle, scformat, ptr ); /* string "outflow____________:" */
   fscanf( bndryfle, scformat, ptr ); /* long integer string */
   bnd.nof = strtol( ptr, endp, DEC );

   if ( BNDOF < bnd.nof )
   {
      fprintf( display, "\n\n Message from function %s:", __func__ );
      fprintf( display, "\n Too many outflow boundary faces defined "
         "in DSC mesh %s !!!", bnd.name );
      fprintf( display, "\n [ Maximum number is %ld = macro BNDOF "
         "in SOLVER.CONF or %s.", ( long ) BNDOF, DO_BOUNDR );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), DO_BOUNDR, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->bnd ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_outfl_faces_",
      ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   fscanf( bndryfle, scformat, ptr ); /* string, "CELL&FACE_LABELS_...", e.g. */

   if ( pp == ONE )
   {
      for ( ii=null; ii<bnd.n; ii++ )
      {
         fscanf( bndryfle, scformat, ptr );
         bnd.m[ii] = strtol( ptr, endp, DEC ); 
         fscanf( bndryfle, scformat, ptr );
         bnd.f[ii] = strtol( ptr, endp, DEC );

         fscanf( bndryfle, scformat, ptr );
         bnd.r00[ii] = strtod( ptr, endp );
         fscanf( bndryfle, scformat, ptr );
         bnd.r01[ii] = strtod( ptr, endp );
         fscanf( bndryfle, scformat,ptr );
         bnd.r10[ii] = strtod( ptr,endp );
         fscanf( bndryfle, scformat, ptr );
         bnd.r11[ii] = strtod( ptr, endp );

         bnd.i00[ii] = ZERO;
         bnd.i01[ii] = ZERO;
         bnd.i10[ii] = ZERO;
         bnd.i11[ii] = ZERO;
      };
   }
   else /* if pp == TWO: complex domain */
   {
      for ( ii=null; ii<bnd.n ; ii++ )
      {
         fscanf( bndryfle, scformat, ptr );
         bnd.m[ii] = strtol( ptr, endp, DEC );
         fscanf( bndryfle, scformat, ptr );
         bnd.f[ii] = strtol( ptr, endp, DEC );

         fscanf( bndryfle, scformat, ptr ); /* real part of s00=( r00+j*i00 ) */
         bnd.r00[ii] = strtod( ptr, endp );  
         fscanf( bndryfle, scformat, ptr ); /* real part of s01=( r01+j*i01 ) */
         bnd.r01[ii] = strtod( ptr, endp );
         fscanf( bndryfle, scformat, ptr ); /* imaginary part ... */
         bnd.i00[ii] = strtod( ptr, endp ); /* [ and so forth ] */
         fscanf( bndryfle, scformat, ptr ); 
         bnd.i01[ii] = strtod( ptr, endp ); /* reading by lines, the order */
         fscanf( bndryfle, scformat,ptr );  /* follows this matrix storage */
         bnd.r10[ii] = strtod( ptr,endp );  /* scheme: */
         fscanf( bndryfle, scformat, ptr ); /* | ( r00 +    |   ( r01 +   | */ 
         bnd.r11[ii] = strtod( ptr, endp ); /* |   i00*j )  |     i01*j ) | */
         fscanf( bndryfle, scformat,ptr );  /* | -------------------------| */
         bnd.i10[ii] = strtod( ptr,endp );  /* | ( r10 +    |   ( r11 +   | */
         fscanf( bndryfle, scformat, ptr ); /* |   i10*j )  |     i11*j ) | */	
         bnd.i11[ii] = strtod( ptr, endp );
      };
   };

   if ( null < bnd.p )
   { 
      fscanf( bndryfle, scformat, ptr );
      phi = strtod( ptr, endp ); /* phase shift between periodic sections */

      bnd.r = cos(  phi );
      bnd.i = sin( -phi );
   };

   for ( ii=null; ii<bnd.p; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.cb[ii] = strtol( ptr, endp, DEC );  /* periodic boundary cell */
      fscanf( bndryfle, scformat, ptr );   
      bnd.cp[ii] = strtol( ptr, endp, DEC );  /* pertinent periodicity cell */
      fscanf( bndryfle, scformat, ptr );
      bnd.p0[ii] = strtol( ptr, endp, DEC );  /* 1st. periodic boundary port */
      fscanf( bndryfle, scformat, ptr );
      bnd.pp0[ii] = strtol( ptr, endp, DEC ); /* pertinent periodicity port */
      fscanf( bndryfle, scformat, ptr );
      bnd.p1[ii] = strtol( ptr, endp, DEC );  /* 2nd. periodic boundary port */
      fscanf( bndryfle, scformat, ptr );
      bnd.pp1[ii] = strtol( ptr, endp, DEC ); /* pertinent periodicity port */
   };
/*............................................................................*/
# if DSC_HCRMDE != 0
   
   for ( ii=null; ii<bnd.nsk; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.msk[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.fsk[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.rs[ii] = strtod( ptr, endp ); /* surface [skin] resisance */

      bnd.hs[ii] = ZERO; /* initialize smoothed heat sources */

      kk = null; do
      {
         qq = null; do
         {
            fscanf( bndryfle, scformat, ptr );        /* skin effect */
            bnd.sk[ii][kk][qq] = strtod( ptr, endp ); /* penetration depth */
         } while(( ++qq ) < TWO );
      } while(( ++kk ) < TWO );
   };

   for ( ii=null; ii<bnd.nhc; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.mhc[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.fhc[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.hc[ii] = strtod( ptr, endp );

      fscanf( bndryfle, scformat, ptr );
      bnd.rd[ii] = strtod( ptr, endp );
   };

   for ( ii=null; ii<bnd.nsc; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.msc[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.fsc[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.sc[ii] = strtod( ptr, endp ); /* surface heat conductce [W/(K*m^2)] */

      fscanf( bndryfle, scformat, ptr );
      bnd.tr[ii] = strtod( ptr, endp ); /* reference tempertature [C] */
   };

   for ( ii=null; ii<bnd.ntf; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.mtf[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.ftf[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.tf[ii] = strtod( ptr, endp );
   };

   for ( ii=null; ii<bnd.ntn; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.mtn[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.tn[ii] = strtod( ptr, endp );
   };
/*............................................................................*/
# if DSC_FLDMDE != 0
/* no-slip boundary faces: */

   for ( ii=null; ii<bnd.nns; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.mns[ii] = strtol( ptr, endp, DEC ); /* mesh cell index */

      fscanf( bndryfle, scformat, ptr );
      bnd.fns[ii] = strtol( ptr, endp, DEC ); /* face index */
/*............................................................................*/
# if NUSSELT != 0
      fscanf( bndryfle, scformat, ptr );
      bnd.nus[ii] = strtod( ptr, endp ); /* Nusselt number at no-slip face */

      if ( fabs( bnd.nus[ii] ) < 1.e-277 ) /* correct vanishing Nusselt number*/
         bnd.nus[ii] = 1.;                 /* to trivial [ 1.] */ 
# endif
/*............................................................................*/
   };

/* free slip boundary faces: */

   for ( ii=null; ii<bnd.nsl; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.msl[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.fsl[ii] = strtol( ptr, endp, DEC );

/* the outer face normal vector: */
      kk = null; do
      {
         fscanf( bndryfle, scformat, ptr );
         bnd.nf[ii][kk] = strtod( ptr, endp );
      } while(( ++kk ) < THREE );
   };

/* inflow boundary faces: */

   for ( ii=null; ii<bnd.nif; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.mif[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.fif[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.ti[ii] = strtod( ptr, endp ); /* inflow temperature */

      kk = null; do
      {
         fscanf( bndryfle, scformat, ptr );
         bnd.uf[ii][kk] = strtod( ptr, endp ); /* inflow fluid velocity [m/s] */
      } while(( ++kk ) < THREE );
   };

/* outflow boundary faces: */

   for ( ii=null; ii<bnd.nof; ii++ )
   {
      fscanf( bndryfle, scformat, ptr );
      bnd.mof[ii] = strtol( ptr, endp, DEC );

      fscanf( bndryfle, scformat, ptr );
      bnd.fof[ii] = strtol( ptr, endp, DEC );

      kk = null; do
      {
         fscanf( bndryfle, scformat, ptr );
         bnd.vf[ii][kk] = strtod( ptr, endp ); /* outflow fluid velocity [m/s]*/
      } while(( ++kk ) < THREE );
   };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/

   fclose( bndryfle );

   if (( state->rmfle ) != ONE )
   {
      fprintf( display, "\n entered: boundary file %s", ( state->bnd ));
   }
   else if (( state->rmfle ) == ONE )
   {
      for ( kk=jj+ONE; kk<( state->nbrjbs ); kk++ )
      {
         if ( bndlbl[kk] == bndlbl[jj] ) 
         {
            fprintf( display, "\n entered: boundary file %s", ( state->bnd ));

            goto neighb;
         };
      };
      ind = remove( state->bnd );

      if ( ind == null ) 
         fprintf( display, "\n entered and removed: "
            "boundary file %s", ( state->bnd ));
   }; 
/*............................................................................*/
   neighb: /* set neighbouring mesh indices of boundary face [index i] */
           /* and of adjacent face to null */

   for ( ii=null; ii<bnd.n; ii++ )
   {
# if DSC_ADJCLL == 0

# if DSC_HCRMDE == 0
/* [ in the heat propagation algorithm, this part conflicts   */
/*   with [skin effect] lossy boundaries - so made inactive ] */

      m = top.m[bnd.m[ii]][(int)bnd.f[ii]];
      if ( null < m )
      {
         pp = abs( top.p[bnd.m[ii]][(int)prt1[(int)bnd.f[ii]]] );
         if ( null < pp )
            top.m[m][(int)fce[pp-ONE]] = null;

         pp = abs( top.p[bnd.m[ii]][(int)prt2[(int)bnd.f[ii]]] );
         if ( null < pp )
            top.m[m][(int)fce[pp-ONE]] = null;
      };
      top.m[bnd.m[ii]][(int)bnd.f[ii]] = null;
# endif /* DSC_HCRMDE == 0 */
# elif DSC_ADJCLL == 1

# if DSC_HCRMDE == 0
/* [ in the heat propagation algorithm, this part conflicts   */
/*   with [skin effect] lossy boundaries - so made inactive ] */

      m = top.m[bnd.m[ii]][(int)prt1[(int)bnd.f[ii]]];
      if ( null < m )
      {
         pp = abs( top.p[bnd.m[ii]][(int)prt1[(int)bnd.f[ii]]] );
         if ( null < pp )
            top.m[m][pp-ONE] = null;
      };

      m = top.m[bnd.m[ii]][(int)prt2[(int)bnd.f[ii]]];
      if ( null < m )
      {
         pp = abs( top.p[bnd.m[ii]][(int)prt2[(int)bnd.f[ii]]] );
         if ( null < pp )
            top.m[m][pp-ONE] = null;
      };
      top.m[bnd.m[ii]][(int)prt1[(int)bnd.f[ii]]] = null;
      top.m[bnd.m[ii]][(int)prt2[(int)bnd.f[ii]]] = null;
# endif /* DSC_HCRMDE == 0 */
# endif /* DSC_ADJCLL == 1 */
   };

   for ( ii=null; ii<bnd.p; ii++ )
   {
      m=bnd.cb[ii];
      if ( null < m )
      {
         pp = abs( bnd.p0[ii] );
         if ( null < pp )
# if DSC_ADJCLL == 0
            top.m[m][(int)fce[pp-ONE]] = null;
# elif DSC_ADJCLL == 1
            top.m[m][pp-ONE] = null;
# endif
         pp = abs( bnd.p1[ii] );
         if ( null < pp ) 
# if DSC_ADJCLL == 0
            top.m[m][(int)fce[pp-ONE]] = null;
# elif DSC_ADJCLL == 1
            top.m[m][pp-ONE] = null;
# endif
         n = bnd.cp[ii];
         if ( null < n )
         {
            pp = abs( bnd.pp0[ii] );
            if ( null < pp ) 
# if DSC_ADJCLL == 0
               top.m[n][(int)fce[pp-ONE]] = null;
# elif DSC_ADJCLL == 1
               top.m[n][pp-ONE] = null;
# endif
            pp = abs( bnd.pp1[ii] );
            if ( null < pp ) 
# if DSC_ADJCLL == 0
               top.m[n][(int)fce[pp-ONE]] = null;
# elif DSC_ADJCLL == 1
               top.m[n][pp-ONE] = null;
# endif
         };
      };
   };

# if BND_DISP == 1

   fprintf( display, "\n\n text ..................................: " );

   ind = strlen( bnd.text );

   if ( 35 < ind )  
      fprintf( display, "\n %s\n", bnd.text ); 
   else
      fprintf( display, "%s", bnd.text );

   if ( bnd.n > null )
   {  
      fprintf( display, "\n type ..................................:"
                                                             " APERIODIC" );
      fprintf( display, "\n\n aperiodic boundary, %4ld faces: >>---"
                         "---------------------------------------->", bnd.n );
      llns = bnd.n/clns + ONE;

  /* cell1: */

      for ( kk=null; kk<llns; kk++ )
      {
         fprintf( display, "\n -> cell:" ); 
         qq = kk*clns;

         for ( ii=ONE; ii<=clns; ii++ )
         {
            if (( lbnd < kk )&&( cbnd < ii ))
            {
               fprintf( display, "   ***%c", bslsh );
               fprintf( display, "%6ld%c", bnd.m[bnd.n-ONE], bslsh);
               goto faces;
            };    
            if ( qq < bnd.n )   
               fprintf( display, "%6ld%c", bnd.m[qq], bslsh );  
            qq++; 
            if ( bnd.n <= qq  )  
               goto faces; 
         };

        faces:   

         fprintf( display, "\n -> face:" );
         qq = kk*clns;

         for ( ii=ONE; ii<=clns; ii++ ) 
         { 
            if (( lbnd < kk )&&( cbnd < ii))
            {
               fprintf( display, "   ***%c", slsh );
               fprintf( display, "%6d%c", (int)bnd.f[bnd.n-ONE], slsh );
               goto ready1;
            };
            if ( qq < bnd.n )  
	       fprintf( display, "%6d%c", (int)bnd.f[qq], slsh );  
            qq++;
            if ( bnd.n <= qq ) 
	       goto ready1; 
         };
      };            
   };

  ready1:

   if ( null < bnd.p )
   {
      fprintf( display, "\n type ..................................:"
                                                                " PERIODIC" );
      fprintf( display, "\n phi = beta*P ..........................: %.15e",
                                                                        phi );
      fprintf( display, "\n\n periodic boundary, %4ld faces: >>----"
                         "---------------------------------------->", bnd.p );
      llns = ( bnd.p/clns ) + ONE;

  /* cell2: */

      for ( kk=null; kk<llns; kk++ )
      {
         fprintf( display, "\n -> cell:" );
         qq = kk*clns;

         for ( ii=ONE; ii<=clns; ii++ )
         {
            if (( lbnd < kk )&&( cbnd < ii ))
            {
               fprintf( display, "   ***%c", bslsh );
               fprintf( display, "%6ld%c", bnd.cb[bnd.p-ONE], bslsh );
               goto periodic;
            };
            if ( qq < bnd.p )  
	       fprintf( display, "%6ld%c", bnd.cb[qq], bslsh );
            qq++;
            if ( bnd.p <= qq ) 
	       goto periodic;
         };

        periodic:

         fprintf( display, "\n -> prdc:" );
         qq = kk*clns;

         for ( ii=ONE; ii<=clns; ii++ )
         {
            if (( lbnd < kk )&&( cbnd < ii))
            {
               fprintf( display, "   ***%c", slsh );
               fprintf( display, "%6ld%c", bnd.cp[bnd.p-ONE], slsh );
               goto ready2;
            };
            if ( qq < bnd.p )  
	       fprintf( display, "%6ld%c", bnd.cp[qq], slsh );
            qq++;
            if ( bnd.p <= qq ) 
	       goto ready2;
         };
      };
   };

  ready2:

   fprintf( display, "\n --------------------------------------"
                      "----------------------------------------");
# endif

   return ONE; /* normal return */
}
/*============================================================================*/
# undef BND_DISP
/***************** end of boundary input function boundr(*) *******************/

