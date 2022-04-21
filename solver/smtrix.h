/* [ file: smtrix.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function smtrix(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   This subroutine enters nodal S-parameters of a DSC system,                 *
*   named top.name<n>, from file DSC_PRFX.S_MATRIX_FILE<n>.                    *
*   The topology of the referenced system is stored in file                    *
*   DNS_PRFX.TOP_FILE<n> and transferred to the main program                   *
*   by function systpl(*).                                                     *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define SMX_DISP 1
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
# include "../math/expmtx.h"
# include "../math/gssjtp.h"
# include "../math/jacbtp.h"
/*----------------------------------------------------------------------------*/
# define SMX_GYMGCRR 0       /* SMX_GYMGCRR 1: evaluate average g-current     */
                             /* density [ integral mean value ] over time     */
                             /* interval I[k] = [ k*dt , (k+1)*dt ); k=0,1,...*/
                             /* SMX_GYMGCRR 0: evaluate left limit of g-curr. */
                             /* at times t[k] = dt*(2*k+1)/2 minus.           */
                             /* [ For entire consistency, this macro should   */
                             /*   coincede with same macro in subroutine      */
                             /*  'smatrx(*)' of program 'elsy.c'. ]           */

# define SMX_GYELCRR 0       /* SMX_GYELCRR 0: evaluate g-current             */
                             /* at times t[k] = dt*(2*k+1)/2.                 */
                             /* [ For entire consistency, this macro should   */
                             /*   coincede with same macro in subroutine      */
                             /*  'smatrx(*)' of program 'elsy.c'. ]           */

/* computational constants [ to be changed with care - only by experts ]:.....*/
/*----------------------------------------------------------------------------*/
# ifndef PRECISION
   # define PRECISION ( double )( 1.000e-14 )
# endif
/*----------------------------------------------------------------------------*/
# define SMX_ONE_PLUS  ( 1.+ 33.*PRECISION )
# define SMX_ONE_MINUS ( 1.- 33.*PRECISION )
# define SMX_SCMPR   ( 1.0000000000000 ) /* S-parameter regularization factor */
# define SMX_GCMPR   ( 1.0000000000000 ) /* gyromagn. S-parameter regl.factor */
# define SMX_GYROBND ( 1.0e-277 )
# define SMX_SNGLBND ( 1.0e-15 )         /* singularity bound for matrix U    */
# define SMX_TRIVIAL ( 1.0e-277 )        /* triviality bound                  */
/*----------------------------------------------------------------------------*/
/* natural constants - to be changed only with the help of God:               */
/* [ all units are international units (mks), if not specified otherwise.]    */

# define EPS_VAC      ( 8.8541878170e-12 ) /* vac. permittivity [A*sec/(V*m)] */
# define MY_VAC_      ( 1.2566370614e-06 ) /* "    permeability [V*sec/(A*m)] */
# define ELECTR_CHRGE (-1.6021773349e-19 ) /* electron charge [Coulomb=A*sec] */
# define ELECTR_MASS_ ( 9.1093897540e-31 ) /* electron mass [kg]              */
# define BOHRs_MAGNTN ( 1.1654071500e-29 ) /* Bohr's magneton [Volt*sec*m]    */
# define PLANCKs_CNST ( 6.6260755400e-34 ) /* Planck's constant [Joule*sec]   */
/*----------------------------------------------------------------------------*/
/* macros: */

# define CLSMX(NN) \
{ \
   fclose( smatrix ); \
\
   ii = (NN); \
   if ( ii != ONE ) \
   { \
      fprintf( display, "\n entered: S-parameters %s", ( state->smx )); \
   } \
   else if ( ii == ONE ) \
   { \
      for ( kk=(jj+ONE); kk<state->nbrjbs; kk++ ) \
      { \
         if ( smxlbl[kk] == smxlbl[jj] ) \
         { \
            fprintf( display, "\n entered: S-parameters %s", ( state->smx )); \
            goto disp; \
         }; \
      }; \
      ind = remove ( state->smx ); \
      if ( ind == null ) \
         fprintf( display, "\n entered and removed:" \
            " S-parameters %s", ( state->smx )); \
   }; \
}
/*----------------------------------------------------------------------------*/
typedef struct
{
   signed char
      rtn;

   double 
    no,tr,ec,ld,tg,mc,
    mi[DIMNS],
    ms[DIMNS],
    hg[DIMNS],
    ue[DIMNS][DIMNS],
    um[DIMNS][DIMNS],
    ve[DIMNS][DIMNS],
    vm[DIMNS][DIMNS],
    qe[DIMNS][DIMNS],
    qm[DIMNS][DIMNS],
   eue[DIMNS][DIMNS],
   eum[DIMNS][DIMNS];
}  GYRMTX; 
static GYRMTX gyr = {null};

typedef struct
{
   signed char
      rtn;

   double aa[DIMNS][DIMNS],
          bb[DIMNS][DIMNS];
}  FORMTX;
static FORMTX frm = {null};
/*----------------------------------------------------------------------------*/
static short smxlbl[DSC_JOBS] = {null};
static GAUSS_JRD gss = {null};
static JACOBI_EV jac = {null};
/*============================================================================*/

short smtrix( const short jj )
{
/*............................................................................*/
/* allusions: */
/*
   extern struct solverstat solver;
   extern struct topology top;
   extern struct tlmsmx smx;
   extern JACOBI_EV jac;
   extern FORMTX frm;
   extern GYRMTX gyr;

# if DSC_DOMAIN == 0
   extern DSC_DEFLCT dfl;
# elif DSC_DOMAIN == 1
   extern DSC_DEFLCT dfl;
# endif

   extern short smxlbl[];
*/
/*............................................................................*/
/* declarations: */

   static struct solverstat
     *state = &solver;

   static FORMTX
     *fmp = &frm;
/*
   static GYRMTX
     *gyp = &gyr;
*/
   static JACOBI_EV
     *jev = &jac;

   static FILE
     *smatrix = NULL;

   static signed char
      qq = null,
      pp = null;

   static short
      kk = null,
      ll = null,
     ind = null;

   static long
      hh = null,
      ii = null,
      me = null,
      mm = null,
      nn = null,
      icll = null,
      fcll = null;

   static double
      hnorm = ZERO,
      hnrmi = ZERO;

   static char 
      cmp[STS_SIZE] = {null},
      ptr[STS_SIZE] = {null},
    **endp = NULL;
/*............................................................................*/
# if DSC_HCRMDE != 0
   static char 
      ttyp[STS_SIZE] = {null};

   static struct hcrsmx
     *hsp = &hcs;
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* the following character strings MUST coincede with equally named strings */
/* in function elsydrv(*) [ cf. directory elsy ] */

   static const char
     *smtrix = "-smx",
     *fields = "Maxwell_field",
     *smxpar= ">>-S-PARAMETERS->>", 
     *smxidx= ">>-S-PARAMETER_LABELS->>",
     *smxptr = S_MATRIX_FILE;
/*............................................................................*/
# if DSC_HCRMDE != 0
   static const char
     *thermal = "Thermal";
# endif
/*............................................................................*/
   static const char
     *scformat = "%80s";

/* prototypes: */

   double sqrt( double x );

   JACOBI_EV 
     *jacobi( JACOBI_EV *jcp );

   int expmtx( short rank );
   short exptue( );
   short exptum( );
   short vebinv( long h ); 
   short vmbinv( long h );

   char 
      *lotos( long h, char qq );
/*----------------------------------------------------------------------------*/
/* allocate memory: - */
/*............................................................................*/
   gyr.ec = ELECTR_CHRGE/ELECTR_MASS_;                       /* [ Coulomb/kg] */
   gyr.mc = -2.*PI*BOHRs_MAGNTN/PLANCKs_CNST;                /* [ m/Coulomb ] */
   smx.adm = sqrt( EPS_VAC/MY_VAC_ );
/*............................................................................*/
# if DSC_LNGNAMES == 1
   strcpy( state->smx, state->prfx );
   strcat( state->smx, smxptr );
# else
   strncpy( state->smx, state->prfx, VSS_SIZE );
   strncat( state->smx, smxptr, ( SHS_SIZE - VSS_SIZE - THREE ));
# endif
/*............................................................................*/
# if DSC_HCRMDE != 0
/*............................................................................*/
# if DSC_FLDMDE != 0
/* volumes and dynamic mean temperatures of  */ 
/* fluid connected components reset to ZERO: */ 

   ii = null; do
   {
      ( hsp->cnn[ii] ) = null;
      ( hsp->volume[ii] ) = ZERO;
      ( hsp->tmean[ii] ) = ZERO;
      ( hsp->crsdt[ii] ) = ZERO;
/*...........................................................................*/
# if CMPRSSBL != 0    
      ( hsp->mass[ii] ) = ZERO;
# endif /* CMPRSSBL != 0 */   
/*...........................................................................*/
   } while (( ii++ ) < NFCNN );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* open s-matrix file: */

   strcat( state->smx, lotos( smxlbl[jj], null ));

   smatrix = fopen( state->smx ,"r");

   if ( smatrix == null )
   {
      fprintf( stderr, "\n\n Error on opening S-parameters file %s\n ",
         state->smx );
      strncpy( state->fcterr, __func__, SHS_SIZE );
      strncpy( state->errmsg, state->smx, SHS_SIZE );
      strncat( state->errmsg, ":unable_to_open_file_", 
         ( LGS_SIZE - SHS_SIZE )); 

      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( smatrix, scformat, ptr ); /* char string: DSC model identifier */
   strncpy( smx.name, ptr, SHS_SIZE - ONE );
/*............................................................................*/
/* file_check: */    

   if( null != strncmp( smx.name, top.name, THREE ))
   {
      CLSMX(null);
      fprintf( stderr, "\n File error: Inconsistent system identifier"
	 "in file %s !!!", state->smx );
      fprintf( stderr, "\n --- Please verify !!! --- \n" );
      strncpy( state->fcterr, __func__, SHS_SIZE );
      strncpy( state->errmsg, state->smx, SHS_SIZE );
      strncat( state->errmsg, ":inconsistent_DSC_model_identifier_",
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( smatrix, scformat, ptr ); /* connected char string [ comment ] */
   strncpy( smx.text, ptr, STS_SIZE - ONE );
   fscanf( smatrix, scformat, ptr ); /* string "_____________________..." */
/*............................................................................*/
   fscanf( smatrix, scformat, ptr ); /* string "first_mesh_cell_index_:" */
   fscanf( smatrix, scformat, ptr ); /* icll = initial mesh index */
   icll = strtol( ptr, endp, DEC );   /* [ should be equal to ONE ] */
                                              
   if ( icll != ONE )
   {
      CLSMX(null);
      fprintf( stderr, "\n Truncated S_matrix file %s: ", state->smx );
      fprintf( stderr, "\n First mesh cell index %ld != ONE", icll );
      strncpy( state->fcterr, __func__, SHS_SIZE );
      strncpy( state->errmsg, state->smx, SHS_SIZE );
      strncat( state->errmsg, ":truncated_cell_labels_", 
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( smatrix, scformat, ptr ); /* string "last_mesh_cell_index_:" */
   fscanf( smatrix, scformat, ptr );
   fcll = strtol( ptr, endp, DEC );

   if ( fcll != top.n )
   {
      fclose( smatrix );

      fprintf( stderr, "\n Error on reading s-matrix file %s:", state->smx );
      fprintf( stderr, "\n Inconsistent number of cells !!!\n" );
      strncpy( state->fcterr, __func__, SHS_SIZE );
      strncpy( state->errmsg, state->smx, SHS_SIZE );
      strncat( state->errmsg, ":inconsistent_cell_number_",
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( smatrix, scformat, ptr ); /* frequency/time domain identifier  */
                                      /* [ only Maxwell field algorithm ] */
   if (( *ptr == 'F' )
     ||( *ptr == 'f' ))
      state->dmn = 'f'; /* [ 'f'requency domain ] */
   else if (( *ptr == 'T' )
          ||( *ptr == 't' ))
      state->dmn = 't'; /* [ 't'ime domain ] */
   else
   {
      fclose( smatrix );

      fprintf( stderr, "\n Error in S_matrix file %s:", state->smx );
      fprintf( stderr, "\n Illeagal domain identifier !!!\n" );
      strncpy( state->fcterr, __func__, SHS_SIZE );
      strncpy( state->errmsg, state->smx, SHS_SIZE );
      strncat( state->errmsg, ":illeagal_domain_",
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
/*............................................................................*/
/* parameters, Maxwell field algorithm: */

   if( state->dmn == 'f' )
   {
      fscanf( smatrix, scformat, ptr ); /* string "frequency_[Hz]_____:" */
      fscanf( smatrix, scformat, ptr );
      state->fr = strtod( ptr, endp );

      fscanf( smatrix, scformat, ptr ); /* string "phase_shift_[rad]__:" */
      fscanf( smatrix, scformat, ptr );
      state->ph = strtod( ptr, endp );

      smx.shr = cos( state->ph );
      smx.shi = sin( state->ph );
      state->dt = ( double) ONE;
   }
   else if ( state->dmn == 't' )
   {
      fscanf( smatrix, scformat, ptr ); /* string "timestep_[s]_______:" */
      fscanf( smatrix, scformat, ptr ); 
      state->dt = strtod( ptr, endp );

      fscanf( smatrix, scformat, ptr ); /* string "mean_frequency_[Hz]:" */
      fscanf( smatrix, scformat, ptr );
      state->fr = strtod( ptr, endp );

      state->ph = ZERO;

      hnorm = ZERO;
   };
/*............................................................................*/
# if DSC_HCRMDE != 0
   fscanf( smatrix, scformat, ptr ); /* string "HCR_timestep_[s]:", e.g. */
   fscanf( smatrix, scformat, ptr ); /* heat current time step [ seconds ] */
   state->hcdt = strtod( ptr, endp );
   fscanf( smatrix, scformat, ptr ); /* string "thermal_start[s]:" */
   fscanf( smatrix, scformat, ptr ); /* start thermal routines with that  */
   state->hcrstart = strtod( ptr, endp ); /* delay [ seconds ] */
/*............................................................................*/
# if DSC_FLDMDE != 0
   fscanf( smatrix, scformat, ptr ); /* string "fluid_start__[s]:" */
   fscanf( smatrix, scformat, ptr ); /* start fluid routines with that  */
   state->fldstart = strtod( ptr, endp ); /* delay [ seconds ] */
   fscanf( smatrix, scformat, ptr ); /* string "connected_fluid_cpts" */
   fscanf( smatrix, scformat, ptr ); /* number of fluid connected compts */
   ( hsp->cnn[null] ) = strtol( ptr, endp, DEC );

   ( hsp->crsdt[null] ) = 4.*( state->hcdt ); /* [ default ] */

   if ( null < ( hsp->cnn[null] ))
   {
      fscanf( smatrix, scformat, ptr ); /* string "COARSENING_PERIOD[S]" */
      
      kk = ONE; do
      {
         fscanf( smatrix, scformat, ptr ); /* string "component..." */
         fscanf( smatrix, scformat, ptr ); /* coarsening time period [s] */
         ( hsp->crsdt[kk] ) = strtod( ptr, endp );

	 if (( hsp->crsdt[kk] ) < SMX_TRIVIAL )
            ( hsp->crsdt[kk] ) = ( hsp->crsdt[( kk - ONE )] ); /* [ default ] */
      } while(( kk++ ) < ( hsp->cnn[null] ));
   }; /* end if ... */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   smx.ge = null;
   smx.gm = null; 

   me = null;
   mm = null;
/*............................................................................*/
/* enter E/H field s-parameters: */

   fscanf( smatrix, scformat, ptr ); /* string ">>-S-PARAMETERS->>", e.g. */

   if ( null != strncmp( ptr, smxpar, SIX ))
   {
      fprintf( stderr, "\n\n Error message from function %s :", __func__ );
      fprintf( stderr, "\n Corrupted s-matrix file %s,", state->smx );
      fprintf( stderr, "\n different FORMER.C and SOLVER.C releases (?)," );
      fprintf( stderr, "\n or incompatible compiler options (?)." );
      strncpy( state->fcterr, __func__, SHS_SIZE ); 
      strncpy( state->errmsg, state->smx, SHS_SIZE );
      strncat( state->errmsg, ":corrupted_file",
         ( LGS_SIZE - SHS_SIZE - TEN ));

      return null; /* abnormal return */
   };
/*...........................................................................*/
   fscanf( smatrix, scformat, ptr );

   strcpy( cmp, fields );
   strcat( cmp, smtrix );

   hh = null;
   while ( null == strncmp( cmp, ptr, SIX )) /* cmp = "s-matrix", e.g.*/
   {                             
      hh++;
      if ( NRSMX < hh )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many s-matrices defined "
            "in DSC system %s !!!", smx.name );
         fprintf( stderr, "\n [ Maximum number is %ld = macro NRSMX "
            "fixed in file '%s'.", ( long ) NRSMX, "SOLVER.CONF" );
         fprintf( stderr, "\n - Change macro in compliance with "
            "memory resources.]\n " );
         strncpy( state->fcterr, __func__, SHS_SIZE );
         strncpy( state->errmsg, state->smx, SHS_SIZE );
         strncat( state->errmsg, ":too_many_s-matrices_defined_",
            ( LGS_SIZE - SHS_SIZE ));

         return null; /* abnormal return */
      };
/*............................................................................*/
      fscanf( smatrix, scformat, ptr ); /* type of node or s-matrix */
      smx.etyp[hh] = ptr[null];

      if ( smx.etyp[hh] != 't' )         /* non_trivial e_node */
      {
         if ( smx.etyp[hh] == 'g' )      /* gyroelectric node */
         {
            if ( GENDS <= me )
            {
               fprintf( stderr, "\n\n Error message from function %s :",
                   __func__);
               fprintf( stderr, "\n Too many gyroelectric nodes defined "
                  "in DSC system %s !!!", smx.name );
               fprintf( stderr, "\n [ Maximum number is %ld = macro GENDS "
                  "fixed in file '%s'.", ( long ) GENDS, "SOLVER.CONF" );
               fprintf( stderr, "\n - Change macro in compliance with "
                  "memory resources.]\n " );
               strncpy( state->fcterr, __func__, SHS_SIZE );
               strncpy( state->errmsg, state->smx, SHS_SIZE );
               strncat( state->errmsg, ":too_many_gyroelectric_nodes_",
                  ( LGS_SIZE - SHS_SIZE ));

               return null; /* abnormal return */
            };
/*............................................................................*/
# if DSC_DOMAIN == 0
            ll = null; do
            {
               dfl.ge[me][ll] = ZERO; /* initialize g-current vct. = ZERO */
               dfl.de[me][ll] = ZERO; /* initialize deflection vector  "" */
            } while(( ++ll ) < DIMNS );
# elif DSC_DOMAIN == 1
            ll = null; do
            {
               dfl.ge[me][ll] = ZERO; /* initialize g-current vct. = ZERO */
               dfl.de[me][ll] = ZERO; /* initialize deflection vector  "" */
            } while(( ++ll ) < DIMNS );
# endif
            me++; 

            if ( GESMX <= smx.ge )
            {
               fprintf( stderr, "\n\n Error message from function %s :",
                  __func__);
               fprintf( stderr, "\n Too many gyroelectric s-matrices defined "
                  "in DSC system %s !!!", smx.name );
               fprintf( stderr, "\n [ Maximum number is %ld = macro GESMX "
                  "fixed in file '%s'.", ( long ) GESMX, "SOLVER.CONF" );
               fprintf( stderr, "\n - Change macro in compliance with "
                  "memory resources.]\n " );
               strncpy( state->fcterr, __func__, SHS_SIZE );
               strncpy( state->errmsg, state->smx, SHS_SIZE ); 
               strncat( state->errmsg, ":too_many_gyroelectric_s-matrices",
                  ( LGS_SIZE - SHS_SIZE ));

               return null; /* abnormal return */
            };
/*............................................................................*/
            if ( smx.ge == null )
            {
               ll = null; do                        /* internal magnetic flux */
               {                                    /* density vector [Tesla] */
                  fscanf( smatrix, scformat, ptr );
                  gyr.mi[ll] = strtod( ptr, endp );
               } while(( ++ll ) < DIMNS );

               fscanf( smatrix, scformat, ptr );/* plasma electron density*/
               gyr.no = strtod( ptr, endp );    /* [ 1 / m^3 ]            */
               fscanf( smatrix, scformat, ptr );/* relaxation time [sec]  */
               gyr.tr = strtod( ptr, endp );
            };
/*............................................................................*/
/* enter area matrix A: */

            ll = null; do
            {
               kk = null; do
               {
                  fscanf( smatrix, scformat, ptr );
                  ( fmp->aa[ll][kk] ) = strtod( ptr, endp );
               } while(( ++kk ) < DIMNS );
            } while(( ++ll ) < DIMNS );
           
            fscanf( smatrix, scformat, ptr ); /* e_mode S-matrix type     */

         }; /* end if smx.etyp == 'g'yroelectric */
/*............................................................................*/
/* enter e-mode S-parameters: */             
                                       /* pp=0: asymmetric                  */
         pp  = ( *ptr == 'd' );        /* pp=1: diagonal                    */
         pp += ( TWO*( *ptr == 's' )); /* pp=2: symmetric non-diagonal      */
         ll = null; do
         {
            if ( state->dmn == 't' )
            {
               switch( pp )
               {
                 case 1:
                  fscanf( smatrix, scformat, ptr ); /* diagonal block KE */
                  smx.se[hh][ll][ll] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  fscanf( smatrix, scformat, ptr ); /* diagonal block ME */
                  smx.se[hh][ll+DIMNS][ll] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  fscanf( smatrix, scformat, ptr ); /* diagonal block NE */
                  smx.se[hh][ll+DIMNS][ll+DIMNS] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  smx.se[hh][ll][ll+DIMNS] =\
                     SMX_SCMPR;                      /* diagonal block LE */

                  kk = null;                   /* all off-diagonal */
                  while( kk < ll )             /* elements set to ZERO */
                  {                           
                     smx.se[hh][ll][kk] = ZERO;
                     smx.se[hh][kk][ll] = ZERO;
                     smx.se[hh][ll+DIMNS][kk] = ZERO;
                     smx.se[hh][kk+DIMNS][ll] = ZERO;
                     smx.se[hh][ll][kk+DIMNS] = ZERO;
                     smx.se[hh][kk][ll+DIMNS] = ZERO;
                     smx.se[hh][ll+DIMNS][kk+DIMNS] = ZERO;
                     smx.se[hh][kk+DIMNS][ll+DIMNS] = ZERO;
                     kk++ ;
                  };
                  break;
                
                 case 2:
                  kk = null; do 
                  {
                     fscanf( smatrix, scformat, ptr );  /* block KE */
                     smx.se[hh][ll][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr );  /* block ME */
                     smx.se[hh][ll+DIMNS][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr );  /* block NE */
                     smx.se[hh][ll+DIMNS][kk+DIMNS] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     smx.se[hh][ll][kk+DIMNS] =\
                        SMX_SCMPR*( double )( ll == kk ); /* block LE */

                     if ( ll != kk )
                     {
                        smx.se[hh][kk][ll] = smx.se[hh][ll][kk];
                        smx.se[hh][kk+DIMNS][ll] = smx.se[hh][ll+DIMNS][kk];
                        smx.se[hh][kk][ll+DIMNS] = smx.se[hh][ll][kk+DIMNS];
                        smx.se[hh][kk+DIMNS][ll+DIMNS] =\
                           smx.se[hh][ll+DIMNS][kk+DIMNS];
                     };
                  } while (( ++kk ) <= ll );
                  break;

                 default:
                  kk = null; do
                  {
                     fscanf( smatrix, scformat, ptr );  /* block KE */
                     smx.se[hh][ll][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr );  /* block ME */
                     smx.se[hh][ll+DIMNS][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr );  /* block NE */
                     smx.se[hh][ll+DIMNS][kk+DIMNS] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     smx.se[hh][ll][kk+DIMNS] =\
                        SMX_SCMPR*( double )( ll == kk ); /* block LE */
                  } while (( ++kk ) < DIMNS );
                  break;
               }; /* end: switch( pp ) */
            } /* end if state->dmn == 't' */
            else if ( state->dmn == 'f' )
            {
               switch( pp )
               {
                 case 1:
                  fscanf( smatrix, scformat, ptr ); /* real part */
                  smx.se[hh][ll][ll] =\
                     SMX_SCMPR*strtod( ptr, endp );
                                      
                  fscanf( smatrix, scformat, ptr ); /* imaginary part */
                  smx.se[hh][ll][ll+DIMNS] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  kk = null;                   /* off-diagonal elements */
                  while( kk < ll )             /* set to ZERO */
                  {
                     smx.se[hh][ll][kk] = ZERO;
                     smx.se[hh][kk][ll] = ZERO;
                     smx.se[hh][ll][kk+DIMNS] = ZERO;
                     smx.se[hh][kk][ll+DIMNS] = ZERO;
                     kk++ ;
                  };
                  break;

                 case 2:
                  kk = null; do
                  {
                     fscanf( smatrix, scformat, ptr ); /* real part */
                     smx.se[hh][ll][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr ); /* imaginary part */
                     smx.se[hh][ll][kk+DIMNS] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     if( kk != ll )
                     {
                        smx.se[hh][kk][ll] = smx.se[hh][ll][kk];
                        smx.se[hh][kk][ll+DIMNS] = smx.se[hh][ll][kk+DIMNS];
                     };
                  } while (( ++kk ) <= ll );
                  break;

                 default:
                  kk = null; do
                  {
                     fscanf( smatrix, scformat, ptr ); /* real part */
                     smx.se[hh][ll][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr ); /* imaginary part */
                     smx.se[hh][ll][kk+DIMNS] =\
                        SMX_SCMPR*strtod( ptr, endp );
                  } while (( ++kk ) < DIMNS );
                  break;
               };/* end switch ... */
            }; /* end if state->dmn == 'f' */
         } while (( ++ll ) < DIMNS );  
/*............................................................................*/
         if ( smx.etyp[hh] == 'g' )
         {
/* processing gyroelectric scattering parameters: */

            smx.me[hh] = smx.ge;

            if ( smx.ge == null )
            {
/*............................................................................*/
/* The following function exptue(*) conditionally returnes matrices exu.er[][], 
   ... in a way dependent on g-current integration macro SMX_GYELCRR :            
   If one [ cell independent ] set of gyroelectric parameters, gyr.mi, gyr.no,
   gyr.tr, .... , is given only one initial call of function 'exptue' is      
   sufficient.                                                                
   Otherwise, for mesh dependent parameters, this must be modified by calling 
   'exptue' each time new mesh parameters are entered.                        */
/*............................................................................*/
/* returns matrices gyr.ue[][], gyr.ve[][], gyr.eue[][], gyr.qe[][]           */
                                               /*.............................*/
               ind = exptue( );               /*                              */
/*..........................................*/                           
            };
/*............................................................................*/
/* returns matrices smx.eue[][], smx.de[][] and smx.sge[][]                   */
                                           /*.................................*/
            ind = vebinv( hh );           /*                                  */
/*......................................*/
            smx.ge++;
         };/* end if ...etyp == 'g'yroelectric */
      }; /* end if ( ...etyp != 't'rivial ) */
/*............................................................................*/
      fscanf( smatrix, scformat, ptr );
      smx.mtyp[hh] = ptr[null];

      if ( smx.mtyp[hh] != 't' )
      {
         if ( smx.mtyp[hh] == 'g' )       /* gyromagnetic node !              */
         {
            if ( GMNDS <= mm )
            {
               fprintf( stderr, "\n\n Error message from function %s :",
                  __func__);
               fprintf( stderr, "\n Too many gyromagnetic nodes defined "
                  "in DSC system %s !!!", smx.name );
               fprintf( stderr, "\n [ Maximum number is %ld = macro GMNDS "
                  "fixed in file '%s'.", ( long ) GMNDS, "SOLVER.CONF" );
               fprintf( stderr, "\n - Change macro in compliance with "
                  "memory resources.]\n " );
               strncpy( state->fcterr, __func__, SHS_SIZE );
               strncpy( state->errmsg, state->smx, SHS_SIZE );
               strncat( state->errmsg, ":too_many_gyromagnetic_nodes_",
                  ( LGS_SIZE - SHS_SIZE ));

               return null; /* abnormal return */
            };
/*............................................................................*/
# if DSC_DOMAIN == 0
            ll = null; do
            {
               dfl.gm[mm][ll] = ZERO; /* initialize g-current vct. = ZERO */
               dfl.dm[mm][ll] = ZERO; /* initialize deflection vector  "" */
            } while(( ++ll ) < DIMNS );
# elif DSC_DOMAIN == 1
            ll = null; do
            {
               dfl.gm[mm][ll] = ZERO; /* initialize g-current vct. = ZERO */
               dfl.dm[mm][ll] = ZERO; /* initialize deflection vector  "" */
            } while(( ++ll ) < DIMNS );
# endif
            mm++; 

            if ( GMSMX <= smx.gm )
            {
               fprintf( stderr, "\n\n Error message from function %s :",
                  __func__);
               fprintf( stderr, "\n Too many gyromagnetic s-matrices defined "
                  "in DSC system %s !!!", smx.name );
               fprintf( stderr, "\n [ Maximum number is %ld = macro GMSMX "
                  "fixed in file '%s'.", ( long ) GMSMX, "SOLVER.CONF" );
               fprintf( stderr, "\n - Change macro in compliance with "
                  "memory resources.]\n " );
               strncpy( state->fcterr, __func__, SHS_SIZE );
               strncpy( state->errmsg, state->smx, SHS_SIZE );
               strncat( state->errmsg, ":too_many_gyromagnetic_s-matices",
                  ( LGS_SIZE - SHS_SIZE ));

               return null; /* abnormal return */
            };
/*............................................................................*/
            if ( smx.gm == null )
            {
               ll = null; do
               {                             
/* saturation magnetization [Tesla]: */
                  fscanf( smatrix, scformat, ptr );
                  gyr.ms[ll] = strtod( ptr, endp ); 
/* internal magnetic field [A/m]: */
                  fscanf( smatrix, scformat, ptr );
                  gyr.hg[ll] = strtod( ptr, endp ); 
               } while(( ++ll ) < DIMNS );

               fscanf( smatrix, scformat, ptr );
               gyr.ld  = strtod( ptr, endp );
               fscanf( smatrix, scformat, ptr );
               gyr.tg  = strtod( ptr, endp );
            };
/*............................................................................*/
/* enter area matrix A: */

            ll = null; do 
            {
               kk = null; do
               {
                  fscanf( smatrix, scformat, ptr );
                  fmp->aa[ll][kk] = strtod( ptr, endp );
               } while(( ++kk ) < DIMNS );
            } while(( ++ll ) < DIMNS );
           
            fscanf( smatrix, scformat, ptr ); /* h_mode S-matrix type  */

         }; /* end if smx.mtyp == 'g'yromagnetic */
/*............................................................................*/
/* enter h-mode S-parameters: */
                                       /* pp=0: asymmetric */
         pp = ( *ptr == 'd' );         /* pp=1: diagonal */
         pp += ( TWO*( *ptr == 's' )); /* pp=2: symmetric [ non-diagonal ] */
                                           
         ll = null; do 
         {
            if ( state->dmn == 't' ) 
            {
               switch( pp )
               {
                 case 1:
                  fscanf( smatrix, scformat, ptr ); /* diagonal, block KH */
                  smx.sh[hh][ll][ll] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  fscanf( smatrix, scformat, ptr ); /* diagonal, block MH */ 
                  smx.sh[hh][ll+DIMNS][ll] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  fscanf( smatrix, scformat, ptr ); /* diagonal, block NH */
                  smx.sh[hh][ll+DIMNS][ll+DIMNS] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  smx.sh[hh][ll][ll+DIMNS] =\
                     SMX_SCMPR;                      /* unit matrix block LH */

                  kk = null;                   /* off-diagonal elements */
                  while( kk < ll )             /* set to ZERO */
                  {                           
                     smx.sh[hh][ll][kk] = ZERO;
                     smx.sh[hh][kk][ll] = ZERO; 
                     smx.sh[hh][ll+DIMNS][kk] = ZERO;
                     smx.sh[hh][kk+DIMNS][ll] = ZERO;
                     smx.sh[hh][ll][kk+DIMNS] = ZERO;
                     smx.sh[hh][kk][ll+DIMNS] = ZERO;
                     smx.sh[hh][ll+DIMNS][kk+DIMNS] = ZERO;
                     smx.sh[hh][kk+DIMNS][ll+DIMNS] = ZERO;
                     kk++ ;
                  };
                  break;
            
                 case 2:
                  kk = null; do
                  {
                     fscanf( smatrix, scformat, ptr );
                     smx.sh[hh][ll][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr );
                     smx.sh[hh][ll+DIMNS][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr );
                     smx.sh[hh][ll+DIMNS][kk+DIMNS] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     smx.sh[hh][ll][kk+DIMNS] =\
                        SMX_SCMPR*( double )( ll == kk );

                     if( kk != ll ) 
                     {
                        smx.sh[hh][kk][ll] = smx.sh[hh][ll][kk];
                        smx.sh[hh][kk+DIMNS][ll] = smx.sh[hh][ll+DIMNS][kk];
                        smx.sh[hh][kk][ll+DIMNS] = smx.sh[hh][ll][kk+DIMNS];
                        smx.sh[hh][kk+DIMNS][ll+DIMNS] =\
                           smx.sh[hh][ll+DIMNS][kk+DIMNS];
                     };
                  } while (( ++kk ) <= ll );
                  break;

                 default:
                  kk = null; do
                  {
                     fscanf( smatrix, scformat, ptr );
                     smx.sh[hh][ll][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr );
                     smx.sh[hh][ll+DIMNS][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr );
                     smx.sh[hh][ll+DIMNS][kk+DIMNS] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     smx.sh[hh][ll][kk+DIMNS] =\
                        SMX_SCMPR*( double)( ll == kk );
                  } while (( ++kk ) < DIMNS );
                  break;
               };/* end switch (pp) */
            }/* end if state->dmn == 't' */
            else if ( state->dmn == 'f' )
            {
               switch (pp)
               {
                 case 1:
                  fscanf( smatrix, scformat, ptr ); /* real part             */
                  smx.sh[hh][ll][ll] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  fscanf( smatrix, scformat, ptr ); /* imaginary part        */
                  smx.sh[hh][ll][ll+DIMNS] =\
                     SMX_SCMPR*strtod( ptr, endp );

                  kk = null;                   /* off-diagonal elements */
                  while( kk < ll )             /* set to ZERO */
                  {
                     smx.sh[hh][ll][kk] = ZERO;
                     smx.sh[hh][kk][ll] = ZERO;
                     smx.sh[hh][ll][kk+DIMNS] = ZERO;
                     smx.sh[hh][kk][ll+DIMNS] = ZERO;
                     kk++ ;
                  };
                  break;

                 case 2:
                  kk = null; do
                  {
                     fscanf( smatrix, scformat, ptr ); /* real part */
                     smx.sh[hh][ll][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr ); /* imaginary part */
                     smx.sh[hh][ll][kk+DIMNS] =\
                        SMX_SCMPR*strtod( ptr, endp );
                     if( kk != ll )
                     {
                        smx.sh[hh][kk][ll] =\
                           smx.sh[hh][ll][kk];
                        smx.sh[hh][kk][ll+DIMNS] =\
                           smx.sh[hh][ll][kk+DIMNS];
                     };
                  } while (( ++kk ) <= ll );
                  break;

                 default:
                  kk = null ; do
                  {
                     fscanf( smatrix, scformat, ptr ); /* real part */
                     smx.sh[hh][ll][kk] =\
                        SMX_SCMPR*strtod( ptr, endp );

                     fscanf( smatrix, scformat, ptr ); /* imaginary part */
                     smx.sh[hh][ll][kk+DIMNS] =\
                        SMX_SCMPR*strtod( ptr, endp );
                  } while (( ++kk ) < DIMNS );
                  break;
               }; /* end switch ... */
            }; /* end if ( state->dmn == 'f' ) */
         } while (( ++ll ) < DIMNS );
/*............................................................................*/
         if ( smx.mtyp[hh] == 'g' )
         {
/* processing gyromagnetic scattering parameters [ if smx.mtyp == 'g' ]*/

            smx.mm[hh] = smx.gm;

            if ( smx.gm == null )
            {
/*............................................................................*/
/* The following function exptum(*) conditionally returnes matrices exu.er[][],
   ... in a way dependent on g-current integration macro SMX_GYMGCRR :   

   If one [ mesh independent ] set of gyromagnetic parameters, gyr.hg, gyr.ms,
   gyr.tg, .... , is given only one initial call of function 'exptum' is     
   sufficient.                                                              
   Otherwise, for mesh dependent parameters, this must be modified by calling 
   'exptum' each time new mesh parameters are entered.                        */
/*............................................................................*/
/* returns matrices gyr.um[][], gyr.vm[][], gyr.eum[][], gyr.qm[][]           */
                                               /*.............................*/
               ind = exptum( );               /*                              */
/*..........................................*/
            };
/*............................................................................*/
/* returns matrices smx.eum[][], smx.dm[][] and smx.sgm[][]                   */
                                           /*.................................*/
            ind = vmbinv( hh );           /*                                  */
/*......................................*/
            smx.gm++;
         };/* end if ...mtyp == 'g'yromagnetic */
      };/* end if ( ...mtyp != 't'rivial ) */

/* check Hilbert norm of smx.se...  */

/* Hilbert_norm: */

      hnrmi = ZERO;

      if (( state->dmn == 't' )&&( smx.etyp[hh] != 't' ))
      { 
         for ( ll=null; ll<DIMNS; ll++ )
         {
            for ( kk=null; kk<DIMNS; kk++ )
            {
               jac.h.r[ll][kk] = ZERO;
               jac.h.i[ll][kk] = ZERO;
               for ( qq=null; qq<DIMNS; qq++ )
               {
                  jac.h.r[ll][kk] += ( smx.se[hh][ll+THREE][qq+THREE] *\
                                       smx.se[hh][kk+THREE][qq+THREE] );
               };
            };
         };
/*............................................................................*/
/* compute Hilbert [spectral ] norm jac.jn: */ 

         ( jev->rank ) = DIMNS;
/*............................................................................*/
         jev = jacobi( jev );          /*                                     */
/*...................................*/
         hnrmi = jac.hn;

      };/* end if ( ...etyp != 't'rivial ) */
        
      if (( state->dmn == 't' )
        &&( smx.mtyp[hh] != 't' ))
      {
         for ( ll=null; ll<DIMNS; ll++ )
         {
            for ( kk=null; kk<DIMNS; kk++ )
            {
               jac.h.r[ll][kk] = ZERO;
               jac.h.i[ll][kk] = ZERO;

               for ( qq=null; qq<DIMNS; qq++ )
               {
                  jac.h.r[ll][kk] += ( smx.sh[hh][ll+THREE][qq+THREE] *\
                                       smx.sh[hh][kk+THREE][qq+THREE] );
               };
            };
         };
/*............................................................................*/
/* compute Hilbert [spectral ] norm jac.jn: */ 

         ( jev->rank ) = DIMNS;
/*............................................................................*/
         jev = jacobi( jev );          /*                                     */
/*...................................*/
         if ( hnrmi < jac.hn ) 
            hnrmi = jac.hn;

      };/* end if ( ...mtyp != 't' ) */
/*............................................................................*/
/* Violated stability condition [ Spectral/Hilbert norm exceeds ONE ]: */ 

      hnrmi = sqrt( hnrmi );

      if ( SMX_ONE_PLUS < hnrmi)
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Pathological S-matrix no. %ld !!!", hh );
         fprintf( stderr, "\n [ Hilbert norm %.15e > 1 may yield unstable "
            "process.]\n ", hnrmi );
         strncpy( state->fcterr, __func__, SHS_SIZE );
         strncpy( state->errmsg, state->smx, SHS_SIZE );
         strncat( state->errmsg, ":unstable_S_matrix[Hilbert_norm>1] ",
            ( LGS_SIZE - SHS_SIZE ));
         CLSMX( null );

         return null; /* abnormal return */
      };
/*............................................................................*/
      if ( hnorm < hnrmi )  
         hnorm = hnrmi;

      fscanf( smatrix, scformat, ptr ); /* string "Maxwll-field-smx<nn>" */
                                         /* or "S-PARAMETER_LABELS...", e.g.*/
   }; /* while ( null == strncmp( cmp, ptr, SIX )); */

   smx.hh[null] = hh;
/*............................................................................*/
# if DSC_HCRMDE != 0 /* enter thermal or fluid 's-parameters' */

   strcpy( cmp, thermal );
   strcat( cmp, smtrix );

   hh = null;
   while ( null == strncmp( ptr, cmp, SIX )) /* cmp = "Thermal-smx", e.g.*/
   {  
      hh++;
      if ( HCSMX < hh )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many thermal s-parameter sets defined "
            "in DSC system %s !!!", smx.name );
         fprintf( stderr, "\n [ Maximum number is %ld = macro HCSMX "
            "fixed in file %s.", ( long ) HCSMX, "SOLVER.CONF" );
         fprintf( stderr, "\n - Change macro in compliance with "
            "memory resources.]\n " );
         strncpy( state->fcterr, __func__, SHS_SIZE );
         strncpy( state->errmsg, state->smx, SHS_SIZE );
         strncat( state->errmsg, ":too_many_s-matrices_defined_",
            ( LGS_SIZE - SHS_SIZE ));

         return null; /* abnormal return */
      };
/*............................................................................*/
      fscanf( smatrix, scformat, ttyp ); /* thermal or fluid type */
      ( hsp->ttyp[hh] ) = *ttyp;

      if ( *ttyp != 't' ) /* non-trivial thermal and/or fluid type */
      { 
         ( hsp->scs[hh] ) = null; /* electric [magnetic] losses: s=1[2] */

/* temperature updating coefficient */
         fscanf( smatrix, scformat, ptr );
         ( hsp->ct[hh] ) = strtod( ptr, endp );

/* cell volume [m^3]: */
         fscanf( smatrix, scformat, ptr );
         ( hsp->vol[hh] ) = strtod( ptr, endp );

/* heat conductivity [W/(m*K)]: */
         fscanf( smatrix, scformat, ptr );
         ( hsp->kh[hh] ) = strtod( ptr, endp );

/* heat capacity per VOLUME [J/(m^3*K) !!!]: */
         fscanf( smatrix, scformat, ptr );
         ( hsp->cv[hh] ) = strtod( ptr, endp );
/*............................................................................*/
/* [ exterior ] cell face vectors: */

         kk = null; do
         {
            ll = null; do
            {
                fscanf( smatrix, scformat, ptr );
                ( hsp->f[hh][kk][ll] ) = strtod( ptr, endp );
            } while(( ++ll ) < DIMNS );
         } while(( ++kk ) < FACES );
/*............................................................................*/
/* cell surface [m^2]: */
         ( hsp->sf[hh] ) = ZERO; 
         kk = null; do
         {
/* area of face kk [m^2]: */
            fscanf( smatrix, scformat, ptr );
            ( hsp->fm[hh][kk] ) = strtod( ptr, endp );
            ( hsp->sf[hh] ) += ( hsp->fm[hh][kk] );
         } while(( ++kk ) < FACES );
/*............................................................................*/
# if DSC_FLDMDE != 0
         if ( *ttyp == 'f' )
         { /* fluid type */
/* normalized flow updating coefficient [sec*m^3/Kg] */
            fscanf( smatrix, scformat, ptr );
            ( hsp->ft[hh] ) = strtod( ptr, endp );

/* mean density: */
            fscanf( smatrix, scformat, ptr );
            ( hsp->rm[hh] ) = strtod( ptr, endp );

/* mean temperature: */
            fscanf( smatrix, scformat, ptr );
            ( hsp->tm[hh] ) = strtod( ptr, endp );

/* mean expansion coefficient [1/K]: */
            fscanf( smatrix, scformat, ptr );
            ( hsp->bm[hh] ) = strtod( ptr, endp );

/* mean adiabatic compression coefficient [1/Pa]: */
            fscanf( smatrix, scformat, ptr );
            ( hsp->cm[hh] ) = strtod( ptr, endp );

/* dynamic viscosity [Kg/(m*sec)]: */
            fscanf( smatrix, scformat, ptr );
            ( hsp->ny[hh] ) = strtod( ptr, endp );

/* constant q1 = Cp/Cv - 1: */
            fscanf( smatrix, scformat, ptr );
            ( hsp->q1[hh] ) = strtod( ptr, endp );

/* dissipation time constant td: */
            fscanf( smatrix, scformat, ptr );
            ( hsp->td[hh] ) = strtod( ptr, endp );

	    if (( SMX_TRIVIAL < ( hsp->td[hh] ))
               &&(( hsp->td[hh] ) < HUGE_VALF ))
/* dissipation factor */
	       ( hsp->dc[hh] ) = 1. - ( state->hcdt )/( hsp->td[hh] );
            else
	       ( hsp->dc[hh] ) = 1.;
/* no dissipation [ may yield unstable proc ] */
/*............................................................................*/
# if (( TURBMOD == 1 )\
    ||( TURBMOD == 2 )) /* Prandtl turbulence models */
            fscanf( smatrix, scformat, ptr );
            ( hsp->rl2[hh] ) = strtod( ptr, endp );
            ( hsp->rl2[hh] ) *= ( hsp->rl2[hh] );
            ( hsp->rl2[hh] ) *= ( hsp->rm[hh] );
# endif
/*............................................................................*/
/* gravitation: */
            ll = null; do
            {
               fscanf( smatrix, scformat, ptr );
               ( hsp->gr[hh][ll] ) = strtod( ptr, endp );
            } while (( ++ll ) < DIMNS );
	    
/* pressure gradient: */
            ll = null; do
            {
               fscanf( smatrix, scformat, ptr );
               ( hsp->gp[hh][ll] ) = strtod( ptr, endp );
            } while (( ++ll ) < DIMNS );
/*............................................................................*/
/* [ exterior ] cell face vectors: */
/* (1) *//*
            kk = null; do
            {
               ll = null; do
               {
                   fscanf( smatrix, scformat, ptr );
                   ( hsp->f[hh][kk][ll] ) = strtod( ptr, endp );
               } while(( ++ll ) < DIMNS );
            } while(( ++kk ) < FACES );
*/
/*............................................................................*/
/* matrix adj(B): */
/*
            kk = null; do
            {
               ll = null; do
               {
                  fscanf( smatrix, scformat, ptr );
                  ( hsp->b[hh][kk][ll] ) = strtod( ptr, endp );
               } while (( ++ll ) < DIMNS );
            } while (( ++kk ) < DIMNS );
*/
         }; /* endif ( ttyp -> "fluid" ) */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
/* matrix adj(B^-1): */

         kk = null; do
         {
            ll = null; do
            {
               fscanf( smatrix, scformat, ptr );
               ( hsp->bi[hh][kk][ll] ) = strtod( ptr, endp );
/*............................................................................*/
# if DSC_FLDMDE != 0
/*
               ii = null; do
               {
                  ( hsp->pb[hh][ll][kk][ii] ) =\
                     ( hsp->bi[hh][kk][ll] )*( hsp->b[hh][ll][ii] );
               } while (( ++ii ) < DIMNS );
*/
# endif /* DSC_FLDMDE != 0 */
            } while (( ++ll ) < DIMNS );
         } while (( ++kk ) < DIMNS );
/*............................................................................*/
/* form vectors;  ( F )^T * adj(B^-1 )[kk]: */

         kk = null; do
         {
            ll = null; do
            {
               fscanf( smatrix, scformat, ptr );
               ( hsp->s[hh][kk][ll] ) = strtod( ptr, endp );
            } while (( ++ll ) < DIMNS );
         } while (( ++kk ) < FACES );
/*............................................................................*/
         if (( null != strrchr( ttyp, (int) 'E' ))
           &&( null != strrchr( ttyp, (int) 'M' )))
         { /* lossy Maxwell field -> heat sources */

            ( hsp->scs[hh] ) = THREE;

            fscanf( smatrix, scformat, ptr );
            ( hsp->ke[hh] ) = strtod( ptr, endp );
            fscanf( smatrix, scformat, ptr );
            ( hsp->km[hh] ) = strtod( ptr, endp );
         }
	 else if ( null != strrchr( ttyp, (int) 'M' ))
         {
            ( hsp->scs[hh] ) = TWO;

            fscanf( smatrix, scformat, ptr );
            ( hsp->km[hh] ) = strtod( ptr, endp );
         }
	 else if ( null != strrchr( ttyp, (int) 'E' ))
         {
            ( hsp->scs[hh] ) = ONE;

            fscanf( smatrix, scformat, ptr );
            ( hsp->ke[hh] ) = strtod( ptr, endp );
         }
         else
            ( hsp->scs[hh] ) = null; /* no electric or magnetic losses */
      };/* end if not ( ttyp -> "trivial" ) */
/*............................................................................*/
      fscanf( smatrix, scformat, ptr ); /* string "crr-smx" */
                                         /* or ">>-S-PARAMETER_LABELS->>"  */
   };
   ( hsp->hh[null] ) = hh;
# endif /* DSC_HCRMDE != 0 */
/*...........................................................................*/
/* enter s-matrix labels */

   if ( null != strncmp( ptr, smxidx, THREE ))
   {
      fprintf( stderr, "\n\n Error message from function %s :", __func__ );
      fprintf( stderr, "\n Corrupted s-matrix file %s,", state->smx );
      fprintf( stderr, "\n different FORMER.C and SOLVER.C releases," );
      fprintf( stderr, "\n or incompatible compiler options !!! " );
      strncpy( state->fcterr, __func__, SHS_SIZE ); 
      strncpy( state->errmsg, state->smx, SHS_SIZE );
      strncat( state->errmsg, ":corrupted_file",
         ( LGS_SIZE - SHS_SIZE - TEN ));

      return null; /* abnormal return */
   };
/*...........................................................................*/
   fscanf( smatrix, scformat, ptr );

   if ( null != strncmp( ptr, "cell", THREE ))
   {
      fprintf( stderr, "\n\n Error message from function %s :", __func__ );
      fprintf( stderr, "\n Corrupted s-matrix file %s,", state->smx );
      fprintf( stderr, "\n different FORMER.C and SOLVER.C releases," );
      fprintf( stderr, "\n or incompatible compiler options !!! " );
      strncpy( state->fcterr, __func__, SHS_SIZE ); 
      strncpy( state->errmsg, state->smx, SHS_SIZE );
      strncat( state->errmsg, ":corrupted_file",
         ( LGS_SIZE - SHS_SIZE - TEN ));

      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( smatrix, scformat, ptr ); /* string "electric" */
/*............................................................................*/
# if DSC_HCRMDE != 0 
   fscanf( smatrix, scformat, ptr ); /* string "thermal/fluid" */
# if DSC_FLDMDE != 0
   fscanf( smatrix, scformat, ptr ); /* string "[fluid_connected_component]" */
/* hsp->cnn[m]: connected component [index] of cell m */
   ( hsp->cnn[null] ) = null;
# endif
# endif
/*............................................................................*/
/* Now we are running through the cells [ indexed ii ]: */

   ii = icll;
   while( ii <= fcll )
   {
      if ( NODES < ii )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many cells defined "
            "in DSC system %s !!!", smx.name );
         fprintf( stderr, "\n [ Maximum number is %ld = macro NODES "
            "fixed in file '%s'.", ( long ) NODES, "SOLVER.CONF" );
         fprintf( stderr, "\n - Change macro in compliance with "
            "memory resources.]\n " );
         strncpy( state->fcterr, __func__, SHS_SIZE );
         strncpy( state->errmsg, state->smx, SHS_SIZE );
         strncat( state->errmsg, ":too_many_nodes_defined_",
            ( LGS_SIZE - SHS_SIZE ));  

         return null; /* abnormal return */
      };
/*............................................................................*/
      fscanf( smatrix, scformat, ptr );
      nn = strtol( ptr, endp, DEC );

      if ( nn != ii )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n S-parameters file %s corrupted in "
            "cell number %ld !!!\n ", state->smx, ii );
         strncpy( state->fcterr, __func__, SHS_SIZE ); 
         strncpy( state->errmsg, state->smx, SHS_SIZE );
         strncat( state->errmsg, ":corrupted_file/cell",
            ( LGS_SIZE - SHS_SIZE - TEN ));

         return null; /* abnormal return */
      };
/*............................................................................*/
      fscanf( smatrix, scformat, ptr );
      smx.hh[ii] = strtol( ptr, endp, DEC );
/*...........................................................................*/
# if DSC_HCRMDE != 0 /* enter thermal [ or fluid ] node s-matrix label */
      fscanf( smatrix, scformat, ptr );
/* s-parameter index, cell ii */
      ( hsp->hh[ii] ) = strtol( ptr, endp, DEC );
/*...........................................................................*/
# if DSC_FLDMDE != 0 /* fluid connected component */

      if (( hsp->ttyp[( hsp->hh[ii] )] ) == 'f' )
      {
         fscanf( smatrix, scformat, ptr );
/* fluid [ connected ] component index */
         ( hsp->cnn[ii] ) = strtol( ptr, endp, DEC );

         if ( NFCNN < ( hsp->cnn[ii] ))
         {
            fprintf( stderr, "\n\n Error message from function %s:", __func__ );
            fprintf( stderr, "\n Too many fluid connected components defined "
               "in DSC system %s !!!", smx.name );
            fprintf( stderr, "\n [ Maximum number is %ld = macro NFCNN "
               "fixed in file '%s'.", ( long ) NFCNN, "SOLVER.CONF" );
            fprintf( stderr, "\n - Change macro in compliance with "
               "memory resources.]\n " );
            strncpy( state->fcterr, __func__, SHS_SIZE );
            strncpy( state->errmsg, state->smx, SHS_SIZE );
            strncat( state->errmsg, ":too_many_fluid_components_",
               ( LGS_SIZE - SHS_SIZE ));

            return null; /* abnormal return */
         }
         else if (( hsp->cnn[ii] ) <= null )
         {
            fprintf( stderr, "\n\n Error message from function %s:",
               __func__ );
            fprintf( stderr, "\n Illegal fluid connected component "
               "index cnn = %d <= null", ( hsp->cnn[ii] ));
            fprintf( stderr, "\n defined in cell %ld, s-parameter set %ld !!!",
               ii, ( hsp->hh[ii] ));
            fprintf( stderr, "\n [ legal: 0 < cnn <= %d "
               "( = macro NFCNN in file SOLVER.CONF )]\n\n ", NFCNN );
            return null;
         };

         ( hsp->volume[( hsp->cnn[ii] )] ) += ( hsp->vol[( hsp->hh[ii] )] );
         ( hsp->tmean[( hsp->cnn[ii] )] ) += \
	 (( hsp->vol[( hsp->hh[ii] )] )*( hsp->tm[( hsp->hh[ii] )] ));

/* eventually adjust number of */
/* fluid [connected] components: */

         if (( hsp->cnn[null] ) < ( hsp->cnn[ii] ))
            ( hsp->cnn[null] ) = ( hsp->cnn[ii] );
/*...........................................................................*/
# if CMPRSSBL != 0    
/* set fluid densities of cell[ii] to initial */
/* mean density in pertinent s-parameter set: */

         ( hsp->mass[( hsp->cnn[ii] )] ) += (( hsp->vol[( hsp->hh[ii] )] )*\
                                             ( hsp->rm[( hsp->hh[ii] )] ));
         qq = null; do
         {
           hcr[qq].rn[ii] = ( hsp->rm[hsp->hh[ii]] );
           pp = null; do
           {
              hcr[qq].rf[ii][pp] = hcr[qq].rn[ii];
           } while(( ++pp ) < FACES );
         } while(( ++qq ) < DSC_HCRMDE );
# endif
/*...........................................................................*/
      } /* end if (( hsp->ttyp[smx.hh[ii]] ) == 'f'luid type ) */
      else /* if (( hsp->ttyp[hsp->hh[ii]] ) != 'f' ) */
         ( hsp->cnn[ii] ) = null;
# endif /* DSC_FLDMDE != 0 */
/*...........................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*...........................................................................*/
      ii++ ;
   }; /* while( ii <= fcll ); [ next cell index ] */
/*...........................................................................*/
# if DSC_HCRMDE != 0
/*...........................................................................*/
# if DSC_FLDMDE != 0
/* initial mean temperatures: */ 

   ii = null; do
   {
      if ( ZERO < ( hsp->volume[ii] ))
         ( hsp->tmean[ii] ) /= ( hsp->volume[ii] );
   } while (( ii++ ) < ( hsp->cnn[null] ));
# endif
/*...........................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*...........................................................................*/
/* close s-matrix file: */

   CLSMX( state->rmfle );

   disp: ;
/*............................................................................*/
# if SMX_DISP == 1

   if ( state->dmn == 't' )
   {
      fprintf( display, "\n [ stability check o.k. :"
                        " Hilbert norm = %.13e <= 1 ]", hnorm );
      fprintf( display, "\n\n DSC time step..........................:"
                        " %.15e sec", state->dt );
   }
   else
   {
      fprintf( display, "\n\n DSC phase shift........................:"
                        " %.15e rad", state->ph );
   };
   ii = fcll - icll + ONE;
   fprintf( display, "\n number of nodes........................:"
		     " %ld", ii );

   fprintf( display, "\n --------------------------------------"
                      "----------------------------------------" );
# endif 

   return ONE; /* normal return */
}
/*============================================================================*/

/*******************************************************************************
*                                                                              *
*   The following function 'exptue(*)' returns matrices:                       *
*                                                                              *
*                gyr.ue[][] = U     [ sec^-1 ]                                 *
*                gyr.ve[][] = V,    [ m^-1   ]                                 *
*               gyr.eue[][] = E(1)  [ dimensionless ] , E(x) = exp(x*dt*U)     *
*                                                                              *
*   and , conditionally,                                                       *
*                                                                              *
*                gyr.qe[][] = (U^-1)*( E(1) - Id )/4.  [ sec ]                 *
*   or           gyr.qe[][] = dt * Id / 4.             [ sec ]                 *
*                                                                              *
*   if  SMX_GYELCRR = null or  SMX_GYELCRR != null,  respectively.             *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

short exptue( void )
{
/* allusions: */
/*
   extern struct tlmsmx smx;
   extern struct expmtx exu;
   extern GAUSS_JRD gss;
   extern GYRMTX gyr;
*/
/* declarations: */ 

   static struct solverstat
     *state = &solver;

   static GAUSS_JRD
     *gjp = &gss;

   static double  
      xx = ZERO,
      norm = ZERO;

   static short 
      ii = null, 
      jj = null,
      kk = null,
      nn = null;
/*
   static char 
      *ptr;
*/
   static signed char 
      nsg = null,
      sgn = null;

   double sqrt( double x );

   int expmtx( short rank );

   GAUSS_JRD 
      *gssjrd( GAUSS_JRD *gjp );
/*............................................................................*/

   if ( SMX_SNGLBND * state->dt < gyr.tr ) /* non-singular U */
      nsg = ONE;
   else
      nsg = null;

   norm = ELECTR_CHRGE*gyr.no;
   norm *= norm;

   for ( ii=null ; ii<DIMNS ; ii++ )
   {
      norm += ( gyr.mi[ii] * gyr.mi[ii] );
   };
   norm = sqrt(norm);
                                          
   for ( ii=null; ii<DIMNS; ii++ )
   {                                        
      if ( nsg == ONE ) /* ( SMX_SNGLBND * state->dt ) < gyr.tr: non singular U */
      {
         gyr.ue[ii][ii] = -1. / gyr.tr;
      }
      else 
      {
         gyr.ue[ii][ii] = ZERO;
      };

      gyr.ve[ii][ii] = ELECTR_CHRGE * gyr.no * gyr.ec;

# if SMX_GYELCRR == null
      ( gjp->mr[ii][ii] ) = state->dt * gyr.ue[ii][ii];
      ( gjp->mi[ii][ii] ) = ZERO; 
      exu.ur[ii][ii] = state->dt * gyr.ue[ii][ii];
# else /* if SMX_GYELCRR != null */
      exu.ur[ii][ii] = state->dt * gyr.ue[ii][ii];
# endif
      exu.ui[ii][ii] = ZERO;

      sgn = ONE;                           
      for ( jj=ONE; jj<DIMNS; jj++ )       
      {                              
         sgn *= -ONE;

         kk = (( ii+jj ) % THREE );
         nn = (( ii+TWO*jj ) % THREE );

         gyr.ue[ii][kk] = - sgn * gyr.ec * gyr.mi[nn]; /* [ sec^-1 ] */

         gyr.ve[ii][kk] = ZERO; 
          
         ( gjp->mr[ii][kk] ) = state->dt * gyr.ue[ii][kk]; /* [ dimensionless ] */
         ( gjp->mi[ii][kk] ) = ZERO;
         exu.ur[ii][kk] = state->dt * gyr.ue[ii][kk];      /* [ dimensionless ] */
         exu.ui[ii][kk] = ZERO;
      };
   };
/*............................................................................*/
   expmtx( DIMNS ); /* returns exponential matrix exp(exu.ur+i*exu.ui) */
/*......................*/

# if SMX_GYELCRR == null

   if ( nsg == ONE ) /* ( SMX_SNGLBND*state->dt ) < gyr.tr : non-singular U */
   {
      ( gjp ->rank ) = DIMNS;
      ( gjp->opt ) = 'i'; /* option: 'i'nverse */
/*............................................................................*/
      gjp = gssjrd( gjp );    /* matrix inversion ( dt*U )^-1                 */
/*..........................*/

      for ( ii=null; ii<DIMNS; ii++ )
      {
         for ( jj=null; jj<DIMNS; jj++ )
         {
            gyr.eue[ii][jj] = exu.er[ii][jj];
            gyr.qe[ii][jj]  = ZERO; 

            for ( kk=null; kk<DIMNS; kk++ ) 
            {
               xx = ( exu.er[ii][kk] - ( double )( ii == kk )); 
               gyr.qe[ii][jj] += ( xx * ( gjp->zr[kk][jj] ));
            };

            gyr.qe[ii][jj] *= ( QUART * state->dt );
         }; 
      };
   }   
   else /* singular U: gyr.tr <= ( SMX_SNGLBND * state->dt ) */
   {
      for ( ii=null; ii<DIMNS; ii++ )
      {
         for ( jj=null; jj<DIMNS; jj++ )
         {
            gyr.eue[ii][jj] = exu.er[ii][jj];
            gyr.qe[ii][jj]  = ( double )( ii == jj ) * QUART * state->dt; 
         };
      };
   }; 

# else /* if SMX_GYELCRR != null */

   for ( ii=null; ii<DIMNS; ii++ )
   {
      for ( jj=null; jj<DIMNS; jj++ )
      {
         gyr.eue[ii][jj] = exu.er[ii][jj];
         gyr.qe[ii][jj]  = ( double )( ii == jj ) * QUART * state->dt;
      };
   };

# endif

   return ONE;
}
/*============================================================================*/

/*******************************************************************************
*                                                                              *
*   The following function 'exptum(*)' returns matrices:                       *
*                                                                              *
*                gyr.um[][] = U     [ sec^-1 ]                                 *
*                gyr.vm[][] = c1*V, [ m^-1   ]        , c1   = MY_VAC_*smx.adm *
*               gyr.eum[][] = E(1)  [ dimensionless ] , E(x) = exp(x*dt*U)     *
*                                                                              *
*   and , conditionally,                                                       *
*                                                                              *
*                gyr.qm[][] = E(1)  [ dimensionless ] ,                        *
*   or           gyr.qm[][] = ((dt*U)^-1)*( E(1) - Id )*E(1/2)                 *
*                                                                              *
*   if  SMX_GYMGCRR != ONE  or  SMX_GYMGCRR == ONE, respectively.              *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

short exptum( void )
{
/* allusions: */
/*
   extern struct tlmsmx smx;
   extern struct expmtx exu;
   extern GYRMTX gyr;

# if SMX_GYMGCRR == ONE
   extern GAUSS_JRD gss;
# endif
 */
/* declarations: */ 

   static struct solverstat
     *state = &solver;

# if SMX_GYMGCRR == ONE
   static GAUSS_JRD *gjp = &gss;
# endif

   static double 
      xx = ZERO,
      norm = ZERO;

# if SMX_GYMGCRR == ONE
   static double 
      yy = ZERO;
# endif

   static short 
      ii = null, 
      jj = null,
      kk = null,
      nn = null,
   /*   ind = null, */
      ratio = null;
/*     
   static char 
      *ptr; 
*/
   static signed char 
      sgn = null;

   double sqrt( double x );

   int expmtx( short rank );

   GAUSS_JRD 
      *gssjrd( GAUSS_JRD *gjp );
/*............................................................................*/

   if (( SMX_SNGLBND*state->dt ) < gyr.tg ) /* case: non-singular U */
      ratio = TWO;
   else
      ratio = ONE;

   norm = ZERO;

   for ( ii=null; ii<DIMNS; ii++ )
   {
      norm += ( gyr.ms[ii] * gyr.ms[ii] ) ;
   };
   norm = sqrt(norm);

   xx = gyr.ld * gyr.mc;    /* gyromagnetic factor 'gamma' */
                            /* dimension:   [ m/(A*sec) ]  */
                                          
   for ( ii=null; ii<DIMNS; ii++ )
   {                                        
      if ( ratio == TWO ) /* ( SMX_SNGLBND*state->dt ) < gyr.tg: non singular U */
      {
         gyr.um[ii][ii] = -1. / gyr.tg;
      }
      else 
      {
         gyr.um[ii][ii] = ZERO;
      };

      gyr.vm[ii][ii] = ZERO;

# if SMX_GYMGCRR == ONE
      ( gjp->mr[ii][ii] ) = state->dt * gyr.um[ii][ii];
      ( gjp->mi[ii][ii] ) = ZERO; 
      exu.ur[ii][ii] = state->dt * gyr.um[ii][ii] / ratio;
# else /* if SMX_GYMGCRR != ONE */
      exu.ur[ii][ii] = state->dt * gyr.um[ii][ii];
# endif
      exu.ui[ii][ii] = ZERO;

      sgn = ONE;                           
      for ( jj=ONE; jj<DIMNS ; jj++ )       
      {                              
         sgn *= -ONE;

         kk = (( ii+jj ) % THREE );
         nn = (( ii+TWO*jj ) % THREE );

         gyr.um[ii][kk] = - sgn * xx * gyr.hg[nn]; /* [ sec^-1 ] */

         if ( SMX_GYROBND < norm )
         {
            gyr.vm[ii][kk] = sgn * smx.adm * xx * gyr.ms[nn]; /* [ m^-1 ] */
         }
         else
         {
            gyr.vm[ii][kk] = ZERO;
         };

# if SMX_GYMGCRR == ONE
         ( gjp->mr[ii][kk] ) = state->dt * gyr.um[ii][kk]; /*[dimensionless]*/
         ( gjp->mi[ii][kk] ) = ZERO;
         exu.ur[ii][kk] = state->dt * gyr.um[ii][kk] / ratio; /*[dimensionless]*/
# else /* if SMX_GYMGCRR != ONE */
         exu.ur[ii][kk] = state->dt * gyr.um[ii][kk];      
# endif
         exu.ui[ii][kk] = ZERO;
      };
   };
/*............................................................................*/
   expmtx( DIMNS ); /* returns exponential matrix exp(exu.ur+i*exu.ui) */
/*............................................................................*/
# if SMX_GYMGCRR == ONE

   if ( ratio == TWO ) /* ( SMX_SNGLBND*state->dt ) < gyr.tg: non-singular U */
   {
      ( gjp ->rank ) = DIMNS;
      ( gjp->opt ) = 'i'; /* option: 'i'nverse */
/*............................................................................*/
      gjp = gssjrd( gjp );    /* matrix inversion ( dt*U )^-1                 */
/*..........................*/

      for ( ii=null; ii<DIMNS; ii++ )
      {
         for ( jj=null; jj<DIMNS; jj++ )
         {
            xx = ZERO; 
            for ( kk=null; kk<DIMNS; kk++ )
            {
               xx += ( exu.er[ii][kk] * exu.er[kk][jj] );
            };
            gyr.eum[ii][jj] = xx;
         };
         for ( jj=null; jj<DIMNS; jj++ )
         {
            gyr.qm[ii][jj] = ZERO; 
            for ( kk=null; kk<DIMNS; kk++ ) 
            {
               xx = ZERO;
               for ( nn=null; nn<DIMNS; nn++ ) 
               {
                  xx +=  ( (gjp->zr[kk][nn]) * exu.er[nn][jj] );
               };
               yy = ( gyr.eum[ii][kk] - ( double )( ii == kk ) );
               gyr.qm[ii][jj] +=  ( xx * yy );
            };
         }; 
      };
   }   
   else /* singular U:  gyr.tg <= ( SMX_SNGLBND * state->dt ) */
   {
      for ( ii=null; ii<DIMNS; ii++ )
      {
         for ( jj=null; jj<DIMNS; jj++ )
         {
            xx = exu.er[ii][jj];
            gyr.eum[ii][jj] = xx;
            gyr.qm[ii][jj] = xx;
         };
      };
   }; 

# else /* if SMX_GYMGCRR != ONE */

   for ( ii=null; ii<DIMNS; ii++ )
   {
      for ( jj=null; jj<DIMNS; jj++ )
      {
         xx = exu.er[ii][jj];
         gyr.eum[ii][jj] = xx;
         gyr.qm[ii][jj] = xx; 
      };
   };

# endif

   return ONE;
}

/*******************************************************************************
*                                                                              *
*   The following function 'vebinv(*)' returns matrices:                       *
*                                                                              *
*    smx.eue[][] = exp(dt*U)     [ exp(dt*U) = gyr.eue[][]                   ] *
*     smx.de[][] = -((T)^-1)*A*Q [ A = smx.a[][], Q = gyr.qe[][]             ] *
*    smx.sge[][] = ( V*B^-1 )/4. [ B^-1 = c1*(A*) ; c1 = 1./sqrt(det(A))     ] *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

short vebinv( long hh )
{
/*............................................................................*/
/* allusions: */
/*
   extern struct tlmsmx smx;
   extern GAUSS_JRD gss;
   extern GYRMTX gyr;
   extern FORMTX frm;
*/
/* declarations: */ 

   static GAUSS_JRD
     *gjp = &gss;
/*
   static char 
      *ptr;
*/
   static short 
      ii = null,
      jj = null,
      nn = null;

   static double 
      xx = ZERO,
      yy = ZERO,
      qtn = ZERO;

   static double 
      mat[DIMNS][DIMNS] = {{ZERO}}; 
    
   double sqrt( double x );

   GAUSS_JRD
     *gssjrd( GAUSS_JRD *gjp );
/*----------------------------------------------------------------------------*/
/* - here starts the job: ----------------------------------------------------*/

/*  computes matrix  B^-1  from matrix A . In `parcel-twines' scheme:         */
/*  B^-1  =  A'*(sqrt(det(A)))^-1   ;  A = frm.aa[j][k] , A' = frm.aa[k][j].  */
/*  Also, returns smx.eum[smx.gm] = gyr.eum = exp(t*U)                        */

   for ( ii=null; ii<DIMNS; ii++ )
   {
      for ( jj=null; jj<DIMNS; jj++ )
      {
         ( gjp->mr[ii][jj] ) = frm.aa[ii][jj];
         ( gjp->mi[ii][jj] ) = ZERO;
      };
   };

   ( gjp->rank ) = DIMNS;
   ( gjp->opt ) = 'd'; /* option: 'd'eterminant */
/*............................................................................*/
   gjp = gssjrd( gjp );    /* returns ( gjp->dtr ) = determinant of A         */
/*.......................*/
    
   qtn = 4.*sqrt(( gjp->dtr ));

   for ( ii=null; ii< DIMNS; ii++ )
   {
      for ( jj=null; jj<DIMNS; jj++ )
      {
         smx.eue[smx.ge][ii][jj] = gyr.eue[ii][jj];

         xx = ZERO;
         yy = ZERO;

         for ( nn=null; nn<DIMNS; nn++ )
         {
            xx += ( frm.aa[ii][nn] * gyr.qe[nn][jj] );
            yy += ( gyr.ve[ii][nn] * frm.aa[jj][nn] );
         };
         mat[ii][jj] = xx;
         smx.sge[smx.ge][ii][jj] = SMX_GCMPR * yy / qtn;  /* = ( V*B-1 )/4. */
      };
   };

   for ( ii=null ; ii< DIMNS ; ii++ ) /* smx.de = - T^-1 * A * Q */
   {                                  /* [ (T)^-1 = (KE + Id3)/smx.adm ] */
      for ( jj=null ; jj<DIMNS ; jj++ )      
      {
         xx = ZERO;
         for ( nn=null; nn<DIMNS; nn++ )     
         {
            xx += ( (( double )( ii == nn )+smx.se[hh][ii][nn])*mat[nn][jj] );
         };
         smx.de[smx.ge][ii][jj] = - SMX_GCMPR * xx / smx.adm;
      };
   };

   return ONE;
}
/************************** end of function vebinv(*) *************************/

/*******************************************************************************
*                                                                              *
*   The following function 'vmbinv(*)' returns matrices:                       *
*                                                                              *
*    smx.eum[][] = exp(dt*U)     [ exp(dt*U) = gyr.eum[][]                   ] *
*     smx.dm[][] = ((yT)^-1)*A*Q [ A = smx.a[][], Q = gyr.qm[][] = Y/MY_VAC_ ] *
*    smx.sgm[][] = c1*V*B^-1     [ B^-1 = c2*(A*) ; c1 = MY_VAC_*smx.adm/4., ] *
*                                [                  c2 = 1./sqrt(det(A))     ] *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

short vmbinv( long hh )
{
/* allusions: */
/*
   extern struct tlmsmx smx;
   extern GAUSS_JRD gss; 
   extern GYRMTX gyr;
   extern FORMTX frm;
*/
/* declarations: */ 

   static GAUSS_JRD *gjp = &gss;
/*
   static char 
      *ptr;
*/
   static short 
      ii = null,
      jj = null,
      nn = null;

   static double 
      xx = ZERO,
      yy = ZERO,
      qtn = ZERO;

   static double 
      mat[DIMNS][DIMNS] = {{ZERO}}; 
    
   double sqrt( double x );

   GAUSS_JRD
     *gssjrd( GAUSS_JRD *gjp );
/*----------------------------------------------------------------------------*/
/* - here starts the job: ----------------------------------------------------*/

/*  computes matrix  B^-1  from matrix A . In `parcel-twines' scheme:         */
/*  B^-1  =  A'*(sqrt(det(A)))^-1   ;  A = frm.aa[j][k] , A' = frm.aa[k][j].  */
/*  Also, returns smx.eum[smx.gm] = gyr.eum = exp(t*U)                        */

   for ( ii=null; ii<DIMNS; ii++ )
   {
      for ( jj=null; jj<DIMNS; jj++ )
      {
         (gjp->mr[ii][jj]) = frm.aa[ii][jj];
         (gjp->mi[ii][jj]) = ZERO;
      };
   };

   ( gjp->rank ) = DIMNS;
   ( gjp->opt ) = 'd'; /* option: 'd'eterminant */
/*............................................................................*/
   gjp = gssjrd( gjp );    /* returns ( gjp->dtr ) = determinant of A         */
/*.......................*/

   qtn = 4.*sqrt(( gjp->dtr ));

   for ( ii=null; ii< DIMNS; ii++ )
   {
      for ( jj=null; jj<DIMNS; jj++ )
      {
         smx.eum[smx.gm][ii][jj] = gyr.eum[ii][jj];

         xx = ZERO;
         yy = ZERO;

         for ( nn=null; nn<DIMNS; nn++ )
         {
            xx += ( frm.aa[ii][nn] * gyr.qm[nn][jj] );
            yy += ( gyr.vm[ii][nn] * frm.aa[jj][nn] );
         };
         mat[ii][jj] = xx;
         smx.sgm[smx.gm][ii][jj] = SMX_GCMPR * yy / qtn;
      };
   };

   for ( ii=null; ii< DIMNS; ii++ )    /* smx.dm = ((yT)^-1 * A * Q ) */
   {                                   /* [ (yT)^-1 = (Id3 - KH) ]    */
      for ( jj=null; jj<DIMNS; jj++ )      
      {
         xx = ZERO;
         for ( nn=null; nn<DIMNS; nn++ )     
         {
            xx += ((( double )( ii == nn ) - smx.sh[hh][ii][nn] ) * \
               mat[nn][jj] );
         };
         smx.dm[smx.gm][ii][jj] = SMX_GCMPR * xx;
      };
   };

   return ONE;
}
/************************** end of function vmbinv(*) *************************/
# undef SMX_DISP
# undef SMX_ONE_PLUS
# undef SMX_ONE_MINUS
# undef SMX_GYMGCRR
# undef SMX_SCMPR
# undef SMX_GYROBND      
# undef SMX_SNGLBND
# undef SMX_TRIVIAL
# undef EPS_VAC
# undef MY_VAC_
# undef ELECTR_CHRGE
# undef ELECTR_MASS_
# undef BOHRs_MAGNTN
# undef PLANCKs_CNST  
# undef CLSMX
/**************** end of DSC S-matrix input function smtrix(*) ****************/
