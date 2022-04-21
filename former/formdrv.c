/* [ file: formdrv.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function formdrv(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   This function generates the DSC system input and output files              *
*   dsc.top[nn], dsc.smx[nn], dsc.bnd[nn], dsc.exc[nn], dsc.val[nn],           *
*   and dsc.batch, for the parameter transfer to program solver.c.             *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 22, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*----------------------------------------------------------------------------*/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
# include <unistd.h> /* system specification header, cf. sysconf( ), etc. */
# include <time.h>   /* cf. time( ), ctime( ), asctime( ), localtime( ), etc. */
/*----------------------------------------------------------------------------*/
# if _ISOC99_SOURCE == 1
   # include <fenv.h>
   # include <iso646.h>
# endif
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize this header for FORMER configuration: */
# include "./FORMER.CONF" /* Edit and customize this FORMER */
                                  /* configuration header */
/*----------------------------------------------------------------------------*/
/* A macro that gives a conveniently readable form to the string returned */
/* with [ the pointer of ] function ctime(&nseconds) */
# include "../tools/TIMEFORM.M"
/*----------------------------------------------------------------------------*/
# ifndef USE_NCURSES
   # define USE_NCURSES 1 /* =1 on BSD or Linux [GNU] systems, e.g.*/
# endif
/*----------------------------------------------------------------------------*/
/* the following [two] macros should be defined in "../CONFIG.H"              */
# ifndef DSC_ADJCLL      /* assign neighboring cell index top.mn[][k] so:    */
   # define DSC_ADJCLL 0 /* 0: to faces [ k is a face index; 0 <= k < 6 ]    */
# endif                   /* 1: to ports [ k is a port index; 0 <= k < 12 ]   */
/*----------------------------------------------------------------------------*/
# ifndef DSC_FCELBL      /* label neighboring face index top.fn[][] so:      */
   # define DSC_FCELBL 0 /* 0: unsigned, started with index null             */
# endif                   /* 1: unsigned, started with index ONE              */
                          /* 2: signed, started with index ONE                */
                          /* [ the negative/positive sign indicating opposed/ */
                          /*   aligned orientation to neighboring cell face ] */
/*----------------------------------------------------------------------------*/
# ifndef LNELNGTH
   # define LNELNGTH 78 /* maximum line length in structure files dsc.xyz */
# endif
/*----------------------------------------------------------------------------*/
/* Disable/enable electric boundary completeness check [0/1/2], 2 with return */
# ifndef FORM_BNDCMPL
   # define FORM_BNDCMPL 1
# endif
/*----------------------------------------------------------------------------*/
# ifndef NFCES
   # define NFCES 0       /* number of labelled faces [ can be 0 ] */
# endif
/*----------------------------------------------------------------------------*/
/* file identifiers [ names, pre- and suffixes etc. ]:                        */
# ifndef DSC_PRFX
   # define DSC_PRFX "dsc."/* common file prefix for SOLVER.C program files */
# endif
# define TOPOLOGY_FILE   "top"
# define S_MATRIX_FILE   "smx"
# define BOUNDARY_FILE   "bnd"
# define EXCITATION_FILE "exc"
# define EVALUATION_FILE "val"
# define SOLVER_BATCH_FILE "batch"
/*----------------------------------------------------------------------------*/
/* batch file limits [ for danse core program ]: */
# define DSC_JOBS  100 /* maximum number of DSC jobs per runtime */
			/* should coincede with same macro in */
			/* program solver [ and not excede 100 ] */
# define DSC_MXJBL 999 /* maximum job [ index ] */
			/* should coincede with same macro in */
			/* program solver [ and not excede 999 ] */
/*----------------------------------------------------------------------------*/
/* surface macros [ display instructions etc., without algorithmic impact ]   */

# define FORM_ITEMS  7 /* number of items in Elsy_menu [ except  0 = end ]    */
# define FORM_REPEAT 1 /* 1: check repeating s-parameters, set repeater flags */
# define FORM_DISP   1 /*  additional monitoring of intermediate results      */
                       /*  if defined CELL_POINTS = n:                        */
                       /*  print vertex points of mesh n                      */
# define LBL_CELLS   1 /* =1 : write cell labels in structure files           */
/*----------------------------------------------------------------------------*/
# ifndef VAL_INITLBL
   # define VAL_INITLBL 1 /* if val.ni < VAL_INITLBL the val.ni = VAL_INITLBL */
# endif
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
   # define BND_CSHAPE 1 /* call cshape(*) function for [ skin effect ] heat */ 
                         /* source definition on boundary faces */ 
   # ifndef DSC_INTLCE
      # define DSC_INTLCE 1 /* 0/1: separated/interlaced internal Maxwell */
   # endif                   /* field and heat current computation loops */
# endif /* DSC_HCRMDE != 0 */
/*----------------------------------------------------------------------------*/
/* precision in formdrv: */
# ifndef PRECISION
   # define PRECISION ( double )( 1.0e-15 )
# endif
/*----------------------------------------------------------------------------*/
# if FORM_REPEAT == 1
/* relative precision with that repeated parameters shall be distinguished */
# define RPT_PRECISION ( double )( 3.3e+01*PRECISION )
# endif
/*----------------------------------------------------------------------------*/
/* very small double in formdrv: */
# ifndef SMALL_VAL
   # define SMALL_VAL ( 1.0e-277 )
# endif
/*----------------------------------------------------------------------------*/
/* giant double in formdrv: */
# ifndef GIANT_VAL
   # define GIANT_VAL ( 1.0e+277 )
# endif
/*----------------------------------------------------------------------------*/
/* structures typedefs, etc.: */
# include "./formtp.h"
# include "./smxtyp.h"
# include "../tools/txctyp.h"
# include "../solver/dsptyp.h"
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
   # include "./hcrstp.h"
   # include "../tools/cshptp.h"
# endif
/*----------------------------------------------------------------------------*/
# ifndef EPS_VAC
   # define EPS_VAC  ( 8.8541878170e-12 ) /* vac. permittivity [A*sec/(V*m)] */
# endif
# ifndef MY_VAC_
   # define MY_VAC_  ( 1.2566370614e-06 ) /* "    permeability [V*sec/(A*m)] */
# endif
/*----------------------------------------------------------------------------*/
static FORMSTATE form = {null};
/*----------------------------------------------------------------------------*/
static TOPOLOGY top = {null};
static PARAMETERS par = {null};
static BOUNDARIES bnd = {null};
static EXCITATION exc = {null};
static EVALUATION val = {null};
/*----------------------------------------------------------------------------*/
static struct coordinates cor = {null};
static struct media med = {null};
static struct repeater rpf = {{null}};
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
/*
static CSHAPE csh = {null};
*/
static HCRSMX hsx = {null};
static struct repeater rph = {{null}};
# if DSC_FLDMDE != 0
static struct fconnect fcn = {{null}};
# endif
# endif /* DSC_HCRMDE != 0 */
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* 'my_terminal' configuration: */

   # include <ncurses.h>
   # include <curses.h>
   # include <term.h> /* terminal type header */

   static char *term; /* terminal type string */ 

   # define CLSCREEN {\
     printf( "%s", tgetstr( "cl", null )); \
   }

   # define PRBLDCLR(a) {\
     printf( "%s%s", tgetstr( "md", null ), (a)); /* bold clear output */ \
   }

   # define PRINVERS(a) {\
     printf( "%s%s", tgetstr( "mr", null ), (a)); /* inverse */ \
   }

   # define PRNORMAL(a) {\
     printf( "%s%s", tgetstr( "me", null ), (a)); /* back to normal output */ \
   }
# else
   # define CLSCREEN { ;\
   }

   # define PRBLDCLR(a) {\
     printf( "%s", (a));\
   }

   # define PRINVERS(a) {\
     printf( "%s", (a));\
   }

   # define PRNORMAL(a) {\
     printf( "%s", (a));\
   }
# endif
/*============================================================================*/

short formdrv( char *err )
{
/* allusions: */

/* extern TOPOLOGY top;   *//* topology structure: mesh gridding scheme */
/* extern PARAMETERS par; *//* parameter structure: coordnts. and media */
/* extern BOUNDARIES bnd; *//* boundary configuration structure */
/* extern EXCITATION exc; *//* excitation file configuration structure */
/* extern EVALUATION val; *//* evaluation file configuration structure */
/* extern FORMSTATE form; *//* [formdrv] menu state transfer structure */
/* extern struct coordinates cor; */
/* extern struct media med; */
/* extern struct repeater rpf; */
/* # if DSC_HCRMDE != 0 */
/* extern struct repeater rph; */
/* # if DSC_FLDMDE != 0 */
/* extern struct fconnect fcn; */
/* # endif */
/* # endif */

/* extern FORMSTATE form; *//* [formdrv] menu state transfer structure */

/* declarations: */

   static FORMSTATE *state = &form;

   static TOPOLOGY *tpt = &top;
   static PARAMETERS *ppt = &par;
   static BOUNDARIES *bpt = &bnd;
   static EXCITATION *ept = &exc;
   static EVALUATION *vpt = &val;

   static struct coordinates *cpt = &cor;
   static struct media *mpt = &med;
   static struct repeater *rep = &rpf;
   
   static TXCNSL *csp; /* cf. txcnsl(*) */
   static DSPLAY *dsp; /* cf. dsplay(*) */ 

   static S_MATRIX 
     *spt;
/*............................................................................*/
# if DSC_HCRMDE != 0
   static struct repeater *rhp;

# if DSC_FLDMDE != 0
   static struct fconnect *fcp;
# endif

   static HCRSMX 
      *hsp;

# if BND_CSHAPE == 1

   static CSHAPE
      *clp = NULL;

   static double
      admce = ZERO,
      ss = ZERO;

# endif /* BND_CSHAPE == 1 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* DSC system files */
/*
   static FILE
     *keyboard = stdin,
     *display = stdout;

   [ This didn't work with Redhat's version of egcc C compiler ]
*/
   static FILE
     *topology,
     *s_matrix, 
     *bndryfle,
     *excitfle,
     *evalfle,
     *batchfle;

/* log files: */

# if VERTX_LOG == 1
   FILE *vtxlog;
   static char *vtxfle = "vtx.log";
# endif
# if UNLNK_LOG == 1
   FILE *linklog;
   static char *lnkfle = "lnk.log";
# endif
# if EMWLL_LOG == 1
   FILE *wallog;
   static char *wllfle = "wll.log";
# endif
# if COORD_LOG == 1
   FILE *coordlog;
   static char *corfle = "crd.log";
# endif
# if CSIZE_LOG != 0
   FILE *csizelog;
   static char *cszfle = "csz.log";
# endif
# if MEDIA_LOG == 1
   FILE *medialog; 
   static char *medfle = "med.log";
# endif
# if SMTRX_LOG == 1
   FILE *smtrxlog;
   static char *smxfle = "smx.log";
# endif

# if CSIZE_LOG != 0
   static double 
      xx1 = ZERO,
      xx2 = ZERO,
      yy1 = ZERO,
      yy2 = ZERO,
      zz1 = ZERO,
      zz2 = ZERO,

      qdist = ZERO,
      csize = ZERO,
      mmsze = ZERO;
   static double
      mxsze[CSIZE_LOG] = {null};
   static long
      mxcll[CSIZE_LOG] = {null};
   static char
      vertx[CSIZE_LOG][TWO] = {{null}};
# endif /* CSIZE_LOG != 0 */

   static char 
      ptr[STS_SIZE] = {null},           
      fleptr[SHS_SIZE] = {null},
    **endp = null;

   static const char
      slsh = 47,  /* ASCII 'slash'      */
      bslsh = 92;  /* ASCII 'back slash' */
    
   static short
      llns = ONE,
      clns = 10,
      lbnd = 0,
      cbnd = 8,  /* <~ cbnd = clns - 2 */
      cc = null,
      fc = null,
      pp = null,
      qq = null,
      jj = null,
      kk = null,
      ll = null,
      idx = null,  
      idx1 = null,
      idx2 = null;

   static long 
      hh = null,
      ii = null,
      mm = null,
      mn = null,
      nn = null,
      ecl = null, /* electric current cell label */
      mcl = null, /* magnetic current cell label */
      node = null;

   static char 
      extyp = null;

   static const char 
/* 1st and 2nd ports pertinent to faces 0,...,5: */
      prt1[FACES] = { 7, 5, 11, 9, 3, 1 }, /* 1st port */
      prt2[FACES] = { 10, 8, 2, 0, 6, 4 }, /* 2nd port */
/* faces pertinent to ports 0,1,...,11: */
      fce[PORTS] = { 3, 5, 2, 4, 5, 1, 4, 0, 1, 3, 0, 2 };
/*............................................................................*/
# if DSC_HCRMDE != 0

   static short
      ff = null,
      fn = null;

   static char 
      exctp = null;
# endif
/*............................................................................*/
/* the following character strings MUST coincede with equally named strings */ 
/* in file smtrix.h [ cf. directory solver ] */ 

   static const char
     *smtrix = "-smx",
     *fields = "Maxwell_field",
     *smxptr = S_MATRIX_FILE,
     *smxidx = ">>-S-PARAMETER_LABELS->>",
     *smxpar = ">>-S-PARAMETERS->>",
     *fclbls = ">>-CELL_and_FACE_LABELS,_parameters,_etc->>";
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE == 0
   static const char
     *thermal = "Thermal";
# else
   static const char
     *thermal = "Thermal-fluid";
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* the DSC system file names [ and other strings ]: */

   static const char
     *prefix = DSC_PRFX,
     *sysptr = TOPOLOGY_FILE, 
     *bndptr = BOUNDARY_FILE, 
     *excptr = EXCITATION_FILE,
     *evlptr = EVALUATION_FILE,
     *btchpt = SOLVER_BATCH_FILE,
     *rmvfle = "remove_files",
     *rtnfle = "retain_files", 
     *spformat = "%-27s\n",
     *i_format = "%ld\n",
     *d_format = "%+.16e\n",
     *time_frm = "created: %24s";
/*............................................................................*/
/* batch file labels: */

   static short
      joblbl[DSC_JOBS] = {null},
      toplbl[DSC_JOBS] = {null},
      smxlbl[DSC_JOBS] = {null},
      bndlbl[DSC_JOBS] = {null},
      exclbl[DSC_JOBS] = {null};
/*............................................................................*/
/* time_t types: */

   time_t nseconds = null;
   time_t   *timer = null;
   time_t time( time_t *timer );

   static char 
      ctmptr[STS_SIZE] = {null}, 
      tmestr[STS_SIZE] = {null}; 
/*............................................................................*/
# if FORM_REPEAT == 1
/* smx mirror structure: */

   static S_MATRIX
      rpt = {null};
/*............................................................................*/
# if DSC_HCRMDE != 0
/* hcrsmx mirror structure: */

   static HCRSMX
      rhc = {null};
# endif
/*............................................................................*/

   static const double /* parameter distinguishing bounds for: */
      se_bound = RPT_PRECISION, /* E_field s-parameters */
      sm_bound = RPT_PRECISION, /* M_field s-parameters */
      no_bound = RPT_PRECISION,
      tr_bound = RPT_PRECISION,
      ld_bound = RPT_PRECISION,
      tg_bound = RPT_PRECISION,
      mi_bound = RPT_PRECISION,
      ms_bound = RPT_PRECISION,
      hg_bound = RPT_PRECISION,
      aa_bound = RPT_PRECISION; /* area matrix (A) */
/*............................................................................*/
# if DSC_HCRMDE != 0
   static const double /* relative parameter variation bounds: */
      hct_bound = RPT_PRECISION, /* heating coefficient ( hsp->ct ) */
      hke_bound = RPT_PRECISION, /* electr [loss] conductivity ( hsp->ke ) */
      hkm_bound = RPT_PRECISION, /* magnet [loss] conductivity ( hsp->km ) */
      hvl_bound = RPT_PRECISION, /* cell volume ( hsp->vol ) */
      hcf_bound = RPT_PRECISION, /* face vector bound */
      hcv_bound = RPT_PRECISION, /* heat capacity */
      hkh_bound = RPT_PRECISION, /* heat conductivity */
      hbi_bound = RPT_PRECISION, /* bound for matrix adj(B^-1) */
      hcs_bound = RPT_PRECISION; /* form vector bound */
/*............................................................................*/
# if DSC_FLDMDE != 0
   static const double /* parameter distinguishing bounds for: */
   /* hcb_bound = RPT_PRECISION, *//* bound for adjoint node matrix adj(B) */
      hft_bound = RPT_PRECISION, /* normalized fluid upd. coeff. ( hsp->ft ) */
      hrm_bound = RPT_PRECISION, /* mean density ( hsp->rm ) */
      htm_bound = RPT_PRECISION, /* mean temperature ( hsp->tm ) */
      hbm_bound = RPT_PRECISION, /* expansion coeff ( hsp->bm ) */
      hcm_bound = RPT_PRECISION, /* compression coeff ( hsp->cm ) */
      hny_bound = RPT_PRECISION, /* dynamic viscosity ( hsp->ny ) */
      hq1_bound = RPT_PRECISION, /* Cp/Cv - 1         ( hsp->q1 ) */
      htd_bound = RPT_PRECISION, /* dissipation time constant ( hsp->td ) */
/*............................................................................*/
/* Prandtl turbulence models */
# if (( TURBMOD == 1 ) \
    ||( TURBMOD == 2 ))
      hLL_bound = RPT_PRECISION, /* characteristic length ( hsp->LL ) */
# endif /* TURBMOD != null */
/*............................................................................*/
      hgr_bound = RPT_PRECISION, /* gravity vector ( hsp->gr ) */
      hgp_bound = RPT_PRECISION; /* pressure gradient ( hsp->gp ) */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
# endif /* FORM_REPEAT == 1 */
/*............................................................................*/
/* prototyping: */

# ifndef _CCBUG
   char 
      *strcpy( char *ptr1, const char *ptr2 ),  
      *strcat( char *ptr1, const char *ptr2 ),
      *strncat( char *ptr1, const char *ptr2, size_t n );
# endif

   char
      *lotos( long mm, char cc );

   char
      *dotos( double x, char precision, char *format );

   DSPLAY 
     *dsplay( DSPLAY *dsp );

   TXCNSL 
     *txcnsl( TXCNSL *csp );

   S_MATRIX 
      *smatrx( S_MATRIX *spt );
/*............................................................................*/
# if DSC_HCRMDE != 0
# if BND_CSHAPE == 1
   CSHAPE
      *cshape( CSHAPE *clp );
# endif

   HCRSMX
      *hcrsmx( HCRSMX *hsp );
# endif
/*............................................................................*/
   TOPOLOGY *
      linker( TOPOLOGY *tpt );

/* model dependent fuctions: */

   TOPOLOGY *
      systop( FORMSTATE *state );
   PARAMETERS *
      syssmx( FORMSTATE *state );
   BOUNDARIES *
      sysbnd( FORMSTATE *state );
   EXCITATION *
      sysexc( FORMSTATE *state );
   EVALUATION *
      sysval( FORMSTATE *state );
/*----------------------------------------------------------------------------*/
# if USE_NCURSES
/* get the terminal info: */

   term = ( char *) getenv( "TERM" ); /* get the terminal */

   kk = tgetent( null, term );

   if ( ONE != kk )
   {
      fprintf( stderr, "Error on getting the termcap info\n" ); 
      exit( EXIT_FAILURE );
   };
# endif
/*............................................................................*/
/* set buffer length = null: */
/*
   kk = setvbuf( keyboard, null, _IONBF, null );
   kk = setvbuf( display, null, _IONBF, null );
   [ This didn't work with Redhat's version of egcc C compiler ]
*/
   kk = setvbuf( stdin, null, _IONBF, null );
   kk = setvbuf( stdout, null, _IONBF, null ); 
/*............................................................................*/
/* memory allocations: */   
/* -- only static variables used -- */  
/*...........................................................................*/
/* initialize: */

   csp = txcnsl( null ); /* initialize the text console */
   dsp = dsplay( null ); /* initialize display [ running cursor function ] */

   strncpy(( ppt->domain ), "time_domain", 11 );

   ( ppt->cpt ) = cpt;
   ( ppt->mpt ) = mpt;
   ( ppt->rep ) = rep;
/*............................................................................*/
# if DSC_HCRMDE != 0
   rhp = &rph;
   ( ppt->rhp ) = rhp;
# if DSC_FLDMDE != 0
   fcp = &fcn;
   ( ppt->fcp ) = fcp;
# endif
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   ( state->tpt ) = tpt;
   ( state->ppt ) = ppt;
   ( state->bpt ) = bpt;
   ( state->ept ) = ept;
   ( state->vpt ) = vpt;
/*............................................................................*/
/* initialize struct S_MATRIX [ defined in smatrx(*) ] and pointer spt: */
   
   spt = smatrx( NULL );
/*............................................................................*/
# if DSC_HCRMDE != 0
/*............................................................................*/
# if BND_CSHAPE == 1
/* initialize struct CSHAPE [ defined in cshape(*) ] and pointer csh: */

   clp = cshape( clp ); 

   admce = sqrt( EPS_VAC / MY_VAC_ ); /* char field admittance */

# endif
/*............................................................................*/
/* initialize struct hsx of type HCRSMX [ defined in hcrsmx(*) ] and pointer */
/* hsp->hsx [ thermal and fluid s-matrix structure]: */

   hsp = &hsx;
   ( hsp->opt ) = null;
/*............................................................................*/
   hsp = hcrsmx( hsp );    /*                                     */
/*.......................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* initialize repeater structure: */

   hh = null; do
   {
      ( ppt->rep->smx[hh] ) = null;
/*............................................................................*/
# if DSC_HCRMDE != 0
      ( ppt->rhp->smx[hh] ) = null;
# if DSC_FLDMDE != 0
      ( ppt->fcp->cnn[hh] ) = null;
# endif
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   } while(( hh++ ) < NODES );
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE != 0
   kk = null; do
   {
      ( ppt->fcp->crsdt[kk] ) = ZERO;
   } while(( kk++ ) < NFCNN );
# endif
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* display: */

   strcpy(( csp->cmmnt ), "Welcome to DSC-FORMER !" );
   strcpy(( csp->cnfrm ), "Nothing done! Do you really want to quit ?" );

   ( csp->clscr ) = -ONE;
   ( csp->dfopt ) = 7; /* the initial default menu option */

  menu:

   ( csp->items ) = 7;
   ( csp->dflnf ) = ( csp->dfopt ); /* set line feed before default option */
                                    /* line */
/*
   if (( csp->clscr ) == null )
       ( csp->clscr ) = ONE;
*/
   strcpy(( csp->envmt ), "DSC-FORMER" );
   strcpy(( csp->tasks ), "Create DSC system file ...:" );

   strcpy(( csp->mline[1] ), "* topology file >------> [ dsc.top<n> ]" );
   strcpy(( csp->mline[2] ), "* S-matrix file >------> [ dsc.smx<n> ]" );
   strcpy(( csp->mline[3] ), "* boundary file >------> [ dsc.bnd<n> ]" );
   strcpy(( csp->mline[4] ), "* excitation file >----> [ dsc.exc<n> ]" );
   strcpy(( csp->mline[5] ), "* evaluation file >----> [ dsc.val<n> ]" );
   strcpy(( csp->mline[6] ), "* SOLVER batch file >--> [ dsc.batch ]" );
   strcpy(( csp->mline[7] ), "* Create all DSC system files" );

   strcpy(( csp->escpe ), "End of program / escape" );
/*............................................................................*/
   csp = txcnsl( csp );   /* build the [ start ] menu                         */
/*......................*/
   ( state->item ) = ( csp->option );

   ( state->stat ) = null;

   if (( state->item ) == FORM_ITEMS )
   {
      idx1 = ONE;
      idx2 = FORM_ITEMS-TWO;

      printf( "\n" );

      strcpy(( csp->rqlng ),
         "Please enter common file index <N>" );
      strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
      csp = txcnsl( csp );      /* enter common file index from text console  */
/*............................*/
      strcpy( ptr, ( lotos(( csp->inlng ), null )));
      strncpy(( state->flbl ), ptr, THREE );
   } 
   else
   {
      idx1 = ( state->item );
      idx2 = ( state->item );
   };
/*............................................................................*/
   for ( idx=idx1; idx<=idx2; idx++ )
   {
      switch( idx )
      {
        case 0:

         PRBLDCLR( "\r" );
         printf( "\r %*s", 78, "DSC-FORMER" );
         PRNORMAL( "\r " );

         return null;
/*............................................................................*/
        case 1: /* >------- topological structure file generation ----------> */

         strcpy( fleptr, prefix );
	 strcat( fleptr, sysptr );

	 if (( state->item ) < FORM_ITEMS )
         {
            printf( "\n" );

            strcpy(( csp->rqlng ),
               "Please enter index <N> of topology file '" );
            strcat(( csp->rqlng ), fleptr );
            strcat(( csp->rqlng ), "<N>'" );
            strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
            csp = txcnsl( csp );           /* enter index N on text console   */
/*.......................................*/
            strcpy( ptr, ( lotos(( csp->inlng ), null )));
            strncpy(( state->flbl ), ptr, THREE );
         };

         strcat( fleptr, ( state->flbl ));
/*............................................................................*/
# if FORM_DISP == 1
         if (( state->item ) == FORM_ITEMS )
         {
            printf( "\n DSC mesh generation function "
               "%s called.                    ", "systop(*)" );
         }
         else
         {
            printf( "\n DSC mesh topology function "
               "%s called.                  ", "systop(*)" );
         };
         PRBLDCLR( "\r" );
         printf( "\n %*s", 78, "DSC-FORMER" );
         PRNORMAL( "" );
# endif
/*............................................................................*/
         ( tpt->rtn ) = -TWO;
/*............................................................................*/
         tpt = systop( state );            /*  enter DSC mesh topology        */
/*.......................................*/
         ( state->tpt ) = tpt;

         switch( tpt->rtn )
         {
           case -TWO:
            printf( "\n Message from function %s:      ", __func__ );
            printf( "\n Unsatisfied code in function %s"
                    " !!! \n", "systop(*)" );
            ( csp->dfopt ) = null; /* the next default menu option */
            goto menu;

           case -ONE:
            ( csp->dfopt ) = null; /* the next default menu option */
            goto menu;

           case ONE:
            printf( "\n\n Message from function %s:", __func__ );
            printf( "\n Error in function %s !!!\n", "systop(*)" );
            ( csp->dfopt ) = null; /* the next default menu option */
            goto menu;

           case FORM_ITEMS:
/*............................................................................*/
# if FORM_DISP == 1
            printf( "\n ==================================="
               "===========================================" );
            PRBLDCLR( "\r" );
            printf( "\n %*s", 78, "DSC-FORMER" );
            PRNORMAL( "\r" );
            printf( "\n DSC mesh generation function "\
               "%s terminated.                  ", "systop(*)" );
# endif
/*............................................................................*/
            idx2 = FORM_ITEMS - TWO;
            ( state->stat ) = FORM_ITEMS;
            ( state->tpt ) = tpt;

            break;

           default:
            ( state->tpt ) = tpt;
/*............................................................................*/
# if FORM_DISP == 1
            printf( "\n ==================================="
               "===========================================" );
            PRBLDCLR( "\r" );
            printf( "\n %*s", 78, "DSC-FORMER" );
            PRNORMAL( "\r" );
            printf( "\n DSC mesh topology function "\
               "%s terminated.                ", "systop(*)" );
# endif
/*............................................................................*/
            break;
         };
/*............................................................................*/
# if FORM_DISP == 1
         printf( "\n\n DSC system identifier: "
            "%.22s", ( tpt->name ));
         printf( "\n Number of mesh cells: "
            "%ld ", (( tpt->mf ) - ( tpt->mi ) + ONE ));
/*............................................................................*/
# if NFCES != 0
         printf( "\n Number of cell faces: "
            "%ld ", (( tpt->ff ) - ( tpt->fi ) + ONE ));
# endif
/*............................................................................*/
         printf( "\n Comment: %.70s ", tpt->text );
# endif /* FORM_DISP == 1 */
/*............................................................................*/
# if VERTX_LOG == 1
/* create log file of vertex point indices */

         ( tpt->ci ) = CPNTS;
         ( tpt->cf ) = null;

         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            cc = null; do 
            {
               if (( tpt->cm[hh][cc] ) < ( tpt->ci ))
                  ( tpt->ci ) = ( tpt->cm[hh][cc] );

               if (( tpt->cf ) < ( tpt->cm[hh][cc] ))
                  ( tpt->cf ) = ( tpt->cm[hh][cc] );
            } while(( ++cc ) < CRNRS );

            hh++;
         };

         vtxlog = fopen( vtxfle, "w" );

         printf( "\n opened: Vertex points logfile %s", vtxfle );

         fprintf( vtxlog, "%s\n\n", "Log file: Mesh cell vertex points " );
         fprintf( vtxlog, "%s\n", ( tpt->name ));
         fprintf( vtxlog, "%s\n\n", ( tpt->text ));
         fprintf( vtxlog, "initial_cell________: %7ld\n", ( tpt->mi ));
         fprintf( vtxlog, "final___cell________: %7ld\n", ( tpt->mf ));
/*............................................................................*/
# if NFCES != 0
         fprintf( vtxlog, "initial_face________: %7ld\n", ( tpt->fi ));
         fprintf( vtxlog, "final___face________: %7ld\n", ( tpt->ff ));
# endif
/*............................................................................*/
         fprintf( vtxlog, "initial_vertex_point: %7ld\n", ( tpt->ci ));
         fprintf( vtxlog, "final___vertex_point: %7ld\n", ( tpt->cf ));

         kk = null; do /* 70 characters '-' */
         {
            fprintf( vtxlog, "----------" );
         } while(( ++kk ) < SEVEN );
         fprintf( vtxlog, "----------\n" );

         fprintf( vtxlog, "| cell |" );
         fprintf( vtxlog, " point1 |" );
         fprintf( vtxlog, " point2 |" );
         fprintf( vtxlog, " point3 |" );
         fprintf( vtxlog, " point4 |" );
         fprintf( vtxlog, " point5 |" );
         fprintf( vtxlog, " point6 |" );
         fprintf( vtxlog, " point7 |" );
         fprintf( vtxlog, " point8 |" );

         fprintf( vtxlog, "\n" );
         kk = null; do /* 60 characters '-' */
         {
            fprintf( vtxlog, "----------" );
         } while(( ++kk ) < SIX );
         fprintf( vtxlog, "---------- DANSE-%s\n", DSC_RELEASE );

         hh = ( tpt->mi );
         while( hh <= ( tpt->mf )) 
         {
            fprintf( vtxlog, "%7ld", hh );
            jj = null; do
            {
               fprintf( vtxlog, " %8ld", ( tpt->cm[hh][jj] ));
            } while(( ++jj ) < CRNRS );
            fprintf( vtxlog, "\n" );
            hh++;
         };

         fprintf( vtxlog, "DANSE-%s ----------", DSC_RELEASE );
         kk = null; do /* 60 characters '-' */
         {
            fprintf( vtxlog, "----------" );
         } while(( ++kk ) < SIX );
/*............................................................................*/
/* close vertex point log file with date & time of creation: */ 

         nseconds = time( timer );
/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); '*/

         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( vtxlog, "\nLog file %s", vtxfle );
         fprintf( vtxlog, " created:\n%s\n", tmestr );

         fclose( vtxlog );

         printf( CLEAR_LINE );
         printf( "\r Vertex point logfile %s ", vtxfle );
         printf( time_frm, tmestr );
# endif /* VERTX_LOG == 1 */
/*............................................................................*/

         printf("\n");
/*............................................................................*/
         tpt = linker( tpt );  /* DSC mesh linker [interconnect neighb.cells] */
/*...........................*/

         if (( tpt->rtn ) == ONE )
         {
            printf( "\n\n Message from function %s:", __func__ );
            printf( "\n Error in function linker(*) !!!\n" );
            ( csp->dfopt ) = null; /* the next default menu option */
            goto menu;
         };
/*............................................................................*/
# if FORM_DISP == 1
         if (( state->stat ) == FORM_ITEMS )
         {
            printf( "\n DSC system %s entered and"
                    " linked. ", "systop(*)" );
         }
         else if (( state->stat ) != FORM_ITEMS )
         {
            printf( "\n DSC mesh topology %s entered and"
                    " linked. ", "systop(*)" );
         };
# endif
/*............................................................................*/
# if UNLNK_LOG == 1
/* create unliked ports file:*/

         linklog = fopen( lnkfle, "w" );

         printf( "\n opened: Link logfile %s", lnkfle );

         fprintf( linklog, "%s\n\n","Unlinked ports log file" );
         fprintf( linklog, "%s\n", ( tpt->name ));
         fprintf( linklog, "%s\n\n", ( tpt->text ));
         fprintf( linklog, "initial_cell________: %7ld\n", ( tpt->mi ));
         fprintf( linklog, "final___cell________: %7ld\n", ( tpt->mf ));

# if NFCES != 0
         fprintf( linklog, "initial_face________: %7ld\n", ( tpt->fi ));
         fprintf( linklog, "final___face________: %7ld\n", ( tpt->ff ));
# endif

         fprintf( linklog, "initial_vertex_point: %7ld\n", ( tpt->ci ));
         fprintf( linklog, "final___vertex_point: %7ld\n\n", ( tpt->cf ));
         fprintf( linklog, "unlinked_ports______: " );

         mn = null;
         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {

# if DSC_ADJCLL == 0
            jj = null; do
            {
               if (( tpt->mn[hh][jj] ) == null )
                  mn += TWO;
            } while(( ++jj ) < FACES );
# elif DSC_ADJCLL == 1
            jj = null; do
            {
               if (( tpt->mn[hh][jj] ) == null )
                  mn++;
            } while(( ++jj ) < PORTS );
# endif
            hh++;
         };

         if ( null < mn )
         {
            fprintf( linklog, "%7ld\n\n", mn );
            fprintf( linklog, "unlinked_ports:\n\n" );
# if DSC_ADJCLL == 0
            fprintf( linklog, "cell\tface\tport1\tport2\n" );
# elif DSC_ADJCLL == 1
            fprintf( linklog, "cell\tface\tport\n" );
# endif
            hh = ( tpt->mi ); 
            while( hh <= ( tpt->mf ))
            {
# if DSC_ADJCLL == 0
               jj = null; do
               {
                  if (( tpt->mn[hh][jj] ) == null )
                  {
                     fprintf( linklog, "%-7ld\t%1d\t%-2d\t%-2d\n",
# if DSC_FCELBL == 0
                        hh, jj, ( prt1[jj]+ONE ), ( prt2[jj]+ONE ));
# else
                        hh, ( jj+ONE ), ( prt1[jj]+ONE ), ( prt2[jj]+ONE ));
# endif /* DSC_FCELBL != 0 */
                  };
               } while(( ++jj ) < FACES );
# elif DSC_ADJCLL == 1
               jj = null; do
               {
                  if (( tpt->mn[hh][jj] ) == null )
                     fprintf( linklog, "%-7ld\t%1d\t%-2d\n",
# if DSC_FCELBL == 0
                        hh, fce[jj], ( jj+ONE ));
# else
                        hh, ( fce[jj]+ONE ), ( jj+ONE ));
# endif /* DSC_FCELBL != 0 */
               } while(( ++jj ) < PORTS );
# endif
               hh++;
            };
         }
         else if ( mn == null )
            fprintf( linklog, "\t none\n" );

         mn = null;
/*............................................................................*/
/* close linklog file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( linklog, "\nLog file %s ", lnkfle );
         fprintf( linklog, "created:\n%s\n", tmestr );

         fclose( linklog );

         printf( CLEAR_LINE ); 
         printf( "\r Link logfile %s ", lnkfle );
         printf( time_frm, tmestr );
# endif /* UNLNK_LOG == 1 */
/*............................................................................*/
# if EMWLL_LOG == ONE
/* create electric/magnetic walls file:*/

         wallog = fopen( wllfle, "w" );

         printf( "\n opened: Wall logfile %s", wllfle );

         fprintf( wallog, "%s\n\n", "log file: "
            "electric and magnetic walls " );
         fprintf( wallog, "%s\n", ( tpt->name ));
         fprintf( wallog, "%s\n\n", ( tpt->text ));
         fprintf( wallog, "initial_cell________: %7ld\n", ( tpt->mi ));
         fprintf( wallog, "final___cell________: %7ld\n", ( tpt->mf ));

# if NFCES != 0
         fprintf( wallog, "initial_face________: %7ld\n", ( tpt->fi ));
         fprintf( wallog, "final___face________: %7ld\n\n", ( tpt->ff ));
# endif

         fprintf( wallog, "electric_walls______: " ); 

         ii = null;
         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            jj = null; do
            {
# if DSC_ADJCLL == 0
               if (( tpt->mn[hh][jj] ) == ( - ONE ))
                  ii++;
# elif DSC_ADJCLL == 1
               if ((( tpt->mn[hh][(int)prt1[jj]] ) == ( - ONE ))
                 ||(( tpt->mn[hh][(int)prt2[jj]] ) == ( - ONE )))
                  ii++;
# endif
            } while(( ++jj ) < FACES );
            hh++;
         };
         if ( ii == null )
            fprintf( wallog, "\t none\n" );
         else
            fprintf( wallog, "%7ld faces\n", ii );

         mn = ii;

         fprintf( wallog, "magnetic_walls______: " );

         ii = null;
         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            jj = null; do
            {
# if DSC_ADJCLL == 0
               if (( tpt->mn[hh][jj] ) == ( - TWO ))
                  ii++;
# elif DSC_ADJCLL == 1
               if ((( tpt->mn[hh][(int)prt1[jj]] ) == ( - TWO ))
                 ||(( tpt->mn[hh][(int)prt2[jj]] ) == ( - TWO )))
                  ii++;
# endif
            } while(( ++jj ) < FACES );
            hh++;
         };
         if ( ii == null )
            fprintf( wallog, "\t none\n" );
         else
            fprintf( wallog, "%7ld faces\n", ii );

         if ( null < mn  )
         { 
            fprintf( wallog, "\nELECTRIC WALLS:\n\n" );
            fprintf( wallog, "cell\tface\tport1\tport2\n" );

            hh = ( tpt->mi ); 
            while( hh <= ( tpt->mf ))
            {
               jj = null; do
               {
# if DSC_ADJCLL == 0
                  if (( tpt->mn[hh][jj] ) == ( - ONE ))
# elif DSC_ADJCLL == 1
                  if ((( tpt->mn[hh][(int)prt1[jj]] ) == ( - ONE ))
                    ||(( tpt->mn[hh][(int)prt2[jj]] ) == ( - ONE )))
# endif
                  {
# if DSC_FCELBL == 0
                     fprintf( wallog, "%-7ld\t%1d\t%-2d\t%-2d\n",
                        hh, jj, ( prt1[jj]+ONE ), ( prt2[jj]+ONE ));
# else
                     fprintf( wallog, "%-7ld\t%1d\t%-2d\t%-2d\n",
                        hh, ( jj+ONE ), ( prt1[jj]+ONE ), ( prt2[jj]+ONE ));
# endif /* DSC_FCELBL != 0 */
                  };
               } while(( ++jj ) < FACES );
               hh++;
            };
         };

         mn = ii;

         if ( null < mn )
         {
            fprintf( wallog, "\nMAGNETIC WALLS:\n\n" );
            fprintf( wallog, "cell\tface\tport1\tport2\n" );

            hh = ( tpt->mi );
            while( hh <= ( tpt->mf ))
            {
               jj = null; do
               {
# if DSC_ADJCLL == 0
                  if (( tpt->mn[hh][jj] ) == ( - TWO ))
# elif DSC_ADJCLL == 1
                  if ((( tpt->mn[hh][(int)prt1[jj]] ) == ( - TWO ))
                    ||(( tpt->mn[hh][(int)prt2[jj]] ) == ( - TWO )))
# endif
                  {
# if DSC_FCELBL == 0
                     fprintf( wallog, "%-7ld\t%1d\t%-2d\t%-2d\n",
                        hh, jj, ( prt1[jj]+ONE ), ( prt2[jj]+ONE ));
# else
                     fprintf( wallog, "%-7ld\t%1d\t%-2d\t%-2d\n",
                        hh, ( jj+ONE ), ( prt1[jj]+ONE ), ( prt2[jj]+ONE ));
# endif /* DSC_FCELBL != 0 */
                  };
               } while(( ++jj ) < FACES );
               hh++;
            };
         }; /* end if ( null < mn ) */
/*............................................................................*/
/* close wall log file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( wallog, "\nLog file %s ", wllfle );
         fprintf( wallog, "created:\n%s\n", tmestr );

         fclose( wallog );

         printf( CLEAR_LINE );
         printf( "\r Wall logfile %s ", wllfle );
         printf( time_frm, tmestr );
# endif /* EMWLL_LOG == 1 */
/*............................................................................*/
/* create topology file: */

         topology = fopen( fleptr, "w" );

         printf( "\n opened: Topological structure file %s", fleptr );
/*............................................................................*/
/* system identifier [model name] and comment: */

         fprintf( topology, spformat, ( tpt->name )); /* DSC system identifier*/
         fprintf( topology, spformat, ( tpt->text )); /* comment, text etc.   */

         strcpy( ptr, "DSC-MODEL_TOPOLOGY_FILE_" );
         strcat( ptr, fleptr );

         pp = LNELNGTH - strlen( ptr );
         qq = null; do
         {
            fprintf( topology, "%c", 95 );
         } while(( ++qq ) < pp );
         fprintf( topology, "%s\n", ptr );
/*............................................................................*/
         fprintf( topology, "%-27s   % -ld\n",
            "first_mesh_cell_index:", ( tpt->mi ));
         fprintf( topology, "%-27s   % -ld\n",
            "last_mesh_cell_index:", ( tpt->mf ));
/*............................................................................*/
# if NFCES != 0
         fprintf( topology, "%-27s   % -ld\n",
            "first_face_index:", ( tpt->fi ));
         fprintf( topology, "%s-27s  % -ld\n",
            "last_face_index:", ( tpt->ff ));
# endif
/*............................................................................*/
         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
/*............................................................................*/
# if LBL_CELLS == 1
            fprintf( topology, "cell%ld:\n", hh );
# endif
/*............................................................................*/
# if DSC_ADJCLL == 0

            kk = null; do 
            {
               if (( tpt->mn[hh][kk] ) == null )
                  fprintf( topology, "0000000000000\n" ); /* unlinked port ! */
               else if (( tpt->mn[hh][kk] ) < null )
                  fprintf( topology, "%ld\n", ( tpt->mn[hh][kk] ));
               else if ( null < ( tpt->mn[hh][kk] ))
               {
                  fprintf( topology, "%-8ld", ( tpt->mn[hh][kk] ));
                  fprintf( topology, " %-2d", ( tpt->pn[hh][(int)prt1[kk]] ));
                  fprintf( topology, " %-2d", ( tpt->pn[hh][(int)prt2[kk]] ));
                  fprintf( topology, " %-2d\n", ( tpt->fn[hh][kk] ));
               };
            } while(( ++kk ) < FACES );

# elif DSC_ADJCLL == 1

            kk = null; do 
            {
               if (( tpt->mn[hh][kk] ) == null )
                  fprintf( topology, "0000000000000\n" ); /* unlinked port ! */
               else if (( tpt->mn[hh][kk] ) < null )
                  fprintf( topology, "%ld\n", ( tpt->mn[hh][kk] ));
               else if ( null < ( tpt->mn[hh][kk] ))
               {
                  fprintf( topology, "%-8ld", ( tpt->mn[hh][kk] ));
                  fprintf( topology, " %-2d", ( tpt->pn[hh][kk] ));
                  fprintf( topology, " %-2d\n", ( tpt->fn[hh][(int)fce[kk]] ));
               };
            } while(( ++kk ) < PORTS );

# endif /* DSC_ADJCLL == 1 */
/*............................................................................*/
            hh++;
         }; /* end while( hh <= ( tpt->mf )) */
/*............................................................................*/
/* close topology file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( topology, "\nDSC model topology file %s ", fleptr );
         fprintf( topology, "created:\n%s\n", tmestr );

         fclose( topology );

         printf( CLEAR_LINE ); 
         printf( "\r Topology file %s ", fleptr );
         printf( time_frm, tmestr );

         if (( state->item ) == FORM_ITEMS )
            printf( "\n ----------------------------------"
               "--------------------------------------------");
         break; 
/*............................................................................*/
        case 2:/* >----------------- s-matrix file generation --------------> */

         if (( state->item ) == FORM_ITEMS )
            printf( "\n\n S-matrix file generation:" );
     
         strcpy( fleptr, prefix );
         strcat( fleptr, smxptr );

         if (( state->stat ) == FORM_ITEMS )
         {
            strcpy(( ppt->name ), ( tpt->name ));
            strcat( fleptr, ( state->flbl ));
            printf( "\n\n Comment: %.70s \n", ( ppt->text ));
         }
         else if (( state->stat ) != FORM_ITEMS )
         {
            if (( state->item ) < FORM_ITEMS )
            {
               printf( "\n" );

               strcpy(( csp->rqlng ),
                  "Please enter index <N> of S-matrix file '" );
               strcat(( csp->rqlng ), fleptr );
               strcat(( csp->rqlng ), "<N>'" );
               strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
               csp = txcnsl( csp );            /* request/enter index N       */
/*...........................................*//* from console                */
               strcpy( ptr, ( lotos(( csp->inlng ), null )));
               strncpy(( state->flbl ), ptr, THREE );
            }
            else
               printf( "\n" );

            strcat( fleptr, ( state->flbl ));
/*............................................................................*/
# if FORM_DISP == 1
            printf( "\n DSC model coordinates and "
               "media functions called.      " );
            PRBLDCLR( "\r" );
            printf( "\n %*s", 78, "DSC-FORMER" );
            PRNORMAL( "" );
# endif
/*............................................................................*/
            if (( state->item ) < FORM_ITEMS )
            { 
/*............................................................................*/
               tpt = systop( state );             /* enter DSC mesh topology  */
/*..............................................*//* [cf. DSC model function] */
               ( state->tpt ) = tpt;

               switch( tpt->rtn )
               { 
                 case -ONE:
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 case ONE:
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Error in function %s !!!\n", "systop(*)" ); 
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 case TWO: /* = ( state->item ) [ coordinates, media ] */
                  ( state->stat ) = ( state->item );
/*............................................................................*/
# if FORM_DISP == 1
                  printf( "\n ==================================="
                     "===========================================" );
                  PRBLDCLR( "\r" );
                  printf( "\n %*s", 78, "DSC-FORMER" );
                  PRNORMAL( "\r" );
                  printf( "\n DSC model coordinates and media "\
                     "function %s terminated.           ", "systop(*)" );
# endif
/*............................................................................*/
                  goto identf_check2;

                 default:
                  break;
               };
            };

            ( state->tpt ) = tpt;
            ( ppt->rtn ) = -TWO;
/*............................................................................*/
/* clear repeater structures: */
            hh = null; do
            {
               ( ppt->rep->smx[hh] ) = null;
/*............................................................................*/
# if DSC_HCRMDE != 0
               ( ppt->rhp->smx[hh] ) = null;
# if DSC_FLDMDE != 0
               ( ppt->fcp->cnn[hh] ) = null;
# endif
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
            } while(( hh++ ) < ( tpt->mf ));
/*............................................................................*/
/* call mesh cell coordinates & media parameters: */

/*............................................................................*/
            ppt = syssmx( state );             /* mesh coordinates, media etc.*/
/*...........................................*//* [ cf. DSC model function ]  */

            switch( ppt->rtn )
            {
              case -TWO:
               printf( "\n Message from function %s:      ", __func__ );
               printf( "\n Unsatisfied code in function %s !!!\n",
                  "syssmx(*)" );
               ( csp->dfopt ) = null; /* the next default menu option */
               goto menu;

              case -ONE: 
               ( csp->dfopt ) = null; /* the next default menu option */
               goto menu;

              case ONE:
               printf( "\n\n Message from function %s:", __func__ );
               printf( "\n Error in function %s !!!          \n",
                  "syssmx(*)" );
               ( csp->dfopt ) = null; /* the next default menu option */
               goto menu;

              default:
               printf( "\n ==================================="
                  "===========================================" );
               PRBLDCLR( "\r" );
               printf( "\n %*s", 78, "DSC-FORMER" );
               PRNORMAL( "\r" );
               printf( "\n DSC model coordinates and media "\
                  "function %s terminated.          ", "syssmx(*)" );
               break;
            };

           identf_check2:

            if ( null != strncmp(( tpt->name ),( ppt->name ), FOUR ) )
            {
               printf( "\n Message from function %s:", __func__ );
               printf( "\n Incompatible DSC system identifiers "
                  "on functions %s and %s.",
                     "systop(*)", "syssmx(*)" );
               printf( "\n Please check functions and restart program !\n" );
               return ONE;             
            }; 
/*............................................................................*/
# if FORM_DISP == 1
            if (( state->item ) < FORM_ITEMS )
            {
               printf( "\n\n DSC system identifier [ structure name ]: "
                  "%.22s", ( ppt->name ));
               printf( "\n Comment: %.70s ", ( ppt->text ));
            };
# endif
/*............................................................................*/
            ( state->tpt ) = tpt;
            ( state->ppt ) = ppt;
         }; /* end if (( state->stat ) != FORM_ITEMS ) */
/*............................................................................*/
# if COORD_LOG == 1
/* create log file of vertex point coordinates */ 

         ( tpt->ci ) = CPNTS;
         ( tpt->cf ) = null;

         hh = ( tpt->mi ); 
         while( hh <= ( tpt->mf )) 
         {
            cc = null; do 
            {
               if (( tpt->cm[hh][cc] ) < ( tpt->ci )) 
                  ( tpt->ci ) = ( tpt->cm[hh][cc] );

               if (( tpt->cf ) < ( tpt->cm[hh][cc] )) 
                  ( tpt->cf ) = ( tpt->cm[hh][cc] );
            } while(( ++cc ) < CRNRS );
            hh++;
         };

         coordlog = fopen( corfle, "w" );

         printf( "\n opened: Coordinates logfile %s", corfle );

         fprintf( coordlog, "%s\n\n",
            "Log file: Vertex point coordinates [ m ]" );
         fprintf( coordlog, "%s\n", ( ppt->name ));
         fprintf( coordlog, "%s\n\n", ( ppt->text ));
         fprintf( coordlog, "initial_cell________: %7ld\n", ( tpt->mi ));
         fprintf( coordlog, "final___cell________: %7ld\n", ( tpt->mf ));
/*............................................................................*/
# if NFCES != 0
         fprintf( coordlog, "initial_face________: %7ld\n", ( tpt->fi ));
         fprintf( coordlog, "final___face________: %7ld\n", ( tpt->ff ));
# endif
/*............................................................................*/
         fprintf( coordlog, "initial_vertex_point: %7ld\n", ( tpt->ci ));
         fprintf( coordlog, "final___vertex_point: %7ld\n", ( tpt->cf ));

         kk = null; do /* 60 characters '-' */
         {
            fprintf( coordlog, "----------" );
         } while(( ++kk ) < SIX );
         fprintf( coordlog, "---------\n" );

         fprintf( coordlog, "| point |" );
         fprintf( coordlog, "         x         |" );
         fprintf( coordlog, "         y         |" );
         fprintf( coordlog, "         z         |" );

         fprintf( coordlog, "\n" );
         kk = null; do /* 50 characters '-' */
         {
            fprintf( coordlog, "----------" );
         } while(( ++kk ) < FIVE );
         fprintf( coordlog, "--------- DANSE-%s\n", DSC_RELEASE );

         hh = ( tpt->ci );
         while( hh <= ( tpt->cf ))
         {
            fprintf( coordlog, "%8ld % .12e % .12e % .12e\n", hh,
               ( ppt->cpt->c[hh][null] ), 
                  ( ppt->cpt->c[hh][ONE] ), 
                     ( ppt->cpt->c[hh][TWO] ));
            hh++;
         };

         fprintf( coordlog, "DANSE-%s ---------", DSC_RELEASE );
         kk = null; do /* 50 characters '-' */
         {
            fprintf( coordlog, "----------" );
         } while(( ++kk ) < FIVE );
/*............................................................................*/
/* close coordinates log file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( coordlog, "\nCoordinates log file %s ", corfle );
         fprintf( coordlog, "created:\n%s\n", tmestr );

         fclose( coordlog );

         printf( CLEAR_LINE );
         printf( "\r Coordinates logfile %s ", corfle );
         printf( time_frm, tmestr );
# endif /* COORD_LOG == 1 */
/*............................................................................*/
# if CSIZE_LOG != 0
/* create cell size log file */

         csizelog = fopen( cszfle, "w" );

         printf( "\n opened: Cellsize logfile %s", cszfle );

         fprintf( csizelog, "%s\n\n","Log file: cell sizes "
            "[ maximum vertex point distances ]" );
         fprintf( csizelog, "%s\n", ( ppt->name ));
         fprintf( csizelog, "%s\n\n", ( ppt->text ));
         fprintf( csizelog, "initial_cell________: %7ld\n", ( tpt->mi ));
         fprintf( csizelog, "final___cell________: %7ld\n\n", ( tpt->mf ));
         fprintf( csizelog, "cell\tvertex pair\tdistance[m]\n" );

         jj = null;
         while( jj <= CSIZE_LOG )
         {
            mxsze[jj] = ZERO;
            jj++;
         };

         mmsze = 1.e+277;
         mn = null;
         ii = null;
         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            if (( ppt->mpt->idx[hh] ) == null )
            {
               fprintf( csizelog, "%-7ld -----trivial_cell----------\n", hh );
            } 
            else
            {
               csize = ZERO;
               jj = null; do
               {
                  mm = ( tpt->cm[hh][jj] );
                  xx1 = ( ppt->cpt->c[mm][null] );
                  yy1 = ( ppt->cpt->c[mm][ONE] );
                  zz1 = ( ppt->cpt->c[mm][TWO] );

                  kk = jj + ONE;
                  while( kk < EIGHT )
                  {
                     nn = ( tpt->cm[hh][kk] );
                     xx2 = xx1 - ( ppt->cpt->c[nn][null] );
                     yy2 = yy1 - ( ppt->cpt->c[nn][ONE] );
                     zz2 = zz1 - ( ppt->cpt->c[nn][TWO] );

                     qdist = xx2*xx2 + yy2*yy2 + zz2*zz2;

                     if ( csize < qdist )
                     {
                        csize = qdist;
                        pp = jj;
                        qq = kk;
                     };
                     kk++ ;
                  };
               } while(( ++jj ) < EIGHT );
               csize = sqrt( csize );

               jj = null;
               while( jj < CSIZE_LOG )
               {
                  if ( mxsze[jj] < csize )
                  {
                     mxsze[jj] = csize;
                     mxcll[jj] = hh;
                     vertx[jj][null] = pp;
                     vertx[jj][ONE] = qq;

                     ii++;
                     break;
                  };
                  jj++;
               };

               if ( csize < mmsze )
               {
                  mmsze = csize;
                  mn = hh;
               };
            
               fprintf( csizelog, "%-7ld\t%d<->%d\t\t%.5e"
                  "\n", hh, pp, qq, csize );
            };
            hh++;
         };

         fprintf( csizelog, "\nThe %d %s\n",
            CSIZE_LOG, "largest cells:" );

         if ( CSIZE_LOG < ii )
            ii = CSIZE_LOG;

         jj = null;
         while( jj < ii )
         {
            fprintf( csizelog, "%-7ld\t%d<->%d\t\t%.5e"
               "\n", mxcll[jj], vertx[jj][0], vertx[jj][1], mxsze[jj] );
            jj++;
         };
         fprintf( csizelog, "\nThe smallest cell:\n%-7ld\tmax.dist.:\t%.5e"
            "\n", mn, mmsze );
/*............................................................................*/
/* close cell size log file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( csizelog, "\nCell size log file %s ", cszfle );
         fprintf( csizelog, "created:\n%s\n", tmestr );

         fclose( csizelog );

         printf( CLEAR_LINE );
         printf( "\r Cellsize logfile %s ", cszfle );
         printf( time_frm, tmestr );
# endif /* CSIZE_LOG == 1 */
/*............................................................................*/
# if MEDIA_LOG == 1
/* create log file of mesh media [ material parameters ]*/

         medialog = fopen( medfle, "w" );

         printf( "\n opened: Media logfile %s", medfle );

         fprintf( medialog, "%s\n\n","Log file: mesh media" );
         fprintf( medialog, "%s\n", ( ppt->name ));
         fprintf( medialog, "%s\n\n", ( ppt->text ));
         fprintf( medialog, "initial_cell________: %7ld\n", ( tpt->mi ));
         fprintf( medialog, "final___cell________: %7ld\n", ( tpt->mf ));

         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            qq = ( ppt->mpt->idx[hh] );

            fprintf( medialog, "\ncell%ld [ medium %ld ]:", hh, qq );
/*............................................................................*/
/* gyroelectric parameters: */

            if (( ppt->mpt->no[qq] ) != ZERO )
            {
               fprintf( medialog, "\nelectron density [1/meter^3]............: "
                  "%.12e\n", ( ppt->mpt->no[qq] ));
               fprintf( medialog, "plasma relaxation time [seconds]........: "
                  "%.12e\n", ( ppt->mpt->tg[qq] ));
               fprintf( medialog, "internal magnetic field [Amperes/meter].:\n"
                  "Hx = %+.9e\tHy = %+.9e\tHz = %+.9e\n",
                     ( ppt->mpt->mi[qq][0] ),
                        ( ppt->mpt->mi[qq][1] ),
                           ( ppt->mpt->mi[qq][2] )); 
            };
/*............................................................................*/
/* gyromagnetic parameters: */

            if (( ppt->mpt->ld[qq] ) != ZERO ) 
            {
               fprintf( medialog, "\nLANDE' factor...........................: "
                  "%.12e\n", ( ppt->mpt->ld[qq] ));
               fprintf( medialog, "spin-spin relaxation time [seconds].....: "
                  "%.12e\n", ( ppt->mpt->tg[qq] ));
               fprintf( medialog, "saturation magnetization [Tesla]........:\n"
                  "Bx = %+.9e\tBy = %+.9e\tBz = %+.9e\n",
                     ( ppt->mpt->ms[qq][0] ),
                        ( ppt->mpt->ms[qq][1] ),
                           ( ppt->mpt->ms[qq][2] ));
               fprintf( medialog, "internal magnetic field [Amperes/meter].:\n"
                  "Hx = %+.9e\tHy = %+.9e\tHz = %+.9e\n",
                     ( ppt->mpt->ms[qq][0] ),
                        ( ppt->mpt->ms[qq][1] ),
                           ( ppt->mpt->ms[qq][2] ));
            };
/*............................................................................*/
# if DSC_HCRMDE != 0
/* thermal parameters: */

            if (( ppt->mpt->cv[qq] ) != ZERO )
            {
               fprintf( medialog, "\nheat capacity [J/(K*m^3)].:\n"
                  "%.12e\n", ( ppt->mpt->cv[qq] ));
            };

            if (( ppt->mpt->kh[qq] ) != ZERO )
            {
               fprintf( medialog, "\nheat conductivity [W/(K*m)].:\n"
                  "%.12e\n", ( ppt->mpt->kh[qq] ));
            };

# if DSC_FLDMDE != 0

            if (( ppt->mpt->rm[qq] ) != ZERO )
            {
               fprintf( medialog, "\nmean density [Kg/(m^3)].:\n"
                  "%.12e\n", ( ppt->mpt->rm[qq] ));
            };

            if (( ppt->mpt->tm[qq] ) != ZERO )
            {
               fprintf( medialog, "\nmean temperature [Kg/(sec^2*m)].:\n"
                  "%.12e\n", ( ppt->mpt->tm[qq] ));
            };

            if (( ppt->mpt->bm[qq] ) != ZERO )
            {
               fprintf( medialog, "\n expansion coefficient [1/K].:\n"
                  "%.12e\n", ( ppt->mpt->bm[qq] ));
            };

            if (( ppt->mpt->cm[qq] ) != ZERO )
            {
               fprintf( medialog, "\n adiabatic compression coeff [1/Pa].:\n"
                  "%.12e\n", ( ppt->mpt->cm[qq] ));
            };

            if (( ppt->mpt->ny[qq] ) != ZERO )
            {
               fprintf( medialog, "\ndynamic viscosity [Kg/(sec*m)].:\n"
                  "%.12e\n", ( ppt->mpt->ny[qq] ));
            };

            if (( ppt->mpt->q1[qq] ) != ZERO )
            {
               fprintf( medialog, "\nCp/Cv - 1 [dimensionless]:\n"
                  "%.12e\n", ( ppt->mpt->q1[qq] ));
            };

            if (( ppt->mpt->td[qq] ) != ZERO )
            {
               fprintf( medialog, "\ndissipation time constant [sec]:\n"
                  "%.12e\n", ( ppt->mpt->td[qq] ));
            };
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* electric parameters: */

            fprintf( medialog, "\nrelative permittivity:\n" );
            fprintf( medialog, "%.9e %.9e %.9e\n",
               ( ppt->mpt->ep[qq][0] ),
                  ( ppt->mpt->ep[qq][3] ),
                     ( ppt->mpt->ep[qq][5] ));
            fprintf( medialog, "\t*\t%.9e %.9e\n",
               ( ppt->mpt->ep[qq][1] ),
                  ( ppt->mpt->ep[qq][4] ));
            fprintf( medialog, "\t*\t\t*\t%.9e\n",
               ( ppt->mpt->ep[qq][2] ));
            fprintf( medialog, "\nelectric conductivity [ Siemens/meter ]:\n");

            fprintf( medialog, "%.9e %.9e %.9e\n",
               ( ppt->mpt->ke[qq][0] ),
                  ( ppt->mpt->ke[qq][3] ),
                     ( ppt->mpt->ke[qq][5] ));
            fprintf( medialog, "\t*\t%.9e %.9e\n",
               ( ppt->mpt->ke[qq][1] ),
                  ( ppt->mpt->ke[qq][4] ));
            fprintf( medialog, "\t*\t\t*\t%.9e\n",
               ( ppt->mpt->ke[qq][2] ));
            fprintf( medialog, "\nrelative permeability:\n" );
/*............................................................................*/
/* magnetic parameters: */

            fprintf( medialog, "%.9e %.9e %.9e\n",
               ( ppt->mpt->my[qq][0] ),
                  ( ppt->mpt->my[qq][3] ),
                     ( ppt->mpt->my[qq][5] ));
            fprintf( medialog, "\t*\t%.9e %.9e\n",
               ( ppt->mpt->my[qq][1] ),
                  ( ppt->mpt->my[qq][4] ));
            fprintf( medialog, "\t*\t\t*\t%.9e\n",
               ( ppt->mpt->my[qq][2] ));
            fprintf( medialog, "\nmagnetic conductivity [ Ohm/meter ]:\n" );

            fprintf( medialog, "%.9e %.9e %.9e\n",
               ( ppt->mpt->km[qq][0] ),
                  ( ppt->mpt->km[qq][3] ),
                     ( ppt->mpt->km[qq][5] ));
            fprintf( medialog, "\t*\t%.9e %.9e\n",
               ( ppt->mpt->km[qq][1] ),
                  ( ppt->mpt->km[qq][4] ));
            fprintf( medialog, "\t*\t\t*\t%.9e\n",
               ( ppt->mpt->km[qq][2] ));
            hh++;
         };
/*............................................................................*/
/* close media log file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( medialog, "\nMedia log file %s ", medfle );
         fprintf( medialog, "created:\n%s\n", tmestr );

         fclose( medialog );

         printf( CLEAR_LINE );
         printf( "\r Media logfile %s ", medfle );
         printf( time_frm, tmestr );
# endif /* MEDIA_LOG == 1 */
/*............................................................................*/
/* The structure of type S_MATRIX, pointed to by spt, is declared in smatrx(*)*/
/* and initialized or reset as follows: */ 
   
/*............................................................................*/
         spt = smatrx( NULL );         /*                                     */
/*...................................*/

/*............................................................................*/
# if DSC_HCRMDE != 0
/* Initialize or reset the structure of type HCRSMX pointed to by hsp: */

         ( hsp->opt ) = null;
/*............................................................................*/
         hsp = hcrsmx( hsp );        /*                                       */
/*............................................................................*/
# endif
/*............................................................................*/
/* Time step determination [ Maxwell field s-matrix ]: */
/* Confirm any default time step or phase shift != 0, if such are stored */ 

         if ( *( ppt->domain ) == 't' )
         { 
            PRBLDCLR( "" );
            printf( "\n TIME STEP [ Maxwell field ]:" );
            PRNORMAL( "\n" );

            if ( ZERO < ( state->dt ))
            {
               printf( "\n There is a stored time step "
                  "%.16e seconds\n", ( state->dt ));
               strcpy(( csp->rqdbl ), "Accept, modify, or start "
                  "stability check [ enter 0 ]" );
               strcpy(( csp->rqfrm ), "brackets" );
               ( csp->dfdbl ) = ( state->dt ); /* default time step */
/*............................................................................*/
               csp = txcnsl( csp );        /* request on text console */
/*.......................................*/
               if ( ZERO < ( csp->indbl ))
               {
                  ( state->dt ) = ( csp->indbl );
                  ( spt->dt ) = ( state->dt );

                  printf( "\r The new time step stored is "
                     "%.16e seconds\n", ( spt->dt ));

                  goto terminal1;
               };
            }; /* end if ( ZERO < ( state->dt )) */
         } /* end if ... [ frequency domain ] */
         else if ( *( ppt->domain ) == 'f' )
	 { 
            PRBLDCLR( "" );
            printf( "\n PHASE SHIFT [ Maxwell field ]:" );
            PRNORMAL( "\n" );

            if ( ZERO < ( state->dp ))
            {
               printf( "\n There is a stored phase shift "
                  "%.16e radians\n", ( state->dp ));
               strcpy(( csp->rqdbl ), "Accept, modify, or start "
                  "stability check [ enter 0 ]" );
               strcpy(( csp->rqfrm ), "brackets" );
               ( csp->dfdbl ) = ( state->dp ); /* default phase shift */
/*............................................................................*/
               csp = txcnsl( csp );        /* request on text console */
/*.......................................*/
               if ( ZERO < ( csp->indbl ))
               {
                  ( state->dp ) = ( csp->indbl );
                  ( state->dt ) = ( state->dp )/( 2.*PI*( ppt->fr ));
                  ( spt->dp ) = ( state->dp );
                  ( spt->dt ) = ( state->dt );

                  printf( "\r The new phase shift stored is "
                     "%.16e radians\n", ( spt->dp ));

                  goto terminal1;
               };
            }; /* end if ( ZERO < ( state->dp )) */
         }; /* end if ... [ frequency domain ] */
/*............................................................................*/
/* display: */

         printf( "\n " );

/* stability time step determination [ for Maxwell field s-matrix ]: */

         if ( *( ppt->domain ) == 'f' )
            strcpy( dsp->messge,
               "Stability phase shift determination started. " );
         else
            strcpy( dsp->messge,
	       "Stability time step determination started. " );

         ( dsp->option ) = 's'; /* 's'tart message */ 
         dsplay( dsp ); 

         strcpy(( dsp->messge ),
            "[Please wait a moment]" );

         ( dsp->option ) = 'm'; /* 'm'essage under running cursor */
         dsplay( dsp ); 

         ( dsp->option ) = 'c'; /* running 'c'ursor effect */
	 ( dsp->range ) = ( tpt->mf ) - ( tpt->mi ) + ONE;
/*............................................................................*/
         ecl = null;
         mcl = null;

         ( state->dt ) = GIANT_VAL;
         ( state->fr ) = ( ppt->fr );

         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            if ( null == ( ppt->rep->smx[hh] ))
            {  
               qq = ( ppt->mpt->idx[hh] ); /* media index */

               if ( qq == null )
               {
                  ( spt->med ) = null;
                  ( spt->opt ) = null;
/*............................................................................*/
                  spt = smatrx( spt );       /* reset parameters */
/*.........................................*/
               }
               else if (( null ) /* [ hence doesn't apply ] */\
                      &&( null == strncmp(( ppt->mpt->type[qq] ),
                  "trivial", FIVE )))
               {
                  ( spt->med ) = qq;
                  ( spt->opt ) = null;
/*............................................................................*/
                  spt = smatrx( spt );       /* reset parameters */
/*.........................................*/
               }
               else /* [ potential ] non trivial medium */
               {
/* compute stability time step: */

                  ( spt->med ) = qq;

                  if ( *( ppt->domain ) == 'f' )
                  { 
                     ( spt->fr ) = ( state->fr );
                     ( spt->dp ) = ( state->dp );
                  }
                  else
                     ( spt->fr ) = ZERO;

                  ( spt->no ) = ( ppt->mpt->no[qq] );
                  ( spt->tr ) = ( ppt->mpt->tr[qq] );
                  ( spt->ld ) = ( ppt->mpt->ld[qq] );
                  ( spt->tg ) = ( ppt->mpt->tg[qq] );

                  kk = null; do
                  {
                     ( spt->mi[kk] ) = ( ppt->mpt->mi[qq][kk] );
                     ( spt->ms[kk] ) = ( ppt->mpt->ms[qq][kk] );
                     ( spt->hg[kk] ) = ( ppt->mpt->hg[qq][kk] );

	             ll = kk;
                     while( ll < DIMNS )
	             {
                        switch( ll - kk )
                        {
                          case 0:
                           ( spt->epr[kk][ll] ) = ( ppt->mpt->ep[qq][kk] );
                           ( spt->myr[kk][ll] ) = ( ppt->mpt->my[qq][kk] );
                           ( spt->ke[kk][ll] ) = ( ppt->mpt->ke[qq][kk] );
                           ( spt->km[kk][ll] ) = ( ppt->mpt->km[qq][kk] );
                           break;

                          case 1:
                           ( spt->epr[kk][ll] ) = \
                              ( ppt->mpt->ep[qq][kk+THREE] );
                           ( spt->myr[kk][ll] ) = \
                              ( ppt->mpt->my[qq][kk+THREE] );
                           ( spt->ke[kk][ll] ) = \
                              ( ppt->mpt->ke[qq][kk+THREE] );
                           ( spt->km[kk][ll] ) = \
                              ( ppt->mpt->km[qq][kk+THREE] );
                           break;
                       
                          case 2:
                           ( spt->epr[kk][ll] ) = ( ppt->mpt->ep[qq][kk+FIVE] );
                           ( spt->myr[kk][ll] ) = ( ppt->mpt->my[qq][kk+FIVE] );
                           ( spt->ke[kk][ll] ) = ( ppt->mpt->ke[qq][kk+FIVE] );
                           ( spt->km[kk][ll] ) = ( ppt->mpt->km[qq][kk+FIVE] );
                           break;
                        };
                        ( spt->epr[ll][kk] ) = ( spt->epr[kk][ll] );
                        ( spt->myr[ll][kk] ) = ( spt->myr[kk][ll] );
                        ( spt->ke[ll][kk] ) = ( spt->ke[kk][ll] );
                        ( spt->km[ll][kk] ) = ( spt->km[kk][ll] );

                        ll++;
                     };
                  } while(( ++kk ) < DIMNS );
/*............................................................................*/
# if CELL_POINTS
                  if ( hh == CELL_POINTS )
                     printf( "\n\n Corner points of cell no. %ld\n", hh );
# endif
/*............................................................................*/
                  cc = null; do
                  {
                     kk = null; do
                     {
                        ( spt->c[cc][kk] ) = \
                           ( ppt->cpt->c[( tpt->cm[hh][cc] )][kk] );
/*............................................................................*/
# if CELL_POINTS
                        if ( hh == CELL_POINTS )
                        {
                           printf( "\n vertex%1ld:    cpt->c[%06ld][%1d]  =  "
                              "%+.14e", cc, ( tpt->cm[hh][cc] ), kk,
                              ( ppt->cpt->c[( tpt->cm[hh][cc] )][kk] ));
                        };
# endif
/*............................................................................*/
                     } while(( ++kk ) < DIMNS );
                  } while(( ++cc ) < CRNRS );
/*............................................................................*/
# if CELL_POINTS
                  if ( hh == CELL_POINTS )
                  {
                     printf( "\n\n Please acknowledge "
                        "( enter any character ) ...........: " );
                     scanf( "%s", ptr );
                     printf( "\n " );
                  };
# endif
/*............................................................................*/
/* compute stability time step: */

                  ( spt->opt ) = 't';
/*............................................................................*/
                  spt = smatrx( spt );       /*                               */
/*.........................................*/
                  if (( spt->rtn ) == null )
                  {
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in 'smatrx(*)', "
                        "option 't', cell number %ld !!!\n", hh );

                     ( csp->dfopt ) = null; /* the next default menu option */

                     goto menu;
                  };

                  if (( *( spt->etyp ) != 't' )||( *( spt->mtyp ) != 't' ))
                  {
                     if (( spt->dt ) < ( state->dt ))
                     {
	                ( state->dt ) = ( spt->dt );
                        ( state->dp ) = 2.*PI*( state->dt )*( ppt->fr );
	                node = hh;
                     };

                     if ( *( spt->getp ) == 'g' )
                        ecl++;

                     if ( *( spt->gmtp ) == 'g' )
                        mcl++;
                  };
               }; /* end if media index qq != null */
            }; /* end if (( ppt->rep->smx[hh] ) == null ) */
	    ( dsp->state ) = ( ++hh );
            dsplay( dsp );
         }; /* next mesh cell index hh */
/*............................................................................*/
/* clear display: */

         printf( CLEAR_LINE );
/*............................................................................*/
         if ( GENDS < ecl )
         {
            printf( "\n\n Message from function %s:", __func__ );
            printf( "\n\n Too many electric current cells "
               "in DSC mesh %s !!!", ( tpt->name ));
            printf( "\n [ Maximum number is %d = macro GENDS "
               "in %s.", GENDS, __func__ );
            printf( "\n - Change macro only in compliance "
               "with memory resources " );
            printf( "\n   and same macro in function scattr(*) "
               "of program 'DSC.C'.] " );
            return ONE;
         };

         if ( GMNDS < mcl )
         {
            printf( "\n\n Message from function %s:", __func__ );
            printf( "\n\n Too many magnetic current cells "
               "in DSC mesh %s !!!", ( tpt->name ));
            printf( "\n [ Maximum number is %d = macro GMNDS "
               "in %s.", GMNDS, __func__ );
            printf( "\n - Change macro only in compliance "
               "with memory resources" );
            printf( "\n   and same macro in function scattr(*) "
               "of program SOLVER.C.]" );
            return ONE;
         }; 
/*............................................................................*/
/* confirm or chage time step / phase shift: */

         ( spt->dt ) = ( state->dt );

         if ( *( ppt->domain ) == 'f' )
         {
            ( state->dp ) = 2.*PI*( state->dt )*( ppt->fr );
            ( spt->dp ) = ( state->dp );

            printf( "\r DSC phase shift determination terminated:"\
               "                                     " );
            printf( "\n The stability phase shift is %.13e radians ",
               ( spt->dp ));
            printf( "\n [ minimum over all cells found "
               "at cell no.%ld ].\n", node );

            strcpy(( csp->rqdbl ),
               "Enter smaller phase shift or return [ escape: 0 ]" );
            strcpy(( csp->rqfrm ), "brackets" );
            ( csp->dfdbl ) = ( state->dp ); /* default phase shift */
/*............................................................................*/
            csp = txcnsl( csp );           /* request ph.shft on text console */
/*.......................................*/
            if (( csp->indbl ) <= ZERO )
            {  
               ( csp->dfopt ) = null; /* the next default menu option */
               goto menu;
            }
            else
               ( spt->dp ) = ( csp->indbl );

            while(( state->dp ) < ( spt->dp )) /* enter new phase shift */
            {  
               printf( "\r Non-stable phase shift - try again:\n" );

               strcpy( ptr, "Enter new phase shift [ 0 < dp <= " ); 
               strcat( ptr, dotos(( state->dp ), 8, "e" ));
               strcat( ptr, " ]" );
               strcpy(( csp->rqdbl ), ptr );
               strcpy(( csp->rqfrm ), "brackets" );
               ( csp->dfdbl ) = ( state->dp ); /* default phase shift */
/*............................................................................*/
               csp = txcnsl( csp );        /* request phase shift on console  */
/*.......................................*/
               if (( csp->indbl ) <= ZERO )
               {
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;
               }
               else
                  ( spt->dp ) = ( csp->indbl );
            };

            if (( spt->dp ) < ( state->dp ))
               printf( "\r The new phase shift is "\
                  "%.16e radians\n", ( spt->dp ));

            ( state->dp ) = ( spt->dp );
            ( state->dt ) = ( state->dp )/( 2.*PI*( ppt->fr ));
         }
         else /* if ( *( ppt->domain ) == 't' ) */
         {
            printf( "\r Maxwell field time step determination terminated:"\
               "                                       " );
            printf( "\n The stability time step is %.13e seconds",
               ( spt->dt ));
	    printf( "\n [ minimum over all cells found "
               "on cell no.%ld ].\n", node );

            strcpy(( csp->rqdbl ),
               "Enter smaller time step or return [ escape: 0 ]" );
            strcpy(( csp->rqfrm ), "brackets" );
            ( csp->dfdbl ) = ( state->dt ); /* default time step */
/*............................................................................*/
            csp = txcnsl( csp );        /* request time step on text console  */
/*....................................*/
            if (( csp->indbl ) <= ZERO )
            {
               ( csp->dfopt ) = null; /* the next default menu option */
               goto menu;
            }
            else
               ( spt->dt ) = ( csp->indbl );
            
            while(( state->dt ) < ( spt->dt )) /* enter new time step */
            {  
               printf( "\r Non-stable time step - try again:\n" );

               strcpy( ptr, "Enter new time step [ 0 < dt <= " ); 
               strcat( ptr, dotos(( state->dt ), 8, "e" ));
               strcat( ptr, " ]" );
               strcpy(( csp->rqdbl ), ptr );
               strcpy(( csp->rqfrm ), "brackets" );
               ( csp->dfdbl ) = ( state->dt ); /* default time step */
/*............................................................................*/
               csp = txcnsl( csp );            /* request t.step on console   */
/*...........................................*/
               if (( csp->indbl ) <= ZERO )
               {
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;
               }
               else
                  ( spt->dt ) = ( csp->indbl );
            };

            if (( spt->dt ) < ( state->dt ))
               printf( "\r The new time step is "\
                  "%.16e seconds\n", ( spt->dt ));

            ( state->dt ) = ( spt->dt );
            ( state->dp ) = ZERO;
         }; /* end if *( ppt->domain ) == ... */

        terminal1: ;

/* end of time step / phase shift determination [ Maxwell field part ] */
/*............................................................................*/
# if DSC_HCRMDE != 0 
/* stability thermal&fluid time step determination: */
/*............................................................................*/
/* Confirm any default time step ( state->hcdt ) != 0, if such is stored */ 

         PRBLDCLR( "" );
/*............................................................................*/
# if DSC_FLDMDE == 0 
         printf( "\n TIME STEP [ thermal ]:" );
# else
         printf( "\n TIME STEP [ thermal&fluid ]:" );
# endif
/*............................................................................*/
         PRNORMAL( "\n" );

         pp = null;

         if ( ZERO < ( state->hcdt ))
         {
            printf( "\n There is a stored time step "
               "%.16e seconds\n", ( state->hcdt ));
            strcpy(( csp->rqdbl ), "Accept, modify, or start "
               "stability check [ enter 0 ]" );
            strcpy(( csp->rqfrm ), "brackets" );
            ( csp->dfdbl ) = ( state->hcdt ); /* default time step */
/*............................................................................*/
            csp = txcnsl( csp );        /* request on text console */
/*....................................*/
            if ( ZERO < ( csp->indbl ))
            {
               ( state->hcdt ) = ( csp->indbl );

               printf( "\r The new time step stored is "
                  "%.16e seconds\n", ( state->hcdt ));
               pp = ONE;
            }
            else
               pp = null;
         }; /* if ZERO < ( state->hcdt ) */
/*............................................................................*/
/* display: */

         if ( pp == null )
         {
            printf( "\n " );

            strcpy(( dsp->messge ),
               "Thermal stability time step determination started. " );

            ( dsp->option ) = 's'; /* 's'tart message */ 
            dsplay( dsp ); 

            strcpy(( dsp->messge ),
               "[Please wait a moment]" );

            ( dsp->option ) = 'm'; /* 'm'essage under running cursor */
            dsplay( dsp ); 

            ( dsp->option ) = 'c'; /* running 'c'ursor effect */
            ( dsp->range ) = ( tpt->mf ) - ( tpt->mi ) + ONE;

            ( state->hcdt ) = GIANT_VAL;
            ( hsp->dt ) = ( state->hcdt );
         }; /* end if ( pp == null ) */
/*............................................................................*/
# if DSC_FLDMDE != 0
         ( ppt->fcp->cnn[null] ) = null;
# endif
/*............................................................................*/
         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            if ( null == ( ppt->rhp->smx[hh] ))
            {  
               qq = ( ppt->mpt->idx[hh] ); /* media index */

               if ( qq == null )
               {
                  ( hsp->opt ) = null;
/*............................................................................*/
                  hsp = hcrsmx( hsp );            /* reset structure hsp */
/*..............................................*/
               }
               else /* non trivial medium ( qq != null ) */
               {
/* find stability time step [minimum over all mesh cells]: */
/* media index */
                  ( hsp->med ) = qq;

/* thermal properties: */
                  ( hsp->cv ) = ( ppt->mpt->cv[qq] );
                  ( hsp->kh ) = ( ppt->mpt->kh[qq] );

/* electromagnetic properties [lossy medium]: */
                  ( hsp->ke ) = ( ppt->mpt->ke[qq][null] );
                  ( hsp->km ) = ( ppt->mpt->km[qq][null] );
/*............................................................................*/
# if DSC_FLDMDE != 0
/* fluid properties: */
                  ( hsp->rm ) = ( ppt->mpt->rm[qq] );
                  ( hsp->tm ) = ( ppt->mpt->tm[qq] );
                  ( hsp->bm ) = ( ppt->mpt->bm[qq] );
                  ( hsp->cm ) = ( ppt->mpt->cm[qq] );
                  ( hsp->ny ) = ( ppt->mpt->ny[qq] );
                  ( hsp->q1 ) = ( ppt->mpt->q1[qq] );
                  ( hsp->td ) = ( ppt->mpt->td[qq] );
/*............................................................................*/
/* Prandtl turbulence models */
# if (( TURBMOD == 1 ) \
    ||( TURBMOD == 2 ))
                  ( hsp->LL ) = ( ppt->mpt->LL[qq] );
# endif /* TURBMOD != null */
/*............................................................................*/
                  kk = null; do
                  {
                     ( hsp->gr[kk] ) = ( ppt->mpt->gr[qq][kk] );
                     ( hsp->gp[kk] ) = ( ppt->mpt->gp[qq][kk] );
                  } while(( ++kk ) < THREE );

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# if CELL_POINTS
                  if ( hh == CELL_POINTS )
                     printf( "\n\n Corner points of cell no. %ld\n", hh );
# endif
/*............................................................................*/
/* vertex point coordinates: */

                  jj = null; do
                  {
                     kk = null; do
                     {
                        ( hsp->c[jj][kk] ) = \
                           ( ppt->cpt->c[( tpt->cm[hh][jj] )][kk] );
/*............................................................................*/
# if CELL_POINTS
                        if ( hh == CELL_POINTS )
                        {
                           printf( "\n vertex%1ld:    cpt->c[%06ld][%1d]  ="
                              "  %+.14e", jj, ( tpt->cm[hh][jj] ), kk,
                              ( ppt->cpt->c[( tpt->cm[hh][jj] )][kk] ));
                        };
# endif
/*............................................................................*/
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < CRNRS );
/*............................................................................*/
# if CELL_POINTS
                  if ( hh == CELL_POINTS )
                  {
                     printf( "\n\n Please acknowledge "
                        "( enter any character ) ...........: " );
                     scanf( "%s", ptr );
                     printf( "\n " );
                  };
# endif
                  ( hsp->opt ) = 't';
/*............................................................................*/
                  hsp = hcrsmx( hsp );            /* stability time step */
/*..............................................*/
                  if (( hsp->rtn ) == null )
                  {
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in function hcrsmx(*), "
                        "option 't', cell number %ld !!!\n", hh );

                     ( csp->dfopt ) = null; /*the next default menu option */
                     goto menu;
                  };

                  if ( *( hsp->ttyp ) != 't' )
                  {
                     if ( pp == null )
                     {
                        if (( hsp->dt ) < ( state->hcdt ))
                        {
                           ( state->hcdt ) = ( hsp->dt );
                           node = hh;
                        };
                     }; /* end if pp == null */
/*............................................................................*/
# if DSC_FLDMDE != 0
                     if ( *( hsp->ttyp ) == 'f' )
                     {
                        if (( ppt->fcp->cnn[hh] ) <= null )
                        {
                           fprintf( stderr, "\n\n Error message from function "
                              "%s:",  __func__ );
                           fprintf( stderr, "\n Illegal fluid connected "
                              "component index cnn = %d <= null",
                               ( ppt->fcp->cnn[hh] ));
                           fprintf( stderr, "\n defined in cell %ld, "
                              "medium %d !!!", hh, qq );
                           fprintf( stderr, "\n [ define: 0 < cnn <= %d "
                              "( = macro NFCNN in file FORMER.CONF )]\n\n ",
                              NFCNN );
                           return ONE;
                        }
                        else if (( ppt->fcp->cnn[null] )<( ppt->fcp->cnn[hh] ))
                           ( ppt->fcp->cnn[null] ) = ( ppt->fcp->cnn[hh] );
                     }; /* end if ( *( hsp->ttyp ) == 'f' ) */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
                  }; /* end if ( *( hsp->ttyp ) != 't' ) */
               }; /* end if media index ( qq != null ); non trivial medium */
            }; /* end if (( ppt->rhp->smx[hh] ) == null ) */
            hh++;

            if ( pp == null )
            {
               ( dsp->state ) = hh;
               dsplay( dsp );
            };
         }; /* next mesh cell index hh */
/*............................................................................*/
/* clear display: */

         printf( CLEAR_LINE );
/*............................................................................*/
         if ( pp == null )
         {
            ( hsp->dt ) = ( state->hcdt );
         
            printf( "\r Thermal time step determination terminated:"\
               "                                   " );
            printf( "\n The thermal stability time step is %.13e seconds ",
               ( hsp->dt ));
            printf( "\n [ minimum over all cells found "
               "on cell no.%ld ].\n", node );

            strcpy(( csp->rqdbl ),
               "Enter smaller time step or return [ escape: 0 ]" );
            strcpy(( csp->rqfrm ), "brackets" );
            ( csp->dfdbl ) = ( state->hcdt ); /* default time step */
/*............................................................................*/
            csp = txcnsl( csp );        /* request time step on text console  */
/*....................................*/
            if (( csp->indbl ) <= ZERO )
            {
               ( csp->dfopt ) = null; /* the next default menu option */
               goto menu;
            }
            else
               ( hsp->dt ) = ( csp->indbl );
            
            while(( state->hcdt ) < ( hsp->dt )) /* enter new time step */
            {  
               printf( "\r Non-stable thermal time step - try again:\n" );

               strcpy( ptr, "Enter new time step [ 0 < dt <= " ); 
               strcat( ptr, dotos(( state->hcdt ), 8, "e" ));
               strcat( ptr, " ]" );
               strcpy(( csp->rqdbl ), ptr );
               strcpy(( csp->rqfrm ), "brackets" );
               ( csp->dfdbl ) = ( state->hcdt ); /* default time step */
/*............................................................................*/
               csp = txcnsl( csp );            /* request t.step on console   */
/*...........................................*/
               if (( csp->indbl ) <= ZERO )
               {
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;
               }
               else
                  ( hsp->dt ) = ( csp->indbl );
            }; /* end while */

            if (( hsp->dt ) < ( state->hcdt ))
               printf( "\r The new thermal time step is "\
                  "%.16e seconds\n", ( hsp->dt ));
            ( state->hcdt ) = ( hsp->dt );
         }; /* end if ( pp == null ) */
# endif /* DSC_HCRMDE != 0 */ /* end of thermal time step determination */
/*............................................................................*/
/* CREATE S-MATRIX FILE: */

         ecl = null;
         mcl = null;
         
         if ( null == strncmp( ppt->domain, "time_domain", TWO ))
         {
            ( spt->fr ) = ( state->fr );
            ( spt->dt ) = ( state->dt );
         }
         else if ( null == strncmp( ppt->domain, "frequency_domain", TWO ))
         {
            ( spt->fr ) = ( state->fr );
            ( spt->dp ) = ( state->dp );
         };

         s_matrix = fopen( fleptr, "w" );

         printf( "\n opened: S-matrix file %s", fleptr );
/*............................................................................*/
/* system identifier [model name] and comment: */

         fprintf( s_matrix, spformat, ( ppt->name )); /* DSC system identifier*/
         fprintf( s_matrix, spformat, ( ppt->text )); /* comment, text etc.   */

         strcpy( ptr, "DSC-MODEL_S-PARAMETER_FILE_" );
         strcat( ptr, fleptr );

         pp = LNELNGTH - strlen( ptr );
         qq = null; do
         {
            fprintf( s_matrix, "%c", 95 );
         } while(( ++qq ) < pp );
         fprintf( s_matrix, "%s\n", ptr );
/*............................................................................*/
         fprintf( s_matrix, "%-27s   % -ld\n",
            "first_mesh_cell_index:", ( tpt->mi ));
         fprintf( s_matrix, "%-27s   % -ld\n",
            "last_mesh_cell_index:", ( tpt->mf ));
/*............................................................................*/
# if NFCES != 0
         fprintf( s_matrix, "%-27s   % -ld\n",
            "first_face_index:", ( tpt->fi ));
         fprintf( s_matrix, "%-27s   % -ld\n",
            "last_face_index:", ( tpt->ff ));
# endif
/*............................................................................*/
# if SMTRX_LOG == 1
/* create S-matrix logfile: */

         smtrxlog = fopen( smxfle, "w" );

         printf( "\n opened: S-matrix logfile %s", smxfle );

         fprintf( smtrxlog, "\n%s\n\n","Log file: s-parameters " );
         fprintf( smtrxlog, "%s\n", ( tpt->name ));
         fprintf( smtrxlog, "%s\n\n", ( tpt->text ));

         if ( null == strncmp(( ppt->domain ), "time_domain", TWO ))
         {
            fprintf( smtrxlog, spformat, "TIME_DOMAIN___________" );
            fprintf( smtrxlog, "DSC_time_step________: %.12e seconds ",
               ( spt->dt ));
         }
         else if ( null == strncmp(( ppt->domain ), "frequency_domain", TWO ))
         {
            fprintf( smtrxlog, spformat, "FREQUENCY_DOMAIN______" );
            fprintf( smtrxlog, "frequency____________: %.12e Hertz\n",
               ( spt->fr ));
            fprintf( smtrxlog, "DSC_phase_shift______: %.12e radians ",
               ( spt->dp ));
         };
# endif
/*............................................................................*/
         if ( null == strncmp( ppt->domain, "time_domain", TWO ))
         {
            strcpy( dsp->messge, 
               "Real nodal s-parameter determination started. " );

            fprintf( s_matrix, spformat, "TIME_DOMAIN" );
            fprintf( s_matrix, "%-27s   %+.16e\n",
               "Mxw.field_tstep_[seconds]:", ( state->dt ));
            fprintf( s_matrix, "%-27s   %+.16e\n",
               "mean_frequency_[Hz]:", ( state->fr ));
         }
         else if ( null == strncmp(( ppt->domain ), "frequency_domain", TWO ))
         {
            strcpy( dsp->messge, 
               "Complex nodal s-parameter determination started. " );

            fprintf( s_matrix, spformat, "FREQUENCY_DOMAIN" );
            fprintf( s_matrix, "%-27s   %+.16e\n",
               "frequency_[Hz]:", ( state->fr ));
            fprintf( s_matrix, "%-27s   %+.16e\n",
               "phase_shift_[radians]:", ( state->dp ));
         };
/*............................................................................*/
# if DSC_HCRMDE != 0
         fprintf( s_matrix, "%-27s   %+.16e\n",
            "therm_time_step_[seconds]:", ( state->hcdt ));

/* start delay for thermal computations */
         fprintf( s_matrix, "%-27s   %+.16e\n",
            "therm_start_at_[seconds]:", ( state->hcrstart ));
/*............................................................................*/
# if DSC_FLDMDE != 0
/* delay with that fluid dynamical computations start */
/* make sure that ( hcrstart <= fldstart ): */
         if (( state->fldstart ) < ( state->hcrstart ))
            ( state->fldstart ) = ( state->hcrstart );

         fprintf( s_matrix, "%-27s   %+.16e\n",
            "fluid_start_at_[seconds]:", ( state->fldstart ));

         fprintf( s_matrix, "%-27s   % -d\n",
            "connected_fluid_components:", ( ppt->fcp->cnn[null] ));

         if ( null < ( ppt->fcp->cnn[null] ))
         {
            fprintf( s_matrix, "%-27s\n", "coarsening_periods" );

            kk = null;
            while(( kk++ ) < ( ppt->fcp->cnn[null] ))
            {
               strcpy( ptr, "component_no" );
               strcat( ptr, ( lotos( kk, null )));
               strcat( ptr, "_[seconds]" );
               fprintf( s_matrix, "%-27s   %+.16e\n", ptr,
                  ( ppt->fcp->crsdt[kk] ));
            };
         };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* compute & store s-parameters: */
/* initialize running cursor: */

         ( dsp->option ) = 's'; /* 's'tart message */
         dsplay( dsp );

         strcpy(( dsp->messge ), 
            "[Please wait a moment]" );

         ( dsp->option ) = 'm'; /* 'm'essage under running cursor */
         dsplay( dsp );

         ( dsp->option ) = 'c'; /* running 'c'ursor effect */
         ( dsp->range ) = ( tpt->mf ) - ( tpt->mi ) + ONE;
/*............................................................................*/
/* compute and store E/H field s-parameters: */

         fprintf( s_matrix, "\n%s\n\n", smxpar ); /* "S-PARAMETERS", e.g. */

         nn = null;
         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            if ( null < ( ppt->rep->smx[hh] ))
            {
               ( ppt->rep->smx[hh] ) =
                  ( ppt->rep->smx[( ppt->rep->smx[hh] )] );
            }
            else
            {
               qq = ( ppt->mpt->idx[hh] );

               if ( qq == null )
               {
                  ( spt->opt ) = null;
/*............................................................................*/
                  spt = smatrx( spt );        /* reset s-parameters */
/*..........................................*/
               }
               else /* non trivial media */
               {
                  ( spt->med ) = qq;

                  if ( *( ppt->domain ) == 'f' )
                  {
                     ( spt->fr ) = ( state->fr );
                     ( spt->dp ) = ( state->dp );
                  }
                  else
                  {
                     ( spt->fr ) = ZERO;
                     ( spt->dp ) = ZERO;
                  };

                  ( spt->dt ) = ( state->dt );

                  ( spt->no ) = ( ppt->mpt->no[qq] );
                  ( spt->tr ) = ( ppt->mpt->tr[qq] );
                  ( spt->ld ) = ( ppt->mpt->ld[qq] );
                  ( spt->tg ) = ( ppt->mpt->tg[qq] );

	          kk = null; do
                  {                          
                     ( spt->mi[kk] ) = ( ppt->mpt->mi[qq][kk] );
                     ( spt->ms[kk] ) = ( ppt->mpt->ms[qq][kk] );
                     ( spt->hg[kk] ) = ( ppt->mpt->hg[qq][kk] );

                     ll = kk;
                     while( ll < DIMNS )
                     {
                        switch( ll - kk )
                        {
                          case 0:
                           ( spt->epr[kk][ll] ) = ( ppt->mpt->ep[qq][kk] );
                           ( spt->myr[kk][ll] ) = ( ppt->mpt->my[qq][kk] );
                           ( spt->ke[kk][ll] ) = ( ppt->mpt->ke[qq][kk]) ;
                           ( spt->km[kk][ll] ) = ( ppt->mpt->km[qq][kk]) ;
                           break;

                          case 1:
                           ( spt->epr[kk][ll] ) = \
                              ( ppt->mpt->ep[qq][kk+THREE] );
                           ( spt->myr[kk][ll] ) = \
                              ( ppt->mpt->my[qq][kk+THREE] );
                           ( spt->ke[kk][ll] ) = \
                              ( ppt->mpt->ke[qq][kk+THREE] );
                           ( spt->km[kk][ll] ) = \
                              ( ppt->mpt->km[qq][kk+THREE] );
                           break;

                          case 2:
                           ( spt->epr[kk][ll] ) = ( ppt->mpt->ep[qq][kk+FIVE] );
                           ( spt->myr[kk][ll] ) = ( ppt->mpt->my[qq][kk+FIVE] );
                           ( spt->ke[kk][ll] ) = ( ppt->mpt->ke[qq][kk+FIVE] );
                           ( spt->km[kk][ll] ) = ( ppt->mpt->km[qq][kk+FIVE] );
                           break;
                        };
                        ( spt->epr[ll][kk] ) = ( spt->epr[kk][ll] );
                        ( spt->myr[ll][kk] ) = ( spt->myr[kk][ll] );
                        ( spt->ke[ll][kk] ) = ( spt->ke[kk][ll] );
                        ( spt->km[ll][kk] ) = ( spt->km[kk][ll] );

                        ll++;
                     };
                  } while(( ++kk ) < DIMNS );

                  jj = null; do
                  {
                     kk = null; do
                     {
                        ( spt->c[jj][kk] ) = \
                           ( ppt->cpt->c[(tpt->cm[hh][jj])][kk] );
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < CRNRS );

                  ( spt->opt ) = 's';
/*............................................................................*/
                  spt = smatrx( spt );       /* compute s-parameters */
/*.........................................*/
                  if (( spt->rtn ) == null )
                  {
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in smatrx(*), "
                        "option 's', cell number %ld !!!\n", hh ); 

                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;
                  };
               }; /* end if media index != null */
/*............................................................................*/
# if FORM_REPEAT == 1

               if ( null == strncmp(( ppt->domain ), "time_domain", TWO ))
               {
                  if ( hh == ( tpt->mi ))
                     goto copy1;

                  if (( null != strcmp( rpt.etyp, ( spt->etyp )))
                    ||( null != strcmp( rpt.mtyp, ( spt->mtyp ))))
                     goto copy1;
                  else if (( *( spt->etyp ) == 't' )\
                         &&( *( spt->mtyp ) == 't' ))
                     goto lbl_nn; /* write s-matrix label */

                  if ( *rpt.etyp != *( spt->etyp ))
                     goto copy1;

                  if ( *rpt.mtyp != *( spt->mtyp ))
                     goto copy1;

                  if ( *rpt.getp != *( spt->getp ))
                     goto copy1;

                  if ( *rpt.gmtp != *( spt->gmtp ))
                     goto copy1;

                  jj = null; do
                  {
                     kk = null; do
                     {
                        if ( SMALL_VAL < fabs( spt->se[jj][kk] ))
                        {
                           if ( se_bound < fabs( 1. - \
                              (( spt->se[jj][kk] )/( rpt.se[jj][kk] ))))
                              goto copy1;
                        } 
                        else if ( SMALL_VAL < fabs( rpt.se[jj][kk] ))
                        {
                           if ( se_bound < fabs( 1. - \
                              (( rpt.se[jj][kk] )/( spt->se[jj][kk] ))))
                              goto copy1;
                        }; 

                        if ( SMALL_VAL < fabs( spt->sm[jj][kk] ))
                        {
                           if ( sm_bound < fabs( 1. - \
                              (( spt->sm[jj][kk] )/( rpt.sm[jj][kk] ))))
                              goto copy1;
                        } 
                        else if ( SMALL_VAL < fabs( rpt.sm[jj][kk] ))
                        {
                           if ( sm_bound < fabs( 1. - \
                              (( rpt.sm[jj][kk] )/( spt->sm[jj][kk] ))))
                              goto copy1;
                        }; 
                     } while(( ++kk ) < SORDR );
                  } while(( ++jj ) < SORDR );

                  if ( *( spt->getp ) == 'g' )
                  {
                     if ( SMALL_VAL < fabs( spt->no ))
                     {
                        if ( no_bound < fabs( 1.- ( rpt.no/( spt->no ))))
                           goto copy1;
                     } 
                     else if ( SMALL_VAL < fabs( rpt.no ))
                     {
                        if ( no_bound < fabs( 1.- (( spt->no )/ rpt.no )))
                           goto copy1;
                     };

                     if ( SMALL_VAL < fabs( spt->tr ))
                     {
                        if ( tr_bound < fabs( 1.- ( rpt.tr/( spt->tr ))))
                           goto copy1;
                     } 
                     else if ( SMALL_VAL < fabs( rpt.tr ))
                     {
                        if ( tr_bound < fabs( 1.- (( spt->tr )/rpt.tr )))
                           goto copy1;
                     };
  
                     jj = null; do
                     {
                        if ( SMALL_VAL < fabs( spt->mi[jj] ))
                        {
                           if ( mi_bound < fabs( 1. -\
                              ( rpt.mi[jj]/( spt->mi[jj] ))))
                              goto copy1;
                        }
                        else if ( SMALL_VAL < fabs( rpt.mi[jj] ))
                        {
                           if ( mi_bound < fabs( 1. -\
                              (( spt->mi[jj] )/rpt.mi[jj] )))
                              goto copy1;
                        };

                        kk = null; do
                        {
                           if ( SMALL_VAL < fabs( spt->a[jj][kk] ))
                           {
                              if ( aa_bound < fabs( 1. - \
                                 ( rpt.a[jj][kk]/( spt->a[jj][kk] ))))
                                 goto copy1;
                           }
                           else if ( SMALL_VAL < fabs( rpt.a[jj][kk] ))
                           {
                              if ( aa_bound < fabs( 1. - \
                                 (( spt->a[jj][kk] )/rpt.a[jj][kk] )))
                                 goto copy1;
                           };
                        } while(( ++kk ) < DIMNS );
                     } while(( ++jj ) < DIMNS );
                  };

                  if ( *( spt->gmtp ) == 'g' )
                  {
                     if ( SMALL_VAL < fabs( spt->tg ))
                     {
                        if ( tg_bound < fabs( 1. - ( rpt.tg/( spt->tg )))) 
                           goto copy1;
                     }
                     else if ( SMALL_VAL < fabs( rpt.tg ))
                     {
                        if ( tg_bound < fabs( 1. - (( spt->tg )/rpt.tg ))) 
                           goto copy1;
                     };

                     if ( SMALL_VAL < fabs( spt->ld ))
                     {
                        if ( ld_bound < fabs( 1. - ( rpt.ld/( spt->ld )))) 
                           goto copy1;
                     }
                     else if ( SMALL_VAL < fabs( rpt.ld ))
                     {
                        if ( ld_bound < fabs( 1. - (( spt->ld )/rpt.ld ))) 
                           goto copy1;
                     };

                     jj = null; do
                     {
                        if ( SMALL_VAL < fabs( spt->ms[jj] ))
                        {
                           if (  ms_bound < fabs( 1. - \
                              ( rpt.ms[jj]/( spt->ms[jj] ))))
                              goto copy1;
                        }
                        else if ( SMALL_VAL < fabs( rpt.ms[jj] ))
                        {
                           if (  ms_bound < fabs( 1. - \
                              (( spt->ms[jj] )/rpt.ms[jj] )))
                              goto copy1;
                        };

                        if ( SMALL_VAL < fabs( spt->hg[jj] ))
                        {
                           if (  hg_bound < fabs( 1. - \
                              ( rpt.hg[jj]/( spt->hg[jj] ))))
                              goto copy1;
                        }
                        else if ( SMALL_VAL < fabs( rpt.hg[jj] ))
                        {
                           if (  hg_bound < fabs( 1. - \
                              (( spt->hg[jj] )/rpt.hg[jj] )))
                              goto copy1;
                        };

                        kk = null; do
                        {
                           if ( SMALL_VAL < fabs( spt->a[jj][kk] ))
                           {
                              if ( aa_bound < fabs( 1. - \
                                 ( rpt.a[jj][kk]/( spt->a[jj][kk] ))))
                                 goto copy1;
                           }
                           else if ( SMALL_VAL < fabs( rpt.a[jj][kk] ))
                           {
                              if ( aa_bound < fabs( 1. - \
                                 (( spt->a[jj][kk] )/rpt.a[jj][kk] )))
                                 goto copy1;
                           };
                        } while(( ++kk ) < DIMNS );
                     } while(( ++jj ) < DIMNS );
                  };
/*............................................................................*/
/* no changes: former s- and form parameters overtaken with present cell      */

                  if ( *( spt->getp ) == 'g' )
                     ecl++;

                  if ( *( spt->gmtp ) == 'g' )
                     mcl++;

                  goto lbl_nn; /* write s-matrix label */
/*............................................................................*/
/* overtake new s-parameters into shadow structure rpt: */

                 copy1:

                  strcpy( rpt.etyp, ( spt->etyp ));
                  strcpy( rpt.mtyp, ( spt->mtyp ));
                  strcpy( rpt.getp, ( spt->getp ));
                  strcpy( rpt.gmtp, ( spt->gmtp ));

                  jj = null; do
                  {
                     kk = null; do
                     {
                        rpt.se[jj][kk] = ( spt->se[jj][kk] );
                        rpt.sm[jj][kk] = ( spt->sm[jj][kk] );
                     } while(( ++kk ) < SORDR );
                  } while(( ++jj ) < SORDR );

                  if ( *( spt->getp ) == 'g' )
                  {
                     rpt.no = ( spt->no );
                     rpt.tr = ( spt->tr );

                     jj = null; do
                     {
                        rpt.mi[jj] = ( spt->mi[jj] );

                        kk = null; do
                        {
                           rpt.a[jj][kk] = ( spt->a[jj][kk] );
                        } while(( ++kk ) < DIMNS );
                     } while(( ++jj ) < DIMNS );
                  }
                  else /* ( spt->getp ) != 'g' */
                  {
                     rpt.no = ZERO;
                     rpt.tr = ZERO;

                     jj = null; do
                     {
                        rpt.mi[jj] = ZERO;

                        kk = null; do
                        {
                           rpt.a[jj][kk] = ( spt->a[jj][kk] );
                        } while(( ++kk ) < DIMNS );
                     } while(( ++jj ) < DIMNS );
                  };

                  if ( *( spt->gmtp ) == 'g' )
                  {
                     rpt.tg = ( spt->tg );
                     rpt.ld = ( spt->ld );

                     jj = null; do
                     {
                        rpt.ms[jj] = ( spt->ms[jj] );
                        rpt.hg[jj] = ( spt->hg[jj] );

                        kk = null; do
                        {
                           rpt.a[jj][kk] = ( spt->a[jj][kk] );
                        } while(( ++kk ) < DIMNS );
                     } while(( ++jj ) < DIMNS );
                  }
                  else /* if spt->gmtp != 'g' */
                  {
                     rpt.tg = ZERO;
                     rpt.ld = ZERO;

                     jj = null; do
                     {
                        rpt.ms[jj] = ZERO;
                        rpt.hg[jj] = ZERO;

                        kk = null; do
                        {
                           rpt.a[jj][kk] = ( spt->a[jj][kk] );
                        } while(( ++kk ) < DIMNS );
                     } while(( ++jj ) < DIMNS );
                  };
               } /* end if *( ppt->domain ) == 't'ime_domain */
               else /* if *( ppt->domain ) == 'f'requency_domain */
               {
                  if ( hh == ( tpt->mi ))
                     goto copy2;

                  if (( null != strcmp( rpt.etyp, ( spt->etyp )))
                    ||( null != strcmp( rpt.mtyp, ( spt->mtyp ))))
                     goto copy2;
                  else if (( *( spt->etyp ) == 't' )\
                         &&( *( spt->mtyp ) == 't' ))
                     goto lbl_nn; /* write s-matrix label */

                  if ( *rpt.etyp != *( spt->etyp ))
                     goto copy2;

                  if ( *rpt.mtyp != *( spt->mtyp ))
                     goto copy2;

                  if ( *rpt.getp != *( spt->getp ))
                     goto copy2;

                  if ( *rpt.gmtp != *( spt->gmtp ))
                     goto copy2;

                  jj = null; do
                  {
                     kk = null; do
                     {
                        if ( SMALL_VAL < fabs( spt->ser[jj][kk] ))
                        {
                           if ( se_bound < fabs( 1. - \
                              ( rpt.ser[jj][kk]/( spt->ser[jj][kk] ))))
                              goto copy2;
                        }
                        else if ( SMALL_VAL < fabs( rpt.ser[jj][kk] ))
                        {
                           if ( se_bound < fabs ( 1. - \
                              (( spt->ser[jj][kk] )/rpt.ser[jj][kk] )))
                              goto copy2;
                        };

                        if ( SMALL_VAL < fabs( spt->sei[jj][kk] ))
                        {
                           if ( se_bound < fabs( 1. - \
                              ( rpt.sei[jj][kk]/( spt->sei[jj][kk] ))))
                              goto copy2;
                        }
                        else if ( SMALL_VAL < fabs( rpt.sei[jj][kk] ))
                        {
                           if ( se_bound < fabs( 1. - \
                              (( spt->sei[jj][kk] )/rpt.sei[jj][kk] )))
                              goto copy2;
                        };

                        if ( SMALL_VAL < fabs( spt->smr[jj][kk] ))
                        {
                           if ( sm_bound < fabs( 1. - \
                              ( rpt.smr[jj][kk]/( spt->smr[jj][kk] ))))
                              goto copy2;
                        }
                        else if ( SMALL_VAL < fabs( rpt.smr[jj][kk] ))
                        {
                           if ( sm_bound < fabs( 1. - \
                              (( spt->smr[jj][kk] )/rpt.smr[jj][kk] )))
                              goto copy2;
                        };

                        if ( SMALL_VAL < fabs( spt->smi[jj][kk] ))
                        {
                           if ( sm_bound < fabs( 1. - \
                              ( rpt.smi[jj][kk]/( spt->smi[jj][kk] ))))
                              goto copy2;
                        }
                        else if ( SMALL_VAL < fabs( rpt.smi[jj][kk] ))
                        {
                           if ( sm_bound < fabs( 1. - \
                              (( spt->smi[jj][kk] )/rpt.smi[jj][kk] )))
                              goto copy2;
                        };
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < DIMNS );
/*............................................................................*/
/* no changes: former s- and form parameters overtaken with present cell      */

                  goto lbl_nn; /* write s-matrix label */
/*............................................................................*/
/* copy new s-parameters into repeater structure: */

                 copy2:

                  strcpy( rpt.etyp, ( spt->etyp ));
                  strcpy( rpt.mtyp, ( spt->mtyp ));
                  strcpy( rpt.getp, ( spt->getp ));
                  strcpy( rpt.gmtp, ( spt->gmtp ));

                  jj = null; do
                  {
                     kk = null; do
                     {
                        rpt.ser[jj][kk] = ( spt->ser[jj][kk] );
                        rpt.sei[jj][kk] = ( spt->sei[jj][kk] );
                        rpt.smr[jj][kk] = ( spt->smr[jj][kk] );
                        rpt.smi[jj][kk] = ( spt->smi[jj][kk] );
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < DIMNS );
               };

# endif /* end if FORM_REPEAT == 1 */
/*............................................................................*/

/*............................................................................*/
# if SMTRX_LOG == 1

               if ( *( ppt->domain ) == 't' )
               {
                  fprintf( smtrxlog, "\n\nKE[%ld]:\n", hh );
                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .15e% ", ( spt->se[ii][jj] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );

                  fprintf( smtrxlog, "\n\nLE[%ld]:\n", hh );
                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .15e% ",
                              ( spt->se[ii][jj+DIMNS] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );

                  fprintf( smtrxlog, "\n\nME[%ld]:\n", hh );
                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .15e% ",
                              ( spt->se[ii+DIMNS][jj] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );

                  fprintf( smtrxlog, "\n\nNE[%ld]:\n", hh );
                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .15e% ",
                              ( spt->se[ii+DIMNS][jj+DIMNS] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );

                  fprintf( smtrxlog, "\n\nKH[%ld]:\n", hh );
                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .15e% ", ( spt->sm[ii][jj] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );

                  fprintf( smtrxlog, "\n\nLH[%ld]:\n", hh );
                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .15e% ",
                              ( spt->sm[ii][jj+DIMNS] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );

                  fprintf( smtrxlog, "\n\nMH[%ld]:\n", hh );
                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .15e% ",
                              ( spt->sm[ii+DIMNS][jj] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );

                  fprintf( smtrxlog, "\n\nNH[%ld]:\n", hh );
                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .15e% ",
                              ( spt->sm[ii+DIMNS][jj+DIMNS] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );
               } /* end if *(ppt->domain) == 't'ime_domain */
               else /* if *(ppt->domain) == 'f'requency_domain */
               {
                  fprintf( smtrxlog, "\n\nSE[%ld]:\n", hh );
                  ii = null; do 
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .4le%+.4le*j ",
                              ( spt->ser[ii][jj] ), ( spt->sei[ii][jj] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS ); 

                  fprintf( smtrxlog, "\n\nSH[%ld]:\n", hh );

                  ii = null; do
                  {
                     fprintf( smtrxlog, "\n" );
                     jj = null; do
                     {
                        fprintf( smtrxlog, "|% .4le%+.4le*j ",
                              ( spt->smr[ii][jj] ), ( spt->smi[ii][jj] ));
                     } while(( ++jj ) < DIMNS );
                     fprintf( smtrxlog, "|" );
                  } while(( ++ii ) < DIMNS );
               }; /* end if *( ppt->domain ) == 'f'requency_domain */

               if (( FORM_REPEAT == 1 )\
                 &&( hh < ( tpt->mf )))
                  fprintf( smtrxlog, "\n\n*\n*\n*" );

# endif /* SMTRX_LOG == 1 */
/*............................................................................*/
               nn++; /* next s-matrix label */
/*............................................................................*/
# if LBL_CELLS == 1
               fprintf( s_matrix, "%s%s%ld\n", fields, smtrix, nn );
# endif
/*............................................................................*/
               if ( *( spt->getp ) == 'g' )
               {
                  fprintf( s_matrix, spformat, ( spt->getp ));

                  if ( ecl == null )
                  {
                     for ( jj=null; jj<DIMNS; jj++ )
                        fprintf( s_matrix, d_format, ( spt->mi[jj] ));

                     fprintf( s_matrix, d_format, ( spt->no ));
                     fprintf( s_matrix, d_format, ( spt->tr ));
                  };

                  jj = null; do
                  {
                     kk = null; do
                     {
                        fprintf( s_matrix, d_format, ( spt->a[jj][kk] ));
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < DIMNS );

                  ecl++;
               };

               fprintf( s_matrix, spformat, ( spt->etyp ));

               if ( *( spt->etyp ) != 't' )
               {                                  /* ll=0: asymmetric */
                  ll = ( *( spt->etyp ) == 'd' ); /* ll=1: diagonal */
                  ll += 2*( *( spt->etyp ) == 's' ); /* ll=2: symmetric and */
                                                  /* non-diagonal */
                  jj = null; do
                  {
                     if ( null == strncmp(( ppt->domain ),
                        "time_domain", TWO ))
                     {
                        switch( ll )
                        {
                          case 1: 
                           fprintf( s_matrix, d_format, ( spt->se[jj][jj] ));
                           fprintf( s_matrix, d_format,
                              ( spt->se[jj+THREE][jj] ));
                           fprintf( s_matrix, d_format,
                              ( spt->se[jj+THREE][jj+THREE] ));
                           break;

                          case 2:
                           kk = null; do
                           {
                              fprintf( s_matrix, d_format, ( spt->se[jj][kk] ));
                              fprintf( s_matrix, d_format,
                                 ( spt->se[jj+THREE][kk] ));
                              fprintf( s_matrix, d_format,
                                 ( spt->se[jj+THREE][kk+THREE] ));
                           } while(( ++kk ) <= jj );
                           break;
     
                          default:
                           kk = null; do
                           {
                              fprintf( s_matrix, d_format, ( spt->se[jj][kk] ));
                              fprintf( s_matrix, d_format,
                                 ( spt->se[jj+THREE][kk] ));
                              fprintf( s_matrix, d_format,
                                 ( spt->se[jj+THREE][kk+THREE] ));
                           } while(( ++kk ) < DIMNS );
                           break;
                        }; /* end switch(ll) */
                     } /* end if *( ppt->domain ) == 't'ime_domain */
                     else if ( null == strncmp(( ppt->domain ),
                        "frequency_domain", TWO ))
                     {
                        switch( ll )
                        {
                          case 1:
                           fprintf( s_matrix, "%+.16e  ", ( spt->ser[jj][jj] ));
                           fprintf( s_matrix, "%+.16e\n", ( spt->sei[jj][jj] ));
                           break;

                          case 2:
                           kk = null; do
                           {
                              fprintf( s_matrix, "%+.16e  ",
                                 ( spt->ser[jj][kk] ));
                              fprintf( s_matrix, "%+.16e\n",
                                 ( spt->sei[jj][kk] ));
                           } while(( ++kk ) <= jj );
                           break;

                          default:
                           kk = null; do
                           {
                              fprintf( s_matrix, "%+.16e  ",
                                 ( spt->ser[jj][kk] ));
                              fprintf( s_matrix, "%+.16e\n",
                                 ( spt->sei[jj][kk] ));
                           } while(( ++kk ) < DIMNS );
                           break;
                        };
                     }; /* end if *( ppt->domain ) == 'f'requency_domain */
                  } while(( ++jj ) < DIMNS );/* next jj */
               }; /* if ...etyp != 't'rivial */

               if ( *( spt->gmtp ) == 'g' )
               {
                  fprintf( s_matrix, spformat, ( spt->gmtp ));

                  if ( mcl == null )
                  {
                     jj = null; do
                     {
                        fprintf( s_matrix, d_format, ( spt->ms[jj] ));
                        fprintf( s_matrix, d_format, ( spt->hg[jj] ));
                     } while(( ++jj ) < DIMNS );

                     fprintf( s_matrix, d_format, ( spt->ld ));
                     fprintf( s_matrix, d_format, ( spt->tg ));
                  };

                  jj = null; do
                  {
                     kk = null; do
                     {
                        fprintf( s_matrix, d_format, ( spt->a[jj][kk] ));
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < DIMNS );

                  mcl++;
               };

               fprintf( s_matrix, spformat, ( spt->mtyp ));

               if ( *( spt->mtyp ) != 't' )
               {                                  /* ll=0: asymmetric */
                  ll = ( *( spt->mtyp ) == 'd' ); /* ll=1: diagonal */
                  ll += 2*( *( spt->mtyp ) == 's' ); /* ll=2: symmetric and */
                                                  /* non-diagonal */
                  jj = null; do 
                  {
                     if ( null == strncmp(( ppt->domain ),
                        "time_domain", TWO ))
                     {
                        switch( ll )
                        {
                          case 1:
                           fprintf( s_matrix, d_format,
                              ( spt->sm[jj][jj] ));
                           fprintf( s_matrix, d_format,
                              ( spt->sm[jj+THREE][jj] ));
                           fprintf( s_matrix, d_format,
                              ( spt->sm[jj+THREE][jj+THREE] ));
                           break;

                          case 2:
                           kk = null; do
                           {
                              fprintf( s_matrix, d_format,
                                 ( spt->sm[jj][kk] ));
                              fprintf( s_matrix, d_format,
                                 ( spt->sm[jj+THREE][kk] ));
                              fprintf( s_matrix, d_format,
                                 ( spt->sm[jj+THREE][kk+THREE] ));
                           } while(( ++kk ) <= jj );
                           break;

                          default:
                           kk = null; do
                           {
                              fprintf( s_matrix, d_format,
                                 ( spt->sm[jj][kk] ));
                              fprintf( s_matrix, d_format,
                                 ( spt->sm[jj+THREE][kk] ));
                              fprintf( s_matrix, d_format,
                                 ( spt->sm[jj+THREE][kk+THREE] ));
                           } while(( ++kk ) < DIMNS );
                           break;
                        }; /* end switch(ll) */
                     } /* end if *( ppt->domain ) == 't'ime_domain */ 
                     else if ( null == strncmp(( ppt->domain ),
                           "frequency_domain", TWO ))
                     {
                        switch( ll )
                        {
                          case 1:
                           fprintf( s_matrix, "%+.16e  ",
                              ( spt->smr[jj][jj] ));
                           fprintf( s_matrix, "%+.16e\n",
                              ( spt->smi[jj][jj] ));
                           break;

                          case 2:
                           kk = null; do
                           {
                              fprintf( s_matrix, "%+.16e  ",
                                 ( spt->smr[jj][kk] ));
                              fprintf( s_matrix, "%+.16e\n",
                                 ( spt->smi[jj][kk] ));
                           } while(( ++kk ) <= jj );
                           break;

                          default:
                           kk = null; do
                           {
                              fprintf( s_matrix, "%+.16e  ",
                                 ( spt->smr[jj][kk] ));
                              fprintf( s_matrix, "%+.16e\n",
                                 ( spt->smi[jj][kk] ));
                           } while(( ++kk ) < DIMNS );
                           break;
                        };
                     };/* end if *( ppt->domain ) == 'f'requency_domain */
                  } while(( ++jj ) < DIMNS ); /* next jj */
               }; /* end if ( ... htyp != 't'rivial ) */
/*............................................................................*/
# if FORM_REPEAT == 1
              lbl_nn:
# endif
/*............................................................................*/
               ( ppt->rep->smx[hh] ) = nn;
            };
            ( dsp->state ) = ( ++hh );
            dsplay( dsp );
         }; /* next mesh cell index hh */
/*............................................................................*/
         ( ppt->rep->smx[null] ) = nn; /* counter */
/*
         fprintf( s_matrix, "\ns-matrix_number:\t%ld\n",
               ( ppt->rep->smx[null] ));
*/
/*............................................................................*/
/* clear display: */

         printf( CLEAR_LINE );
         printf( "\r Maxwell field s-parameters terminated."\
            "                                        " );
# if DSC_HCRMDE != 0
/*............................................................................*/
/* compute & store current s-parameters: */

         printf( "\n" );

# if LBL_CELLS == 1
         strcpy( ptr, thermal );
         strcat( ptr, smtrix );
# endif
         strcpy(( dsp->messge ),
            "Thermal and fluid s-parameter determination started." );
/*............................................................................*/
/* initialize running cursor: */

         ( dsp->option ) = 's'; /* 's'tart message */
         dsplay( dsp );

         strcpy(( dsp->messge ),
            "[Please wait a moment]" );

         ( dsp->option ) = 'm'; /* 'm'essage under running cursor */
         dsplay( dsp );

         ( dsp->option ) = 'c'; /* running 'c'ursor effect */
         ( dsp->range ) = ( tpt->mf ) - ( tpt->mi ) + ONE;
/*............................................................................*/
/* compute and store thermal - fluid s-parameters: */

         nn = null;
         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            if ( null < ( ppt->rhp->smx[hh] ))
            {
               ( ppt->rhp->smx[hh] ) =
                     ( ppt->rhp->smx[( ppt->rhp->smx[hh] )] );
/*............................................................................*/
# if DSC_FLDMDE != 0
               ( ppt->fcp->cnn[hh] ) =
                     ( ppt->fcp->cnn[( ppt->rhp->smx[hh] )] );
# endif
/*............................................................................*/
            }
            else
            {
               qq = ( ppt->mpt->idx[hh] ); /* media index qq */

               if ( qq == null )
               {
                  ( hsp->opt ) = null;
/*............................................................................*/
                  hsp = hcrsmx( hsp );             /* reset parameters */
/*...............................................*/
               }
               else /* non trivial medium */
               {
/* media index */
                  ( hsp->med ) = qq;
		  
/* thermal properties: */
                  ( hsp->cv ) = ( ppt->mpt->cv[qq] );
                  ( hsp->kh ) = ( ppt->mpt->kh[qq] );

/* electromagnetic properties [lossy medium]: */
                  ( hsp->ke ) = ( ppt->mpt->ke[qq][null] );
                  ( hsp->km ) = ( ppt->mpt->km[qq][null] );
/*............................................................................*/
# if DSC_FLDMDE != 0 
/* fluid properties: */
                  ( hsp->rm ) = ( ppt->mpt->rm[qq] );
                  ( hsp->tm ) = ( ppt->mpt->tm[qq] );
                  ( hsp->bm ) = ( ppt->mpt->bm[qq] );
                  ( hsp->cm ) = ( ppt->mpt->cm[qq] );
                  ( hsp->ny ) = ( ppt->mpt->ny[qq] );
                  ( hsp->q1 ) = ( ppt->mpt->q1[qq] );
                  ( hsp->td ) = ( ppt->mpt->td[qq] );
/*............................................................................*/
/* Prandtl turbulence models */
# if (( TURBMOD == 1 ) \
    ||( TURBMOD == 2 ))
                  ( hsp->LL ) = ( ppt->mpt->LL[qq] );
# endif /* TURBMOD != null */
/*............................................................................*/
                  kk = null; do
                  {
                     ( hsp->gr[kk] ) = ( ppt->mpt->gr[qq][kk] );
                     ( hsp->gp[kk] ) = ( ppt->mpt->gp[qq][kk] );
                  } while(( ++kk ) < THREE );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# if CELL_POINTS
                  if ( hh == CELL_POINTS )
                     printf( "\n\n Corner points of cell no. %ld\n", hh );
# endif
/*............................................................................*/
/* vertex point coordinates: */

                  jj = null; do
                  {
                     kk = null; do
                     {
                        ( hsp->c[jj][kk] ) = \
                           ( ppt->cpt->c[( tpt->cm[hh][jj] )][kk] );
/*............................................................................*/
# if CELL_POINTS
                        if ( hh == CELL_POINTS )
                        {
                           printf( "\n vertex%1ld:    cpt->c[%06ld][%d] =  "
                              "%+.14e", jj, ( tpt->cm[hh][jj] ), kk,
                              ( ppt->cpt->c[( tpt->cm[hh][jj] )][kk] ));
                        };
# endif
/*............................................................................*/
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < CRNRS );
/*............................................................................*/
# if CELL_POINTS
                  if ( hh == CELL_POINTS )
                  {
                     printf( "\n\n Please acknowledge "
                        "( enter any character ) ...........: " );
                     scanf( "%s", ptr );
                     printf( "\n " );
                  };
# endif
/*............................................................................*/
/* call s-parameter function: */

                  ( hsp->dt ) = ( state->hcdt );
                  ( hsp->opt ) = 's'; /* option: s-parameters */
/*............................................................................*/
                  hsp = hcrsmx( hsp );           /* compute s-parameters */
/*.............................................*//* for given time step */
                  if (( hsp->rtn ) == null )
                  {
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in function hcrsmx(*), "
                        "option 't', cell number %ld !!!\n", hh );

                     ( csp->dfopt ) = null; /* next default menu option */

                     goto menu;
                  };
               }; /* end if media index qq != null */
/*............................................................................*/
# if FORM_REPEAT == 1

               if ( ONE ) /* time domain [ here allways ] */
               {
                  if ( hh == ( tpt->mi ))
                     goto copy3;

                  if ( null != strcmp( rhc.ttyp, ( hsp->ttyp )))
                     goto copy3;
                  else if ( *( hsp->ttyp ) == 't' )
                     goto lbl_nn1; /* write s-matrix label */
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( ppt->fcp->cnn[hh] ) != ( ppt->fcp->cnn[hh-ONE] ))
                     goto copy3;
# endif
/*............................................................................*/
                  if ( rhc.med != ( hsp->med ))
                     goto copy3;

                  if ( SMALL_VAL < fabs( hsp->ct ))
                  {
                     if ( hct_bound < fabs( 1. - \
                        (( rhc.ct )/( hsp->ct ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.ct ))
                  {
                     if ( hct_bound < fabs( 1. - \
                        (( hsp->ct )/rhc.ct )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->vol ))
                  {
                     if ( hvl_bound < fabs( 1. - \
                        ( rhc.vol/( hsp->vol ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.vol ))
                  {
                     if ( hvl_bound < fabs( 1. - \
                        (( hsp->vol )/rhc.vol )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->kh ))
                  {
                     if ( hkh_bound < fabs( 1. - \
                        (( rhc.kh )/( hsp->kh ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.kh ))
                  {
                     if ( hkh_bound < fabs( 1. - \
                        (( hsp->kh )/rhc.kh )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->cv ))
                  {
                     if ( hcv_bound < fabs( 1. - \
                        (( rhc.cv )/( hsp->cv ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.cv ))
                  {
                     if ( hcv_bound < fabs( 1. - \
                        (( hsp->cv )/rhc.cv )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->ke ))
                  {
                     if ( hke_bound < fabs( 1. - \
                        (( rhc.ke )/( hsp->ke ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.ke ))
                  {
                     if ( hke_bound < fabs( 1. - \
                        (( hsp->ke )/rhc.ke )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->km ))
                  {
                     if ( hkm_bound < fabs( 1. - \
                        (( rhc.km )/( hsp->km ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.km ))
                  {
                     if ( hkm_bound < fabs( 1. - \
                        (( hsp->km )/rhc.km )))
                        goto copy3;
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if ( SMALL_VAL < fabs( hsp->ft ))
                  {
                     if ( hft_bound < fabs( 1. - \
                        (( rhc.ft )/( hsp->ft ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.ft ))
                  {
                     if ( hft_bound < fabs( 1. - \
                        (( hsp->ft )/rhc.ft )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->rm ))
                  {
                     if ( hrm_bound < fabs( 1. - \
                        (( rhc.rm/hsp->rm ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.rm ))
                  {
                     if ( hrm_bound < fabs( 1. - \
                        (( hsp->rm )/rhc.rm )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->tm ))
                  {
                     if ( htm_bound < fabs( 1. - \
                        (( rhc.tm/hsp->tm ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.tm ))
                  {
                     if ( htm_bound < fabs( 1. - \
                        (( hsp->tm )/rhc.tm )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->bm ))
                  {
                     if ( hbm_bound < fabs( 1. - \
                        (( rhc.bm/hsp->bm ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.bm ))
                  {
                     if ( hbm_bound < fabs( 1. - \
                        (( hsp->bm )/rhc.bm )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->cm ))
                  {
                     if ( hcm_bound < fabs( 1. - \
                        (( rhc.cm/hsp->cm ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.cm ))
                  {
                     if ( hcm_bound < fabs( 1. - \
                        (( hsp->cm )/rhc.cm )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->ny ))
                  {
                     if ( hny_bound < fabs( 1. - \
                        (( rhc.ny/hsp->ny ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.ny ))
                  {
                     if ( hny_bound < fabs( 1. - \
                        (( hsp->ny )/rhc.ny )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->q1 ))
                  {
                     if ( hq1_bound < fabs( 1. - \
                        (( rhc.q1/hsp->q1 ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.q1 ))
                  {
                     if ( hq1_bound < fabs( 1. - \
                        (( hsp->q1 )/rhc.q1 )))
                        goto copy3;
                  };

                  if ( SMALL_VAL < fabs( hsp->td ))
                  {
                     if ( htd_bound < fabs( 1. - \
                        (( rhc.td/hsp->td ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.td ))
                  {
                     if ( htd_bound < fabs( 1. - \
                        (( hsp->td )/rhc.td )))
                        goto copy3;
                  };
/*............................................................................*/
/* Prandtl turbulence models */
# if (( TURBMOD == 1 ) \
    ||( TURBMOD == 2 ))
                  if ( SMALL_VAL < fabs( hsp->LL ))
                  {
                     if ( hLL_bound < fabs( 1. - \
                        (( rhc.LL/hsp->LL ))))
                        goto copy3;
                  }
                  else if ( SMALL_VAL < fabs( rhc.LL ))
                  {
                     if ( hLL_bound < fabs( 1. - \
                        (( hsp->LL )/rhc.LL )))
                        goto copy3;
                  };
# endif /* TURBMOD != null */
/*............................................................................*/
                  kk = null; do
                  {
                     if ( SMALL_VAL < fabs( hsp->gr[kk] ))
                     {
                        if ( hgr_bound < fabs( 1. - \
                           ( rhc.gr[kk]/( hsp->gr[kk] ))))
                           goto copy3;
                     }
                     else if ( SMALL_VAL < fabs( rhc.gr[kk] ))
                     {
                        if ( hgr_bound < fabs( 1. - \
                           (( hsp->gr[kk] )/rhc.gr[kk] )))
                           goto copy3;
                     };

                     if ( SMALL_VAL < fabs( hsp->gp[kk] ))
                     {
                        if ( hgp_bound < fabs( 1. - \
                           ( rhc.gp[kk]/( hsp->gp[kk] ))))
                           goto copy3;
                     }
                     else if ( SMALL_VAL < fabs( rhc.gp[kk] ))
                     {
                        if ( hgp_bound < fabs( 1. - \
                           (( hsp->gp[kk] )/rhc.gp[kk] )))
                           goto copy3;
                     };
                  } while(( ++kk ) < DIMNS );
/* adj(B): */
/*		  
                  jj = null; do
                  {
                     kk = null; do
                     {
                        if ( SMALL_VAL < fabs( hsp->b[jj][kk] ))
                        {
                           if ( hcb_bound < fabs( 1. - \
			      ( rhc.b[jj][kk]/( hsp->b[jj][kk] ))))
                              goto copy3;
                        }
                        else if ( SMALL_VAL < fabs( rhc.b[jj][kk] ))
                        {
                           if ( hcb_bound < fabs( 1. - \
                              (( hsp->b[jj][kk] )/rhc.b[jj][kk] )))
                              goto copy3;
                        };
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < DIMNS );
*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
/* further geometric parameters */
/* area [ geometric size ] of face jj: */

                  jj = null; do
                  {
                     if ( SMALL_VAL < fabs( hsp->fm[jj] ))
                     {
                        if ( hcf_bound < fabs( 1. - \
                           ( rhc.fm[jj]/( hsp->fm[jj] ))))
                           goto copy3;
                     }
                     else if ( SMALL_VAL < fabs( rhc.fm[jj] ))
                     {
                        if ( hcf_bound < fabs( 1. - \
                           (( hsp->fm[jj] )/rhc.fm[jj] )))
                           goto copy3;
                     };
                  } while(( ++jj ) < FACES );

/* [ extrerior ] cell face vectors: */
		  
                  jj = null; do
                  {
                     kk = null; do
                     {
                        if ( SMALL_VAL < fabs( hsp->f[jj][kk] ))
                        {
                           if ( hcf_bound < fabs( 1. - \
                              ( rhc.f[jj][kk]/( hsp->f[jj][kk] ))))
                              goto copy3;
                        }
                        else if ( SMALL_VAL < fabs( rhc.f[jj][kk] ))
                        {
                           if ( hcf_bound < fabs( 1. - \
                              (( hsp->f[jj][kk] )/rhc.f[jj][kk] )))
                              goto copy3;
                        };
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < FACES );
/* adj(B^-1): */
                  jj = null; do
                  {
                     kk = null; do
                     {
                        if ( SMALL_VAL < fabs( hsp->bi[jj][kk] ))
                        {
                           if ( hbi_bound < fabs( 1. - \
                              ( rhc.bi[jj][kk]/( hsp->bi[jj][kk] ))))
                              goto copy3;
                        }
                        else if ( SMALL_VAL < fabs( rhc.bi[jj][kk] ))
                        {
                           if ( hbi_bound < fabs( 1. - \
                              (( hsp->bi[jj][kk] )/rhc.bi[jj][kk] )))
                              goto copy3;
                        };
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < DIMNS );

/* s[jj][] = F[jj] * adj(B^-1): */

                  jj = null; do
                  {
                     kk = null; do
                     {
                        if ( SMALL_VAL < fabs( hsp->s[jj][kk] ))
                        {
                           if ( hcs_bound < fabs( 1. - \
                              ( rhc.s[jj][kk]/( hsp->s[jj][kk] ))))
                              goto copy3;
                        }
                        else if ( SMALL_VAL < fabs( rhc.s[jj][kk] ))
                        {
                           if ( hcs_bound < fabs( 1. - \
                              (( hsp->s[jj][kk] )/rhc.s[jj][kk] )))
                              goto copy3;
                        };
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < FACES );
/*............................................................................*/
/* no changes: former s- and form parameters overtaken for present cell: */

                  goto lbl_nn1; /* write s-matrix label */
/*............................................................................*/
/* copy new s-parameters into shadow structure rhc: */

                 copy3:

                  strcpy( rhc.ttyp, ( hsp->ttyp ));

                  rhc.med = ( hsp->med );
                  rhc.ct = ( hsp->ct );
                  rhc.vol = ( hsp->vol );

                  rhc.cv = ( hsp->cv );
                  rhc.kh = ( hsp->kh );
                  rhc.ke = ( hsp->ke );
                  rhc.km = ( hsp->km );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if ( NFCNN < ( ppt->fcp->cnn[ii] ))
                  {
                     fprintf( stderr, "\n\n Message from function %s:",
                         __func__ );
                     fprintf( stderr, "\n\n Too many fluid connected "
                        "components defined in DSC system %s !!!", 
                           ( tpt->name ));
                     fprintf( stderr, "\n [ Maximum number is %d = macro "
                        "NFCNN fixed in file %s.", NFCNN, "FORMER.CONF" );
                     fprintf( stderr, "\n - Change macro only in compliance "
                        "with memory resources " );
                     fprintf( stderr, "\n   and same macro in configuration "
                        "file 'SOLVER.CONF' of program 'solver(*)'.] " );

                     exit( EXIT_FAILURE ); /* abnormal return */
                  };

                  rhc.ft = ( hsp->ft );
                  rhc.rm = ( hsp->rm );
                  rhc.tm = ( hsp->tm );
                  rhc.bm = ( hsp->bm );
                  rhc.cm = ( hsp->cm );
                  rhc.ny = ( hsp->ny );
                  rhc.q1 = ( hsp->q1 );
                  rhc.td = ( hsp->td );
/*............................................................................*/
/* Prandtl turbulence models */
# if (( TURBMOD == 1 ) \
    ||( TURBMOD == 2 ))
                  rhc.LL = ( hsp->LL );
# endif /* TURBMOD != null */
/*............................................................................*/
                  kk = null; do
                  {
                     rhc.gr[kk] = ( hsp->gr[kk] );
                     rhc.gp[kk] = ( hsp->gp[kk] );
                  } while(( ++kk ) < DIMNS );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
                  jj = null; do
                  {
                     rhc.fm[jj] = ( hsp->fm[jj] );
                  } while(( ++jj ) < FACES );

                  jj = null; do
                  {
                     kk = null; do
                     {
                        rhc.f[jj][kk] = ( hsp->f[jj][kk] );
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < FACES );
/*............................................................................*/
# if DSC_FLDMDE != 0
/*
                  jj = null; do
                  {
                     kk = null; do
                     {
                        rhc.b[jj][kk] = ( hsp->b[jj][kk] );
                     } while(( ++kk ) < DIMNS ); 
                  } while(( ++jj ) < DIMNS ); 
*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
                  jj = null; do
                  {
                     kk = null; do
                     {
                        rhc.bi[jj][kk] = ( hsp->bi[jj][kk] );
                     } while(( ++kk ) < DIMNS ); 
                  } while(( ++jj ) < DIMNS ); 

                  jj = null; do
                  {
                     kk = null; do
                     {
                        rhc.s[jj][kk] = ( hsp->s[jj][kk] );
                     } while(( ++kk ) < DIMNS );
                  } while(( ++jj ) < FACES );
               } /* end if *( domain ) == 't'ime_domain [ void condition ] */
# endif /* end if FORM_REPEAT == 1 */
/*............................................................................*/

               nn++; /* next s-matrix label */

/*............................................................................*/
               if ( HCSMX < nn )
               {
                  fprintf( stderr,
                     "\n\n Message from function %s:", __func__ );
                  fprintf( stderr,
                     "\n Too many thermal&fluid s-parameter sets defined "
                     "in DSC mesh %s !!!", ( ppt->name ));
                  fprintf( stderr,
                     "\n [ Maximum number is %ld = macro HCSMX "
                     "in file %s.", ( long ) HCSMX, "FORMER.CONF" );
                  fprintf( stderr,
                     "\n - Change macro only in compliance "
                     "with memory resources.]\n" );

                  exit( EXIT_FAILURE ); /* abnormal return */
               };
/*............................................................................*/
# if LBL_CELLS == 1
               strcpy( ptr, thermal );
               strcat( ptr, smtrix );
               fprintf( s_matrix, "%s%ld\n", ptr, nn );
# endif
/*............................................................................*/
               fprintf( s_matrix, spformat, ( hsp->ttyp ));

               if ( *( hsp->ttyp ) != 't' ) /* non trivial cell */
               { 
/*............................................................................*/
/* 'normalized' thermal time step */
                  fprintf( s_matrix, d_format, ( hsp->ct ));

/* cell volume [m^3]: */
                  fprintf( s_matrix, d_format, ( hsp->vol ));
/*............................................................................*/
/* heat conductivity [ W/(m*K) ] and heat capacity per volume [ J/(K*m^3) ] */
                  fprintf( s_matrix, d_format, ( hsp->kh ));
                  fprintf( s_matrix, d_format, ( hsp->cv ));
/*............................................................................*/
/* [ exterior ] cell face vectors: */

                  ll = ONE;
                  jj = null; do
                  {
                     ll *= ( -ONE );
                     kk = null; do
                     {
                        fprintf( s_matrix, "%+.16e ",
                           ((( double ) ll )*( hsp->f[jj][kk] )));
                     } while(( ++kk ) < DIMNS );
                     fprintf( s_matrix, "\n" );
                  } while(( ++jj ) < FACES );
/*............................................................................*/
/* cell faces [m^2]: */
                  jj = null; do
                  {
                     fprintf( s_matrix, "%+.16e ", ( hsp->fm[jj] ));
                  } while(( ++jj ) < THREE );
                  fprintf( s_matrix, "\n" );
                  do
                  {
                     fprintf( s_matrix, "%+.16e ", ( hsp->fm[jj] ));
                  } while(( ++jj ) < FACES );
                  fprintf( s_matrix, "\n" );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if ( null == strncmp( hsp->ttyp, "fluid", THREE ))
                  {
                     fprintf( s_matrix, d_format, ( hsp->ft ));
                     fprintf( s_matrix, d_format, ( hsp->rm ));
                     fprintf( s_matrix, d_format, ( hsp->tm ));
                     fprintf( s_matrix, d_format, ( hsp->bm ));
                     fprintf( s_matrix, d_format, ( hsp->cm ));
                     fprintf( s_matrix, d_format, ( hsp->ny ));
                     fprintf( s_matrix, d_format, ( hsp->q1 ));
                     fprintf( s_matrix, d_format, ( hsp->td ));
/*............................................................................*/
/* Prandtl turbulence models */
# if (( TURBMOD == 1 ) \
    ||( TURBMOD == 2 ))
                     fprintf( s_matrix, d_format, ( hsp->LL ));
# endif /* TURBMOD != null */
/*............................................................................*/
                     kk = null; do
                     {
                        fprintf( s_matrix, "%+.16e ", ( hsp->gr[kk] ));
                     } while(( ++kk ) < DIMNS );
                     fprintf( s_matrix, "\n" );

                     kk = null; do
                     {
                        fprintf( s_matrix, "%+.16e ", ( hsp->gp[kk] ));
                     } while(( ++kk ) < DIMNS );
                     fprintf( s_matrix, "\n" );
/*............................................................................*/
/* [ exterior ] cell face vectors: */
/*
                     ll = ONE;
                     jj = null; do
                     {
		        ll *= ( -ONE );
                        kk = null; do
                        {
                           fprintf( s_matrix, "%+.16e ",
                              ((( double ) ll )*( hsp->f[jj][kk] )));
                        } while(( ++kk ) < DIMNS );
                        fprintf( s_matrix, "\n" );
                     } while(( ++jj ) < FACES );
*/
/*............................................................................*/
/* matrix adj( B ): */
/*
                     jj = null; do
                     {
                        kk = null; do
                        {
                           fprintf( s_matrix, "%+.16e ",
                              ( hsp->b[jj][kk] ));
                        } while(( ++kk ) < DIMNS );
                        fprintf( s_matrix, "\n" );
                     } while(( ++jj ) < DIMNS );
*/
                  }; /* end if ( *( hsp->ttyp ) == 'f'luid ) */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
/* matrix adj( B^-1 ): */

                  jj = null; do
                  {
                     kk = null; do
                     {
                        fprintf( s_matrix, "%+.16e ",
                           ( hsp->bi[jj][kk] ));
                     } while(( ++kk ) < DIMNS );
                     fprintf( s_matrix, "\n" );
                  } while(( ++jj ) < DIMNS );
/*............................................................................*/
/* form vectors  s = F * adj(B^-1): */

                  jj = null; do
                  {
                     kk = null; do
                     {
                        fprintf( s_matrix, "%+.16e ",
                           ( hsp->s[jj][kk] ));
                     } while(( ++kk ) < DIMNS );
                     fprintf( s_matrix, "\n" );
                  } while(( ++jj ) < FACES );

                  if (( null != strrchr(( hsp->ttyp ), (int) 'E' ))\
                    &&( null != strrchr(( hsp->ttyp ), (int) 'M' )))
                  { /* electric sources */
                     fprintf( s_matrix, d_format, ( hsp->ke ));
                     fprintf( s_matrix, d_format, ( hsp->km ));
                  }
		  else if ( null != strrchr(( hsp->ttyp ), (int) 'M' ))
                  { /* magnetic sources */
                     fprintf( s_matrix, d_format, ( hsp->km ));
                  }
		  else if ( null != strrchr(( hsp->ttyp ), (int) 'E' ))
                  { /* magnetic sources */
                     fprintf( s_matrix, d_format, ( hsp->ke ));
                  };
               }; /* end if *( hsp->ttyp ) != 't'rivial */
/*............................................................................*/
# if FORM_REPEAT == 1

              lbl_nn1:
# endif
/*............................................................................*/
               ( ppt->rhp->smx[hh] ) = nn;
            };
            ( dsp->state ) = ( ++hh );
            dsplay( dsp );
         }; /* next mesh cell index hh */

         ( ppt->rhp->smx[null] ) = nn; /* counter */
/*............................................................................*/
/* clear display: */

         printf( CLEAR_LINE );
         printf( "\r Thermal s-parameters terminated.   "\
            "                                           " );
/*............................................................................*/
# endif
/* [ smxidx = ">>-S-PARAMETER_LABELS->>", sprmts = "smx-identifiers", e.g.] */

         fprintf( s_matrix, "\n%s\n\n", smxidx );
         fprintf( s_matrix, "%-15s %-15s", "cell", "electric" );
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE == 0
         fprintf( s_matrix, " %-15s", "thermal" );
# elif DSC_FLDMDE != 0 
         fprintf( s_matrix, " %-15s", "thermal_&_fluid" );
         fprintf( s_matrix, " %-s", "[fluid_connected_component]" );
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         fprintf( s_matrix, "\n" );

         hh = ( tpt->mi );
         while( hh <= ( tpt->mf ))
         {
            fprintf( s_matrix, "%-15ld %-15ld", hh, ( ppt->rep->smx[hh] ));
/*............................................................................*/
# if DSC_HCRMDE != 0
/* thermal node s-matrix [ index ] */
            fprintf( s_matrix, " %-15ld", ( ppt->rhp->smx[hh] ));

# if DSC_FLDMDE != 0
/* fluid connection component [ index cnn ] */
            if (( null < ( ppt->fcp->cnn[null] ))\
              &&( null < ( ppt->fcp->cnn[hh] )))
               fprintf( s_matrix, " %-15d", ( ppt->fcp->cnn[hh] ));
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
            fprintf( s_matrix, "\n" );

            hh++ ;
         }; /* next hh */
/*............................................................................*/
/* close s-parameter file with date & time of creation: */ 

         nseconds = time( timer );
/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); '*/

         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( s_matrix, "\nDSC model s-parameters file %s ", fleptr );
         fprintf( s_matrix, "created:\n%s\n", tmestr );

         fclose( s_matrix );

         printf( "\n S-matrix file %s", fleptr );
/*............................................................................*/
# if SMTRX_LOG == 1

         fprintf( smtrxlog, "\n\nLogfile %s ", smxfle );
         fprintf( smtrxlog, "created:\n%s", tmestr );

         fclose( smtrxlog );
         printf( ",\n S-matrix log file %s", smxfle );
# endif
/*............................................................................*/
         printf( time_frm, tmestr );

         if (( state->item ) == FORM_ITEMS )
            printf( "\n ----------------------------------"
               "--------------------------------------------" );
         break; 
/*............................................................................*/
        case 3: /* >---------- boundary file generation -------------------> */

         if (( state->item ) == FORM_ITEMS )
            printf( "\n\n Boundary file generation: " );

         strcpy( fleptr, prefix );
         strcat( fleptr, bndptr );

         if (( state->stat ) == FORM_ITEMS )
         {
            strcpy(( bpt->name ), ( tpt->name ));
            strcat( fleptr, ( state->flbl ));
            printf( "\n\n Comment: %.70s \n", ( bpt->text ));
         }
         else if (( state->stat ) != FORM_ITEMS )
         {
            if (( state->item ) < FORM_ITEMS )
            {
               printf( "\n" );

               strcpy(( csp->rqlng ),
                  "Please enter index <N> of boundary file '" );
               strcat(( csp->rqlng ), fleptr );
               strcat(( csp->rqlng ), "<N>'" );
               strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
               csp = txcnsl( csp );            /* enter index N on console    */
/*...........................................*/
               strcpy( ptr, ( lotos(( csp->inlng ), null )));
               strncpy(( state->flbl ), ptr, THREE );
            }
            else 
               printf( "\n" ); 

            strcat( fleptr, ( state->flbl ));

            strcpy(( csp->rqfrm ), "brackets" );
            strcpy(( csp->rqstr ), "Input on keyboard or subroutine " );
            strcat(( csp->rqstr ), "sysbnd(*) ? [k/s]" );
            strcpy(( csp->dfstr ), "s" ); /* default */
/*............................................................................*/
            csp = txcnsl( csp );      /* enter string on text console         */
/*..................................*/
            strcpy( ptr, ( csp->instr ));

            if (( *ptr == 's' )||( *ptr == 'S' ))
            {
/*............................................................................*/
# if FORM_DISP == 1
               printf( "\n DSC model boundary generation "
                  "functions started.               " );
               PRBLDCLR( "\r" );
               printf( "\n %*s", 78, "DSC-FORMER" );
               PRNORMAL( "" );
# endif
/*............................................................................*/
               if (( state->item ) < FORM_ITEMS )
               {
/*............................................................................*/
                  tpt = systop( state );            /* transfer: check DSC -  */
/*................................................*//* mesh topology          */
                  switch( tpt->rtn )
                  {
                    case -ONE:
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;

                    case ONE:
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in function %s !!!\n", "systop(*)" );
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;
                      
                    case THREE: /* = ( state->item ) [ boundaries ] */
                     ( state->stat ) = ( state->item );
                     ( state->tpt ) = tpt;
/*............................................................................*/
# if FORM_DISP == 1
                     printf( "\n ==================================="
                        "===========================================");
                     PRBLDCLR( "\r" );
                     printf( "\n %*s", 78, "DSC-FORMER" );
                     PRNORMAL( "\r" );
                     printf( "\n DSC model boundary function "\
                        "%s terminated.               \n", "sysbnd(*)" );
# endif
/*............................................................................*/
                     goto identf_check3;
                    
                    default:
                     ( state->tpt ) = tpt;
                     break; 
                  };
/*............................................................................*/
                  ppt = syssmx( state );             /* DSC mesh coordinates, */
/*.................................................*//* media etc.            */
                  switch( ppt->rtn )
                  {
                    case -ONE:
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;

                    case ONE:
                     printf( "\n\n Message from function %s: ", __func__ );
                     printf( "\n Error in function %s !!!\n", "syssmx(*)"  );
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;

                    default:
                     ( state->ppt ) = ppt;
                     break;
                  };
               };

               ( bpt->p ) = null;
               ( bpt->n ) = null;
/*............................................................................*/
# if DSC_HCRMDE != 0
/* reset temperature boundary face counter */
               ( bpt->ntf ) = null;
/* temperature node counter */
               ( bpt->ntn ) = null;
/* heat current boundary faces counter */
               ( bpt->nhc ) = null;
/* surface heat conducting faces counter */
               ( bpt->nsc ) = null;
/* heat source boundary faces counter */
               ( bpt->nsk ) = null;
/*............................................................................*/
# if DSC_FLDMDE != 0
               ( bpt->nns ) = null;
               ( bpt->nsl ) = null;
               ( bpt->nif ) = null;
               ( bpt->nof ) = null;
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
               ( bpt->rtn ) = -TWO;
/*............................................................................*/
               bpt = sysbnd( state );            /* DSC mesh boundary         */
/*.............................................*//* configuration             */

               switch( bpt->rtn )
               {
                 case -TWO:
                  printf( "\n Message from function %s:      ", __func__ );
                  printf( "\n Unsatisfied code in function %s "
                     "!!! \n", "sysbnd(*)" );
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 case -ONE:
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 case ONE:
                  printf( "\n\n Message from function %s:     ", __func__ );
                  printf( "\n Error in function %s !!!\n", "sysbnd(*)" );
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 default:
# if FORM_DISP == 1
                  printf( "\n ==================================="
                     "===========================================" );
                  PRBLDCLR( "\r" );
                  printf( "\n %*s", 78, "DSC-FORMER" );
                  PRNORMAL( "\r" );
                  printf( "\n DSC model boundary function "\
                     "%s terminated.               ", "sysbnd(*)" );
# endif
                  break;
               };

              identf_check3:

               if ( null != strncmp(( tpt->name ), ( bpt->name ), FOUR ))
               {
                  printf( "\n\n Message from function %s: ", __func__ );
                  printf( "\n Incompatible DSC system identifiers "
                     "on functions %s and %s.",
                     "systop(*)", "syssmx(*)" );
                  printf( "\n Please check functions and restart program ! " );
                  return ONE;
               };

# if FORM_DISP == 1
               if (( state->item ) < FORM_ITEMS )
               {
                  printf( "\n\n DSC system identifier [ structure name ]: "
                     "%.22s", ( bpt->name ));
                  printf( "\n Comment: %.70s \n", ( bpt->text ));
               };
# endif
               if ( BNDAP < ( bpt->n ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many aperiodic boundary ports "
                     "in DSC mesh %s !!!", tpt->name );
                  printf( "\n [ Maximum number is %d = macro BNDAP "
                     "in %s. ", BNDAP, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function 'scattr(*)' "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( BNDPR < ( bpt->p ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many periodic boundary ports "
                     "in DSC mesh %s !!!", tpt->name );
                  printf( "\n [ Maximum number is %d = macro BNDPR "
                     "in %s.", BNDPR, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scattr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };
/*............................................................................*/
# if DSC_HCRMDE != 0
               if ( BNDSK < ( bpt->nsk ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many [skin effekt] heat source faces "
                     "defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDSK "
                     "in %s.", BNDSK, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( BNDHC < ( bpt->nhc ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many heat current boundary faces "
                     "defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDHC "
                     "in %s.", BNDHC, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( BNDSC < ( bpt->nsc ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many surface heat conducting "
                     "faces defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDSC "
                     "in %s.", BNDSC, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( BNDTF < ( bpt->ntf ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many temperature boundary faces "
                     "defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDTF "
                     "in %s.", BNDTF, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( BNDTN < ( bpt->ntn ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many temperature boundary nodes "
                     "defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDTN "
                     "in %s.", BNDTN, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };
/*............................................................................*/
# if DSC_FLDMDE != 0

               if ( BNDNS < ( bpt->nns ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many no-slip boundary faces "
                     "defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDNS "
                     "in %s.", BNDNS, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( BNDSL < ( bpt->nsl ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many free slip boundary faces "
                     "defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDSL "
                     "in %s.", BNDSL, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( BNDIF < ( bpt->nif ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many inflow boundary faces "
                     "defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDIF "
                     "in %s.", BNDIF, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( BNDOF < ( bpt->nof ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many outflow boundary faces "
                     "defined in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro BNDOF "
                     "in %s.", BNDOF, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scathcr(*) "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
            }/* end if( *ptr == 's'ubroutine ) */
            else /* if( *ptr != 's'): keyboard input */
            {
               printf( "\n\n Please enter: >----------------------"
                  "------------------------------------>" );
               printf( "\n -> DSC system identifier ( structure name ) "
                  "......: " );
               scanf( "%s", ( bpt->name ));
               printf( " -> Comment ................................"
                  ".......: " );
               scanf( "%s", ( bpt->text ));

           /* boundary_faces: */

               printf( "\n Boundary input - please enter: --------"
                  "---------------------------------->");

              aperiodic:

               printf( "\n number of aperiodic boundary faces"
                       " ...............: " );
               scanf( "%s", ptr );
               ( bpt->n ) = strtol( ptr, endp, DEC );

               if ( BNDAP < ( bpt->n ))
               {
                  printf( "\n\n Too many aperiodic boundary "
                     "ports in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n Maximum number is %d = macro BNDAP "
                     "in %s.", BNDAP, __func__ );
                  printf( "\n [ Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scattr(*) "
                     "of program SOLVER.C.]\n" );
                  goto aperiodic;
               };

               ii = null;
               while( ii < ( bpt->n ))
               {
	          printf( " \n");
                  printf( " mesh cell index of %5ld. face ......."
                     ".............: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->m[ii] ) = strtol( ptr, endp, DEC );
                  printf( " face index of %5ld. face [ 0 <= i <= 5 ] "
                     ".........: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->f[ii] ) = strtol( ptr, endp, DEC );
                  printf( " reflection coefficient r00 at %5ld. face "
                     ".........: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->r00[ii] ) = strtod( ptr, endp );
                  printf( " reflection coefficient r01 at %5ld. face "
                     ".........: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->r01[ii] ) = strtod( ptr, endp );
                  printf( " reflection coefficient r10 at %5ld. face "
                     ".........: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->r10[ii] ) = strtod( ptr, endp );
                  printf( " reflection coefficient r11 at %5ld. face "
                     ".........: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->r11[ii] ) = strtod( ptr, endp );
                  ii++;
               }; 

              periodic:

               printf( "\n\n number of periodic boundary faces "
                  "................: ");
               scanf( "%s", ptr );
               ( bpt->p ) = strtol( ptr, endp, DEC );

               if ( BNDPR < ( bpt->p ))
               {
                  printf( "\n\n Too many periodic boundary ports in "
                     "DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n Maximum number is %d = macro BNDPR "
                     "in %s.", BNDPR, __func__ );
                  printf( "\n [ Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in function scattr(*) "
                     "of program SOLVER.C.]\n" );
                  goto periodic;
               };

               if ( null < ( bpt->p ))
               {
                  printf( "\n propagation phase angle phi [= beta*P ] "
                     "...........: " );
                  scanf( "%s", ptr );
                  ( bpt->phi ) = strtod( ptr, endp );
               };

               ii = null; 
               while( ii < ( bpt->p ))
               {
                  printf( " \n");
                  printf( " index of %5ld. boundary cell ............."
                     "...........: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->cb[ii] ) = strtol( ptr, endp, DEC );
                  printf( " index of pertinent periodic cell .........."
                     "...........: "); 
                  scanf( "%s", ptr );
                  ( bpt->cp[ii] ) = strtol( ptr, endp, DEC );
                  printf( " 1st. port index of %5ld. boundary cell ...."
                     "...........: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->p0[ii] ) = strtol( ptr, endp, DEC );
                  printf(" ( signed ) index of pertinent periodic port "
                     ".........: ");
                  scanf( "%s", ptr );
                  ( bpt->pp0[ii] ) = strtol( ptr, endp, DEC ); 
                  printf( " 2nd. port index of %5ld. boundary cell ...."
                     "...........: ", ( ii+ONE ));
                  scanf( "%s", ptr );
                  ( bpt->p1[ii] ) = strtol( ptr, endp, DEC );
                  printf( " ( signed ) index of pertinent periodic port "
                     ".........: " );
                  scanf( "%s", ptr );
                  ( bpt->pp1[ii] ) = strtol( ptr, endp, DEC );
                  ii++;
               };
            }; /* end if( *ptr != 's' ): keyboard input */
/*............................................................................*/
            /* Completenes check for boundary faces: */

            if (( state->item ) < FORM_ITEMS )
            {
               printf( "\n Boundary completeness check"
                  "                      " );
               printf( "\n" );
               strcpy(( csp->rqfrm ), "brackets" );
               strcpy(( csp->rqstr ), "Search undefined adjacent boundary " );
               strcat(( csp->rqstr ), "faces ? [y/n]" );
               strcpy(( csp->dfstr ), "y" ); /* default */
/*............................................................................*/
               csp = txcnsl( csp );           /* enter on text console */
/*..........................................*/
               strcpy( ptr, ( csp->instr ));

               if (( *ptr != 'y' )\
                 &&( *ptr != 'Y' ))
                  goto disp3;
/*............................................................................*/
               tpt = linker( tpt );           /* link neighbouring cells      */
/*..........................................*/
            };
         }; /* end if (( state->stat ) != FORM_ITEMS ) */
/*............................................................................*/
# if FORM_BNDCMPL != 0 /* electric boundary completeness check: */
/*............................................................................*/
# if FORM_DISP == 1
         printf( "\n Boundary completeness check started           " );
         printf( "\n [ searching undefined adjacent boundary faces ]." );
# endif
/*............................................................................*/
         ii = null; nn = null;
         while( ii < ( bpt->n ))
         {
            mm = ( bpt->m[ii] ); /* boundary cell index */
            hh = ( ppt->mpt->idx[mm] ); /* its media index */

            if ( null != strncmp(( ppt->mpt->type[hh] ),
                 "trivial", FIVE )) /* non trivial electric cell */
            {
               fc = ( bpt->f[ii] );

# if DSC_ADJCLL == 0
               mn = ( tpt->mn[mm][fc] );
# elif DSC_ADJCLL == 1
               mn = ( tpt->mn[mm][(int)prt1[fc]] );

               if ( mn != ( tpt->mn[mm][(int)prt2[fc]] ))
               {
                  printf( "\n\n Error message from function %s:",
                     __func__ );
                  printf( "\n Inconsistent neighbouring mesh "
                     "cell definition on" );
                  printf( "\n cell no %ld, ports no %d and %d \n",
                     mm, ( prt1[fc]+ONE ), ( prt2[fc]+ONE ));

                  return ONE;
               };

# endif /* end if DSC_ADJCLL == 1 */

               if ( null < mn )
               {
                  hh = ( ppt->mpt->idx[mn] );
                  if ( null != strncmp(( ppt->mpt->type[hh] ),
                     "trivial", FIVE )) /* non trv electr neighbouring cell */
                  {
                     pp = abs( tpt->pn[mm][(int)prt1[fc]] ) - ONE;
                     qq = abs( tpt->pn[mm][(int)prt2[fc]] ) - ONE;

                     kk = null;
                     while( kk < ( bpt->n ))
                     {
                        if ( mn == ( bpt->m[kk] ))
                        {
                           if ( pp == prt1[( bpt->f[kk] )] )
                              goto next_boundary;
                           if ( pp == prt2[( bpt->f[kk] )] )
                              goto next_boundary;
                           if ( qq == prt1[( bpt->f[kk] )] )
                              goto next_boundary;
                           if ( qq == prt2[( bpt->f[kk] )] )
                              goto next_boundary;
                        };
                        kk++;
                     };
                     nn++;
/*...........................................................................*/
/* missing_adjacent/opposite_boundary: */

# if FORM_BNDCMPL == 1
                     printf( "\n\n Warning message from function %s:",
                        __func__ );
# elif FORM_BNDCMPL == 2
                     printf( "\n\n Error message from function %s:",
                        __func__ );
# endif
                     printf( "\n Undefined adjacent electric boundaries !!!" );
# if DSC_FCELBL == 0
                     printf( "\n Electric boundary on cell "
                        "%ld, face %d.", mm, fc );
                     printf( "\n adjacent to undefined boundary on cell "
                        "%ld, face %d.", mn, fce[pp] );
# else
                     printf( "\n Electric boundary on cell "
                        "%ld, face %d.", mm, fc );
                     printf( "\n adjacent to undefined boundary on cell "
                        "%ld, face %d.", mn, fce[pp]+ONE );
# endif

# if FORM_BNDCMPL == 2
                     return ONE;
# endif
                  }; /* end if *( ppt->mpt->type[hh] ) != 't'rivial [elect] ) */
               }; /* end if ( null < mn ) */
            }; /* end if ( null != strncmp(( ppt->mpt->type[hh] )... */
      next_boundary: ii++;
         }; /* while( ii < bpt->n ) */

# if FORM_DISP == 1
         printf( "\n\n Electric boundary check terminated."\
              "                                           " );
         if ( nn == null )
            printf( "\n No undefinded adjacent boundaries found." );
         else
         {
            PRBLDCLR( "\n\n " );
            printf( "%ld undefinded adjacent electric boundary faces !!!",
               ( long ) nn );
            PRNORMAL( "\n\n please acknowledge [enter any character]: " );
            scanf( "%s", ptr );
         };
	 printf( "\n" );
# endif /* FORM_DISP == 1 */
# endif /* FORM_BNDCMPL == 1 */
/*............................................................................*/
        disp3:

         ( state->bpt ) = bpt;
/*............................................................................*/
# if FORM_DISP == 1
         if ( null < ( bpt->n )) 
	 {
            printf( "\n %05ld aperiodic boundary faces: >-----------------"
               "--------------------------->\n", ( bpt->n ));
            llns = (( bpt->n ) / clns ) + ONE;

        /* bnd_cells1: */

            for ( kk=null; kk<llns ; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( bpt->m[( bpt->n )-ONE] ), bslsh );
                     goto bnd_faces;
                  };

                  if ( pp < ( bpt->n ))   
                     printf( "%6ld%c", ( bpt->m[pp] ), bslsh );

                  if (( bpt->n ) <= ( ++pp ))  
                     goto bnd_faces;
               };

              bnd_faces:

               printf( "\n -> face:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6d%c\n", ( bpt->f[( bpt->n )-ONE] ), slsh );
                     goto listo1;
                  };

                  if ( pp < ( bpt->n ))   
                     printf( "%6d%c", ( bpt->f[pp] ), slsh );

                  if (( bpt->n ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto listo1;
                  };
               };
            };
            listo1: ;
         }; /* end if null < ( bpt->n ) */

         if ( null < ( bpt->p ))
         {
            printf( "\n phi = beta*P ................................."
               "........: %.12e", bpt->phi );
            printf( "\n %05ld periodic boundary faces: >--------------"
               "------------------------------->\n", ( bpt->p ));
            llns = (( bpt->p ) / clns ) + ONE;

        /* bnd_cells2: */

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( bpt->cb[( bpt->p )-ONE] ), bslsh );
                     goto bnd_periodic;
                  };

                  if ( pp < ( bpt->p ))   
                     printf( "%6ld%c", ( bpt->cb[pp] ), bslsh);

                  if (( bpt->p ) <= ( ++pp ))  
                     goto bnd_periodic;
               };

              bnd_periodic:

               printf( "\n -> ngbr:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6ld%c\n", ( bpt->cp[( bpt->p )-ONE] ), slsh );
                     goto listo2;
                  };

                  if ( pp < ( bpt->p ))   
                     printf( "%6ld%c", ( bpt->cp[pp] ), slsh );

                  if (( bpt->p ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto listo2;
                  };
               };
            };
           listo2: ;
         };
# else /* if FORM_DISP != 1 */
         printf( "\n" );
# endif /* FORM_DISP ... */
/*............................................................................*/
         bndryfle = fopen( fleptr, "w" );

         printf( "\n opened: Boundary file %s", fleptr );
/*............................................................................*/
/* system identifier [model name] and comment: */

         fprintf( bndryfle, spformat, ( bpt->name ));
         fprintf( bndryfle, spformat, ( bpt->text ));

         strcpy( ptr, "DSC-MODEL_BOUNDARY_FILE_" );
         strcat( ptr, fleptr );

         pp = LNELNGTH - strlen( ptr );
         qq = null; do
         {
            fprintf( bndryfle, "%c", 95 );
         } while(( ++qq ) < pp );
         fprintf( bndryfle, "%s\n", ptr );
/*............................................................................*/
         fprintf( bndryfle, spformat, "Maxwell_field_operation" );

         if ( *( ppt->domain ) == 't' )
            fprintf( bndryfle, spformat, "TIME_DOMAIN" );
         else /* if ( ppt->domain ) == 'f'requency_domain */
            fprintf( bndryfle, spformat, "FREQUENCY_DOMAIN" );

         fprintf( bndryfle, spformat, "Maxwell_field_boundaries" );
         fprintf( bndryfle, "%-27s   %-ld\n",  "aperiodic_boundary_faces:",
            ( long ) ( bpt->n ));
         fprintf( bndryfle, "%-27s   %-ld\n",  "periodic_boundary_faces:",
            ( long ) ( bpt->p ));
/*............................................................................*/
# if DSC_HCRMDE != 0
/*............................................................................*/
# if DSC_FLDMDE == 0
         fprintf( bndryfle, "\n%-27s\n", "Thermal_operation" );
# else
         fprintf( bndryfle, "\n%-27s\n", "Thermal-fluid_operation" );
# endif
/*............................................................................*/
         fprintf( bndryfle, "%-27s\n", "TIME_DOMAIN" );
         fprintf( bndryfle, "%-27s   %-s\n", "thermal_boundaries", "number" );

/* skin-effect heat source faces: */
         fprintf( bndryfle, "%-27s   %ld\n", "skin-effect_heat_sources:",
            ( long ) ( bpt->nsk ));
/* fixed [imposed] heat current density: */
         fprintf( bndryfle, "%-27s   %ld\n", "fixed_heat_current_density:",
            ( long ) ( bpt->nhc ));
/* environment-surface heat conductivity [with respect to a given envm.tmp.] */
         fprintf( bndryfle, "%-27s   %ld\n", "environmnt-surface_h-cndct:",
            ( long ) ( bpt->nsc ));
/* fixed temperature faces: */
         fprintf( bndryfle, "%-27s   %ld\n", "fixed_temperature_faces:",
            ( long ) ( bpt->ntf ));
/* fixed temperature nodes: */
         fprintf( bndryfle, "%-27s   %ld\n", "fixed_temperature_nodes:",
            ( long ) ( bpt->ntn ));
/*............................................................................*/
# if DSC_FLDMDE != 0
         fprintf( bndryfle, "\n%-27s   %s\n", "fluid_boundaries",
           "number" );

/* no-slip faces: */
         fprintf( bndryfle, "%-27s   %ld\n", "no-slip:",
            ( long ) ( bpt->nns ));
/* free_slip faces: */
         fprintf( bndryfle, "%-27s   %ld\n", "free_slip:",
            ( long ) ( bpt->nsl ));
/* inflow faces: */
         fprintf( bndryfle, "%-27s   %ld\n", "inflow:",
            ( long ) ( bpt->nif ));
/* outflow faces: */
         fprintf( bndryfle, "%-27s   %ld\n", "outlow:",
            ( long ) ( bpt->nof ));
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* store boundary nodes, faces & s-parameters: */

         fprintf( bndryfle, "\n%s\n\n", fclbls );
            /* string "...LABELS_etc....", e.g. */

         if ( *( ppt->domain ) == 't' ) /* time_domain [ real parameters ] */
         { 
            for ( ii=null; ii<( bpt->n ); ii++ )
            {
               fprintf( bndryfle, i_format, ( bpt->m[ii] ));
               fprintf( bndryfle, i_format, ( bpt->f[ii] ));
            
               fprintf( bndryfle, "%+.16e  ", ( bpt->r00[ii] ));
               fprintf( bndryfle, "%+.16e\n", ( bpt->r01[ii] ));
               fprintf( bndryfle, "%+.16e  ", ( bpt->r10[ii] ));
               fprintf( bndryfle, "%+.16e\n", ( bpt->r11[ii] ));
            };
         }
         else /* if *( ppt->domain ) == 'f'req-domain [complex parameters] */
         {
            for ( ii=null; ii<( bpt->n ); ii++ )
            {
               fprintf( bndryfle, i_format, ( bpt->m[ii] ));
               fprintf( bndryfle, i_format, ( bpt->f[ii] ));

               fprintf( bndryfle, "%+.16e  ", ( bpt->r00[ii] ));
               fprintf( bndryfle, "%+.16e\n", ( bpt->r01[ii] ));
               fprintf( bndryfle, "%+.16e  ", ( bpt->i00[ii] ));
               fprintf( bndryfle, "%+.16e\n", ( bpt->i01[ii] ));
               fprintf( bndryfle, "%+.16e  ", ( bpt->r10[ii] ));
               fprintf( bndryfle, "%+.16e\n", ( bpt->r11[ii] ));
               fprintf( bndryfle, "%+.16e  ", ( bpt->i10[ii] ));
               fprintf( bndryfle, "%+.16e\n", ( bpt->i11[ii] ));
            };
         }; /* end if *( ppt->domain ) == 'f'requency_domain */

         if ( null < ( bpt->p ))
            fprintf( bndryfle, d_format, ( bpt->phi ));

         for ( ii=null; ii<( bpt->p ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->cb[ii] ));
            fprintf( bndryfle, i_format, ( bpt->cp[ii] ));
            fprintf( bndryfle, i_format, ( bpt->p0[ii] ));
            fprintf( bndryfle, i_format, ( bpt->pp0[ii] ));
            fprintf( bndryfle, i_format, ( bpt->p1[ii] ));
            fprintf( bndryfle, i_format, ( bpt->pp1[ii] ));
         };
/*............................................................................*/
# if DSC_HCRMDE != 0
         for( ii=null; ii<( bpt->nsk ); ii++ )
         {
/*............................................................................*/
/* heat sources opposite to [ skin effect ] lossy boundaries: */

            mm = ( bpt->msk[ii] );
            hh = ( ppt->mpt->idx[mm] ); /* media index mm */

/* check if non-trivial thermal cell: */
            if (( null < hh )\
              &&( NULL == strstr(( ppt->mpt->type[hh] ), "trv_hcr" )))
            {
               fc = ( bpt->fsk[ii] );

               if (( NULL == strstr(( ppt->mpt->type[hh] ),
                  "trivial" ))) /* cell mm non-trivial electric */
               {
                  mn = mm;
                  fn = fc;
               }
               else /* trivial electric cell mm [ yet check neighb cell ] */
               {
/*............................................................................*/
# if DSC_ADJCLL == 0
                  mn = ( tpt->mn[mm][fc] );
# elif DSC_ADJCLL == 1
                  mn = ( tpt->mn[mm][(int)prt1[fc]] );

                  if ( mn != ( tpt->mn[mm][(int)prt2[fc]] ))
                  {
                     printf( "\n\n Error message from function %s:",
                        __func__ );
                     printf( "\n Inconsistent neighbouring mesh "
                        "cell definition on" );
                     printf( "\n cell %ld, ports %d and %d !\n",
                        mm, ( prt1[fc]+ONE ), ( prt2[fc]+ONE ));
                     return ONE;
                  };
# endif /* end if DSC_ADJCLL == 1 */
/*............................................................................*/
                  if ( null < mn ) /* existing neigboring mesh cell */
                  {
                     hh = ( ppt->mpt->idx[mn] ); /* [ media index ] */

                     if (( hh != null ) /* non-trivial electric [nghb] cell */\
                       &&( null != strncmp(( ppt->mpt->type[hh] ),
                          "trivial", FIVE )))
                     {
/*............................................................................*/
# if DSC_FCELBL == 1
                        fn = ( tpt->fn[mm][fc] ) - ONE;
# elif DSC_FCELBL == 2
                        fn = ( tpt->fn[mm][fc] );

                        if ( fn < null )
                           fn = -fn;

                        fn -= ONE;
# else /* DSC_FCELBL == 0, e.g. */
                        fn = ( tpt->fn[mm][fc] );
# endif /* DSC_FCELBL == ... */
/*............................................................................*/
                     }; /* end if non-trivial neighb cell */
                  }; /* end if ( null < mn ) */
               }; /* end if cell mm trivial electric */
/*............................................................................*/
# if BND_CSHAPE == 1
/*   transfer vertex point coordinates:
*//* [ if ( clp->cell ) != 0 is given, 
*//*   then the coordinate transfer is done in function cshape(*) ]
*//*
               kk = null; do
               {
                  ll = null; do
                  {
                     ( clp->c[kk][ll] ) = \
                     ( ppt->cpt->c[( tpt->cm[mn][kk] )][ll] );
                  } while(( ++ll ) < DIMNS );
               } while(( ++kk ) < CRNRS );
*/
/*............................................................................*/
/* call cell shape function, option 3, cell mn, face fn: */

               ( clp->opt ) = 'f';
               ( clp->cell ) = mn;
               ( clp->face ) = fn;
/*............................................................................*/
               clp = cshape( clp );
/*..................................................*//* for given time step */
               if (( clp->rtn ) == ONE )
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Error on calling function cshape(*), "
                     "option 3, cell number %ld !!!\n", mn );

                  ( csp->dfopt ) = null; /* next default menu option */
                  goto menu;
               };

               ss = ( clp->fm[fn] )*( bpt->rs[ii] );
               ss = admce*sqrt( ss );

               kk = null; do
               {
                  ll = null; do
                  {
                     ( bpt->sk[ii][kk][ll] ) = ss* \
                        ( double )( TWO*(( ll+fn )%TWO ) - ONE )* \
                           ( clp->vu[kk][( short )(( ll+ONE )%TWO )] );
                  } while(( ++ll ) < TWO );
               } while(( ++kk ) < TWO );
# endif /* BND_CSHAPE == 1 */
/*............................................................................*/
/* check for undefined pertinent [skin effect] lossy boundary: */

               kk = null;
               while( kk < ( bpt->n ))
               {
                  nn = ( bpt->m[kk] );
                  ff = ( bpt->f[kk] );

                  if (( nn == mn )\
                    &&( ff == fn ))
                     goto wrt_bndhsc;
                  kk++;
               };
/*...........................................................................*/
/* missing [ skin effect lossy ] boundary: */

               printf( "\n\n Error message from function %s:", __func__ );
               printf( "\n Boundary definition incomplete  !!!" );
               printf( "\n Missing skin effect lossy boundary on cell "
                  "%ld, face %d.", mn, fn );
               printf( "\n [ pertinent to defined thermal boundary on cell "
                  "%ld, face %d ].\n", mm, fc );

               return ONE;
/*............................................................................*/
            }; /* end if cell mm thermal non-trivial */

            kk = null; do
            {
               ll = null; do
               {
                  ( bpt->sk[ii][kk][ll] ) = ZERO;
               } while(( ++ll ) < TWO );
            } while(( ++kk ) < TWO );
/*............................................................................*/
           wrt_bndhsc:
       
            fprintf( bndryfle, i_format, ( bpt->msk[ii] ));
            fprintf( bndryfle, i_format, ( bpt->fsk[ii] ));
            fprintf( bndryfle, d_format, ( bpt->rs[ii] ));

            kk = null; do
            {
               ll = null; do
               {
                  fprintf( bndryfle, "%+.16e  ", ( bpt->sk[ii][kk][ll] ));
               } while(( ++ll ) < TWO );
               fprintf( bndryfle, "\n" );
            } while(( ++kk ) < TWO );
         }; /* next ii */

/* end of heat source boundary parameter determination */
/*............................................................................*/
/* [ imposed ] heat currents incident on the boundary: */

         for( ii=null; ii<( bpt->nhc ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->mhc[ii] ));
            fprintf( bndryfle, i_format, ( bpt->fhc[ii] ));
            fprintf( bndryfle, d_format, ( bpt->hc[ii] ));
            fprintf( bndryfle, d_format, ( bpt->rd[ii] ));
         };

/* surface heat conducting faces [ approximating convection, e.g.]: */

         for( ii=null; ii<( bpt->nsc ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->msc[ii] ));
            fprintf( bndryfle, i_format, ( bpt->fsc[ii] ));
            fprintf( bndryfle, d_format, ( bpt->sc[ii] )); /* hcd [W/(K*m^2)] */
            fprintf( bndryfle, d_format, ( bpt->tr[ii] )); /* refer. temp [C] */
         };

/* fixed temperature faces: */

         for( ii=null; ii<( bpt->ntf ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->mtf[ii] ));
            fprintf( bndryfle, i_format, ( bpt->ftf[ii] ));
            fprintf( bndryfle, d_format, ( bpt->tf[ii] )); /* face temp [C] */
         };

/* fixed temperature nodes: */

         for( ii=null; ii<( bpt->ntn ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->mtn[ii] ));
            fprintf( bndryfle, d_format, ( bpt->tn[ii] )); /* node temp [C] */
         };
/*............................................................................*/
# if DSC_FLDMDE != 0 
/* no-slip boundary faces: */

         for( ii=null; ii<( bpt->nns ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->mns[ii] )); /* cell index */
            fprintf( bndryfle, i_format, ( bpt->fns[ii] )); /* face index */
/*............................................................................*/
# if NUSSELT != 0
            fprintf( bndryfle, d_format, ( bpt->nus[ii] )); /* Nusselt number */
# endif
/*............................................................................*/
         };

/* free slip boundary faces: */

         for( ii=null; ii<( bpt->nsl ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->msl[ii] ));
            fprintf( bndryfle, i_format, ( bpt->fsl[ii] ));

            kk = null; do
            {
               fprintf( bndryfle, "%+.16e  ", ( bpt->nf[ii][kk] ));
            } while(( ++kk ) < THREE );
            fprintf( bndryfle, "\n" );
         };

/* inflow boundary faces: */

         for( ii=null; ii<( bpt->nif ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->mif[ii] ));
            fprintf( bndryfle, i_format, ( bpt->fif[ii] ));

/* fluid temperature at inflow: */

            fprintf( bndryfle, d_format, ( bpt->ti[ii] ));

/* inflow fluid velocity: */

            kk = null; do
            {
               fprintf( bndryfle, "%+.16e  ", ( bpt->uf[ii][kk] )); /* [m/s] */
            } while(( ++kk ) < THREE );
            fprintf( bndryfle, "\n" );
         };

/* outflow boundary faces: */

         for( ii=null; ii<( bpt->nof ); ii++ )
         {
            fprintf( bndryfle, i_format, ( bpt->mof[ii] ));
            fprintf( bndryfle, i_format, ( bpt->fof[ii] ));

            kk = null; do
            {
               fprintf( bndryfle, "%+.16e  ", ( bpt->vf[ii][kk] )); /* [m/s] */
            } while(( ++kk ) < THREE );
            fprintf( bndryfle, "\n" );
         };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* close boundary file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( bndryfle, "\nDSC model boundary file %s ", fleptr );
         fprintf( bndryfle, "created:\n%s\n", tmestr );

         fclose( bndryfle );

         printf( CLEAR_LINE );
         printf( "\r Boundary file %s ", fleptr );
         printf( time_frm, tmestr );

         if (( state->item ) == FORM_ITEMS )
            printf( "\n -----------------------------------"
               "-------------------------------------------");
         break;
/*............................................................................*/
        case 4: /* >---------- excitation file generation -----------------> */

         if (( state->item ) == FORM_ITEMS )
            printf( "\n\n Excitation file generation: " );

         strcpy( fleptr, prefix );
         strcat( fleptr, excptr );

         if (( state->stat ) == FORM_ITEMS )
         {
            strcpy(( ept->name ), ( tpt->name ));
            strcat( fleptr, ( state->flbl ));
            printf( "\n\n Comment: %.70s \n", ( ept->text ));
         }
         else if (( state->stat ) != FORM_ITEMS )
         {
            if (( state->item ) < FORM_ITEMS )
            {
               printf( "\n" );

               strcpy(( csp->rqlng ),
                  "Please enter index <N> of excitation file '" );
               strcat(( csp->rqlng ), fleptr );
               strcat(( csp->rqlng ), "<N>'" );
               strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
               csp = txcnsl( csp );            /* enter index N on console    */
/*...........................................*/
               strcpy( ptr, ( lotos(( csp->inlng ), null )));
               strncpy(( state->flbl ), ptr, THREE );
            }
            else 
               printf( "\n" );

            strcat( fleptr, ( state->flbl ));

            strcpy(( csp->rqfrm ), "brackets" );
            strcpy(( csp->rqstr ), "Input on keyboard or subroutine " );
            strcat(( csp->rqstr ), "sysbnd(*) ? [k/s]" );
            strcpy(( csp->dfstr ), "s" ); /* default */
/*............................................................................*/
            csp = txcnsl( csp );      /* enter string on text console         */
/*..................................*/
            strcpy( ptr, ( csp->instr ));
          
            if (( *ptr == 's' )||( *ptr == 'S' ))
            {
/*............................................................................*/
# if FORM_DISP == 1
               printf( "\n DSC model excitation functions "
                  "started.                        " );
               PRBLDCLR( "\r" );
               printf( "\n %*s", 78, "DSC-FORMER" );
               PRNORMAL( "" );
# endif
/*............................................................................*/
               if (( state->item ) < FORM_ITEMS )
               {
/*............................................................................*/
                  tpt = systop( state );            /* transfer: check DSC -  */
/*................................................*//* mesh topology          */

                  switch( tpt->rtn )
                  {
                    case -ONE:
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu; 

                    case ONE:
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in function %s !!!\n", "systop(*)" );
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;

                    case FOUR: /* = ( state->item ) [ excitation ] */
                     ( state->tpt ) = tpt;
                     ( state->stat ) = ( state->item );
/*............................................................................*/
# if FORM_DISP == 1
                     printf( "\n ==================================="
                        "===========================================" );
                     PRBLDCLR( "\r" );
                     printf( "\n %*s", 78, "DSC-FORMER" );
                     PRNORMAL( "\r" );
                     printf( "\n DSC model excitation function "\
                        "%s terminated.                 \n", "sysexc(*)" );
# endif
/*............................................................................*/
                     goto identf_check4;

                    default: 
                     ( state->tpt ) = tpt;
                     break;
                  }; 
/*............................................................................*/
                  ppt = syssmx( state );             /* DSC mesh coordinates, */
/*.................................................*//* media etc.            */

                  switch( ppt->rtn )
                  {
                    case -ONE:
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;

                    case ONE:
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in function %s !!!\n", "syssmx(*)" );
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;
 
                    default:
                     ( state->ppt ) = ppt;
                     break;
                  };
               };/* end if (( state->item ) < FORM_ITEMS ) */

               ( ept->rtn ) = -TWO;
/*............................................................................*/
               ept = sysexc( state );            /* DSC mesh excitation file  */
/*.............................................*//* configuration             */

               switch( ept->rtn )
               {
                 case -TWO:
                  printf( "\n Message from function %s:", __func__ );
                  printf( "\n Unsatisfied code in function %s !!!\n",
                     "sysexc(*)" );
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 case -ONE:
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 case ONE:
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Error in function %s !!!\n", "sysexc(*)" );
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 break;

                 default:
/*............................................................................*/
# if FORM_DISP == 1
                  printf( "\n ==================================="
                     "===========================================" );
                  PRBLDCLR( "\r" );
                  printf( "\n %*s", 78, "DSC-FORMER" );
                  PRNORMAL( "\r" );
                  printf( "\n DSC model excitation function "\
                     "%s terminated.                 \n", "sysexc(*)" );
# endif
/*............................................................................*/
                  break;
               };
              
              identf_check4:

               if ( null != strncmp( tpt->name, ept->name, FOUR ) )
               {
                  printf( "\n Message from function %s:", __func__ );
                  printf( "\n Incompatible DSC system identifiers "
                     "on functions %s and %s.", "systop(*)", "syssmx(*)" );
                  printf( "\n Please check functions and restart program ! ");
                  return ONE; 
               };

               if (( EXCEP < ( ept->ne ) )||( EXCHP < ( ept->nh )))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  if ( EXCEP < ( ept->ne ) )
                  {
                     printf( "\n Too many excited E-ports "
                        "in DSC mesh %s !!!", ( tpt->name ));
                     printf( "\n [ Maximum number is %d = macro EXCEP "
                        "in %s.", EXCEP, "FORMER.CONF" );
                  };
                  if ( EXCHP < ( ept->nh ))
                  {
                     printf( "\n Too many excited H-ports "
                        "in DSC mesh %s !!!", ( tpt->name ));
                     printf( "\n [ Maximum number is %d = macro EXCHP "
                        "in %s.", EXCHP, "FORMER.CONF" );
                  };
                  printf( "\n - Change macro only in compliance "
                     "with memory resources" );
                  printf( "\n   and same macro in function excite(*) "
                     "of program SOLVER.C.]\n" );
                  return ONE;
               };
/*............................................................................*/
# if DSC_HCRMDE != 0
               if ( EXCHC < ( ept->nhc ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many heat current faces to be excited "
                     "in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro EXCHC "
                     "in %s.", EXCHC, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in header 'SOLVER.CONF' "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( EXCTF < ( ept->ntf ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many temperature faces to be excited "
                     "in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro EXCTF "
                     "in %s.", EXCTF, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in header 'SOLVER.CONF' "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( EXCTN < ( ept->ntn ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many temperature nodes to be excited "
                     "in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro EXCTN "
                     "in %s.", EXCTN, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in header 'SOLVER.CONF' "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };
# endif
/*............................................................................*/
            }/* end if( *ptr == 's'ubroutine ) */
            else /* keyboard input */
            {
               printf( "\n\n Please enter: --------"
                  "-------------------------------->" );
               printf( "\n -> DSC system identifier [ structure name ] "
                  "..........: " );
               scanf( "%s", ( ept->name ));
               printf( " -> Excitation type [ e.g. 'DIRAC_PULSE' ... ] "
                  "........: " );  
               scanf( "%s", ( ept->type ));

              types1: /* no more than 19 characters !!! */

               if ( null == strncmp (( ept->type ),
                  "ZERO_EM_FIELD", TWO ))
               {
                  strcpy(( ept->type ), "ZERO_EM_FIELD______" );
                  extyp = 0;
               }
               else if (( *( ppt->domain ) == 'f' )||
                  ( null == strncmp (( ept->type ),
                  "STEADY_STATE", TWO )))
               {
                  strcpy(( ept->type ), "STEADY_STATE_______" );
                  extyp = 1;
               }
               else if ( null == strncmp(( ept->type ),
                  "DIRAC_PULSE", TWO ) )
               {
                  strcpy(( ept->type ), "DIRAC_PULSE________" );
                  extyp = 2;
               }
               else if ( null == strncmp(( ept->type ),
                  "HARMONIC_SINUSOIDAL", TWO ) )
               {
                  strcpy(( ept->type ), "HARMONIC_SINUSOIDAL" );
                  extyp = 3;
                  printf( "\n Please enter exciting frequency [ Hz ] "
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->fr[null] ) = strtod( ptr, endp );
               }
               else if ( null == strncmp(( ept->type ),
                  "SMOOTH_HARMONIC", TWO ) )
               {
                  strcpy(( ept->type ), "SMOOTH_HARMONIC____" );
                  extyp = 4;
                  printf( "\n Please enter exciting frequency [ Hz ] "
                     "...............: ");
                  scanf( "%s", ptr );
                  ( ept->fr[null] ) = strtod( ptr, endp );
                  printf( "\n Enter transition time [ seconds ] ....."
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->dt ) = strtod( ptr, endp );
               }
               else if ( null == strncmp(( ept->type ),
                  "MULTIPLE_HARMONIC", TWO ) )
               {
                  strcpy(( ept->type ), "MULTIPLE_HARMONIC__" );
                  extyp = 5;
                  printf( "\n Enter exciting frequencies "
                     " ( Escape: enter ZERO ) >--->\n" );
                  ( ept->nn ) = null; do
                  {
                     printf( " Enter %2d. frequency [ Hz ] "
                        "...........................: ", (( ept->nn ) + ONE ));
                     scanf( "%s", ptr );
                     ( ept->fr[( ept->nn )] ) = strtod( ptr, endp );
                     ( ept->nn )++;
                  }  while(((ept->nn)<EXCFR)\
                         &&((ept->fr[(ept->nn)-ONE])!=ZERO));
                  ( ept->nn ) -= ONE;
               }
               else if ( null == strncmp(( ept->type ),
                  "GAUSS_PULSE", TWO ) )
               {
                  strcpy(( ept->type ), "GAUSS_PULSE________" );
                  extyp = 6;
                  printf( "\n Please enter sigma time [ seconds ] "
                     "..................: " );
                  scanf( "%s", ptr );
                  ( ept->rt ) = strtod( ptr, endp );
                  printf( "\n Enter delay time [ seconds ] "
                     ".........................: " );
                  scanf( "%s", ptr );
                  ( ept->dt ) = strtod( ptr, endp );
               }
               else if ( null == strncmp(( ept->type ),
                  "OSCILLATING_GAUSS", TWO ) )
               {
                  strcpy(( ept->type ), "OSCILLATING_GAUSS__" );
                  extyp = 7;
                  printf( "\n Please enter exciting frequency [ Hz ] "
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->fr[null] ) = strtod( ptr, endp );
                  printf( "\n Enter sigma time [ seconds ] .........."
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->rt ) = strtod( ptr, endp );
                  printf( "\n Enter delay time [ seconds ] .........."
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->dt ) = strtod( ptr, endp );
               }
               else if ( null == strncmp(( ept->type ),
                  "RAMP_PULSE", TWO ) )
               {
                  strcpy(( ept->type ), "RAMP_PULSE_________" );
                  extyp = 8;
                  printf( "\n Please enter rise time [ seconds ] "
                     "...................: " );
                  scanf( "%s", ptr );
                  ( ept->rt ) = strtod( ptr, endp );
               }
               else if ( null == strncmp(( ept->type ),
                  "RECTANGULAR_PERIODIC", TWO ) )
               {
                  strcpy(( ept->type ), "RECTANGULAR_PERIODC" );
                  extyp = 9;
                  printf( "\n Please enter exciting frequency [ Hz ] "
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->fr[null] ) = strtod( ptr, endp );
               }
               else if ( null == strncmp(( ept->type ),
                  "SAW_TOOTH_PERIODIC", TWO ) )
               {
                  strcpy( ept->type, "SAW_TOOTH_PERIODIC_" );
                  extyp = 10;
                  printf( "\n Please enter exciting frequency [ Hz ] "
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->fr[null] ) = strtod( ptr, endp );
               } 
               else if ( null == strncmp(( ept->type ),
                  "HEAVISIDE_STEP", TWO ) )
               {
                  strcpy(( ept->type ), "HEAVISIDE_STEP_____" );
                  extyp = 11;
               }
               else if ( null == strncmp(( ept->type ),
                  "DSC_TIMESTEP_PRDC", TWO ) )
               {
                  strcpy(( ept->type ), "DSC_TIMESTEP_PRDC__" ); 
                  extyp = 12;   
               }
               else if ( null == strncmp(( ept->type ),
                  "NOISE_AT_RANDOM", TWO ) )
               {
                  strcpy(( ept->type ), "NOISE_AT_RANDOM____" );
                  extyp = 13;
               }
               else if ( null == strncmp(( ept->type ),
                  "MORLET_WAVELET", TWO ) )
               {
                  strcpy(( ept->type ), "MORLET_WAVELET_____" );
                  extyp = 14;
                  printf( "\n Please enter exciting frequency [ Hz ] "
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->fr[null] ) = strtod( ptr, endp );
                  printf( "\n Enter sigma time [ seconds ] "
                     ".........................: " );
                  scanf( "%s", ptr );
                  ( ept->rt ) = strtod( ptr, endp );
                  printf( "\n Enter delay time [ seconds ] "
                     ".........................: " );
                  scanf( "%s", ptr );
                  ( ept->dt ) = strtod( ptr, endp );
               }
               else if ( null == strncmp(( ept->type ),
                  "WAVE_PACKET_HARMNC", TWO ) )
               {
                  strcpy(( ept->type ), "WAVE_PACKET_HARMNC_" );
                  extyp = 15;
                  printf( "\n Please enter exciting frequency [ Hz ] "
                     "...............: " );
                  scanf( "%s", ptr );
                  ( ept->fr[null] ) = strtod( ptr, endp );
                  printf( "\n Enter transition time  [ seconds ] "
                     "..................: " );
                  scanf( "%s", ptr );
                  ( ept->dt ) = strtod( ptr, endp );
                  printf( "\n Enter plateau time  [ seconds ] "
                     ".....................: " );
                  scanf( "%s", ptr );
                  ( ept->ht ) = strtod( ptr, endp );
               }
               else
               {
                  printf( "\n\n Unknown excitation type !!! " );
                  printf( "\n\n Select one of the following types: " );  
                  printf( "\n\n ZERO-EM-FIELD______ " );
                  printf( "\n STEADY_STATE_______ " );
                  printf( "\n DIRAC_PULSE________ " );
                  printf( "\n HARMONIC_SINUSOIDAL " );  
                  printf( "\n SMOOTH_HARMONIC____ " );
                  printf( "\n MULTIPLE_HARMONIC__ " );
	          printf( "\n GAUSS_PULSE________ " );
                  printf( "\n OSCILLATORY_GAUSS__ " );
                  printf( "\n RAMP_PULSE_________ " );
                  printf( "\n RECTANGULAR_PERIODC " );
                  printf( "\n SAW_TOOTH_PERIODIC_ " );
                  printf( "\n HEAVISIDE_STEP_____ " ); 
                  printf( "\n DSC_TIMESTEP_PRDC__ " );
                  printf( "\n NOISE_AT_RANDOM____ " );
                  printf( "\n MORLET_WAVELET_____ " );
                  printf( "\n WAVE_PACKET_HARMNC_ " );
                  printf( "\n ZERO-EM-FIELD______ " );
                  printf( "\n\n Please enter excitation type [ at "
                     "least two leading characters ]:" );
                  printf( "\n -> Enter type [ e.g 'DI'RAC_PULSE ] "
                     "..................: " );
                  scanf( "%s", ept->type );

                  goto types1;
               }; 

               printf( "\n Please enter: ------------------------"
                  "---------------->" );
               printf( "\n -> Comment ..........................."
                  "................: " );
               scanf( "%s", ept->text );

              e_input1:

               printf( "\n E-field excitation: " );
               printf( "\n\n Please enter cell and port indices and "
                  "voltages >-----> " );
               printf( "\n ( Escape: enter ZERO ) >----------------"
                  "-------------->\n" );

               ( ept->ne ) = null; do
               {
                  printf( " mesh cell index [ %5d. port ] ?......."
                     "..............: ", (( ept->ne ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto nep;

                  ( ept->me[( ept->ne )] ) = strtol( ptr, endp, DEC );
                  printf( " port index      [ %5d. port ] ?......."
                     "..............: ", (( ept->ne ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto nep;

                  ( ept->pe[( ept->ne )]) = strtol( ptr, endp, DEC );
                  printf( " voltage/V       [ %5d. port ] ?......."
                     "..............: ", (( ept->ne ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto nep;

                  ( ept->er[( ept->ne )] ) = strtod( ptr, endp );
                  ( ept->ne )++;
               } while((( ept->ne ) < EXCEP )\
                      &&( null < ( ept->me[( ept->ne ) - ONE] )));
               ( ept->ne )-- ;

              nep:
               printf( "\n Number of excited <E|p> ports "
                  "........................: %-6d ", ( ept->ne ));

           /* h_input1: */

               printf( "\n\n H-field excitation: " );
               printf( "\n\n Please enter cell and port indices "
                  "and currents >-----> " );
               printf( "\n [ Escape: enter ZERO ] >------------"
                  "------------------> \n" );

               ( ept->nh ) = null; do
               {
                  printf( " mesh cell index [ %5d. port ] ?..."
                     "..................: ", (( ept->nh ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto nhp;

                  ( ept->mh[( ept->nh )] ) = strtol( ptr, endp, DEC );
                  printf( " port index      [ %5d. port ] ?..."
                     "..................: ", (( ept->nh ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto nhp;

                  ( ept->ph[( ept->nh )] ) = strtol( ptr, endp, DEC );
                  printf( " current/A       [ %5d. port ] ?..."
                     "..................: ", (( ept->nh ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto nhp;

                  ( ept->hr[( ept->nh )] ) = strtod( ptr, endp );
                  ( ept->nh )++ ;
               } while((( ept->nh )< EXCHP )\
                      &&( null < ( ept->mh[( ept->nh ) - ONE] )));
               ( ept->nh )-- ;

              nhp:

               printf( "\n Number of excited <H|p> ports "
                  "........................: %-6d ", ( ept->nh ));

               printf( "\n\n                                  "
                  "                                     ) <- ?" );
               printf( "\r Input correct ? >--------------"
                  "----------------------> [ y/n ] >--> (" );
               scanf( "%s", ptr );
               if (( *ptr == 'n' )||( *ptr == 'N' ))
                  goto e_input1;
            }; /* end if( *ptr != 's'): keyboard input */

            ( state->ept ) = ept;

         }; /* end if (( state->stat ) != FORM_ITEMS ) */
/*............................................................................*/
/* Neither time nor frequency domain Maxwell field excitation [ZERO EM field] */
/* - Only thermal&fluid computations, e.g. */

         if ( null == strncmp (( ept->type ), "ZE", TWO ))
         {
            strcpy (( ept->type ), "ZERO_EM_FIELD______" );
            extyp = 0;
            goto disp4;
         }
/*............................................................................*/
/* frequency domain Maxwell field excitation types [steady state]: */

         else if (( *( ppt->domain ) == 'f' )||
             ( null == strncmp (( ept->type ), "ST", TWO )))
         {
            strcpy (( ept->type ), "STEADY_STATE_______" );
            extyp = 1;
            goto disp4;
         }
/*............................................................................*/
/* time domain Maxwell field excitation types: */

         else if ( null == strncmp (( ept->type ),
            "DIRAC_PULSE" , TWO ))
         {
            strcpy (( ept->type ), "DIRAC_PULSE________" );
            extyp = 2;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "HARMONIC_SINUSOIDAL", TWO ))
         {
            strcpy (( ept->type ), "HARMONIC_SINUSOIDAL" );
            extyp = 3;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "SMOOTH_HARMONIC", TWO ))
         {
            strcpy (( ept->type ), "SMOOTH_HARMONIC____" );
            extyp = 4;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "MULTIPLE_HARMONIC", TWO ))
         {
            strcpy (( ept->type ), "MULTIPLE_HARMONIC__" );
            extyp = 5;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "GAUSS_PULSE", TWO ))
	 {
            strcpy (( ept->type ), "GAUSS_PULSE________" );
	    extyp = 6;
	    goto disp4;
	 }
         else if ( null == strncmp(( ept->type ),
            "OSCILLATING_GAUSS", TWO ))
         {
            strcpy (( ept->type ), "OSCILLATING_GAUSS__" );
            extyp = 7;
            goto disp4;
         } 
         else if ( null == strncmp(( ept->type ),
            "RAMP_PULSE", TWO ))
         {
            strcpy (( ept->type ), "RAMP_PULSE_________" );
            extyp = 8;
            goto disp4;
         } 
         else if ( null == strncmp(( ept->type ),
            "RECTANGULAR_PERIODIC", TWO ))
         {
            strcpy (( ept->type ), "RECTANGULAR_PERIODC" );
            extyp = 9;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "SAW_TOOTH_PERIODIC", TWO ))
         {
            strcpy (( ept->type ), "SAW_TOOTH_PERIODIC_" );
            extyp = 10;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "HEAVISIDE_STEP", TWO ))
         {
            strcpy (( ept->type ), "HEAVISIDE_STEP_____" );
            extyp = 11;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "DSC_TIMESTEP_PRDC", TWO ))
         {
            strcpy (( ept->type ), "DSC_TIMESTEP_PRDC__" );
            extyp = 12;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "NOISE_AT_RANDOM", TWO ))
         {
            strcpy (( ept->type ), "NOISE_AT_RANDOM____" );
            extyp = 13;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "MORLET_WAVELET", TWO ))
         {
            strcpy (( ept->type ), "MORLET_WAVELET_____" );
            extyp = 14;
            goto disp4;
         }
         else if ( null == strncmp(( ept->type ),
            "WAVE_PACKET_HARMNC", TWO ))
         {
            strcpy (( ept->type ), "WAVE_PACKET_HARMNC_" );
            extyp = 15;
            goto disp4;
         }
         else 
         {
            printf( "\n\n Message from function %s:", __func__ );
            printf( "\n Unknown excitation type %s !!! ", ( ept->type ));
            return ONE; 
         };

        disp4: ;
/*............................................................................*/
# if DSC_HCRMDE != 0 
/* thermal excitation types: */

         if ( null == strncmp (( ept->hctp ),
            "PASSIVE" , TWO ))
         {
            strcpy (( ept->hctp ), "PASSIVE____________" );
            exctp = null;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "FLOATING" , TWO ))
         {
            strcpy (( ept->hctp ), "FLOATING___________" );
            exctp = ONE;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "STEADY_STATE" , TWO ))
         {
            strcpy (( ept->hctp ), "STEADY_STATE_______" );
            exctp = 2;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "DIRAC_PULSE" , TWO ))
         {
            strcpy (( ept->hctp ), "DIRAC_PULSE________" );
            exctp = 3;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "SMOOTH_PULSE" , TWO ))
         {
            strcpy (( ept->hctp ), "SMOOTH_PULSE_______" );
            exctp = 4;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "GAUSS_PULSE" , TWO ))
         {
            strcpy (( ept->hctp ), "GAUSS_PULSE________" );
            exctp = 5;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "RAMP_PULSE" , TWO ))
         {
            strcpy (( ept->hctp ), "RAMP_PULSE_________" );
            exctp = 6;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "HEAVISIDE_STEP" , TWO ))
         {
            strcpy (( ept->hctp ), "HEAVISIDE_STEP_____" );
            exctp = 7;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "PLATEAU_SMOOTHED" , TWO ))
         {
            strcpy (( ept->hctp ), "PLATEAU_SMOOTHED___" );
            exctp = 8;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "PERIODIC_PULSE" , TWO ))
         {
            strcpy (( ept->hctp ), "PERIODIC_PULSE_____" );
            exctp = 9;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "DOUBLE_PERIODIC" , TWO ))
         {
            strcpy (( ept->hctp ), "DOUBLE_PERIODIC____" );
            exctp = 10;
            goto disp5;
         }
         else if ( null == strncmp (( ept->hctp ),
            "SWITCHED_SOURCES", TWO ))
         {
            strcpy (( ept->hctp ), "SWITCHED_SOURCES___" );
            exctp = 11;
            goto disp5;
         }
         else
         {
            strcpy (( ept->hctp ), "PASSIVE____________" );
            exctp = null;
            goto disp5;
         };

        disp5: ;

# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
# if FORM_DISP == 1
         if (( state->item ) < FORM_ITEMS )
         {
            printf( "\n DSC system identifier [ structure name ]: "
               "%.22s", ( ept->name ));
            printf( "\n Comment: %.70s\n", ( ept->text ));
         };
/*............................................................................*/
         printf( "\n Maxwell field excitation type: %s", ( ept->type ));
        
         switch ( extyp )
         {
           default:
            break;

           case 0:
            printf( "\n ZERO Maxwell field exciation" );
            ( ept->ne ) = null;
            ( ept->nh ) = null;
            ( ept->nn ) = null;
            ( ept->dt ) = ZERO;
            ( ept->ht ) = ZERO;
            ( ept->rt ) = ZERO;
            ( ept->fr[null] ) = ZERO;
            break;

           case 1:

            if ( *( ppt->domain ) == 'f' )
            {
               ( ept->fr[null] ) = ( ppt->fr );
               printf( "\n Exciting frequency ................."
                  "..................: %.12e Hz", ( ppt->fr ));
            };

            if ( ZERO < ( ept->dt ))
            {
               if (( ept->nn ) == null )
                  ( ept->nn ) = TWO;

               printf( "\n Smoothing order ...................."
                  "..................: %d", ( ept->nn ));

               if ( *( ppt->domain ) == 'f' )
                  printf( "\n Number of smoothing cycles ........."
                     "..................: %ld", ( long ) ( ept->dt ));
               else
                  printf( "\n Smoothing time ....................."
                     "..................: %.12e s", ( ept->dt ));
            }
            else
               ( ept->nn ) = null;

            break;

           case 2:
            break;
	    
           case 3:
           case 9:
           case 10:
            printf( "\n Exciting frequency ..............................."
               "....: %.12e Hz", ( ept->fr[null] ));
            break;

           case 4:
            printf( "\n Exciting frequency ..............................."
               "....: %.12e Hz", ( ept->fr[null] ));

            printf( "\n Transition time .................................."
               "....: %.12e s", ( ept->dt ));
            break;

           case 5:
            printf( "\n " );
            for ( ii=null; ii < ( ept->nn ); ii++ )
            {
               printf( "\n %2ld. exciting frequency ........................"
                  ".......: %.12e Hz", ( ii+ONE ), ( ept->fr[ii] ));
            };
            break;

           case 6: 
	    printf( "\n Sigma time ......................................."
	       "....: %.12e s", ( ept->rt ));
            printf( "\n Delay ............................................"
               "....: %.12e s", ( ept->dt ));
            break;	

           case 7:
           case 14:

            printf( "\n Exciting frequency ..............................."
               "....: %.12e Hz", ( ept->fr[null] ));
            printf( "\n Sigma time ......................................."
               "....: %.12e s", ( ept->rt ));
            printf( "\n Delay ............................................"
               "....: %.12e s", ( ept->dt ));
            break; 

           case 8:

            printf( "\n Rise time ........................................"
               "....: %.12e s", ( ept->rt ));
            break; 

           case 15:

            printf( "\n Exciting frequency ..............................."
               "....: %.12e Hz", ( ept->fr[null] ));
            printf( "\n Transition time .................................."
               "....: %.12e s", ( ept->dt ));
            printf( "\n Plateau time ....................................."
               "....: %.12e s", ( ept->ht ));
            break;
         };
         printf( "\n" );
/*............................................................................*/
# if DSC_HCRMDE != 0

         printf( "\n thermal excitation type: %s", ( ept->hctp ));
        
         switch ( exctp )
         {
           default:
           case 0:
           case 1: /* PASSIVE / FLOATING */
            ( ept->hcn ) = null;
            ( ept->hcdt ) = ZERO;
            ( ept->hcht ) = ZERO;
            ( ept->hcrt ) = ZERO;

            ( ept->nhc ) = null;
            ( ept->ntf ) = null;
            ( ept->ntn ) = null;

            if ( exctp == null )
               printf( "\n Thermal exciation inactive" );
            else
               printf( "\n Floating thermal exciation" );

            break;

           case 2: /* STEADY_STATE */
            if (( ept->hcn ) == null )
               ( ept->hcn ) = TWO;

            printf( "\n Smoothing order .................................."
               "....: %d", ( ept->hcn ));
            printf( "\n Smoothing time ..................................."
               "....: %.12e s", ( ept->hcdt ));
            break;

           case 3: /* DIRAC_PULSE */
            break;

           case 4: /* SMOOTH_RECTANGLE */
            printf( "\n Transition time .................................."
               "....: %.12e s", ( ept->hcdt ));
            break;

           case 5: /* GAUSS_PULSE */
            printf( "\n Sigma time ......................................."
               "....: %.12e s", ( ept->hcrt ));
            printf( "\n Delay ............................................"
               "....: %.12e s", ( ept->hcdt ));
            break;

           case 6: /* RAMP_PULSE */
	    printf( "\n Rise time ........................................"
	       "....: %.12e s", ( ept->hcrt ));
            break;

           case 7: /* HEAVISIDE_STEP */
            break;

           case 8: /* PLATEAU_SMOOTHED___ */
            printf( "\n Delay time ......................................."
               "....: %.12e s", ( ept->hcdt ));
            printf( "\n Hold time ........................................"
               "....: %.12e s", ( ept->hcht ));
            printf( "\n Smoothing order .................................."
               "....: %d", ( ept->hcn ));
            break;

           case 9: /* PERIODIC_PULSE_____ */
            printf( "\n Hold time ........................................"
               "....: %.12e s", ( ept->hcht ));
            printf( "\n Delay time ......................................."
               "....: %.12e s", ( ept->hcdt ));
            break;

           case 10: /* DOUBLE_PERIODIC____ */
            printf( "\n 1st hold time ...................................."
               "....: %.12e s", ( ept->hcht ));
            printf( "\n 1st delay time ..................................."
               "....: %.12e s", ( ept->hcdt ));
            printf( "\n 2nd hold time ...................................."
               "....: %.12e s", ( ept->hcht2 ));
            printf( "\n 2nd delay time ..................................."
               "....: %.12e s", ( ept->hcdt2 ));
            break;

           case 11: /* SWITCHED_SOURCES___ */
            printf( "\n 1st hold time ...................................."
               "....: %.12e s", ( ept->hcht ));
            printf( "\n 1st delay time ..................................."
               "....: %.12e s", ( ept->hcdt ));
            printf( "\n 2nd hold time ...................................."
               "....: %.12e s", ( ept->hcht2 ));
            printf( "\n 2nd delay time ..................................."
               "....: %.12e s", ( ept->hcdt2 ));
            break;
         };
         printf( "\n" );

# endif /* DSC_HCRMDE != 0 */

     /* E_field1: */

         if ( null < ( ept->ne ))
         {
            printf( "\n E-field excitation at %05d ports: >--------"
               "--------------------------------->\n", ( ept->ne ));
            llns = (( ept->ne ) / clns) + ONE;

          /* E_cells: */

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( ept->me[( ept->ne )-ONE] ), bslsh );
                     goto E_ports1;
                  };

                  if ( pp < ( ept->ne )) 
                     printf( "%6ld%c", ( ept->me[pp] ), bslsh );

                  if (( ept->ne ) <= ( ++pp ))  
                     goto E_ports1; 
               };

              E_ports1:

               printf( "\n -> port:" );

               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6d%c\n", ( ept->pe[( ept->ne )-ONE] ), slsh );
                     goto H_field;
                  };

                  if ( pp < ( ept->ne ))   
                     printf( "%6d%c", ( ept->pe[pp] ), slsh );

                  if (( ept->ne ) <= ( ++pp )) 
                  {
                     printf( "\n" );
                     goto H_field;
                  };
               };
            };
         };

        H_field:

         if ( null < ( ept->nh ))
         {
            printf( "\n H-field excitation at %05d ports: >--------"
               "--------------------------------->\n", ( ept->nh ));
            llns = (( ept->nh ) / clns ) + ONE;

        /* h_cells: */

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( ept->mh[( ept->nh )-ONE] ), bslsh );
                     goto H_ports1;
                  };

                  if ( pp < ( ept->nh ))   
                     printf( "%6ld%c", ( ept->mh[pp] ), bslsh );

                  if (( ept->nh ) <= ( ++pp ))  
                     goto H_ports1;
               };

              H_ports1:

               printf( "\n -> port:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6d%c\n", ( ept->ph[( ept->nh )-ONE] ), slsh );
                     goto listo_exc;
                  };

                  if ( pp < ( ept->nh ))   
                     printf( "%6d%c", ( ept->ph[pp] ), slsh );  

                  if (( ept->nh ) <= ( ++pp )) 
                  {
                     printf( "\n" );
                     goto listo_exc; 
                  };
               };
            };
         };

        listo_exc:
# else
         printf( "\n" );
# endif
/*............................................................................*/
/* store excitation mode [parameters]: */
     
         excitfle = fopen( fleptr, "w" );

         printf( "\n opened: Excitation file %s", fleptr );
/*............................................................................*/
/* system identifier [model name] and comment: */

         fprintf( excitfle, spformat, ( ept->name ));
         fprintf( excitfle, spformat, ( ept->text ));

         strcpy( ptr, "DSC-MODEL_EXCITATION_FILE_" );
         strcat( ptr, fleptr );

         pp = LNELNGTH - strlen( ptr );
         qq = null; do
         {
            fprintf( excitfle, "%c", 95 );
         } while(( ++qq ) < pp );
         fprintf( excitfle, "%s\n", ptr );
/*............................................................................*/
         fprintf( excitfle, "%-27s\n", "Maxwell_field_excitation" );

         if ( *( ppt->domain ) == 't' )
            fprintf( excitfle, "%s\n", "TIME_DOMAIN" );
         else /* if *( ppt->domain ) == 'f'requency_domain */
            fprintf( excitfle, "%s\n", "FREQUENCY_DOMAIN" );

         fprintf( excitfle, spformat, ( ept->type ));
         fprintf( excitfle, spformat, "parameters:" );

         switch ( extyp )
         {
           case 0: /* ZERO_EM-FIELD______ */
            break;

           case 1: /* STEADY_STATE */
            fprintf( excitfle, i_format, ( ept->nn )); /* smoothing order */
            fprintf( excitfle, d_format, ( ept->dt )); /* smoothing time  */
            break;

           case 2: /* DIRAC_PULSE */
            break;

           case 3: /* HARMONIC_SINUSOIDAL */
            fprintf( excitfle, d_format, ( ept->fr[null] ));
            break;

           case 4: /* SMOOTH_HARMONIC */
            fprintf( excitfle, d_format, ( ept->fr[null] ));
            fprintf( excitfle, d_format, ( ept->dt ));
            break;
         
           case 5: /* MULTIPLE_HARMONIC */ 
            fprintf( excitfle, i_format, ( ept->nn ));
            for ( ii=null; ii<( ept->nn ); ii++ )
            {
                fprintf( excitfle, d_format, ( ept->fr[ii] ));
            };
            break;

	   case 6: /* GAUSS_PULSE */
	    fprintf( excitfle, d_format,  ( ept->rt ));
	    fprintf( excitfle, d_format,  ( ept->dt ));
	    break;

           case 7: /* OSCILLATORY_GAUSS */
            fprintf( excitfle, d_format, ( ept->fr[null] ));
            fprintf( excitfle, d_format, ( ept->rt ));
            fprintf( excitfle, d_format, ( ept->dt ));
            break;
       
	   case 8: /* RAMP_PULSE */
	    fprintf( excitfle, d_format, ( ept->rt ));
            break;

           case 9: /* RECTANGULAR_PERIODC */
            fprintf( excitfle, d_format, ( ept->fr[null] )); 
            break;

           case 10: /* SAW_TOOTH_PERIODIC */
            fprintf( excitfle, d_format, ( ept->fr[null] ));
            break;

           case 11: /* HEAVISIDE_STEP */
            break;

           case 12: /* DSC_TIMESTEP_PRDC */
            break;
 
           case 13: /* NOISE_AT_RANDOM */
            break;

           case 14: /* MORLET_WAVELET */
            fprintf( excitfle, d_format, ( ept->fr[null] ));
            fprintf( excitfle, d_format, ( ept->rt ));
            fprintf( excitfle, d_format, ( ept->dt ));
            break;

           case 15: /* WAVE_PACKET_HARMNC */
            fprintf( excitfle, i_format, ( ept->nn )); /* smoothing order */
            fprintf( excitfle, d_format, ( ept->dt )); /* transition time */
            fprintf( excitfle, d_format, ( ept->ht )); /* plateau time */
            fprintf( excitfle, d_format, ( ept->fr[null] )); /* frequency */
            break;
	 };

         fprintf( excitfle, "%-27s   %-s\n", "excited_objects", "number" );

         fprintf( excitfle, "%-27s   %-ld\n", "E-ports",
            ( long ) ( ept->ne ));
         fprintf( excitfle, "%-27s   %-ld\n", "H-ports",
            ( long ) ( ept->nh ));
/*............................................................................*/
# if DSC_HCRMDE != 0
/*............................................................................*/
# if DSC_FLDMDE == 0
         fprintf( excitfle, "\n%-27s\n", "Thermal_excitation" );
# else
         fprintf( excitfle, "\n%-27s\n", "Thermal-fluid_excitation" );
# endif
/*............................................................................*/
         fprintf( excitfle, "%-27s\n", "TIME_DOMAIN" );
         fprintf( excitfle, "%-27s\n", ( ept->hctp ));
         fprintf( excitfle, "%-27s\n", "parameters:" );

         switch ( exctp )
         {
           default: /* PASSIVE, FLOATING */
            break;

           case 2: /* STEADY_STATE */
            fprintf( excitfle, i_format, ( ept->hcn )); /* smoothing order */
            fprintf( excitfle, d_format, ( ept->hcdt )); /* smoothing time */
            break;

           case 3: /* DIRAC_PULSE */
            break;

           case 4: /* SMOOTH_PULSE */
            fprintf( excitfle, d_format, ( ept->hcdt )); /* smoothing time */
            break;

           case 5: /* GAUSS_PULSE */
            fprintf( excitfle, d_format, ( ept->hcrt )); /* sigma time */
            fprintf( excitfle, d_format, ( ept->hcdt )); /* delay time */
            break;

           case 6: /* RAMP_PULSE */
            fprintf( excitfle, d_format, ( ept->hcrt )); /* rise time */
            break;

           case 7: /* HEAVISIDE_STEP */
            break;

           case 8: /* PLATEAU_SMOOTHED___ */
            fprintf( excitfle, i_format, ( ept->hcn )); /* smoothing order */
            fprintf( excitfle, d_format, ( ept->hcht )); /* hold time  */
            fprintf( excitfle, d_format, ( ept->hcdt )); /* decline time */
            break;

           case 9: /* PERIODIC_PULSE_____ */
            fprintf( excitfle, d_format, ( ept->hcht )); /* hold time  */
            fprintf( excitfle, d_format, ( ept->hcdt )); /* decline time */
            break;

           case 10: /* DOUBLE_PERIODIC____ */
            fprintf( excitfle, d_format, ( ept->hcht )); /* hold time  */
            fprintf( excitfle, d_format, ( ept->hcdt )); /* decline time */
            fprintf( excitfle, d_format, ( ept->hcht2 )); /* 2nd hold tme  */
            fprintf( excitfle, d_format, ( ept->hcdt2 )); /* 2nd decl tme */
            break;

           case 11: /* SWITCHED_SOURCES___ */
            fprintf( excitfle, d_format, ( ept->hcht )); /* hold time  */
            fprintf( excitfle, d_format, ( ept->hcdt )); /* decline time */
            fprintf( excitfle, d_format, ( ept->hcht2 )); /* 2nd hold tme  */
            fprintf( excitfle, d_format, ( ept->hcdt2 )); /* 2nd decl tme */
            break;
         };

         fprintf( excitfle, "%-27s   %-s\n", "excited_objects", "number" );

         fprintf( excitfle, "%-27s   %-ld\n", "heat_currents:",
            ( long ) ( ept->nhc ));
         fprintf( excitfle, "%-27s   %-ld\n", "face_temperatures:",
            ( long ) ( ept->ntf ));
         fprintf( excitfle, "%-27s   %-ld\n", "node_temperatures:",
            ( long ) ( ept->ntn ));
# endif
/*............................................................................*/
/* from here on, the excited ports, nodes & excitation parameters are stored: */

         fprintf( excitfle, "\n%s\n\n", fclbls );
            /* string "...LABELS_etc....", e.g. */

         for ( ii=null; ii<( ept->ne ); ii++ )
         {
            fprintf( excitfle, "%-ld\n", ( ept->me[ii] ));
            fprintf( excitfle, "%-d\n", ( ept->pe[ii] ));

            if ( *( ppt->domain ) == 't' ) /* time_domain */
               fprintf( excitfle, "%+.16e\n", ( ept->er[ii] ));
            else /* if *( ppt->domain ) == 'f'requency_domain */
               fprintf( excitfle, "%+.16e  %+.16e\n", 
                  ( ept->er[ii] ), ( ept->ei[ii] ));
         };

         for ( ii=null; ii<( ept->nh ); ii++ )
         {
            fprintf( excitfle, "%-ld\n", ( ept->mh[ii] ));
            fprintf( excitfle, "%-d\n", ( ept->ph[ii] ));

            if ( *( ppt->domain ) == 't' ) /* time_domain */
               fprintf( excitfle, "%+.16e\n", ( ept->hr[ii] ));
            else /* if *( ppt->domain ) == 'f'requency_domain */
               fprintf( excitfle, "%+.16e  %+.16e\n",
                  ( ept->hr[ii] ), ( ept->hi[ii] ));
         };
/*............................................................................*/
# if DSC_HCRMDE != 0
         for ( ii=null; ii<( ept->nhc ); ii++ )
         {
            fprintf( excitfle, "%-ld\n", ( ept->mhc[ii] ));
            fprintf( excitfle, "%-d\n", ( ept->fhc[ii] ));
            fprintf( excitfle, "%+.16e\n", ( ept->hc[ii] ));
         };
         for ( ii=null; ii<( ept->ntf ); ii++ )
         {
            fprintf( excitfle, "%-ld\n", ( ept->mtf[ii] ));
            fprintf( excitfle, "%-d\n", ( ept->ftf[ii] ));
            fprintf( excitfle, "%+.16e\n", ( ept->tf[ii] ));
         };
         for ( ii=null; ii<( ept->ntn ); ii++ )
         {
            fprintf( excitfle, "%-ld\n", ( ept->mtn[ii] ));
            fprintf( excitfle, "%+.16e\n", ( ept->tn[ii] ));
         };
# endif
/*............................................................................*/
/* close excitation file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( excitfle, "\nDSC model excitation file %s ", fleptr );
         fprintf( excitfle, "created:\n%s\n", tmestr );

         fclose( excitfle );

         printf( CLEAR_LINE );
         printf( "\r Excitation file %s ", fleptr );
         printf( time_frm, tmestr );
	      
         if (( state->item ) == FORM_ITEMS )
            printf( "\n -----------------------------------"
               "-------------------------------------------");
         break;
/*............................................................................*/
        case 5: /* >---------- evaluation file generation -----------------> */

         if (( state->item ) == FORM_ITEMS )
            printf( "\n\n Evaluation file generation: " );

         strcpy( fleptr, prefix );
	 strcat( fleptr, evlptr );

         ( vpt->n )  = -ONE; /* the number of iteration cycles */
         ( vpt->ni ) = -ONE; /* the first evaluated cycle [ Maxwell field ] */
         ( vpt->nf ) = -ONE; /* the last evaluated cycle [ '' ] */
         ( vpt->r )  = -ONE; /* int. repetition rate (iterations per cycle) */
                             /* [ Maxwell field computation ] */ 
/*............................................................................*/
# if DSC_HCRMDE != 0
         ( vpt->nj ) = -ONE; /* the first evaluated cycle [ thermal - fluid ] */
         ( vpt->nt ) = -ONE; /* the last evaluated cycle [ dito ] */
         ( vpt->rc ) = -ONE; /* int. repetition rate (iterations per cycle) */
                             /* [ thermal - fluid ] */
# endif
/*............................................................................*/
         if (( state->stat ) == FORM_ITEMS )
         {
            strcpy(( vpt->name ), ( tpt->name ));
            strcat( fleptr, ( state->flbl ));
            printf( "\n\n Comment: %.70s \n", ( vpt->text ));
         }
         else if (( state->stat ) != FORM_ITEMS )
         {
            if (( state->item ) < FORM_ITEMS )
            {
               printf( "\n" );

               strcpy(( csp->rqlng ),
                  "Please enter index <N> of evaluation file '" );
               strcat(( csp->rqlng ), fleptr );
               strcat(( csp->rqlng ), "<N>'" );
               strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
               csp = txcnsl( csp );            /* enter index N on console    */
/*...........................................*/
               strcpy( ptr, ( lotos(( csp->inlng ), null )));
               strncpy(( state->flbl ), ptr, THREE );
            }
            else 
               printf( "\n" );

	    strcat( fleptr, ( state->flbl ));

            strcpy(( csp->rqfrm ), "brackets" );
            strcpy(( csp->rqstr ), "Input on keyboard or subroutine " );
            strcat(( csp->rqstr ), "sysbnd(*) ? [k/s]" );
            strcpy(( csp->dfstr ), "s" ); /* default */
/*............................................................................*/
            csp = txcnsl( csp );      /* enter string on text console         */
/*..................................*/
            strcpy( ptr, ( csp->instr ));

            if (( *ptr == 's' )||( *ptr == 'S' ))
            {
/*............................................................................*/
# if FORM_DISP == 1
               printf( "\n DSC model evaluation functions "
                  "started.                        " );
               PRBLDCLR( "\r" );
               printf( "\n %*s", 78, "DSC-FORMER" );
               PRNORMAL( "" );
# endif
/*............................................................................*/
               if (( state->item ) < FORM_ITEMS )
               {
/*............................................................................*/
                  tpt = systop( state );            /* transfer: check DSC -  */
/*................................................*//* mesh topology          */

                  switch( tpt->rtn )
                  {
                    case -ONE:
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;      

                    case ONE:
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in function %s !!!      \n",
                        "systop(*)" );
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;

                    case FIVE: /* = ( state->item ) [ evaluation ] */
                     ( state->stat ) = ( state->item );
                     ( state->tpt ) = tpt;
/*............................................................................*/
# if FORM_DISP == 1
                     printf( "\n ==================================="
                        "===========================================" );
                     PRBLDCLR( "\r" );
                     printf( "\n %*s", 78, "DSC-FORMER" );
                     PRNORMAL( "\r" );
                     printf( "\n DSC model evaluation function "\
                        "%s terminated.                 \n", "sysval(*)" );
# endif
/*............................................................................*/
                     goto identf_check5;

                    default:
                     ( state->tpt ) = tpt;
                     break;
                  };
/*............................................................................*/
                  ppt = syssmx( state );             /* DSC mesh coordinates, */
/*.................................................*//* media etc.            */

                  switch( ppt->rtn )
                  {
                    case -ONE:
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;

                    case ONE:
                     printf( "\n\n Message from function %s:", __func__ );
                     printf( "\n Error in function %s !!!\n", "syssmx(*)" );
                     ( csp->dfopt ) = null; /* the next default menu option */
                     goto menu;

                    default:
                     ( state->ppt ) = ppt;
                     break;
                  };
               };

               ( vpt->rtn ) = -TWO;
/*............................................................................*/
               vpt = sysval( state );            /* DSC mesh evaluation file  */
/*.............................................*//* configuration             */

               switch( vpt->rtn )
               {
                 case -TWO:
                  printf( "\n Message from function %s:      ", __func__ );
                  printf( "\n Unsatisfied code in function %s"
                          " !!! \n", "sysval(*)" );
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 case -ONE:
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 case ONE:
                  printf( "\n\n Message from function %s:    ",
                     __func__ );
                  printf( "\n Error in function %s !!!     \n",
                     "sysval(*)"  );
                  ( csp->dfopt ) = null; /* the next default menu option */
                  goto menu;

                 default:
/*............................................................................*/
# if FORM_DISP == 1
                  printf( "\n ==================================="
                     "===========================================" );
                  PRBLDCLR( "\r" );
                  printf( "\n %*s", 78, "DSC-FORMER" );
                  PRNORMAL( "\r" );
                  printf( "\n DSC model evaluation function "\
                     "%s terminated.                 \n", "sysval(*)" );
# endif
/*............................................................................*/
                  break;
               };

              identf_check5:

               if ( null != strncmp(( tpt->name ), ( vpt->name ), FOUR ) )
               {
                  printf( "\n Message from function %s:", __func__ );
                  printf( "\n Incompatible DSC system identifiers "
                     "on functions %s and %s.",
                     "systop(*)", "syssmx(*)" );
                  printf( "\n Please check functions and restart program ! ");
                  return ONE;
               };

               if (( EVLEP < ( vpt->nep ))
                 ||( EVLEN < ( vpt->nen ))
                 ||( EVLHP < ( vpt->nhp ))
                 ||( EVLHN < ( vpt->nhn )))
               {
                  printf( "\n\n Message from function %s: ", __func__ );

                  if ( EVLEP < ( vpt->nep ))
                  {
                     printf( "\n Too many evaluated E-ports "
                        "in DSC mesh %s !!!", ( tpt->name ));
                     printf( "\n [ Maximum number is %d = macro EVLEP "
                        "in file %s.", EVLEP, "FORMER.CONF" );
                  };
                  if ( EVLEN < ( vpt->nen ))
                  {
                     printf( "\n Too many evaluated E-nodes "
                        "in DSC mesh %s !!!", ( tpt->name ));
                     printf( "\n [ Maximum number is %d = macro EVLEN "
                        "in file %s.", EVLEN, "FORMER.CONF" );
                  };
                  if ( EVLHP < ( vpt->nhp ))
                  {
                     printf( "\n Too many evaluated H-ports "
                        "in DSC mesh %s !!!", ( tpt->name ));
                     printf( "\n [ Maximum number is %d = macro EVLHP "
                        "in file %s.", EVLHP, "FORMER.CONF" );
                  };
                  if ( EVLHN < ( vpt->nhn ))
                  {
                     printf( "\n Too many evaluated H-nodes "
                        "in DSC mesh %s !!!", ( tpt->name ));
                     printf( "\n [ Maximum number is %d = macro EVLHN "
                        "in file %s.", EVLHN, "FORMER.CONF" );
                  };
                  printf( "\n - Change macro only in compliance "
                     "with memory resources " );
                  printf( "\n   and same macro in function values(*) "
                     "of program SOLVER.C.]\n" );
                  return ONE;
               };
/*............................................................................*/
# if DSC_HCRMDE != 0
               if ( EVLHC < ( vpt->nhc ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many heat current faces to be evaluated "
                     "in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro EVLHC "
                     "in %s.", EVLHC, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in header 'SOLVER.CONF' "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( EVLTF < ( vpt->ntf ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many temperature faces to be evaluated "
                     "in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro EVLTF "
                     "in %s.", EVLTF, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in header 'SOLVER.CONF' "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };

               if ( EVLTN < ( vpt->ntn ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many temperature nodes to be evaluated "
                     "in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro EVLTN "
                     "in %s.", EVLTN, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in header 'SOLVER.CONF' "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };
/*............................................................................*/
# if DSC_FLDMDE != 0
               if ( EVLUN < ( vpt->nun ))
               {
                  printf( "\n\n Message from function %s:", __func__ );
                  printf( "\n Too many nodal velocities to be evaluated "
                     "in DSC mesh %s !!!", ( tpt->name ));
                  printf( "\n [ Maximum number is %d = macro EVLUN "
                     "in %s.", EVLUN, __func__ );
                  printf( "\n - Change macro only in compliance with "
                     "memory resources" );
                  printf( "\n   and same macro in header 'SOLVER.CONF' "
                     "of program SOLVER.C.]\n");
                  return ONE;
               };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
            } /* end if( *ptr == 's'ubroutine ) */
            else /* keyboard input */
            {
               printf( "\n Please enter: >--------------------"
                  "------------------->" );
               printf( "\n -> DSC system identifier [ structure name ]"
                  " ..........: " );
               scanf( "%s", vpt->name );
               printf( " -> Comment .........................."
                  ".................: " );
               scanf( "%s", vpt->text );

               printf( "\n" );

               do
               {
                  strcpy( csp->rqlng, "Iteration cycles ? >--"
                     "----------> [ 0 < N, please enter N ]" );
                  strcpy( csp->rqfrm, "points" );
/*............................................................................*/
                  csp = txcnsl( csp );           /* request N on text console */
/*.............................................*/
                  ( vpt->n ) = ( csp->inlng );

               } while(( vpt->n ) <= null );
               do
               {
                  strcpy(( csp->rqlng ), "First evaluated cycle [Mxwf] ? "
                     ">> [ " );
                  strcat(( csp->rqlng ), ( lotos(( long ) VAL_INITLBL, null )));
		  strcat(( csp->rqlng )," <= I <= " );
                  strcat(( csp->rqlng ), ( lotos(( vpt->n ), null )));
		  strcat(( csp->rqlng ),", enter I ]" );
                  strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
                  csp = txcnsl( csp );           /* request I on text console */
/*.............................................*/
                  ( vpt->ni ) = ( csp->inlng );

                  if (( vpt->ni ) < VAL_INITLBL )
                     ( vpt->ni ) = VAL_INITLBL;

               } while((( vpt->ni ) < VAL_INITLBL )
		     ||(( vpt->n ) < ( vpt->ni )));
               do
               {
                  strcpy(( csp->rqlng ), "Last evaluated cycle [Mxwf] ? "
                     ">-> [ " );
                  strcat(( csp->rqlng ), ( lotos(( long ) VAL_INITLBL, null )));
		  strcat(( csp->rqlng )," <= I <= " );
                  strcat(( csp->rqlng ), ( lotos(( vpt->n ), null )));
		  strcat(( csp->rqlng ),", enter I ]" );
                  strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
                  csp = txcnsl( csp );           /* request I on text console */
/*.............................................*/
                  ( vpt->nf ) = ( csp->inlng );
               } while(( vpt->n ) < ( vpt->nf ));

               if (( vpt->r ) == null )
               {
                  strcpy(( csp->rqlng ),
                     "Repetition rate [Maxwell field] ? "
		     ">---> [ 0 <= R, enter R ]" );
                  strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
                  csp = txcnsl( csp );           /* request R on text console */
/*.............................................*/
                  ( vpt->r ) = ( csp->inlng );
               };               
/*............................................................................*/
# if DSC_HCRMDE != 0
/* repetition rates [ thermal - fluid ]: */
# if DSC_INTLCE == 1 /* interlaced internal Maxwell field and thermal */
                      /* computation loops */
	       
               ( vpt->nj ) = ( vpt->ni );
               ( vpt->nt ) = ( vpt->nf );
               ( vpt->rc ) = ( vpt->r );

# elif DSC_INTLCE == 0 /* separated Maxwell field and thermal loops */

               do
               {
                  strcpy(( csp->rqlng ), "First evaluated Th&Fld cycle ? "
                     ">> [ " );
                  strcat(( csp->rqlng ), ( lotos(( long ) VAL_INITLBL, null )));
		  strcat(( csp->rqlng )," <= I <= " );
                  strcat(( csp->rqlng ), ( lotos(( vpt->n ), null )));
		  strcat(( csp->rqlng ),", enter I ]" );
                  strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
                  csp = txcnsl( csp );           /* request I on text console */
/*.............................................*/
                  ( vpt->nj ) = ( csp->inlng );

                  if (( vpt->nj ) < VAL_INITLBL )
                     ( vpt->nj ) = VAL_INITLBL;

               } while((( vpt->nj ) < VAL_INITLBL )
		     ||(( vpt->n ) < ( vpt->nj )));
               do
               {
                  strcpy(( csp->rqlng ), "Last evaluated thermal cycle ? "
                     ">> [ " );
                  strcat(( csp->rqlng ), ( lotos(( long ) VAL_INITLBL, null )));
		  strcat(( csp->rqlng )," <= I <= " );
                  strcat(( csp->rqlng ), ( lotos(( vpt->n ), null )));
		  strcat(( csp->rqlng ),", enter I ]" );
                  strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
                  csp = txcnsl( csp );           /* request I on text console */
/*.............................................*/
                  ( vpt->nt ) = ( csp->inlng );
               } while(( vpt->n ) < ( vpt->nt ));

# endif /* DSC_INTLCE == 0 */

               if (( vpt->rc ) == null )
               {
                  strcpy(( csp->rqlng ),
                     "Repetition rate [ thermal & fluid ] ? "
		     ">-> [ 0 <= R, enter R ]" );
                  strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
                  csp = txcnsl( csp );           /* request R on text console */
/*.............................................*/
                  ( vpt->rc ) = ( csp->inlng );
               };               
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
           /* e_input2: */

               printf( "\n\n E-field evaluation at PORTS: " );
               printf( "\n\n Please enter cell and port indices "
                  ">------------------> " );
               printf( "\n [ Escape: enter ZERO ] >-"
                  "-----------------------------> \n" );

               ( vpt->nep ) = null; do
               {
                  printf( " mesh cell index [ %5d. port ] ?..."
                     "..................: ", (( vpt->nep ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto vep;

                  ( vpt->mep[( vpt->nep )] ) = strtol( ptr, endp, DEC );
                  printf( " port index      [ %5d. port ] ?..."
                     "..................: ", (( vpt->nep ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto vep;

                  ( vpt->pep[( vpt->nep )] ) = strtol( ptr, endp, DEC );
                  ( vpt->nep )++ ;
               } while((( vpt->nep ) < EVLEP )\
                      &&( null < ( vpt->mep[( vpt->nep ) - ONE] )));

               ( vpt->nep )-- ;

              vep:

               printf( "\n Number of evaluated <E|p> ports .."
                  "....................: %-6d ", ( vpt->nep ));

               if ( ONE < ( vpt->nep ))
               {
                  printf( "\n\n individual evaluation / averaging "
                     "over ports ? ( i/a ): " );
                  scanf( "%s", ptr );

                  if (( *ptr == 'a' )||( *ptr == 'A' ))
                     strcpy( vpt->mode_ep, "avarage" );
                  else
                     strcpy( vpt->mode_ep, "individual" );
               };
             
               printf( "\n\n E-field evaluation in NODES: " ); 
               printf( "\n\n Please enter cell indices and components "
                  ">------------> " );
               printf( "\n [ Escape: enter ZERO ] >------------------"
                  "------------> \n" );

               ( vpt->nen ) = null; do
               {
                  printf( " mesh cell index [ %5d. node ] ?........."
                     "............: ", (( vpt->nen ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto ven;

                  ( vpt->men[( vpt->nen )] ) = strtol( ptr, endp, DEC );
                  printf( " component       [ %5d. node ] ?  "
                     "enter u,v or w ....: ", (( vpt->nen ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto ven;

                  ( vpt->cen[( vpt->nen )] ) = ptr[null];
                  ( vpt->nen )++;
               } while((( vpt->nen ) < EVLEN )\
                      &&( null < ( vpt->men[( vpt->nen ) - ONE] )));

               ( vpt->nen )-- ;

              ven:

               printf( "\n Number of evaluated <E|n> nodes ."
                  ".....................: %-6d ", ( vpt->nen ));

               if ( ONE < ( vpt->nen ))
               {
                  printf( "\n\n individual evaluation / averaging "
                     "over nodes ? ( i/a ): " );
                  scanf( "%s", ptr );

                  if (( *ptr == 'a' )||( *ptr == 'A' ))
                     strcpy( vpt->mode_en, "average" );
                  else
                     strcpy( vpt->mode_en, "individual" );
               };
            
           /* h_input2: */

               printf( "\n\n H-field evaluation at PORTS: " );
               printf( "\n\n Please enter cell and port indices "
                  ">------------------> " );
               printf( "\n [ Escape: enter ZERO ] >------------"
                  "------------------> \n" );

               ( vpt->nhp ) = null; do
               {
                  printf( " mesh cell index [ %5d. port ] ?..."
                     "..................: ", (( vpt->nhp ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto vhp;

                  ( vpt->mhp[( vpt->nhp )] ) = strtol( ptr, endp, DEC );
                  printf( " port index      [ %5d. port ] ?...."
                     ".................: ", (( vpt->nhp ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto vhp;

                  ( vpt->php[( vpt->nhp )] ) = strtol( ptr, endp, DEC );
                  ( vpt->nhp )++ ;
               } while((( vpt->nhp ) < EVLHP )\
                      &&( null <( vpt->mhp[( vpt->nhp ) - ONE] )));

               ( vpt->nhp ) -= ONE;

              vhp:

               printf( "\n Number of evaluated <H|p> ports ..."
                  "...................: %-6d ", ( vpt->nhp ));

               if ( ONE < ( vpt->nhp ))
               {
                  printf( "\n\n individual evaluation / averaging "
                     "over ports ? ( i/a ): " );
                  scanf( "%s", ptr );

                  if (( *ptr == 'a' )||( *ptr == 'A' ))
                     strcpy( vpt->mode_hp, "average" );
                  else
                     strcpy( vpt->mode_hp, "individual" );
               };

               printf( "\n\n H-field evaluation in NODES: " );
               printf( "\n\n Please enter cell indices and components "
                  ">------------> " );
               printf( "\n [ Escape: enter ZERO ] >------------------"
                  "------------> \n" );

               ( vpt->nhn ) = null; do
               {
                  printf(" mesh cell index [ %5d. node ] ?........."
                     "............: ", (( vpt->nhn ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto vhn;

                  ( vpt->mhn[( vpt->nhn )] ) = strtol( ptr, endp, DEC );

                  printf( " component       [ %5d. node ] ?  "
                     "enter u,v or w ....: ", (( vpt->nhn ) + ONE ));
                  scanf( "%s", ptr );

                  if ( *ptr == '0' ) 
                     goto vhn;

                  ( vpt->chn[( vpt->nhn )] ) = ptr[null];
                  ( vpt->nhn )++;
               } while((( vpt->nhn ) < EVLHN )\
                      &&( null < ( vpt->mhn[( vpt->nhn ) - ONE] )));

               ( vpt->nhn )-- ;

              vhn:

               printf( "\n Number of evaluated <H|n> nodes .."
                  "....................: %-6d ", ( vpt->nhn ));

               if ( ONE < ( vpt->nhn ))
               {
                  printf( "\n\n individual evaluation / averaging "
                     "over nodes ? ( i/a ): " );
                  scanf( "%s", ptr );

                  if (( *ptr == 'a' )||( *ptr == 'A' ))
                     strcpy( vpt->mode_hn, "average" );
                  else
                     strcpy( vpt->mode_hn, "individual" );
               };
            }; /* end if( *ptr != 's'): keyboard input */ 
         };/* end if (( state->stat ) != FORM_ITEMS ) */
/*............................................................................*/
/* evaluation mode terminated */
/*............................................................................*/
         printf( "\n" );

         while(( vpt->n ) <= null )
         {
            strcpy(( csp->rqlng ), "Iteration cycles ? "
               ">------------> [ 0 < N, please enter N ]" );
            strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
            csp = txcnsl( csp );           /* request N on text console       */
/*.......................................*/
            ( vpt->n ) = ( csp->inlng );
         };

         while((( vpt->ni ) < null )
             ||(( vpt->n ) < ( vpt->ni )))
         {
            strcpy(( csp->rqlng ), "First evaluated cycle ? "
               ">--------> [ 0 <= I <= N, enter I ]" );
            strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
            csp = txcnsl( csp );           /* request I on text console       */
/*.......................................*/
            ( vpt->ni ) = ( csp->inlng );
         };

         if (( vpt->r ) == null )               
         {
            strcpy(( csp->rqlng ),
               "Repetition rate [Maxwell field] ? "
               ">---> [ 0 <= R, enter R ]" );
            strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
            csp = txcnsl( csp );           /* request R on text console */
/*.......................................*/
            ( vpt->r ) = ( csp->inlng );
         };               
/*............................................................................*/
# if DSC_HCRMDE != 0

# if DSC_INTLCE == 1 /* interlaced internal Maxwell field and heat current */
                      /* computation loops */
         ( vpt->nj ) = ( vpt->ni );
         ( vpt->nt ) = ( vpt->nf );
         ( vpt->rc ) = ( vpt->r );
# endif
         if (( vpt->rc ) == null )               
         {
            strcpy(( csp->rqlng ),
               "Repetition rate [heat currents] ? "
               ">---> [ 0 <= R, enter R ]" );
            strcpy(( csp->rqfrm ), "points" );
/*............................................................................*/
            csp = txcnsl( csp );           /* request R on text console */
/*.......................................*/
            ( vpt->rc ) = ( csp->inlng );
         };               
# endif
/*............................................................................*/
# if FORM_DISP == 1
         if (( state->item ) < FORM_ITEMS )
         {
            printf( "\n DSC system identifier [ structure name ]: "
               "%.22s", ( vpt->name ));
            printf( "\n Comment: %.70s \n", ( vpt->text ));
         };
# endif
/*............................................................................*/
         ( state->vpt ) = vpt;

/*............................................................................*/
# if FORM_DISP == 1 /* display evaluated E/H field cells & ports: */

     /* E_field2: */

         if ( null < ( vpt->nep ))
         {
            if (( ONE < ( vpt->nep ))\
              &&(( vpt->mode_ep[null] ) == 'a' ))
               printf( "\n E-field averaging over %5d ports: >-------"
                  "--------------------------------->\n", ( vpt->nep ));
            else 
               printf( "\n E-field evaluation at %5d ports: >--------"
                  "--------------------------------->\n", ( vpt->nep ));

            llns = (( vpt->nep ) / clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( vpt->mep[( vpt->nep )-ONE] ), bslsh );
                     goto E_ports2;
                  };

                  if ( pp < ( vpt->nep ))   
                     printf( "%6ld%c", vpt->mep[pp], bslsh ); 

                  if (( vpt->nep ) <= ( ++pp ))  
                     goto E_ports2; 
               };

              E_ports2:

               printf( "\n -> port:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6d%c\n", ( vpt->pep[( vpt->nep )-ONE] ), slsh );
                     goto E_nodes;
                  };

                  if ( pp < ( vpt->nep ))   
                     printf( "%6d%c", vpt->pep[pp], slsh ); 

                  if (( vpt->nep ) <= ( ++pp )) 
                  {
                     printf( "\n" );
                     goto E_nodes;
                  };
               };
            };
         };/* end if null < ( vpt->nep ) */

        E_nodes:

         if ( null < ( vpt->nen ))
         {
            if (( ONE < ( vpt->nen ))\
              &&(( vpt->mode_en[null] ) == 'a' ))
               printf( "\n E-field averaging over %5d nodes: >-------"
                  "--------------------------------->\n", ( vpt->nen ));
            else
               printf( "\n E-field evaluation at %5d nodes: >--------"
                  "--------------------------------->\n", ( vpt->nen ));

            llns = (( vpt->nen ) / clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> node:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", vpt->men[( vpt->nen )-ONE], bslsh );
                     goto E_comps;
                  };

                  if ( pp < ( vpt->nen ))   
                     printf( "%6ld%c", vpt->men[pp], bslsh ); 

                  if (( vpt->nen ) <= ( ++pp ))  
                     goto E_comps; 
               };

              E_comps: 

               printf( "\n -> vect:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6c%c\n",
                        (( vpt->cen[( vpt->nen )-ONE] )+117 ), slsh );
                     goto H_field2;
                  };

                  if ( pp < ( vpt->nen ))   
                     printf( "%6c%c",
                        (( vpt->cen[pp] )+117 ), slsh );

                  if (( vpt->nen ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto H_field2;
                  };
               };
            };
         };/* end if null< ( vpt->nen ) */

        H_field2:

         if ( null < ( vpt->nhp ))
         {
            if (( ONE < ( vpt->nhp ))\
              &&(( vpt->mode_hp[null] ) == 'a' ))
               printf( "\n H-field averaging over %5d ports: >-------"
                  "--------------------------------->\n", ( vpt->nhp ));
            else
               printf( "\n H-field evaluation at %5d ports: >--------"
                  "--------------------------------->\n", ( vpt->nhp ));

            llns = (( vpt->nhp ) / clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( vpt->mhp[( vpt->nhp )-ONE] ), bslsh );
                     goto H_ports2;
                  };

                  if ( pp < ( vpt->nhp ))   
                     printf( "%6ld%c", vpt->mhp[pp], bslsh ); 

                  if (( vpt->nhp ) <= ( ++pp ))  
                     goto H_ports2; 
               };

              H_ports2:

               printf( "\n -> port:" );

               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6d%c\n", ( vpt->php[( vpt->nhp )-ONE] ), slsh );
                     goto H_nodes; 
                  };

                  if ( pp < ( vpt->nhp ))   
                     printf( "%6d%c", vpt->php[pp], slsh ); 

                  if (( vpt->nhp ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto H_nodes;
                  };
               };
            };
         };/* end if null < ( vpt->nhp ) */

        H_nodes:

         if ( null < ( vpt->nhn ))
         {
            if (( ONE < ( vpt->nhn ))\
              &&(( vpt->mode_hn[null] ) == 'a' ))
               printf( "\n H-field averaging over %05d nodes: >-------"
                  "--------------------------------->\n", ( vpt->nhn ));
            else
               printf( "\n H-field evaluation at %05d nodes: >--------"
                  "--------------------------------->\n", ( vpt->nhn ));

            llns = (( vpt->nhn ) / clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> node:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", vpt->mhn[( vpt->nhn )-ONE], bslsh );
                     goto H_comps;
                  };

                  if ( pp < ( vpt->nhn ))   
                     printf( "%6ld%c", vpt->mhn[pp], bslsh ); 

                  if (( vpt->nhn ) <= ( ++pp ))  
                     goto H_comps; 
               };

              H_comps:

               printf( "\n -> vect:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6c%c\n",
                        (( vpt->chn[( vpt->nhn )-ONE] )+117 ), slsh );
                     goto Maxwfld_listo;
                  };

                  if ( pp < ( vpt->nhn ))   
                     printf( "%6c%c",
                        (( vpt->chn[pp] )+117 ), slsh );

                  if (( vpt->nhn ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto Maxwfld_listo;
                  };
               };
            };
         };/* end if null < ( vpt->nhn ) */

        Maxwfld_listo:
/*............................................................................*/
# if DSC_HCRMDE != 0 /* evaluated thermal objects */
/*............................................................................*/
# if DSC_FCTEMP == 1 /* display face temperatures to be evaluated: */

         if ( null < ( vpt->ntf ))
         {
            printf( "\n temperature evaluation at %05d faces: >---------"
               "---------------------------->\n", ( vpt->ntf ));

            llns = (( vpt->ntf ) / clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c",
                        ( vpt->mtf[( vpt->ntf )-ONE] ), bslsh );
                     goto TF_faces;
                  };

                  if ( pp < ( vpt->ntf ))   
                     printf( "%6ld%c", ( vpt->mtf[pp] ), bslsh ); 

                  if (( vpt->ntf ) <= ( ++pp ))  
                     goto TF_faces; 
               };

              TF_faces:

               printf( "\n -> face:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6ld%c\n",
                        ( long ) ( vpt->ftf[( vpt->ntf )-ONE] ), slsh );

                     goto HC_crrnts;
                  };

                  if ( pp < ( vpt->ntf ))   
                     printf( "%6ld%c", ( long ) ( vpt->ftf[pp] ), slsh );

                  if (( vpt->ntf ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto HC_crrnts;
                  };
               }; /* end for ( ii=ONE; ii<=clns; ii++ ) */
            }; /* end for ( kk=null; kk<llns; ii++ ) */
         }; /* end if null < ( vpt->ntf ) */

        HC_crrnts:
# endif /* DSC_FCTEMP == 1 */
/*............................................................................*/
/* display face currents to be evaluated: */

         if ( null < ( vpt->nhc ))
         {
            printf( "\n heat current evaluation at %05d faces: >--------"
               "---------------------------->\n", ( vpt->nhc ));

            llns = (( vpt->nhc ) / clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c",
                        ( vpt->mhc[( vpt->nhc )-ONE] ), bslsh );
                     goto HC_faces;
                  };

                  if ( pp < ( vpt->nhc ))   
                     printf( "%6ld%c", ( vpt->mhc[pp] ), bslsh ); 

                  if (( vpt->nhc ) <= ( ++pp ))  
                     goto HC_faces; 
               };

              HC_faces:

               printf( "\n -> face:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6ld%c\n", ( long ) \
                        ( vpt->fhc[( vpt->nhc )-ONE] ), slsh );

                     goto Nde_temps;
                  };

                  if ( pp < ( vpt->nhc ))   
                     printf( "%6ld%c", ( long ) ( vpt->fhc[pp] ), slsh );

                  if (( vpt->nhc ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto Nde_temps;
                  };
               }; /* end for ( ii=ONE; ii<=clns; ii++ ) */
            }; /* end for ( kk=null; kk<llns; ii++ ) */
         }; /* end if null < ( vpt->nhc ) */
/*............................................................................*/
/* display node temperatures to be evaluated: */

        Nde_temps:

         if ( null < ( vpt->ntn ))
         {
            printf( "\n temperature evaluation at %05d nodes: >---------"
               "---------------------------->\n", ( vpt->ntn ));

            llns = (( vpt->ntn ) / clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( long ) \
                        ( vpt->mtn[( vpt->ntn )-ONE] ), bslsh );
                     goto T_ports;
                  };

                  if ( pp < ( vpt->ntn ))   
                     printf( "%6ld%c", ( vpt->mtn[pp] ), bslsh ); 

                  if (( vpt->ntn ) <= ( ++pp ))  
                     goto T_ports; 
               };

              T_ports:

               printf( "\n -> port:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6c%c\n", '-', slsh );

                     goto hcrr_listo;
                  };

                  if ( pp < ( vpt->ntn ))   
                     printf( "%6c%c", '-', slsh );

                  if (( vpt->ntn ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto hcrr_listo;
                  };
               }; /* end for ( ii=ONE; ii<=clns; ii++ ) */
            }; /* end for ( kk=null; kk<llns; ii++ ) */
         }; /* end if null < ( vpt->ntn ) */

        hcrr_listo:;
/*............................................................................*/
# if DSC_FLDMDE != 0
/* display nodal velocities to be evaluated: */

         if ( null < ( vpt->nun ))
         {
            printf( "\n flow velocity evaluation at %05d nodes: >-------"
               "---------------------------->\n", ( vpt->nun ));
            llns = (( vpt->nun )/clns )+ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               printf( "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( long ) \
                        ( vpt->mun[(( vpt->nun )-ONE )] ), bslsh );
                     goto U_compts;
                  };

                  if ( pp < ( vpt->nun ))   
                     printf( "%6ld%c", ( vpt->mun[pp] ), bslsh ); 

                  if (( vpt->nun ) <= ( ++pp ))
                     goto U_compts; 
               };

              U_compts:

               printf( "\n -> cmpt:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )\
                    &&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6c%c\n",
                        (( vpt->cun[( vpt->nun )-ONE] )+120 ), slsh );
                     goto flow_listo;
                  };

                  if ( pp < ( vpt->nun ))
                     printf( "%6c%c",
                        (( vpt->cun[pp] )+120 ), slsh );

                  if (( vpt->nun ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto flow_listo;
                  };
               }; /* end for ( ii=ONE; ii<=clns; ii++ ) */
            }; /* end for ( kk=null; kk<llns; ii++ ) */
         }; /* end if null < ( vpt->nun ) */

        flow_listo:;

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
# else /* if DSC_DISP != 1 */
         printf( "\n" );

# endif /* FORM_DISP != 1 */
/*............................................................................*/
/* store operational instructions [ after further consistence checks ]: */

         if (( vpt->ni ) < VAL_INITLBL )
            ( vpt->ni ) = VAL_INITLBL;

         if (( vpt->n ) < ( vpt->ni ))
            ( vpt->n ) = ( vpt->ni );
	 
         if (( vpt->nf ) < VAL_INITLBL )
            ( vpt->nf ) = ( vpt->n );
	 else if (( vpt->n ) < ( vpt->nf ))
            ( vpt->n ) = ( vpt->nf );
/*............................................................................*/
# if DSC_HCRMDE != 0
         if (( vpt->nj ) < VAL_INITLBL )
            ( vpt->nj ) = VAL_INITLBL;

         if (( vpt->n ) < ( vpt->nj ))
            ( vpt->n ) = ( vpt->nj );

         if (( vpt->nt ) < VAL_INITLBL )
            ( vpt->nt ) = ( vpt->n );
	 else if (( vpt->n ) < ( vpt->nt ))
            ( vpt->n ) = ( vpt->nt );
# endif
/*............................................................................*/
         evalfle = fopen( fleptr, "w" );

         printf( "\n opened: Evaluation file %s", fleptr );
/*............................................................................*/
/* system identifier [model name] and comment: */

         fprintf( evalfle, spformat, ( vpt->name ));
         fprintf( evalfle, spformat, ( vpt->text ));

         strcpy( ptr, "DSC-MODEL_EVALUATION_FILE_" );
         strcat( ptr, fleptr );

         pp = LNELNGTH - strlen( ptr );
         qq = null; do
         {
            fprintf( evalfle, "%c", 95 );
         } while(( ++qq ) < pp );
         fprintf( evalfle, "%s\n", ptr );
/*............................................................................*/
/* the total number of iteration cycles: */
         fprintf( evalfle, "%-27s   %-ld\n", "iteration_cycles:",
            ( long ) ( vpt->n ));
/* the first evaluated Maxwell field cycle: */
         fprintf( evalfle, "%-27s   %-ld\n", "1st_evlt'd_cycle,_Mxwfld:",
            ( long ) ( vpt->ni ));
/* the last evaluated Maxwell field cycle: */
         fprintf( evalfle, "%-27s   %-ld\n", "Lst_evlt'd_cycle,_Mxwfld:",
            ( long ) ( vpt->nf ));
/* the internal repetition rate [Maxwell field comput]: */
         fprintf( evalfle, "%-27s   %-ld\n", "iterations/cycle,_Mxwfld:",
            ( long ) ( vpt->r ));
/*............................................................................*/
# if DSC_HCRMDE != 0
/* the first evaluated cycle, thermal&fluid computation: */
         fprintf( evalfle, "%-27s   %-ld\n", "1st_evlt'd_cycle,_Ht&Fld:",
            ( long ) ( vpt->nj ));
/* the last evaluated cycle, thermal&fluid computation: */
         fprintf( evalfle, "%-27s   %-ld\n", "Lst_evlt'd_cycle,_Ht&Fld:",
            ( long ) ( vpt->nt ));
/* the internal repetition rate, thermal&fluid computation: */
         fprintf( evalfle, "%-27s   %-ld\n", "iterations/cycle,_Ht&Fld:",
            ( long ) ( vpt->rc ));
# endif
/*............................................................................*/
/* the following would require running non-interlaced Mxw.-field and thermal */
/* processes [to be separately evaluated] - which is not the present concept.*/
/*
# if DSC_HCRMDE != 0
         fprintf( evalfle, i_format, ( vpt->nc ));
         fprintf( evalfle, i_format, ( vpt->nic ));
         fprintf( evalfle, i_format, ( vpt->rc ));
# endif
*/
/*............................................................................*/
/* evaluated E-ports & nodes: */
         fprintf( evalfle, spformat, "\nMaxwell_field_evaluation" );
         fprintf( evalfle, "%-27s   %-s\n", "evaluated_objects", "number" );

         fprintf( evalfle, "%-27s   %-10ld", "E-ports",
            ( long ) ( vpt->nep ));

         if ( ONE < ( vpt->nep )) 
         {
            if ( null == strstr(( vpt->mode_ep ), "average" ))
               fprintf( evalfle, " %s", "[individual]" );
            else
               fprintf( evalfle, " %s", "[average]" );
         };

         fprintf( evalfle, "\n%-27s   %-10ld", "E-nodes",
            ( long ) ( vpt->nen ));

         if ( ONE < ( vpt->nen ))
         {
            if ( null == strstr(( vpt->mode_en ), "average" ))
               fprintf( evalfle, " %s", "[individual]" );
            else
               fprintf( evalfle, " %s", "[average]" );
         };

         fprintf( evalfle, "\n%-27s   %-10ld", "H-ports",
            ( long ) ( vpt->nhp ));

         if ( ONE < ( vpt->nhp )) 
         {
            if ( null == strstr(( vpt->mode_hp ), "average" ))
               fprintf( evalfle, " %s", "[individual]" );
            else
               fprintf( evalfle, " %s", "[average]" );
         };

         fprintf( evalfle, "\n%-27s   %-10ld", "H-nodes",
            ( long ) ( vpt->nhn ));

         if ( ONE < ( vpt->nhn )) 
         {
            if ( null == strstr(( vpt->mode_hn ), "average" ))
               fprintf( evalfle, " %s", "[individual]" );
            else
               fprintf( evalfle, " %s", "[average]" );
         };

         fprintf( evalfle, "\n" );
/*...........................................................................*/
# if DSC_HCRMDE !=0
/*...........................................................................*/
# if DSC_FLDMDE ==0
         fprintf( evalfle, "\n%-27s\n", "Thermal_evaluation" );
# else
         fprintf( evalfle, "\n%-27s\n", "Thermal-fluid_evaluation" );
# endif /* DSC_FLDMDE ... */
/*...........................................................................*/
         fprintf( evalfle, "%-27s   %-s\n", "evaluated_objects", "number" );
         fprintf( evalfle, "%-27s   %-ld\n", "heat_current:",
            ( long )( vpt->nhc ));
/*...........................................................................*/
# if DSC_FCTEMP == 1
         fprintf( evalfle, "%-27s   %-ld\n", "face_temperatures:",
            ( long ) ( vpt->ntf ));
# endif
/*...........................................................................*/
         fprintf( evalfle, "%-27s   %-ld\n", "node_temperatures:",
            ( long ) ( vpt->ntn ));
/*...........................................................................*/
# if DSC_FLDMDE !=0
         fprintf( evalfle, "%-27s   %-ld\n", "nodal_velocities:",
            ( long ) ( vpt->nun ));
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* store evaluated node & port indices: */

         fprintf( evalfle, "\n%s\n\n", fclbls );
            /* string "...LABELS_etc....", e.g. */
          
         for ( ii=null; ii<( vpt->nep ); ii++ )
            fprintf( evalfle, "%-15ld   %-d\n",
               ( vpt->mep[ii] ), ( vpt->pep[ii] ));

         for ( ii=null; ii<( vpt->nen ); ii++ )
            fprintf( evalfle, "%-15ld   %-d\n",
               ( vpt->men[ii] ), ( vpt->cen[ii] ));

         for ( ii=null; ii<( vpt->nhp ); ii++ )
            fprintf( evalfle, "%-15ld   %-d\n",
               ( vpt->mhp[ii] ), ( vpt->php[ii] ));

         for ( ii=null; ii<( vpt->nhn ); ii++ )
            fprintf( evalfle, "%-15ld   %-d\n",
               ( vpt->mhn[ii] ), ( vpt->chn[ii] ));
/*...........................................................................*/
# if DSC_HCRMDE != 0
         for ( ii=null; ii<( vpt->nhc ); ii++ )
            fprintf( evalfle, "%-15ld   %-d\n",
               ( vpt->mhc[ii] ), ( vpt->fhc[ii] ));
/*...........................................................................*/
# if DSC_FCTEMP == 1
         for ( ii=null; ii<( vpt->ntf ); ii++ )
            fprintf( evalfle, "%-15ld   %-d\n",
               ( vpt->mtf[ii] ), ( vpt->ftf[ii] ));
# endif
/*...........................................................................*/
         for ( ii=null; ii<( vpt->ntn ); ii++ )
            fprintf( evalfle, i_format, ( vpt->mtn[ii] ));
/*...........................................................................*/
# if DSC_FLDMDE !=0
         for ( ii=null; ii<( vpt->nun ); ii++ )
            fprintf( evalfle, "%-15ld   %-d\n",
               ( vpt->mun[ii] ), ( vpt->cun[ii] ));
# endif /* DSC_FLDMDE != 0 */
/*...........................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* close evaluation file with date & time of creation: */ 

         nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( evalfle, "\nDSC model evaluation file %s ", fleptr );
         fprintf( evalfle, "created:\n%s\n", tmestr );

         fclose( evalfle );

         printf( CLEAR_LINE );
         printf( "\r Evalutation file %s ", fleptr );
         printf( time_frm, tmestr ); 

         break;
/*............................................................................*/
        case 6: /* >---------- dsc batch file generation ------------------> */
    
         printf( "\n ======================================"
                  "========================================" );
         PRBLDCLR( "\r" );
         printf( "\n %*s", 78, "DSC-FORMER" );
         PRNORMAL( "\r" );
         printf( "\n DSC batch file creation: " );

         strcpy( fleptr, prefix );
         strcat( fleptr, btchpt );

        number_of_jobs:

         printf( "\n" );
         strcpy(( csp->rqfrm ), "brackets" );
         strcpy( ptr, "Please enter number of jobs " );
         strcat( ptr, "[ 0 <= n <=" );
         strcat( ptr, ( lotos( DSC_JOBS, null )));
         strcat( ptr, " ]" );
         strcpy(( csp->rqlng ), ptr );
         ( csp->dflng ) = ONE;
/*............................................................................*/
         csp = txcnsl( csp );        /*                                       */
/*.................................*/
         ii = ( csp->inlng );

         if ( ii <= null )
            ii = null;
         else if ( DSC_JOBS < ii )
         {
            printf( " Too many jobs !!!  [ maximum number is "
               "%3d ]\n", DSC_JOBS );
            goto number_of_jobs;
         }
         else if ( null < ii )
         {
            printf( "\n" );
            strcpy(( csp->rqfrm ), "brackets" );
            strcpy(( csp->rqstr ), "Relabel DSC system files ? [y/n]" );
            strcpy(( csp->dfstr ), "n" ); /* default */
/*............................................................................*/
            csp = txcnsl( csp );       /*                                     */
/*...................................*/
            strcpy( ptr, ( csp->instr ));

            if (( *ptr == 'y' )||( *ptr == 'Y' ))
            {
              relabel:

               for ( kk=null; kk<ii; kk++ )
               {
                  printf( "\n" );
                  strcpy(( csp->rqfrm ), "brackets" );

                  strcpy( ptr, "Enter " );

		  if ( kk == null )
                     ( csp->dflng ) = kk;
                  else
		     ( csp->dflng ) = joblbl[kk-ONE] + ONE;

                  strcat( ptr, ( lotos(( kk+ONE ), null )));
                  strcat( ptr, ". job index" );
                  strcpy(( csp->rqlng ), ptr );
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  joblbl[kk] = ( csp->inlng );

		  strcpy( ptr, "Enter topology file index of job no." ); 
		  strcat( ptr, ( lotos( joblbl[kk], null )));
                  strcpy(( csp->rqlng ), ptr );
                  ( csp->dflng ) = ( csp->inlng );
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  toplbl[kk] = ( csp->inlng );

                  strcpy( ptr, "Enter s-matrix file index of job no." ); 
		  strcat( ptr, ( lotos( joblbl[kk], null )));
                  strcpy(( csp->rqlng ), ptr );
                  ( csp->dflng ) = ( csp->inlng );
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  smxlbl[kk] = ( csp->inlng );

                  strcpy( ptr, "Enter boundary file index of job no." ); 
                  strcat( ptr, ( lotos( joblbl[kk], null )));
                  strcpy(( csp->rqlng ), ptr );
                  ( csp->dflng ) = ( csp->inlng );
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  bndlbl[kk] = ( csp->inlng );

                  strcpy( ptr, "Enter excitation file index of job no." ); 
                  strcat( ptr, ( lotos( joblbl[kk], null )));
                  strcpy(( csp->rqlng ), ptr );
                  ( csp->dflng ) = ( csp->inlng );
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  exclbl[kk]= ( csp->inlng );
               };

               printf( "\n" );
               strcpy(( csp->rqfrm ), "brackets" );
               strcpy(( csp->rqstr ), "Input correct ? [y/n] " );
               strcat(( csp->dfstr ), "y" ); /* default */
/*............................................................................*/
               csp = txcnsl( csp );       /*                                  */
/*......................................*/
               strcpy( ptr, ( csp->instr ));

               if (( *ptr == 'n' )||( *ptr == 'N' ))
                  goto relabel;
            }
            else /* natural labels */
            {
               kk = DSC_MXJBL - ii;

              first_job:

               strcpy( ptr, "Please enter index of first job " );
               strcat( ptr, "[ 0 <= i <=" );
               strcat( ptr, ( lotos( kk , null )));
               strcat( ptr, " ]" );
               strcpy(( csp->rqlng ), ptr );
               ( csp->dflng ) = null;
/*............................................................................*/
               csp = txcnsl( csp );        /*                                 */
/*.......................................*/
               jj = ( csp->inlng );

               if (( jj < null )||( kk < jj ))
               {
                  printf( " illegal job index !!! "
                     "[ 0 <= i <= %3d ]", kk );
                  goto first_job;
               };

               for ( kk=null; kk<ii; kk++ )
               {
                  joblbl[kk] = jj;
                  toplbl[kk] = jj;
                  smxlbl[kk] = jj;
                  bndlbl[kk] = jj;
                  exclbl[kk] = jj;
                  jj++ ;
               };
            };
/*............................................................................*/
/* open batch file, type DSC_PRFX.batch: */

            batchfle = fopen( fleptr, "w" );

            printf( "\n opened: Batch file %s", fleptr );

            fprintf( batchfle, "%s", "DSC_batchfile" );

            if (( *ptr == 'y' )||( *ptr == 'Y' )) /* still ptr = csp->instr */
               fprintf( batchfle, "\n%s", rmvfle ); /* cf.last call: txcnsl(*)*/
            else
               fprintf( batchfle, "\n%s", rtnfle );

            fprintf( batchfle, "\n%s %3d", "number_of_jobs:", ( short ) ii );

            fprintf( batchfle, "\n\n%s", "job " );

            fprintf( batchfle, " %s", DSC_PRFX );
            fprintf( batchfle, "%s%s", TOPOLOGY_FILE, "[_]" );
            fprintf( batchfle, " %s", DSC_PRFX );
            fprintf( batchfle, "%s%s", S_MATRIX_FILE, "[_]" );
            fprintf( batchfle, " %s", DSC_PRFX );
            fprintf( batchfle, "%s%s", BOUNDARY_FILE, "[_]" );
            fprintf( batchfle, " %s", DSC_PRFX );
            fprintf( batchfle, "%s%s", EXCITATION_FILE, "[_]" );
            fprintf( batchfle, " %s", DSC_PRFX );
            fprintf( batchfle, "%s%s", EVALUATION_FILE, "[_]" );

            fprintf( batchfle, " %s", "state" );
            fprintf( batchfle, " %s", "CPU_time/err" );
               
            for ( kk=null; kk<ii; kk++ )
            {
               fprintf( batchfle, "\n%3d", joblbl[kk] );
               fprintf( batchfle, " %10d", toplbl[kk] );
               fprintf( batchfle, " %10d", smxlbl[kk] );
               fprintf( batchfle, " %10d", bndlbl[kk] );
               fprintf( batchfle, " %10d", exclbl[kk] );
               fprintf( batchfle, " %10d", joblbl[kk] );

               fprintf( batchfle, " %6s", "--" );
               fprintf( batchfle, " %12s", "------------" );
            };
/*............................................................................*/
/* close batch file with date & time of creation: */ 

            nseconds = time( timer );

/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */
            strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

            TIMEFORM( tmestr, ctmptr );

            fprintf( batchfle, "\n\n%s %s %s\n%s\n", "DSC batch file",
               fleptr, "created:", tmestr );

            fclose( batchfle );

            printf( CLEAR_LINE );
            printf( "\r DSC batch file %s ", fleptr );
            printf( time_frm, tmestr );
            ( csp->dfopt ) = null; /* usually, here terminates input */
            goto menu;
         }
         else /* case ii=null [ number of jobs ] */
         {
            ( csp->dfopt ) = null; /* the next default menu option */
            goto menu;
         };

         break;
/*............................................................................*/

        default:
         break;
      }; /* end switch ( idx ) */
   }; /* next idx */

   strcpy(( csp->cmmnt ), "Welcome back to DSC-FORMER !" );
   ( csp->dfopt ) = 7; /* the next default menu option */
   goto menu;
}  
/*============================================================================*/
# undef VERTX_LOG
# undef UNLNK_LOG
# undef EMWLL_LOG
# undef COORD_LOG
# undef CSIZE_LOG
# undef MEDIA_LOG
# undef SMTRX_LOG
# undef FORM_DISP
# undef CELL_POINTS
# undef PRECISION
# undef SMALL_VAL
# undef GIANT_VAL
# undef VAL_INITLBL
# undef BND_HCRSMX
# undef BND_CSHAPE
# undef LNELNGTH
/*********************** end of function 'formdrv(*)' *************************/
