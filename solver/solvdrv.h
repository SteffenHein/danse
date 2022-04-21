/* [ file: solvdrv.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function solvdrv(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   This function controles the DSC process, joblbl[jj] on the basis of all    *
*   batch file instructions  [ file elf.batch ], i.e. reading the DSC system   *
*   files [ dsc.top<*>, dsc.smx<*>, dsc.bnd<*>, dsc.exc<*>, dsc.val<*> ], as   *
*   these are relabeled by main program SOLVER.C                               *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* [ yet included in solver.c ]
# include "../CONFIG.H"
# include "./SOLVER.CONF"
*/
/*----------------------------------------------------------------------------*/
/* A macro that gives a conveniently readable form to the string returned */
/* with [ the pointer of ] function ctime(&nseconds) */
# include "../tools/TIMEFORM.M"
/*----------------------------------------------------------------------------*/
# define POLAK_RIBIERE__ 1
# define FLETCHER_REEVES 0
/*----------------------------------------------------------------------------*/
/* the following macros should be defined in "../CONFIG.H": */

# ifndef DSC_ADJCLL      /* assign neighbouring cell index top.mn[][k] so:    */
   # define DSC_ADJCLL 0 /* 0: to faces [ k is a face index; 0 <= k < 6 ]     */
# endif                  /* 1: to ports [ k is a port index; 0 <= k < 12 ]    */
/*----------------------------------------------------------------------------*/
# ifndef DSC_FCELBL      /* label neighbouring face index top.fn[][] so:      */
   # define DSC_FCELBL 0 /* 0: unsigned, started with index null              */
# endif                  /* 1: unsigned, started with index ONE               */
                         /* 2: signed, started with index ONE                 */
                         /* [ the negative/positive sign indicating opposed/  */
                         /*   aligned orientation to neighbouring cell face ] */
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
   # ifndef DSC_INTLCE        /* 1: interlaced internal Maxwell field and */
      # define DSC_INTLCE 1   /* heat current computation loops */
   # endif                    /* 0: separate internal loops [allows different */
                              /* repetition rates for Maxw and hcurr comput] */ 
/*----------------------------------------------------------------------------*/
   # ifndef DSC_ORTHO        /* 0: [ precise ] non-orthogonal scheme */
      # define DSC_ORTHO 0   /* 1/2: [ approximate ] quasi-orthogonal scheme */
   # endif                 
/*----------------------------------------------------------------------------*/
   # ifndef DSC_HCGAUGE      /* 2: separate gauge heat currents */
      # define DSC_HCGAUGE 0 /* 1: separate incident heat currents */
   # endif                    /* 0: most efficient code */
/*----------------------------------------------------------------------------*/
   # ifndef DSC_FCTEMP
      # define DSC_FCTEMP == 1 /* 1: store face temperatures */
   # endif
# endif
/*----------------------------------------------------------------------------*/
/* time/frequency domain compiler flags:                                      */

# ifndef DSC_DOMAIN       
   # define DSC_DOMAIN 0 /* 1: time domain, 2: frequ.domain, 
                                          0: either [ time or frequ. domain ] */
# endif
/*----------------------------------------------------------------------------*/
/* 'soft' dimensional constants [ may be modified, depending on application   */
/*  - in compliance with memory resources ]:                                  */

# ifndef NODES
   # define NODES 30000 /* maximum number of mesh cells                       */
# endif
# ifndef NRSMX
   # define NRSMX 30000 /* maximum number of Maxwell field s-parameter sets   */
# endif
# ifndef GENDS
   # define GENDS     1 /* maximum number of gyroelectric nodes               */
# endif
# ifndef GESMX
   # define GESMX     1 /* maximum number of gyroelectric s-parameter sets    */
# endif
# ifndef GMNDS
   # define GMNDS     1 /* maximum number of gyromagnetic nodes               */
# endif
# ifndef GMSMX
   # define GMSMX     1 /* maximum number of gyromagnetic s-parameter sets    */
# endif
# ifndef BNDAP
   # define BNDAP  1000 /* maximum number of aperiodic boundary faces         */
# endif
# ifndef BNDPR
   # define BNDPR     1 /* maximum number of periodic boundary faces          */
# endif
# ifndef EXCEP
   # define EXCEP  1000 /* maximum number of excited E ports                  */
# endif
# ifndef EXCHP
   # define EXCHP  1000 /* maximum number of excited H ports                  */
# endif
# ifndef EXCFR
   # define EXCFR     5 /* maximum number of exc. frequencies [type:MULT._H.] */
# endif
# ifndef EVLEP
   # define EVLEP  1000 /* maximum number of evaluated E ports                */
# endif
# ifndef EVLEN
   # define EVLEN     1 /* maximum number of evaluated E nodes                */
# endif
# ifndef EVLHP
   # define EVLHP     1 /* maximum number of evaluated H ports                */
# endif
# ifndef EVLHN
   # define EVLHN     1 /* maximum number of evaluated H nodes                */
# endif
/*............................................................................*/
/* define thermal extensions, if required [ and not yet done in CONFIG.H ]: */
# ifndef DSC_HCRMDE
# define DSC_HCRMDE  1
   # define HCNDS NODES
   # define HCSMX NODES /* maximum number of heat & fluid s-parameter sets    */
   # define BNDHC     1 /* maximum number of heat current boundary faces      */
   # define BNDTF BNDAP /* maximum number of temperature boundary faces       */
   # define BNDTN     1 /* maximum number of imposed node temperatures        */
   # define EXCHC     1 /* maximum number of excited heat current nodes       */
   # define EXCTF     1 /* maximum number of excited temperature faces        */
   # define EXCTN     1 /* maximum number of excited temperature nodes        */
   # define EVLHC     1 /* maximum number of evaluated heat current faces     */
   # define EVLTF     1 /* maximum number of evaluated temperature faces      */
   # define EVLTN     2 /* maximum number of evaluated temperature nodes      */
   # ifndef DSC_FCTEMP
      # define DSC_FCTEMP 1 /* 1: store face temperatures */
   # endif
/*............................................................................*/
/* define fluid extensions, if required [ and not yet done in CONFIG.H ]: */
# ifndef DSC_FLDMDE
# define DSC_FLDMDE  1
   # define NFCNN     1 /* maximum number of fluids [ connected components ] */
   # define BNDSL BNDAP /* maximum number of free slip boundary faces */
   # define BNDOF     1 /* maximum number of outflow boundary faces */
   # define BNDIF     1 /* maximum number of imposed boundary face flows */
   # define EVLUN     1 /* maximum number of evaluated nodal flows */
   # define EVLUF     1 /* maximum number of evaluated cell face flows */
# endif
# endif
/*----------------------------------------------------------------------------*/
# ifndef SLV_ALERT
   # define SLV_ALERT 1 /* 1: produce audible signal at crucial program steps*/
# endif
/*----------------------------------------------------------------------------*/
/* 'hard' dimensional constants [ usually def in solvtp.h - don't modify !!! ]
# define DIMNS  3 
# define PORTS 12 *//* number of external ports per cell *//*
# define FACES  6 *//* number of cell faces *//*
# define SENTR  6 *//* number of S-matrix enries *//*
# define HDIMS  9
# define STUBS  6 *//* number of internal ports ['stubs'] per cell *//*       */
/*----------------------------------------------------------------------------*/
/* 'hard' operational constants [ don't modify !!! ]:                         */

# define ELECTRIC_WALL (-1)
# define MAGNETIC_WALL (-2)
/*----------------------------------------------------------------------------*/
/* structures:                                                                */

# include "./solvtp.h" /* solv typedefs */
# include "./dsptyp.h" /* typedef of dsplay(*) [messg/running cursor function]*/
/*----------------------------------------------------------------------------*/
static struct solverstat solver = {null};
static struct topology top = {null};
static struct tlmsmx smx = {null};
static struct boundary bnd = {null};
static struct excitation exc = {null};
static struct evaluation val = {null};
/*............................................................................*/
# if DSC_HCRMDE != 0
static struct hcrsmx hcs = {null};
static DSC_HCRRTS hcr[DSC_HCRMDE] = {{null}};
static CLUSTER temperature = {null};
/*............................................................................*/
# if DSC_FLDMDE != 0
static CLUSTER pressure = {null};
static CLUSTER velocity[THREE] = {{null}};
/*............................................................................*/
# if CMPRSSBL != 0
static CLUSTER density = {null};
static CLUSTER flow[THREE] = {{null}};
# endif /* CMPRSSBL != 0 */
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# if DSC_HCGAUGE == 2
   static DSC_GGEHCR ggc[DSC_HCRMDE] = {{null}};
# endif
/*............................................................................*/
# endif
/*----------------------------------------------------------------------------*/
# if DSC_DOMAIN == 1

static DSC_FIELDS fld = {null};
static DSC_GGEFLD ggf = {null};
static DSC_DEFLCT dfl = {null};
/*----------------------------------------------------------------------------*/
# elif DSC_DOMAIN == 2

static DSC_FIELDS fld[TWO] = {{null}};
/*----------------------------------------------------------------------------*/
# else /* default: DSC_DOMAIN == 0 */

static DSC_FIELDS fld[TWO] = {{null}};
static DSC_DEFLCT dfl = {null};

# endif /* DSC_DOMAIN == ... */
/*----------------------------------------------------------------------------*/
# ifdef _Include
   # include "./scattr.h"
# endif
# include "./toplgy.h"
# include "./smtrix.h"
# include "./boundr.h"
# include "./excite.h"
# include "./trnsnt.h"
# include "./values.h"
# include "./clearv.h"
/*----------------------------------------------------------------------------*/
# ifndef PRECISION
   # define PRECISION ( double )( 1.000e-14 )
# endif
/*============================================================================*/

short solvdrv( short jj )
{
/*............................................................................*/
/* allusions: */
/*
   extern struct solverstat solver;

# if DSC_DOMAIN == 1
   extern DSC_FIELDS fld;
   extern DSC_GGEFLD ggf;
   extern DSC_DEFLCT dfl;
# else
   extern DSC_FIELDS fld[];
   extern DSC_DEFLCT dfl;
# endif
*/
/*............................................................................*/
/* declarations [types]: */

   static long 
      kk = null,
      mm = null,
      nn = null,
      rr = null;

   static short
      hh = null,
      ii = null,
      pp = null,
      sgn = null,
      ind = null;
/*............................................................................*/
# if DSC_HCRMDE != 0
   static short
      fc = null,
      cc = null;
# endif
/*............................................................................*/
   static char
      tmestr[STS_SIZE] = {null};
/*
   static char
      ptr[STS_SIZE] = {null};
*/
   static const char
      start  = 's',
      restrt = 'm',
      cursor = 'c',
      percent = 'p';
   
   static const signed char
      eport[THREE][FOUR] = {{ 0, 1, 2, 3}, { 4, 5, 6, 7}, { 8, 9, 10, 11}},
      hport[THREE][FOUR] = {{-10, 5, 12,-7}, {-2, 9, 4,-11}, {-6, 1, 8,-3}};
/*............................................................................*/
# if DSC_DOMAIN != 2
   static long 
      mp = null;
# endif
/*............................................................................*/
# if DSC_STOPTIME == 1
   time_t 
      nseconds = null,
     *timer = null;

   static char
      ctmptr[STS_SIZE] = {null};
# endif
/*............................................................................*/
/* static FILE *display = stdout; */
/*............................................................................*/
/* declarations [structures]: */

   static struct solverstat
     *state = &solver;

   static DSPLAY
     *dsp = NULL;

   static DSC_FIELDS
     *inc = NULL,
     *out = NULL;
/*............................................................................*/
# if DSC_DOMAIN == 0
   static DSC_FIELDS
     *gfp = &fld[ONE];

   static DSC_DEFLCT
     *dfp = &dfl;
# elif DSC_DOMAIN == 1
   static DSC_GGEFLD
     *gfp = &ggf;

   static DSC_DEFLCT
     *dfp = &dfl;
# endif
/*............................................................................*/
# if DSC_HCRMDE != 0
   static struct hcrsmx
      *hsp = NULL;
   static DSC_HCRRTS
      *hci[DSC_HCRMDE] = {NULL};
   static DSC_HCRRTS
      *hre[DSC_HCRMDE] = {NULL};
   static CLUSTER 
      *tmp = NULL;
/*............................................................................*/
# if DSC_FLDMDE != 0
   static CLUSTER
      *prs = NULL,
      *vel[THREE] = { NULL };
/*............................................................................*/
# if CMPRSSBL != 0
   static CLUSTER 
      *dns = NULL,
      *flw[THREE] = { NULL };
# endif /* if CMPRSSBL != 0 */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# if DSC_HCGAUGE == 2
   static DSC_GGEHCR
      *gcp[DSC_HCRMDE] = {NULL};
# endif
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* system prototypes: */
# ifndef _CCBUG
   char 
     *strcat( char *ptr1, const char *ptr2 ),
     *strncat( char *ptr1, const char *ptr2, size_t n );
# endif

# if DSC_STOPTIME == 1   
   time_t time( time_t *timer );
   char *ctime( const time_t *timer );
# endif
/*............................................................................*/
/* user prototypes: */

   char 
     *lotos ( long lngint, char length );

   DSPLAY
     *dsplay( DSPLAY *dsp );

   short 
      trnfld( const double ttime ),
      trnhcr( const double hctme ),
      values( const short jj );
/*............................................................................*/
/* include scattering/connection maps [if not compiled separately as objects] */
# ifdef _Include
   DSC_FIELDS
     *excfld( void ),
     *scatfld( void ),
     *cnctfld( void );
/*............................................................................*/
# if DSC_HCRMDE != 0
   DSC_HCRRTS
     *exchcr( const short cc ),
     *scathcr( const short cc ),
     *cncthcr( const short cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
   DSC_HCRRTS
     *sctflow( const short cc ),
     *conflow( const short cc );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
# else /* if _Include not defined */
   DSC_FIELDS
     *excfld( struct solverstat *epp ),
     *scatfld( struct solverstat *epp ),
     *cnctfld( struct solverstat *epp );
/*............................................................................*/
# if DSC_HCRMDE != 0
   DSC_HCRRTS
     *exchcr( struct solverstat *epp ),
     *scathcr( struct solverstat *epp ),
     *cncthcr( struct solverstat *epp );
/*............................................................................*/
# if DSC_FLDMDE != 0
   DSC_HCRRTS
     *sctflow( struct solverstat *epp ),
     *conflow( struct solverstat *epp );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*----------------------------------------------------------------------------*/
/* initialize struct solverstat: */

   ( state->mxwnn ) = -ONE; /* reset Maxwell field iterations counter */

   ( state->tpt ) = &top;
   ( state->spt ) = &smx;
   ( state->bpt ) = &bnd;
   ( state->ept ) = &exc;
   ( state->vpt ) = &val;
/*............................................................................*/
# if (( DSC_DOMAIN == 0 )\
    ||( DSC_DOMAIN == 1 ))
   ( state->gfp ) = gfp;
   ( state->dfp ) = dfp;
# endif
/*............................................................................*/
# if DSC_HCRMDE != 0
   ( state->hcrnn ) = null; /* reset thermal iterations counter */

   hsp = &hcs;
   ( state->hsp ) = hsp;
   tmp = &temperature;
   ( tmp->hsp ) = &hcs;
   ( tmp->par ) = 't'; /* [ parameter: 't'emperature ] */
   ( state->tmp ) = tmp;

   cc = null; do
   {
      hci[cc] = &hcr[cc];
      ( state->hci[cc] ) = hci[cc];
      hre[cc] = &hcr[cc];
      ( state->hre[cc] ) = hre[cc];

# if DSC_HCGAUGE == 2
      gcp[cc] = &ggc[cc];
      ( state->gcp[cc] ) = gcp[cc];
# endif

   } while(( ++cc ) < DSC_HCRMDE );
/*............................................................................*/
# if DSC_FLDMDE != 0

   ( state->fldnn ) = null; /* reset fluid iterations counter */

   prs = &pressure;
   ( prs->hsp ) = &hcs;
   ( prs->par ) = 'p'; /* [ parameter: 'p'ressure ] */
   ( state->prs ) = prs;

   cc = null; do
   {
     vel[cc] = &velocity[cc];
     ( vel[cc]->hsp ) = &hcs;
     ( vel[cc]->par ) = 'v'; /* [ parameter: 'v'elocity ] */
     ( state->vel[cc] ) = vel[cc];
   } while(( ++cc ) < THREE );
/*............................................................................*/
# if CMPRSSBL != 0

   dns = &density;
   ( dns->hsp ) = &hcs;
   ( dns->par ) = 'd'; /* [ parameter: 'd'ensity ] */
   ( state->dns ) = dns;

   cc = null; do
   {
     flw[cc] = &flow[cc];
     ( flw[cc]->hsp ) = &hcs;
     ( flw[cc]->par ) = 'f'; /* [ parameter: 'f'low ] */
     ( state->flw[cc] ) = flw[cc];
   } while(( ++cc ) < THREE );

# endif /* CMPRSSBL != 0 */
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* initialize display function: */

   dsp = dsplay( null );

   if( *( state->opt ) != null )
   {
      display = fopen(( state->logfle ), "r+" );
      fseek( display, ( state->logps ), SEEK_SET );
   };

   setvbuf( display, null, unbuff, bufsize );
   ( dsp->display ) = display;
/*............................................................................*/
/* launch the DSC process:                                                    */

   strcpy( syslbl, lotos( joblbl[jj], null ));
   strcpy(( dsp->messge ), "DSC process no " );
   strcat(( dsp->messge ), syslbl );
   strncat(( dsp->messge ), " launched", 9 );
/*............................................................................*/
# if DSC_STOPTIME == 1
   nseconds = time( timer );
   strcpy( tmestr, ctime( &nseconds ) + 11 );
   strncat(( dsp->messge ), " at ", 4 );
   strncat(( dsp->messge ), tmestr, 8 );
# endif
/*............................................................................*/
   if( *( state->opt ) == null )
      strncat(( dsp->messge ), " - please don't disturb !!!", 28 );
   
   ( state->ttime ) = ZERO;
   ( state->swing ) = ZERO;
/*............................................................................*/
# if DSC_HCRMDE != 0 
   ( state->sws ) = null;
   ( state->hctme ) = ZERO;
   ( state->hcswg ) = ZERO;
   ( state->hsk ) = ZERO;
/*............................................................................*/
# if DSC_FLDMDE != 0
/* the process log file [ only for fluid computations; */
/* stores flow parameters computed in function conflow(*) ] */

   strcpy(( state->dsclog ), "dsc.log" );
   strcat(( state->dsclog ), syslbl );

   ( state->dscstr ) = fopen(( state->dsclog ), "w+" );

   fprintf(( state->dscstr ), "SOR_logfile" );
   fprintf(( state->dscstr ), "\nDSC_process_no_%s", syslbl );
   fprintf(( state->dscstr ), "\n" );

   kk = null; do /* 74 characters '-' */
   {
      fprintf(( state->dscstr ), "----------" );
   } while(( ++kk ) < SEVEN );

   fprintf(( state->dscstr ), "-----\n" );
   fprintf(( state->dscstr ), "| iteration " );
   fprintf(( state->dscstr ), "| time_[sec] " );
   fprintf(( state->dscstr ), "| SOR_steps " );
   fprintf(( state->dscstr ), "| sum_dw*dw  " );
   fprintf(( state->dscstr ), "| max_|df|   " );
   fprintf(( state->dscstr ), "| cell_no  |" );
   fprintf(( state->dscstr ), "\n" );

   kk = null; do /* fill line up with '-' to 78 characters */
   {
      fprintf(( state->dscstr ), "----------" );
   } while(( ++kk ) < SIX );
   fprintf(( state->dscstr ), "--- DANSE-%s\n", DSC_RELEASE );

   fflush( state->dscstr );
/*............................................................................*/
/* the pressure log file [ only for fluid computations; */
/* stores extrema of pressure computed in conflow(*) ] */

   strcpy(( state->prslog ), "prs.log" );
   strcat(( state->prslog ), syslbl );

   ( state->prsstr ) = fopen(( state->prslog ), "w+" );

   fprintf(( state->prsstr ), "PRESSURE_logfile" );
   fprintf(( state->prsstr ), "\nDSC_process_no_%s", syslbl );
   fprintf(( state->prsstr ), "\n" );

   kk = null; do /* up to 78 characters '-' */
   {
      fprintf(( state->prsstr ), "----------" );
   } while(( ++kk ) < SEVEN );

   fprintf(( state->prsstr ), "--------\n" );
   fprintf(( state->prsstr ), "| iteration " );
   fprintf(( state->prsstr ), "| time_[sec] " );
   fprintf(( state->prsstr ), "| minimum_[Pa] " );
   fprintf(( state->prsstr ), "| cell_no  " );
   fprintf(( state->prsstr ), "| maximum_[Pa] " );
   fprintf(( state->prsstr ), "| cell_no  |" );
   fprintf(( state->prsstr ), "\n" );

   kk = null; do /* fill line up with '-' to 78 characters */
   {
      fprintf(( state->prsstr ), "----------" );
   } while(( ++kk ) < SIX );
   fprintf(( state->prsstr ), "------ DANSE-%s\n", DSC_RELEASE );

   fflush( state->prsstr );
/*............................................................................*/
# if SLV_TMPDSP == 1
/* store temperature field [nodal fluid temperature - as computed in function */
/* sctflow(*)] */

   strcpy(( state->tmpdsp ), "tmp.dsp" );
/*............................................................................*/
# if SLV_JOBLBLS != 0
   strcat(( state->tmpdsp ), syslbl );
# endif
/*............................................................................*/
   ( state->tmpfst ) = fopen(( state->tmpdsp ), "w+" );

   fprintf(( state->tmpfst ), "TEMPERATURE_field" );
   fprintf(( state->tmpfst ), "\nDSC_process_no_%s", syslbl );
   fprintf(( state->tmpfst ), "\n" );

   kk = null; do
   {
      fprintf(( state->tmpfst ), "----------" );
   } while(( ++kk ) < THREE );

   fprintf(( state->tmpfst ), "---\n" );
   fprintf(( state->tmpfst ), "| cell index |" );

   if ( TEMPGGE < 1.0e-77 )
      fprintf(( state->tmpfst ), "  temperature [K] |" );
   else
      fprintf(( state->tmpfst ), "  temperature [C] |" );

   fprintf(( state->tmpfst ), "\n" );

   kk = null; do
   {
      fprintf(( state->tmpfst ), "----------" );
   } while(( ++kk ) < TWO );
   fprintf(( state->tmpfst ), "- DANSE-%s", DSC_RELEASE );

   ( state->tmpfp ) = ftell( state->tmpfst );

   fclose( state->tmpfst );

# endif /* SLV_TMPDSP == 1 */
/*............................................................................*/
# if SLV_PRSDSP == 1
/* store pressure field [nodal fluid pressure - as computed in function */
/* sctflow(*)] */

   strcpy(( state->prsdsp ), "prs.dsp" );
/*............................................................................*/
# if SLV_JOBLBLS != 0
   strcat(( state->prsdsp ), syslbl );
# endif
/*............................................................................*/
   ( state->prsfst ) = fopen(( state->prsdsp ), "w+" );

   fprintf(( state->prsfst ), "PRESSURE_field___" );
   fprintf(( state->prsfst ), "\nDSC_process_no_%s", syslbl );
   fprintf(( state->prsfst ), "\n" );

   kk = null; do
   {
      fprintf(( state->prsfst ), "----------" );
   } while(( ++kk ) < THREE );

   fprintf(( state->prsfst ), "---\n" );
   fprintf(( state->prsfst ), "| cell index |" );
   fprintf(( state->prsfst ), "   pressure [Pa]  |" );
   fprintf(( state->prsfst ), "\n" );

   kk = null; do
   {
      fprintf(( state->prsfst ), "----------" );
   } while(( ++kk ) < TWO );
   fprintf(( state->prsfst ), "- DANSE-%s", DSC_RELEASE );

   ( state->prsfp ) = ftell( state->prsfst );

   fclose( state->prsfst );
# endif /* SLV_PRSDSP == 1 */
/*............................................................................*/
# if CMPRSSBL != 0
# if SLV_DNSDSP == 1
/* store density field [nodal fluid densuty - as computed in function */
/* sctflow(*)] */

   strcpy(( state->dnsdsp ), "dns.dsp" );
/*............................................................................*/
# if SLV_JOBLBLS != 0
   strcat(( state->dnsdsp ), syslbl );
# endif
/*............................................................................*/
   ( state->dnsfst ) = fopen(( state->dnsdsp ), "w+" );

   fprintf(( state->dnsfst ), "DENSITY_field____" );
   fprintf(( state->dnsfst ), "\nDSC_process_no_%s", syslbl );
   fprintf(( state->dnsfst ), "\n" );

   kk = null; do
   {
      fprintf(( state->dnsfst ), "----------" );
   } while(( ++kk ) < THREE );

   fprintf(( state->dnsfst ), "---\n" );
   fprintf(( state->dnsfst ), "| cell index |" );
   fprintf(( state->dnsfst ), " density [kg/m^3] |" );
   fprintf(( state->dnsfst ), "\n" );

   kk = null; do
   {
      fprintf(( state->dnsfst ), "----------" );
   } while(( ++kk ) < TWO );
   fprintf(( state->dnsfst ), "- DANSE-%s", DSC_RELEASE );

   ( state->dnsfp ) = ftell( state->dnsfst );

   fclose( state->dnsfst );
# endif /* SLV_DNSDSP == 1 */
# endif /* CMPRSSBL != 0 */
/*............................................................................*/
# if SLV_FLWDSP == 1
/* store fluid velocity field [nodal fluid velocity - as computed in function */
/* sctflow(*)] */

   strcpy(( state->veldsp ), "flw.dsp" );
/*............................................................................*/
# if SLV_JOBLBLS != 0
   strcat(( state->veldsp ), syslbl );
# endif
/*............................................................................*/
   ( state->velfst ) = fopen(( state->veldsp ), "w+" );

   fprintf(( state->velfst ), "VELOCITY_field___" );
   fprintf(( state->velfst ), "\nDSC_process_no_%s", syslbl );
   fprintf(( state->velfst ), "\n" );

   kk = null; do
   {
      fprintf(( state->velfst ), "----------" );
   } while(( ++kk ) < FIVE );

   fprintf(( state->velfst ), "------\n" );
   fprintf(( state->velfst ), "| cell index |" );
   fprintf(( state->velfst ), "  u_x [m/s]  |" );
   fprintf(( state->velfst ), "  u_y [m/s]  |" );
   fprintf(( state->velfst ), "  u_z [m/s]  |" );

   fprintf(( state->velfst ), "\n" );
   kk = null; do
   {
      fprintf(( state->velfst ), "----------" );
   } while(( ++kk ) < FOUR );
   fprintf(( state->velfst ), "---- DANSE-%s", DSC_RELEASE );

   ( state->velfp ) = ftell( state->velfst );

   fclose( state->velfst );

# endif /* SLV_FLWDSP == 1 */
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   if(( state->dmn ) == 't' ) /* [ 't'ime domain process ] */
   {
/*............................................................................*/
# if DSC_DOMAIN == 2 
      fprintf( display, "\n Error message from function %s :", __func__ );
      fprintf( display, "\n\n Solver is compiled in frequency "
         "domain mode 'DSC_DOMAIN = 2'" );
      fprintf( display, "\n - but started in time domain mode 1 !!!" );
      fprintf( display, "\n\n [ Change time/frequency domain macro DSC_DOMAIN "
         "in SOLVER.CONF from 2 to" );
      fprintf( display, "\n   1: time domain, or" );
      fprintf( display, "\n   0: time & frequency domain." );
      fprintf( display, "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[solver_compiled_"
         "in_mode_DSC_DOMAIN=2] ", LGS_SIZE );
      return null;
/*............................................................................*/
# else /* DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
   # if DSC_DOMAIN == 1
      inc = &fld;
      out = &fld;
   # else
      inc = &fld[null];
      out = &fld[null];
   # endif
/*............................................................................*/
      ( state->inc ) = inc;
      ( state->out ) = out;
/*............................................................................*/
/* start message - display on screen:                                         */

      if( *( state->opt ) == null )
      {
         ( dsp->option ) = start;
         dsplay( dsp );

         ( dsp->option ) = restrt;
         strcpy(( dsp->messge ), "[Time_domain_DSC_process]" );
         dsplay( dsp );
         ( dsp->option ) = cursor;
      }
      else
      {
         fprintf( display, " %s", ( dsp->messge ));
         ( state->logps ) = ftell( display );
         ( dsp->fleps ) = ftell( display );
         ( dsp->option ) = percent;
      };
/*............................................................................*/
/* here starts the time domain process: */
/* 
*//* external loop [ val.n iterations ]
*/
      ( dsp->range ) = val.n;
      ( dsp->messge[null] ) = null;

      nn = ONE;
      while( nn <= val.n )
      {
         ( state->nn ) = nn;
/*............................................................................*/
/*
*//* internal loop, Maxwell field computation computation 
*//* [ averaging over val.r repetitions ]
*/
/*............................................................................*/
         kk = null;
         while(( ++kk ) <= val.r )
         { 
/*............................................................................*/
/* excitation cycle:                                                          */

            if ( nn <= val.nf )
            {
/*............................................................................*/
               ind = trnfld( state->ttime );         /* compute swing */
/*.................................................*/
# ifdef _Include      
               if ( ind == ONE )
                  inc = excfld( );
# else /* if _Include not defined */
               if ( ind == ONE )
                  inc = excfld( state );
# endif /* _Include not defined */
/*............................................................................*/
/* evaluate ports [Maxwell field]: */

               if ( val.ni <= nn )
               {
                  rr = null; 
                  while( rr < val.nep )
                  {
                     mm = val.mep[rr];
                     pp = val.pep[rr]-ONE;
                     val.epr[rr] += (( inc->i[mm][pp] )+( out->r[mm][pp] ));

                     if ( null < bnd.p )
                     {
                        mm += top.n;
                        val.epi[rr] += (( inc->i[mm][pp] )+( out->r[mm][pp] ));
                     };
                     rr++ ;
                  }; /* next rr */

                  rr = null; 
                  while( rr < val.nhp )
                  {
                     mm = val.mhp[rr];
                     pp = val.php[rr]-ONE;
                     val.hpr[rr] += (( inc->i[mm][pp] )-( out->r[mm][pp] ));

                     if ( null < bnd.p )
                     {
                        mm += top.n;
                        val.hpi[rr] += (( inc->i[mm][pp] )-( out->r[mm][pp] ));
                     };
                     rr++ ;
                  }; /* next rr */
               }; /* end if ( val.ni <= nn) */
            }; /* end if ( nn <= val.nf ) */
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_INTLCE == 1
/*............................................................................*/
            if ( nn <= val.nt )
            {
               cc = null; do
               {
/*............................................................................*/
                  ind = trnhcr( state->hctme );      /* compute current swing */
/*.................................................*/
# ifdef _Include       
                  if ( ind == ONE )
                     hci[cc] = exchcr( cc );
# else /* if _Include not defined */
                  ( state->hclbl ) = cc;

                  if ( ind == ONE )
                     hci[cc] = exchcr( state );
# endif /* _Include not defined */
               } while(( ++cc ) < DSC_HCRMDE );
/*............................................................................*/
/* evaluate ports [ heat and fluids ]: */

               if ( val.nj <= nn )
               {
                  cc = null; do
                  {
                     rr = null;
                     while( rr < val.nhc[cc] )
                     {
                        mm = val.mhc[cc][rr];
                        fc = val.fhc[cc][rr];
/*............................................................................*/
# if DSC_HCGAUGE == 0
                        val.hc[cc][rr] += ( hci[cc]->ic[mm][fc] );
# elif DSC_HCGAUGE == 1
                        val.hc[cc][rr] += (( hci[cc]->ic[mm][fc] ) - \
                                              ( hci[cc]->rc[mm][fc] ));
# elif DSC_HCGAUGE == 2 
                        val.hc[cc][rr] += (( hci[cc]->ic[mm][fc] ) - \
                                           ( hci[cc]->rc[mm][fc] ) + \
                                           ( gcp[cc]->gc[mm][fc] ));
# endif
/*............................................................................*/
                        rr++ ;
                     }; /* end while... */
/*............................................................................*/
# if DSC_FCTEMP == 1 /* face temperature evaluation */
                     rr = null;
                     while( rr < val.ntf[cc] )
                     {
                        mm = val.mtf[cc][rr];
                        fc = val.ftf[cc][rr];
                        val.tf[cc][rr] += ( hre[cc]->tf[mm][fc] );
                        rr++ ;
                     };
# endif /* DSC_FCETEMP == 1 */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
               }; /* end if ( val.nj <= nn ) */
            }; /* end if ( nn <= val.nt ) */
# endif /* DSC_INTLCE == 1 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* reflection cycle [ Maxwell field ]: */

            if ( nn <= val.nf )
            {
               ( state->ttime ) += (( state->dt )/2.);
/*............................................................................*/
# ifdef _Include       
               out = scatfld( );
# else /* _Include not defined */
               out = scatfld( state );
# endif /* _Include not defined */
/*............................................................................*/
/* evaluate nodes [Maxwell field]: */

               if ( val.ni <= nn )
               {
                  rr = null;
                  while( rr < val.nen )
                  {
                     mm = val.men[rr];
                     ii = val.cen[rr];

                     hh = null; do
                     {
                        pp = eport[ii][hh];
                        val.enr[rr] += \
                           (( inc->i[mm][pp] )+( out->r[mm][pp] ));

                        if ( null < bnd.p )
                        {
                           mp = mm+top.n;
                           val.eni[rr] += \
                              (( inc->i[mp][pp] )+( out->r[mp][pp] ));
                        }; 
                     } while(( ++hh ) < FOUR );
                     rr++ ;
                  };

                  rr = null;
                  while( rr < val.nhn )
                  {
                     mm = val.mhn[rr];
                     ii = val.chn[rr];

                     hh = null; do
                     {
                        pp = abs( hport[ii][hh] );
                        sgn = hport[ii][hh] / pp;
                        pp -= ONE;

                        val.hnr[rr] += \
                           sgn*(( inc->i[mm][pp] )-( out->r[mm][pp] ));

                        if ( null < bnd.p )
                        {
                           mp = mm+top.n;
                           val.hni[rr] += \
                              sgn*(( inc->i[mp][pp] )-( out->r[mp][pp] ));
                        }; 
                     } while(( ++hh ) < FOUR );
                     rr++ ;
                  };
	       }; /* end if ( val.ni <= nn ) */
	    }; /* end if ( nn <= val.nf ) */
/*............................................................................*/
# if DSC_HCRMDE != 0     
# if DSC_INTLCE == 1 
/*............................................................................*/
/* reflection cycle [ heat and fluids ]: */

            if ( nn <= val.nt )
            {
               ( state->hctme ) += (( state->hcdt )/2.);

               cc = null; do
               {
/*............................................................................*/
# ifdef _Include       
                  if (( state->hcrstart ) <= ( state->hctme ))
                     hre[cc] = scathcr( cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                     hre[cc] = sctflow( cc );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# else /* _Include not defined */
                  if (( state->hcrstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hre[cc] = scathcr( state );
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hre[cc] = sctflow( state );
                  };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
/*............................................................................*/
/* evaluate nodes [ heat and fluids ]: */

               if ( val.nj <= nn )
               {
                  cc = null; do
                  {
                     rr = null;
                     while( rr < val.ntn[cc] )
                     {
                        mm = val.mtn[cc][rr];
                        val.tn[cc][rr] += ( hre[cc]->tn[mm] );
                        rr++ ;
                     };
/*............................................................................*/
# if DSC_FLDMDE != 0 
                     rr = null;
                     while( rr < val.nun[cc] )
                     {
                        mm = val.mun[cc][rr];
                        pp = val.cun[cc][rr];
                        val.un[cc][rr] += ( hre[cc]->un[mm][pp] );
                        rr++ ;
                     };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
               }; /* end if ( val.nj <= nn ) */
            }; /* end if ( nn <= val.nt ) */
# endif /* DSC_INTLCE == 1 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* connection cycle [ Maxwell field ]: */

            if ( nn <= val.nf )
            {
               ( state->ttime ) += (( state->dt )/2.);
/*............................................................................*/
# ifdef _Include       
               inc = cnctfld( );
# else /* _Include not defined */
               inc = cnctfld( state );
# endif /* _Include not defined */
/*............................................................................*/
	    }; /* end if ( nn <= val.nf ) */
/*............................................................................*/
# if DSC_HCRMDE != 0 
# if DSC_INTLCE == 1 
/*............................................................................*/
/* connection cycle [ heat and fluids ]: */

            if ( nn <= val.nt )
            {
               ( state->hctme ) += (( state->hcdt )/2.);

               cc = null; do
               {
/*............................................................................*/
# ifdef _Include       
                  if (( state->hcrstart ) <= ( state->hctme ))
                     hci[cc] = cncthcr( cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                     hci[cc] = conflow( cc );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# else /* _Include not defined */
                  if (( state->hcrstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hci[cc] = cncthcr( state );
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hci[cc] = conflow( state );
                  };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
	    }; /* end if ( nn <= val.nt ) */
/*............................................................................*/
# endif /* DSC_INTLCE == 1 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         }; /* end while(( ++kk ) <= val.r ) */
/*............................................................................*/
/*
*//* end of internal loop [Maxwell field computation]
*/
/*............................................................................*/
# if DSC_HCRMDE != 0 
# if DSC_INTLCE == 0
/*
*//* internal loop, heat current computation 
*//* [ averaging over val.rc repetitions ]
*/
/*............................................................................*/
         kk = null;
         while(( ++kk ) <= val.rc )
         { 
/*............................................................................*/
/* excitation [ heat and fluids ]: */

            if ( nn <= val.nt )
            {
               cc = null; do
               {
/*............................................................................*/
                  ind = trnhcr( state->hctme );      /* compute current swing */
/*.................................................*/
# ifdef _Include       
                  if ( ind == ONE )
                     hci[cc] = exchcr( cc );
# else /* if _Include not defined */
                  ( state->hclbl ) = cc;

                  if ( ind == ONE )
                     hci[cc] = exchcr( state );
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
/*............................................................................*/
/* face evaluation [ heat and fluids ]: */

               if ( val.nj <= nn )
               {
                  cc = null; do
                  {
                     rr = null;
                     while( rr < val.nhc[cc] )
                     {
                        mm = val.mhc[cc][rr];
                        fc = val.fhc[cc][rr];
/*............................................................................*/
# if DSC_HCGAUGE == 0
                        val.hc[cc][rr] += ( hci[cc]->ic[mm][fc] );
# elif DSC_HCGAUGE == 1
                        val.hc[cc][rr] += (( hci[cc]->ic[mm][fc] ) - \
                                           ( hci[cc]->rc[mm][fc] ));
# elif DSC_HCGAUGE == 2 
                        val.hc[cc][rr] += (( hci[cc]->ic[mm][fc] ) - \
                                           ( hci[cc]->rc[mm][fc] ) + \
                                           ( gcp[cc]->gc[mm][fc] ));
# endif
/*............................................................................*/
                        rr++ ;
                     }; /* end while... */
/*............................................................................*/
# if DSC_FCTEMP == 1 /* face temperature evaluation */
                     rr = null;
                     while( rr < val.ntf[cc] )
                     {
                        mm = val.mtf[cc][rr];
                        fc = val.ftf[cc][rr];
                        val.tf[cc][rr] += ( hre[cc]->tf[mm][fc] );
                        rr++ ;
                     };
# endif /* DSC_FCETEMP == 1 */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
               }; /* end if ( val.nj <= nn ) */
/*............................................................................*/
/* reflection cycle [ heat and fluids ]: */
	       
               ( state->hctme ) += (( state->hcdt )/2.);

               cc = null; do
               {
/*............................................................................*/
# ifdef _Include       
                  if (( state->hcrstart ) <= ( state->hctme ))
                     hre[cc] = scathcr( cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                     hre[cc] = sctflow( cc );   
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# else /* _Include not defined */
                  if (( state->hcrstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hre[cc] = scathcr( state );
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hre[cc] = sctflow( state );
                  };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
/*............................................................................*/
/* evaluate nodes [ heat and fluids ]: */

               if ( val.nj <= nn )
               {
                  cc = null; do
                  {
                     rr = null;
                     while( rr < val.ntn[cc] )
                     {
                        mm = val.mtn[cc][rr];
                        val.tn[cc][rr] += ( hre[cc]->tn[mm] );
                        rr++ ;
                     };
/*............................................................................*/
# if DSC_FLDMDE != 0 
                     rr = null;
                     while( rr < val.nun[cc] )
                     {
                        mm = val.mun[cc][rr];
                        pp = val.cun[cc][rr];
                        val.un[cc][rr] += ( hre[cc]->un[mm][pp] );
                        rr++ ;
                     };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
               }; /* end if ( val.nj <= nn ) */
/*............................................................................*/
/* connection cycle [ heat and fluids ]: */

               ( state->hctme ) += (( state->hcdt )/2.);

               cc = null; do
               {
/*............................................................................*/
# ifdef _Include       
                  if (( state->hcrstart ) <= ( state->hctme ))
                     hci[cc] = cncthcr( cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                     hci[cc] = conflow( cc );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# else /* _Include not defined */
                  if (( state->hcrstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hci[cc] = cncthcr( state );
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hci[cc] = conflow( state );
                  };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
            }; /* end if ( nn <= val.nt ) */
         }; /* end while(( ++kk ) <= val.rc ); internal hcrr loop terminated */
/*............................................................................*/
/*
*//* end of internal loop [heat current]
*/
/*............................................................................*/
# endif /* DSC_INTLCE == 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         ind = values( -ONE );        /* store values                         */
/*..................................*/
/* display process state [ running cursor or percent value ]: */

         ( dsp->state ) = nn;
         dsplay( dsp );

         nn++;
      }; /* while( nn <= val.n ) [ val.n iterations terminated ] */
/*............................................................................*/
/* 
*//* end of external loop - job terminated
*/
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
   } /* end of time domain DSC process */
/*............................................................................*/
   else /* if ( *state->dmn == 'f'requency domain ) */
   {
/*............................................................................*/
# if DSC_DOMAIN == 1
      fprintf( display, "\n Error message from function %s :", __func__ );
      fprintf( display, "\n\n Solver is compiled in time "
         "domain mode 'DSC_DOMAIN = 1'" );
      fprintf( display, "\n - but started in frequency domain mode !!!" );
      fprintf( display, "\n\n [ Change time/frequency domain macro DSC_DOMAIN "
         "in SOLVER.CONF from 1 to" );
      fprintf( display, "\n   0: time & frequency domain "
         "( relaxation scheme ) or " );
      fprintf( display, "\n   2: frequency domain "
         "( relaxation scheme )." );
      fprintf( display, "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[solver_compiled_"
         "in_mode_DSC_DOMAIN=1] ", LGS_SIZE );
      return null;
/*............................................................................*/
# else /* DSC_DOMAIN = 0 or 2 */
/*............................................................................*/
/* frequency domain [ relaxation scheme ]: */

      inc = &fld[null];
      out = &fld[ONE];

      ( state->inc ) = inc;
      ( state->out ) = out;
/*............................................................................*/
/* display start message on screen: */

      if( *( state->opt ) == null )
      {
         ( dsp->option ) = start;
         dsplay( dsp );
         ( dsp->option ) = restrt;
         strcpy( dsp->messge, "[Frequency_domain_DSC_iteration]" );
         dsplay( dsp );
         ( dsp->option ) = cursor;
      }
      else
      {
         fprintf( display, " %s", ( dsp->messge ));
         ( state->logps ) = ftell( display );
         ( dsp->fleps ) = ftell( display );
         ( dsp->option ) = percent;
      };
/*............................................................................*/
/* frequency domain DSC algorithm: */
/* 
*//* external loop [ val.n + ONE iterations ]
*/
      ( dsp->range ) = val.n;
      ( dsp->messge[null] ) = null;

      nn = null;
      while( nn <= val.n )
      {
         ( state->nn ) = nn;
/*............................................................................*/
/* 
*//* internal loop, Maxwell field computation computation 
*//* [ averaging over val.r repetitions ]
*/
/*............................................................................*/
         kk = null;
         while(( ++kk ) <= val.r )
         {
/*............................................................................*/
/* excitation [ Maxwell field ]: */

            if ( nn <= val.nf )
            {
/*............................................................................*/
               ind = trnfld( state->ttime );           /* compute swing */
/*...................................................*/
# ifdef _Include       
               if ( ind == ONE )
                  inc = excfld( );
# else /* _Include not defined */
               if ( ind == ONE )
                  inc = excfld( state );
# endif /* _Include not defined */
/*............................................................................*/
/* evaluate ports [ Maxwell field ]: */

               if ( val.ni <= nn )
               {
                  rr = null;
                  while( rr < val.nep )
                  {
                     mm = val.mep[rr];
                     pp = val.pep[rr]-ONE;

                     val.epr[rr] += \
                        (( inc->r[mm][pp] ) + ( out->r[mm][pp] ));
                     val.epi[rr] += \
                        (( inc->i[mm][pp] ) + ( out->i[mm][pp] ));

                     rr++ ;
                  }; /* next rr */

                  rr = null;
                  while( rr < val.nhp )
                  {
                     mm = val.mhp[rr];
                     pp = val.php[rr]-ONE;

                     val.hpr[rr] += \
                        (( inc->r[mm][pp] ) - ( out->r[mm][pp] ));
                     val.hpi[rr] += \
                        (( inc->i[mm][pp] ) - ( out->i[mm][pp] ));
                     rr++ ;
                  }; /* next rr */
	       }; /* end if ( val.ni <= nn ) */
            }; /* end if ( nn <= val.nf ) */
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_INTLCE == 1 
/* excitation [ heat and fluids ]: */

            if ( nn <= val.nt )
            {
/*............................................................................*/
               ind = trnhcr( state->hctme );         /* compute current swing */
/*.................................................*/
               if ( ind == ONE )
               {
                  cc = null; do
                  {
/*............................................................................*/
# ifdef _Include       
                     hci[cc] = exchcr( cc );
# else /* _Include not defined */
                     ( state->hclbl ) = cc;
                     hci[cc] = exchcr( state );
# endif /* _Include not defined */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
               }; /* end if ( ind == ONE ) */
/*............................................................................*/
/* evaluate ports [ heat and fluids ]: */

               if ( val.nj <= nn )
               {
                  cc = null; do
                  {
                     rr = null; 
                     while( rr < val.nhc[cc] )
                     {
                        mm = val.mhc[cc][rr];
                        fc = val.fhc[cc][rr];
/*............................................................................*/
# if DSC_HCGAUGE == 0
                        val.hc[cc][rr] += ( hci[cc]->ic[mm][fc] );
# elif DSC_HCGAUGE == 1
                        val.hc[cc][rr] += (( hci[cc]->ic[mm][fc] ) - \
                                           ( hci[cc]->rc[mm][fc] ));
# elif DSC_HCGAUGE == 2 
                        val.hc[cc][rr] += (( hci[cc]->ic[mm][fc] ) - \
                                           ( hci[cc]->rc[mm][fc] ) + \
                                           ( gcp[cc]->gc[mm][fc] ));
# endif
/*............................................................................*/
                        rr++ ;
                     };
/*............................................................................*/
# if DSC_FCTEMP == 1
                     rr = null; 
                     while( rr < val.ntf[cc] )
                     {
                        mm = val.mtf[cc][rr];
                        fc = val.ftf[cc][rr];
                        val.tf[cc][rr] += ( hre[cc]->tf[mm][fc] );
                        rr++ ;
                     };
# endif /* DSC_FCETEMP == 1 */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
	       }; /* end if ( val.nj <= nn ) */
	    }; /* end if ( nn <= val.nt ) */
# endif /* DSC_INTLCE == 1 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* reflection cycle [ Maxwell field ]: */

	    if ( nn <= val.nf ) 
            {
               ( state->ttime ) += (( state->dt )/2.);
/*............................................................................*/
# ifdef _Include       
               out = scatfld( );     
# else /* _Include not defined */
               out = scatfld( state );
# endif /* _Include not defined */
/*............................................................................*/
/* evaluate nodes [ Maxwell field ]: */

               if ( val.ni <= nn )
               {
                  rr = null;
                  while( rr < val.nen )
                  {
                     mm = val.men[rr];
                     ii = val.cen[rr];

                     hh = null; do
                     {
                        pp = eport[ii][hh];

                        val.enr[rr] += \
                           (( inc->r[mm][pp] ) + ( out->r[mm][pp] ));
                        val.eni[rr] += \
                           (( inc->i[mm][pp] ) + ( out->i[mm][pp] )); 
                     } while(( ++hh ) < FOUR );
                     rr++ ;
                  };

                  rr = null;
                  while( rr < val.nhn )
                  {
                     mm = val.mhn[rr];
                     ii = val.chn[rr];

                     hh = null; do
                     {
                        pp = abs( hport[ii][hh] );
                        sgn = hport[ii][hh] / pp;
                        pp -= ONE;

                        val.hnr[rr] += \
                           sgn*(( inc->r[mm][pp] )-( out->r[mm][pp] ));
                        val.hnr[rr] += \
                           sgn*(( inc->i[mm][pp] )-( out->i[mm][pp] ));
                     } while(( ++hh ) < FOUR );
                     rr++ ;
                  };
	       }; /* end if ( val.ni <= nn ) */
	    }; /* end if ( nn <= val.nf ) */
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_INTLCE == 1
/* reflection cycle [ heat and fluids ]: */

	    if ( nn <= val.nt ) 
            {
               ( state->hctme ) += (( state->hcdt )/2.);

               cc = null; do
               {
/*............................................................................*/
# ifdef _Include       
                  if (( state->hcrstart ) <= ( state->hctme ))
                     hre[cc] = scathcr( cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                     hre[cc] = sctflow( cc );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# else /* _Include not defined */
                  if (( state->hcrstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hre[cc] = scathcr( state );
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hre[cc] = sctflow( state );
                  };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
/*............................................................................*/
/* evaluate nodes [ heat and fluids ]: */

	       if ( val.nj <= nn ) 
               {
                  cc = null; do
                  {
                     rr = null;
                     while( rr < val.ntn[cc] )
                     {
                        mm = val.mtn[cc][rr];
                        val.tn[cc][rr] += ( hre[cc]->tn[mm] );
                        rr++ ;
                     };
/*............................................................................*/
# if DSC_FLDMDE != 0 
                     rr = null;
                     while( rr < val.nun[cc] )
                     {
                        mm = val.mun[cc][rr];
                        pp = val.cun[cc][rr];
                        val.un[cc][rr] += ( hre[cc]->un[mm][pp] );
                        rr++ ;
                     };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
	       }; /* end if ( val.nj <= nn ) */
	    }; /* end if ( nn <= val.nt ) */
# endif /* DSC_INTLCE == 1 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* connection cycle [ Maxwell field ]: */

            if ( nn <= val.nf ) 
            {
               ( state->ttime ) += (( state->dt )/2.);
/*............................................................................*/
# ifdef _Include       
               inc = cnctfld( );
# else /* _Include not defined */
               inc = cnctfld( state );
# endif /* _Include not defined */
/*............................................................................*/
	    }; /* end if ( nn <= val.nf ) */
/*............................................................................*/
# if DSC_HCRMDE != 0 
# if DSC_INTLCE == 1
/* connection cycle [ heat and fluids ]: */

            if ( nn <= val.nt ) 
            {
               ( state->hctme ) += (( state->hcdt )/2.);

               cc = null; do
               {
/*............................................................................*/
# ifdef _Include       
                  if (( state->hcrstart ) <= ( state->hctme ))
                     hci[cc] = cncthcr( cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                     hci[cc] = conflow( cc );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# else /* _Include not defined */
                  if (( state->hcrstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hci[cc] = cncthcr( state );
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hci[cc] = conflow( state );
                  };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
	    }; /* end if ( nn <= val.nt ) */
# endif /* DSC_INTLCE == 1 */
# endif /* endif DSC_HCRMDE != 0 */
/*............................................................................*/
         }; /* end while(( ++kk ) <= val.rc ); internal hcrr loop terminated */
/*............................................................................*/
/*
*//* end of internal loop
*/
/*............................................................................*/
# if DSC_HCRMDE != 0     
# if DSC_INTLCE == 0
/* 
*//* internal loop, heat current computation 
*//* [ averaging over val.rc repetitions ]
*/
/*............................................................................*/
         kk = null;
	 while(( ++kk ) <= val.rc )
         {
/*............................................................................*/
/* excitation [ heat and fluids ]: */

            if ( nn <= val.nt )
            {
/*............................................................................*/
               ind = trnhcr( state->hctme );         /* compute current swing */
/*.................................................*/
               if ( ind == ONE )
               {
                  cc = null; do
                  {
/*............................................................................*/
# ifdef _Include       
                     hci[cc] = exchcr( cc );
# else /* _Include not defined */
                     ( state->hclbl ) = cc;
                     hci[cc] = exchcr( state );
# endif /* _Include not defined */
                  } while(( ++cc ) < DSC_HCRMDE );
               }; /* end if ( ind == ONE ) */
/*............................................................................*/
/* evaluate ports [ heat and fluids ]: */

               if ( val.nj <= nn )
               {
                  cc = null; do
                  {
                     rr = null; 
                     while( rr < val.nhc[cc] )
                     {
                        mm = val.mhc[cc][rr];
                        fc = val.fhc[cc][rr];
/*............................................................................*/
# if DSC_HCGAUGE == 0
                        val.hc[cc][rr] += ( hci[cc]->ic[mm][fc] );
# elif DSC_HCGAUGE == 1
                        val.hc[cc][rr] += (( hci[cc]->ic[mm][fc] ) - \
                                           ( hci[cc]->rc[mm][fc] ));
# elif DSC_HCGAUGE == 2 
                        val.hc[cc][rr] += (( hci[cc]->ic[mm][fc] ) - \
                                           ( hci[cc]->rc[mm][fc] ) + \
                                           ( gcp[cc]->gc[mm][fc] ));
# endif
/*............................................................................*/
                        rr++ ;
                     };
/*............................................................................*/
# if DSC_FCTEMP == 1
                     rr = null; 
                     while( rr < val.ntf[cc] )
                     {
                        mm = val.mtf[cc][rr];
                        fc = val.ftf[cc][rr];
                        val.tf[cc][rr] += ( hre[cc]->tf[mm][fc] );
                        rr++ ;
                     };
# endif /* DSC_FCETEMP == 1 */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
               }; /* end if ( val.nj <= nn ) */
/*............................................................................*/
/* reflection cycle [ heat and fluids ]: */

               ( state->hctme ) += (( state->hcdt )/2.);

               cc = null; do
               {
/*............................................................................*/
# ifdef _Include       
                  if (( state->hcrstart ) <= ( state->hctme ))
                     hre[cc] = scathcr( cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                     hre[cc] = sctflow( cc );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# else /* _Include not defined */
                  if (( state->hcrstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hre[cc] = scathcr( state );
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hre[cc] = sctflow( state );
                  };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
/*............................................................................*/
/* evaluate nodes [ heat and fluids ]: */

               if ( val.nj <= nn )
               {
                  cc = null; do
                  {
                     rr = null;
                     while( rr < val.ntn[cc] )
                     {
                        mm = val.mtn[cc][rr];
                        val.tn[cc][rr] += ( hre[cc]->tn[mm] );
                        rr++ ;
                     };
/*............................................................................*/
# if DSC_FLDMDE != 0 
                     rr = null;
                     while( rr < val.nun[cc] )
                     {
                        mm = val.mun[cc][rr];
                        pp = val.cun[cc][rr];
                        val.un[cc][rr] += ( hre[cc]->un[mm][pp] );
                        rr++ ;
                     };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
                  } while(( ++cc ) < DSC_HCRMDE );
               }; /* end if ( val.nj <= nn ) */
/*............................................................................*/
/* connection cycle [ heat and fluids ]: */
	    
               ( state->hctme ) += (( state->hcdt )/2.);

               cc = null; do
               {
/*............................................................................*/
# ifdef _Include       
                  if (( state->hcrstart ) <= ( state->hctme ))
                     hci[cc] = cncthcr( cc );
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                     hci[cc] = conflow( cc );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# else /* _Include not defined */
                  if (( state->hcrstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hci[cc] = cncthcr( state );
                  };
/*............................................................................*/
# if DSC_FLDMDE != 0
                  if (( state->fldstart ) <= ( state->hctme ))
                  {
                     ( state->hclbl ) = cc;
                     hci[cc] = conflow( state );
                  };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* _Include not defined */
/*............................................................................*/
               } while(( ++cc ) < DSC_HCRMDE );
            }; /* end if ( nn <= val.nt ) */
         }; /* end while(( ++kk ) <= val.rc ); internal hcrr loop terminated */
/*............................................................................*/
/*
*//* end of internal loop [ heat current computation ]
*/
/*............................................................................*/
# endif /* DSC_INTLCE == 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         ind = values( -ONE );        /* store values                         */
/*..................................*/
/* display process state [ running cursor or percent value ]: */

         ( dsp->state ) = nn;
         dsplay( dsp );

         nn++ ;
      }; /* while( nn <= val.n ) [ val.n iterations terminated ] */
/*............................................................................*/
/*
*//* end of external loop - job terminated
*/
# endif /* DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
   };
/*............................................................................*/
/* close the job: */
/*............................................................................*/
   ind = values( -TWO );         /* close evaluation file                     */
/*.............................*/
# if DSC_STOPTIME == 1
   nseconds = time( timer );
   strncpy( ctmptr, ctime(&nseconds), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

   TIMEFORM( tmestr, ctmptr );
# endif
/*............................................................................*/
# if DSC_HCRMDE != 0 /* Fluid flow processes only: */
/*............................................................................*/
# if DSC_FLDMDE != 0 /* close flow parameters log files */

   fprintf(( state->dscstr ), "DANSE-" );
   fprintf(( state->dscstr ), "%s ---", DSC_RELEASE );
   kk = null; do /* fill line up with '-' to 78 characters */
   {
      fprintf(( state->dscstr ), "----------" );
   } while(( ++kk ) < SIX );
   fprintf(( state->dscstr ), "\nDSC process no %s terminated:\n", syslbl );

# if DSC_STOPTIME == 1
   fprintf(( state->dscstr ), "%s", tmestr );
# endif

   fclose( state->dscstr );
   
   fprintf(( state->prsstr ), "DANSE-" );
   fprintf(( state->prsstr ), "%s ------", DSC_RELEASE );
   kk = null; do /* fill line up with '-' to 78 characters */
   {
      fprintf(( state->prsstr ), "----------" );
   } while(( ++kk ) < SIX );
   fprintf(( state->prsstr ), "\nDSC process no %s terminated:\n", syslbl );

# if DSC_STOPTIME == 1
   fprintf(( state->prsstr ), "%s", tmestr );
# endif

   fclose( state->prsstr );
/*............................................................................*/
# if SLV_TMPDSP == 1
   
   ( state->tmpfst ) = fopen(( state->tmpdsp ), "r+" );
   fseek(( state->tmpfst ), ( state->tmpef ), SEEK_SET );

   fprintf(( state->tmpfst ), "----------" );
   kk = null; do
   {
      fprintf(( state->tmpfst ), "----------" );
   } while(( ++kk ) < SEVEN );
   fprintf(( state->tmpfst ), "\nDSC process no %s terminated:\n", syslbl );

# if DSC_STOPTIME == 1
   fprintf(( state->tmpfst ), "%s", tmestr );
# endif

   fclose( state->tmpfst );

# endif /* SLV_TMPDSP == 1 */
/*............................................................................*/
# if SLV_PRSDSP == 1

   ( state->prsfst ) = fopen(( state->prsdsp ), "r+" );
   fseek(( state->prsfst ), ( state->prsef ), SEEK_SET );

   fprintf(( state->prsfst ), "----------" );

   kk = null; do
   {
      fprintf(( state->prsfst ), "----------" );
   } while(( ++kk ) < SEVEN );

   fprintf(( state->prsfst ), "\nDSC process no %s terminated:\n", syslbl );

# if DSC_STOPTIME == 1
   fprintf(( state->prsfst ), "%s", tmestr );
# endif

   fclose( state->prsfst );

# endif /* SLV_PRSDSP == 1 */
/*............................................................................*/
# if CMPRSSBL != 0
# if SLV_DNSDSP == 1

   ( state->dnsfst ) = fopen(( state->dnsdsp ), "r+" );
   fseek(( state->dnsfst ), ( state->dnsef ), SEEK_SET );

   fprintf(( state->dnsfst ), "----------" );

   kk = null; do
   {
      fprintf(( state->dnsfst ), "----------" );
   } while(( ++kk ) < SEVEN );

   fprintf(( state->dnsfst ), "\nDSC process no %s terminated:\n", syslbl );

# if DSC_STOPTIME == 1
   fprintf(( state->dnsfst ), "%s", tmestr );
# endif

   fclose( state->dnsfst );

# endif /* SLV_DNSDSP == 1 */
# endif /* CMPRSSBL != 0 */
/*............................................................................*/
# if SLV_FLWDSP == 1

   ( state->velfst ) = fopen(( state->veldsp ), "r+" );
   fseek(( state->velfst ), ( state->velef ), SEEK_SET );

   fprintf(( state->velfst ), "------" );

   kk = null; do
   {
      fprintf(( state->velfst ), "----------" );
   } while(( ++kk ) < FIVE );

   fprintf(( state->velfst ), "\nDSC process no %s terminated:\n", syslbl );

# if DSC_STOPTIME == 1
   fprintf(( state->velfst ), "%s", tmestr );
# endif

   fclose( state->velfst );

# endif /* SLV_FLWDSP != 0 */
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   strcpy(( dsp->messge ), "DSC process no " );
   strcat(( dsp->messge ), syslbl );
   strncat(( dsp->messge ), " terminated", 11 );

# if DSC_STOPTIME == 1
   strncat(( dsp->messge ), " at ", 4 );
   strncat(( dsp->messge ), tmestr, 8 );
# endif

   if( *( state->opt ) == null )
   {
      strncat(( dsp->messge ), "\n", 2 );

      ( dsp->option ) = null;
      dsplay( dsp ); /* clear display */

# if SLV_ALERT == 1
   fprintf( display, "\a" ); /* produce audible signal */
# endif                      /* [ only while running in shell ] */
   }
   else
   {
      fprintf( display, "\n %s", ( dsp->messge ));
      ( state->logps ) = ftell( display );
      fclose( display );
   };

   return ONE;
}   
/*============================================================================*/
# undef POLAK_RIBIERE__ 
# undef FLETCHER_REEVES 
/************************* end of function solvdrv(*) *************************/
