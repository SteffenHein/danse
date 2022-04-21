/* [ file: solver.c ] */
# define PROGRAM "solver"
/*******************************************************************************
*                                                                              *
*   ANSI C program solver.c                                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   ELectromagnetic Field Evaluation program                                   *
*   based on the time and frequency domain DSC method                          *
*   of numerical approximation to Maxwell's equations                          *
*                                                                              *
*   [ Version supporting non-orthogonal mesh cell and                          *
*   gyrotropic material characteristics ]                                      *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 22, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
# include <time.h>	/* cf. time(*), ctime(*), asctime(*), localtime(*)    */
# include <unistd.h>	/* UNIX type systems standard header, cf. sysconf(*)  */
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
/* Edit and customize this header for SOLVER configuration: */
# include "./SOLVER.CONF"
/*----------------------------------------------------------------------------*/
/* Don't mix different releases: */
# ifndef SLV_RELEASE
   # define SLV_RELEASE DSC_RELEASE
# endif
# ifndef SLV_UPDATE
   # define SLV_UPDATE DSC_UPDATE
# endif
/*----------------------------------------------------------------------------*/
# ifndef USE_NCURSES
   # define USE_NCURSES 1
# endif
/*----------------------------------------------------------------------------*/
/* Fix these constants, e.g.[ sample excerpt of SOLVER.CONF ]: */
/* [ */
/* # define NODES 50000 *//* maximum number of mesh cells                     */
/* # define NRSMX NODES *//* maximum number of [ different ] s-matrices       */
/* # define GENDS     1 *//* maximum number of gyroelectric nodes             */
/* # define GESMX     1 *//* maximum number of [ diff.] gyroelct. s-matrices  */
/* # define GMNDS     1 *//* maximum number of gyromagnetic nodes             */
/* # define GMSMX     1 *//* maximum number of [ diff.] gyromagn. s-matrices  */
/* # define BNDAP  2000 *//* maximum number of aperiodic boundary faces       */
/* # define BNDPR     1 *//* maximum number of periodic boundary faces        */
/* # define EXCEP  1000 *//* maximum number of excited E ports                */
/* # define EXCHP  1000 *//* maximum number of excited H ports                */
/* # define EVLEP   100 *//* maximum number of evaluated E ports              */
/* # define EVLEN     1 *//* maximum number of evaluated E nodes              */
/* # define EVLHP   100 *//* maximum number of evaluated H ports              */
/* # define EVLHN     1 *//* maximum number of evaluated H nodes              */
/*............................................................................*/
/* define additional currents, if needed [better defined in CONFIG.H, though] */
/* # ifndef DSC_HCRMDE                                                       */
/*    # define DSC_HCRMDE 1                                                  */
/* # endif                                                                    */
/*............................................................................*/
/* # if DSC_HCRMDE != 0                                                      */
/*    # define HCNDS NODES                                                    */
/*    # define BNDHC     1 *//* maximum number of heat current boundary faces */
/*    # define BNDTF BNDAP *//* maximum number of temperature boundary faces  */
/*    # define BNDTN     1 *//* maximum number of imposed node temperatures   */
/*    # define EXCHC     1 *//* maximum number of excited current nodes       */
/*    # define EXCTF     1 *//* maximum number of excited temperature faces   */
/*    # define EXCTN     1 *//* maximum number of excited temperature nodes   */
/*    # define EVLHC     1 *//* maximum number of evaluated heat current faces*/
/*    # define EVLTF     1 *//* maximum number of evaluated temperature faces */
/*    # define EVLTN     1 *//* maximum number of evaluated temperature nodes */
/*............................................................................*/
/*    # if DSC_FLDMDE != 0 */
/*       # define further constants [ cf. SOLVER.CONF ] */
/*    # endif */
/* # endif ] */
/*----------------------------------------------------------------------------*/
/* This may yet be defined in the general configuration header CONFIG.H:      */
# ifndef SYS_TIMES
   # define SYS_TIMES 0 /* 2: use function times(*)                           */
                        /* 1: use function difftime(*) [ ANSI C standard ]    */
                        /* 0: use function clock(*)    [ ANSI C standard ]    */
# endif
/*----------------------------------------------------------------------------*/
/* system time function header, defining times(*), e.g.: */
# if SYS_TIMES == 2
   # include <sys/times.h>
# endif
/*----------------------------------------------------------------------------*/
/* the solver text console header [ typedefs of function txcnsl(*) etc.]: */
# include "../tools/txctyp.h"
/*----------------------------------------------------------------------------*/
# ifndef DSC_DOMAIN
   # define DSC_DOMAIN 0 /* 1: time domain,                                   */
		/* 2: frequency domain [ steady state relax. algorithm ]      */
		/* 0: time&frequency domains [ steady state relax. alg.]      */
# endif
/*----------------------------------------------------------------------------*/
/* operational constants: [ may be modified as suitable - without computatio- */
/*                          nal impact ]                                      */
# define DSC_STOPTIME 1 /* 1: print start and stopping time in solvdrv(*)     */

# ifndef DSC_DAEMON
   # define DSC_DAEMON 1 /* 1: allows solver to be run in daemon option "-d"  */
# endif
# if DSC_DAEMON == 1
/* where the daemon finds the job instructions [ dsc.batch ] and DSC system   */
/* files [ dsc.top<N>,..., dsc.val<N> ]:                                      */
   # ifndef DSC_WRK
      # define DSC_WRK "/users/dsc/wrk/"
   # endif
   # define DSC_LNGNAMES 1 /* must be 1 with DSC_DAEMON set to 1 !!! */
   # ifndef DSC_LOG
      # define DSC_LOG "/users/dsc/log/dsc.log"
   # endif
   # ifndef DSC_PID
      # define DSC_PID "/users/dsc/run/dsc.pid"
   # endif
/* the daemon initializer function dminit(*) requires this header: */
   # include "../math/unxsys.h"
# else
   # define DSC_LNGNAMES 0 /* 0/1: allow 20/80 char file names [with care] */
# endif

# define DSC_DISP 1 /* 1: various display functions in toplgy(*) etc.  */
# define DSC_INIT 2 /* explicit parameter initialization:              */
                    /* 1: before, 2: after -- every job executation    */
# define DSC_BATCH 'b' /* 'b': batch file program mode [ default ]     */
                       /*  n : single job program mode for input files */
                       /*      indexed n, 0 <= n < DSC_JOBS            */
# if DSC_BATCH == 'b'
   # define DSC_BTCHFLE "batch" /*extension for [jobs definition] batch files*/
   # define DSC_PRTCTFL 1 /* 1: protect DSC system files from removal if */
# endif                /*    job can't be executed [ for some reason]   */

# define DSC_JOBS  100 /* maximum number of jobs in batch mode */
                       /* [ don't fix more than 100 ]          */
# define DSC_MXJBL 999 /* maximum job label                    */
                       /* [ don't fix more than 999 ]          */
/*----------------------------------------------------------------------------*/
/* DSC system input and output file names, pre- and suffixes: >-------------> */

# define DSC_PRFX    "dsc." /* common prefix of SOLVER.C operation files    */

# define TOPOLOGY_FILE   "top"  /* topology file                   [ suffix ] */
# define S_MATRIX_FILE   "smx"  /* s-matrix file                   [ "      ] */
# define BOUNDARY_FILE   "bnd"  /* boundary conditions file        [ "      ] */
# define EXCITATION_FILE "exc"  /* excitation                      [ "      ] */
# define EVALUATION_FILE "val"  /* and evaluation file             [ "      ] */
# define STDOUT_FILE "out"
/*----------------------------------------------------------------------------*/
/* system oriented declarations, allusions and function prototypes: >-------> */
                                        
/* system function prototypes: >--------------------------------------------> */

time_t time( time_t *timer );                
char *ctime( const time_t *timer );          

struct tm *localtime( const time_t *timer ); 
char *asctime( const struct tm *lt_ptr );

# if SYS_TIMES == 2
   clock_t times( struct tms *cl_ptr );
# elif SYS_TIMES == 1
   double difftime( time_t time1, time_t time0 );
# else
   clock_t clock( void );
# endif

long sysconf( int s_name ); 

# define BUFFER 1024
int setvbuf( FILE *display, char *buff, int unbuff, size_t bufsize ); 

/*-- time function structures: ---------------------------------------------->*/

static struct tm loct = {null};
static struct tm *lt_ptr= &loct;

# if SYS_TIMES == 2
   struct tms cpu = {null};
   struct tms *cl_ptr = &cpu;
# endif

static const char *wday[] = { "Sun","Mon","Tue","Wed","Thu","Fri","Sat" };

/*-- global variables:------------------------------------------------------->*/

static int unbuff     = _IONBF;
static size_t bufsize =  null;

/*-- user oriented declarations, allusions and function prototypes: --------->*/
   static char 
      syslbl[VSS_SIZE] = {null};
/*----------------------------------------------------------------------------*/
   static FILE 
     *display = NULL;
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* enhance 'my_terminal' window: */

   # include <ncurses.h>
   # include <curses.h>
   # include <term.h> /* the terminal type header */

   static char *term;        /* terminal type string */ 

   # define CLSCREEN /* clear screen */ \
   { \
     fprintf( display, "%s", tgetstr( "cl", null )); \
   }

   # define PRBLDCLR(a) /* bold clear output */ \
   { \
      fprintf( display, "%s%s", tgetstr( "md", null ), (a)); \
   }

   # define PRINVERS(a) /* print inverse ( black on white ) */ \
   { \
     fprintf( display, "%s%s", tgetstr( "mr", null ), (a)); \
   }

   # define PRNORMAL(a) /* back to normal output */ \
   {\
     fprintf( display, "%s%s", tgetstr( "me", null ), (a)); \
   }
# else
   # define CLSCREEN { ;\
   }

   # define PRBLDCLR(a) {\
     fprintf( display, "%s", (a));\
   }

   # define PRINVERS(a) {\
     fprintf( display, "%s", (a));\
   }

   # define PRNORMAL(a) {\
     fprintf( display, "%s", (a));\
   }
# endif
/*----------------------------------------------------------------------------*/
static short joblbl[DSC_JOBS] = {null}; /* job indices array */
/*----------------------------------------------------------------------------*/
/* DANSE library  subroutines and functions: */
# include "./solvdrv.h"
/*============================================================================*/

int main( int argc, char **argv )
{
                          /* -- system oriented declarations and allusions -->*/
   static short
      hrs = null,
      min = null,
      sec = null;

   static time_t
      nseconds = null,
      *timer = null;

   static double
      scs = ZERO;

# if SYS_TIMES == 2
   static clock_t
      ticks = null,
      job_time = null,
      usr_time = null;
# elif SYS_TIMES == 1
   static time_t 
      time0 = null,
      time1 = null,
      time2 = null,
      time3 = null;

   static double
      usr_time = ZERO,
      job_time = ZERO;
# else
   static clock_t
      job_time = null,
      usr_time = null;
# endif

/*-- system function prototypes: -------------------------------------------->*/

   time_t time( time_t *timer );
   char *ctime( const time_t *timer );

   struct tm *localtime( const time_t *timer );
   char *asctime ( const struct tm *lt_ptr );

   long sysconf( int s_name );

   clock_t clock( void );

   int setvbuf( FILE *display, char *buff, int mode, size_t bfsize );

   long strtol( const char *ptr, char **endp, int n );
   double strtod( const char *ptr, char **endp );

# if USE_NCURSES == 1
   int setupterm( NCURSES_CONST char *term, int fildes, int *errret );
# endif

# ifndef _CCBUG
   char 
     *strcat( char *pointer1, const char *pointer2 ),
     *strncat( char *pointer1, const char *pointer2, size_t n );
# endif

/*-- user oriented declarations and allusions ---->*/

   static FILE 
     *batchfle,
     *checkfle;

/* allusions: */
/*
   extern struct solverstat solver;

   extern struct topology top;
   extern struct boundary bnd;

   extern char joblbl[];
   extern char toplbl[];
   extern char smxlbl[];
   extern char bndlbl[];
   extern char exclbl[];
*/
   static char
      ptr[STS_SIZE] = {null};

   static char
     *tmeptr,
     *bchfle,
     *fleptr,
     *status,
    **endp = null;

   static const char
     *astrx = "***",
     *btchp = DSC_BTCHFLE,
     *rmfle = "remove_files",
     *rtfle = "retain_files";

# if DSC_PRTCTFL == 1
/* file protection arrays: they store all DSC model file names */
/* that pertain to unterminated jobs, thus protecting these from */
/* premature removal [ option solver.rmfle = TWO ]. */

   static char 
      topfls[FIVE*DSC_JOBS+ONE] = {null},
      smxfls[FIVE*DSC_JOBS+ONE] = {null},
      bndfls[FIVE*DSC_JOBS+ONE] = {null},
      excfls[FIVE*DSC_JOBS+ONE] = {null};
# endif

   static short 
      ii = null, 
      jj = null,
      kk = null,
     ind = null;

   static long 
      mm = null;

   static signed char 
      hh = null,
      ll = null,
      pp = null,
      qq = null;

   static TXCNSL
     *csp;

# if DSC_DAEMON == 1
   static DMINIT
     *dpt;
# endif

/*-- math. function prototypes: --------------------------------------------->*/

   double fabs ( double x );
   int abs( int n );

/*-- user defined function prototypes: -------------------------------------->*/

   short toplgy( const short nn );
   short smtrix( const short nn );
   short boundr( const short nn );
   short excite( const short nn );
   short values( const short nn );
   short trnsnt( const double ttime );         /* excitation 'swing; function */
   short clearv( long, long, long, long, long, long, long, 
                 long, long, long, long, long, long, long );
   char
     *lotos ( long lngint, char length );

   TXCNSL
     *txcnsl( TXCNSL *csp );
/*----------------------------------------------------------------------------*/
/* initialize display on stdout */

   display = stdout;
/*............................................................................*/
/* the selected clock may require: */

# if SYS_TIMES == 2
   ticks = sysconf( _SC_CLK_TCK );      /* ticks = clocks per second          */
# endif
/*............................................................................*/
/* initial memory allocations: */

   status  = ( char * ) calloc( TWO, ONE );
   tmeptr  = ( char * ) calloc( STS_SIZE, ONE );
   bchfle  = ( char * ) calloc( STS_SIZE, ONE );
   fleptr  = ( char * ) calloc( STS_SIZE, ONE );

   solver.opt = ( char * ) calloc( SHS_SIZE, ONE );

# if DSC_LNGNAMES == 1
   solver.top = ( char * ) calloc( STS_SIZE, ONE );
   solver.smx = ( char * ) calloc( STS_SIZE, ONE );
   solver.bnd = ( char * ) calloc( STS_SIZE, ONE );
   solver.exc = ( char * ) calloc( STS_SIZE, ONE );
   solver.val = ( char * ) calloc( STS_SIZE, ONE );
   solver.prfx = ( char * ) calloc( STS_SIZE, ONE );
# else
   solver.top = ( char * ) calloc( SHS_SIZE, ONE );
   solver.smx = ( char * ) calloc( SHS_SIZE, ONE );
   solver.bnd = ( char * ) calloc( SHS_SIZE, ONE );
   solver.exc = ( char * ) calloc( SHS_SIZE, ONE );
   solver.val = ( char * ) calloc( SHS_SIZE, ONE );
   solver.prfx = ( char * ) calloc( VSS_SIZE, ONE );
# endif

   solver.logfle = ( char * ) calloc( SHS_SIZE, ONE );
   solver.fcterr = ( char * ) calloc( SHS_SIZE, ONE );
   solver.errmsg = ( char * ) calloc( LGS_SIZE, ONE );
/*............................................................................*/
/* initialize time/frequency domain indicator: */

   solver.dmn = 't'; /* default: 't'ime domain process */
/*............................................................................*/
/* initialize strings and job labels: */

# if DSC_LNGNAMES == 1
   strcpy ( solver.prfx, DSC_PRFX );
# else
   strncpy ( solver.prfx, DSC_PRFX, VSS_SIZE );
# endif

# if DSC_BATCH == 'b'
/* batch file operation mode */

   # if DSC_LNGNAMES == 1
      strcpy( bchfle, solver.prfx );
      strcat( bchfle, btchp );
   # else
      strncpy( bchfle, solver.prfx, VSS_SIZE );
      strncat( bchfle, btchp, SHS_SIZE - VSS_SIZE );
   # endif
# endif
/*............................................................................*/
/* initialize DSC system file indicators: */
/* [ natural labels; may be changed with batch file input ] */

   jj = null;
   while( jj < DSC_JOBS )
   {
      joblbl[jj] = jj; 
      toplbl[jj] = jj;
      smxlbl[jj] = jj;
      bndlbl[jj] = jj;
      exclbl[jj] = jj;
      jj++ ;
   };  
/*............................................................................*/
/* initialize [ clear ] parameter arrays: */

   ind = clearv( null, null, null, null, null, null, null,
                 null, null, null, null, null, null, null );
/*............................................................................*/
/* stop program start time: */

   nseconds = time( timer );
   tmeptr = ctime( &nseconds );
/* ['tmeptr = ctime(*)' same as 'tmeptr = asctime( localtime( &nseconds ));'] */
/*............................................................................*/
/* extract options solver.opt from the command line: */

   if( null < --argc )
   {
      do
      {
         strcpy( ptr, *++argv );
         
         if (( null == strcmp( ptr, "v" ))||
             ( null == strcmp( ptr, "-v" )))
         {
            fprintf( display, "\rsolver release %s;", SLV_RELEASE );
            fprintf( display, " %s\n", SLV_UPDATE );
            exit( EXIT_SUCCESS );
         }
         else if (( null == strcmp( ptr, "d" ))||
                  ( null == strcmp( ptr, "-d" ))) /* daemon mode */
         {
            strcat(( solver.opt+ONE ), "d" );
            continue;
         }
         else if (( null == strcmp( ptr, "b" ))||
                  ( null == strcmp( ptr, "-b" ))) /* background running mode */
         {
            strcat(( solver.opt+ONE ), "b" );
            continue;
         }
         else if (( null == strcmp( ptr, "m" ))||
                  ( null == strcmp( ptr, "-m" ))) /* menu = terminal mode */
         {
            strcat(( solver.opt+ONE ), "m" );
            continue;
         }
	 else
	 {
            fprintf( display, "\runknown program option solver %s\n", ptr );
            exit( EXIT_FAILURE );
         };
      } while(( null < ( --argc ))&&( strlen( solver.opt ) < SHS_SIZE ));
   }; /* end if ( null < argc ) */
/*............................................................................*/
/* evaluate the command line options solver.opt by priorities: */ 
/*............................................................................*/
# if DSC_DAEMON == 1

   if ( null != strchr(( solver.opt + ONE ), 'd' ))
   {
/* initialize structure dmn [ type DMINIT ]: */
/*............................................................................*/
      dpt = dminit( null );   /*                                           */
/*..........................*/
   
      strcpy ( solver.prfx, DSC_WRK );
      strcat ( solver.prfx, DSC_PRFX );

      strcpy( dpt->dmndir, DSC_WRK );
      strcpy( dpt->logfle, DSC_LOG );
      strcpy( dpt->pidfle, DSC_PID );

      ( dpt->uid ) = null; /* require root permissions to start the daemon */
/*............................................................................*/
      dpt = dminit( dpt );   /*                                               */
/*.........................*/
		
      strcpy( solver.logfle, ( dpt->logfle ));
      *solver.opt = TWO;
   }
   else if ( null != strchr(( solver.opt+ONE ), 'b' ))
   {
/* initialize daemon structure *dpt = &dmn [ of type DMINIT ]: */
/*............................................................................*/
      dpt = dminit( null );   /*                                              */
/*..........................*/

      strcpy ( solver.prfx, DSC_PRFX );

      strcpy(( dpt->dmndir ), "./" );
      strcpy(( dpt->logfle ), solver.prfx );
      strcat(( dpt->logfle ), "log" );
      strcpy(( dpt->pidfle ), solver.prfx );
      strcat(( dpt->pidfle ), "pid" );

      fprintf( display, "\rdsc-solver started in background mode," );
      fprintf( display, "\nterminal output is directed to file %s\n",
         ( dpt->logfle ));
      ( dpt->uid ) = ONE; /* any user may start the daemon */
/*............................................................................*/
      dpt = dminit( dpt );   /*                                               */
/*.........................*/
      strcpy( solver.logfle, ( dpt->logfle ));
      *solver.opt = ONE;
   }
   else
      strcpy ( solver.prfx, DSC_PRFX );

# else /* depending on allowed length of filenames: */

   # if DSC_LNGNAMES == 0
      strncpy ( solver.prfx, DSC_PRFX, VSS_SIZE );
   # else
      strcpy ( solver.prfx, DSC_PRFX );
   # endif

   if ( null != strchr(( solver.opt+ONE ), 'd' ))
   {
      strcpy( ptr, DSC_PRFX );
      strcat( ptr, "log" );

      fprintf( display, "\rdsc-solver not compiled for daemon operation," );
      fprintf( display, "\nstarted with terminal set to '%s'\n", ptr );

      strcpy( solver.logfle, ptr );
      *solver.opt = ONE;
   }
   else if ( null != strchr(( solver.opt+ONE ), 'b' ))
   {
      strcpy( ptr, DSC_PRFX );
      strcat( ptr, "log" );

      fprintf( display, "\rdsc-solver started in background mode," );
      fprintf( display, "\nterminal output is directed to '%s'\n", ptr );

      strcpy( solver.logfle, ptr );
      *solver.opt = ONE;
   };
# endif /* DSC_DAEMON ==... */
/*............................................................................*/

   if (( *solver.opt != TWO )&&( null != strchr(( solver.opt+ONE ), 'm' )))
      *solver.opt = null;

/*............................................................................*/
/* terminal mode [ keybord input ] */

   if ( *solver.opt == null ) /* initialize program in terminal mode */
   {
/*............................................................................*/
# if USE_NCURSES == 1
/* set the terminal type: */

# if defined ( _BSD )
      setupterm( "cons25", 1, (int *) 0 );
# elif defined ( _GNU_Linux )
      setupterm( "linux", 1, (int *) 0 );
# elif defined ( _Linux )
      setupterm( "linux", 1, (int *) 0 );
# endif
/*............................................................................*/
/* get the terminal info: */

      term = ( char *) getenv( "TERM" ); /* get the terminal */

      kk = tgetent( null, term );

      if( ONE != kk )
      {
         fprintf( stderr, "Error on getting the termcap info\n" );
         exit( EXIT_FAILURE );
      };
# endif
/*............................................................................*/
/* set buffer length to null: */

      ind = setvbuf( stdin, null, unbuff, bufsize );
      ind = setvbuf( stdout, null, unbuff, bufsize ); 
/*............................................................................*/
      csp = txcnsl( null );    /* initialize the text console                 */
/*...........................*/
      ( csp->clscr ) = 2; /* clear screen and scroll that number of lines */

      strcpy( csp->title, "DSC program " ); /* start message: */
      strcat( csp->title, PROGRAM );        /* DSC program ...*/
      strcat( csp->title, " started:\n " ); /* started:<time> */
      strncat( csp->title, tmeptr, 24 );
      
      strcpy( csp->cmmnt, "Welcome to DSC SOLVER !" );
      strcpy( csp->envmt, "DSC SOLVER" );
/*............................................................................*/
# if DSC_BATCH != 'b'
/*............................................................................*/
      csp = txcnsl( csp );    /* display the start message [ title ]          */
/*..........................*/
# elif DSC_BATCH == 'b' /* batch file mode: */

      ( csp->dfopt ) = ONE; /* the initial menu option */
      strcpy( csp->cnfrm, "Nothing done! Do you really want to quit ?" );

     menu:  /* mode */

      ( csp->clscr ) = 2; /* clear screen; scroll that number of lines */
      ( csp->items ) = 2;
      strcpy( csp->tasks, "Enter job instructions ..." );
      strcpy( csp->mline[1], "* from batch file [ 'dsc.batch' ]" );
      strcpy( csp->mline[2], "* on keyboard" );
      strcpy( csp->escpe, "End of program / escape:" );
/*............................................................................*/
      csp = txcnsl( csp );    /* build the start menu                         */
/*..........................*/
      ind = ( csp->option );

      if( ind == TWO ) /* create new batch file: */
      {
         strcpy( csp->rqfrm, "brackets" );

        number_of_jobs:

         fprintf( display, "\n" );
         strcpy( ptr, "Please enter number of jobs " );
         strcat( ptr, "[ 0 <= n <=" );
         strcat( ptr, lotos( DSC_JOBS, null ));
         strcat( ptr, " ]" );
         strcpy( csp->rqlng, ptr );
         ( csp->dflng ) = ONE;
/*............................................................................*/
         csp = txcnsl( csp );        /*                                       */
/*.................................*/
         solver.nbrjbs = csp->inlng;

         if ( solver.nbrjbs <= null )
         {
            solver.nbrjbs = null;
            ind = null;
         }
         else if ( DSC_JOBS < solver.nbrjbs )
         {
            fprintf( display, " Too many jobs !!!  [ maximum number is"
                              " %3d ]\n", DSC_JOBS );
            goto number_of_jobs;
         }
         else
         {
            fprintf( display, "\n" );
            strcpy( csp->rqfrm, "brackets" );
            strcpy( csp->rqstr, "Relabel DSC system files ? [y/n]" );
            strcpy( csp->dfstr, "n" ); /* default */
/*............................................................................*/
            csp = txcnsl( csp );       /*                                     */
/*...................................*/
            strcpy( ptr, csp->instr );

            if (( *ptr == 'y' )||( *ptr == 'Y' ))
            {
              relabel:

               for ( kk=null; kk<solver.nbrjbs; kk++ )
               {
                  fprintf( display, "\n" );
                  strcpy( csp->rqfrm, "brackets" );

                  if ( kk == null )
                     ( csp->dflng ) = null;
                  else
                     ( csp->dflng) = joblbl[kk-ONE] + ONE;

                  strcpy( ptr, "Enter " ); 
                  strcat( ptr, lotos(( kk+ONE ), null ));
                  strcat( ptr, ". job index" );
                  strcpy( csp->rqlng, ptr );
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  joblbl[kk] = csp->inlng;

		  strcpy( ptr, "Enter topology file index of job no." ); 
		  strcat( ptr, lotos( joblbl[kk], null ));
                  strcpy( csp->rqlng, ptr );
                  ( csp->dflng ) = csp->inlng;
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  toplbl[kk] = csp->inlng;

                  strcpy( ptr, "Enter s-matrix file index of job no." ); 
		  strcat( ptr, lotos( joblbl[kk], null ));
                  strcpy( csp->rqlng, ptr );
                  ( csp->dflng ) = csp->inlng;
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  smxlbl[kk]= csp->inlng;

                  strcpy( ptr, "Enter boundary file index of job no." ); 
                  strcat( ptr, lotos( joblbl[kk], null ));
                  strcpy( csp->rqlng, ptr );
                  ( csp->dflng ) = csp->inlng;
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  bndlbl[kk]= csp->inlng;

                  strcpy( ptr, "Enter excitation file index of job no." ); 
                  strcat( ptr, lotos( joblbl[kk], null ));
                  strcpy( csp->rqlng, ptr );
                  ( csp->dflng ) = csp->inlng;
/*............................................................................*/
                  csp = txcnsl( csp );        /*                              */
/*..........................................*/
                  exclbl[kk]= csp->inlng;
               };

               fprintf( display, "\n" );
               strcpy( csp->rqstr, "Input correct ? [y/n] " );
               strcat( csp->dfstr, "y" ); /* default */
/*............................................................................*/
               csp = txcnsl( csp );       /*                                  */
/*......................................*/
               strcpy( ptr, csp->instr );

               if (( *ptr == 'n' )||( *ptr == 'N' ))
                  goto relabel;
            }
            else /* natural labels */
            {
               kk = DSC_MXJBL - solver.nbrjbs;

              first_job:

               strcpy( ptr, "Please enter index of first job " );
               strcat( ptr, "[ 0 <= i <=" );
               strcat( ptr, lotos( kk, null ));
               strcat( ptr, " ]" );
               strcpy( csp->rqlng, ptr );
               ( csp->dflng ) = null;
/*............................................................................*/
               csp = txcnsl( csp );        /*                                 */
/*.......................................*/
               jj = csp->inlng;

               if (( jj < null )||( kk < jj ))
               {
                  fprintf( display, " illegal job index !!!"
                                    " [ 0 <= i <= %3d ]", kk );
                  goto first_job;
               };

               for ( kk=null; kk<solver.nbrjbs; kk++ )
               {
                  joblbl[kk] = jj;
                  toplbl[kk] = jj;
                  smxlbl[kk] = jj;
                  bndlbl[kk] = jj;
                  exclbl[kk] = jj;
                  jj++;
               };
            };

            fprintf( display, "\n" );
            strcpy( csp->rqfrm, "brackets" );
            strcat( csp->rqstr, "Remove DSC system files on program " );
            strcat( csp->rqstr, "execution ? [y/n]" );
            strcpy( csp->dfstr, "n" ); /* default */
/*............................................................................*/
            csp = txcnsl( csp );       /*                                     */
/*...................................*/
            strcpy( ptr, csp->instr );

            if (( *ptr == 'y' )||( *ptr == 'Y' ))
            {
/*............................................................................*/
# if DSC_PRTCTFL == 1
               solver.rmfle = TWO;
# else
               solver.rmfle = ONE;
# endif
/*............................................................................*/
            }
            else
               solver.rmfle = null;

            batchfle = fopen( bchfle, "w" );

            if ( null == batchfle )
            {
               fprintf( stderr,
                  "\n\n Error message from program %s :", PROGRAM );
               fprintf( stderr,
                  "\n Can't open job instructions file %s !", bchfle );
               fprintf( stderr,
                  "\n [ Check file location and/or permissions.]" );
               fprintf( stderr,
                  "\n Program %s stopped.\n", PROGRAM );

               exit( EXIT_FAILURE );
            };

            fprintf( batchfle, "%s", "DSC_batchfile" );

            if ( solver.rmfle != null )
               fprintf( batchfle, "\n%s", rmfle );
            else
               fprintf( batchfle, "\n%s", rtfle );

            fprintf( batchfle, "\n%s %3d", "number_of_jobs:", solver.nbrjbs );
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

            fprintf( batchfle," %s", "state" );
            fprintf( batchfle," %s", "CPU_time/err" );

            for ( kk=null; kk<solver.nbrjbs; kk++ )
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

            nseconds = time( timer );
            tmeptr   = ctime( &nseconds );

            fprintf( batchfle, "\n\n%s%s%s\n%.24s", "DSC batch file '",
                                 bchfle, "' created:   ", tmeptr );
            fprintf( batchfle, "\nEOF" );

            fclose( batchfle );

            printf( "\n %s%s%s\n %.24s\n", "DSC batch file '",
                    bchfle, "' created:", tmeptr );

         };/* end if null < solver.nbrjbs [ = number of jobs ] */
      };/* end if ind == TWO: create new batch file */
   } /* end if ( solver.opt ) == null */
/*............................................................................*/
   else /* if(*solver.opt != 0): non-terminal [daemon or background] operation: */
   {
      freopen( solver.logfle, "a+", display );
      setvbuf( display, null, unbuff, bufsize );

      fprintf( display, "\n DSC program " );        /* start message: */
      fprintf( display, "%s ", PROGRAM );           /* DSC program ...*/
      fprintf( display, "started: %.24s", tmeptr ); /* started:<time> */
                                                    
      ind = ONE;
   };
/*............................................................................*/
/* enter job instructions from file DSC_PRFX.batch: */

  job_instruct:

   if ( ind == THREE )
      fclose( batchfle );
     
   if ( null < ind )
   {    
      batchfle = fopen( bchfle, "r+" );

      if ( batchfle == null )
      {
         if( *solver.opt == null )
	 {
            CLSCREEN; /* clear screen */

            fprintf( display, "\n ======================================="
                              "=======================================" );
            fprintf( display, "\n Error on opening DSC batch file"
                    " '%s' !", bchfle );
            fprintf( display, "\n Check directory, if file '%s' is"
                    " present.", bchfle ); 
            fprintf( display, "\n [ Otherwise you can directly enter job"
                    " instructions on keyboard, cf. option (2)" );
            fprintf( display, "\n - batch file will be simultaneously"
                              " created.]" );
            fprintf( display, "\n ======================================="
                              "=======================================" );
            PRBLDCLR( "\r" );
            fprintf( display, "\n %*s", 78, "DSC SOLVER" );
            PRNORMAL( "\r" );

            ( csp->dfopt ) = 2; /* the next menu default option */
	 
            goto menu; /* back to menu */
         }
	 else /* background or daemon mode */
	 { /* can't open or read batch file; stop program */
            fprintf( display,
               "\n\n Error message from program %s :", PROGRAM );
            fprintf( display,
               "\n Can't open job instructions file %s !", bchfle );
            fprintf( display,
               "\n [ Check file location and/or permissions.]" );
            fprintf( display,
               "\n Program %s stopped.\n", PROGRAM );
/*............................................................................*/
# if DSC_DAEMON == 1
            remove( dpt->pidfle ); /* process identity number */
	                           /* in daemon or background operation mode */
# endif
/*............................................................................*/
            exit( EXIT_FAILURE );
         };
      } 
      else if ( ind == ONE )
         fprintf( display, "\n opened: %s'%s'", "DSC batch file ", bchfle );
      
      fscanf( batchfle, "%s", ptr ); /* comment */

/* the following depends on old/new batch file format: */

      fscanf( batchfle, "%s", ptr ); /* new format: 'remove/retain_files' */
                                     /* old format: 'number_of_jobs' */
                                     /* in the latter case ...: */

      if ( null == strncmp( ptr, "number_of_jobs", SIX )) /* [ old format ] */
      {
         fscanf( batchfle, "%s", ptr ); /* job number */
         solver.nbrjbs = strtol( ptr, endp, DEC );

         fscanf( batchfle, "%s", ptr ); /* string: 'first_job_index' */
         fscanf( batchfle, "%s", ptr ); /* index */
         fscanf( batchfle, "%s", ptr ); /* string: 'remove/retain_files' */
      };

      if ( null == strncmp( rmfle, ptr, SIX )) /* case: remove_files_ */
      {
/*............................................................................*/
# if DSC_PRTCTFL == 1
         solver.rmfle = TWO;

	 jj = null;
	 while( jj <= FIVE*DSC_JOBS )
	 {
            topfls[jj] = null; 
            smxfls[jj] = null; 
            bndfls[jj] = null; 
            excfls[jj] = null; 
           jj++ ;
         };
# else
         solver.rmfle = ONE;
# endif
/*............................................................................*/
      }
      else /* retain files */
      {
         solver.rmfle = null;
      };
      
      fscanf( batchfle, "%s", ptr ); /* new format: 'number_of_jobs' */
                                     /* old format: 'job' */
                                      
      if ( null == strncmp( ptr, "number_of_jobs", SIX )) /* [ new format ] */
      {
         fscanf( batchfle, "%s", ptr ); /* number */
         solver.nbrjbs = strtol( ptr, endp, DEC );
         fscanf( batchfle,"%s", ptr );  /* string: 'job' */
      };

      for ( kk=ONE; kk<EIGHT; kk++ )
      {
         fscanf( batchfle,"%s", ptr ); /* enter rest of line [ past 'job'] */
      };
         
      if ( DSC_JOBS < solver.nbrjbs )
      {
         fprintf( display, "\n Too many jobs defined in batch file !!! " );
         fprintf( display, "\n [ reduced to maximum number %3d ]",
                                                                  DSC_JOBS );
         solver.nbrjbs = DSC_JOBS; /* maximum number of jobs */
      };

      solver.rstart = solver.nbrjbs;
      solver.stofs = null;

      ii = null; jj = null;
      while ( jj < solver.nbrjbs )
      {
         fscanf( batchfle, "%s", ptr );	/* jj_th job index */
         joblbl[jj] = strtol( ptr, endp, DEC );
         fscanf( batchfle, "%s", ptr ); /* pertinent topology file index */
         toplbl[jj] = strtol( ptr, endp, DEC );
         fscanf( batchfle, "%s", ptr ); /* pertinent s-matrix file index */
         smxlbl[jj] = strtol( ptr, endp, DEC );
         fscanf( batchfle, "%s", ptr ); /* pertinent boundary file index */
         bndlbl[jj] = strtol( ptr, endp, DEC );
         fscanf( batchfle, "%s", ptr ); /* pertinent excitation file idx */
         exclbl[jj] = strtol( ptr, endp, DEC );
         fscanf( batchfle, "%s", ptr ); /* evaluation file index = job index */

         mm = ftell( batchfle );

         fscanf( batchfle, "%s", status );  /* states: "--", "ok", "??", e.g. */
/*............................................................................*/
/* check presence of DSC system files if status != ok: */

         if ( null != strncmp( status, "ok", TWO ))
         {
# if DSC_LNGNAMES == 1
            strcpy( solver.top, DSC_PRFX );
            strcat( solver.top, TOPOLOGY_FILE );

            strcpy( solver.smx, DSC_PRFX );
            strcat( solver.smx, S_MATRIX_FILE );
            
            strcpy( solver.bnd, DSC_PRFX );
            strcat( solver.bnd, BOUNDARY_FILE );

            strcpy( solver.exc, DSC_PRFX );
            strcat( solver.exc, EXCITATION_FILE );

            strcpy( solver.val, DSC_PRFX );
            strcat( solver.val, EVALUATION_FILE );
# else
            strncpy( solver.top, DSC_PRFX, VSS_SIZE );
            strncat( solver.top, TOPOLOGY_FILE, ( SHS_SIZE - VSS_SIZE - TWO ));

            strncpy( solver.smx, DSC_PRFX, VSS_SIZE );
            strncat( solver.smx, S_MATRIX_FILE, ( SHS_SIZE - VSS_SIZE - TWO ));
            
            strncpy( solver.bnd, DSC_PRFX, VSS_SIZE );
            strncat( solver.bnd, BOUNDARY_FILE, ( SHS_SIZE - VSS_SIZE - TWO ));

            strncpy( solver.exc, DSC_PRFX, VSS_SIZE );
            strncat( solver.exc, EXCITATION_FILE, ( SHS_SIZE - VSS_SIZE - TWO ));

            strncpy( solver.val, DSC_PRFX, VSS_SIZE );
            strncat( solver.val, EVALUATION_FILE, ( SHS_SIZE - VSS_SIZE - TWO ));
# endif
            strcat( solver.top, lotos( toplbl[jj], null ));
            strcat( solver.smx, lotos( smxlbl[jj], null ));
            strcat( solver.bnd, lotos( bndlbl[jj], null ));
            strcat( solver.exc, lotos( exclbl[jj], null ));
            strcat( solver.val, lotos( joblbl[jj], null ));

            strcpy( fleptr, solver.top );
            checkfle = fopen( fleptr, "r" );

            if ( checkfle != null ) /* check OK: topology file present */
            {
               fclose ( checkfle );
               strcpy( fleptr, solver.smx );
               checkfle = fopen( fleptr, "r" );

               if ( checkfle != null ) /* check OK: s-matrix file present */
               {
                  fclose ( checkfle );
                  strcpy( fleptr, solver.bnd );
                  checkfle = fopen( fleptr, "r" );

                  if ( checkfle != null ) /* check OK: boundary file present */
                  {
                     fclose ( checkfle );
                     strcpy( fleptr, solver.exc );
                     checkfle = fopen( fleptr, "r" );

                     if ( checkfle != null ) /* check OK: excit. file present */
                     {
                        fclose ( checkfle );
                        strcpy( fleptr, solver.val );
                        checkfle = fopen( fleptr, "r" );
                     };
                  };
               };
            };

            if ( checkfle == null ) /* check FAILDED: absent file(s) */
            {   
/*............................................................................*/
# if DSC_PRTCTFL == 1

               if ( solver.rmfle == TWO ) /* protect files against removal */
               {
                  strcpy( ptr, lotos( toplbl[jj], null ));
                  if (( null == strstr( topfls, ptr ))&&
                      ( strlen( topfls ) < ( FOUR*DSC_JOBS )))
                  {
                     strcat( topfls, ptr );
                     strncat( topfls, astrx, ONE );
                  };

                  strcpy( ptr, lotos( smxlbl[jj], null ));
                  if (( null == strstr( smxfls, ptr ))&&
                      ( strlen( smxfls ) < ( FOUR*DSC_JOBS )))
                  {
                     strcat( smxfls, ptr );
                     strncat( smxfls, astrx, ONE );
                  };

                  strcpy( ptr, lotos( bndlbl[jj], null ));
                  if (( null == strstr( bndfls, ptr ))&&
                      ( strlen( bndfls ) < ( FOUR*DSC_JOBS )))
                  {
                     strcat( bndfls, ptr );
                     strncat( bndfls, astrx, ONE );
                  };

                  strcpy( ptr, lotos( exclbl[jj], null ));
                  if (( null == strstr( excfls, ptr ))&&
                      ( strlen( excfls ) < ( FOUR*DSC_JOBS )))
                  {
                     strcat( excfls, ptr );
                     strncat( excfls, astrx, ONE );
                  };
               };
# endif /* DSC_PRTCTFL == 1 */
/*............................................................................*/
               strncpy( solver.fcterr, "solver.c", SHS_SIZE );
               strncpy( solver.errmsg, fleptr, SHS_SIZE );
               strncat( solver.errmsg, "...:_unable_to_open_file_",
                                     LGS_SIZE - SHS_SIZE );
               fseek( batchfle, mm, SEEK_SET );

               fprintf( batchfle," %6s", "??" ); /* write "??" state */
               fprintf( batchfle," %.9s%.3s\n", solver.errmsg, "..?" );

               fseek( batchfle, mm, SEEK_SET );

               fscanf( batchfle, "%s", status );
               
               if ( ii == null )
               {
                  ii = ONE;
               };

               fprintf( display, "\n\n Job no.%d: Can't open file(s)"
                  " %.9s ... ! ", joblbl[jj], fleptr );
               fprintf( display, "[ overriding.] " );
/*............................................................................*/
               values( jj );              /* write error message              */
/*......................................*//* into evaluation file elf.val[]   */
            }
            else if ( checkfle != null )
            {
               fclose( checkfle );

               fseek( batchfle, mm, SEEK_SET );

               fprintf( batchfle, " %6s", "--" );      /* write '--' state */
               fprintf( batchfle, " %12s", "------------" );

               fseek( batchfle, mm, SEEK_SET );

               fscanf( batchfle, "%s", ptr );
               
               if ( solver.stofs == null )
               {
                  solver.stofs = mm;
                  solver.rstart = jj;
               };
            };
         };
       
         fscanf( batchfle, "%s", ptr );   /* dummy: CPU_time or err messg */

         jj++;
      };

      solver.endfle = ftell( batchfle );
/*...........................................................................*/
/* all jobs terminated: */

      if (( solver.stofs == null )\
        &&( ii == null ))
      {
         fseek ( batchfle, solver.endfle, SEEK_SET );

         kk = null;
         while( kk < 5 )
         {
            fscanf( batchfle, "%s", ptr );
            if ( null != strstr( ptr, "terminated" ))
            {   
               fprintf( display, "\n Message on '%s':\n", bchfle );
               while( kk < 10 )
               {
                  fprintf( display, " %s", ptr );
                  fscanf( batchfle, "%s", ptr );
                  kk++;
               };
            };/* end if .. */
            kk++;
         };/* end while ... */
         fprintf( display, "\n All jobs marked %cok%c"
                            " [terminated].", 34, 34 );
      };
      fclose( batchfle );

   };/* end if null < ind */
/*............................................................................*/
   if ( *solver.opt == null )
   {
      PRBLDCLR( "\r" );
      fprintf( display, "\r %*s", 78, "DSC SOLVER" );
      PRNORMAL( "\r " );
   };

   fprintf( display, "\n ======================================"
                      "========================================" );
   if( ind <= null )
      goto stop_user_time;
   
   if( ind == THREE ) /* i.e. restart from jobs loop, chk2_nbrjbs     */
         goto start;  /* return to jobs loop without setting user tme */
#endif
/*............................................................................*/
/* set user time: */

# if SYS_TIMES == 2
   cl_ptr   = &cpu;
   ind      = times( cl_ptr );
   usr_time = ( cl_ptr -> tms_utime );
# elif SYS_TIMES == 1
   time0 = time( timer );
# else
   usr_time = clock( );
# endif

   if ( solver.nbrjbs <= null )
      goto stop_user_time;
/*............................................................................*/
/* here start the jobs */

  start:

   jj = solver.rstart;
   while( jj < solver.nbrjbs )
   {
/*............................................................................*/
# if DSC_INIT == 1
      ind = clearv( NODES, NRSMX, GENDS, GESMX, GMNDS, GMSMX, EXCEP,
                    EXCHP, BNDPR, BNDAP, EVLEP, EVLEN, EVLHP, EVLHN );
# endif
/*............................................................................*/
      nseconds = time( timer );
        lt_ptr = localtime( &nseconds );
           hrs = ( short )( lt_ptr -> tm_hour );
           min = ( short )( lt_ptr -> tm_min  );
           sec = ( short )( lt_ptr -> tm_sec  );
           ind = ( short )( lt_ptr -> tm_wday );

      fprintf( display, "\n\n -------------------- DSC job no %-3d started: %s"
         " %02d:%02d:%02d --------------------",
         joblbl[jj], wday[ind], hrs, min, sec );
/*............................................................................*/
/* set job time: */  

# if SYS_TIMES == 2
      ind = times( cl_ptr );
      job_time = ( cl_ptr -> tms_utime );
# elif SYS_TIMES == 1
      time2 = time( timer );
# else
      job_time = clock( );
# endif
/*............................................................................*/
/* reset error message to '***' [ void ]: */

      strncpy( solver.errmsg, astrx, LGS_SIZE );
/*............................................................................*/
/* enter DSC system files: */

/* enter DSC mesh topology from file elf.top<*>: */
      ind = toplgy( jj );

      if( ind == null )	/* ind = null: abnormal return [ error ], */
         goto valuefle;	/* write error message into evaluation file */

/* enter s-parameters from file elf.smx<*>: */
      ind = smtrix( jj );

      if( ind == null ) /* abnormal return [ error ] */           
         goto valuefle;	/* write error message into evaluation file */

/* enter boundary conditions from file elf.bnd<*>: */
      ind = boundr( jj );	

      if( ind == null ) /* abnormal return [ error ] */
         goto valuefle;	/* write error message into evaluation file */

/* enter excitation mode from file elf.exc<*>: */
      ind = excite( jj );

      if( ind == null ) /* abnormal return [ error ] */
         goto valuefle;	/* write error message into evaluation file */

/* enter evaluation mode */
/* [ evaluated cells and ports etc.] from file elf.val<*> and */
/* store operation parameters or error messages in that file: */

     valuefle:
      ind = values( jj );
/*............................................................................*/
/* DSC process:  */

      if ( ind == ONE ) /* ind == ONE: normal return from function values(*)  */
      {
         if( *solver.opt != null )
         {
            solver.logps = ftell( display );
            fclose( display );
         };
/*............................................................................*/
         ind = solvdrv( jj );         /*                                      */
/*..................................*/

         if( *solver.opt != null )
         {
            display = fopen( solver.logfle, "r+" );
            fseek( display, solver.logps, SEEK_SET );
         };

         if( ind == ONE ) /* ind == ONE: normal return from solvdrv(*)        */
         {
/* stop job time: */

# if DSC_STOPTIME == 0
            nseconds = time( timer );
            lt_ptr = localtime( &nseconds );
               hrs = ( short )( lt_ptr -> tm_hour );
               min = ( short )( lt_ptr -> tm_min  );
               sec = ( short )( lt_ptr -> tm_sec  );
               ind = ( short )( lt_ptr -> tm_wday );

            if( *solver.opt == null )
            {
               fprintf( display, "\r Job no%3d terminated: %s %02d:%02d:%02d\n",
                  joblbl[jj], wday[ind], hrs, min, sec );
            }
            else
            {
               fprintf( display, "\n Job no%3d terminated: %s %02d:%02d:%02d",
                  joblbl[jj], wday[ind], hrs, min, sec );
            };
# endif
/*............................................................................*/
/* stop and print job time: */

# if SYS_TIMES == 2
            ind = times( cl_ptr );
            job_time = ( cl_ptr -> tms_utime ) - job_time;
            job_time /= ticks;
# elif SYS_TIMES == 1
            time3 = time( timer ); 
            job_time = difftime( time3, time2 );
# else
            job_time = clock( ) - job_time;
            job_time /= CLOCKS_PER_SEC;
# endif
            scs = fabs(( double ) job_time );
            hrs = ( short ) floor ( scs/3600. );
            min = ( short ) floor (( scs - hrs*3600.)/60.);
            scs = ( double )( scs - hrs*3600. - min*60. );
            sec = ( short ) floor( scs );

            if ( *solver.opt == null )
            {
               fprintf( display, " CPU time: %03d hours, %02d minutes "
	          "and %.3e seconds", hrs, min, scs );
            }
            else
            {
               fprintf( display, "\n CPU time: %03d hours, %02d minutes "
	          "and %.3e seconds", hrs, min, scs );
            };
/*............................................................................*/
/* normal end of DSC process: */

            goto clear; /* job normally terminated */
         }
         else /* ind == null: unnormal return from solvdrv [ error ]  */
            ind = values( jj );       /* write error message                  */
      };                              /* into evaluation file                 */
/*............................................................................*/
/* on error [ from values(*) or solvdrv(*) ]: */

      fprintf( display, "\n ---------------------------- overriding "
         "job no.%3d ----------------------------\n", joblbl[jj] );

/*............................................................................*/
/* clear last job parameters: */

     clear: 
/*............................................................................*/
# if DSC_INIT == 2
      ind = clearv( top.n, smx.hh[null], smx.gm, smx.gm, smx.ge, smx.ge,
         exc.ne, exc.nh, bnd.p, bnd.n, val.nep, val.nen, val.nhp, val.nhn );
# endif
/*............................................................................*/
/* check state: bill "ok" [ or "??" if errors occured during DSC process ]:  */

/* check_state: */ 

      ind = THREE;/* op.mark for return to job_instruct [ perturbed batch f.] */
      ii = jj;
      
# if DSC_BATCH == 'b'

      batchfle = fopen( bchfle, "r+" );

      if ( batchfle == null )
      {
         fprintf( stderr,
            "\n\n Error message from program %s :", PROGRAM );
         fprintf( stderr,
            "\n Can't open job instructions file %s !", bchfle );
         fprintf( stderr,
            "\n [ Check file location and/or permissions.]" );
         fprintf( stderr,
            "\n Program %s stopped.\n", PROGRAM );

         exit( EXIT_FAILURE );
      }
      else
      {
         fscanf( batchfle, "%s", ptr ); /* comment */

         fscanf( batchfle, "%s", ptr ); /* new format: 'remove/retain_ files'*/
	 if ( null == strncmp( ptr, "number_of_jobs", SIX )) /* [ old format ]*/
	 {
            fscanf( batchfle, "%s", ptr ); /* number of jobs */
            kk = strtol( ptr, endp, DEC ); /* [ which can be modified */
                                           /* during program run time ] */
            if ( solver.nbrjbs != kk )
               solver.nbrjbs = kk;

            fscanf( batchfle, "%s", ptr ); /* string: 'first_job_index' */
            fscanf( batchfle, "%s", ptr ); /* index */
            fscanf( batchfle, "%s", ptr ); /* string: 'remove/retain_files' */
         };

         if ( null == strncmp( rmfle, ptr, SIX )) /* case: remove__files_ */
         { 
# if DSC_PRTCTFL == 1
            solver.rmfle = TWO;
# else
            solver.rmfle = ONE;
# endif
         }
         else
         {
            solver.rmfle = null;
         };

         fscanf( batchfle, "%s", ptr );
	 if ( null == strncmp( ptr, "number_of_jobs", SIX )) /* [new format] */
	 {
            fscanf( batchfle, "%s", ptr ); /* number of jobs */
            kk = strtol( ptr, endp, DEC ); /* [ which can be modified */
                                           /* during program run time ] */
            if ( solver.nbrjbs != kk )
               solver.nbrjbs = kk;
         };

         fseek( batchfle, solver.stofs, SEEK_SET );
         fscanf( batchfle, "%s", status ); /* state should be '--' */

         if ( null != strncmp( status , "--" , TWO )) /* perturbed batch file:*/
            goto job_instruct;                        /* restart program      */

/*............................................................................*/
/* bill state "ok" or "??" [ error ]: */

         fseek( batchfle, solver.stofs, SEEK_SET );

         if ( null == strncmp( solver.errmsg, astrx, THREE ))
         {
            fprintf( batchfle, " %6s", "ok" );   /* bill 'ok' state */
            fprintf( batchfle, " %03dh:%02dm:%02ds", hrs, min, sec );
         }                                       /* bill CPU_time */
         else
         {
            fprintf( batchfle, " %6s", "??" );
            fprintf( batchfle, " %.9s%.3s", solver.errmsg, "..." );

# if DSC_PRTCTFL == 1
            if ( solver.rmfle == TWO ) /* protect files from removal */
            {
               strcpy( ptr, lotos( toplbl[ii], null ));
               if (( null == strstr( topfls, ptr ))&&
                   ( strlen( topfls ) < ( FOUR*DSC_JOBS )))
               {
                  strcat( topfls, ptr );
                  strncat( topfls, astrx, ONE );
               };

               strcpy( ptr, lotos( smxlbl[ii], null ));
               if (( null == strstr( smxfls, ptr ))&&
                   ( strlen( smxfls ) < ( FOUR*DSC_JOBS )))
               {
                  strcat( smxfls, ptr );
                  strncat( smxfls, astrx, ONE );
               };

               strcpy( ptr, lotos( bndlbl[ii], null ));
               if (( null == strstr( bndfls, ptr ))&&
                   ( strlen( bndfls ) < ( FOUR*DSC_JOBS )))
               {
                  strcat( bndfls, ptr );
                  strncat( bndfls, astrx, ONE );
               };

               strcpy( ptr, lotos( exclbl[ii], null ));
               if (( null == strstr( excfls, ptr ))&&
                   ( strlen( excfls ) < ( FOUR*DSC_JOBS )))
               {
                  strcat( excfls, ptr );
                  strncat( excfls, astrx, ONE );
               };
            };
# endif
         };

         fseek( batchfle, solver.stofs, SEEK_SET );

         fscanf( batchfle, "%s", ptr ); /* state */
         fscanf( batchfle, "%s", ptr ); /* CPU_time/err */

         hh = null;
         ll = null;
         pp = null;
         qq = null;
	 
         if ( DSC_JOBS < solver.nbrjbs )
         {
            fprintf( display, "\n Too many jobs defined in batch file !!! " );
            fprintf( display, "\n [ reduced to maximum number %3d ]",
                                                                  DSC_JOBS );
            solver.nbrjbs = DSC_JOBS; /* maximum number of jobs */
         };

         solver.stofs = null;
         kk = jj+ONE;
	 jj = solver.nbrjbs;
         while ( kk < solver.nbrjbs )
         {
            fscanf( batchfle, "%s", ptr ); /* joblbl[kk] */

            if ( joblbl[kk] != strtol( ptr, endp, DEC ))
               goto job_instruct;          /* [ re-enter modified file ] */

            fscanf( batchfle, "%s", ptr ); /* topology file index */

            if ( toplbl[kk] != strtol( ptr, endp, DEC ))
               goto job_instruct;               

            fscanf( batchfle, "%s", ptr ); /* s-matrix file index */

            if ( smxlbl[kk] != strtol( ptr, endp, DEC ))
               goto job_instruct;

            fscanf( batchfle, "%s", ptr ); /* boundary conf. file */

            if ( bndlbl[kk] != strtol( ptr, endp, DEC ))
               goto job_instruct;

            fscanf( batchfle, "%s", ptr ); /* excitation file */
                     
            if ( exclbl[kk] != strtol( ptr, endp, DEC ))
               goto job_instruct;

            fscanf( batchfle, "%s", ptr ); /* evaluation file */

            if ( joblbl[kk] != strtol( ptr, endp, DEC ))
               goto job_instruct;

            mm = ftell( batchfle );

            fscanf( batchfle, "%s", status );

            if ( null == strncmp( status, "--" , TWO ))
            {
               if ( solver.stofs == null )
               {
                  solver.stofs = mm;
                  jj = kk;
               };

               if ( solver.rmfle != null ) /* check if files 'DSC_PRFX...[ii]' */
               {                         /* will be used in further jobs      */
                  if ( toplbl[kk] == toplbl[ii] )
                     hh = ONE;

                  if ( smxlbl[kk] == smxlbl[ii] )
                     pp = ONE;
                
                  if ( bndlbl[kk] == bndlbl[ii] )
                     qq = ONE;

                  if ( exclbl[kk] == exclbl[ii] )
                     ll = ONE;
               };
            };
            kk++;

            fscanf( batchfle, "%s", ptr ); /* '------------' [ CPU_time/err ] */
         };

         solver.endfle = ftell( batchfle );
          
         fclose( batchfle );

      };/* end if batchfle != null */
/*............................................................................*/
/* remove files if solver.rmfle == TWO: */

      if ( solver.rmfle == TWO ) /* remove used files 'DSC_PRFX...'   */
      {                        /* [ ii = job index just terminated ] */
         if ( hh == null )
         {
/*............................................................................*/
# if DSC_PRTCTFL == 1
            strcpy( ptr, lotos( toplbl[ii], null ));
            if ( null == strstr( topfls, ptr ))
            {
               ind = remove ( solver.top );

               if ( ind == null ) 
                  fprintf( display, "\n removed: topology file %s",
                     solver.top );
            }
            else
               fprintf( display, "\n retained: topology file %s",
                  solver.top );
# else
            ind = remove ( solver.top );

            if ( ind == null ) 
               fprintf( display, "\n removed: topology file %s",
                  solver.top );
# endif
/*............................................................................*/
         };

         if ( pp == null )
         {
/*............................................................................*/
# if DSC_PRTCTFL == 1
            strcpy( ptr, lotos( smxlbl[ii], null ));
            if ( null == strstr( smxfls, ptr ))
            {
               ind = remove ( solver.smx );

               if ( ind == null ) 
                  fprintf( display, "\n removed: s-matrix file %s",
                     solver.smx );
            }
            else
               fprintf( display, "\n retained: s-matrix file %s",
                  solver.smx );
# else
            ind = remove ( solver.smx );

            if ( ind == null ) 
               fprintf( display, "\n removed: s-matrix file %s",
                  solver.smx );
# endif
/*............................................................................*/
         };

         if ( qq == null )
         { 
/*............................................................................*/
# if DSC_PRTCTFL == 1
            strcpy( ptr, lotos( bndlbl[ii], null ));
            if ( null == strstr( bndfls, ptr ))
            {
               ind = remove ( solver.bnd );

               if ( ind == null ) 
                  fprintf( display, "\n removed: boundary file %s",
                     solver.bnd );
            } 
            else
               fprintf( display, "\n retained: boundary file %s",
                  solver.bnd );
# else
            ind = remove ( solver.bnd );

            if ( ind == null ) 
               fprintf( display, "\n removed: boundary file %s",
                  solver.bnd );
# endif
/*............................................................................*/
         };

         if ( ll == null )
         {
/*............................................................................*/
# if DSC_PRTCTFL == 1
            strcpy( ptr, lotos( exclbl[ii], null ));
            if ( null == strstr( excfls, ptr ))
            {
               ind = remove ( solver.exc );

               if ( ind == null ) 
                  fprintf( display, "\n removed: excitation file %s",
                     solver.exc );
            }
            else
               fprintf( display, "\n retained: excitation file %s",
                  solver.exc );
# else
            ind = remove ( solver.exc );

            if ( ind == null ) 
               fprintf( display, "\n removed: excitation file %s",
                  solver.exc );
# endif
/*............................................................................*/
         };
      };

# endif
/*............................................................................*/
/* next_job: */ 

      fprintf( display, "\n ======================================="
                        "=======================================" );

   }; /* while ( jj < solver.nbrjbs ) */

/*....................... program terminated:<date> ..........................*/
# if DSC_BATCH == 'b'

   if ( solver.rstart < solver.nbrjbs )
   { 
      nseconds = time( timer );
      tmeptr = ctime( &nseconds );
      batchfle = fopen( bchfle ,"r+" );
      fseek( batchfle , solver.endfle , SEEK_SET );
      fprintf( batchfle, "\n\n%s%s%s\n%.24s", "DSC batch jobs '",
         bchfle,"' terminated:", tmeptr );
      fclose ( batchfle );
   };

# endif
/*................................ user time ................................*/
/* stop and print user time: */

  stop_user_time:
         
   nseconds = time( timer );
   tmeptr = ctime( &nseconds );
   fprintf( display, "\n DSC program %s terminated: %.24s",
                                                            PROGRAM, tmeptr );

# if SYS_TIMES == 2
   ind = times( cl_ptr );
   usr_time = ( cl_ptr -> tms_utime ) - usr_time;
   usr_time /= ticks;
# elif SYS_TIMES == 1
   time1 = time( timer );
   usr_time = difftime( time1, time0 );
# else
   usr_time = clock( ) - usr_time;
   usr_time /= CLOCKS_PER_SEC;
# endif

   scs = fabs(( double ) usr_time );
   hrs = ( short ) floor ( scs/3600. );
   min = ( short ) floor (( scs - hrs*3600.)/60.);
   scs = ( double )( scs - hrs*3600. - min*60. );

   fprintf( display, "\n The CPU time was: %02d hours, %02d minutes and"
           " %1.3e seconds.", hrs, min, scs );

   if( *solver.opt == null ) /* terminal mode */
   {
      PRINVERS( "" );
      fprintf( display, "\t  Bye - bye !!!" );
      PRNORMAL( "\n " );
   }
   else
      fprintf( display, "\t  Bye - bye !!!" );
/*............................................................................*/
/* normal program end: */

# if DSC_DAEMON == 1
   remove( dpt->pidfle );  /* process identity number */
# endif

   exit( EXIT_SUCCESS );
}  
/*============================================================================*/
/**************** end of DSC analysis main program solver.c *******************/
