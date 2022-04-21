/* [ file: former.c ] */
# define PROGRAM "FORMER"
/*******************************************************************************
*                                                                              *
*   ANSI C program former.c                                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Danse preprocessing [ model and process definition ] routine,              *
*   creates input files for the DSC analysis program solver.c                  *
*                                                                              *
*   This  program  creates  structure  topology  file  'DSC_PRFX.top',         *
*   S-matrix  file  'DSC_PRFX.smx' ,  boundary  file  'DSC_PRFX.bnd',          *
*   as well as  field excitation and evaluation files  'DSC_PRFX.exc',         *
*   'DSC_PRFX.val' - used as batch files in DSC main program solver.c.         *
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
# include "FORMER.CONF"
/*----------------------------------------------------------------------------*/
/* Don't mix different releases: */
# ifndef FORMER_RELEASE
   # define FORMER_RELEASE DSC_RELEASE
# endif
# ifndef FORMER_UPDATE
   # define FORMER_UPDATE DSC_UPDATE
# endif
/*----------------------------------------------------------------------------*/
# ifndef USE_NCURSES
   # define USE_NCURSES 1
# endif
/*----------------------------------------------------------------------------*/
/* This may yet be defined in the general configuration header CONFIG.H: */
/*
*//* # define SYS_TIMES 0 *//* 2: use function times(*)                       */
                            /* 1: use function difftime(*) [ ANSI C standard ]*/
                            /* 0: use function clock(*)    [ ANSI C standard ]*/
# if SYS_TIMES == 2
   # include <sys/times.h>  /* system time: cf. times(*)                      */
# endif
/*----------------------------------------------------------------------------*/
# include "../tools/txctyp.h" /* the text console structure type */
# include "../former/formtp.h" /* some structure type definitions */
/*----------------------------------------------------------------------------*/
/* some further operational constants */

# define DISP 1              /* set DISP 1 for displaying intermediate results*/
# define EXTENDED 0          /* EXTENDED (0) 1 for (reduced) extended files   */
# define CLEAR 0             /* CLEAR    (0) 1 for parameter clearing         */
# define BUFFER 1024
/*----------------------------------------------------------------------------*/
                                       /*=== system orientated declarations =>*/
                                       /*--- file pointers ------------------>*/
                                       /*--- system functions --------------->*/
time_t time( time_t *timer );                
char *ctime( const time_t *timer );          

struct tm *localtime( const time_t *timer ); 
char *asctime( const struct tm *lt_ptr );

# if SYS_TIMES == 2
   long sysconf(int s_name); 
   clock_t times( struct tms *cl_ptr );
# elif SYS_TIMES == 1
   double difftime( time_t time1, time_t time0 );
# else
   clock_t clock( void );
# endif

int setvbuf( FILE *display, char *buff, int unbuff, size_t bufsize );
 
                                       /*--- global structures -------------->*/
static struct tm loct    = {null};
static struct tm *lt_ptr = &loct;

# if SYS_TIMES == 2
static struct tms cpu    = {null};
static struct tms *cl_ptr= &cpu;
# endif
                                       /*--- global variables --------------->*/
static int unbuff     = _IONBF;
static size_t bufsize =  null;

static char *err = "job_error";
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* 'my_terminal' configuration: */

   # include <ncurses.h>
   # include <curses.h>
   # include <term.h> /* the terminal type header */

   static char *term;        /* terminal type string */ 

   # define CLSCREEN { \
     fprintf( stdout, "%s", tgetstr( "cl", null )); \
   }

   # define PRBLDCLR(a) { \
     fprintf( stdout, \
        "%s%s", tgetstr( "md", null ), (a)); /* bold clear output */ \
   }

   # define PRINVERS(a) { \
     fprintf( stdout, "%s%s", tgetstr( "mr", null ), (a)); /* inverse */ \
   }

   # define PRNORMAL(a) { \
     fprintf( stdout, \
        "%s%s", tgetstr( "me", null ), (a)); /* back to normal output */ \
   }
# else
   # define CLSCREEN { \
     fprintf( stdout, "\f" );\
   }

   # define PRBLDCLR(a) { \
     fprintf( stdout, "%s", (a));\
   }

   # define PRINVERS(a) { \
     fprintf( stdout, "%s", (a));\
   }

   # define PRNORMAL(a) { \
     fprintf( stdout, "%s", (a));\
   }
# endif
/*============================================================================*/

int main( int argc, char **argv )
{
/*--- system function prototypes -------------------------------------------->*/

   long   strtol( const char *ptr, char **endp, int base );
   double strtod( const char *ptr, char **endp );

   time_t time( time_t *timer );
   char *ctime( const time_t *timer );

   struct tm *localtime( const time_t *timer );
   char *asctime ( const struct tm *lt_ptr );

# if SYS_TIMES == 2
   long sysconf(int s_name); 
   clock_t times( struct tms *cl_ptr );
# elif SYS_TIMES == 1
   double difftime( time_t t1, time_t t0 );
# else
   clock_t clock( void );
# endif

   int setvbuf( FILE *display, char *buff, int mode, size_t bfsize );

# if USE_NCURSES == 1
   int setupterm( NCURSES_CONST char *term, int fildes, int *errret );
# endif

/*--- system  orientated declarations and allusions ------------------------->*/

   static short 
      hrs = null,
      min = null;

   static time_t
      nseconds = null,
     *timer    = null;

   static double 
      seconds = ZERO;

# if SYS_TIMES == 2
   static clock_t
      ticks = null,
      usr_time = null;
# elif SYS_TIMES == 1
   static time_t
      time0 = null,
      time1 = null;

   static double
      usr_time = ZERO;
# else
   static clock_t
      usr_time = null;
# endif

/*--- user orientated declarations and allusions -------------------------->*/

   static TXCNSL
     *csp;

   static char 
     ptr[STS_SIZE] = {null},
    *tmeptr;

   static short 
      job_number = null,
      first_job  = null,
      idx        = null;         
  
   static unsigned char 
      ind = null;

/*--- user function prototypes: -------------------------------------------->*/

   short 
      formdrv( char *err );   

   TXCNSL 
      *txcnsl( TXCNSL *csp );

/*---------------------- end of declaration part -----------------------------*/
/* read the command line: */

   if( null < --argc )
   {
      do
      {
         strcpy( ptr, *++argv );

         if (( null == strcmp( ptr, "v" ))||
             ( null == strcmp( ptr, "-v" )))
         {
            fprintf( stderr, "\rDANSE release %s;", FORMER_RELEASE );
            fprintf( stderr, " %s\n", FORMER_UPDATE );

            exit( EXIT_SUCCESS );
         }
	 else
	 {
            fprintf( stderr, "\runknown program option 'former %s'\n", ptr );
            exit( EXIT_FAILURE );
         };
      } while ( null < ( --argc ));
   }; /* end if ( null < argc ) */
/*............................................................................*/
# if USE_NCURSES == 1
/* set the terminal type: */

# if defined ( _BSD )
   setupterm( "cons25", 1, ( int *) 0 );
# elif defined ( _GNU_Linux )
   setupterm( "linux", 1, ( int *) 0 );
# elif defined ( _Linux )
   setupterm( "linux", 1, ( int *) 0 );
# endif

/* get the terminal info: */

   term = ( char *) getenv( "TERM" ); /* get the terminal */

   ind = tgetent( null, term );

   if( ONE != ind )
   {
      fprintf( stderr, "Error on getting the termcap info\n " ); 
      exit( EXIT_FAILURE );
   };
# endif
/*............................................................................*/
/* set buffer length = null: */

   ind = setvbuf( stdin, null, unbuff, bufsize );/*set buffer length = null*/
   ind = setvbuf( stdout, null, unbuff, bufsize ); 

/*............................................................................*/
/* allocate memory: */

   tmeptr = ( char * ) calloc( STS_SIZE, ONE );
   
/*............................................................................*/
/* stop the startup time: */

   nseconds = time( timer );
   tmeptr = ctime( &nseconds );
/* [ same as 'tmeptr = asctime(localtime(&nseconds));' ] */
/*............................................................................*/
/* display | menu: */

   csp = txcnsl( null ); /* initialize the text console */
/*
   CLSCREEN;
   PRNORMAL( "\f\n\n\n" );
*/
   ( csp->clscr ) = 2; /* clear screen; scroll that number of lines */
   strcpy(( csp->title ), "DSC program " );
   strcat(( csp->title ), PROGRAM );
   strcat(( csp->title ), " started: " );
   strcat(( csp->title ), ptr );
   strcat(( csp->title ), " The DSC system function is " );
   strcat(( csp->title ), DSC_SYSTEM );
/*............................................................................*/
   csp = txcnsl( csp );   /*                                                  */
/*......................*/
/* the function shoud do something like that: */
/*
   fprintf( stdout, "DSC program %s started: ", PROGRAM );
   PRNORMAL( "" );
   fprintf( stdout, " %s ", tmeptr );
   fprintf( stdout, "The DSC system function is '%s'.", DSC_SYSTEM );
   fprintf( stdout, "\n ======================================="
      "=======================================" );
*/
/*............................................................................*/
   job_number = ONE;
   first_job  = null;
/*............................................................................*/
/* set_time: */

# if SYS_TIMES == 2
   cl_ptr   = &cpu;
   ind      = times( cl_ptr );
   usr_time = ( cl_ptr->tms_utime );
# elif SYS_TIMES == 1
   time0 = time( timer );
# else
   usr_time = clock( );
# endif

/*............................................................................*/
/* iterate: */

# if CLEAR == 1
   fprintf( stdout, "\n initializing parameters" );
   ind = clear( ); 
   fprintf( stdout, "\r                         \r" );
# endif 

   for ( idx = first_job; idx < (first_job + job_number); idx++ )
   {
      nseconds = time(timer);
      lt_ptr = localtime(&nseconds);
      ind = ( short ) ( lt_ptr->tm_wday );
      hrs = ( short ) ( lt_ptr->tm_hour );
      min = ( short ) ( lt_ptr->tm_min  );
      seconds = ( double ) ( lt_ptr->tm_sec  );
/*============================================================================*/
/* here starts the job [number idx]: */   
            
/*
      write job specific instructions between the following braces: 
*/
      {
         ind = formdrv( err );            /*  call job  'formdrv'             */
      };

/* and here the job [number idx] is terminated */

# if DISP == 1
      if ( ind == ONE ) 
      {
         fprintf( stdout, " \n Abnormal job termination.  " );
         goto job_err;   
      };
# endif
                    
/*----------------------------------------------------------------------------*/
      nseconds = time(timer); 
      lt_ptr = localtime(&nseconds);
      ind = ( short ) ( lt_ptr->tm_wday );
      hrs = ( short ) ( lt_ptr->tm_hour );
      min = ( short ) ( lt_ptr->tm_min  );
      seconds = ( double ) ( lt_ptr->tm_sec  ); 

      goto next_job; 
/*............................................................................*/
     job_err:

# if CLEAR == 1
      ind = clear ( );           
# endif

      fprintf( stdout, "\n -------------------------- error: job overridden"
         " ----------------------------\n" );

/*......................... next job: index i+1 ..............................*/

     next_job:;

# if CLEAR == 1 
      ind = clear( );
# endif

      fprintf( stdout, "\n ======================================"
              "========================================\n" );

      nseconds = time( timer );
      tmeptr = ctime( &nseconds );

      fprintf( stdout, " Program %s terminated:  %s ", PROGRAM, tmeptr );
   };
/*............................... user time ..................................*/
/* stop_time: */

# if SYS_TIMES == 2
   ind = times( cl_ptr );
   usr_time = ( cl_ptr->tms_utime ) - usr_time;
   ticks = sysconf( _SC_CLK_TCK ); /* ticks = clocks per second -> */
   usr_time /= ticks;              /* -> cf. header <unistd.h>     */
# elif SYS_TIMES == 1
   time1 = time( timer );
   usr_time = difftime( time1, time0 );
# else 
   usr_time = clock( ) - usr_time;
   usr_time /= CLOCKS_PER_SEC;
# endif

   hrs     = ( short ) floor( usr_time/3600.);
   min     = ( short ) floor(( usr_time - hrs*3600.)/60.);
   seconds = ( double ) usr_time - hrs*3600. -min*60.;

   fprintf( stdout, "The CPU time was: %.2d hours, %.2d minutes and "
      "%1.3e seconds. ", hrs, min, seconds ); 

   fprintf( stdout, "\t  " );
   PRINVERS( "" );
   fprintf( stdout, "Bye - bye !!!" );
   PRNORMAL( "\n " );
/*............................................................................*/
/* exit main program FORMER: */

   exit( EXIT_SUCCESS );
}
/*============================================================================*/
/* here, the DSC model defining functions
   sygrid(*), sycord(*), sybndr(*), sysexc(*), sysval(*)
   may be included.

   Alternatively, these  functions  can be included  at the end of
   the  formdrv.h  function  [ which, however, requires  some care
   on defining static variables for these model functions ].                  */
/*
# include DSC_SYSTEM
*/
/*============================================================================*/
/*************** end of DSC system definition program former.c ****************/
