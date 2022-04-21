/* [ file: poster.c ] */
# define PROGRAM "POSTER"
/*******************************************************************************
*                                                                              *
*   ANSI C program poster.c                                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   This program evaluates and processes solver result files values<N>,        *
*   created by DSC analysis programs former.c and solver.c.                    *
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
# include <float.h>
# include <math.h>
# include <limits.h>
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
/* Edit and customize this header for POSTER configuration: */
# include "POSTER.CONF"
/*----------------------------------------------------------------------------*/
/* Don't mix different releases: */
# ifndef POSTER_RELEASE
   # define POSTER_RELEASE DSC_RELEASE
# endif
# ifndef POSTER_UPDATE
   # define POSTER_UPDATE DSC_UPDATE
# endif
/*----------------------------------------------------------------------------*/
# ifndef USE_NCURSES
   # define USE_NCURSES 1 /* =1 on BSD or Linux [GNU] systems, e.g.*/
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
# include "../poster/posttp.h" /* some structure type definitions */
/*----------------------------------------------------------------------------*/
/* operational constants, flags: */

# define DISP 1              /* set DISP 1 for displaying intermediate results*/
# define EXTENDED 0          /* EXTENDED (0) 1 for (reduced) extended files   */
# define BATCH 0             /* BATCH    (0) 1 for batch file operation       */
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
   long sysconf( int s_name ); 
   clock_t times( struct tms *cl_ptr );
# elif SYS_TIMES == 1
   double difftime( time_t time1, time_t time0 );
# else
   clock_t clock( void );
# endif

int setvbuf( FILE *display, char *buff, int unbuff, size_t bufsize );
 
/*--- global structures ----------------------------------------------------->*/
static struct tm loct    = {null};
static struct tm *lt_ptr = &loct;

# if SYS_TIMES == 2
static struct tms cpu    = {null};
static struct tms *cl_ptr= &cpu;
# endif

/*--- global variables ------------------------------------------------------>*/
static int unbuff     = _IONBF;
static size_t bufsize =  null;

# if BATCH == 1
static char *wday[] = { "Sun","Mon","Tue","Wed","Thu","Fri","Sat" };
# endif

static char *err = "job_error";
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* 'my_terminal' configuration: */

   # include <ncurses.h>
   # include <curses.h>
   # include <term.h> /* terminal type header */

   static char *term; /* terminal type string */ 

   # define CLSCREEN { \
     fprintf( stdout, "%s", tgetstr( "cl", null )); \
   }

   # define PRBLDCLR(a) { \
     fprintf( stdout, \
        "%s%s", tgetstr( "md", null ), (a)); /* bold clear output */ \
   }

   # define PRINVERS(a) { \
     fprintf( stdout, \
        "%s%s", tgetstr( "mr", null ), (a)); /* inverse */ \
   }

   # define PRNORMAL(a) { \
     fprintf( stdout, \
        "%s%s", tgetstr( "me", null ), (a)); /* back to normal output */ \
   }
# else
   # define CLSCREEN { \
     fprintf( stdout, "\f" ); \
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

   long strtol( const char *ptr, char **endp, int base );
   double strtod( const char *ptr, char **endp );

   time_t time( time_t *timer );
   char *ctime( const time_t *timer );

   struct tm *localtime( const time_t *timer );
   char *asctime ( const struct tm *lt_ptr );

# if SYS_TIMES == 2
   long sysconf( int s_name );
   clock_t times( struct tms *cl_ptr );
# elif SYS_TIMES == 1
   double difftime( time_t time1, time_t time0 );
# else
   clock_t clock( void );
# endif

   int setvbuf( FILE *display, char *buff, int mode, size_t bfsize );

# if USE_NCURSES == 1
   int setupterm( NCURSES_CONST char *term, int fildes, int *errret );
# endif

/*--- system orientated declarations and allusions -------------------------->*/

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

/*--- user orientated declarations and allusions ---------------------------->*/

   static TXCNSL
      *csp;

   static char 
      ptr[STS_SIZE] = {null},
     *tmeptr;

# if BATCH == 1
    **endp = null; 
# endif

   static short 
      job_number  = null,
      first_job   = null,
      ind         = null,
      idx         = null;         

/*--- user function prototypes: --------------------------------------------->*/

   short 
      postdrv( char *err );   

   TXCNSL 
      *txcnsl( TXCNSL *csp );

/* int clear( int i, int j, int k, int h,int n, int m ); */

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
            fprintf( stdout, "\rposter release %s;", POSTER_RELEASE );
            fprintf( stdout, " %s\n", POSTER_UPDATE );
            exit( EXIT_SUCCESS );
         }
	 else
	 {
            fprintf( stderr, "\runknown program option 'poster %s'\n", ptr );
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
      fprintf( stderr, "Error on getting the termcap info\n" ); 
      exit( EXIT_FAILURE );
   };
# endif
/*............................................................................*/
/* memory allocations: */

   tmeptr = ( char *) calloc( STS_SIZE , ONE );
 
/*............................................................................*/
/*set buffer length = null:*/

   setbuf( stdin, null );
   ind = setvbuf( stdin, null, unbuff, bufsize );
   ind = setvbuf( stdout, null, unbuff, bufsize ); 

/*............................................................................*/
/* set startup time: */

   nseconds = time( timer );
   tmeptr = ctime( &nseconds );
/* [ same as 'ptr = asctime(localtime(&nseconds));' ] */
/*............................................................................*/
/* display | menu: */

   csp = txcnsl( null ); /* initialize text console */
/*
   CLSCREEN;
   PRNORMAL( "\f\n\n\n" );
*/
   ( csp->clscr ) = 2; /* clear screen; scroll that number of lines */
   strcpy( csp->title, "DSC program " );
   strcat( csp->title, PROGRAM );
   strcat( csp->title, " started: " );
   strcat( csp->title, tmeptr );

# ifdef DSC_SYSTEM
   strcat( csp->title, " The DSC system function is " );
   strcat( csp->title, DSC_SYSTEM );
# endif
/*............................................................................*/
   csp = txcnsl( csp );   /*                                                  */
/*......................*/

/*............................................................................*/
/* This function shoud do somthing like that: */
/*
   fprintf( stdout, "DSC program %s started: ", PROGRAM );
   PRNORMAL( "" );
   fprintf( stdout, " %s ", tmeptr );

# ifdef DSC_SYSTEM
   fprintf( stdout, "The DSC system function is '%s'.", DSC_SYSTEM );
   fprintf( stdout, "\n ======================================="
# endif
   fprintf( stdout, "\n ==============================="
       "===============================================" );
*/
/*............................................................................*/
# if BATCH == 1
   fprintf( stdout, "\n\n please enter number of jobs...................: " );
   scanf( "%s", ptr);
   job_number = strtol( ptr, endp, DEC );
   fprintf( stdout, "\n please enter index of first job ( n >= 0 )....: " );
   scanf( "%s", ptr );
   first_job = strtol( ptr, endp, DEC );
   fprintf( stdout, "\n ==============================="
       "===============================================" );
# else
   job_number = ONE;
   first_job  = null;
# endif
/*............................................................................*/
/* set_time: */ 

# if SYS_TIMES == 2
   cl_ptr = &cpu;
   ind = times( cl_ptr );
   usr_time = ( cl_ptr -> tms_utime );
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

   for( idx = first_job; idx< (first_job + job_number); idx++ )
   {
      nseconds = time(timer);
      lt_ptr = localtime(&nseconds);
      ind = ( short ) ( lt_ptr -> tm_wday );
      hrs = ( short ) ( lt_ptr -> tm_hour );
      min = ( short ) ( lt_ptr -> tm_min  );
      seconds = ( double ) ( lt_ptr -> tm_sec  );
# if BATCH == 1
      fprintf( stdout,
         "\n\n --------------- job no. %.2d  started:    %s %.2d:%.2d:%.2f"
         " --------------------------", 
	 idx, wday[ind], hrs, min, seconds );
# endif
/*----------------------------------------------------------------------------*/

/* job: */
            
/* job specific instructions: */

      {
         ind = postdrv( err ); /* call model configuration function postdrv(*)*/
      };

/* end_of_job: */

# if DISP == 1
      if ( ind == ONE )
      {
         fprintf( stdout, " \n abnormal job termination.  " );
         goto job_err;   
      };
# endif
                    
/*----------------------------------------------------------------------------*/
      nseconds = time(timer); 
      lt_ptr = localtime(&nseconds);
      ind = ( short ) ( lt_ptr -> tm_wday );
      hrs = ( short ) ( lt_ptr -> tm_hour );
      min = ( short ) ( lt_ptr -> tm_min  );
      seconds = ( double ) ( lt_ptr -> tm_sec  );
# if BATCH == 1
      fprintf( stdout,
         "\n --------------- job no. %.2d  terminated: %s %.2d:%.2e:"
         " %.2d -------------------------\n", 
	 idx, wday[ind], hrs, min, seconds);
# endif

      goto next_job; 
/*............................................................................*/
     job_err:

# if BATCH == 1
      fprintf( stdout,
         "\n -------------------------- job no. %2d overridden "
         "-----------------------------\n", idx );  
# else
      fprintf( stdout,
         "\n --------------------------- error: job overriden "
         "-----------------------------\n" );
# endif

/*......................... next job: index i+1 ..............................*/

     next_job: ;

# if CLEAR == 1 
      ind = clear( );
# endif 

      fprintf( stdout, "\n ==============================="
          "===============================================\n" );

      nseconds = time(timer);
      tmeptr = ctime(&nseconds);

      fprintf( stdout, " Program %s terminated:  %s ", PROGRAM, tmeptr);
   };
/*............................... user time ..................................*/
/* stop_time: */

# if SYS_TIMES == 2
   ind = times( cl_ptr );
   usr_time = ( cl_ptr -> tms_utime ) - usr_time;
   ticks = sysconf( _SC_CLK_TCK ); /* ticks = clocks per second -> ... */
   usr_time /= ticks;
# elif SYS_TIMES == 1
   time1 = time( timer );
   usr_time = difftime( time1, time0 );
# else                             /* ... -> cf. header <sysconf.h> */
   usr_time = clock( ) - usr_time;
   usr_time /= CLOCKS_PER_SEC;
# endif

   hrs     = ( short ) floor( usr_time/3600.);
   min     = ( short ) floor(( usr_time - hrs*3600.)/60.);
   seconds = ( double ) usr_time - hrs*3600. - min*60.;

   fprintf( stdout, "The CPU time was: %.2d hours, %.2d minutes and "
      "%.2f seconds. ", hrs, min, seconds ); 

   fprintf( stdout, "\t  " );
   PRINVERS( "" );
   fprintf( stdout, "Bye - bye !!!" );
   PRNORMAL( "\n " );
/*............................................................................*/
/* end: */

   exit( EXIT_SUCCESS );
}
/*============================================================================*/
/************************** end of main program poster.c **********************/
