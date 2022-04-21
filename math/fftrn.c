/* [ file: fftrn.c ] */
# define FFT_NAME "fftrn.c"
/*******************************************************************************
*                                                                              *
*   ANSI C Program fftrn(*)                                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Input/output routine for fast Fourier transformation program fftrf.c       *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define _POSIX_SOURCE 1 /* some headers of the POSIX.1 standard will be used */
/*----------------------------------------------------------------------------*/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
# include <time.h>      /* cf. time( ), ctime( ),asctime( ), localtime( )     */
# include <sys/times.h> /* system time: cf. times( ), etc.                    */
# include <unistd.h>    /* system specification header, cf. sysconf( ), etc.  */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# include "./maths.h"  /* 'my' computation environment headers */
# include "./consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H"
/*----------------------------------------------------------------------------*/
# ifndef
   # define LINES 63
# endif
# ifndef
   # define WIDTH 80
# endif
/*----------------------------------------------------------------------------*/
/* double format [ "E" or "e" ]:                                              */

# define DBL_FORMAT "E"
/*----------------------------------------------------------------------------*/
/* operation marks:                                                           */

# define LIST_FILE 0
/*----------------------------------------------------------------------------*/
/* file_names: */

# define outfle "results"
# define spcfle "spectrum"
# define disfle "distribution"
/*----------------------------------------------------------------------------*/
/* system function prototypes: */
double cos ( double );
double sin ( double );
clock_t times( struct tms *cl_ptr );
long szsconf( int s_name );
/*----------------------------------------------------------------------------*/
/* user defined headers: */
# include "../math/FFT.CONF"
# include "../math/fftyp.h"
/*----------------------------------------------------------------------------*/
static struct tms cpu = {null};
static struct tms *cl_ptr = &cpu;
/*============================================================================*/

int main( )
{
/* allusions: */
/*
   extern struct tms cpu;
   extern struct tms *cl_ptr;
*/
/* declarations: */

   static FFT 
      fft = {null},
     *fpt = &fft;

   FILE *eingabe;
   FILE *ausgabe;

   char 
     *s_format = "%s\n",
     *i_format = "%ld\n",
     *d_format,
     *d_format_lb;

   char 
      opt[SHS_SIZE] = {null}, /* option: for_backward transformation */ 
     *ptr,
     *ptr1,
     *mode,
     *text,
     *unit_abs,
     *unit_ord,
     *f_unit = "Hz",
     *t_unit = "seconds";

   char **endp = null; 

   char 
      string[TWO][30] = {{ null }},
      file[FTR_NMBR+ONE][50] = {{ null}};

   static short
      ll = null,
      hh = null,
      hrs = null,
      min = null;

   static long 
      ii = null,
      jj = null,
      kk = null,
      mm = null,
      nn = null,
      pp = null,
      qq = null,
      init  = null,
      final = null;  
 
   static long ticks = null;

   time_t nseconds = null;
   time_t *timer   = null;

   static double 
      usr_time = ZERO,
           sec = ZERO;
        
   static double 
      x1 = ZERO,
      x2 = ZERO,
      xx = ZERO,
      dx = ZERO,
      xlower = ZERO,
      xupper = ZERO;

# if READ_REVERSE_ == 1
   static double 
      yy = ZERO;
# endif
/*----------------------------------------------------------------------------*/
/* function prototypes */
/* (1) system functions: */

   long strtol( const char *ptr, char **endp, int base );
   double strtod( const char *ptr, char **endp );

   time_t time( time_t *timer );
   char *ctime( const time_t *timer ); 

   long sysconf( int s_name );
   clock_t times( struct tms *cl_ptr );

/* (2) mathematical functions: */

   double ceil( double x );
   double floor( double x );

/* (3) user defined functions: */

   FFT
     *fftrf( FFT *fpt ),
     *ffnrm( FFT *fpt );
/*----------------------------------------------------------------------------*/
/* memory allocations */

   ptr  = calloc ( STS_SIZE , ONE );
   ptr1 = calloc ( STS_SIZE , ONE ); 
   mode = calloc ( STS_SIZE , ONE );
   text = calloc ( SHS_SIZE , ONE );

   d_format    = calloc( VSS_SIZE, ONE );
   d_format_lb = calloc( VSS_SIZE, ONE );

   unit_abs = calloc ( STS_SIZE , ONE );
   unit_ord = calloc ( STS_SIZE , ONE );

/*----------------------------------------------------------------------------*/

   ticks = sysconf( _SC_CLK_TCK );           /* ticks = clocks per second     */

   setvbuf( stdout, null, _IONBF, null );    /* set buffer length = null      */

   if ( null == strncmp( DBL_FORMAT, "E", ONE ))
   {
      strcpy( d_format,    "%+.15lE   " );
      strcpy( d_format_lb, "%+.15lE \n" );
   }
   else
   {
      strcpy( d_format,    "%+.15le   " );
      strcpy( d_format_lb, "%+.15le \n" );
   };

/*----------------------------- date -----------------------------------------*/
   nseconds = time( timer );
   ptr = ctime( &nseconds );
   printf( "\n\n Program %s started:     %s", FFT_NAME, ptr );

   printf( "====================================="
      "=========================================" );

  menu:  /* mode */

   printf( "\n\n Computation mode:" );
   printf( "\n\n Fourier Transform.........."
      "( time domain to frequency domain ) ---> (1) " );
   printf( "\n Inverse Fourier Transform..( "
      "frequency domain to time domain ) ---> (2) " );
   printf( "\n\n                                  "
      "                                     ) <- ?" ); 
   printf( "\r Please select mode of computation............."
      "[ enter 1 or 2 ] ---> (" );   
   scanf( "%s", ptr );
   ii = strtol( ptr, endp, DEC );

   printf( "\n ------------------------------------"
      "------------------------------------------\n" );
   {
      if ( ii == ONE )   /* time domain ---> freqency domain  */
      {
         strcpy( opt, "forward" );
	 strcpy( mode, "Fourier_transform" );
	 strcpy( unit_abs, f_unit );
         printf( "\n                                    "
            "                                   ) <- ? " );
         printf( "\r compute convolutions ? [ y / n ] "
            "---------------------------------> (" );
	 scanf( "%s", ptr );

	 if ( *ptr == 'y' )
	 {
            printf( "                                     "
               "                                 ......? " );
            printf( "\r enter number of  d i s t i n c t  distributions "
               "[ 1 <= n <= 10 ]--> ");
	    scanf( "%s", ptr1 );  
	    hh = strtol( ptr1, endp, DEC );

            for ( ll=ONE; ll<=hh; ll++ )
	    {
               printf( "                                    "
                  "                                  ......?" );
               printf( "\r please enter filename of distribution %2d "
                  "-------------------------> ", ll );
	       scanf( "%s", file[ll] );
               printf( "                                    "
                  "                                  ......?" );
               printf( "\r please enter multiplicity of distribution %2d "
                  "---------------------> ", ll );
               scanf( "%s", ptr1 );
	       fpt->mult[ll] = strtol( ptr1, endp, DEC );
	    };

            printf( "                                    "
               "                                   ) <- ?" );
            printf( "\r normalize distribution functions ? "
               "[ y / n ] ---------------------> (" );
            scanf( "%s", ptr1 );

	    if ( *ptr1 == 'y' ) 
	    {
               printf( "                                    "
                  "                                  ......?" );
               printf( "\r please enter normalization constant "
                  "( non-neg. integral ) --------> " );
	       scanf( "%s", ptr1 );
	       ( fpt->nor ) = strtod( ptr1, endp );
            }
            else
            {
               ( fpt->nor ) = -1. ;
            };
	 } /* end if ( *ptr == 'y'); compute convolutions */
         else if ( *ptr != 'y' ) /* no convolutions */
	 {
	    hh = ONE;
	    ( fpt->mult[ONE] ) = ONE;
            printf( "\n                                   "
               "                                   ......?" );
            printf( "\r please enter filename of distribution "
               "----------------------------> " );
            scanf( "%s", file[ONE] );
            printf( "                                    "
               "                                   ) <- ?" );
            printf( "\r normalize distribution function ? [ y / n ] "
               "----------------------> (" );
	    scanf( "%s", ptr1 );

	    if ( *ptr1 == 'y' ) 
	    {
               printf( "                                    "
                  "                                  ......?" );
               printf( "\r please enter normalization constant "
                  "( non-neg. integral ) --------> " );
	       scanf( "%s", ptr1 );
	       ( fpt->nor ) = strtod( ptr1, endp );
            }
            else
            {
               ( fpt->nor ) = -1.;
            };
         }; /* end of case *ptr != 'y' [ no convolutions ] */
/*............................................................................*/

         printf( "\n please enter frequency interval "
            "[ s1, s2 ] ( unit: %s ) -------> ", unit_abs );
         printf( "\n ------------------------------"
            "--------------------->  enter s1....: " );
	 scanf( "%s", ptr );
	 xlower = strtod( ptr, endp );
	 printf( " -------------------------------"
            "-------------------->  enter s2....: " );
	 scanf( "%s", ptr );
	 xupper = strtod(ptr,endp);
	 printf( " frequency interval: [ %+.5e , %+.5e ] "
            "( unit: %s )", xlower, xupper, unit_abs );
      } /* end of case ii == ONE: time domain ---> frequency domain */
      else if ( ii != ONE ) /* frequency domain ---> time domain */ 
      {
         strcpy( opt, "backwrd" );
	 strcpy( mode, "inverse_Fourier_transform" );
	 strcpy( unit_abs, t_unit );

	 hh = ONE;
	 ( fpt->mult[ONE] ) = ONE;

         printf( "\n                                   "
            "                                   ......?" );
         printf( "\r please enter filename of spectrum "
            "--------------------------------> " );
	 scanf( "%s", file[ONE] );
         printf( "                                    "
            "                                   ) <- ?" );
         printf( "\r normalize spectrum ?  [ y / n ] "
            "----------------------------------> (" );
	 scanf( "%s", ptr1 );

         if ( *ptr1 == 'y' )
	 {
            printf( "                                     "
               "                                  ......?" );
            printf( "\r please enter normalization constant "
               "( non-neg. integral ) --------> " );
            scanf( "%s", ptr1 );
            ( fpt->nor ) = strtod( ptr1, endp );
         }
         else
         {
            ( fpt->nor ) = -ONE;
         };

	 printf( "\n please enter time interval  [ t1, t2 ] "
            "( unit: %s ) ------------->  ", unit_abs );
         printf( "\n ------------------------------"
            "---------------------->  enter t1...: " );
	 scanf( "%s", ptr );
	 xlower = strtod( ptr, endp );
	 printf( " --------------------------------"
            "-------------------->  enter t2...: " );
	 scanf( "%s", ptr );
	 xupper = strtod( ptr, endp );
	 printf( " time interval:       [ %+.5e , %+.5e ] "
            "( unit: %s ) ", xlower, xupper, unit_abs );
      };
      printf( "\n\n                                  "
         "                                    ......?");
      printf( "\r please enter comment ( max. 30 characters "
         ") ----------------------> " );
      scanf( "%s", text );
   };

   printf( "\n                                   "
      "                                    ) <- ?" );
   printf( "\r input correct ?  [ y / n ] ----"
      "-----------------------------------> (" );
   scanf( "%s", ptr );

   printf( "\n ==================================="
      "===========================================" );

   if (( *ptr != 'y' )&&( *ptr != 'Y' ))
      goto menu;

/*............................................................................*/
/* enter distribution or spectrum: */

   if (( *opt != 'b' )&&( *opt != 'B' ))
      printf( "\n\n entering distributions --> \n" );
   else if (( *opt == 'b' )||( *opt == 'B' )) 
      printf( "\n\n entering spectral distributions --> \n" );

   for ( ll=ONE ; ll<=hh ; ll++)
   {
      printf( "\n --> %d.) %s: \n\n", ll, file[ll] );
      eingabe = fopen( file[ll], "r" );

      fscanf( eingabe, s_format, string[0] );
      printf( "        %s\n", string[0] );
      fscanf( eingabe, s_format, string[1] );
      printf( "        %s\n", string[1] );
      fscanf( eingabe, s_format, unit_abs );
      printf( "\n        abscisse_unit.................:  "
         "%s\n", unit_abs );
      fscanf( eingabe, s_format, unit_ord );
      printf( "        ordinate_unit.................:  "
         "%s\n", unit_ord );

      nn = TWO;

      if (( *opt != 'b' )&&( *opt != 'B' ))
      {
         fscanf( eingabe, s_format, ptr );
   	 ( fpt->t[ll] ) = strtod( ptr, endp );
	 fscanf( eingabe, s_format, ptr );
	 ( fpt->tt[ll] ) = strtod( ptr, endp );
	 fscanf( eingabe, s_format, ptr );
	 ( fpt->dt[ll] ) = strtod( ptr, endp );
	 fscanf( eingabe, s_format, ptr );
         kk = strtol( ptr, endp, DEC );

         printf( "\n        Time interval ....:  [ %+.5e , "
            "%+.5e ]  %s\n", ( fpt->t[ll] ), ( fpt->tt[ll] ), t_unit );
         printf( "        Delta t ..........:  %+.15e     "
            "      %s\n", ( fpt->dt[ll] ), t_unit );
         printf( "        Divisions ........:  %ld\n", kk );

         while (( nn < kk )&&( nn < FTR_SIZE ))
            nn *= TWO;
         
         ( fpt->ttlg[ll] ) = nn;

         strcat( unit_ord, "*" );
         strcat( unit_ord, f_unit );
         strcpy( unit_abs, t_unit );
      }
      else if (( *opt == 'b' )||( *opt == 'B' ))
      {
         fscanf( eingabe, s_format, ptr );
         ( fpt->s[ll] ) = strtod( ptr, endp );
         fscanf( eingabe, s_format, ptr );
         ( fpt->ss[ll] ) = strtod( ptr, endp );
         fscanf( eingabe, s_format, ptr );
         ( fpt->ds[ll] ) = strtod( ptr, endp );
         fscanf( eingabe, s_format, ptr );
         kk = strtol( ptr, endp, DEC );

         printf( "\n        Frequency interval:  [ %+.5e , "
            "%+.5e ]  %s\n", ( fpt->s[ll] ), ( fpt->ss[ll] ), f_unit );
         printf( "        Delta f ..........:  %+.15e     "
            "      %s\n", fpt->ds[ll], f_unit );
         printf( "        Divisions ........:  %ld\n", kk );

         while (( nn < kk )&&( nn < FTR_SIZE ))
            nn *= TWO;

         ( fpt->stlg[ll] ) = nn;
   
         strcat( unit_ord, "*" );
         strcat( unit_ord, t_unit );
         strcpy( unit_abs, f_unit );
      };

      mm = nn / TWO;
      pp = ( nn - kk ) / TWO;
      qq = pp + kk ;

      ii = null;
      while ( ii < pp )
      {
         ( fpt->r[ll][ii] ) = ZERO;
         ( fpt->i[ll][ii] ) = ZERO;
         ii++ ;
      };
      while ( ii < qq )
      {
         fscanf( eingabe, s_format, ptr );
         ( fpt->r[ll][ii] ) = strtod( ptr, endp );
         fscanf( eingabe, s_format, ptr );
         ( fpt->i[ll][ii] ) = strtod( ptr, endp );

# if LIST_FILE == 1
         printf( " %d:   % .16e %+.16e*j \n",
            ii, ( fpt->r[ll][ii] ), ( fpt->i[ll][ii] ));
# endif
         ii++;
      };
      while ( ii < nn )
      {
         ( fpt->r[ll][ii] ) = ZERO;
         ( fpt->i[ll][ii] ) = ZERO;
         ii++ ;
      };

# if READ_REVERSE_ == 1
      ii = pp; do
      {
         jj = ii + mm;
         xx = ( fpt->r[ll][jj] );
         yy = ( fpt->i[ll][jj] );
         ( fpt->r[ll][jj] ) = ( fpt->r[ll][ii] );
         ( fpt->i[ll][jj] ) = ( fpt->i[ll][ii] );
         ( fpt->r[ll][ii] ) = xx;
         ( fpt->i[ll][ii] ) = yy;
      } while (( ++ii ) < mm );
# endif

      fclose( eingabe );
   }; /*  next ll ( next file )  */ 
/*............................................................................*/

   cl_ptr = &cpu;
   times( cl_ptr );
   usr_time = ( cl_ptr -> tms_utime );
/*............................................................................*/

   if ( null <= ( fpt->nor )) 
   {
      ( fpt->p ) = ONE;
      ( fpt->q ) = hh;
      strcpy(( fpt->opt ), opt );
/*............................................................................*/
      fpt = ffnrm( fpt );    /* normalization function */
/*.........................*/

      if (( fpt->rtn ) == ONE )
      {
         printf( "\n\n Error on calling normalization function %s, ",
            DO_FFNRM );
         printf( "\n abnormal program termination.\n " );

         if (( *opt != 'b' )&&( *opt != 'B' ))
            ausgabe  = fopen( spcfle, "w" );
         else if (( *opt == 'b' )||( *opt == 'B' ))
            ausgabe  = fopen( disfle, "w" );

         fprintf( ausgabe, "***** error in subroutine %s ***** ", DO_FFNRM );

         nseconds = time( timer );
         ptr = ctime( &nseconds );

         fprintf( ausgabe, "\n abnormal end of program %s:\n%s",
                  FFT_NAME, ptr );

         fclose( ausgabe );

         goto ende;
      };
   };

   ( fpt->p ) = ONE;
   ( fpt->q ) = hh;
   strcpy(( fpt->opt ), opt );
/*............................................................................*/
   fpt = fftrf( fpt );     /* fast Fourier transformation function */
/*.......................*/

   if (( fpt->rtn ) == ONE )
   {
      printf( "\n\n Error in program %s:", DO_FFTRF );
      printf( "\n abnormal program termination." ); 

      if (( *opt != 'b' )&&( *opt != 'B' ))
         ausgabe  = fopen( spcfle, "w" );
      else if (( *opt == 'b' )||( *opt == 'B' ))
         ausgabe  = fopen( disfle, "w" );

      fprintf( ausgabe, "***** error in subroutine %s ***** ", DO_FFTRF );

      nseconds = time( timer );
      ptr = ctime( &nseconds );

      fprintf( ausgabe, "\n abnormal end of program %s:\n%s", 
         FFT_NAME, ptr );                           
      fclose( ausgabe );

      goto ende;
   };
/*----------------------------------------------------------------------------*/

   if (( *opt != 'b' )&&( *opt != 'B' ))
      ausgabe  = fopen(spcfle,"w");
   else if (( *opt == 'b' )||( *opt == 'B' ))
      ausgabe  = fopen( disfle, "w" );

   fprintf( ausgabe, s_format, mode );
   fprintf( ausgabe, s_format, text );
   fprintf( ausgabe, s_format, unit_abs );
   fprintf( ausgabe, s_format, unit_ord );

   if (( *opt != 'b' )&&( *opt != 'B' ))
   { 
      fprintf( ausgabe, d_format_lb, ( fpt->s[null] ));
      fprintf( ausgabe, d_format_lb, ( fpt->ss[null] ));
      fprintf( ausgabe, d_format_lb, ( fpt->ds[null] ));
      fprintf( ausgabe, i_format, ( fpt->stlg[null] ));

      final = ( fpt->stlg[null] );
   }
   else if (( *opt == 'b')||( *opt == 'B' ))
   {
      fprintf( ausgabe, d_format_lb, ( fpt->t[null] ));
      fprintf( ausgabe, d_format_lb, ( fpt->tt[null] ));
      fprintf( ausgabe, d_format_lb, ( fpt->dt[null] ));
      fprintf( ausgabe, i_format, ( fpt->ttlg[null] ));

      final = ( fpt->ttlg[null] );
   };
   
   ii = null;
   while ( ii < final ) 
   {
      fprintf( ausgabe, d_format, ( fpt->r[null][ii] ));
      fprintf( ausgabe, d_format_lb, ( fpt->i[null][ii] ));

      ii++ ;
   }; 

   nseconds = time( timer );
   ptr = ctime( &nseconds );
   fprintf( ausgabe, "\nProgram %s terminated:\n%s", FFT_NAME, ptr );

   fclose( ausgabe );
/*............................................................................*/
/* create sectional file */

   ausgabe = fopen( outfle ,"w+" );

   fprintf( ausgabe, s_format, mode );
   fprintf( ausgabe, s_format, text );
   fprintf( ausgabe, s_format, unit_abs );
   fprintf( ausgabe, s_format, unit_ord );

   if (( *opt != 'b' )&&( *opt != 'B' ))
   {
       x1 = ( fpt->s[null] );
       x2 = ( fpt->ss[null] );
       dx = ( fpt->ds[null] );
       nn = ( fpt->stlg[null] );
   }
   else if (( *opt == 'b')||( *opt == 'B' ))
   {
       x1 = ( fpt->t[null] );
       x2 = ( fpt->tt[null] );
       dx = ( fpt->dt[null] );
       nn = ( fpt->ttlg[null] );
   };

# if WRITE_REVERSE == 1
   ii = mm;
   xx = x1;
# else
   ii = null;
   xx = ZERO;
# endif

   if( xx < xlower )
      kk = ( long ) floor (( xlower - xx ) / dx );
   else
      kk = null;
   if( nn < kk )
      kk = nn;
      
   init = ii + kk;
   x1 = xx + kk*dx;

   if( xx < xupper )
      kk = ( long ) ceil (( xupper - xx ) / dx );
   else
      kk = null;
   if( nn < kk )
      kk = nn;

   final = ii + kk + ONE;
   x2 = xx + kk*dx;

   kk = final - init;

/* [ Fourier transformations terminated ] */
/*............................................................................*/

   fprintf( ausgabe, d_format_lb, x1 ); 
   fprintf( ausgabe, d_format_lb, x2 );
   fprintf( ausgabe, d_format_lb, dx );
   fprintf( ausgabe, i_format, kk );

   ii = init;
   jj = init;
   while ( ii < final )
   {
      if ( nn <= jj )
          jj -= nn;
      
      fprintf( ausgabe, d_format, ( fpt->r[null][jj] ));
      fprintf( ausgabe, d_format_lb , ( fpt->i[null][jj] ));

      ii++ ;
      jj++ ;
   };

   nseconds = time( timer );
   ptr = ctime( &nseconds );
   fprintf( ausgabe, "\nProgram %s terminated:\n%s", FFT_NAME, ptr );

   fclose( ausgabe );

   printf( "\n\n results on file '%s' ", outfle );

   if (( *opt != 'b' )&&( *opt != 'B' ))
      printf( "\n [ complete spectrum on file %s ] ", spcfle );
   else if (( *opt == 'b')||( *opt == 'B' ))
      printf( "\n [ complete pulse on file %s ] ", disfle );

  ende:

   printf( "\n\n ====================================="
           "=========================================" );
/*............................................................................*/
/* display time and date: */

   nseconds = time( timer );
   ptr = ctime( &nseconds );
   printf( "\n Program %s terminated:\n %s", FFT_NAME, ptr );

/*............................................................................*/
/* system and user time: */

   times( cl_ptr );
   usr_time = ( cl_ptr -> tms_utime ) - usr_time;
   usr_time /= ticks;
   hrs = usr_time/3600;
   min = ( usr_time - hrs*3600 )/60;
   sec = usr_time - hrs*3600 -min*60;

   printf( "The CPU time was: %.2d hours, %.2d "
      "minutes and %+.3e seconds. ", hrs, min, sec ); 
   printf( "\t  Bye - bye !!!\n" );

   exit( EXIT_SUCCESS );
}  
/*============================================================================*/
# undef READ_REVERSE_
# undef WRITE_REVERSE
# undef LIST_FILE
# undef DISP
/************************ end of main program fftrn.c *************************/
