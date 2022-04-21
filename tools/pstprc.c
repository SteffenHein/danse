/* [ file: pstprc.c ] */
# define DO_PSTPRC "pstprc(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function pstprc(*), DANSE release 1.0.                              *
*   [ in DANSE program package ]                                               *
*                                                                              *
*   Computes the S-parameters of a symmetric coupler by evaluating the         *
*   reflection coefficients of all the four quarter systems obtained by        *
*   putting magnetic and electric walls into the two symmetry axes of          *
*   the coupler.                                                               *
*   [ It is thus taken advantage from the equivalence of these quarter         *
*     systems with the coupler under symmetric/antisymmetric [even/odd mode]   *
*     operation conditions.]                                                   *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 27, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
# include <unistd.h>         /* system specification header, cf. sysconf()    */
# include <time.h>           /* cf. time( ),ctime( ),asctime( ),localtime( )  */
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
# include "../poster/POSTER.CONF" /* The general configuration header */
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
   # include <ncurses.h>
# endif 
/*----------------------------------------------------------------------------*/
# include "../poster/posttp.h"
# include "../tools/txctyp.h"
# include "../tools/smithchrt.h"
/*----------------------------------------------------------------------------*/
# define CPL_DIMNS 500           /* maximum number of frequency points */
# define CPL_PRECISION (1.e-12 ) /* relative precision */
# ifndef CPL_CPHASES
   # define CPL_CPHASES  1       /* 1: compute phases */
# endif
# if CPL_CPHASES == 1
   # ifndef CPL_PBRNCH
      # define CPL_PBRNCH 0
   # endif
   # ifndef CPL_PHUNIT
      # define CPL_PHUNIT "DEG"  /* phase unit: "DEG" or "rad" */
   # endif                       /* ["DEG" any char / "rad" lower case char] */
# endif
# ifdef EVL_WRTFREQ
   # define CPL_SVFRQ EVL_WRTFREQ /* 1: read/write frequencies */
# else                            /*    from/into 'spc.xyz' files */ 
   # define CPL_SVFRQ 0
# endif
/*----------------------------------------------------------------------------*/
SPLINES spl = {null};
GRAPHICS gph = {null};
/*----------------------------------------------------------------------------*/
# define CPRODUCT( a, b, c, d, z, w ) \
{ \
   (z) = (a)*(c) - (b)*(d); \
   (w) = (b)*(c) + (a)*(d); \
}
# define CQUOTIENT( a, b, c, d, z, w ) \
{ \
   (w) = (c)*(c) + (d)*(d); \
   (z) = ((a)*(c) + (b)*(d))/(w); \
   (w) = ((b)*(c) - (a)*(d))/(w); \
}
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* 'my_terminal' configuration: */

   # include <termcap.h>     /* terminal type header */
   static char *term;        /* terminal type string */ 

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
/*----------------------------------------------------------------------------*/
# define FLE_CHCK1( ) \
{ \
   if( null == stream ) \
   { \
      printf( "\n\n Error message from function pstprc(*) :" ); \
      printf( "\n Cannot open file %s !!!", fleptr ); \
      printf( "\n [ Please verify if file exists " \
         "in presently working directory.]\n " ); \
      exit( EXIT_FAILURE ); \
   }; \
}
/*----------------------------------------------------------------------------*/
# define FLE_CHCK2( ) \
{ \
   ind = null; \
   if( null == strstr( ptr, eop )) \
   { \
      printf( "\n Sorry, there are no %s mode S-parameters " \
         "in this file !!!", eop ); \
      printf( "\n [ Please re-enter file name.]\n\n " ); \
      ind = ONE; \
   }; \
}
/*----------------------------------------------------------------------------*/
# define CPL_SPLCPY( opt, ipl ) \
{ \
   spt->vct[ii][null] = ss; \
   if ((opt) == 'p' ) \
      spt->vct[ii][ONE] = 100.*sqrt( real1*real1 + imag1*imag1 ); \
   else if ((opt) == 'd' ) \
      spt->vct[ii][ONE] = 10.*log10( real1*real1 + imag1*imag1 ); \
   else if ((opt) == 'r' ) \
      spt->vct[ii][ONE] = real1; \
   ind = null; \
   do \
   { \
      spt->dmn[nn] = ss; \
      ss += ( dx/((ipl) + ONE )); \
      ind++ ; \
      nn++ ; \
   } while(( ind <= (ipl))&&( ii < ( mm - ONE ))); \
} \
/*----------------------------------------------------------------------------*/
# if CPL_CPHASES == 1
   # include "../math/argc.h"
# endif
/*============================================================================*/

POSTSTATE *\
pstprc( POSTSTATE *stp )
{
/* declarations: */

   static POSTSTATE
     *state;

   static FILE 
     *stream;

   static GRAPHICS
      *gpt = &gph;

   static SPLINES
      *spt = &spl;

   static TXCNSL
     *csp;

   static time_t 
      nseconds = null;

   static short
      ind = null;

   static long
      ii = null,
      nn = null,
      mm = null;

   static double
      xx1 = ZERO,
      xx2 = ZERO,
      dx = ZERO,
      real1 = ZERO,
      real2 = ZERO,
      imag1 = ZERO,
      imag2 = ZERO,
      ss = ZERO,
      tt = ZERO;

   static const char
      items = SEVEN,
     *hrzev = "horz.even",
     *hrzod = "horz.odd_",
     *vrtev = "vert.even",
     *vrtod = "vert.odd_";

# if CPL_CPHASES == 1
   static const char
     *phunit = "radians";
# endif

   static char
      option = null,
      fleptr[STS_SIZE] = {null},
      tlmmod[STS_SIZE] = {null},
      xunit[STS_SIZE] = {null},
      zunit[STS_SIZE] = {null},
      idx[SHS_SIZE] = {null},
      ptr[STS_SIZE] = {null},
      eop[STS_SIZE] = {null},
     *tmeptr,
    **endp = null;

   static double
      scpl11[CPL_DIMNS][TWO] = {{ZERO}},
      scpl21[CPL_DIMNS][TWO] = {{ZERO}},
      scpl31[CPL_DIMNS][TWO] = {{ZERO}},
      scpl41[CPL_DIMNS][TWO] = {{ZERO}},
      direct[CPL_DIMNS][TWO] = {{ZERO}},
      see11[CPL_DIMNS][TWO] = {{ZERO}},
      seo11[CPL_DIMNS][TWO] = {{ZERO}},
      soe11[CPL_DIMNS][TWO] = {{ZERO}},
      soo11[CPL_DIMNS][TWO] = {{ZERO}};

# if CPL_CPHASES == 1
   static COMPLEX cc = {null},
                *cpt = &cc;

   static double
      phs11[CPL_DIMNS] = {ZERO},
      phs21[CPL_DIMNS] = {ZERO},
      phs31[CPL_DIMNS] = {ZERO},
      phs41[CPL_DIMNS] = {ZERO};
# endif

/* prototypes: */

# ifndef _CCBUG
   int strcmp( const char *s1, const char *s2 );
   char *strstr( const char *s1, const char *s2 );
# endif

   time_t time( time_t *timer );
   char *ctime( const time_t *timer );

   SPLINES
     *spline( SPLINES *spt );

   int
     graphp( GRAPHICS *gpt );

   TXCNSL 
     *txcnsl( TXCNSL *csp );

   char 
     *lotos( long n, char m );

# if CPL_CPHASES == 1
   COMPLEX
     *argc( COMPLEX *ipt, short nn );
# endif

   POSTSTATE
     *smithchrt( POSTSTATE *stp );
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
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
   state = stp;

   tmeptr = ( char * ) calloc( STS_SIZE, ONE );

/*............................................................................*/
/* initialize menu: */

   csp = txcnsl( null );

   ( csp->dfopt ) = items - ONE;

  menu:

   ( csp->items ) = items;
   ( csp->dflnf ) = 1; /* set line feed before and after default menu option */

   strcpy( csp->tasks, "COMPUTE ... " );
   strcpy( csp->cmmnt, "Even/odd mode postprocessing:" );
   strcpy( csp->envmt, "PSTPRC" );

   strcpy( csp->mline[1], "* reflection coefficient >--> [ s11 ]" );
   strcpy( csp->mline[2], "* coupling coefficient >----> [ s21 ]" );
   strcpy( csp->mline[3], "* coupling coefficient >----> [ s31 ]" );
   strcpy( csp->mline[4], "* coupling coefficient >----> [ s41 ]" );
   strcpy( csp->mline[5], "* directivity >---------> [ s41/s21 ]" );
   strcpy( csp->mline[6], "* Compute all S-parameters" );
   strcpy( csp->mline[7], "* Generate Smith chart" );

   strcpy( csp->escpe, "End of program / escape:" );

/*............................................................................*/
   csp = txcnsl( csp );    /*                                                 */
/*.......................*/
   
   option = ( csp->option );

/*............................................................................*/
   if ( option == items )
   {
      PRBLDCLR( "\r" );
      printf( " \t\t\t\t\t\t\t\t\t PSTPRC" );
      PRNORMAL( "" );

      state = smithchrt( state );
      goto menu;
   }
/*............................................................................*/
/* s_[horizontal_even/vertical_even]: */

   else if(( null < option )&&( option < items ))
   {
      printf( "\n" );

     even_even:

      strcpy( eop, hrzev ); 
      strcat( eop, "/" );
      strcat( eop, vrtev );

      strcpy( csp->rqlng, "Enter job index of '" );
      strcat( csp->rqlng, eop );
      strcat( csp->rqlng, "' mode S-prmtrs" );
      strcpy( csp->rqfrm, "brackets" );
      ( csp->dflng ) = 1;
/*............................................................................*/
      csp = txcnsl( csp );           /* enter index N on text console   */
/*.................................*/
      strcpy( idx, lotos( csp->inlng, null ));
      strcpy ( fleptr, "spc.s11_" );
      strcat ( fleptr, idx );

      stream = fopen( fleptr, "r" );

      FLE_CHCK1( );

      fscanf( stream, "%s", tlmmod ); /* DSC system identifier */
      fscanf( stream, "%s", ptr );    /* comment */

      FLE_CHCK2( ); /* check even/even mode identifier */

      if ( ind == ONE )
         goto even_even;

      fscanf( stream, "%s", xunit ); 

      if ( null == strcmp( xunit, "1/seconds" ))
         strcpy( xunit, "Hz" );

      fscanf( stream, "%s", zunit );

      fscanf( stream, "%s", ptr );
      xx1 = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      xx2 = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      dx = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );

      mm = strtol( ptr, endp, DEC );

      if( CPL_DIMNS < mm )
      {
         printf( "\n\n Error message from function %s:", DO_PSTPRC );
         printf( "\n Too many frequency points !!!" );
         printf( "\n [ The maximum number is %ld = macro CPL_DIMNS",
            (long) CPL_DIMNS );
         printf( "\n   - change macro in SCPLCONF.H ]\n " );
         exit( EXIT_FAILURE );
      };

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fscanf( stream, "%s", ptr );
# endif
         fscanf( stream, "%s", ptr );
         see11[ii][null] = strtod( ptr, endp );
         fscanf( stream, "%s", ptr );
         see11[ii][ONE] = strtod( ptr, endp );
         ii++ ;
      };

      fclose( stream );
   
/*............................................................................*/
/* s_[horizontal_even/vertical_odd]: */

     even_odd:

      strcpy( eop, hrzev ); 
      strcat( eop, "/" );
      strcat( eop, vrtod );

      strcpy( csp->rqlng, "Enter job index of '" );
      strcat( csp->rqlng, eop );
      strcat( csp->rqlng, "' mode S-prmtrs" );
      strcpy( csp->rqfrm, "brackets" );
      ( csp->dflng ) = 2;
/*............................................................................*/
      csp = txcnsl( csp );           /* enter index N on text console   */
/*.................................*/
      strcpy( idx, lotos( csp->inlng, null ));
      strcpy ( fleptr, "spc.s11_" );
      strcat ( fleptr, idx );
   
      stream = fopen( fleptr, "r" );

      FLE_CHCK1( );

      fscanf( stream, "%s", ptr ); /* DSC system identifier */
      fscanf( stream, "%s", ptr ); /* comment */

      FLE_CHCK2( );
 
      if( ind == ONE )
         goto even_odd;

      fscanf( stream, "%s", ptr );
      fscanf( stream, "%s", ptr );
      fscanf( stream, "%s", ptr );

      xx1 = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      xx2 = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      dx = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );

      nn = strtol( ptr, endp, DEC );
      if( nn != mm )
      {
         printf( "\n\n Error message from function %s:", DO_PSTPRC );
         printf( "\n Incompatible S-parameter files s11[even/even], "
            "s11[even/odd] !!!" );
         printf( "\n [ The number of frequency points differs: %ld != %ld ]",
            mm, nn );
         printf( "\n " );
         exit( EXIT_FAILURE );
      };

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fscanf( stream, "%s", ptr );
# endif
         fscanf( stream, "%s", ptr );
         seo11[ii][null] = strtod( ptr, endp );
         fscanf( stream, "%s", ptr );
         seo11[ii][ONE] = strtod( ptr, endp );
         ii++ ;
      };

      fclose( stream );

/*............................................................................*/
/* s_[horizontal_odd/vertical_even]: */

     odd_even:

      strcpy( eop, hrzod ); 
      strcat( eop, "/" );
      strcat( eop, vrtev );

      strcpy( csp->rqlng, "Enter job index of '" );
      strcat( csp->rqlng, eop );
      strcat( csp->rqlng, "' mode S-prmtrs" );
      strcpy( csp->rqfrm, "brackets" );
      ( csp->dflng ) = 3;
/*............................................................................*/
      csp = txcnsl( csp );           /* enter index N on text console   */
/*.................................*/
      strcpy( idx, lotos( csp->inlng, null ));
      strcpy ( fleptr, "spc.s11_" );
      strcat ( fleptr, idx );
   
      stream = fopen( fleptr, "r" );

      FLE_CHCK1( );

      fscanf( stream, "%s", ptr ); /* DSC system identifier */
      fscanf( stream, "%s", ptr ); /* comment */

      FLE_CHCK2( );

      if( ind == ONE )
         goto odd_even;

      fscanf( stream, "%s", ptr );
      fscanf( stream, "%s", ptr );
      fscanf( stream, "%s", ptr );

      xx1 = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      xx2 = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      dx = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );

      nn = strtol( ptr, endp, DEC );
      if( nn != mm )
      {
         printf( "\n\n Error message from function %s:", DO_PSTPRC );
         printf( "\n Incompatible S-parameter files s11[even/even], "
            "s11[odd/even] !!!" );
         printf( "\n [ The number of frequency points differs: %ld != %ld ]",
            mm, nn );
         printf( "\n " );
         exit( EXIT_FAILURE );
      };

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fscanf( stream, "%s", ptr );
# endif
         fscanf( stream, "%s", ptr );
         soe11[ii][null] = strtod( ptr, endp );
         fscanf( stream, "%s", ptr );
         soe11[ii][ONE] = strtod( ptr, endp );
         ii++ ;
      };

      fclose( stream );

/*............................................................................*/
/* s_[horizontal_odd/vertical_even]: */

     odd_odd:

      strcpy( eop, hrzod ); 
      strcat( eop, "/" );
      strcat( eop, vrtod );

      strcpy( csp->rqlng, "Enter job index of '" );
      strcat( csp->rqlng, eop );
      strcat( csp->rqlng, "' mode S-prmtrs" );
      strcpy( csp->rqfrm, "brackets" );
      ( csp->dflng ) = 4;
/*............................................................................*/
      csp = txcnsl( csp );           /* enter index N on text console   */
/*.................................*/
      strcpy( idx, lotos( csp->inlng, null ));
      strcpy ( fleptr, "spc.s11_" );
      strcat ( fleptr, idx );
   
      stream = fopen( fleptr, "r" );

      FLE_CHCK1( );

      fscanf( stream, "%s", ptr ); /* DSC system identifier */
      fscanf( stream, "%s", ptr ); /* comment */

      FLE_CHCK2( );

      if( ind == ONE )
         goto odd_odd;

      fscanf( stream, "%s", ptr );
      fscanf( stream, "%s", ptr );
      fscanf( stream, "%s", ptr );

      xx1 = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      xx2 = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      dx = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );

      nn = strtol( ptr, endp, DEC );
      if( nn != mm )
      {
         printf( "\n\n Error message from function %s:", DO_PSTPRC );
         printf( "\n Incompatible S-parameter files s11[even/even], "
            "s11[odd/odd] !!!" );
         printf( "\n [ The number of frequency points differs: %ld != %ld ]",
            mm, nn );
         printf( "\n " );
         exit( EXIT_FAILURE );
      };

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fscanf( stream, "%s", ptr );
# endif
         fscanf( stream, "%s", ptr );
         soo11[ii][null] = strtod( ptr, endp );
         fscanf( stream, "%s", ptr );
         soo11[ii][ONE] = strtod( ptr, endp );
         ii++ ;
      };

      fclose( stream );
   }; /* end if (( null < option )&&( option < items )) */

/*............................................................................*/
   if(( option == 1 )||( option == ( items-ONE )))
   {
/* reflection: */

      ii = null;
      while( ii < mm )
      {
         scpl11[ii][null] = ( see11[ii][null] + seo11[ii][null] +\
                              soe11[ii][null] + soo11[ii][null] )/4.;

         scpl11[ii][ONE] = ( see11[ii][ONE] + seo11[ii][ONE] +\
                             soe11[ii][ONE] + soo11[ii][ONE] )/4.;
# if CPL_CPHASES == 1
         ( cpt->r ) = scpl11[ii][null];
         ( cpt->i ) = scpl11[ii][ONE];
/*............................................................................*/
         cpt = argc( cpt, CPL_PBRNCH );      /* function arg(z) */
/*.........................................*//* [ branch CPL_PBRNCH ] */
         if ( null == strstr( phunit, CPL_PHUNIT ))
            phs11[ii] = 180.*( cpt->arg )/PI;
         else
            phs11[ii] = cpt->arg;
# endif
         ii++ ;
      };

      strcpy( fleptr, "cpl.ref" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "s11" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", zunit );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

# if CPL_SVFRQ == 1
      ss = xx1;
# endif

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fprintf( stream, "\n%+.15e", ss );
	 ss += dx;
         fprintf( stream, "\t %+.15e", scpl11[ii][null] );
# else
         fprintf( stream, "\n%+.15e", scpl11[ii][null] );
# endif
         fprintf( stream, "\t %+.15e", scpl11[ii][ONE] );
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                       fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "EOF\n");
      fclose( stream );

# if CPL_CPHASES == 1
      strcpy( fleptr, "cpl.phs11" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "phase_s11" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", CPL_PHUNIT );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

      ss = xx1;
      ii = null;
      while( ii < mm )
      {
         fprintf( stream, "\n%+.15e", ss );
         fprintf( stream, "   %+.15e", phs11[ii] );
         ss += dx;
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                       fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "EOF\n" );
      fclose( stream );
# endif

/*............................................................................*/
/* graphics, reflection [percent]: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = scpl11[ii][null];
         imag1 = scpl11[ii][ONE];

         CPL_SPLCPY( 'p', EVL_GPHINTP ); /* input: spline(*) function */
                                          /* [ opt, 'p'ercent ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for interpolating reflection !!!" );
         printf( "\n [ program stopped.]\n " );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "ref" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "reflection_[0/0]" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, "0/0" );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   
/*............................................................................*/
/* graphics, reflection [dB]: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = scpl11[ii][null];
         imag1 = scpl11[ii][ONE];

         CPL_SPLCPY( 'd', EVL_GPHINTP ); /* input: spline(*) function */
                                          /* [ opt, 'd'ecibel ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for s11 interpolation !!!" );
         printf( "\n [ program stopped.]\n " );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "s11" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "s11" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, "dB" );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   

# if CPL_CPHASES == 1 
/*............................................................................*/
/* graphics, phase s11: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = phs11[ii];

         CPL_SPLCPY( 'r', null ); /* input: spline(*) function */
                                  /* [ opt, 'r'eal ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for s11 phase interpolation !!!" );
         printf( "\n [ program stopped.]\n " );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "phs11" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "phase_s11" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, CPL_PHUNIT );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   
# endif

   };/* end if ( option == ... ) */
/*............................................................................*/
   if(( option == 2 )||( option == ( items-ONE )))
   {
/* transmission s21: */

      ii = null;
      while( ii < mm )
      {
         scpl21[ii][null] = ( see11[ii][null] - seo11[ii][null] +\
                              soe11[ii][null] - soo11[ii][null] )/4.;

         scpl21[ii][ONE] = ( see11[ii][ONE] - seo11[ii][ONE] +\
                             soe11[ii][ONE] - soo11[ii][ONE] )/4.;
# if CPL_CPHASES == 1
         ( cpt->r ) = scpl21[ii][null];
         ( cpt->i ) = scpl21[ii][ONE];
/*............................................................................*/
         cpt = argc( cpt, CPL_PBRNCH );      /* function arg(z) */
/*.........................................*//* [ branch CPL_PBRNCH ] */
         if ( null == strstr( phunit, CPL_PHUNIT ))
            phs21[ii] = 180.*( cpt->arg )/PI;
         else
            phs21[ii] = cpt->arg;
# endif
         ii++ ;
      };

      strcpy( fleptr, "cpl.s21" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "s21" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", zunit );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

# if CPL_SVFRQ == 1
      ss = xx1;
# endif 

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fprintf( stream, "\n%+.15e", ss );
	 ss += dx;
         fprintf( stream, "\t %+.15e", scpl21[ii][null] );
# else
         fprintf( stream, "\n%+.15e", scpl21[ii][null] );
# endif
         fprintf( stream, "\t %+.15e", scpl21[ii][ONE] );
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                        fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "EOF\n" );

      fclose( stream );

# if CPL_CPHASES == 1
      strcpy( fleptr, "cpl.phs21" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "phase_s21" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", CPL_PHUNIT );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

      ss = xx1;
      ii = null;
      while( ii < mm )
      {
         fprintf( stream, "\n%+.15e", ss );
         fprintf( stream, "   %+.15e", phs21[ii] );
         ss += dx;
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                       fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "EOF\n" );
      fclose( stream );
# endif

/*............................................................................*/
/* graphics s21 [dB]: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = scpl21[ii][null];
         imag1 = scpl21[ii][ONE];

         CPL_SPLCPY( 'd', EVL_GPHINTP ); /* input: spline(*) function */
                                          /* [ opt, 'd'ecibel ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for s21 interpolation !!!" );
         printf( "\n [ program stopped.]\n" );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "s21" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "coupling_s21" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, "dB" );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   

# if CPL_CPHASES == 1
/*............................................................................*/
/* graphics, phase s21: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = phs21[ii];

         CPL_SPLCPY( 'r', null ); /* input: spline(*) function */
                                  /* [ opt, 'r'eal ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for s21 phase interpolation !!!" );
         printf( "\n [ program stopped.]\n" );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "phs21" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "phase_s21" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, CPL_PHUNIT );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   
# endif

   };/* end if ( option == ... ) */
/*............................................................................*/
   if(( option == 3 )||( option == ( items-ONE )))
   {
/* coupling s31: */

      ii = null;
      while( ii < mm )
      {
         scpl31[ii][null] = ( see11[ii][null] + seo11[ii][null] -\
                              soe11[ii][null] - soo11[ii][null] )/4.;

         scpl31[ii][ONE] = ( see11[ii][ONE] + seo11[ii][ONE] -\
                             soe11[ii][ONE] - soo11[ii][ONE] )/4.;
# if CPL_CPHASES == 1
         ( cpt->r ) = scpl31[ii][null];
         ( cpt->i ) = scpl31[ii][ONE];
/*............................................................................*/
         cpt = argc( cpt, CPL_PBRNCH );      /* function arg(z) */
/*.........................................*//* [ branch CPL_PBRNCH ] */
         if ( null == strstr( phunit, CPL_PHUNIT ))
            phs31[ii] = 180.*( cpt->arg )/PI;
         else
            phs31[ii] = cpt->arg;
# endif
         ii++ ;
      };

      strcpy( fleptr, "cpl.s31" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "s31" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", zunit );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

# if CPL_SVFRQ == 1
      ss = xx1;
# endif

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fprintf( stream, "\n%+.15e", ss );
	 ss += dx;
         fprintf( stream, "\t %+.15e", scpl31[ii][null] );
# else
         fprintf( stream, "\n%+.15e", scpl31[ii][null] );
# endif
         fprintf( stream, "\t %+.15e", scpl31[ii][ONE] );
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                        fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "EOF\n" );

      fclose( stream );

# if CPL_CPHASES == 1
      strcpy( fleptr, "cpl.phs31" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "phase_s31" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", CPL_PHUNIT );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

      ss = xx1;
      ii = null;
      while( ii < mm )
      {
         fprintf( stream, "\n%+.15e", ss );
         fprintf( stream, "   %+.15e", phs31[ii] );
         ss += dx;
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                       fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "EOF\n" );
      fclose( stream );
# endif

/*............................................................................*/
/* graphics s31 [dB]: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = scpl31[ii][null];
         imag1 = scpl31[ii][ONE];

         CPL_SPLCPY( 'd', EVL_GPHINTP ); /* input: spline(*) function */
                                          /* [ opt, 'd'ecibel ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for s31 interpolation !!!" );
         printf( "\n [ program stopped.]\n" );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "s31" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "coupling_s31" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, "dB" );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   

# if CPL_CPHASES == 1
/*............................................................................*/
/* graphics, phase s31: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = phs31[ii];

         CPL_SPLCPY( 'r', null ); /* input: spline(*) function */
                                  /* [ opt, 'r'eal ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for s31 phase interpolation !!!" );
         printf( "\n [ program stopped.]\n" );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "phs31" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "phase_s31" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, CPL_PHUNIT );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   
# endif

   };/* end if ( option == ... ) */
/*............................................................................*/
   if(( option == 4 )||( option == ( items-ONE )))
   {
/* coupling s41: */

      ii = null;
      while( ii < mm )
      {
         scpl41[ii][null] = ( see11[ii][null] - seo11[ii][null] -\
                              soe11[ii][null] + soo11[ii][null] )/4.;

         scpl41[ii][ONE] = ( see11[ii][ONE] - seo11[ii][ONE] -\
                             soe11[ii][ONE] + soo11[ii][ONE] )/4.;
# if CPL_CPHASES == 1
         ( cpt->r ) = scpl41[ii][null];
         ( cpt->i ) = scpl41[ii][ONE];
/*............................................................................*/
         cpt = argc( cpt, CPL_PBRNCH );      /* function arg(z) */
/*.........................................*//* [ branch CPL_PBRNCH ] */
         if ( null == strstr( phunit, CPL_PHUNIT ))
            phs41[ii] = 180.*( cpt->arg )/PI;
         else
            phs41[ii] = cpt->arg;
# endif
         ii++ ;
      };

      strcpy( fleptr, "cpl.s41" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "s41" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", zunit );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

# if CPL_SVFRQ == 1
      ss = xx1;
# endif

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fprintf( stream, "\n%+.15e", ss );
	 ss += dx;
         fprintf( stream, "\t %+.15e", scpl41[ii][null] );
# else
         fprintf( stream, "\n%+.15e", scpl41[ii][null] );
# endif
         fprintf( stream, "\t %+.15e", scpl41[ii][ONE] );
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                        fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "EOF\n" );

      fclose( stream );

# if CPL_CPHASES == 1
      strcpy( fleptr, "cpl.phs41" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "phase_s41" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", CPL_PHUNIT );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

      ss = xx1;
      ii = null;
      while( ii < mm )
      {
         fprintf( stream, "\n%+.15e", ss );
         fprintf( stream, "   %+.15e", phs41[ii] );
         ss += dx;
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                       fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "EOF\n" );
      fclose( stream );
# endif

/*............................................................................*/
/* graphics s41 [dB]: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = scpl41[ii][null];
         imag1 = scpl41[ii][ONE];

         CPL_SPLCPY( 'd', EVL_GPHINTP ); /* input: spline(*) function */
                                          /* [ opt, 'd'ecibel ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for interpolating reflection !!!" );
         printf( "\n [ program stopped.]\n" );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "s41" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "coupling_s41" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, "dB" );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   

# if CPL_CPHASES == 1
/*............................................................................*/
/* graphics, phase s41: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = phs41[ii];

         CPL_SPLCPY( 'r', null ); /* input: spline(*) function */
                                  /* [ opt, 'r'eal ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for s41 phase interpolation !!!" );
         printf( "\n [ program stopped.]\n" );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "phs41" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "phase_s41" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, CPL_PHUNIT );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   
# endif

   };/* end if ( option == ... ) */
/*............................................................................*/
   if(( option == 5 )||( option == ( items-ONE )))
   {
/* directivity: */

      ii = null;
      while( ii < mm )
      {
         real1 = ( see11[ii][null] + seo11[ii][null] -\
                   soe11[ii][null] - soo11[ii][null] )/4.;

         imag1 = ( see11[ii][ONE] + seo11[ii][ONE] -\
                   soe11[ii][ONE] - soo11[ii][ONE] )/4.;

         real2 = ( see11[ii][null] - seo11[ii][null] -\
                   soe11[ii][null] + soo11[ii][null] )/4.;

         imag2 = ( see11[ii][ONE] - seo11[ii][ONE] -\
                   soe11[ii][ONE] + soo11[ii][ONE] )/4.;

         CQUOTIENT( real1, imag1, real2, imag2, ss, tt ); /* s31/s41 */

         direct[ii][null] = ss;
         direct[ii][ONE] = tt;
         ii++ ;
      };

      strcpy( fleptr, "cpl.dir" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "s31/s41" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", zunit );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

# if CPL_SVFRQ == 1
      ss = xx1;
# endif

      ii = null;
      while( ii < mm )
      {

# if CPL_SVFRQ == 1
         fprintf( stream, "\n%+.15e", ss );
	 ss += dx;
         fprintf( stream, "\t %+.15e", direct[ii][null] );
# else
         fprintf( stream, "\n%+.15e", direct[ii][null] );
# endif
         fprintf( stream, "\t %+.15e", direct[ii][ONE] );
         ii++ ;
      };

      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n%s%s%s\n%s", "ASCII file '",
                       fleptr, "' created: ", tmeptr );
      fprintf( stream, "\n%s", "                                     \n");

      fclose( stream );

/*............................................................................*/
/* graphics, directivity [dB]: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = direct[ii][null];
         imag1 = direct[ii][ONE];

         CPL_SPLCPY( 'd', EVL_GPHINTP ); /* input: spline(*) function */
                                          /* [ opt, 'd'ecibel ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Graphics require too many points !!!" );
         printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
            "in funtion %s.", (long) nn, (long) GPH_POINTS, "graphp(*)" );
         printf( "\n [ Change macro in compliance "
            "with memory resources.]\n " );
         exit( EXIT_FAILURE );
      };

      ( spt->mm ) = mm; /* number of support points */
      ( spt->nn ) = nn; /* number of interpolated points */
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: reflection           */
/*...........................*//*                                             */
      if ( ( spt->rtn ) != null )                                      
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n Error on calling spline function" );
         printf( "\n for interpolating reflection !!!" );
         printf( "\n [ program stopped.]\n" );
         exit( EXIT_FAILURE );
      }
      else
      {
         ii = null;
         while ( ii < nn )
         {
            gpt->vct[ii][null] = spt->dmn[ii];
            gpt->vct[ii][ONE]  = spt->fct[ii];
            ii++ ;
         };

         strcpy( gpt->format, PLOT_FORMAT );
         strcpy( gpt->file, "dir" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "directivity_s41/s21" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, "dB" );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   
   };/* end if ( option == ... ) */

   if ( option != null )
      goto menu;

   PRBLDCLR( "\r" );
   printf( " \t\t\t\t\t\t\t\t\t PSTPRC" );
   PRNORMAL( "" );

   ( state->rtn ) = null;
   return state;
}
/*============================================================================*/
# undef CPL_PRECISION
# undef CPL_CPHASES
# undef CPL_PBRNCH
# undef CPL_PHUNIT
# undef CPL_DIMNS
# undef CPL_SVFRQ
/***************** end of postprocessing function pstprc(*) *******************/
