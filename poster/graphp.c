/* [ file: graphp.c ] */
# define DO_GRAPHP "graphp(*)"
/*******************************************************************************
*                                                                              *
*   Graphics files creation function graphp(*)                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 08, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <time.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# include "../CONFIG.H"
# include "./POSTER.CONF"
/*----------------------------------------------------------------------------*/
# include "./posttp.h"
/*--------------------------------------------------------------------------*/
# ifndef PLOT_FORMAT
   # define PLOT_FORMAT "gnuplot_2D" /*plot file format: "SPLINE" or "gnuplot"*/
# endif
# ifndef GPH_USCALE
   # define GPH_USCALE 0    /* 1: take uniform scales on xy axes */
# endif
# ifndef GPH_CURVES
   # define GPH_CURVES 20   /* maximum number of simultaneously plotted */
# endif                     /* graphics */
# ifndef GPH_POINTS
   # define GPH_POINTS 5000 /* maximum number of graphics points */
# endif
# ifndef GPH_DATFLE
   # define GPH_DATFLE "dat"
# endif
# ifndef GPH_RNDOFF         /* roundoff to that number of floating pnt digits */
   # define GPH_RNDOFF 12
# endif
/* The unscaled gnuplot range: abs( log10(x) ) < GPH_UNSCALED [ i.e. beyond */
/* the interval ( 1.e-GPH_UNSCALED, 1.e+GPH_UNSCALED ) values are rescaled ] */
# ifndef GPH_UNSCALED
   # define GPH_UNSCALED 9
# endif                  
/*----------------------------------------------------------------------------*/
/* macros: */

# define GPH_FOPEN( pfx, fle ) \
{ \
   strcpy( pptr, (pfx) ); \
   strcat( pptr, (fle) ); \
 \
   pltfle = fopen( pptr, "w" ); \
 \
   if (( pltfle == null ) \
     &&( null != ( gpt->dsp ))) \
   { \
      fprintf( stderr, "\n\n Message from function %s:", DO_GRAPHP ); \
      fprintf( stderr, "\n\n Error on opening file %s !!!", pptr ); \
      fprintf( stderr, "\n [ overrides and returns to calling program ].\n "); \
      return null; \
   }; \
}
/*--------------------------------------------------------------------------*/
# ifndef TYPE_GRAPHICS
# define TYPE_GRAPHICS 1
/* graphics data transfer structure [ used in function graphp(*), e.g.]: */
typedef struct
{
   signed char
      rtn, /* return operation mark: 0: returm with error */
      dsp; /* display operation mark: 1 display some file saving messages */

   char
      name[STS_SIZE],
      text[STS_SIZE];

   char
      file[STS_SIZE],
      curve[GPH_CURVES][STS_SIZE],
      format[SHS_SIZE];

   char
      xunit[SHS_SIZE],
      yunit[SHS_SIZE];

   short
      nc; /* nc = number of graphics [ 'curves' ] */

   long
      nn,
      np[GPH_CURVES]; /* np = number of sample points */

   double
      xmin,
      xmax,
      ymin,
      ymax,
      vct[GPH_POINTS][GPH_CURVES+ONE];

} GRAPHICS;
# endif
/*============================================================================*/

int graphp( GRAPHICS *gpt )
{
/* allusions: */

/* declarations; */

   time_t nseconds = null;
   time_t   *timer = null;

   static short 
      ii = null,
      lx = null,
      ly = null;

   static long 
      kk = null;

   static double
      x_min   = ZERO,
      y_min   = ZERO,
      x_max   = ZERO,
      y_max   = ZERO,
      x_mean  = ZERO, 
      y_mean  = ZERO,
      x_shift = ZERO,
      y_shift = ZERO,
      xx      = ZERO,
      yy      = ZERO,
      uu      = ZERO,
      vv      = ZERO,
      rscale  = .95;  /* range scaling factor */

   static char
     *timefrm = " created: %.24s ",
     *timeptr,
     *pptr,
     *dptr,
     *prefix;

# if GPH_USCALE == 1
   static char 
      **endp = null; 
# endif

/* prototypes: */

   double fabs( double x );
   double pow( double x, double y );
   double log( double x );
   double log10( double x );
   double floor( double x );

   time_t time( time_t *timer );

/* streams: */
       
   FILE 
     *pltfle,
     *datfle;

   char *lotos ( long ll, char digits );

# ifdef GPH_RNDOFF
   double rndoff( double xx, short nn );
# endif
/*----------------------------------------------------------------------------*/
/* memory allocations: */

   pptr    = ( char *) calloc( STS_SIZE, ONE );
   dptr    = ( char *) calloc( STS_SIZE, ONE );
   prefix  = ( char *) calloc( VSS_SIZE, ONE );
   timeptr = ( char *) calloc( STS_SIZE, ONE );

   x_shift = ZERO;
   y_shift = ZERO;
   rscale  = .95;
/*............................................................................*/
   if (( gpt->nc ) == null ) /* no curve number transferred */
      ( gpt->nc ) = ONE;
   else if ( GPH_CURVES < ( gpt->nc ))
   {
      fprintf( stderr, "\n\n Error in function %s:", DO_GRAPHP );
      fprintf( stderr, "\n\n Too many [=%d] curves transferred to function !!!",
         ( gpt->nc ));
      fprintf( stderr, "\n [ The maximum number is %ld = macro GPH_CURVES "
         "in function %s.", ( long ) GPH_CURVES, DO_GRAPHP );
      fprintf( stderr, "\n - Change macro only in compliance "
         "with memory resources.]\n" );

      exit( EXIT_FAILURE );
   };
/*............................................................................*/
   ii = null;
   while( ii < ( gpt->nc ))
   {
      if (( gpt->nn ) < ( gpt->np[ii] ))
         ( gpt->nn ) = ( gpt->np[ii] );
      ii++;
   };

   ii = null;
   while( ii < ( gpt->nc ))
   {
      if (( gpt->np[ii] ) == null ) /* no individual point number transferred */
         ( gpt->np[ii] ) = ( gpt->nn );
      ii++;
   };    
/*............................................................................*/
   ii = null;
   while ( ii < ( gpt->nc ))
   {
      if ( GPH_POINTS < ( gpt->np[ii] ))
      {
         fprintf( stderr, "\n\n Error in function %s:", DO_GRAPHP );
         fprintf( stderr,
            "\n\n Too many [=%ld] points transferred in array vct[][%d] !!!",
            ( long )( gpt->np[ii] ), ( ii + ONE ));
         fprintf( stderr, "\n [ The maximum number is %ld = macro GPH_POINTS "
            "in function %s.", (long) GPH_POINTS, DO_GRAPHP );
         fprintf( stderr, "\n - Change macro only in compliance "
            "with memory resources.]\n" );

         exit( EXIT_FAILURE );
      };
      ii++;
   };
/*............................................................................*/
# if GPH_USCALE == 1

   fprintf( stdout, "\n\t\t\t\t\t\t\t\t       ) <- ?" );
   fprintf( stdout, "\r Rescale coordinates ? >----------"
      "-------------------> [ y/n ] >---> (" );
   scanf( "%s", pptr );

   if ( *pptr == 'y' )
   {
     scaling:
      fprintf( stdout, " Please enter scaling factor [ ZERO "
         "< factor ] ....................: " );
      scanf( "%s", pptr );
      rscale = strtod( pptr, endp );

      if ( rscale < 1.e-277 )
         goto scaling;

      fprintf( stdout, " Enter x-shift [ -1 <= shift <= 1 ] "
         "...............................: ");
      scanf( "%s", pptr );
      x_shift = strtod( pptr , endp );

      fprintf( stdout, " Enter y-shift [ -1 <= shift <= 1 ] "
         "...............................: ");
      scanf( "%s", pptr );
      y_shift = strtod( pptr , endp );
   }; 
# endif

   ( gpt->xmin ) = 1.e+277;
   ( gpt->ymin ) = 1.e+277;

   ( gpt->xmax ) = -( gpt->xmin );
   ( gpt->ymax ) = -( gpt->ymin );

   kk = null;
   while( kk < ( gpt->nn ))
   {
      xx = ( gpt->vct[kk][null] );

      if (( gpt->xmax ) < xx )
          ( gpt->xmax ) = xx;

      if( xx < ( gpt->xmin ))
         ( gpt->xmin ) = xx;
      kk++;
   };

   ii = null;
   while( ii < ( gpt->nc ))
   {
      kk = null;
      while( kk < ( gpt->np[ii] ))
      {
         yy = ( gpt->vct[kk][ii+ONE] );

         if(( gpt->ymax ) < yy )
            ( gpt->ymax ) = yy;

         if( yy < ( gpt->ymin )) 
            ( gpt->ymin ) = yy;
         kk++;
      };
      ii++;
   };

   xx = (( gpt->xmax ) - ( gpt->xmin ))/2.;
   yy = (( gpt->ymax ) - ( gpt->ymin ))/2.;

# if GPH_USCALE == 1
   if( xx < yy ) 
      xx = yy;

   yy = xx;
# endif

   x_mean = (( gpt->xmin ) + ( gpt->xmax ))/2.;
   x_min = x_mean + ( x_shift - 1./rscale )*xx;
   x_max = x_mean + ( x_shift + 1./rscale )*xx;

   y_mean = (( gpt->ymin ) + ( gpt->ymax ))/2.;
   y_min = y_mean + ( y_shift - 1./rscale )*yy;
   y_max = y_mean + ( y_shift + 1./rscale )*yy;
/*............................................................................*/
/* plot system file */

   if ( null == strncmp(( gpt->format ), "gnuplot_2D", THREE ))
   {
/* renormalize absolute x and y values to match interval 
                                       [ 1.e-GPH_UNSCALED, 1.e+GPH_UNSCALED ] */
      xx = fabs( gpt->xmax );

      if ( xx < fabs( gpt->xmin ))
         xx = fabs( gpt->xmin );

      if ( ZERO < xx ) 
         lx = floor( log10( xx ));
      else
         lx = null;

      if ( abs( lx ) < GPH_UNSCALED )
      {
         xx = 1.;
         lx = null;
      }
      else
         xx = pow( 10., -lx ); /* = 10^(-lx ) = exp( -lx*log( 10.)) */

      x_min *= xx;
      x_max *= xx;

      yy = fabs( gpt->ymax );

      if ( yy < fabs( gpt->ymin ))
         yy = fabs( gpt->ymin );

      if ( ZERO < yy )
         ly = floor( log10( yy ));
      else
         ly = null;

      if ( abs( ly ) < GPH_UNSCALED )
      {
         yy = 1.;
         ly = null;
      }
      else
         yy = pow( 10., -ly );     /* = 10^(-ly ) = exp( -ly*log( 10.)) */

      y_min *= yy;
      y_max *= yy;

      if ( y_min == y_max )
      {
          y_min *= 0.9999;
          y_max *= 1.0001;
      };

      strcpy( prefix , "gnu." );

      GPH_FOPEN( prefix, ( gpt->file ));

      if ( null != ( gpt->dsp ))
         fprintf( stdout, "\n opened: gnuplot system file %s", pptr );

      fprintf( pltfle, "set title '%s %s'\n", gpt->name, gpt->text );
      fprintf( pltfle, "set xrange [%.15e:%.15e]\n", x_min, x_max );
      fprintf( pltfle, "set yrange [%.15e:%.15e]\n", y_min, y_max );

      if ( lx == null )
         fprintf( pltfle, "set xlabel '%s'\n", ( gpt->xunit ));
      else
         fprintf( pltfle, "set xlabel '1.E%+d %s'\n", lx, ( gpt->xunit ));

      if ( ly == null )
         fprintf( pltfle, "set ylabel '%s'\n", ( gpt->yunit ));
      else
         fprintf( pltfle, "set ylabel '1.E%+d %s'\n", ly, ( gpt->yunit ));

      fprintf( pltfle, "set grid\n" );
      fprintf( pltfle, "set border\n" );
      fprintf( pltfle, "%s %c\n", "plot", 92 );

# ifdef USER_PATH
      fprintf( pltfle, "%s/", USER_PATH );
# endif

      ii = null;
      while ( ii < (( gpt->nc )-ONE ))
      {
         if ( null < strlen( gpt->curve[ii] ))
            strcpy( dptr, ( gpt->curve[ii] ));
         else
         {
            strcpy( dptr, prefix );
            strcat( dptr, ( gpt->file ));
            strcat( dptr, GPH_DATFLE );
            strcat( dptr, lotos( ii, null ));
         };

         fprintf( pltfle, "%c%s%c with lines, %c\n", 39, dptr, 39, 92 );
	 ii++;            /* 39 = ASCII semicolon */ 
      };                  /* 92 = ASCII backslash */

      if ( ii < ( gpt->nc ))
      {
         if ( null < strlen( gpt->curve[ii] ))
            strcpy( dptr, ( gpt->curve[ii] ));
         else
         {
            strcpy( dptr, prefix );
            strcat( dptr, ( gpt->file ));
            strcat( dptr, GPH_DATFLE );
            strcat( dptr, lotos( ii, null ));
         };

         fprintf( pltfle, "%c%s%c with lines\n", 39, dptr, 39 );
      };

      fprintf( pltfle, "pause -1 '[ hit return to continue ]'\n" );

      nseconds = time( timer );
      timeptr = ctime( &nseconds );

      fprintf( pltfle, "\n# 2D-gnuplot system file '%s'\n", pptr );
      fprintf( pltfle, "# created:%.24s", timeptr );
      fprintf( pltfle, "\n# EOF" );

      fclose( pltfle );

      if ( null != ( gpt->dsp ))
      {
         fprintf( stdout, CLEAR_LINE );
         fprintf( stdout, "\r gnuplot system file %s", pptr );
         fprintf( stdout, timefrm, timeptr );
      };
/*............................................................................*/
/* 2D-plot points file: */ 

      ii = null;
      while ( ii < ( gpt->nc ))
      {
         if ( null < strlen( gpt->curve[ii] ))
            strcpy( dptr, ( gpt->curve[ii] ));
         else
         {
            strcpy( dptr, prefix );
            strcat( dptr, ( gpt->file ));
            strcat( dptr, GPH_DATFLE );
            strcat( dptr, lotos( ii, null ));
         };

         datfle = fopen( dptr, "w" );

         if (( datfle == null )
           &&( null != ( gpt->dsp )))
         {
            fprintf( stderr,
               "\n\n Message from function %s:", DO_GRAPHP );
            fprintf( stderr,
               "\n\n Unknown error on opening plot data file %s !!!", dptr );
            fprintf( stderr,
               "\n [ overrides and returns to calling program ].\n" );
            return null;
         };

         if ( null != ( gpt->dsp ))
            fprintf( stdout,
               "\n opened: gnuplot data file %s ", dptr );

         fprintf( datfle, "# %s", ( gpt->name ));

         if ( lx == null )
            fprintf( datfle, " [x-unit: %s, ", ( gpt->xunit ));
         else
            fprintf( datfle, " [x-unit: 1.E%+d %s, ", lx, ( gpt->xunit ));

         if ( ly == null )
            fprintf( datfle, "y-unit: %s]\n", ( gpt->yunit ));
         else
            fprintf( datfle, "y-unit: 1.E%+d %s]\n", ly, ( gpt->yunit ));

         kk=null;
	 while ( kk < ( gpt->np[ii] ))
         {
            uu = xx*( gpt->vct[kk][null] );
            vv = yy*( gpt->vct[kk][ii+ONE] );

# ifdef GPH_RNDOFF
            vv = rndoff( vv, ( GPH_RNDOFF+ONE ));
            fprintf( datfle, "%+.*e %+.*e\n",
               ( int ) GPH_RNDOFF, uu, ( int ) GPH_RNDOFF, vv );
# else
            fprintf( datfle, "%+.15e %+.15e\n", uu, vv );
# endif
            kk++;
         };

         fclose( datfle );

         if ( null != ( gpt->dsp ))
         {
            fprintf( stdout, CLEAR_LINE );
            fprintf( stdout, "\r gnuplot data file %s", dptr );
            fprintf( stdout, timefrm, timeptr );
         };
     
         ( gpt->np[(ii++)] ) = null; /* clear */
      };
   } /* end if ( gpt->format ) == "gnuplot_2D" */
   else if ( null == strncmp(( gpt->format ), "SPLINE", THREE ))
   {
      strcpy( prefix , "SPL." );

      GPH_FOPEN( prefix, ( gpt->file ));

      if ( null != ( gpt->dsp ))
         fprintf( stdout, "\n opened: SPLINE file '%s' ", pptr );

      fprintf( pltfle, "%s\n" , ( gpt->name ));
      fprintf( pltfle, "%s\n" , ( gpt->text ));
      fprintf( pltfle, "%s\n" , ( gpt->xunit ));
      fprintf( pltfle, "%s\n" , ( gpt->yunit ));
      fprintf( pltfle, "%ld\n", ( gpt->nn ));

      ii = null;
      while( ii< ( gpt->nc ))
      {
         kk = null;
         while( kk< ( gpt->np[ii] ))
         {
            uu = ( gpt->vct[kk][null] );
            vv = ( gpt->vct[kk][ii+ONE] );

# ifdef GPH_RNDOFF
            vv = rndoff( vv, ( GPH_RNDOFF+ONE ));
            fprintf( pltfle, "%+.*E   %+.*E\n",
               ( int ) GPH_RNDOFF, uu, ( int ) GPH_RNDOFF, vv );
# else
            fprintf( pltfle, "%+.15E   %+.15E\n", uu, vv );
# endif
            kk++;
         };
         fprintf( pltfle, "\n" );

         ( gpt->np[(ii++)] ) = null; /* clear */
      };

      nseconds = time( timer );
      timeptr = ctime( &nseconds );

      fprintf( pltfle, "\n# SPLINE file '%s'\n", pptr );
      fprintf( pltfle, "# created:%.24s", timeptr );
      fprintf( pltfle, "\n# EOF" );

      fclose( pltfle );

      if ( null != ( gpt->dsp ))
      {
         fprintf( stdout, CLEAR_LINE );
         fprintf( stdout, "\r SPLINE file '%s'", pptr );
         fprintf( stdout, timefrm, timeptr );
      };
   }; /* end if ( gpt->format ) == "SPLINE" */

   return ONE;
}
/*===========================================================================*/
# undef GPH_USCALE
# undef GPH_UNSCALED
# undef GPH_FOPEN
# undef GPH_DATFLE
/*************** end of plot file generation function graphp(*) ***************/
