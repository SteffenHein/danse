/* [ file: smithchrt.c ] */
# define DO_SMITHCHRT "smithchrt(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function smithchrt(*)                                               *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Creates smith chart for S-parameters computed in function modval(*).       *
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
# include <time.h>           /* cf. time( ),ctime( ),asctime( ),localtime( )  */
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* some user defined constants */
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
# include "../poster/POSTER.CONF" /* The general configuration header */
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
   # include <ncurses.h>
# endif 
/*----------------------------------------------------------------------------*/
# ifndef DO_ARGC
   # include "../math/argc.h"
# endif
# include "../tools/txctyp.h"
# include "../poster/posttp.h"
/*----------------------------------------------------------------------------*/
# define PST_SMITH 1             /* 1: plot in smith chart */
# define PST_POLAR 1             /* 1: plot in polar coardinates */
# define PST_CRTSN 1             /* 1: plot in Cartesian coordinates */
# define PST_SCALE 1             /* 0: crt. plot range [ -1., 1] X [-1., 1.] */
                                 /* 1: "... [x_min, x_max] X [y_min, y_max] */
                                 /* 2: "... [-|z_max|, |z_max| ]^2 */
# ifdef EVL_WRTFREQ
   # define PST_SVFRQ EVL_WRTFREQ /* 1: read/write frequencies */
# else                            /*    from/into 'spc.xyz' files */ 
   # define PST_SVFRQ 0
# endif

# define PST_DIMNS 500           /* maximum number of frequency points        */
# define PST_PRECISION (1.e-12 ) /* relative precision                        */
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
/*============================================================================*/

POSTSTATE
*smithchrt( POSTSTATE *stp )
{
/* declarations: */

   static POSTSTATE
     *state;

   static FILE 
     *stream;

   static COMPLEX 
      cc = {null},
     *cpt = &cc;

   static TXCNSL
     *csp;

   static time_t 
      nseconds = null,
     *timer = null;

   static short
      ind = null;

   static long
      ii = null,
      jj = null,
      nn = null,
      mm = null;

   static double
      xx1 = ZERO,
      xx2 = ZERO,
      xx3 = ZERO,
      dx = ZERO,
      uu1 = ZERO,
      uu2 = ZERO,
      uu3 = ZERO,
      vv1 = ZERO,
      vv2 = ZERO,
      vv3 = ZERO,
      phi1 = ZERO,
      phi2 = ZERO,
      phi3 = ZERO;

   static const char
      items = NINE;

   static char
      option = null,
      mflag[THREE] = {null},
      fleptr[STS_SIZE] = {null},
      tlmmod[STS_SIZE] = {null},
      xunit[STS_SIZE] = {null},
      yunit[STS_SIZE] = {null},
      ptr[STS_SIZE] = {null},
      plot[STS_SIZE] = {null},
      dat1[STS_SIZE] = {null},
      dat2[STS_SIZE] = {null},
     *tmeptr,
    **endp = null;

   static double
      ssn1[PST_DIMNS][TWO] = {{ZERO}};

/* prototypes: */

   time_t time( time_t *timer );
   char *ctime( const time_t *timer );

   COMPLEX
     *argc( COMPLEX *cpt, short nn );

   char *
      lotos( long n, char m );

   SPLINES
     *spline( SPLINES *spt );

   int
     graphp( GRAPHICS *gpt );

   TXCNSL 
     *txcnsl( TXCNSL *csp );
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

/*............................................................................*/
   csp = txcnsl( null );   /* initialize text console                         */
/*.......................*/

   strcpy( csp->cnfrm, "Nothing done, do you really want to quit ?" );
   ( csp->dfopt ) = 1; /* initial default menu option */
   mflag[0] = ONE;
   mflag[1] = -ONE;
   mflag[2] = null;

  menu:

   ( csp->mflag[8] ) = mflag[0];
   ( csp->mflag[9] ) = mflag[1];
   ( csp->mflag[10] ) = mflag[2];

   strcpy( csp->cmmnt, "Plot file builder" );

# if PST_SMITH == 1
   strcpy( csp->envmt, "SMITH" );
# else
   strcpy( csp->envmt, "PLOT" );
# endif

   strcpy( csp->tasks, "Create plot files of ..." );
   strcpy( csp->flags, "using ..." );
   strcpy( csp->escpe, "End of program / escape:" );

   ( csp->items ) = items;

   strcpy( csp->mline[1], "* reflection coefficient >---> [ s11 ]" );
   strcpy( csp->mline[2], "* transmission coefficient >-> [ s21 ]" );
   strcpy( csp->mline[3], "* coupling coefficient >-----> [ s31 ]" );
   strcpy( csp->mline[4], "* coupling coefficient >-----> [ s41 ]" );
   strcpy( csp->mline[5], "* directivity >----------> [ s41/s21 ]" );
   strcpy( csp->mline[6], "* any S-parameters specified by file name" );
   strcpy( csp->mline[7], "* Plot all S-parameters" );
   strcpy( csp->mline[8], "* cartesian coordinates [ Gauss plane ]" );
   strcpy( csp->mline[9], "* Smith chart" );
   strcpy( csp->mline[10],"* polar coordinates" );
/*............................................................................*/
   csp = txcnsl( csp );    /* build the menu [ on text console ]              */
/*.......................*/
   option = ( csp->option );

   mflag[0] = ( csp->mflin[8] );
   mflag[1] = ( csp->mflin[9] );
   mflag[2] = ( csp->mflin[10] );
   
   if ( option == null )
   {
      PRBLDCLR( "\r" );
# if PST_SMITH == 1
      printf( "\t\t\t\t\t\t\t\t\t  SMITH" );
# else
      printf( "\t\t\t\t\t\t\t\t\t   PLOT" );
# endif
      PRNORMAL( "" );

      ( state->rtn ) = null;
      return state;
   }
   else if ( option == 6 )
   {
      strcpy( csp->rqstr, "Please enter file name" );
      strcpy( csp->rqfrm, "points" );
/*............................................................................*/
      csp = txcnsl( csp );    /*                                              */
/*..........................*/
      strncpy ( fleptr, csp->instr, SHS_SIZE );
      jj = null;
   }
   else if ( option == 7 )
   {
      printf( "\n Option %d is presently not implemented.", option );
      ( csp->dfopt ) = null; /* next default menu option */
      goto menu;
   }
   else if(( 0 < option )
         &&( option < 6 ))
   {
      jj = option;

      strcpy( ptr, lotos(( long) jj, null ));

      strcpy( csp->rqlng,
         "Please enter job index N of s-parameter file spc.s" );
      strcat( csp->rqlng, ptr );
      strcat( csp->rqlng, "1_<N>" );
/*
      strcpy( csp->rqfrm, "points" );
*/
/*............................................................................*/
      csp = txcnsl( csp );    /*                                              */
/*..........................*/
      ind = ( csp->inlng );
      printf( "\n" );

      strcpy ( fleptr, "spc.s" );
      strcat ( fleptr, ptr );
      strcat ( fleptr, "1_" );
      strcat( fleptr, lotos( ind, null ));
   }
   else
   {
      PRBLDCLR( "\r" );
# if PST_SMITH == 1
      printf( "\t\t\t\t\t\t\t\t\t  SMITH" );
# else
      printf( "\t\t\t\t\t\t\t\t\t   PLOT" );
# endif
      PRNORMAL( "" );

      ( state->rtn ) = null;
      return state;
   };
/*............................................................................*/
/* enter s-parameters from file: */

   stream = fopen( fleptr, "r+" );

   if( null == stream )
   {
      printf( "\n Error message from function '%s':", DO_SMITHCHRT );
      printf( "\n Cannot open file '%s' !!! ", fleptr );
      printf( "\n [ Please verify if file is actually"
              " present in chosen directory.] " );
      printf( "\n " );
      exit( EXIT_FAILURE );
   };

   fscanf( stream, "%s", tlmmod ); /* DSC system identifier */
   fscanf( stream, "%s", ptr );    /* comment */
   fscanf( stream, "%s", xunit );
   fscanf( stream, "%s", yunit );

   fscanf( stream, "%s", ptr );
   xx1 = strtod( ptr, endp );
   fscanf( stream, "%s", ptr );
   xx3 = strtod( ptr, endp );
   fscanf( stream, "%s", ptr );
   dx = strtod( ptr, endp );
   fscanf( stream, "%s", ptr );

   mm = strtol( ptr, endp, DEC ) - ONE;

   if( PST_DIMNS <= mm )
   {
      printf( "\r Error message from function '%s':", DO_SMITHCHRT );
      printf( "\n Too many frequency points !!! " );
      printf( "\n [ The maximum number is %ld = macro PST_DIMNS ",
                 (long) PST_DIMNS );
      printf( "\n   - change macro in '%s' ]", DO_SMITHCHRT );
      printf( "\n " );
      exit( EXIT_FAILURE );
   };

   nn = ( mm/2 );
   xx2 = xx1 + nn*dx;

   ii = null;

# if PST_SVFRQ == 1
   fscanf( stream, "%s", ptr );
# endif

   fscanf( stream, "%s", ptr );
   ssn1[ii][null] = strtod( ptr, endp );
   fscanf( stream, "%s", ptr );
   ssn1[ii][ONE] = strtod( ptr, endp );
   uu1 = ssn1[null][null];
   vv1 = ssn1[null][ONE];
   while( ii < nn )
   {
      ii++ ;

# if PST_SVFRQ == 1
   fscanf( stream, "%s", ptr );
# endif

      fscanf( stream, "%s", ptr );
      ssn1[ii][null] = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      ssn1[ii][ONE] = strtod( ptr, endp );
   };
   uu2 = ssn1[nn][null];
   vv2 = ssn1[nn][ONE];
   while( ii < mm )
   {
      ii++ ;

# if PST_SVFRQ == 1
   fscanf( stream, "%s", ptr );
# endif

      fscanf( stream, "%s", ptr );
      ssn1[ii][null] = strtod( ptr, endp );
      fscanf( stream, "%s", ptr );
      ssn1[ii][ONE] = strtod( ptr, endp );
   };
   uu3 = ssn1[mm][null];
   vv3 = ssn1[mm][ONE];

   fclose( stream );

   ( cpt->r ) = uu1;
   ( cpt->i ) = vv1;
/*............................................................................*/
   cpt = argc( cpt, ONE );    /* computes = arg(z) [ phase range phi1 ]    */
/*.............................*/
   phi1 = 180.*( cpt->arg )/PI;

   ( cpt->r ) = uu2;
   ( cpt->i ) = vv2;
/*............................................................................*/
   cpt = argc( cpt, ONE );    /* computes = arg(z) [ phase range phi2 ]    */
/*.............................*/
   phi2 = 180.*( cpt->arg )/PI;

   ( cpt->r ) = uu3;
   ( cpt->i ) = vv3;
/*............................................................................*/
   cpt = argc( cpt, ONE );    /* computes = arg(z) [ phase range phi3 ]    */
/*.............................*/
   phi3 = 180.*( cpt->arg )/PI;
      
   if( null < jj )
   {
      strcpy( ptr, "s" );
      strcat( ptr, lotos(( long ) jj, null ));
      strcat( ptr, "1_" );
      strcat( ptr, lotos( ind, null ));
   }
   else if ( null == jj )
   {
      strcpy( ptr, fleptr );
   };

   strcpy( dat1, "gnu." );
   strcat( dat1, ptr );
   strcpy( dat2, dat1 );
   strcat( dat1, ".lower" );
   strcat( dat2, ".upper" );

/*............................................................................*/
# if PST_SMITH == 1
/*............................................................................*/
/* gnuplot smith chart: */

   if( mflag[1] == ONE )
   {
      strcpy( plot, "smc." );
      strcat( plot, ptr );

      stream = fopen( plot, "w" );
/*............................................................................*/
/*
#!/user7/vonhagen/bin/gnuplot -persist
# adjust previous line to local conditions
#
# you need a file with the following ordering
# freq Re(S) Im(S)
# example: (put this into a file called "vdns11.dat"
#----------------snip-------------------
## remove the 1st # in each line
## van der Neut's patch antenna computed by HP ADS
## It's in a technical report, basically
## a off-centre fed
## don't get ffoled, the following data is normalized
## to 100 \ohm, not 50 \ohm
## here comes the data:
##
## f/GHz   Re(S_{11})  Im(S_{11})
#4.20      -0.181	 0.525
#4.22      -0.121	 0.478
#4.24      -0.067	 0.418
#4.26      -0.021	 0.347
#4.28       0.012	 0.265
#4.30       0.030	 0.176
#4.32       0.031	 0.084
#4.34       0.014	-0.006
#4.36      -0.018	-0.089
#4.38      -0.065	-0.163
#4.40      -0.120	-0.244
##-------------snip
# 
# in gnuplot
# gnuplot> load "smith.gpl" # plots smith chart
# gnuplot> rep "vdns11.dat" u 2:3
#
# released under GPPL (c) Juergen v.Hagen 19990506
*/
      fprintf( stream, "set title '%s'", ptr );

      fprintf( stream, "\nset angle degree" );
      fprintf( stream, "\nset grid polar 360.000000" );
      fprintf( stream, "\nset grid xtics noytics noztics nox2tics "
                       "noy2tics %c", 92 );
      fprintf( stream, "\nnomxtics nomytics nomztics nomx2tics "
                       "nomy2tics %c", 92 );
      fprintf( stream, "\nlt 0 lw 1.000, lt 0 lw 1.000" );
      fprintf( stream, "\nset samples 100, 100" );
      fprintf( stream, "\nset size ratio 1 1,1" );
      fprintf( stream, "\nset origin 0,0" );
      fprintf( stream, "\nset data style points" );
      fprintf( stream, "\nset function style linespoints" );
      fprintf( stream, "\nset tics in" );
      fprintf( stream, "\nset ticslevel 0.5" );
      fprintf( stream, "\nset noxtics #axis nomirror norotate autofreq" ); 
      fprintf( stream, "\nset noytics #axis nomirror norotate autofreq" ); 
      fprintf( stream, "\nset rrange [ 0 : 1 ] noreverse nowriteback" );
      fprintf( stream, "\nset trange [ 0 : 360 ] noreverse nowriteback" );
      fprintf( stream, "\nset xrange [ -1.05 : 1.05 ] noreverse "
                       "nowriteback" );
      fprintf( stream, "\nset yrange [ -1.05 : 1.05 ] noreverse "
                       "nowriteback" );
      fprintf( stream, "\nset zero 1e-08" );

      fprintf( stream, "\nset nopolar" );
      fprintf( stream, "\nset para" );
/*
# stop if circle leaves the |S| = 1 condition.
# we only consider here passive S parameters
*/
      fprintf( stream, "\nxblend(t,x,max) = %c", 92 );
      fprintf( stream, "\n((1./x*cos(t)+1.)**2 + (sin(t)+1.)**2/x**2 %c",
                            92 );
      fprintf( stream, "\n>= max**2 ?  1./0. : 1.)" );
/*
# impedance functions
# crx = constant resistance, x coordinate
# cry = constant resistance, y coordinate
# cxx = constant reactance, x coordinate
# cxy = constant reactance, y coordinate
*/
      fprintf( stream, "\ncrx(t,r) = (cos(t)+r)/(r+1.)" );
      fprintf( stream, "\ncry(t,r) = 1.*sin(t)/(r+1.)" );
      fprintf( stream, "\ncxx(t,x) = (1./x*cos(t) + 1.)*xblend(t,x,1.01)" );
      fprintf( stream, "\ncxy(t,x) = (sin(t) + 1.)*xblend(t,x,1.01)/x" );
/*
# admittance functions
# cgx = constant G, x coordinate
# cgy = constant G, y coordinate
# cbx = constant B, x coordinate
# cby = constant B, y coordinate
*/
      fprintf( stream, "\ncgx(t,r) = (cos(t)-r)/(r+1.)" );
      fprintf( stream, "\ncgy(t,r) = 1.*sin(t)/(r+1.)" );
      fprintf( stream, "\ncbx(t,x) = "
                       "- (1./x*cos(t) + 1.)*xblend(t,x,1.01)" );
      fprintf( stream, "\ncby(t,x) = (sin(t) + 1.)*xblend(t,x,1.01)/x" );
      fprintf( stream, "\nset samples 300" );

      fprintf( stream, "\npl t-180, 0 noti w lines 3, %c", 92 );
      fprintf( stream, "\ncxx(t,1), cxy(t,1) noti w lines 3, %c", 92 );
      fprintf( stream, "\ncxx(t,-1), cxy(t,-1) noti w lines 3, %c", 92 );
      fprintf( stream, "\ncxx(t,0.5), cxy(t,0.5) noti w lines 3, %c", 92 );
      fprintf( stream, "\ncxx(t,-0.5), cxy(t,-0.5) noti w lines 3, %c", 92 );
      fprintf( stream, "\ncrx(t,0), cry(t,0) noti w lines 3, %c", 92 );
      fprintf( stream, "\ncrx(t,1), cry(t,1) noti w lines 3, %c", 92 );
      fprintf( stream, "\ncrx(t,0.3333), cry(t,0.3333) noti w lines 3, %c",
                        92 );
      fprintf( stream, "\ncrx(t,3), cry(t,3) noti w lines 3" );
      fprintf( stream, "\nset key Left reverse 1.25,0.0" );
/*
#
# the above command is for impedance form.
# for admittance form use:
#
      fprintf( stream, "\npl t-180, 0 noti w lines 1, %c", 92 );
      fprintf( stream, "\ncbx(t,1), cby(t,1) noti w lines 1, %c", 92 );
      fprintf( stream, "\ncbx(t,-1), cby(t,-1) noti w lines 1, %c", 92 );
      fprintf( stream, "\ncbx(t,0.5), cby(t,0.5) noti w lines 1, %c", 92 );
      fprintf( stream, "\ncbx(t,-0.5), cby(t,-0.5) noti w lines 1, %c", 92 );
      fprintf( stream, "\ncgx(t,0), cgy(t,0) noti w lines 1, %c", 92 );
      fprintf( stream, "\ncgx(t,1), cgy(t,1) noti	w lines 1, %c", 92 );
      fprintf( stream, "\ncgx(t,0.3333), cgy(t,0.3333) noti w lines 1, %c",
                        92 );
      fprintf( stream, "\ncgx(t,3), cgy(t,3) noti w lines 1" );
#
#  EOF
*/
      fprintf( stream, "\nset border" );
      fprintf( stream, "\nreplot %c", 92 );

      fprintf( stream, "\n'%s',%c", dat1, 92 );
      fprintf( stream, "\n'%s'", dat2 );
/*
      fprintf( stream, "\n'%s' with lines ,%c", dat1, 92 );
      fprintf( stream, "\n'%s' with lines", dat2 );
*/
      fprintf( stream, "\npause -1 '[ hit return to continue ]'" );

      nseconds = time( timer );
      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n#%s%s%s\n#%s", "ASCII file '",
                 plot, "' created:", tmeptr );
      fprintf( stream, "#EOF" );
      fclose( stream );

      printf( "\r %s%s%s\n %s\r", "ASCII file '",
                  plot, "' created:", tmeptr );
   }; /* end if ( mflag[1] == ONE ) */
# endif
/*............................................................................*/
# if PST_CRTSN == 1
/*............................................................................*/
/* gnuplot, Cartesian coordinates: */

   if( mflag[0] == ONE )
   {

# if PST_SCALE == 0
      uu1 = -1.03;
      uu2 = 1.03;
      vv1 = -1.03;
      vv2 = 1.03;
# else
      uu1 = 1.0e+77;
      uu2 = -uu1;
      vv1 = uu1;
      vv2 = -vv1;

      ii = null;
      while( ii <= mm )
      {
	 uu3 = ssn1[ii][null];
	 vv3 = ssn1[ii][ONE];

# if PST_SCALE == 1
         if ( uu3 < uu1 )
            uu1 = uu3;
         if ( uu2 < uu3 )
            uu2 = uu3;
         if ( vv3 < vv1 )
            vv1 = vv3;
         if ( vv2 < vv3 )
            vv2 = vv3;
# elif PST_SCALE == 2
         uu1 = uu3*uu3 + vv3*vv3;

	 if ( uu2 < uu1 )
	    uu2 = uu1;
# endif
         ii++ ;
      };

# if PST_SCALE == 1
      uu3 = .7*( uu2 - uu1 );
      vv3 = .7*( vv2 - vv1 );
      if ( vv3 < uu3 )
      {
         uu1 = .5*( uu1 + uu2 ) - uu3;
         uu2 = uu1 + 2.*uu3;
         vv1 = .5*( vv1 + vv2 ) - uu3;
         vv2 = vv1 + 2.*uu3;
      }
      else
      {
         uu1 = .5*( uu1 + uu2 ) - vv3;
         uu2 = uu1 + 2.*vv3;
         vv1 = .5*( vv1 + vv2 ) - vv3;
         vv2 = vv1 + 2.*vv3;
      };
      
# elif PST_SCALE == 2
      uu2 = sqrt( uu2 );
      uu1 = - uu2;
      vv1 = uu1;
      vv2 = uu2;
# endif

      if ( uu1 < -1.03 )
         uu1 = -1.03;
      if ( 1.03 < uu2 )
         uu2 = 1.03;
      if ( vv1 < -1.03 )
         vv1 = -1.03;
      if ( 1.03 < vv2 )
         vv2 = 1.03;
# endif

      strcpy( plot, "crt." );
      strcat( plot, ptr );

/*............................................................................*/

      stream = fopen( plot, "w" );

      fprintf( stream, "set title '%s'", ptr );
      fprintf( stream, "\nset xrange [%.7e:%.7e]", uu1, uu2 );
      fprintf( stream, "\nset yrange [%.7e:%.7e]", vv1, vv2 );
      fprintf( stream, "\nset xlabel 'Re(%s)'", ptr );
      fprintf( stream, "\nset ylabel 'Im(%s)'", ptr );
/*
      fprintf( stream, "\nset size square" );
*/
      fprintf( stream, "\nset size ratio 1 1,1" );
      fprintf( stream, "\nset grid" );
      fprintf( stream, "\nset border" );
      fprintf( stream, "\nplot %c", 92 );
/*
      fprintf( stream, "\n'%s',%c", dat1, 92 );
      fprintf( stream, "\n'%s'", dat2 );
*/
      fprintf( stream, "\n'%s' with lines ,%c", dat1, 92 );
      fprintf( stream, "\n'%s' with lines", dat2 );

      fprintf( stream, "\npause -1 '[ hit return to continue ]'" );

      nseconds = time( timer );
      tmeptr = ctime( &nseconds );
      fprintf( stream, "\n\n#%s%s%s\n#%s", "ASCII file '",
                         plot, "' created:", tmeptr );
      fprintf( stream, "#EOF" );
      fclose( stream );

      printf( "\r %s%s%s\n %s\r", "ASCII file '",
                  plot, "' created:", tmeptr );
   };/* end if ( mflag[0] == ONE ) */
# endif
/*............................................................................*/
/* create data files: */
/*............................................................................*/
/* lower frequency domain: */

   stream = fopen( dat1, "w" );

   fprintf( stream, "#%s", tlmmod );
   fprintf( stream, "\n#%s_[lower_frequency_domain]", ptr );
   fprintf( stream, "\n#%s", xunit );
   fprintf( stream, "\n#%s", yunit );
   fprintf( stream, "\n#%+.15e", xx1 );
   fprintf( stream, "\n#%+.15e", xx2 );
   fprintf( stream, "\n#%+.15e", dx );
   fprintf( stream, "\n#%ld", ( nn + ONE ));

   ii = null;
   while( ii <= nn )
   {
      fprintf( stream, "\n%+.15e", ssn1[ii][null] );
      fprintf( stream, "   %+.15e", ssn1[ii][ONE] );
      ii++ ;
   };

   fprintf( stream, "\n\n#lower_phase_shift: %+.15e deg", phi1-phi2 );
   fprintf( stream, "\n#total_phase_shift: %+.15e deg", phi1-phi3 );

   nseconds = time( timer );
   tmeptr = ctime( &nseconds );
   fprintf( stream, "\n\n#%s%s%s\n#%s", "ASCII file '",
                          dat1, "' created:", tmeptr );
   fprintf( stream, "#EOF" );
   fclose( stream );

   printf( "\r %s%s%s\n %s\r", "ASCII file '",
               dat1, "' created:", tmeptr );
/*............................................................................*/
/* upper frequency domain: */

   stream = fopen( dat2, "w" );

   fprintf( stream, "#%s", tlmmod );
   fprintf( stream, "\n#%s_[upper_frequency_domain]", ptr );
   fprintf( stream, "\n#%s", xunit );
   fprintf( stream, "\n#%s", yunit );
   fprintf( stream, "\n#%+.15e", xx2 );
   fprintf( stream, "\n#%+.15e", xx3 );
   fprintf( stream, "\n#%+.15e", dx );
   fprintf( stream, "\n#%ld", ( mm - nn + ONE ));

   ii = nn;
   while( ii <= mm )
   {
      fprintf( stream, "\n%+.15e", ssn1[ii][null] );
      fprintf( stream, "   %+.15e", ssn1[ii][ONE] );
      ii++ ;
   };

   fprintf( stream, "\n\n#upper_phase_shift: %+.15e deg", phi2-phi3 );
   fprintf( stream, "\n#total_phase_shift: %+.15e deg", phi1-phi3 );

   nseconds = time( timer );
   tmeptr = ctime( &nseconds );
   fprintf( stream, "\n\n#%s%s%s\n#%s", "ASCII file '",
                             dat2, "' created:", tmeptr );
   fprintf( stream, "#EOF" );
   fclose( stream );

   printf( "\r %s%s%s\n %.24s\r", "ASCII file '",
               dat2, "' created:", tmeptr );
/*............................................................................*/

   ( csp->dfopt ) = null; /* next default menu option */
   goto menu;

/*............................................................................*/
   return state;
}
/*============================================================================*/
# undef PST_SMITH
# undef PST_POLAR
# undef PST_CRTSN
# undef PST_SCALE
# undef PST_SVFRQ
/*********************** end of function smithchrt(*) *************************/
