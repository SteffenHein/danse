/* [ file: pstprc2.h ] */
# define DO_SMITHCHRT "smithchrt(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function smithchrt.c                                                *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Creates smith chart for S-parameters computed with 'modval(*)'             *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*----------------------------------------------------------------------------*/
# define PST_SMITH 1             /* 1: plot in smith chart */
# define PST_POLAR 1             /* 1: plot in polar coardinates */
# define PST_CRTSN 1             /* 1: plot in Cartesian coordinates */
# define PST_SCALE 1             /* 0: crt. plot range [ -1., 1] X [-1., 1.] */
                                 /* 1: "... [x_min, x_max] X [y_min, y_max] */
                                 /* 2: "... [-|z_max|, |z_max| ]^2 */
# ifdef EVL_WRTFREQ
   # define PST_SVFRQ EVL_WRTFREQ /* 1: save frequencies in files cpl.sN1 */
# else
   # define PST_SVFRQ 0
# endif
/*............................................................................*/
# define PST_DIMNS 2000          /* maximum number of frequency points        */
# define PST_PRECISION ( 1.000e-12 ) /* relative precision                    */
/*----------------------------------------------------------------------------*/
# include "../math/argc.h"
/*
# include "../math/txctyp.h"
*/
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

POSTSTATE *\
smithchrt( POSTSTATE *stp )
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
      items = SIX;

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

   ( csp->clscr ) = null; /* clear screen, then that number of line feeds */
   ( csp->dfopt ) = ONE; /* initial default menu option */
   ( csp->dflnf ) = null; /* no line feed before default option line */

   mflag[0] = ONE;
   mflag[1] = -ONE;
   mflag[2] = null;

  menu:

   strcpy( csp->cmmnt, "Plot file builder" );

# if PST_SMITH == 1
   strcpy( csp->envmt, "SMITH" );
# else
   strcpy( csp->envmt, "PLOT" );
# endif

   ( csp->mflag[5] ) = mflag[0];
   ( csp->mflag[6] ) = mflag[1];
   ( csp->mflag[7] ) = mflag[2];

   strcpy( csp->tasks, "Create plot files of ..." );
   strcpy( csp->flags, "using ..." );
   strcpy( csp->escpe, "End of program / escape:" );

   ( csp->items ) = items;

   strcpy( csp->mline[1], "* reflection coefficient >---> [ s11 ]" );
   strcpy( csp->mline[2], "* transmission coefficient >-> [ s21 ]" );
/*
   strcpy( csp->mline[3], "* coupling coefficient >-----> [ s31 ]" );
   strcpy( csp->mline[4], "* coupling coefficient >-----> [ s41 ]" );
   strcpy( csp->mline[5], "* directivity >----------> [ s41/s21 ]" );
*/
   strcpy( csp->mline[3], "* any S-parameters specified by file name" );
   strcpy( csp->mline[4], "* Plot all S-parameters" );
   strcpy( csp->mline[5], "* cartesian coordinates [ Gauss plane ]" );
   strcpy( csp->mline[6], "* Smith chart" );
   strcpy( csp->mline[7], "* polar coordinates" );
/*............................................................................*/
   csp = txcnsl( csp );    /* build the menu [ on text console ]              */
/*.......................*/
   option = ( csp->option );

   mflag[0] = ( csp->mflin[5] );
   mflag[1] = ( csp->mflin[6] );
   mflag[2] = ( csp->mflin[7] );
   
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
   else if ( option == THREE )
   {
      strcpy( csp->rqstr,
               "Please enter file name" );
      strcpy( csp->rqfrm, "points" );
/*............................................................................*/
      csp = txcnsl( csp );    /*                                              */
/*..........................*/
      strncpy ( fleptr, csp->instr, SHS_SIZE );
      jj = null;
   }
   else if ( option == FOUR )
   {
      PRBLDCLR( "\r This option is presently not implemented." );
# if PST_SMITH == 1
      printf( "\t\t\t\t  SMITH" );
# else
      printf( "\t\t\t\t   PLOT" );
# endif
      PRNORMAL( "" );

      ( csp->dfopt ) = null; /* next default menu option */

      goto menu;
   }
   else if(( null < option )&&( option < 3 ))
   {
      jj = option;

      strcpy( ptr, lotos(( long) jj, null ));
      strcpy( csp->rqstr,
               "Please enter name of s" );
      strcat(( csp->rqstr ), ptr );
      strcat(( csp->rqstr ), "1-parameter file" );

      strcpy(( csp->dfstr ), "cpl.s" );
      strcat(( csp->dfstr ), ptr );
      strcat(( csp->dfstr ), "1" );

      strcpy( csp->rqfrm, "points" );

/*............................................................................*/
      csp = txcnsl( csp );    /*                                              */
/*..........................*/
      printf( "\n" );
      strcpy( fleptr, ( csp->instr ));
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
      printf( "\n Cannot open file '%s' !!!", fleptr );
      printf( "\n [ Please verify if file exists" \
         " in presently working directory.]" ); \
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
      printf( "\r Error message from function %s:", DO_SMITHCHRT );
      printf( "\n Too many frequency points !!!" );
      printf( "\n [ The maximum number is %ld = macro PST_DIMNS ",
         (long) PST_DIMNS );
      printf( "\n   - change macro in %s ]", DO_SMITHCHRT );
      printf( "\n " );
      exit( EXIT_FAILURE );
   };

   nn = ( mm/2 );
   xx2 = xx1 + nn*dx;

   ii = null;

# if PST_SVFRQ == 1
   fscanf( stream, "%s", ptr ); /* frequency */
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
      fscanf( stream, "%s", ptr ); /* frequency */
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
      fscanf( stream, "%s", ptr ); /* frequency */
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
      fprintf( stream, "\nset style data points" );
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
      fprintf( stream, "\ncgx(t,r) = ( cos(t)-r)/(r+1.)" );
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

   if ( option < TWO )
      ( csp->dfopt ) = option + ONE; /* next default menu option */
   else
      ( csp->dfopt ) = null;

   goto menu;

/*............................................................................*/
   return state;
}
/*============================================================================*/
# undef PST_SMITH
# undef PST_POLAR
# undef PST_CRTSN
/*********************** end of function smithchrt.c **************************/










# define DO_PSTPRC "pstprc(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function pstprc(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Computes the S-parameters of a symmetric 2 port by evaluating the          *
*   reflection coefficients of the two systems obtained with setting           *
*   magnetic and electric walls into the symmetry axis.                        *
*   [ Obviously, this is nothing but the well-known open/short circuit         *
*     network analysis for symmetric devices.]                                 *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: April 14, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*----------------------------------------------------------------------------*/
# ifndef PST_DIMNS
   # define PST_DIMNS 2000          /* maximum number of frequency points */
# endif
# ifndef PST_PRECISION
   # define PST_PRECISION (1.e-12 ) /* relative precision */
# endif
# ifndef PST_PHASES
   # define PST_PHASES  1       /* 1: compute phases */
# endif
# if PST_PHASES == 1
   # ifndef PST_PBRNCH
      # define PST_PBRNCH 0
   # endif
   # ifndef PST_PHUNIT
      # define PST_PHUNIT "DEG"  /* phase unit: "DEG" or "rad" */
   # endif                       /* ["DEG" any char / "rad" lower case char ] */
# endif
# ifndef PST_SVFRQ
   # define PST_SVFRQ 1          /* 1: save frequencies in cpl.s11,... files */
# endif
/*----------------------------------------------------------------------------*/
/*
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
*/
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
   if( null == stream ) \
   { \
      printf( "\n\n Error message from function pstprc(*):" ); \
      printf( "\n Cannot open file '%s' !!!", fleptr ); \
      printf( "\n [ Please verify if file is actually" \
         " present in chosen directory.]" ); \
      printf( "\n " ); \
      exit( EXIT_FAILURE ); \
   }; 
/*----------------------------------------------------------------------------*/
# define FLE_CHCK2( ) \
   ind = null; \
   if( null == strstr( ptr, eop )) \
   { \
      printf( "\n Sorry, there are no %s mode S-parameters" \
              " in this file !!!", eop ); \
      printf( "\n [ Please re-enter file name.]\n " ); \
      ind = ONE; \
   };
/*----------------------------------------------------------------------------*/
# define PST_SPLCPY( OPT, IPL ) \
   { \
      spt->vct[ii][null] = ss; \
      if ((OPT) == 'p' ) \
         spt->vct[ii][ONE] = 100.*sqrt( real1*real1 + imag1*imag1 ); \
      else if ((OPT) == 'd' ) \
         spt->vct[ii][ONE] = 10.*log10( real1*real1 + imag1*imag1 ); \
      else if ((OPT) == 'r' ) \
         spt->vct[ii][ONE] = real1; \
      ind = null; \
      do \
      { \
         spt->dmn[nn] = ss; \
         ss += ( dx/((IPL) + ONE )); \
         ind++ ; \
         nn++ ; \
      } while(( ind <= (IPL))&&( ii < ( mm - ONE ))); \
   } \
/*----------------------------------------------------------------------------*/
/*
# if PST_PHASES == 1
   # include "math/argc.h"
# endif
*/
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
      imag1 = ZERO,
      ss = ZERO;

   static const char
      items = FOUR,
     *evemde = "even",
     *oddmde = "odd";

# if PST_PHASES == 1
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
      scpl11[PST_DIMNS][TWO] = {{ZERO}},
      scpl21[PST_DIMNS][TWO] = {{ZERO}},
      sev11[PST_DIMNS][TWO] = {{ZERO}},
      sod11[PST_DIMNS][TWO] = {{ZERO}};

# if PST_PHASES == 1
   static COMPLEX cc = {null},
                *cpt = &cc;

   static double
      phs11[PST_DIMNS] = {ZERO},
      phs21[PST_DIMNS] = {ZERO};
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

   POSTSTATE
     *smithchrt( POSTSTATE *stp );

   char 
     *lotos( long n, char m );

# if PST_PHASES == 1

   COMPLEX
     *argc( COMPLEX *ipt, short nn );

# endif
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

   ( csp->clscr ) = null; /* clear screen, then that number of line feeds */
   ( csp->dfopt ) = THREE;
   ( csp->dflnf ) = null; /* no line feed before default option line */

   strcpy( csp->tasks, "COMPUTE ... " );
   strcpy( csp->cmmnt, "Even/odd mode postprocessing:" );
   strcpy( csp->envmt, "PSTPRC" );

  menu:

   ( csp->items ) = items;

   strcpy( csp->mline[1], "* reflection coefficient >--> [ s11 ]" );
   strcpy( csp->mline[2], "* coupling coefficient >----> [ s21 ]" );
   strcpy( csp->mline[3], "* compute all S-parameters" );
   strcpy( csp->mline[4], "* generate Smith chart" );

   strcpy( csp->escpe, "Return to main program / escape:" );

/*............................................................................*/
   csp = txcnsl( csp );    /*                                                 */
/*.......................*/
   
   option = ( csp->option );

   if(( null < option )&&( option < FOUR ))
      printf( "\n" );

/*............................................................................*/
/* s_even: */

   if(( null < option )&&( option <= THREE ))
   {
     even_mode:

      strcpy( eop, evemde ); 

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
         goto even_mode;

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

      if( PST_DIMNS < mm )
      {
         printf( "\n\n Error message from function %s:", DO_PSTPRC );
         printf( "\n Too many frequency points !!!" );
         printf( "\n [ The maximum number is %ld = macro PST_DIMNS",
            (long) PST_DIMNS );
         printf( "\n   - change macro in SCPLCONF.H ]" );
         printf( "\n\n " );
         exit( EXIT_FAILURE );
      };

      ii = null;
      while( ii < mm )
      {
# if PST_SVFRQ == 1
         fscanf( stream, "%s", ptr );
# endif
         fscanf( stream, "%s", ptr );
         sev11[ii][null] = strtod( ptr, endp );
         fscanf( stream, "%s", ptr );
         sev11[ii][ONE] = strtod( ptr, endp );
         ii++ ;
      };

      fclose( stream );
   
/*............................................................................*/
/* s_odd: */

     odd_mode:

      strcpy( eop, oddmde ); 

      strcpy( csp->rqlng, "Enter job index of '" );
      strcat( csp->rqlng, eop );
      strcat( csp->rqlng, "' mode S-prmtrs" );
      strcpy( csp->rqfrm, "brackets" );
      ( csp->dflng ) = ( csp->inlng ) + ONE;
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
         goto odd_mode;

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
         printf( "\n Incompatible S-parameter files s11[even], "
            "s11[odd] !!!" );
         printf( "\n [ The number of frequency points differs: %ld != %ld ]",
            mm, nn );
         printf( "\n\n " );
         exit( EXIT_FAILURE );
      };

      ii = null;
      while( ii < mm )
      {
# if PST_SVFRQ == 1
         fscanf( stream, "%s", ptr );
# endif
         fscanf( stream, "%s", ptr );
         sod11[ii][null] = strtod( ptr, endp );
         fscanf( stream, "%s", ptr );
         sod11[ii][ONE] = strtod( ptr, endp );
         ii++ ;
      };

      fclose( stream );
   };
/*............................................................................*/
   if(( option == ONE )||( option == THREE ))
   {
/* reflection: */

      ii = null;
      while( ii < mm )
      {
         scpl11[ii][null] = ( sev11[ii][null] + sod11[ii][null] )/2.;
         scpl11[ii][ONE] = ( sev11[ii][ONE] + sod11[ii][ONE] )/2.;

# if PST_PHASES == 1
         ( cpt->r ) = scpl11[ii][null];
         ( cpt->i ) = scpl11[ii][ONE];
/*............................................................................*/
         cpt = argc( cpt, PST_PBRNCH );      /* function arg(z) */
/*.........................................*//* [ branch PST_PBRNCH ] */
         if ( null == strstr( phunit, PST_PHUNIT ))
            phs11[ii] = 180.*( cpt->arg )/PI;
         else
            phs11[ii] = cpt->arg;
# endif
         ii++ ;
      };

      strcpy( fleptr, "cpl.s11" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "s11" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", zunit );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

# if PST_SVFRQ == 1
      ss = xx1;
# endif

      ii = null;
      while( ii < mm )
      {
# if PST_SVFRQ == 1
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

# if PST_PHASES == 1
      strcpy( fleptr, "cpl.phs11" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "phase_s11" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", PST_PHUNIT );
      fprintf( stream, "\n%+.15e", xx1 );
      fprintf( stream, "\n%+.15e", xx2 );
      fprintf( stream, "\n%+.15e", dx );
      fprintf( stream, "\n%ld", mm );

      ss = xx1;
      ii = null;
      while( ii < mm )
      {
         fprintf( stream, "\n%+.15e", ss );
         fprintf( stream, "\t%+.15e", phs11[ii] );
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

         PST_SPLCPY( 'p', EVL_GPHINTP ); /* input: spline(*) function */
                                         /* [ opt, 'p'ercent ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n\n Graphics require too many points !!!" );
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

         PST_SPLCPY( 'd', EVL_GPHINTP ); /* input: spline(*) function */
                                         /* [ opt, 'd'ecibel ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n\n Graphics require too many points !!!" );
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

# if PST_PHASES == 1 
/*............................................................................*/
/* graphics, phase s11: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = phs11[ii];

         PST_SPLCPY( 'r', null ); /* input: spline(*) function */
                                  /* [ opt, 'r'eal ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n\n Graphics require too many points !!!" );
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
         strcpy( gpt->yunit, PST_PHUNIT );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   
# endif

   };/* end if ( option == ... ) */
/*............................................................................*/
   if(( option == TWO )||( option == THREE ))
   {
/* transmission s21: */

      ii = null;
      while( ii < mm )
      {
         scpl21[ii][null] = ( sev11[ii][null] - sod11[ii][null] )/2.;
         scpl21[ii][ONE] = ( sev11[ii][ONE] - sod11[ii][ONE] )/2.;

# if PST_PHASES == 1
         ( cpt->r ) = scpl21[ii][null];
         ( cpt->i ) = scpl21[ii][ONE];
/*............................................................................*/
         cpt = argc( cpt, PST_PBRNCH );      /* function arg(z) */
/*.........................................*//* [ branch PST_PBRNCH ] */
         if ( null == strstr( phunit, PST_PHUNIT ))
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

# if PST_SVFRQ == 1
      ss = xx1;
# endif

      ii = null;
      while( ii < mm )
      {
# if PST_SVFRQ == 1
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

# if PST_PHASES == 1
      strcpy( fleptr, "cpl.phs21" );
      stream = fopen( fleptr, "w" );

      fprintf( stream, "%s", tlmmod );
      fprintf( stream, "\n%s", "phase_s21" );
      fprintf( stream, "\n%s", xunit );
      fprintf( stream, "\n%s", PST_PHUNIT );
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

         PST_SPLCPY( 'd', EVL_GPHINTP ); /* input: spline(*) function */
                                         /* [ opt, 'd'ecibel ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n\n Graphics require too many points !!!" );
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
         printf( "\n [ program stopped .]\n " );
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

# if PST_PHASES == 1
/*............................................................................*/
/* graphics, phase s21: */

      ss = xx1;
      nn = null;
      ii = null;
      while ( ii < mm )
      {
         real1 = phs21[ii];

         PST_SPLCPY( 'r', null ); /* input: spline(*) function */
                                  /* [ opt, 'r'eal ] */
         ii++ ;
      };

      if( GPH_POINTS <= nn )
      {
         printf( "\n\n Message from function %s:", DO_PSTPRC );
         printf( "\n\n Graphics require too many points !!!" );
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
         strcpy( gpt->file, "phs21" );
         strcpy( gpt->name, tlmmod );
         strcpy( gpt->text, "phase_s21" );
         strcpy( gpt->xunit, xunit );
         strcpy( gpt->yunit, PST_PHUNIT );

         ( gpt->nn ) = ( spt->nn );
/*............................................................................*/
         ind = graphp( gpt );         /* create graphics file:                */
/*.................................*//*                                       */
      };   
# endif
      ( csp->dfopt ) = FOUR;
   } /* end if ( option == ... ) */
   else if( option == FOUR )
   {
      PRBLDCLR( "\r" );
      printf( "\r %*s", 78, "PSTPRC" );
      PRNORMAL( "" );

      state = smithchrt( state );
      ( csp->dfopt ) = null;
   };
/*............................................................................*/

   if ( option != null )
   {
      strcpy( csp->tasks, "COMPUTE ... " );
      strcpy( csp->cmmnt, "Postprocessing, generate Smith chart etc.:" );
      strcpy( csp->envmt, "PSTPRC" );

      goto menu;
   };
      
   PRBLDCLR( "\r" );
   printf( "\r %*s", 78, "PSTPRC" );
   PRNORMAL( "" );

   ( state->rtn ) = null;
   return state;
}
/*============================================================================*/
# undef PST_PRECISION
# undef PST_PHASES
# undef PST_PBRNCH
# undef PST_PHUNIT
# undef PST_DIMNS
# undef PST_SCALE
# undef PST_SVFRQ
/****************** end of postprocessing function pstprc(*) ******************/
