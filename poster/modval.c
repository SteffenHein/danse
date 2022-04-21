# ifdef _POST_MODEL
/* SECTION COMPILED IN OPTION -D_POST_MODEL */
/*____________________________________________________________________________*/
/* [ function modval(*) ] */
# define DSC_MODEL "mod_ccoax.G"
/*******************************************************************************
*                                                                              *
*   ANSI-C function modval(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0.]                                                            *
*                                                                              *
*   This function determines special [ i.e. DSC model dependent ] evaluation   *
*   modes for files dsc.val<n>                                                 *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 22, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define MOD_DEFLT 7 /* default models */
/*----------------------------------------------------------------------------*/
# define MOD_FLOWS 1 /* evaluate flows - none:0, only vertical [y]:1, all:2 */
/*----------------------------------------------------------------------------*/
# define MOD_ENBLSKN 1 /* 1: enable [ skin effect ] heat source bnd faces */
# define MOD_ENBLESC 1 /* 1: enable envrmnt-surface heat conductivity */
/*----------------------------------------------------------------------------*/
/* excitation type: EXC_MAXW = "DIRAC_PULSE", "HARMONIC", "MULTIPLE_HARMONIC" */
/*
# define EXC_MAXW  "DIRAC_PULSE"
# define EXC_MAXW  "GAUSS_PULSE"
# define EXC_MAXW  "HARMONIC"
# define EXC_MAXW  "SMOOTH_HARMONIC"
# define EXC_MAXW  "OSCILLATING_GAUSS__" [e.g.]

# define EXC_HCRR  "PASSIVE____________"
# define EXC_HCRR  "PERIODIC_PULSE_____"
# define EXC_HCRR  "DOUBLE_PERIODIC____" [e.g.]
*/
# define EXC_MAXW  "OSCILLATING_GAUSS__"
# define EXC_HCRR  "PASSIVE____________"
/*----------------------------------------------------------------------------*/
# define MDV_PLOT_SWR 1 /*1: create standing waves plots at frequencies frq[j]*/
# define MDV_SPCLPRTS 3 /* number of evalt'd special ports: 0<= n <MAX_PERIOD */
# define MDV_AUTOEVAL 0 /*1: automatic port labeling for standing wave eval.  */
# define MDV_SPLINTPL 1 /*>0: spline interpolate spc.format [references etc.] */
# define MDV_GPHINTPL 5 /*>0: spline interpl. gph.format [gnu.graphics,e.g.]  */
# define MDV_NORMALZE 0 /*1: normalize S-parmts to abs.value 1 [ cautious ! ] */
/*----------------------------------------------------------------------------*/
# define _POSIX_SOURCE 1 /* some headers of the POSIX.1 standard will be used */
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
# include "./math/maths.h"  /* 'my' computation environment headers */
# include "./math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "./CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize this header for POSTER configuration: */
# include "./poster/POSTER.CONF" 
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
   # include <ncurses.h>
# endif 
/*----------------------------------------------------------------------------*/
# include "./tools/txctyp.h"  
# include "./poster/posttp.h"  /* typedefs: POSTSTATE, EVALUATE, SPLINES, etc.*/
/*----------------------------------------------------------------------------*/
# define LARGE_LOG_VAL ( 1.e+3 )
/*----------------------------------------------------------------------------*/
static SPLINES spl = {null};
static GRAPHICS gph = {null};
/*----------------------------------------------------------------------------*/
/* 'my_terminal' configuration: */

# if USE_NCURSES == 1
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
   # define CLSCREEN {\
     printf( "\f" );\
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
/* configure function pstprc(*): */
/*...........................................................................*/
# define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X' ] */
/*...........................................................................*/
# define PST_MESH  0      /* 1: plot mesh cell system */
# define PST_WALLS 2      /* 1: plot E-walls, 2: plot E- and M-walls */
# define PST_ISOSCLE 1    /* 1: equally scaled axes */
/*...........................................................................*/
# define PST_TMPPLT  4    /* Disable/enable temperature plot mode [0/1/...] */
# define PST_PRSPLT  3    /* Disable/enable pressure plot mode [0/1/...] */
# define PST_FLWPLT  4    /* Disable/enable velocity plot mode [0/1/...] */
# define PST_PRSINIT 1
/*...........................................................................*/
# if CMPRSSBL == 0
/* Disable density plot mode */
   # define PST_DNSPLT  0
   # define PST_DNSINIT 0
# else /* if CMPRSSBL != 0 */
/* Disable/enable density plot mode [0/1/...] */
   # define PST_DNSPLT  3
   # define PST_DNSINIT 1  
# endif /* CMPRSSBL != 0 */
/*...........................................................................*/
/* scale, vectors: */
/* PST_SCALE1: ratio of PST_SCLFLW or, if this is ZERO, of max | flow | to */
/* the geometric mesh extension */
# define PST_SCALE1 ( 1.000e-01 )

/* scale, vector pointers: */ 
/* PST_SCALE2: the ratio of the vector pointer size to the mesh extension */
# define PST_SCALE2 ( 1.000e-02 )
/*...........................................................................*/
# ifndef PST_SCLFLW /* scale of fluid flow [ maximum fluid velocity , e.g. ] */
# if MOD_DEFLT == 1
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 2
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 3
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 4
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 5
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 6
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 7
   # undef PST_PAUSE
   # define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X'] */
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 8
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 9
   # undef PST_PAUSE
   # define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X'] */
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 10
   # undef PST_PAUSE
   # define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X'] */
   # define PST_SCLFLW ( 1.000e-01 )
# else
   # define PST_SCLFLW ZERO
# endif /* MOD_DEFLT == ... */
# endif /* not defined PST_SCLFLW */
/*...........................................................................*/
# include "./tools/pstprc.fld"
/*============================================================================*/

POSTSTATE *
modval( POSTSTATE *state )
{
/* allusions: */
/*
   extern SPLINES spl;
   extern GRAPHICS gph;
*/
/* declarations: */

   static FILE
      *evalfle = NULL,
      *timefle[MAX_PERIOD],
      *specfle[MAX_PERIOD];

   static FFT *fpt;
   static TXCNSL *csp;
   static SPLINES *spt = &spl;
   static GRAPHICS *gpt = &gph;
   static EVALUATE *vpt;

   static char     /* operation marks                                         */
      opm1 = null, /* opm1 = [0] 1: reference file [not] charged              */
      opm2 = null; /* opm2 = [0] 1: reference is [not] mean inc. wave         */

   static const char
     *timefrm = "created: %.24s ";

   static short
      j_ = null,
      k_ = null,
      ll = null,
      ind = null,
      lbl = null,
      frqlbl = null,

   envl[MAX_PERIOD] = {null};

   static long 
      hh = null,
      h_ = null,
      ii = null,
      i_ = null,
      jj = null,
      kk = null,
      mm = null,
      nn = null,
      pp = null,
      qq = null,
      init = null,
      final = null;

   static short 
     *idx = NULL;

   static const char 
     *spcpfx = "spc.",
     *option = "forward", /* Fourier transformation option */
     *spformat = "%s\n",
     *ipformat = "%ld\n",
     *scformat = "%80s";

   static char 
      ptr[STS_SIZE] = {null},
      tmeptr[STS_SIZE] = {null},
      spcfle[STS_SIZE] = {null},
      dpformat[STS_SIZE] = {null},
      absc_unit[STS_SIZE] = {null},
      ordn_unit[STS_SIZE] = {null},
      t_type[MAX_PERIOD][STS_SIZE] = {{null}},
      t_text[MAX_PERIOD][STS_SIZE] = {{null}},
      s_type[MAX_PERIOD][STS_SIZE] = {{null}},
      s_text[MAX_PERIOD][STS_SIZE] = {{null}},
      fleptr[MAX_PERIOD][STS_SIZE] = {{null}},
    **endp = null;

   static double
      xx = ZERO,
      yy = ZERO,
      zz = ZERO,
      x1 = ZERO,
      x2 = ZERO,
      dx = ZERO,
      df = ZERO,
      norm = ZERO,
      real = ZERO,
      imag = ZERO,
      evl0 = ZERO,
      evl1 = ZERO,
      dlt0 = ZERO,
      dlt1 = ZERO,
      minm = ZERO,
      maxm = ZERO,
      xlower = ZERO,
      xupper = ZERO,
      frequency = ZERO;

   static double 
      frq[EVL_SPF] = {ZERO},
      vswr[EVL_SPF] = {ZERO},
      refl[EVL_SPF] = {ZERO},
      mean[EVL_SPF] = {ZERO},
      rref[EVL_SPF] = {ZERO},
      iref[EVL_SPF] = {ZERO},
      absmaxm[MAX_PERIOD] = {ZERO},
      timemax[MAX_PERIOD] = {ZERO},
      absmean[MAX_PERIOD] = {ZERO},
      meanenv[MAX_PERIOD] = {ZERO},
      maxspec[MAX_PERIOD] = {ZERO},
      freqmax[MAX_PERIOD] = {ZERO},
      modulation[MAX_PERIOD] = {ZERO},
      rspc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      ispc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      aspc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      rcpl[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      icpl[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      acpl[MAX_PERIOD][EVL_SPF] = {{ZERO}};

   static time_t
      nseconds = null,
     *timer = null;

/* system parameters and function prototypes */

   time_t
      time( time_t *timer );

   char
     *ctime( const time_t *timer );

# ifndef _CCBUG
   char
     *strcpy( char *ptr1, const char *ptr2 ),
     *strcat( char *ptr1, const char *ptr2 ),
     *strncat(char *ptr1, const char *ptr2, size_t n );
# endif

/* mathematical function prototypes */

   double 
      sqrt( double x ),
      log10( double x ),
      ceil( double x ),
      floor( double x );

/* user defined function prototypes: */

   FFT *
      fftrf( FFT *fpt );

   int 
      dspval( void );

   SPLINES
     *spline( SPLINES *spt );

   int 
      graphp( GRAPHICS *gpt );

   POSTSTATE 
     *pstprc( POSTSTATE *stp );

   char
     *lotos( long, char );

   TXCNSL 
     *txcnsl( TXCNSL *csp );

   EVALUATE
     *readval( POSTSTATE *state, char option );
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
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
/* string assignments [ initializations ], memory allocation etc. */
   
/* double format: */

   if ( null == strncmp( PLOT_FORMAT, "SPLINE" ,THREE ))
      strncpy( dpformat ,"%+.15E%s", TEN );
   else
      strncpy( dpformat , "%+.15e%s", TEN );

   strcpy( gph.format, PLOT_FORMAT );
/*............................................................................*/
/* initialize text console: */

   csp = txcnsl( null );
/*............................................................................*/
/* Error messages: */

# if LINUX_C_SNTX != 1
   hh = strlen ( LNGSTR );
   if ( hh < MAX_PERIOD * ( DEC + ONE ) ) /* <- length of long string LNGSTR  */
   {                                      /*    initializer in 'consts.h'     */
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD > %d !!!", 250/DEC );
      printf( "\n [ Filename pointer array 'fleptr[]' may perturb memory.]" );
   };
# endif

   mm = MAX_PERIOD;
   if ( FTR_NMBR < mm )                /* <- FTR_NMBR: macro in POSTER.CONF */
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD = %d !!!", MAX_PERIOD );
      printf( "\n [ Required is MAX_PERIOD <= %d = FTR_NMBR "
        "<- macro defined in FORMER.CONF.]", FTR_NMBR );
   };

   if ( FTR_NMBR < mm )                /* <- FTR_NMBR: macro in POSTER.CONF */
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD = %d !!!", MAX_PERIOD );
      printf( "\n [ Required is MAX_PERIOD <= %d = FTR_NMBR "
              "<- macro defined in 'dmnsnd.h'.]", FTR_NMBR );
   };
/*............................................................................*/
/* initialize pointers to structures [ of type EVALUATE and FFT ]: */

   vpt = ( state->vpt );
   fpt = ( vpt->fpt );
   idx = ( state->idx );
/*............................................................................*/
/* evaluation mode lbl = 1,2,... */

   if ( null == strncmp( vpt->exctyp, "HARMONIC", FIVE ))
      lbl = 4;
   else if ( null == strncmp( vpt->exctyp, "SMOOTH_HARMONIC", FIVE ))
      lbl = 4;
   else if ( null == strncmp( vpt->exctyp, "STEADY_STATE", FIVE ))
      lbl = 5;
   else
      lbl = 6;
/*............................................................................*/
   PRBLDCLR( "\n" );
   printf( " %*s", 78, "MODVAL" );
   PRNORMAL( "" );

   if (( ONE < lbl )
     &&( lbl < SEVEN ))
   {
/* enter special ports 1,2, ... */

      ii = null;
      while(( ii < ( state->fldprd ))
          &&( ii < MDV_SPCLPRTS )
          &&( ii < MAX_PERIOD ))
      {
         ii++;
         idx[ii] = ii;
      };
      idx[null] = ii;

      while (( ii < ( state->fldprd ))
           &&( ii < MAX_PERIOD ))
      {
         ii++;
         hh = ii - MDV_SPCLPRTS;
         h_ = hh - ONE;
/*............................................................................*/
# if MDV_AUTOEVAL == 1
         idx[ii] = ii;
# else
/*............................................................................*/
         if ( hh == ONE )
            printf( "\n VSWR evaluation: Please enter indices "
               "[ Escape: enter null ] >---->\n" );
        next_index:

         printf( " >----> enter %6ld. index ..................."
            "....................: ", hh );
         scanf( "%s", ptr );
         idx[ii] = ( short ) strtol( ptr, endp, DEC );
/*............................................................................*/
/* check if index is in domain: */

         if (( state->fldprd ) < idx[ii] )
         {
            printf( "\n This is not a Maxwell field port "
               "in selected file %s !!!", ( state->file ));
            if ( idx[ii] <= ( state->period ))
               printf( "\n [ It is a heat or fluid port that can be "
                  "separately evaluated.]" );
            printf( "\n\n >----> Enter new index >------"
               "------------------------------------>\n");
            goto next_index;
         };
# endif /* if MDV_AUTOEVAL != 1 */
/*............................................................................*/
         if ( null < idx[ii] )
         {
            if (( ONE < lbl )
              &&( lbl < FOUR ))
            {
               printf( " >----> enter name of file to be created "
                  "for index %-6d..........: ", idx[ii] );
               scanf( "%s", ptr );
               strcpy( fleptr[h_], ptr );
            };
            idx[null] = ii;

            if ( MAX_PERIOD <= idx[null] )
            {
               printf( "\n The last evaluation index %d has been "
                  "accepted.", idx[idx[null]] );
               printf( "\n [ Maximum number is %ld ",
                  ( long ) MAX_PERIOD );
               printf( "= macro MAX_PERIOD "
                  "in header '%s'.]\n", "POSTER.CONF" );

              goto read_values;
            };
         }
         else /* if ( idx[ii] <= null ) */
            goto read_values;
      }; 

     read_values:
/*............................................................................*/
      vpt = readval( state, 'r' );       /*                                   */
/*.....................................*/
   }; /* end if (( ONE < lbl )&&( lbl < SEVEN )) */
/*........... values dsc.val<n> read & copied for selected ports ...........*/
/* evaluation [ options labelled lbl ]: */

   switch ( lbl )
   {
     case null:
      return state;

     case 1:

      printf( "\n" );
      break;

     case 2: /* >------ save time response ---------------------------------> */

      printf( "\n\n Please wait a moment !" );
      printf( "\n [ Writing data into files %s ... ]\n",
         fleptr[null] );

      for ( hh=ONE; hh<=idx[null]; hh++ )
      {
         h_ = hh - ONE;

         strncpy( t_type[h_], ( vpt->name ), STS_SIZE );
         strncpy( t_text[h_], ( vpt->text ), STS_SIZE );

         timefle[h_] = fopen( fleptr[h_], "w+" );

         fprintf( timefle[h_], spformat, t_type[h_] );
         fprintf( timefle[h_], "%s%s%d%s\n", t_text[h_], "_[evaluated"
            "_index", idx[hh], "]" );
         fprintf( timefle[h_], spformat, absc_unit );
         fprintf( timefle[h_], spformat, ordn_unit );

         fprintf( timefle[h_], dpformat, ( fpt->t[null] ), "\n" );
         fprintf( timefle[h_], dpformat, ( fpt->tt[null] ), "\n" );
         fprintf( timefle[h_], dpformat, ( fpt->dt[null] ), "\n" );
         fprintf( timefle[h_], ipformat, ( fpt->ttlg[null] ));

         for ( i_=null ; i_<( fpt->ttlg[null] ); i_++ )
         {
            fprintf( timefle[h_], dpformat, ( fpt->r[hh][i_] ), "  " );
            fprintf( timefle[h_], dpformat, ( fpt->i[hh][i_] ), "\n" );

# if LIST_FILE == 1
            printf( "  %ld:   %1.16e  + ( %1.16e )*j \n",
                    i_, ( fpt->r[hh][i_] ), ( fpt->i[hh][i_] ));
# endif
         };

         nseconds = time( timer );
         strcpy( tmeptr, ctime( &nseconds ));
         fprintf( timefle[h_], "\n%s %s %s %.24s\n", "DSC time response file",
            fleptr[h_], "terminated:", tmeptr );
         printf( " %s %s %s %.24s", "DSC time response file",
            fleptr[h_], "terminated:", tmeptr );

         fclose( timefle[h_] );
      };/* next hh */

      break;
/*...................... end case2 [ save time response ] ....................*/

     case 3: /* >------------ save spectral response -----------------------> */
     case 6: /* >---------- evaluate spectral response ---------------------> */

      if ( null < idx[null]  )
      {
         if (( lbl == SIX ) 
           &&( null != strncmp(( vpt->text ), "reference", NINE )))
         {                    /* enter external reference spectrum  'spc.ref' */
            strcpy( spcfle, spcpfx );
            strcat( spcfle, "ref" );

            printf( "\n" );

           open_reference1:

            evalfle = fopen( spcfle, "r+" );

            if ( evalfle == null )
            {
               printf( "\n Reference file %s not found "
                  "in present directory:", spcfle );
               printf( "\n Please re-enter filename [ Escape: "
                  "enter null ] >----> " );
               scanf( "%s", spcfle );

               if ( *spcfle == '0' )
                  goto frq_interval;
               else
                  goto open_reference1;
            }
            else if ( evalfle != null )
            {
               printf( "\n opened: reference spectrum file %s", spcfle ); 
               fscanf( evalfle, scformat, ptr );

               if ( null != strncmp( ptr, ( vpt->name ), THREE ))
               {
                  printf( "\n\n Error message from function %s :", __func__ );
                  printf( "\n\n Incompatible system identifier '%s'", ptr );
                  printf( "\n in reference spectrum, file %s !!!", spcfle );
                  printf( " [ overriding.]\n" );

                  fclose( evalfle );

                  goto frq_interval;
               };
       
               fscanf( evalfle, scformat, ptr );
               fscanf( evalfle, scformat, absc_unit );
               fscanf( evalfle, scformat, ordn_unit );

/* the domain lower bound: */
               fscanf( evalfle, scformat, ptr );
               xlower = strtod( ptr, endp ); /* the lower frequency bound */
	        
/* the domain upper bound: */
               fscanf( evalfle, scformat, ptr );
               xupper = strtod( ptr, endp ); /* the upper frequency bound */

/* the increment: */
               fscanf( evalfle, scformat, ptr );
/*............................................................................*/
# if MDV_SPLINTPL == 0
               df = strtod( ptr, endp );
# else
               dx = strtod( ptr, endp );
# endif
/*............................................................................*/
/* the number of sample points: */

               fscanf( evalfle, scformat, ptr );
               frqlbl = strtol( ptr, endp, DEC );

               if ( EVL_SPF < frqlbl )
                  frqlbl = EVL_SPF;

               if ( SPL_INTPL < frqlbl )
                  frqlbl = SPL_INTPL;
               
               jj = null;
               while( jj < frqlbl )
               {
/*............................................................................*/
# if EVL_WRTFREQ == 1
                  fscanf( evalfle, scformat, ptr ); /* read the frequency */
# endif
/*............................................................................*/
                  fscanf( evalfle, scformat, ptr ); /* the real part */
                  rref[jj] = strtod( ptr, endp );
                  fscanf( evalfle, scformat, ptr ); /* the imaginary part */
                  iref[jj] = strtod( ptr, endp );

                  jj++ ;   
               };
               printf( "\r entered: reference spectrum %s .      \n",
                  spcfle );
               fclose( evalfle );
               opm1 = ONE;
            }
            else /* case: no reference spectrum found or lbl != SIX */
               goto frq_interval;
         }
         else /* case: create reference spectrum */
              /* specify frequency domain */
         {
           frq_interval:

            opm1 = null;

            strncpy( ordn_unit, ( vpt->yunit ), STS_SIZE );
            strcat( ordn_unit, "*" );
            strncat( ordn_unit, ( vpt->xunit ), STS_SIZE );
            strcpy( absc_unit, "1/" );
            strcat( absc_unit, ( vpt->xunit ));

            printf( "\n\n Please enter frequency interval "
               "[ f1, f2 ] / GHz\n\n" );

            strcpy(( csp->rqfrm ), "points" );
            strcpy(( csp->rqdbl ), "Enter lower frequeny " );
            strcat(( csp->rqdbl ), ">>---------------------------> f1 / GHz" );
/*............................................................................*/
            csp = txcnsl( csp );        /* text console                       */
/*....................................*/
            xlower = ( csp->indbl );
            xlower *= ( 1.0e+09 );

            strcpy(( csp->rqfrm ), "points" );
            strcpy(( csp->rqdbl ), "Enter upper frequeny " );
            strcat(( csp->rqdbl ), ">>---------------------------> f2 / GHz" );
/*............................................................................*/
            csp = txcnsl( csp );        /* text console                       */
/*....................................*/
            xupper = ( csp->indbl );
            xupper *= ( 1.0e+09 );

            printf( "\n frequency interval: [ %1.5e , %1.5e ] "
               "GHz\n\n", ( 1.0e-09*xlower ), ( 1.0e-09*xupper ));

            strcpy(( csp->rqfrm ), "bracket" );
            strcpy(( csp->rqstr ), "Input correct >>---------" );
            strcat(( csp->rqstr ), "------------------------> [ Y/n ] ?" );
            strcpy(( csp->dfstr ), "y" );
/*............................................................................*/
            csp = txcnsl( csp );        /* text console                    */
/*....................................*/
            strcpy( ptr, ( csp->instr ));

            if (( *ptr == 'n' )||( *ptr == 'N' )) 
               goto frq_interval;
         };
      };
/* evaluation parameters [ reference spectrum, frequency domain ] specified.  */
/*............................................................................*/










/*............................................................................*/
/* Fourier transformations: */

      printf( "\n Please wait a moment !" );

      if ( lbl == THREE )
         printf( "\n [ Fourier transforms "
            "- writing data on files '%s' ... ]\n", fleptr[null] );
      if ( lbl == SIX )
         printf( "\n [ Fourier transforms ]\n" );

      kk = ( fpt->ttlg[null] );

      nn = TWO;
      while (( nn < kk )&&( nn < FTR_SIZE ))
         nn *= TWO;

      mm = nn / TWO;
      pp = ( nn - kk ) / TWO;
      qq = pp + kk;

      hh = ONE; 
      while( hh <= idx[null] )
      {
         h_ = hh - ONE;

         ii = nn - ONE;
         while ( qq <= ii )
         {
            ( fpt->r[hh][ii] ) = ZERO;
            ( fpt->i[hh][ii] ) = ZERO;
            ii-- ;
         };
         while ( pp <= ii )
         {
            ( fpt->r[hh][ii] ) = ( fpt->r[hh][ii-pp] );
            ( fpt->i[hh][ii] ) = ( fpt->i[hh][ii-pp] );
            ii-- ;
         };
         while ( null <= ii )
         {
            ( fpt->r[hh][ii] ) = ZERO;
            ( fpt->i[hh][ii] ) = ZERO;
            ii-- ;
         };
/*............................................................................*/
# if READ_REVERSE_ == 1 /* interchange upper [negative] with lower [positive] */
                        /* frequencies: */
         ii = pp; do
         {
            jj = ii + mm;

            xx = ( fpt->r[hh][jj] );
            yy = ( fpt->i[hh][jj] );
	    
            ( fpt->r[hh][jj] ) = ( fpt->r[hh][ii] );
            ( fpt->i[hh][jj] ) = ( fpt->i[hh][ii] );

            ( fpt->r[hh][ii] ) = xx;
            ( fpt->i[hh][ii] ) = yy;

            ii++ ;
         }  while ( ii < mm );

         ( fpt->t[hh] ) = - mm*( fpt->dt[null] );
# else
         ( fpt->t[hh] ) = ZERO;
# endif /* READ_REVERSE */
/*............................................................................*/
         ( fpt->tt[hh] ) = ( fpt->t[hh] ) + nn*( fpt->dt[null] );
         ( fpt->dt[hh] ) = ( fpt->dt[null] );
         ( fpt->ttlg[hh] ) = nn;
         ( fpt->mult[hh] ) = ONE;
         ( fpt->p ) = hh;
         ( fpt->q ) = hh;
         strcpy( fpt->opt, option);
/*............................................................................*/
         fpt = fftrf( fpt );      /* Fast Fourier transform                   */
/*..............................*/
         hh++ ;
      }; /* next hh */

      dx = ( fpt->ds[null] );
/*............................................................................*/
# if WRITE_REVERSE == 1
      ii = mm;
      xx = ( fpt->s[null] );
# else
      ii = null;
      xx = ZERO;
# endif
/*............................................................................*/
      if (( MDV_SPLINTPL == 0 )
        &&( opm1 == ONE ))
      {
         x1 = xlower + dx/2.;
         x2 = xupper - dx/2.;
      }
      else
      {
         x1 = xlower;
         x2 = xupper;
      };

      if( xx < x1 )
         kk = ( long ) floor (( x1 - xx ) / dx ); /* kk = largest integer not */
      else                                        /* greater than argument */
         kk = null;

      if( nn < kk )
         kk = nn;
      
      init = ii + kk;
      x1 = xx + kk*dx; /* [1]; cf. remark [2], below */

      if( xx < x2 )
         kk = ( long ) ceil (( x2 - xx ) / dx ); /* kk = smallest integer not */
      else                                       /* less than argument */
         kk = null;

      if( nn < kk )
         kk = nn;

      final = ii + kk + ONE;
      x2 = xx + kk*dx;

      kk = final - init;

/* Fourier transformations terminated */
/*............................................................................*/










/*............................................................................*/
      if ( lbl == THREE ) /* save spectral response */ 
      {
         hh = ONE;
         while( hh <= idx[null] )
         {
            h_ = hh - ONE;

            strncpy( s_type[h_], ( vpt->name ), STS_SIZE );
            strncpy( s_text[h_], ( vpt->text ), STS_SIZE );

            specfle[h_] = fopen( fleptr[h_], "w+" );

            fprintf( specfle[h_], spformat, s_type[h_] );
            fprintf( specfle[h_], "%s%s%d%s\n", s_text[h_],"_[evaluated"
               "_index", idx[hh], "]" );
            fprintf( specfle[h_], spformat, absc_unit );
            fprintf( specfle[h_], spformat, ordn_unit );
            fprintf( specfle[h_], dpformat, x1, "\n" );
            fprintf( specfle[h_], dpformat, x2, "\n" );
            fprintf( specfle[h_], dpformat, dx, "\n" );
            fprintf( specfle[h_], ipformat, kk );

# if EVL_WRTFREQ == 1
            frequency = x1;
# endif
            jj = init;
            ii = init;
            while ( ii < final )
            {
               if ( nn <= jj )
                  jj -= nn;
 
# if EVL_WRTFREQ == 1
               fprintf( specfle[h_], dpformat, frequency, "  " );
	       frequency += dx;
# endif
               fprintf( specfle[h_], dpformat, ( fpt->r[hh][jj] ), "  " );
               fprintf( specfle[h_], dpformat, ( fpt->i[hh][jj] ), "\n" );

# if LIST_FILE == 1
               printf( "  %d:   %1.16le  + ( %1.16le )*j \n",
                  jj, ( fpt->r[hh][jj] ), ( fpt->i[hh][jj] ));
# endif
               jj++ ;
               ii++ ;
            };

            nseconds = time( timer );
            strcpy( tmeptr, ctime( &nseconds ));
            fprintf( specfle[h_], "\n%s %s %s %.24s\n", "DSC spectrum file",
               fleptr[h_], "terminated:", tmeptr );

            printf( "\n\n %s %s %s %.24s", "DSC spectrum file",
               fleptr[h_], "terminated:", tmeptr );

            fclose( specfle[h_] );

            hh++ ;
         };
         break;
      };/* end if lbl == THREE */
/*............................................................................*/
/* spectral response saved */
/*............................................................................*/
/* still case 6 [ cont'd ] - spectral response; spline interpolation:         */
/*............................................................................*/
# if MDV_SPLINTPL == 0
/*
      printf( "\n frqlbl = %.7d, kk = %.7d", frqlbl, kk );
      printf( "\n xlower = %.15le , x1 = %.15le", xlower, x1 );
      printf( "\n xupper = %.15le , x2 = %.15le", xupper, x2 );
      printf( "\n df     = %.15le , dx = %.15le", df, dx );
      scanf( "%s", ptr );
*/
      if( opm1 == ONE )
      {
         if(( frqlbl != kk )||
            ( 1.e-7 < fabs( ( df - dx ) / dx ))||
            ( 1.e-7 < fabs( ( xlower - x1 ) / dx ))||
            ( 1.e-7 < fabs( ( xupper - x2 ) / dx )))
         {
            printf( "\n\n Error message from function %s :", __func__ );
            printf( "\n\n Evaluated file does not match reference file" );
            printf( " !!!\n " );
            break;
         };
      };

      frqlbl = kk;
      xlower = x1;
      xupper = x2;

# elif MDV_SPLINTPL < 10

      if( opm1 == null ) /* no reference file entered [ yet generated, e.g. ] */
      { 
         frqlbl = kk*MDV_SPLINTPL;
         if ( EVL_SPF < frqlbl )
            frqlbl = EVL_SPF;
      };

      dx = ( xupper - xlower ) / ( frqlbl - ONE );

# elif MDV_SPLINTPL >= 30

      if( opm1 == null ) /* no reference file entered [ yet generated, e.g. ] */
      {
         frqlbl = MDV_SPLINTPL;
         if ( EVL_SPF < frqlbl )
            frqlbl = EVL_SPF;
      };

      dx = ( xupper - xlower ) / ( frqlbl - ONE );
# endif
/*............................................................................*/
      if ( SPL_SPNTS < kk )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many supp. points in frequency domain !!!" );
         fprintf( stderr, "\n [ Number %ld exceeds maximum number of spline "
            "interpolation points %ld", kk, (long) SPL_SPNTS );
         fprintf( stderr, "\n   = macro SPL_INTPL in function %s.", "spline(*)" );
         fprintf( stderr, "\n   - Change macro in compliance with memory "
            "resources.]\n" );

         exit( EXIT_FAILURE );
      }
      else if ( SPL_INTPL < frqlbl )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many number of frequency points !!!" );
         fprintf( stderr, "\n [ Number %ld exceeds maximum number of spline "
            "interpolation points %ld", ( long ) frqlbl, (long) SPL_INTPL );
         fprintf( stderr, "\n   = macro SPL_INTPL in function %s.", "spline(*)" );
         fprintf( stderr, "\n   - Change macro in compliance with memory "
            "resources.]\n" );

         exit( EXIT_FAILURE );
      }
      else if ( EVL_SPF < frqlbl )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many number of frequency points !!!" );
         fprintf( stderr, "\n [ Number %ld exceeds maximum number %ld",
            ( long ) frqlbl, (long) EVL_SPF );
         fprintf( stderr, "\n   = macro EVL_SPF in function %s.", __func__ );
         fprintf( stderr, "\n   - Change macro in compliance with memory "
            "resources.]\n" );

         exit( EXIT_FAILURE );
      };

      if( xlower < x1 ) /* [2]: this can only happen due to roundoff error in */
         xlower = x1;   /* formula [1].                                       */
      if( x2 < xupper ) /* [ analogous statement ]                            */
         xupper = x2;

      df = ( xupper - xlower )/( frqlbl - ONE );
      frequency = xlower;

      ii = ONE;
      jj = frqlbl-ONE;
      frq[null] = xlower;
      spl.dmn[null] = xlower;
      while ( ii < jj ) 
      {
         frequency += df;
         frq[ii] = frequency;
         spl.dmn[ii] = frequency;
         ii++;
      };                    /* [3]; spl.dmn[] must lie strictly within the    */
      frq[ii] = xupper;     /* closed interval                                */
      spl.dmn[ii] = xupper; /* [ spl.vct[0][0], spl.vct[frqlbl-ONE][0] ],     */
                            /* which ist granted by virtue of [2] above       */
                            /* and condition [4], below.                      */
/*............................................................................*/
      hh = ONE;
      while( hh <= idx[null] )
      {
         h_ = hh - ONE;
         
/* spectrum, port hh [real part]: */

         frequency = x1;
         jj = init;
         ii = init;
         while ( ii < final )
         {
            if ( nn <= jj )
               jj -= nn;

            spl.vct[ii-init][null] = frequency;
            spl.vct[ii-init][ONE] = ( fpt->r[hh][jj] );

            frequency += ( fpt->ds[hh] );

            jj++ ;
            ii++ ;
         };
         spl.vct[kk-ONE][null] = x2; /* [4]: cf. comment [3] above */
/*............................................................................*/
# if MDV_SPLINTPL == 0

         j_ = null; 
         while( j_ < frqlbl )
         {
            rspc[h_][j_] = spl.vct[j_][ONE];
            j_++ ;
         };
# else
         spl.nn = frqlbl;
         spl.mm = kk;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation:                 */
/*.................................*//* real part                             */
         if (( spt->rtn ) == ONE )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on calling spline function" );
            printf( "\n for computing rspc[%ld][] !!!", h_ );
            printf( "\n [ program stopped.]\n" );

            exit( EXIT_FAILURE );
         };

         j_ = null; 
         while( j_ < frqlbl )
         {
            rspc[h_][j_] = spl.fct[j_];
            j_++ ;
         };
# endif
/*............................................................................*/
/* spectrum, port hh [maginary part]: */

         jj = init;
         ii = init;
         while ( ii < final )
         {
            if ( nn <= jj )
               jj -= nn;
               
            spl.vct[ii-init][ONE] = ( fpt->i[hh][jj] );

            jj++ ;
            ii++ ;
         };
/*............................................................................*/
# if MDV_SPLINTPL == 0

         j_ = null;
         while( j_ < frqlbl )
         {
            ispc[h_][j_] = spl.vct[j_][ONE];
            j_++ ;
         };
# else
         spl.nn = frqlbl;
         spl.mm = kk;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation:                 */
/*.................................*//* imaginary part                        */
         if (( spt->rtn ) == ONE )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on calling spline function" );
            printf( "\n for computing ispc[%ld][] !!!", h_ );
            printf( "\n [ program stopped.]\n" );

            exit( EXIT_FAILURE );
         };

         j_ = null; 
         while ( j_ < frqlbl )
         {
            ispc[h_][j_] = spl.fct[j_];
            j_++ ;
         };
# endif /* MDV_SPLINTPL != 0 */
/*............................................................................*/
/* spectrum, port hh [absolute value]: */

         jj = init;
         ii = init;
         while ( ii < final )
         {
            if ( nn <= jj )
               jj -= nn;

            real = ( fpt->r[hh][jj] );
            imag = ( fpt->i[hh][jj] );
            spl.vct[ii-init][ONE] = sqrt( real*real + imag*imag );

            jj++ ;
            ii++ ;
         };

         maxspec[h_] = ZERO;
         freqmax[h_] = ZERO;
/*............................................................................*/
# if MDV_SPLINTPL == 0

         j_ = null;
         while ( j_ < frqlbl )
         {
            norm = fabs( spl.vct[j_][ONE] );
            aspc[h_][j_] = norm;

            if ( maxspec[h_] < norm )
            {
               maxspec[h_] = norm;
               freqmax[h_] = spl.vct[j_][null];
            };
            j_++ ;
         };
# else
         spl.nn = frqlbl;
         spl.mm = kk;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation:                 */
/*.................................*//* absolute value                        */
         if (( spt->rtn ) == ONE )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on calling spline function" );
            printf( "\n for computing aspc[%ld][] !!!", h_ );
            printf( "\n [ program stopped.]\n" );

            exit( EXIT_FAILURE );
         };

         j_ = null; 
         while ( j_ < frqlbl )
         {
            norm = fabs( spl.fct[j_] );
            aspc[h_][j_] = norm;

            if ( maxspec[h_] < norm )
            {
               maxspec[h_] = norm;
               freqmax[h_] = spl.dmn[j_];
            };
            j_++ ;
         };
# endif /* MDV_SPLINTPL != 0 */
/*............................................................................*/
         freqmax[h_] /= 1.e+09; /* in GHz */

         hh++ ;

      };/* next hh */
/*......................... spectral response ready ..........................*/










/*............................................................................*/
/* generate reference spectrum: */
/*
      printf( "\n\n                                   "
              "                                  ) <- ?" );
      printf( "\r Generate reference spectrum ? >----"
              "----------------> [ y/n ] >--> (" );
      scanf( "%s", ptr );

      if (( *ptr == 'y' )||( *ptr == 'Y' ))
      {
*/ /* } */
      printf( "\n" );

      if ( null == strncmp(( vpt->text ), "reference", NINE ))
      {
         strcpy( spcfle, spcpfx );
         strcat( spcfle, "ref" );
         strcat( spcfle, ( state->flbl ));

         evalfle = fopen( spcfle, "w+" );

         if ( evalfle == null )
         {
            printf( "\n Error on opening reference spectrum file"
               " %s !", spcfle );
            printf( "\n [overriding: spectrum not saved.]\n " );
            break;
         };

         fprintf( evalfle, spformat, ( vpt->name ));
         fprintf( evalfle, spformat, "reference_spectrum" );
         fprintf( evalfle, spformat, absc_unit );
         fprintf( evalfle, spformat, ordn_unit );
         fprintf( evalfle, dpformat, xlower, "\n" );
         fprintf( evalfle, dpformat, xupper, "\n" );
         fprintf( evalfle, dpformat, dx, "\n" );
         fprintf( evalfle, ipformat, frqlbl );

# if EVL_WRTFREQ == 1
	 frequency = xlower;
# endif
         j_ = null;
         while ( j_ < frqlbl )
         {
# if EVL_WRTFREQ == 1
            fprintf( evalfle, dpformat, frequency, "  " );
	    frequency += dx;
# endif
            fprintf( evalfle, dpformat, rspc[null][j_], "  " );
            fprintf( evalfle, dpformat, ispc[null][j_], "\n" );
            j_++ ;
         };

         nseconds = time(timer);
         strcpy( tmeptr, ctime( &nseconds ));

         fprintf( evalfle, "\n# SPECTRUM file %s\n", spcfle );
         fprintf( evalfle, "# created:%s", tmeptr );

         fclose( evalfle );

         printf( "\n" );
         printf( CLEAR_LINE );
         printf( "\r REFSPC.file %s ", spcfle );
         printf( timefrm, tmeptr );

         if ( idx[null] <= MDV_SPCLPRTS )
            break;
      }; /* end if ( null == strncmp( vpt->text, "reference", NINE )) */
/*............................................................................*/
/* reference spectrum ready [ and stored ] */
/*............................................................................*/











/*............................................................................*/
/* display on screen [ console ]: */

      printf( "\n\n Special evaluated spectral parameters "
         "of file %s :", ( state->file ));
      printf( "\n\n  index|abs. maximum :at frequency " );

      for ( j_= null; (( j_< frqlbl )
        &&( j_< THREE )); j_++ )
         printf( "|ampl.at GHz->" );

      printf( "\n       |[%.13s]       [%.3s]", ordn_unit, "GHz" );

      for ( j_= null; (( j_< frqlbl )
        &&( j_< THREE )); j_++ )
         printf( "|%.7e", frq[j_]/1.e+09 );

      printf( "\n --------------------------------------"
              "----------------------------------------" );

      hh = ONE;
      while( hh <= idx[null] )
      {
         h_ = hh - ONE;
         printf( "\n %6d|%.7e:%.7e", idx[hh], maxspec[h_], freqmax[h_] );

         for ( j_=null ; (( j_< frqlbl )&&( j_< THREE )) ; j_++ )
            printf( "|%.7e", aspc[h_][j_] );

         hh++ ;
      };/* next hh */

      jj = ( int ) (( frqlbl - THREE + ( FIVE - ONE )) / FIVE );

      ii = null;
      while ( ii < jj )
      {
         j_ = THREE + ii*FIVE;
         if ( j_ < frqlbl )
         {
            printf( "\n\n  index" );
            while (( j_ < THREE + ( ii + ONE )*FIVE )&&( j_ < frqlbl ))
            {
               printf( "|ampl.at GHz->" );
               j_++;
            };
         };
         j_ = THREE + ii*FIVE;
         if ( j_ < frqlbl )
         {
            printf( "\n       " );
            while (( j_ < THREE + ( ii + ONE )*FIVE )&&( j_ < frqlbl ))
            {
               printf( "|%.7e", ( 1.e-09*frq[j_] ));
               j_++;
            };
         };
         printf( "\n -----------------------------------"
            "-------------------------------------------" );

         hh = ONE;
         while( hh <= idx[null] )
         {
            h_ = hh - ONE;
            j_ = THREE + ii*FIVE;
            if ( j_ < frqlbl )
            {
               printf( "\n %6d", idx[hh] );
               while (( j_ < THREE + ( ii + ONE )*FIVE )&&( j_ < frqlbl ))
               {
                  printf( "|%.7e", aspc[h_][j_] );
                  j_++;
               };
            };
            hh++ ;
         };
         ii++;
      };/* end while ii < jj */

      printf( "\n -----------------------------------"
         "-------------------------------------------" );
/*............................ display ready ................................*/










/*............................................................................*/
      opm2 = ONE;

      if (( opm1 == null )
        &&( MDV_SPCLPRTS < idx[null] ))  
      {
         printf( "\n\n normalizing S-parameters to absolute value"
            "\n of incident wave amplitude.\n " );
      }
      else if (( opm1 == null )
             &&( idx[null] <= MDV_SPCLPRTS ))
      {
         opm2 = null;
         printf( "\n\n normalizing S-parameters relative to port1.\n " ); 
      };
         
/* = number of vswr evaluation points: */

      kk = idx[null] - MDV_SPCLPRTS;

/* number of spline interpolataion points between support points */
/* spl.vct[][null]: */

      if ( null < kk )
      {
         i_ = MDV_GPHINTPL;
   
         if ( GPH_POINTS <= kk * ( i_ + ONE ) ) 
            i_ = ( short ) (( double)( GPH_POINTS / kk ) - ONE );
      };

      j_ = null; 
      while ( j_ < frqlbl )
      {
         if ( null < kk ) /* i.e. if ( MDV_SPCLPRTS < idx[null] ) */ 
         {
            ii = null;
            hh = ONE;
            h_ = null;
            while( hh < kk )
            {
               mm = hh + MDV_SPCLPRTS;
               nn = hh + MDV_SPCLPRTS + ONE;
               spl.vct[h_][null] = idx[mm];
               spl.vct[h_][ONE]  = aspc[mm-ONE][j_];
               spl.dmn[ii] = idx[mm];
               ii++;
               for ( k_= ONE; k_<= i_ ; k_++ )
               {
                  x2 = ( double ) k_ / ( i_ + ONE ) ;
                  x1 = 1. - x2;
                  spl.dmn[ii]  = x1*idx[mm] + x2*idx[nn];
                  ii++;
               };
               h_ = hh;
               hh++ ;
            };/* next hh */
            spl.vct[h_][null] = idx[idx[null]];
            spl.vct[h_][ONE]  = aspc[idx[null]-ONE][j_];
            spl.dmn[ii] = idx[idx[null]];
            ii++ ;

            spl.nn = ii;
            spl.mm = kk;
/*............................................................................*/
            spt = spline( spt );        /* spline interpolation: stand. wave  */
/*....................................*//* over pnts.idx, at frequency frq[j_]*/
            if (( spt->rtn ) == ONE )
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling spline function !!!" );
               printf( "\n [ program stopped.]\n" );
	       
               exit( EXIT_FAILURE );
            };

            minm = +1.e+277;
            maxm = - minm;

            jj = null;
            while( jj < ii )
            {
               if ( fabs( spl.fct[jj] ) < minm )
                  minm = fabs( spl.fct[jj] );
               if ( maxm < fabs( spl.fct[jj] ) )
                  maxm = fabs( spl.fct[jj] );
               jj++ ;
            };

            mean[j_] = .5*( minm + maxm ); /* mean incident wave amplitude */
/*............................................................................*/
# if MDV_PLOT_SWR == 1

            if ( j_ == null )
            {
               jj = null;
               while ( jj < ii )
               {
                  gph.vct[jj][null] = spl.dmn[jj];
                  gph.vct[jj][ONE]  = spl.fct[jj];
                  jj++ ;
               };
               strcpy( gph.file, "stw" );
               strcat( gph.file, lotos( j_, null ));
               strcat( gph.file, state->flbl );
               strcpy( gph.name, state->name );
               strcpy( gph.text, "standing_wave" );
               strcpy( gph.xunit, "index" );
               strcpy( gph.yunit, "volts" );
               gph.nn = ii;
/*............................................................................*/
               ind = graphp( gpt );               /* graphics file: std.wave  */
/*..............................................*//* at frequency frq[j_]     */
            };
# endif /* MDV_PLOT_SWR == 1 */
/*............................................................................*/
         }
         else /* if ( idx[null] <= MDV_SPCLPRTS ) *//* normalization to port1 */
            mean[j_] = aspc[null][j_];

         if ( opm1 == null ) /* absent reference spectrum:  */
                             /* normalize to incident wave  */
                             /* [ mean absolute amplitude ] */
         {
            real = mean[j_];
            imag = ZERO;
            norm = real*real + imag*imag; 

            if ( null < kk ) /* i.e standing wave ports to be evaluated       */
            {                
               if ( 1.e-277 < minm )
               {
                  yy = maxm/minm;
                  vswr[j_] = yy;
                  zz = ( yy - 1. )/( yy + 1. );
                  refl[j_] = 100.*zz;
               }
               else
               {
                   vswr[j_] = HUGE_VALF;
                   refl[j_] = 100.; 
               };
/*  
               rspc[null][j_] -= real;             ???  - inspection !!!    
               ispc[null][j_] -= imag;             
*/
            };
         }
         else if ( opm1 == ONE )  /* normalization to reference spectrum      */
         {
            real = rref[j_];
            imag = iref[j_];
            norm = real*real + imag*imag;

            rspc[null][j_] -= real;
            ispc[null][j_] -= imag;
            xx = ( real*rspc[null][j_] + imag*ispc[null][j_] ) / norm;
            yy = ( real*ispc[null][j_] - imag*rspc[null][j_] ) / norm;
            zz = sqrt( xx*xx + yy*yy );

            if ( 1.e-277 < fabs( 1. - zz ) )
            {
               vswr[j_] = ( 1. + zz )/( 1. - zz );
               refl[j_] = 100.*zz;
            } 
            else
            {
               vswr[j_] = HUGE_VALF;
               refl[j_] = 100.;
            };
         };

         ii = MDV_SPCLPRTS;

         if ( idx[null] < ii )
            ii = idx[null]; 

         h_= null;
         while( h_< ii )
         {
            xx = ( real*rspc[h_][j_] + imag*ispc[h_][j_] ) / norm;
            yy = ( real*ispc[h_][j_] - imag*rspc[h_][j_] ) / norm;

            rcpl[h_][j_] = xx;
            icpl[h_][j_] = yy;

            zz = xx*xx + yy*yy; 

            if ( 1.e-277 < zz )
            {
/*............................................................................*/
# if MDV_NORMALZE == 1 /* normalize S-parameters to abs.value one */

               zz = sqrt( zz ); 
               rcpl[h_][j_] /= zz;
               icpl[h_][j_] /= zz;

               acpl[h_][j_] = 0.;
# else
               acpl[h_][j_] = 10.*log10( zz ); 
# endif
/*............................................................................*/
            }
            else
               acpl[h_][j_] = - LARGE_LOG_VAL;
            h_++ ;
         };
         j_++ ;
      };/* while ( j_ < frqlbl ) */
      
/* s-parameters ready */
/*............................................................................*/
/* save S-parameters: */

      ii = MDV_SPCLPRTS;

      if ( idx[null] < ii )
         ii = idx[null];
       
      h_ = null;
      while( h_ < ii )
      {
         strcpy( ptr, "s" );
         strcat( ptr, lotos ( h_+ONE, null ) );
         strcat( ptr, "1_" );
         strcat( ptr, state->flbl );
         
         strcpy( spcfle, spcpfx );
         strcat( spcfle, ptr );

         evalfle = fopen( spcfle, "w+" );

         jj = strlen( ptr );
         strcat( ptr, "_<<" );
         strncat( ptr, vpt->text, ( STS_SIZE - SIX - jj ));
         strcat( ptr, ">>" );
              
         fprintf( evalfle, spformat, ( vpt->name ));
         fprintf( evalfle, spformat, ptr );
         fprintf( evalfle, spformat, absc_unit );
         fprintf( evalfle, spformat, "---" );
         fprintf( evalfle, dpformat, xlower, "\n" );
         fprintf( evalfle, dpformat, xupper, "\n" );
         fprintf( evalfle, dpformat, dx, "\n" );
         fprintf( evalfle, ipformat, frqlbl );
/*............................................................................*/
# if EVL_WRTFREQ == 1
	 frequency = xlower;
# endif
/*............................................................................*/
         j_ = null;
         while ( j_ < frqlbl )
         {
/*............................................................................*/
# if EVL_WRTFREQ == 1
            fprintf( evalfle, dpformat, frequency, "  " );
	    frequency += dx;
# endif
/*............................................................................*/
            fprintf( evalfle, dpformat, rcpl[h_][j_], "  " );
            fprintf( evalfle, dpformat, icpl[h_][j_], "\n" );
            j_++ ;
         };

         nseconds = time(timer);
         strcpy( tmeptr, ctime( &nseconds ));

         fprintf( evalfle, "\n# SPECTRUM file %s\n", spcfle );
         fprintf( evalfle, "# created:%s", tmeptr );

         fclose( evalfle );

         printf( "\n" );
         printf( CLEAR_LINE );
         printf( "\r SPECTR.file %s ", spcfle );
         printf( timefrm, tmeptr );

         h_++ ;
      };
/*............................................................................*/
      i_ = MDV_GPHINTPL;         /* number of interpoltation points           */
                                 /* between frequencies frq                   */
      if ( GPH_POINTS <= frqlbl * ( i_ + ONE ) )
             i_ = ( short )(( double )( GPH_POINTS / frqlbl ) - ONE );

      ii = null;
      j_ = null;
      while( j_< frqlbl-ONE )
      {
         spl.vct[j_][null] = frq[j_];

         if ( opm2 == ONE )
            spl.vct[j_][ONE] = vswr[j_];

         spl.dmn[ii] = frq[j_];
         ii++;
         for ( k_= ONE ; k_<= i_ ; k_++ )
         {
            x2 = ( double ) k_ / ( i_ + ONE );
            x1 = 1. - x2;
            spl.dmn[ii]  = x1*frq[j_] + x2*frq[j_+ONE];
            ii++ ;
         };
         j_++ ;
      }; /* until j_ = frqlbl - ONE */
      spl.vct[j_][null] = frq[j_];

      if ( opm2 == ONE )
         spl.vct[j_][ONE] = vswr[j_];

      spl.dmn[ii] = frq[j_];
      ii++;

      if ( opm2 == ONE )
      {
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation: VSWR            */
/*.................................*/
         ii = null;
         j_ = null;
         while( j_ < frqlbl-ONE )
         {
            spl.vct[j_][ONE] = refl[j_];
            gph.vct[ii][null] = spl.dmn[ii];
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][null] = spl.dmn[ii];
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++;
         }; /* until j_ = frqlbl - ONE */
         spl.vct[j_][ONE] = refl[j_];
         gph.vct[ii][null] = spl.dmn[ii];
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++ ;

         strcpy( gph.file, "swr" );
         strcat( gph.file, state->flbl );
         strcpy( gph.name, state->name );
         strcpy( gph.text, "voltage_standing_wave_ratio " );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "---" );
         gph.nn = ii;

         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );              /* graphics file: VSWR             */
                                          /*                                  */
         spt = spline( spt );            /* spline interpolation: reflexion   */
/*.....................................*/
         ii = null;
         j_ = null; 
         while( j_ < frqlbl-ONE )
         {
            zz = fabs( refl[j_]/100. );

            if ( 0. < zz )
               spl.vct[j_][ONE] = 20.*log10( zz );
            else                                  /* Replacing LARGE_LOG_VAL  */
               spl.vct[j_][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                                  /* made spline(*) function  */
                                                  /* inoperative              */
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++;
         }; /* until j_ = frqlbl - ONE */

         zz = fabs( refl[j_]/100. );

         if ( 0. < zz )
            spl.vct[j_][ONE] = 20.*log10( zz );
         else                                  /* Replacing LARGE_LOG_VAL  */
            spl.vct[j_][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                               /* made spline(*) function  */
                                               /* inoperative              */
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++ ;

         strcpy( gph.file, "ref" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "reflexion " );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "0/0" );

         gph.nn = ii;
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );              /* graphics file: reflex.[%], port1*/
                                          /*                                  */
         spt = spline( spt );            /* spline interpolation: refl.[dB],p1*/
/*.....................................*/
         ii = null;
         j_ = null;
         while( j_ < frqlbl-ONE )
         {
            zz = fabs( refl[j_]/100. );

            if( zz < 1. )
               spl.vct[j_][ONE] = 10.*log10( 1. - zz*zz );
            else                                  /* Replacing LARGE_LOG_VAL  */
               spl.vct[j_][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                                  /* made spline(*) function  */
                                                  /* inoperative              */
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++;
         }; /* until j_ = frqlbl - ONE */

         zz = fabs( refl[j_]/100. );

         if ( zz < 1. )
            spl.vct[j_][ONE] = 10.*log10( 1. - zz*zz );
         else                                  /* Replacing LARGE_LOG_VAL  */
            spl.vct[j_][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                               /* makes spline(*) function */
                                               /* inoperative              */
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++ ;

         strcpy( gph.file, "rtl" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "return_loss" );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "dB" );

         gph.nn = ii;
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );              /* graphics file: return loss [dB] */
                                          /*                                  */
         spt = spline( spt );            /* spline intpl.: insert.loss [dB],p1*/
/*.....................................*/
      }; /* end if ( opm2 == ONE ) */

      ii = null;
      j_ = null; 
      while( j_ < frqlbl-ONE )
      {
         if ( null < MDV_SPCLPRTS )
            spl.vct[j_][ONE] = acpl[null][j_];

         if ( opm2 == ONE ) 
         {
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ; 
            };
         };
         j_++ ;
      }; /* until j_ = frqlbl - ONE */

      if ( null < MDV_SPCLPRTS )
         spl.vct[j_][ONE] = acpl[null][j_];

      if ( opm2 == ONE ) 
      {
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++ ;

         strcpy( gph.file, "isl" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "insertion_loss " );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "dB" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );              /* graphics file: insert.loss [dB] */
/*.......................................*/
      }; /* end if ( opm2 == ONE ) */
/* reflexion, VSWR,..., insertion loss etc. have been saved */
/*............................................................................*/
/* plot s-parameters, special ports: [ in PLOT_FORMAT, e.g. "SPLINE" ] */

      ll = null;
      while (( ll < MDV_SPCLPRTS )
           &&( ll < idx[null] ))
      {
         strcpy( ptr, lotos( ll+ONE, null ));

         spl.mm = frqlbl;
         spl.nn = ii;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation: |sn1|           */
/*.................................*/ 
         ii = null;
         j_ = null;
         while ( j_ < frqlbl-ONE )
         {
            spl.vct[j_][ONE] = rcpl[ll][j_];
            gph.vct[ii][null] = spl.dmn[ii];
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++ ;
         }; /* until j_ = frqlbl-ONE */

         spl.vct[j_][ONE] = rcpl[ll][j_];
         gph.vct[ii][null] = spl.dmn[ii];
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++;

         strcpy( gph.file, "s" );
         strcat( gph.file, ptr ); 
         strcat( gph.file, "1_" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "|S" );
         strcat( gph.text, ptr );
         strcat( gph.text, "1|" );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "dB" );

         gph.nn = ii;
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );                /* graphics file: |sn1|          */
                                            /*                                */
         spt = spline( spt );              /* spline interpolation: real(sn1) */
/*.......................................*/
         ii = null;
         j_ = null;
         while ( j_ < frqlbl-ONE )
         {
            spl.vct[j_][ONE] = icpl[ll][j_];
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_ ++;
         }; /* until j_ = frqlbl-ONE */
         spl.vct[j_][ONE] = icpl[ll][j_];
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++;

         strcpy( gph.file, "r" );
         strcat( gph.file, ptr );
         strcat( gph.file, "1_" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "real(s" );
         strcat( gph.text, ptr );
         strcat( gph.text, "1)" );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "--" );

         gph.nn = ii;
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );                /* graphics file: real(sn1)      */
                                            /*                                */
         spt = spline( spt );              /* spline interpolation: imag(sn1) */
/*.......................................*/
         ll++;

         ii = null;
         j_ = null;
         while ( j_ < frqlbl-ONE )
         {
            if (( ll < MDV_SPCLPRTS  )
              &&( ll < idx[null] ))
               spl.vct[j_][ONE] = acpl[ll][j_];

            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++ ;
         }; /* until j_ = frqlbl-ONE */

         if (( ll < MDV_SPCLPRTS  )
           &&( ll < idx[null] ))
            spl.vct[j_][ONE] = acpl[ll][j_];

         gph.vct[ii][ONE] = spl.fct[ii];
         ii++;

         strcpy( gph.file, "i" );
         strcat( gph.file, ptr );
         strcat( gph.file, "1_" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "imag(s" );
         strcat( gph.text, ptr );
         strcat( gph.text, "1)" );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "--" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );                /* graphics file: imag(sn1)      */
/*.........................................*/

      };/* end while ( ll < SPECIAL_ ... ) */
      break;
/*.................... s-parameter computation finished ......................*/

     case 4:

      hh = ONE;
      while ( hh <= idx[null] )
      {
         h_ = hh - ONE;

         evl1 = ZERO;
         dlt1 = ZERO;

         absmaxm[h_] = ZERO;
         absmean[h_] = ZERO;
         meanenv[h_] = ZERO;
         timemax[h_] = ZERO;
            envl[h_] = null;

         i_ = null;
         while ( i_ < ( fpt->ttlg[null] ))
         {
            evl0 = ( fpt->r[hh][i_] );
            norm = fabs( evl0 );
            dlt0 = evl0 - evl1;

            absmean[h_] += norm;

            if ( absmaxm[h_] < norm )
            {
               absmaxm[h_] = norm;
               timemax[h_] = ( fpt->t[null] ) + i_*( fpt->dt[null] );
            };

            if ( dlt0*dlt1 < ZERO )/* rel.extremum: compute mean.envlp.*/
            {
               envl[h_]++;
               meanenv[h_] += fabs( evl1 );
            };

            if ( null < i_ )
               dlt1 = dlt0;

            evl1 = evl0;

            i_++ ;
         }; /* next i_ */

         absmean[h_] /= ( fpt->ttlg[null] );

         if ( null < envl[h_] )
         {
            meanenv[h_] /= envl[h_];
            modulation[h_] = 100.*(fabs( absmaxm[h_] - meanenv[h_] ));

            if ( 1.e-277 < fabs(meanenv[h_] ) )
               modulation[h_] /= meanenv[h_];
            else
               modulation[h_] = HUGE_VALF;
         }
         else
         {
            meanenv[h_] = ZERO;
            modulation[h_] = ZERO;
         };
         hh++ ;
      }; /* next hh */

/* tabl4: */

      printf( "\n\n Special parameters evaluated in file %s :",
         ( state->file ));
      printf( "\n\n index | abs. maximum : found at time|    abs. mean|  "
         "mean envlp.|   modulation" );
      printf( "\n       | [%10s] : [ %10s]|    [%7s]|  [%9s]| %12s",
         ( vpt->yunit ), ( vpt->xunit ), ( vpt->yunit ), ( vpt->yunit ),
         "[ percent]" );
      printf( "\n -----------------------------------"
         "-------------------------------------------" );

      hh = ONE;
      while ( hh <= idx[null] )
      {
         h_ = hh - ONE;
         printf( "\n %6d| %.7e: %.7e|%.7e|%.7e|%.7e", idx[hh], absmaxm[h_],
         timemax[h_], absmean[h_], meanenv[h_], modulation[h_] );
         hh++ ;
      };
      printf( "\n -----------------------------------"
         "-------------------------------------------" );

/* menu4: */

      kk = idx[null] - MDV_SPCLPRTS;

      i_ = MDV_GPHINTPL; /* number of interpolation points */

      ii = null;
      h_ = null;
      hh = ONE;
      while ( hh < kk )
      {
         mm = hh + MDV_SPCLPRTS;
         nn = hh + MDV_SPCLPRTS + ONE;

         spl.vct[h_][null] = idx[mm];
         spl.vct[h_][ONE]  = meanenv[mm-ONE];
         spl.dmn[ii] = spl.vct[h_][null];

         k_ = ONE;
         while ( k_ <= i_ )
         {
            x2 = ( double ) k_ / ( i_ + ONE );
            x1 = 1. - x2;
            spl.dmn[ii]  = x1*idx[mm] + x2*idx[nn];
            ii++ ;
            k_++ ;
         };
         h_ = hh;
         hh++ ;
      }; /* until hh = kk; h_ = kk - ONE */
      spl.vct[h_][null] = idx[idx[null]];
      spl.vct[h_][ONE]  = meanenv[idx[null]-ONE];
      spl.dmn[ii] = idx[idx[null]];
      ii++ ;

      printf( "\n %s: VSWR spline interpolation started.", __func__ );

      spl.mm = kk;
      spl.nn = ii;
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: VSWR                 */
/*............................*/

      printf( "\r %s: VSWR spline interpolation terminated.", __func__ );

      spl.fmax = fabs( spl.fmax );
      spl.fmin = fabs( spl.fmin );

      yy = spl.fmax/spl.fmin;               /* voltage standing wave ratio    */
      vswr[null] = yy;

      xx =  ( yy - 1. )/( yy + 1. );        /* reflexion factor               */
      refl[null] = xx;

      if ( fabs( xx ) < .75 )
         printf( "\n\n Source line matching:" );
      else
         printf( "\n\n Source line matching [ approximate values ]:" );

      printf( "\n\n VSWRatio  = % .12e ", yy );

      zz = 20.*log10( xx );                 /* reflexion [dB]                 */
      printf( "\n\n REFLEXION = % .12e dB [ %.12e percent ] ", zz, 100.*xx );

      yy = sqrt( 1. - xx*xx );              /* insertion factor               */
      zz = 20.*log10( yy );                 /* insertion loss [dB]            */
      printf( "\n INSERTION = % .12e dB [ %.12e percent ] ", zz, 100.*yy );

      printf( "\n\n S-parameters: \n" );

      zz = 20.*log10( xx );                 /* S11 = reflexion [dB]           */
      printf( "\n S[1,1] = % .12e dB ", zz );

      xx = .5*( spl.fmax + spl.fmin );      /* incident wave amplitude        */
      xx *= xx;                              

      hh = ONE;
      while ( hh < MDV_SPCLPRTS )
      {
         yy = meanenv[ hh ];
         zz = 10.*log10( yy*yy / xx );
         printf( "\n S[%ld,1] = % .12e dB ", hh+ONE, zz );
         hh++ ;
      };

      printf( "\n\n Please acknowledge [ enter any character ]: " );
      scanf( "%s", ptr );

      strcpy( gph.file, "stw" );
      strcat( gph.file, ( state->flbl ));
      strcpy( gph.name, ( state->name ));
      strcpy( gph.text, "standing_wave" );
      strcpy( gph.xunit, "index" );
      strcpy( gph.yunit, "volts" );

      ii = null;
      hh = ONE;
      while ( hh < kk )
      {
         gph.vct[ii][null] = spl.dmn[ii];
         gph.vct[ii][ONE]  = spl.fct[ii];

         k_ = ONE;
         while ( k_ <= i_ )
         {
            gph.vct[ii][null] = spl.dmn[ii];
            gph.vct[ii][ONE]  = spl.fct[ii];
            ii++ ;
            k_++ ;
         };
         hh++ ;
      }; /* next hh */
      gph.vct[ii][null] = spl.dmn[ii];
      gph.vct[ii][ONE]  = spl.fct[ii];
      ii++;

      gph.nn = ii;
/*............................................................................*/
      ind = graphp( gpt );      /* graphics file: standing wave               */
/*............................*/
      break;
/*............................................................................*/

     case 5:

      if ( null == strncmp( vpt->text, "reference", NINE ))
      {
         strcpy( spcfle, spcpfx );
         strcat( spcfle, "ref" );
         strcat( spcfle, state->flbl );

         evalfle = fopen( spcfle, "w+" );

         if ( evalfle == null )
         {
            printf( "\n Error on opening reference file "
               "%s !", spcfle );
            printf( "\n [overriding: values not saved.]\n " );
            break;
         };

         fprintf( evalfle, "%s\n", ( vpt->name ));
         fprintf( evalfle, "%s\n", "reference_spectrum" );
         fprintf( evalfle, "%s\n", ( vpt->xunit ));
         fprintf( evalfle, "%s\n", ( vpt->yunit ));
	 
         fprintf( evalfle, dpformat, ( double )( fpt->t[null] ), "\n" );
         fprintf( evalfle, dpformat, ( double )( fpt->tt[null] ), "\n" );
         fprintf( evalfle, dpformat, ( double )( fpt->dt[null] ), "\n" );
         fprintf( evalfle, "%ld \n", ( fpt->ttlg[null] ));
/*............................................................................*/
# if EVL_WRTFREQ == 1
	 frequency = ( fpt->t[null] );
# endif
/*............................................................................*/
         jj = null;
         while ( jj < ( fpt->ttlg[null] ))
         {
/*............................................................................*/
# if EVL_WRTFREQ == 1
            fprintf( evalfle, dpformat, frequency, "  " );
	    frequency += ( fpt->dt[null] );
# endif
/*............................................................................*/
            fprintf( evalfle, dpformat, ( fpt->r[ONE][jj] ), "  " );
            fprintf( evalfle, dpformat, ( fpt->i[ONE][jj] ), "\n" );
            jj++ ;
         };

         nseconds = time( timer );
         strcpy( tmeptr, ctime( &nseconds ));

         fprintf( evalfle, "\n# DSC reference %s\n", spcfle );
         fprintf( evalfle, "# created:%s", tmeptr );

         fclose( evalfle );

	 printf( "\n" );
         printf( CLEAR_LINE );
         printf( "\n Reference file %s ", spcfle );
         printf( "created:\n %.24s\n ", tmeptr );

      } /* end if ( null == strncmp( vpt->text, "reference", NINE )) */
      else if ( null < idx[null]  )
      {
/* enter external reference <spcpfx>.ref: */

         strcpy( spcfle, spcpfx );
         strcat( spcfle, "ref" );

         printf( "\n" );

        open_reference2:

         evalfle = fopen( spcfle, "r+" );

         if ( evalfle == null )
         {
            printf( "\n Reference file %s not found "
               "in present directory:", spcfle );
            printf( "\n Please re-enter filename [ Escape: "
               "enter null ] >----> " );
            scanf( "%s", spcfle );

            if ( *spcfle == '0' )
               break;
            else
               goto open_reference2;
         }
         else if ( evalfle != null )
         {
            printf( "\n opened: reference file %s ", spcfle );

            fscanf( evalfle, scformat, ptr );

            if ( null != strncmp( ptr, vpt->name, THREE ))
            {
               printf( "\n\n Error message from function %s :", __func__ );
               printf( "\n\n Incompatible system identifier '%s'", ptr );
               printf( "\n on reference spectrum, file %s !!!", spcfle );
               printf( " [ overriding. ]\n" );

               fclose( evalfle );

               break;
            };

            fscanf( evalfle, scformat, ptr );
            fscanf( evalfle, scformat, absc_unit );
            fscanf( evalfle, scformat, ordn_unit );
            fscanf( evalfle, scformat, ptr );/* xlower */
            xlower = strtod( ptr, endp );

            fscanf( evalfle, scformat, ptr );/* xupper */
            xupper = strtod( ptr, endp );

            fscanf( evalfle, scformat, ptr ); /* dx */
            dx = strtod( ptr, endp );

            fscanf( evalfle, scformat, ptr );
            ii = strtol( ptr, endp, DEC );

            xx = xlower;
            jj = null;
            while( jj < ii )
            {
/*............................................................................*/
# if EVL_WRTFREQ == 1
               fscanf( evalfle, scformat, ptr );
# endif
/*............................................................................*/
               fscanf( evalfle, scformat, ptr );
               rref[jj] = strtod( ptr, endp );
               fscanf( evalfle, scformat, ptr );
               iref[jj] = strtod( ptr, endp );
               xx += dx;

               jj++ ;
            };
            xupper = xx - dx;
            printf( "\r entered: reference %s .      \n",
               spcfle );
            fclose( evalfle );
            opm1 = ONE;
         };

         jj = null; 
         while( jj < ii )
         {
            real = rref[jj];
            imag = iref[jj];
            norm = real*real + imag*imag;

            ( fpt->r[ONE][jj] ) -= real;
            ( fpt->i[ONE][jj] ) -= imag;

            rspc[null][jj] = ( fpt->r[ONE][jj] );
            ispc[null][jj] = ( fpt->i[ONE][jj] );

            if ( 1.e-277 < norm )
            {
               xx = ( real*rspc[null][jj] + imag*ispc[null][jj] ) / norm;
               yy = ( real*ispc[null][jj] - imag*rspc[null][jj] ) / norm;
               zz = sqrt( xx*xx + yy*yy );

               if ( 1.e-277 < fabs( 1. - zz ) )
               {
                  vswr[jj] = ( 1. + zz )/( 1. - zz );
                  refl[jj] = 100.*zz;
               } 
               else
               {
                  vswr[jj] = HUGE_VALF;
                  refl[jj] = 100.;
               };
            }
            else
            {
               xx = HUGE_VALF;
               yy = HUGE_VALF;
               zz = HUGE_VALF;
               vswr[jj] = HUGE_VALF;
               refl[jj] = HUGE_VALF;
            };

            kk = MDV_SPCLPRTS;

            if ( idx[null] < kk )
               kk = idx[null]; 

            hh = ONE;
            while( hh <= kk )
            {
               h_ = hh - ONE;

               real = rref[jj];
               imag = iref[jj];
               norm = real*real + imag*imag;

               rspc[h_][jj] = ( fpt->r[hh][jj] );
               ispc[h_][jj] = ( fpt->i[hh][jj] );

               if ( 1.e-277 < norm )
               {
                  xx = ( real*rspc[h_][jj] + imag*ispc[h_][jj] ) / norm;
                  yy = ( real*ispc[h_][jj] - imag*rspc[h_][jj] ) / norm;
                  zz = xx*xx + yy*yy;

                  if ( 1.e-277 < zz )
                  {
/*............................................................................*/
# if MDV_NORMALZE == 1

                     zz = sqrt( zz );
                     xx /= zz;
                     yy /= zz;
                  
                     acpl[h_][jj] = 0.;
# else
                     acpl[h_][jj] = 10.*log10( zz );
# endif
/*............................................................................*/
                  }
                  else
                     acpl[h_][jj] = - LARGE_LOG_VAL;
               }
               else
               {
                  xx = HUGE_VALF;
                  yy = HUGE_VALF;
                  zz = HUGE_VALF; 

                  acpl[h_][jj] = LARGE_LOG_VAL;
               };

               rcpl[h_][jj] = xx;
               icpl[h_][jj] = yy;

               hh++ ;
            };
            jj++ ;

         };/* next jj */

         jj = null; xx = xlower; 
         while( jj < ii )
         {
            gph.vct[jj][ONE] = vswr[jj];
            gph.vct[jj][null] = xx;

            xx += dx;
            jj++ ;
         };/* next jj */

         strcpy( gph.file, "swr" );
         strcat( gph.file, state->flbl );
         strcpy( gph.name, state->name );
         strcpy( gph.text, "voltage_standing_wave_ratio" );
         strcpy( gph.xunit, absc_unit );
         strcpy( gph.yunit, "---" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );        /* graphics file: VSWR                   */
/*.................................*/
         jj = null; xx = xlower;
         while( jj < ii )
         {
            gph.vct[jj][ONE] = refl[jj];
            gph.vct[jj][null] = xx;

            xx += dx ;
            jj++ ;
         };/* next jj */

         strcpy( gph.file, "ref" );
         strcat( gph.file, state->flbl );
         strcpy( gph.name, state->name );
         strcpy( gph.text, "reflexion" );
         strcpy( gph.xunit, absc_unit );
         strcpy( gph.yunit, "0/0" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );        /* graphics file: reflexion [0/0]        */
/*.................................*/
         jj = null; xx = xlower;
         while( jj < ii )
         {
            zz = refl[jj];

            if( 1.e-277 < zz  )
               gph.vct[jj][ONE] = 20.*log10( zz/100. );
            else                                  /* Replacing LARGE_LOG_VAL  */
               gph.vct[jj][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                                  /* made graphp(*) function  */
                                                  /* inoperative              */
            gph.vct[jj][null] = xx;

            xx += dx;
            jj++ ;
         };/* next jj */

         strcpy( gph.file, "rtl" );
         strcat( gph.file, ( state->flbl ));
         strcpy( gph.text, "return_loss" );
         strcpy( gph.xunit, absc_unit );
         strcpy( gph.yunit, "dB" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );        /* graphics file: return loss [dB]       */
/*.................................*/
         jj = null; xx = xlower;
         while( jj < ii )
         {
            zz = refl[jj]/100.;
       
            if( zz < 1. )
               gph.vct[jj][ONE] = 10.*log10( 1. - zz*zz );
            else                                  /* Replacing LARGE_LOG_VAL  */
               gph.vct[jj][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                                  /* made graphp(*) function  */
                                                  /* inoperative              */
            gph.vct[jj][null] = xx;

            xx += dx;
            jj++ ;
         };/* next jj */

         strcpy( gph.file, "isl" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "insertion_loss" );
         strcpy( gph.xunit, absc_unit );
         strcpy( gph.yunit, "dB" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );        /* graphics file: insert.loss [dB]       */
/*.................................*/
/* still as above:
         kk = MDV_SPCLPRTS;
	 
         if ( idx[null] < kk )
            kk = idx[null]; 
*/
         hh = ONE;
         while( hh <= kk )
         {
            h_ = hh - ONE;
            jj = null; xx = xlower;
            while( jj < ii )
            {
               gph.vct[jj][ONE] = acpl[h_][jj];
               gph.vct[jj][null] = xx;

               xx += dx;
               jj++ ;
            };/* next jj */
            
            strcpy( ptr, "s" );
            strcat( ptr, lotos( hh, null ) );
            strcat( ptr , "1_" );
            strcat( ptr, state->flbl );
            strcpy( gph.file, ptr );
            strcpy( gph.text, ptr );
/*
            strcat( gph.text, "[steady_state]" );
*/
            strcat( gph.text, "_" );
            strcat( gph.text, vpt->text );
            strcpy( gph.xunit, absc_unit );
            strcpy( gph.yunit, "dB" );
            gph.nn = ii;
/*............................................................................*/
            ind = graphp( gpt );           /* graphics file: s[k,1] [dB]      */
/*.......................................*/
/* complex s-parameters : */

            strcpy( spcfle, spcpfx );
            strcat( spcfle, ptr );

            evalfle = fopen( spcfle, "w+" );

            if ( evalfle == null )
            {
               printf( "\n Error on opening file "
                  "%s !", spcfle );
               printf( "\n [overriding: values not saved.]\n " );
            }
            else
            {
               fprintf( evalfle, "%s\n", ( vpt->name ));
               fprintf( evalfle, "%s_%s\n", ptr, ( vpt->text ));
               fprintf( evalfle, "%s\n", ( vpt->xunit ));
               fprintf( evalfle, "%s\n", "---" );
               fprintf( evalfle, dpformat, xlower, "\n" );
               fprintf( evalfle, dpformat, xupper, "\n" );
               fprintf( evalfle, dpformat, dx, "\n" );
               fprintf( evalfle, "%ld\n", ii );
/*............................................................................*/
# if EVL_WRTFREQ == 1
	       frequency = xlower;
# endif
/*............................................................................*/
               jj = null;
               while ( jj < ii )
               {
/*............................................................................*/
# if EVL_WRTFREQ == 1
                  fprintf( evalfle, dpformat, frequency, "  " );
	          frequency += dx;
# endif
/*............................................................................*/
                  fprintf( evalfle, dpformat, rcpl[h_][jj], "  " );
                  fprintf( evalfle, dpformat, icpl[h_][jj], "\n" );
                  jj++ ;
               };

               nseconds = time( timer );
               strcpy( tmeptr, ctime( &nseconds ));

               fprintf( evalfle, "\n# DSC file %s\n", spcfle );
               fprintf( evalfle, "# created:%s", tmeptr );

               fclose( evalfle );

            }; /* end if evalfle != null */  

            hh++ ;
         }; /* end while ( hh <= kk [ kk <= MDV_SPCLPRTS ) ] */
      }; /* end if ( null < idx[null]  ) */

      break;
/*............................................................................*/

     default:
      break;
   }; /* end switch ( lbl ) */

   PRBLDCLR( "\r" );
   printf( " %*s", 78, "MODVAL" );
   PRNORMAL( "" );
   ( state->rtn ) = null;

   return state;
}
/*=================== end of function body modval(*) =========================*/
# undef LARGE_LOG_VAL
/************************ end of function modval(*) ***************************/


/*----------------------------------------------------------------------------*/
# endif /* [ end of section compiled with option -D_POST_MODEL ] */
/*----------------------------------------------------------------------------*/
# undef EXC_MAXW
# undef EXC_HCRR
/*************************** end of file model.c ******************************/
