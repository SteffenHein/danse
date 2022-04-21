/* [ file: postdrv.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function postdrv(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   This function evaluates the output files of type dsc.val<N> wherein        *
*   values computed with the DSC solver SOLVER.C are stored                    *
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
# ifdef OPTIMIZE
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 3
# endif
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize this header for POSTER.C configuration: */
# include "./POSTER.CONF"
/*----------------------------------------------------------------------------*/
/* A macro that gives a conveniently readable form to the string returned */
/* with [ the pointer of ] function ctime(&nseconds) */
# include "../tools/TIMEFORM.M"
/*----------------------------------------------------------------------------*/
# if USE_NCURSES
   # include <ncurses.h>
# endif 
/*----------------------------------------------------------------------------*/
# ifndef VAL_INITLBL
   # define VAL_INITLBL 1 /* if val.ni < VAL_INITLBL the val.ni = VAL_INITLBL */
# endif
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
   # ifndef DSC_FCTEMP
      # define DSC_FCTEMP 1 /* 1: read face temperatures */
   # endif

   # ifndef DSC_INTLCE
      # define DSC_INTLCE 1 /* 0/1: separate/interlaced internal Maxwell */
   # endif                   /* field and heat current computation loops */
# endif
/*----------------------------------------------------------------------------*/
/* surface macros [ e.g. display instructions without influence on internal   */
/* algorithms ]:                                                              */

# define EVL_ITEMS 7 /* number of items in Eval.menu [ except  0 = end ] */
# define EVL_DISP  1 /* additional monitoring of intermediate results */
/*----------------------------------------------------------------------------*/
# ifndef EVL_GPHINTP
   # define EVL_GPHINTP 3
# endif
/*----------------------------------------------------------------------------*/
/* file macros  [ names, format instructions etc. ]:                          */

# define DSC_PRFX        "dsc." /* common file prefix for DSC model files */
# define EVALUATION_FILE "val"  /* evaluation file identifier [ with index ] */
# define RESPONSE_FILE  "time"
# define SPECTRUM_FILE  "spec"
/*----------------------------------------------------------------------------*/
/* data alignment modes [ depending on Fast Fourier transform conventions
   followed - usually defined in POSTER.CONF ] 
   READ_REVERSE_ : read distributions starting with upper half of domain
   WRITE_REVERSE : save spectra starting with upper half of domain           */

# ifndef READ_REVERSE_
   # define READ_REVERSE_ 0
# endif
# ifndef WRITE_REVERSE
   # define WRITE_REVERSE 1
# endif
/*----------------------------------------------------------------------------*/
# include "../poster/posttp.h" /* typedefs: EVALUATE etc. */
# include "../math/gssjtp.h"
# include "../tools/txctyp.h"
/*----------------------------------------------------------------------------*/
static EVALUATE evl = {null};
static GAUSS_JRD gss = {null};
static GRAPHICS gph = {null};
static SPLINES spl = {null};
static FFT fft = {null};
static POSTSTATE eval = {ZERO};
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* 'my_terminal' configuration: */

   # include <termcap.h> /* terminal type header */
   static char *term;    /* terminal type string */ 

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
# include "readval.h"
/*============================================================================*/

short postdrv( char *err )
{
/* allusions: */
/*
   extern POSTSTATE eval;
   extern GAUSS_JRD gss;
   extern GRAPHICS gph;
   extern SPLINES spl;
   extern FFT fft;
   extern EVALUATE evl;
*/
/* declarations: */

   static FILE
     *timefle[MAX_PERIOD],
     *specfle[MAX_PERIOD],
     *evalfle = NULL,
     *pltfle = NULL,
     *datfle = NULL;

   static char 
      ptr[STS_SIZE] = {null},             
      preptr[STS_SIZE] = {null},             
      fleptr[MAX_PERIOD][20] = {{null}};

   static POSTSTATE *state = &eval;
   static GAUSS_JRD *gjp = &gss;
   static GRAPHICS *gpt = &gph;
   static SPLINES *spt = &spl;
   static FFT *fpt = &fft;
   static EVALUATE *vpt = &evl;

   static short 
      lbl = null,
      ind = null,
      lgx = null,
      lgy = null,
      frqlbl = null,

      extr[MAX_PERIOD+ONE] = {null};

   static long
      h_ = null,
      hh = null,
      ii = null,
      jj = null,
      kk = null,
      ll = null,
      mm = null,
      nn = null,
      pp = null,
      qq = null,
      init = null,
      final = null;

   static long flb[EVL_SPF+ONE] = {null};

   static const char
     *clsmsg = "DSC",
     *outmsg = "output",
     *evlptr = EVALUATION_FILE,
     *option = "forward", /* fourier transform option */
     *i_format = "%ld\n",
     *scformat = "%80s",
     *spformat = "%s\n";

   static char 
      fpformat[VSS_SIZE] = {null},
      unit_abs[SHS_SIZE] = {null},
      unit_ord[SHS_SIZE] = {null},
      pltptr[SHS_SIZE] = {null},
      datptr[SHS_SIZE] = {null},
      ctmptr[STS_SIZE] = {null},
      tmestr[STS_SIZE] = {null};

   static char 
     *t_type[MAX_PERIOD],
     *t_text[MAX_PERIOD],
     *s_type[MAX_PERIOD],
     *s_text[MAX_PERIOD],
     *timefrm = " created: %.24s ",
     *prefix = "gnu.",

    **endp = null;

   static double 
      xx = ZERO,
      yy = ZERO,
      zz = ZERO,
      x1 = ZERO,
      x2 = ZERO,
      dx = ZERO,
      norm = ZERO,
      real = ZERO,
      imag = ZERO,
      evl0 = ZERO,
      evl1 = ZERO,
      dlt0 = ZERO,
      dlt1 = ZERO,
      xlower = ZERO,
      xupper = ZERO,
      frequency = ZERO;

   static double
      x_min   = ZERO,
      y_min   = ZERO,
      z_min   = ZERO,
      x_max   = ZERO,
      y_max   = ZERO,
      z_max   = ZERO, 
      x_mean  = ZERO, 
      y_mean  = ZERO,
      z_mean  = ZERO,
      uu      = ZERO,
      vv      = ZERO,
      ww      = ZERO;

   static double
      rspc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      ispc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      aspc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      absmaxm[MAX_PERIOD] = {ZERO},
      timemax[MAX_PERIOD] = {ZERO},
      absmean[MAX_PERIOD] = {ZERO},
      meanenv[MAX_PERIOD] = {ZERO},
      modulation[MAX_PERIOD] = {ZERO},
      maxspec[MAX_PERIOD] = {ZERO},
      freqmax[MAX_PERIOD] = {ZERO},
      frq[EVL_SPF] = {ZERO};

# if DSC_HCRMDE != 0
   static short
      cc = null;
# endif

/* system parameters and function prototypes */

   time_t nseconds = null;
   time_t   *timer = null;
   time_t time(time_t *timer);

   SPLINES *
      spline( SPLINES *spt );

   int 
      graphp( GRAPHICS *gpt );

   char
      *lotos( long mm, char cc );

# ifndef _CCBUG
   char
      *strcpy( char *ptr1, const char *ptr2 ),  
      *strcat( char *ptr1, const char *ptr2 ),
      *strncat( char *ptr1, const char *ptr2, size_t n );
# endif

   static TXCNSL
     *csp;

/* mathematical function prototypes */

   double 
      sqrt( double x ),
      ceil( double x ),
      floor( double x );

/* user defined function prototypes: */ 

   FFT
     *fftrf( FFT *fpt );

   GAUSS_JRD
     *gssjrd( GAUSS_JRD *gjp );

   POSTSTATE
     *modval( POSTSTATE *state ),
     *pstprc( POSTSTATE *state );

   EVALUATE
     *readval( POSTSTATE *state, char option );

   TXCNSL
     *txcnsl( TXCNSL *csp );

# ifdef GPH_RNDOFF
   double rndoff( double xx, short nn );    
# endif
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* get the terminal type: */

   term = ( char *) getenv( "TERM" ); /* get the terminal */

   ind = tgetent( null, term );

   if( ONE != ind )
   {
      fprintf( stderr, "Error on getting the termcap info\n" ); 
      exit( EXIT_FAILURE );
   };
# endif
/*............................................................................*/
/* set buffer length = null: */

   ind = setvbuf( stdin, null, _IONBF, null );
   ind = setvbuf( stdout, null, _IONBF, null );
/*............................................................................*/
/* intialize text console */

   csp = txcnsl( null );
/*...........................................................................*/
/* port evaluation [ storage ] periods */
/* assignments, plot formats [initializations]: */

   ( state->vpt ) = vpt;
   ( vpt->fpt ) = fpt;
   
   strcpy( gpt->format, PLOT_FORMAT );
/*............................................................................*/
/* memory allocations: */

   hh = null; do
   {
      t_type[hh] = ( char *) calloc( SHS_SIZE, ONE );
      t_text[hh] = ( char *) calloc( STS_SIZE, ONE );
      s_type[hh] = ( char *) calloc( SHS_SIZE, ONE );
      s_text[hh] = ( char *) calloc( STS_SIZE, ONE );
   } while (( ++hh ) < MAX_PERIOD );

   if ( null == strncmp( PLOT_FORMAT, "SPLINE", THREE ))
      strncpy( fpformat, "%+.15E%s", TEN );
   else
      strncpy( fpformat, "%+.15e%s", TEN );

   strcpy(( vpt->mode_ep ), "individual" );
   strcpy(( vpt->mode_en ), "individual" );
   strcpy(( vpt->mode_hp ), "individual" );
   strcpy(( vpt->mode_hn ), "individual" );
/*............................................................................*/
/* Error messages: */

# if LINUX_C_SNTX != 1 
   hh = strlen (LNGSTR);
   if ( hh < ( MAX_PERIOD*( DEC+ONE ))) /* <- length of long string LNGSTR */
   {                                    /*    initializer in 'consts.h'    */
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD > %d !!!", 250/DEC );
      printf( "\n [ Filename pointer array 'fleptr[]' may perturb memory.]" ); 
   };
# endif
   mm = MAX_PERIOD;
   if ( FTR_NMBR < mm )                 /* <- FTR_NMBR: macro in POSTER.CONF */
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD = %d !!!", MAX_PERIOD );
      printf( "\n [ Required is MAX_PERIOD <= %d = FTR_NMBR"
              " <- macro defined in POSTER.CONF. ]", FTR_NMBR );
   };
   if ( FTR_NMBR < mm )                  /* <- FTR_NMBR: macro in POSTER.CONF */
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD = %d !!!", MAX_PERIOD );
      printf( "\n [ Required is MAX_PERIOD <= %d = FTR_NMBR"
              " <- macro defined in POSTER.CONF.]", FTR_NMBR );
   };
/*...........................................................................*/
/* initialize structure EVALUATE: */

   ( vpt->rtn ) = null;
   ( vpt->n ) = null;
   ( vpt->ni ) = null;
   ( vpt->nf ) = null;
   ( vpt->ti ) = ZERO;
   ( vpt->tf ) = ZERO;
   ( vpt->tc ) = ZERO;

# if DSC_HCRMDE != 0
   ( vpt->nj ) = null;
   ( vpt->nt ) = null;
   ( vpt->ctj ) = ZERO;
   ( vpt->ctf ) = ZERO;
   ( vpt->ctc ) = ZERO;
# endif
/*...........................................................................*/
/* port evaluation [ storage ] periods */

   hh = null; do
   {
      ( state->idx[hh] ) = null;
   } while (( hh++ ) < MAX_PERIOD );

   ( state->fldprd ) = null;
   ( state->period ) = null;

# if DSC_HCRMDE != 0
   ( state->hcrprd ) = null;
# endif
/*............................................................................*/
   strcpy( csp->cmmnt, "Welcome to DSC-POSTER !" );

# ifdef DSC_SYSTEM
   ( csp->dfopt ) = 7; /* initial default menu option */
# else
   ( csp->dfopt ) = 6;
# endif

   strcpy(( csp->cnfrm ), "Nothing done! Do you really want to quit ?" );
   ( csp->clscr ) = -ONE;

  menu:

   ( csp->dflnf ) = ( csp->dfopt ) ; /* set special line feed before default */
/*
   if (( csp->clscr ) == null )
      ( csp->clscr ) = ONE;
*/
   strcpy(( state->file ), DSC_PRFX );
   strcat(( state->file ), evlptr );

   if ( null == strlen( csp->cmmnt ))
      strcpy(( csp->cmmnt ), "Welcome back to DSC-POSTER !" );

   strcpy(( csp->envmt ), "DSC-POSTER" );
   strcpy(( csp->tasks ), "EVALUATION options ...:" );
   strcpy(( csp->mline[1] ), "* display evaluation scheme "
      ">-------> [ PREFX.val<n> ]" );
   strcpy(( csp->mline[2] ), "* create time response files "
      ">------> [ type time<n> ]" );
   strcpy(( csp->mline[3] ), "* create spectrum files "
      ">-----------> [ type spec<n> ]" );
   strcpy(( csp->mline[4] ), "* evaluate special parameters "
      ">-----> [ abs.max, -mean ...]" );
   strcpy(( csp->mline[5] ), "* evaluate special parameters on "
      "Fourier transforms" );
   strcpy(( csp->mline[6] ), "* any charged postprocessing function "
      "[ smithchrt(*), e.g.]" );

# ifdef DSC_SYSTEM
   strcpy(( csp->mline[7] ), "* DSC model dependent job evaluation"
      " >>[ function modval(*)]" );
   ( csp->items ) = 7;
# else
   ( csp->items ) = 6;
# endif

   strcpy(( csp->escpe ), "End of program / escape:" );
/*............................................................................*/
   csp = txcnsl( csp );   /* call menu building function                      */
/*......................*/
   lbl = ( csp->option );
/*............................................................................*/
   if ( lbl == null ) 
   {
      PRBLDCLR( "\r" );
      printf( "\r %*s", 78, "DSC-POSTER" );
      PRNORMAL( "\r " );
      return null;
   }
   else if ( lbl == ( EVL_ITEMS - ONE ))
   {
      PRBLDCLR( "\r" );
      printf( "\r %*s", 78, "DSC-POSTER" );
      PRNORMAL( "\r" );
/*............................................................................*/
      state = pstprc( state );       /* any postprocessing especially         */
/*...............................*//*   implemented in that function          */
      if (( state->rtn ) == ONE )
      {
         printf( "\n\n Message from function %s :", __func__ );
         printf( "\n Error on calling function pstprc(*) !!!" );
         printf( "\n [ Please verify code.]");
      };
      ( csp->dfopt ) = null; /* next default menu option: End of program  */
      goto menu;
   };
   printf( "\n" );

   strcpy( ptr, "Please enter index N of values file " );
   strcat( ptr, DSC_PRFX );
   strcat( ptr, EVALUATION_FILE );
   strcat( ptr, "<N>" );
   strcpy( csp->rqlng, ptr );
/*............................................................................*/
   csp = txcnsl( csp );     /* input on text console                          */
/*........................*/
   strcpy( ptr, lotos(( csp->inlng ), null ));
   strncpy(( state->flbl ), ptr, THREE ); 
   strcat(( state->file ), ( state->flbl ));

  open_eval:

   evalfle = fopen(( state->file ), "r+" );

   if ( evalfle == null )
   {
      printf( "\n Error on opening values file %s\n ", ( state->file ));
      strcpy( csp->rqstr, "Please re-enter filename [ Escape: enter null ]" );
/*............................................................................*/
      csp = txcnsl( csp );    /* input on text console                        */
/*..........................*/
      strcpy(( state->file ), ( csp->instr ));

      if ( *( state->file ) == '0' )
      {
         ( csp->dfopt ) = null; /* next default menu option: End of program  */
         goto menu;
      }
      else
         goto open_eval;
   } 
   else  /* file DSC_PRFX.val<n> reading section */
   {  
      fscanf( evalfle, scformat, ptr );
      strncpy( vpt->name, ptr, SHS_SIZE );

      if ( *( vpt->name ) == ' ' )
      {
         printf( "\n Message from function %s :", __func__ );
         printf( "\n Illegal format of file %s !!!", state->file );
         printf( "\n [ Cannot read DSC system specification in line 1.]" );

         fclose( evalfle );
         
         ( csp->dfopt ) = null; /* next default menu option: End of program  */

         goto menu;
      };
/*............................................................................*/
/* operational instructions [ evaluated ports, computation modes, etc.]: */

      fscanf( evalfle, scformat, ptr );
      strncpy(( vpt->text ), ptr, STS_SIZE );

      fscanf( evalfle, scformat, ptr ); /* string "______________________..." */

      fscanf( evalfle, scformat, ptr ); /* string "iteration_cycles" */
      fscanf( evalfle, scformat, ptr ); /* long integer string */
      ( vpt->n ) = strtol( ptr, endp, DEC );
      
/* the first evaluated cycle, Maxwell field: */

      fscanf( evalfle, scformat, ptr ); /* string "1st_evlt'd_cycle_Mxwfld" */
      fscanf( evalfle, scformat, ptr );
      ( vpt->ni ) = strtol( ptr, endp, DEC );

      if (( vpt->ni ) < VAL_INITLBL )
          ( vpt->ni ) = VAL_INITLBL;

/* the last evaluated [ and computed ] cycle, Maxwell field: */

      fscanf( evalfle, scformat, ptr ); /* string "Lst_cmpt'd_cycle,_Mxwfld" */
      fscanf( evalfle, scformat, ptr );
      ( vpt->nf ) = strtol( ptr, endp, DEC );
      
/* internal repetition rate, Maxwell field: */

      fscanf( evalfle, scformat, ptr ); /* string "iterations/cycle,_Mxwfld" */
      fscanf( evalfle, scformat, ptr );
      ( vpt->r ) = strtol( ptr, endp, DEC );

      if (( vpt->r ) < null )
         ( vpt->r ) = null;
/*............................................................................*/
# if DSC_HCRMDE != 0

/* the first evaluated cycle, heat and fluids: */
      fscanf( evalfle, scformat, ptr ); /* string "1st_evlt'd_cycle,_Ht&Fld" */
      fscanf( evalfle, scformat, ptr );
      ( vpt->nj ) = strtol( ptr, endp, DEC );
      
      if (( vpt->nj ) < VAL_INITLBL )
          ( vpt->nj ) = VAL_INITLBL;

/* the last evaluated [ and computed ] cycle, heat and fluids: */

      fscanf( evalfle, scformat, ptr ); /* string "Lst_cmpt'd_cycle,_Ht&Fld" */
      fscanf( evalfle, scformat, ptr ); 
      ( vpt->nt ) = strtol( ptr, endp, DEC );
      
/* internal repetition rate [ heat and fluids]: */
      fscanf( evalfle, scformat, ptr ); /* string "iterations/cycle,_Ht&Fld" */
      fscanf( evalfle, scformat, ptr );
/*............................................................................*/
# if DSC_INTLCE == 1 /* interlaced Maxwell field and thermal computation */
                      /* within repetition loop */
      ( vpt->rc ) = ( vpt->r );
# else
      ( vpt->rc ) = strtol( ptr, endp, DEC );

      if (( vpt->rc ) < null )
         ( vpt->rc ) = null;
# endif
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* the following needs decoupled E/H field and current storage modes */
/* [ in different evaluation files, e.g.] - This option is not used, here */
/*
# if DSC_HCRMDE !=0
      cc = null; do
      {
         fscanf( evalfle, scformat, ptr );
         ( vpt->nc[cc] ) = strtol( ptr, endp, DEC );
         fscanf( evalfle, scformat, ptr );
         ( vpt->nic[cc] ) = strtol( ptr, endp, DEC );
         fscanf( evalfle, scformat, ptr );
         ( vpt->rc[cc] ) = strtol( ptr, endp, DEC );
         cc++;
      } while ( cc < DSC_HCRMDE );
# endif
*/
/*............................................................................*/
      fscanf( evalfle, scformat, ptr ); /* string "Maxwell_field_evaluation" */
      fscanf( evalfle, scformat, ptr ); /* string "evaluated_ports" */
      fscanf( evalfle, scformat, ptr ); /* string "number..." */

      fscanf( evalfle, scformat, ptr ); /* string "E-ports:" */
      fscanf( evalfle, scformat, ptr ); /* the number [ long integer string ] */
      ( vpt->nep ) = strtol( ptr, endp, DEC );

      if ( ONE < ( vpt->nep ))
      {
         fscanf( evalfle, scformat, ptr ); /* any string that contains */
                                       /* "average" or "individual" */
         if ( null != strstr( ptr, "average" ))
            strcpy(( vpt->mode_ep ), "average" );
         else 
            strcpy(( vpt->mode_ep ), "individual" );
      };

      fscanf( evalfle, scformat, ptr ); /* string "E-nodes:" */
      fscanf( evalfle, scformat, ptr );
      ( vpt->nen ) = strtol( ptr, endp, DEC);

      if ( ONE < ( vpt->nen )) 
      {
         fscanf( evalfle, scformat, ptr ); /* any string that contains */
                                       /* "average" or "individual" */
         if ( null != strstr( ptr, "average" ))
            strcpy(( vpt->mode_en ), "average" );
         else 
            strcpy(( vpt->mode_en ), "individual" );
      };

      fscanf( evalfle, scformat, ptr ); /* string "H-ports:" */
      fscanf( evalfle, scformat, ptr );
      ( vpt->nhp ) = strtol( ptr, endp, DEC );

      if ( ONE < ( vpt->nhp )) 
      {
         fscanf( evalfle, scformat, ptr ); /* any string that contains */
                                       /* "average" or "individual" */
         if ( null != strstr( ptr, "average" ))
            strcpy(( vpt->mode_hp ), "average" );
         else 
            strcpy(( vpt->mode_hp ), "individual" );
      };

      fscanf( evalfle, scformat, ptr ); /* string "H-nodes:" */
      fscanf( evalfle, scformat, ptr );
      ( vpt->nhn ) = strtol( ptr, endp, DEC );

      if ( ONE < ( vpt->nhn )) 
      {
         fscanf( evalfle, scformat, ptr ); /* any string that contains */
                                       /* "average" or "individual" */
         if ( null != strstr( ptr, "average" ))
            strcpy(( vpt->mode_hn ), "average" );
         else 
            strcpy(( vpt->mode_hn ), "individual" );
      };
/*...........................................................................*/
# if DSC_HCRMDE !=0
/*...........................................................................*/
# if DSC_FLDMDE == 0
      fscanf( evalfle, scformat, ptr ); /* string "Thermal_evaluation" */
# else
      fscanf( evalfle, scformat, ptr ); /* string "Heat_&_fluid_evaluation" */
# endif /* DSC_FLDMDE ... */
/*...........................................................................*/
      fscanf( evalfle, scformat, ptr ); /* string "evaluated_ports" */
      fscanf( evalfle, scformat, ptr ); /* string "number..." */

      cc = null; do
      {
/*............................................................................*/
/* evaluated heat current faces: */

         fscanf( evalfle, scformat, ptr ); /* string "heat_current" */
         fscanf( evalfle, scformat, ptr );
         ( vpt->nhc[cc] ) = strtol( ptr, endp, DEC );

         if ( EVLHC < ( vpt->nhc[cc] ))
         {
            printf( "\n\n Message from function %s :",
               __func__ );
            printf( "\n Too many heat current faces to be evaluated "
               "in DSC mesh %s !!!", vpt->name );
            printf( "\n [ Maximum number is %ld = macro EVLHC "
               "in %s.", ( long ) EVLHC, "POSTER.CONF" );
            printf( "\n - Change macro only in compliance "
               "with memory resources !\n " );
            return null;
         };
/*............................................................................*/
# if DSC_FCTEMP == 1
/* evaluated temperature faces: */

         fscanf( evalfle, scformat, ptr ); /* string "face_temperature" */
         fscanf( evalfle, scformat, ptr );
         ( vpt->ntf[cc] ) = strtol( ptr, endp, DEC );

         if ( EVLTF < ( vpt->ntf[cc] ))
         {
            printf( "\n\n Message from function %s :",
               __func__ );
            printf( "\n Too many temperature faces to be evaluated "
               "in DSC mesh %s !!!", vpt->name );
            printf( "\n [ Maximum number is %ld = macro EVLTF "
               "in %s.", (long) EVLTF, "POSTER.CONF" );
            printf( "\n - Change macro only in compliance "
               "with memory resources !\n " );
            return null;
         };
# endif /* DSC_FCTEMP == 1 */
/*............................................................................*/
/* evaluated temperature nodes: */

         fscanf( evalfle, scformat, ptr ); /* string "node_temperature" */
         fscanf( evalfle, scformat, ptr );
         ( vpt->ntn[cc] ) = strtol( ptr, endp, DEC );

         if ( EVLTN < ( vpt->ntn[cc] ))
         {
            printf( "\n\n Message from function %s :",
               __func__ );
            printf( "\n Too many temperature nodes to be evaluated "
               "in DSC mesh %s !!!", vpt->name );
            printf( "\n [ Maximum number is %ld = macro EVLTN "
               "in %s.", (long) EVLTN, "POSTER.CONF" );
            printf( "\n - Change macro only in compliance "
               "with memory resources !\n " );
            return null;
         };
/*............................................................................*/
# if DSC_FLDMDE != 0 
/* evaluated nodal flows: */

         fscanf( evalfle, scformat, ptr ); /* string "nodal_velocity" */
         fscanf( evalfle, scformat, ptr );
         ( vpt->nun[cc] ) = strtol( ptr, endp, DEC );

         if ( EVLUN < ( vpt->nun[cc] ))
         {
            printf( "\n\n Message from function %s :",
               __func__ );
            printf( "\n Too many nodal velocities to be evaluated "
               "in DSC mesh %s !!!", vpt->name );
            printf( "\n [ Maximum number is %ld = macro EVLUN "
               "in %s.", (long) EVLUN, "POSTER.CONF" );
            printf( "\n - Change macro only in compliance "
               "with memory resources !\n " );
            return null;
         };
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
         cc++;
      } while ( cc < DSC_HCRMDE );
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      fscanf( evalfle, scformat, ptr ); /* string ">>LABELS---...>>" */

      for ( ii=null; ii<( vpt->nep ); ii++ )
      {
         fscanf( evalfle, scformat, ptr );
         ( vpt->mep[ii] ) = strtol( ptr, endp, DEC );
         fscanf( evalfle, scformat, ptr );
         ( vpt->pep[ii] ) = strtol( ptr, endp, DEC );
      };

      for ( ii=null; ii<( vpt->nen ); ii++ )
      {
         fscanf( evalfle, scformat, ptr );
         ( vpt->men[ii] ) = strtol( ptr, endp, DEC );
         fscanf( evalfle, scformat, ptr );
         ( vpt->cen[ii] ) = ptr[null];
      };

      for ( ii=null; ii<( vpt->nhp ); ii++ )
      {
         fscanf( evalfle, scformat, ptr );
         ( vpt->mhp[ii] ) = strtol( ptr, endp, DEC );
         fscanf( evalfle, scformat, ptr );
         ( vpt->php[ii] ) = strtol( ptr, endp, DEC );
      };

      for ( ii=null; ii<( vpt->nhn ); ii++ )
      {
         fscanf( evalfle, scformat, ptr );
         ( vpt->mhn[ii] ) = strtol( ptr, endp, DEC );
         fscanf( evalfle, scformat, ptr );
         ( vpt->chn[ii] ) = ptr[null];
      };
/*............................................................................*/
# if DSC_HCRMDE != 0
      cc = null; do
      {
         for ( ii=null; ii<( vpt->nhc[cc] ); ii++ )
         {
            fscanf( evalfle, scformat, ptr );
            ( vpt->mhc[cc][ii] ) = strtol( ptr, endp, DEC );
            fscanf( evalfle, scformat, ptr );
            ( vpt->fhc[cc][ii] ) = strtol( ptr, endp, DEC );
         };
/*............................................................................*/
# if DSC_FCTEMP == 1
         for ( ii=null; ii<(vpt->ntf[cc]); ii++ )
         {
            fscanf( evalfle, scformat, ptr );
            ( vpt->mtf[cc][ii] ) = strtol( ptr, endp, DEC );
            fscanf( evalfle, scformat, ptr );
            ( vpt->ftf[cc][ii] ) = strtol( ptr, endp, DEC );
         };
# endif /* DSC_FCTEMP == 1 */
/*............................................................................*/
         for ( ii=null; ii<( vpt->ntn[cc] ); ii++ )
         {
            fscanf( evalfle, scformat, ptr );
            ( vpt->mtn[cc][ii] ) = strtol( ptr, endp, DEC );
         };
/*............................................................................*/
# if DSC_FLDMDE != 0 
         for ( ii=null; ii<(vpt->nun[cc]); ii++ )
         {
            fscanf( evalfle, scformat, ptr );
            ( vpt->mun[cc][ii] ) = strtol( ptr, endp, DEC );
            fscanf( evalfle, scformat, ptr );
            ( vpt->cun[cc][ii] ) = strtol( ptr, endp, DEC );
         };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
         cc++ ;
      } while( cc < DSC_HCRMDE );
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* end of evaluated cell, face and port labels input */

      fscanf( evalfle, scformat, ptr ); /* "DSC" [ file 'dsc.val' created: etc.] */

      if ( null != strncmp( ptr, clsmsg, THREE )) /* check string *clsmsg = "DSC" */
      {
         printf( "\n\n Message from function %s :", __func__ );
         printf( "\n Illegal format of file %s !!!", ( state->file ));
         printf( "\n [ Cannot read closure message 'DSC ...' ] " );

         fclose( evalfle );

         ( csp->dfopt ) = null; /* next default menu option: End of program  */
         goto menu;
      };

      printf( "\n " );
      printf( "%s ", ptr );

      ll=null; do
      {
         fscanf( evalfle, scformat, ptr );
         printf( "%s ", ptr );
      } while(( ++ll ) < 10 );

      printf( "\n " );

      fscanf( evalfle, scformat, ptr ); /* [ string "__________..." ] */
      fscanf( evalfle, scformat, ptr ); /* [ string "output" ] */

/* check string *outmsg="output" */ 

      if ( null != strncmp( ptr, outmsg, FIVE )) 
      {
         strcpy( preptr, ptr );
         fscanf( evalfle, scformat, ptr );

         if ( null == strcpy( ptr, preptr ))
         {
            printf( "\n Message on file %s: ", ( state->file ));
            printf( "\n %s", preptr );

            ll=null; do
            {
               printf( " %s", ptr );
               strcpy( preptr, ptr );
               fscanf( evalfle, scformat, ptr );
               ll++;
            } while (( ll < 5 )\
                   &&( null != strcmp( ptr, preptr )));
         };

         fclose( evalfle );

         ( csp->dfopt ) = null; /* next default menu option: End of program  */
         goto menu;
      };

      printf( "\n %s ", ptr ); /* print "output of program SOLVER.C ... " */

      ll=null; do
      {
         fscanf( evalfle, scformat, ptr );
         printf( "%s ", ptr );
      } while(( ++ll ) < 10 );

      fscanf( evalfle, scformat, ptr ); /* [ closure string "__________..." ] */

      fscanf( evalfle, scformat, ptr ); /* boundary_type: */
      strncpy(( vpt->bndtyp ), ptr, SHS_SIZE ); /* APERIODIC_STRUCTURE ,      */
      printf( "\n\n %s ", ( vpt->bndtyp ));     /* PERIODIC_STRUCTURE_        */

      if ( null == strncmp(( vpt->bndtyp ), "APERIODIC", SIX ))
         ( vpt->read ) = ONE;
      else if ( null == strncmp(( vpt->bndtyp ), "PERIODIC", SIX ))
         ( vpt->read ) = TWO;
      else
      {
         printf( "\n\n Error on values file %s :", ( state->file ));
         printf( "\n Can't read boundary type identifier '%s'.",
            ( vpt->domain ));
         printf( "\n [ Legal is `APERIODIC' or `PERIODIC' ] " );

         fclose( evalfle );

         ( csp->dfopt ) = null; /* next default menu option: End of program  */
         goto menu;
      };
     
      fscanf( evalfle, scformat, ptr ); /* str "Maxwell_field_computation:" */
      printf( "\n %s", ptr );

      fscanf( evalfle, scformat, ptr ); /* domain [ identifier ] */
      strncpy(( vpt->domain ), ptr, SHS_SIZE ); /* TIME_DOMAIN________ , or */
      printf( "\n\n %s", ( vpt->domain ));      /* FREQUENCY_DOMAIN___ */

      if ( null == strncmp(( vpt->domain ), "TIME_DOMAIN", TWO ))
      { 
         ( vpt->read ) = ONE; 
                                            
         fscanf( evalfle, scformat, ptr ); /* left bracket */
         fscanf( evalfle, scformat, ptr ); /* string "Timestep" */
         fscanf( evalfle, scformat, ptr ); /* '=' */

         fscanf( evalfle, scformat, ptr ); /* time step */
         ( vpt->dt ) = strtod( ptr, endp );

         ( vpt->tc ) = ( vpt->r )*( vpt->dt ); /* cycle time */
         ( vpt->ti ) = ( vpt->ni )*( vpt->tc );
         ( vpt->tf ) = ( vpt->ti ) +\
             (( vpt->nf )-( vpt->ni ))*( vpt->tc );
                                          
         fscanf( evalfle, scformat, ptr ); /* string "seconds" */
         fscanf( evalfle, scformat, ptr ); /* right bracket */

         printf( "\n [ Timestep = %.15e seconds ]", ( vpt->dt ));
      }
      else if ( null == strncmp( vpt->domain, "FREQUENCY_DOMAIN", TWO ))
      {
         ( vpt->read ) = TWO;
                                             
         fscanf( evalfle, scformat, ptr ); /* left bracket */
         fscanf( evalfle, scformat, ptr ); /* string "Frequency" */
         fscanf( evalfle, scformat, ptr ); /* '=' */

         fscanf( evalfle, scformat, ptr ); /* frequency */
         ( vpt->fr ) = strtod( ptr, endp );

         ( vpt->tc ) = ( double )( vpt->r ); /* cycle step */
         ( vpt->ti ) = ( double )( vpt->ni )*( vpt->tc );
         ( vpt->tf ) = ( vpt->ti ) +\
             (( vpt->nf )-( vpt->ni ))*( vpt->tc );
                                              
         fscanf( evalfle, scformat, ptr ); /* string "Hertz" */
         fscanf( evalfle, scformat, ptr ); /* right bracket */

         printf( "\n [ Frequency = %.15e Hertz ]", ( vpt->fr ));
      }
      else
      {
         printf( "\n\n Error on values file %s :", ( state->file ));
         printf( "\n Can't read TIME/FREQUENCY domain identifier `%s'.",
            ( vpt->domain ));
         printf( "\n [ Legal is 'TIME_DOMAIN' or 'FREQUENCY_DOMAIN' ]" );

         fclose( evalfle );

         ( csp->dfopt ) = null; /* next default menu option: End of program  */
         goto menu;
      };

      fscanf( evalfle, scformat, ptr ); /* strng "Maxwell_field_excitation:" */
      printf( "\n %s", ptr );
      fscanf( evalfle, scformat, ptr ); /* EM-field exc type [ identifier] */
      strncpy(( vpt->exctyp ), ptr, SHS_SIZE ); /* STEADY_STATE, DIRAC_PULSE,.*/
      printf( "\n %s", ( vpt->exctyp ));
/*...........................................................................*/
# if DSC_HCRMDE != 0
      fscanf( evalfle, scformat, ptr ); /* string "thermal_computation:", e.g.*/
      printf( "\n %s", ptr );

      fscanf( evalfle, scformat, ptr ); /* string "TIME_DOMAIN", e.g. */

      printf( "\n %s", ptr );

      fscanf( evalfle, scformat, ptr ); /* left bracket */
      fscanf( evalfle, scformat, ptr ); /* string "Timestep" */
      fscanf( evalfle, scformat, ptr ); /* '=' */
      fscanf( evalfle, scformat, ptr ); /* time step */

      ( vpt->cdt ) = strtod( ptr, endp );

      if ( null < ( vpt->rc ))
      {
         ( vpt->ctc ) = ( vpt->rc )*( vpt->cdt );
         ( vpt->ctj ) = ( vpt->nj )*( vpt->ctc );
         ( vpt->ctf ) = ( vpt->ctj ) +\
             (( vpt->nt )-( vpt->nj ))*( vpt->ctc );

         printf( "\n [ Time step  = %.15e seconds", ( vpt->cdt ));
         printf( "\n   Cycle time = %.15e seconds ]", ( vpt->ctc ));
      };

      fscanf( evalfle, scformat, ptr ); /* string "seconds" */
      fscanf( evalfle, scformat, ptr ); /* right bracket */

      fscanf( evalfle, scformat, ptr ); /* "therm-fluid_excitation:", e.g.*/
      printf( "\n %s", ptr );

      fscanf( evalfle, scformat, ptr ); /* therm-fl exc type [identifier] */
      strncpy(( vpt->excctp ), ptr, SHS_SIZE );
      printf( "\n %s", ( vpt->excctp ));
# endif /* DSC_HCRMDE != 0 */
/*...........................................................................*/
      fscanf( evalfle, scformat, ptr ); /* system name */
      strncpy(( state->name ), ptr, SHS_SIZE );
      printf( "\n\n %s", state->name );

      if ( null != strncmp( state->name, vpt->name, TEN ))
      {
         printf( "\n\n Error on values file %s :", state->file );
         printf( "\n DSC process output identifier '%s'", state->name );
         printf( "\n differs from system identifier '%s'"
            "at head of file !!!", vpt->name );

         fclose( evalfle );

         ( csp->dfopt ) = null; /* next default menu option: End of program  */
         goto menu;
      };

      fscanf( evalfle, scformat, ptr ); /* text string [ any comment ] */
      strncpy(( state->text ), ptr, STS_SIZE );
      printf( "\n %s", ( state->text ));

      printf( "\n\n Units, EM-field computation" );

      fscanf( evalfle, scformat, ptr );
      strncpy(( vpt->xunit ), ptr, SHS_SIZE );
      printf( "\n    abscisse: %s", ( vpt->xunit ));

      fscanf( evalfle, scformat, ptr );
      strncpy(( vpt->yunit ), ptr, SHS_SIZE );
      printf( "\n    ordinate: %s", ( vpt->yunit ));
/*............................................................................*/
# if DSC_HCRMDE != 0
/*............................................................................*/
# if DSC_FLDMDE == 0
      printf( "\n\n Units, thermal computation" );
# else /* DSC_FLDMDE != 0 */
      printf( "\n\n Units, thermal-fluid computation" );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      printf( "\n    abscisse: %s", "seconds" );

      if ( TEMPGGE < 1.00e-77 )
         printf( "\n    ordinate: %s", "Kelvin" );
      else
         printf( "\n    ordinate: %s", "deg_Celsius" );
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      fscanf( evalfle, scformat, ptr );
      fscanf( evalfle, scformat, ptr );
      fscanf( evalfle, scformat, ptr );
      fscanf( evalfle, scformat, ptr );

      ( state->fleofs ) = ftell( evalfle );

      fclose( evalfle );
/*............................................................................*/
/* check number of sample points: */
      
      ii = ( vpt->nf ) - ( vpt->ni ) + ONE;
/*............................................................................*/
# if DSC_HCRMDE != 0
      if ( ii < (( vpt->nt ) - ( vpt->nj ) + ONE ))
         ii = ( vpt->nt ) - ( vpt->nj ) + ONE;
# endif
/*............................................................................*/
      if ( FTR_SIZE < ii )
      {
         fprintf( stderr, "\n\n Warning message from function %s :", __func__ );
         fprintf( stderr, "\n Number of sample points = %ld exceeds maximum "
            "= %ld,", ( long ) ii, ( long ) FTR_SIZE );
         fprintf( stderr, "\n fixed by macro FTR_SIZE in configuration header "
            "POSTER.CONF !!!" );
         fprintf( stderr, "\n [ Change macro in compliance with "
            "memory resources.]\n " );

         exit( EXIT_FAILURE );
      };
/*............................................................................*/
/* DSC_PRFX.val<n> file displaying section: */

      if ( null == strncmp(( vpt->domain ), "TIME_DOMAIN", TWO ))
      {
         printf( "\n\n Time domain, Maxwell field.............: "
            "[ %.5e, %.5e ] %s", ( vpt->ti ), ( vpt->tf ), ( vpt->xunit ));
         printf( "\n Cycle time.............................:   %.15e %s",
            ( vpt->tc ), ( vpt->xunit ));
      }
      else
      {
         printf( "\n\n Computational domain, Maxwell field....: "
            "[ %.5e, %.5e ] %s", ( vpt->ti ), ( vpt->tf ), ( vpt->xunit ));
         printf( "\n Cycle length ..........................:   %.15e %s",
            ( vpt->tc ), ( vpt->xunit ));
      };
/*............................................................................*/
# if DSC_HCRMDE != 0
      if (( vpt->cdt ) < 1.0e+77 )
      {
/*............................................................................*/
# if DSC_FLDMDE == 0
         printf( "\n\n Time domain, thermal computation.......: "
            "[ %.5e, %.5e ] %s", ( vpt->ctj ), ( vpt->ctf ), "seconds" );
# else /* if DSC_FLDMDE != 0 */
         printf( "\n\n Time domain, thermal-fluid computation.: "
            "[ %.5e, %.5e ] %s", ( vpt->ctj ), ( vpt->ctf ), "seconds" );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
         printf( "\n Cycle time ............................:   %.15e %s",
            ( vpt->ctc ), "seconds" );
         if ( ONE < ( vpt->rc ))
            printf( "\n Internal time step ....................:   %.15e %s",
            ( vpt->cdt ), "seconds" );
      };
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      printf( "\n\n Total number of iterations ( cycles )..:   %ld",
         ( long ) ( vpt->n ));
      if ( ONE < ( vpt->ni ))
         printf( "\n First computed cycle, Maxwell field....:   %ld",
            ( long ) ( vpt->ni ));
      if (( vpt->nf ) < ( vpt->n ))
         printf( "\n Last computed cycle, Maxwell field.....:   %ld",
            ( long ) ( vpt->nf ));
      if ( ONE < ( vpt->r ))
         printf( "\n Internal repetition rate, Maxwell field:   %ld",
            ( long ) ( vpt->r ));
/*............................................................................*/
# if DSC_HCRMDE != 0
      if (( vpt->cdt ) < 1.0e+77 )
      {
/*............................................................................*/
# if DSC_FLDMDE == 0
         if ( ONE < ( vpt->nj ))
            printf( "\n First computed thermal cycle...........:   %ld",
               ( long ) ( vpt->nj ));
         if (( vpt->nt ) < ( vpt->n ))
            printf( "\n Last computed thermal cycle............:   %ld",
            ( long ) ( vpt->nf ));
         if ( ONE < ( vpt->rc ))
            printf( "\n Internal repetition rate, thermal cycle:   %ld",
            ( long ) ( vpt->rc ));
# else /* DSC_FLDMDE != 0 */
         if ( ONE < ( vpt->nj ))
            printf( "\n First computed thermal-fluid cycle.....:   %ld",
               ( long ) ( vpt->nj ));
         if (( vpt->nt ) < ( vpt->n ))
            printf( "\n Last computed thermal-fluid cycle......:   %ld",
            ( long ) ( vpt->nf ));
         if ( ONE < ( vpt->rc ))
            printf( "\n Internal repetition rate, t-fluid cycle:   %ld",
            ( long ) ( vpt->rc ));
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      };
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* determine periods [ of evaluated port storage ]: */
      
      ( state->fldprd ) = null;

      if ( *( vpt->mode_ep ) != 'a' )
         ( state->fldprd ) += ( vpt->nep ); /* number of evaluated E-ports */
      else if ( null < ( vpt->nep ))
         ( state->fldprd )++;

      if ( *( vpt->mode_en ) != 'a' )
         ( state->fldprd ) += ( vpt->nen ); /* number of evaluated E-nodes */
      else if ( null < ( vpt->nen ))
         ( state->fldprd )++;

      if ( *( vpt->mode_hp ) != 'a' )
         ( state->fldprd ) += ( vpt->nhp ); /* number of evaluated H-ports */
      else if ( null < ( vpt->nhp ))
         ( state->fldprd )++;

      if ( *( vpt->mode_hn ) != 'a' )
         ( state->fldprd ) += ( vpt->nhn ); /* number of evaluated H-nodes */
      else if ( null < ( vpt->nhn ))
         ( state->fldprd )++;

      ( state->period ) = ( state->fldprd );
/*............................................................................*/
# if DSC_HCRMDE != 0 /* thermal evaluation */
      ( state->hcrprd ) = ( state->fldprd );
      if ( null < ( vpt->rc ))
      {
         cc = null; do
         {                                         /* number of evaluated ... */
            ( state->hcrprd ) += ( vpt->nhc[cc] ); /* ... nodal heat currents */
/*............................................................................*/
# if DSC_FCTEMP == 1
            ( state->hcrprd ) += ( vpt->ntf[cc] ); /* ... face temperatures */
# endif
/*............................................................................*/
            ( state->hcrprd ) += ( vpt->ntn[cc] ); /* ... nodal temperatures */
         } while (( ++cc ) < DSC_HCRMDE );
         ( state->period ) = ( state->hcrprd );
/*............................................................................*/
# if DSC_FLDMDE != 0 /* fluid velocities */
         cc = null; do
         {                                         /* number of evaluated ... */
            ( state->period ) += ( vpt->nun[cc] ); /* ... nodal velocities */
         } while (( ++cc ) < DSC_HCRMDE );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      }; /* end if ( null < ( vpt->rc )) */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* evaluated cycles: [ Maxwell field ]: */

      ii = ( vpt->nf ) - ( vpt->ni ) + ONE;

/* initialize: */
      hh = null; do
      {
         ( fpt->t[hh] ) = ( vpt->ti );
         ( fpt->tt[hh] ) = ( vpt->tf );
         ( fpt->dt[hh] ) = ( vpt->tc );
         ( fpt->ttlg[hh] ) = ii;
      } while(( hh++ ) < ( state->period ));
/*............................................................................*/
      printf( "\n" );
/*............................................................................*/
      vpt = readval( state, 'd' );     /* display evaluation scheme           */
/*...................................*/
      if ( lbl == EVL_ITEMS )
      {
         PRBLDCLR( "\r" );
         printf( "\r %*s", 78, "DSC-POSTER" );
         PRNORMAL( "\r" );
         printf( "\n ======================================="
                    "=======================================" );
/*............................................................................*/
         state = modval( state );             /* model dependent evaluation   */
/*..........................................*//* function                     */
         if (( state->rtn ) == ONE )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on calling function modval(*) !!!" );
            printf( "\n [ Please verify code. ]");
         };
    
         ( csp->dfopt ) = 6; /* next default menu option: any postprocessing */
         goto menu;
      };
/*............................................................................*/
/* input: port indices to be evaluated: */

      if (( ONE < lbl )
        &&( lbl < EVL_ITEMS ))
      {
         if ((( vpt->nf ) < ( vpt->ni ))
           ||(( vpt->r ) <= null ))
            hh = ( state->fldprd );
         else
            hh = null;

         ( state->idx[null] ) = null;
         while ((( state->idx[null] ) + hh ) < ( state->period ))
         {
            if ( state->idx[null] == null )
            {
               printf( "\n Please enter indices to be evaluated "
                  "[ Escape: enter null ] >--->\n" );
            };

           next_index:

            printf( " >----> enter %6d. index ................"
               ".....................: ", ( state->idx[null]+ONE ));
            scanf( "%s", ptr );

            if ( ptr[0] == '0' )
               break;

            ( state->idx[( state->idx[null]+ONE )] ) = hh +\
               strtol( ptr, endp, DEC );

/* index out of domain ??? */

            if (( state->period ) < ( state->idx[( state->idx[null]+ONE )] ))
            { 
               printf( "\n This port is not evaluated in "
                  "selected file '%s' !!!", ( state->file ));
               printf( "\n\n >----> Enter another index >-"
                  "----------------------------------->\n" );
               goto next_index;
            };

            if ( null < ( state->idx[(( state->idx[null] )+ONE )] ))
            {
               if (( ONE < lbl )
                 &&( lbl < FOUR ))
               {
                  printf( " >----> enter name of file to be created "
                     "for index %-6d........: ", (( state->idx[null] )+ONE ));
                  scanf( "%s", ptr );
                  strcpy( fleptr[( state->idx[null] )], ptr );
               };

               ( state->idx[null] )++;

               if ( MAX_PERIOD <= ( state->idx[null] ))
               {
                  printf( "\n\n The last evaluation index %ld has been accep"
                     "ted.", ( long ) ( state->idx[( state->idx[null] )] ));
                  printf( "\n [ Maximum number is %ld ", ( long ) MAX_PERIOD );
                  printf( "\r\t\t\t\t\b\b\b\b\b= macro MAX_PERIOD "
                     "in '%s'.]\n", "POSTER.CONF" );

                  goto read_values;
               };
            }
            else if (( state->idx[(( state->idx[null] )+ONE )] ) <= null )
               goto read_values;
         }; /* while (( state->idx[null] ) < ( state->period )) */
      }; /* end if (( ONE < lbl )&&( lbl < EVL_ITEMS )) */ 

     read_values:
/*............................................................................*/
      vpt = readval( state, 'r' );      /*                                    */
/*....................................*/
   }; /* end else if ( evalfle != null ) */
      /* end of DSC_PRFX.val<n> file displaying and reading section */

   if (( state->idx[null] ) < null )
   {
      ( csp->dfopt ) = null; /* next default menu option: End of program */
      goto menu;
   };
/*----------- values DSC_PRFX.val<n> entered, for specified ports ------------*/









/*----------------------------------------------------------------------------*/
/* evaluation [ options labelled lbl ]:                                       */

   switch( lbl )
   {
     case null:
      return ONE;

     case 1:
      break;

     case 2: /* >------ time domain response file generation -----------> */ 

      printf( "\n\n Please wait a moment !" );
      printf( "\n [ Writing data on files %s ... ] ", fleptr[null] );
        
      hh = ONE;
      while ( hh <= ( state->idx[null] ))
      {
         h_ = hh - ONE;

         strcpy( t_type[h_], vpt->name );
         strcpy( t_text[h_], vpt->text );

         timefle[h_] = fopen( fleptr[h_], "w+" );
         fprintf( timefle[h_], spformat, t_type[h_] );
/*............................................................................*/
         if ( state->idx[hh] <= ( state->fldprd )) /* Maxwell field port */
         {
            ii = ( vpt->nf ) - ( vpt->ni ) + ONE;

            dx = ( vpt->tc );
            xx = ( vpt->ti );

	    ( fpt->dt[hh] ) = ( vpt->tc );
	    ( fpt->t[hh] ) = ( vpt->ti );
	    ( fpt->tt[hh] ) = ( vpt->tf );
	    ( fpt->ttlg[hh] ) = ii;

            strcpy( unit_abs, ( vpt->xunit ));
            strcpy( unit_ord, ( vpt->yunit ));

            fprintf( timefle[h_], "%s_%s%ld%s\n", t_text[h_],
               "[evaluated_Maxwfld_port", ( long ) state->idx[hh], "]" );
            fprintf( timefle[h_], spformat, unit_abs );
            fprintf( timefle[h_], spformat, unit_ord );

            fprintf( timefle[h_], fpformat, ( vpt->ti ), "\n" );
            fprintf( timefle[h_], fpformat, ( vpt->tf ), "\n" );
            fprintf( timefle[h_], fpformat, ( vpt->tc ), "\n" );
	 }; /* end if ( state->idx[hh] <= state->fldprd ) */
/*............................................................................*/
# if DSC_HCRMDE != 0 
         if ((( state->fldprd ) < ( state->idx[hh] )) /* heat and fluid port */
           &&(( state->idx[hh] ) <= ( state->hcrprd )))
         {
            ii = ( vpt->nt ) - ( vpt->nj ) + ONE;

            dx = ( vpt->ctc );
            xx = ( vpt->ctj );
	    
	    ( fpt->dt[hh] ) = ( vpt->ctc );
	    ( fpt->t[hh] ) = ( vpt->ctj );
	    ( fpt->tt[hh] ) = ( vpt->ctf );
	    ( fpt->ttlg[hh] ) = ii;

            strcpy( unit_abs, "seconds" );

            if ( TEMPGGE < 1.00e-77 )
               strcpy( unit_ord, "Kelvin" );
            else
               strcpy( unit_ord, "deg_Celsius" );

            fprintf( timefle[h_], "%s_%s%ld%s\n", t_text[h_],
               "[evaluated_thermal_port",
               ( long ) (( state->idx[hh] ) - ( state->fldprd )), "]" );
            fprintf( timefle[h_], spformat, unit_abs );
            fprintf( timefle[h_], spformat, unit_ord );

            fprintf( timefle[h_], fpformat, ( vpt->ctj ), "\n" );
            fprintf( timefle[h_], fpformat, ( vpt->ctf ), "\n" );
            fprintf( timefle[h_], fpformat, ( vpt->ctc ), "\n" );
         };
/*............................................................................*/
# if DSC_FLDMDE != 0 
         if ((( state->hcrprd ) < ( state->idx[hh] ))
           &&(( state->idx[hh] ) <= ( state->period )))
         {
            ii = ( vpt->nt ) - ( vpt->nj ) + ONE;

            dx = ( vpt->ctc );
            xx = ( vpt->ctj );
	    
	    ( fpt->dt[hh] ) = ( vpt->ctc );
	    ( fpt->t[hh] ) = ( vpt->ctj );
	    ( fpt->tt[hh] ) = ( vpt->ctf );
	    ( fpt->ttlg[hh] ) = ii;

            strcpy( unit_abs, "seconds" );
            strcpy( unit_ord, "meter/seconds" );

            fprintf( timefle[h_], "%s_%s%ld%s\n", t_text[h_],
               "[evaluated_fluid_port",
               ( long ) (( state->idx[hh] ) - ( state->hcrprd )), "]" );
            fprintf( timefle[h_], spformat, unit_abs );
            fprintf( timefle[h_], spformat, unit_ord );

            fprintf( timefle[h_], fpformat, ( vpt->ctj ), "\n" );
            fprintf( timefle[h_], fpformat, ( vpt->ctf ), "\n" );
            fprintf( timefle[h_], fpformat, ( vpt->ctc ), "\n" );
         };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* store number of computed cycles: */

         fprintf( timefle[h_], i_format, ( long ) ii );

         jj = null;
         while ( jj < ii )
         {
            yy = ( fpt->r[hh][jj] );
            zz = ( fpt->i[hh][jj] );

# ifdef GPH_RNDOFF
            yy = rndoff( yy, ( GPH_RNDOFF+ONE ));
            zz = rndoff( zz, ( GPH_RNDOFF+ONE ));
            fprintf( timefle[h_], "%+.*e  %+.*e\n",
               ( int ) GPH_RNDOFF, yy, ( int ) GPH_RNDOFF, zz );
# else
            fprintf( timefle[h_], fpformat, yy, "  " );
            fprintf( timefle[h_], fpformat, zz, "\n" );
# endif
            jj++ ;
         };
 
/* close file with date and time of creation: */

         nseconds = time( timer );
/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */

         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( timefle[h_], "\n%s%s%s\n%.24s", "DSC "
            "time response file '", fleptr[h_], "' terminated: ", tmestr );
         fprintf( timefle[h_], "\n" );

         fclose( timefle[h_] );

         printf( "\n %s%s%s%s", "DSC "
            "time response file '", fleptr[h_], "' terminated: \n ", tmestr );

# if EVL_TIMERSP == 1
/*----------------------------------------------------------------------------*/
/* graphics, time domain response: */

         x_min = 1.e+277;
         y_min = 1.e+277;
         z_min = 1.e+277;

         x_max = -x_min;
         y_max = -y_min;
         z_max = -z_min;
/*
*//* [ still, the following applies: ]
*//* # if DSC_HCRMDE == 0
*//*         xx = ( fpt->t[null] );
*//*         dx = ( fpt->dt[null] );
*//* # elif DSC_HCRMDE != 0
*//*         if( state->idx[hh] <= state->fldprd )
*//*         {
*//*            xx = ( fpt->t[null] );
*//*            dx = ( fpt->dt[null] );
*//*         }
*//*         elif( state->fldprd < state->idx[hh] )
*//*         {
*//*            xx = ( vpt->ctj );
*//*            dx = ( vpt->ctc );
*//*         };
*//* # endif
*/
/* determine mimnimum and maximum of range and domain, */
/* then eventually rescale: */ 

         /* still ii as above ! */

         jj = null; kk = null;
         while ( jj < ii )
         {
            yy = ( fpt->r[hh][jj] );
            zz = ( fpt->i[hh][jj] );

            if( x_max < xx )
               x_max = xx; 

            if( xx < x_min )
               x_min = xx; 

            if( y_max < yy )
               y_max = yy;

            if( yy < y_min )
               y_min = yy;

            if( z_max < zz )
               z_max = zz;

            if( zz < z_min )
               z_min = zz;
           
            xx += dx;
            jj++ ;
         };

         if ( y_max < z_max )
            y_max = z_max;

         if ( z_min < y_min )
            y_min = z_min;

         xx = ( x_max - x_min )/2.;
         yy = ( y_max - y_min )/2.;
         zz = ( z_max - z_min )/2.;

         x_mean = ( x_min + x_max )/2.;
         x_min = x_mean - 1.01*xx;
         x_max = x_mean + 1.01*xx;

         y_mean = ( y_min + y_max )/2.;
         y_min = y_mean - 1.01*yy;
         y_max = y_mean + 1.01*yy;

         z_mean = ( z_min + z_max )/2.;
         z_min = z_mean - 1.01*zz;
         z_max = z_mean + 1.01*zz;

         xx = fabs( x_max );

         if ( xx < fabs( x_min ))
            xx = fabs( x_min );

         if ( ZERO < xx ) /* floor := largest integer not geater than argument*/
            lgx = floor( log10( xx ));
         else
            lgx = null;

         if ( abs( lgx ) < EVL_UNSCALED )
         {
            uu = 1.;
            lgx = null;
         }
         else
            uu = pow( 10., -lgx ); /* = 10^(-lgx ) = exp( -lgx*log( 10.)) */

         x_min *= uu;
         x_max *= uu;

         yy = fabs( y_max );

         if ( yy < fabs( y_min ))
            yy = fabs( y_min );

         if ( ZERO < yy )
            lgy = floor( log10( yy ));
         else
            lgy = null;

         if ( abs( lgy ) < EVL_UNSCALED )
         {
            vv = 1.;
            lgy = null;
         }
         else
            vv = pow( 10., -lgy ); /* = 10^(-lgy ) = exp( -lgy*log( 10.)) */

         y_min *= vv;
         y_max *= vv;
/*
         zz = fabs( z_max );

         if ( zz < fabs( z_min ))
            zz = fabs( z_min );

         if ( ZERO < zz )
            lgz = floor( log10( zz ));
         else
            lgz = ZERO;

         if ( fabs( lgz ) < EVL_UNSCALED )
         {
            ww = 1.;
            lgz = ZERO;
         }
         else
            ww = pow( 10., -lgz );

         z_min *= ww;
         z_max *= ww;
*/
         ww = vv;

         strcpy( pltptr, prefix );
         strncat( pltptr, fleptr[h_], ( SHS_SIZE - strlen( prefix )));

         pltfle = fopen( pltptr, "w+" );

         printf( "\n opened: gnuplot system file '%s' ", pltptr );

         fprintf( pltfle, "set title '%s %s'\n", vpt->name, vpt->text );
         fprintf( pltfle, "set xrange [%.15e:%.15e]\n", x_min, x_max );
         fprintf( pltfle, "set yrange [%.15e:%.15e]\n", y_min, y_max );

         if ( lgx == null )
            fprintf( pltfle, "set xlabel '%s'\n", unit_abs );
         else
            fprintf( pltfle, "set xlabel '1.E%+d %s'\n",
               lgx, unit_abs );

         if ( lgy == null )
            fprintf( pltfle, "set ylabel '%s'\n", unit_ord );
         else
            fprintf( pltfle, "set ylabel '1.E%+d %s'\n",
               lgy, unit_ord );

         fprintf( pltfle, "set nogrid\n" );
         fprintf( pltfle, "set border\n" );
         fprintf( pltfle, "plot %c\n", 92 );

         strcpy( datptr, prefix );
         strcat( datptr, fleptr[h_] );
         strcat( datptr, "real" );

         if(( vpt->read ) == ONE )
         {
            fprintf( pltfle, "'%s' with lines\n", datptr );
         }
         else
         {
            fprintf( pltfle, "'%s' with lines,%c\n", datptr, 92 );

            strcpy( datptr, prefix );
            strcat( datptr, fleptr[h_] );
            strcat( datptr, "imag" );

            fprintf( pltfle, "'%s' with lines\n", datptr );
         };

         fprintf( pltfle, "pause -1 '[ hit return to continue ]'\n" );

/* close plotfile with date and time of creation: */

         nseconds = time( timer );
/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */

         strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

         TIMEFORM( tmestr, ctmptr );

         fprintf( pltfle, "\n# 2D-gnuplot system file '%s'\n", pltptr );
         fprintf( pltfle, "# created:%.24s", tmestr );
         fprintf( pltfle, "\n#\n" );

         fclose( pltfle );

         printf( CLEAR_LINE );
         printf( "\r gnuplot system file '%s'", pltptr );
         printf( timefrm, tmestr );
/*............................................................................*/
/* plot real part: */

         strcpy( datptr, prefix );
         strcat( datptr, fleptr[h_] );
         strcat( datptr, "real" );

         datfle = fopen( datptr, "w+" );

         if ( datfle == null )
         {
            printf( "\n\n Message from function %s : ", __func__ );
            printf( "\n Unknown error on opening file %s !!! ", datptr );
            printf( "\n [ overrides and returns to calling program ].\n" );
            return ONE;
         };
         printf( "\n opened: gnuplot data - file '%s' ", datptr );

         fprintf( datfle, "# %s\n", vpt->name );
/*............................................................................*/
         if( state->idx[hh] <= state->fldprd )
         {
            xx = uu*( vpt->ti );
            dx = uu*( vpt->tc );
         };
/*............................................................................*/
# if DSC_HCRMDE != 0
         if ( state->fldprd < state->idx[hh] )
         {
            xx = uu*( vpt->ctj );
            dx = uu*( vpt->ctc );
         };
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         /* still ii as above ! */

         jj = null;
         while ( jj < ii )
         {
            yy = vv*( fpt->r[hh][jj] );
# ifdef GPH_RNDOFF
            zz = rndoff( xx, ( GPH_RNDOFF+ONE ));
            yy = rndoff( yy, ( GPH_RNDOFF+ONE ));
            fprintf( timefle[h_], "%+.*e %+.*e\n",
               ( int ) GPH_RNDOFF, zz, ( int ) GPH_RNDOFF, yy );
# else
            fprintf( timefle[h_], fpformat, xx, " " );
            fprintf( timefle[h_], fpformat, yy, "\n" );
# endif
            xx += dx;
            jj++ ;
         };

         fclose( datfle );

         printf( CLEAR_LINE );
         printf( "\r gnuplot data - file '%s' ", datptr );
         printf( timefrm, tmestr );
/*............................................................................*/
/* plot imaginary part: */

         if(( vpt->read ) == TWO )
         {
            strcpy( datptr, prefix );
            strcat( datptr, fleptr[h_] );
            strcat( datptr, "imag" );

            datfle = fopen( datptr, "w+" );

            if ( datfle == null )
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Unknown error on opening file %s !!!", datptr );
               printf( "\n [ overrides and returns to calling program ]. \n" );
               return ONE;
            };

            printf( "\n opened: gnuplot data - file '%s' ", datptr );

            fprintf( datfle, "# %s\n", vpt->name );
/*............................................................................*/
            if ( state->idx[hh] <= state->fldprd )
            {
               dx = uu*( vpt->tc );
               xx = uu*( vpt->ti );
            };
/*............................................................................*/
# if DSC_HCRMDE != 0
            if ( state->fldprd < state->idx[hh] )
            {
               dx = uu*( vpt->ctc );
               xx = uu*( vpt->ctj );
            };
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
            /* still ii as above ! */

            jj = null;
            while ( jj < ii )
            {
               zz = ww*( fpt->i[hh][jj] );
# ifdef GPH_RNDOFF
               yy = rndoff( xx, ( GPH_RNDOFF+ONE ));
               zz = rndoff( zz, ( GPH_RNDOFF+ONE ));
               fprintf( datfle, "%+.*e %+.*e\n",
                  ( int ) GPH_RNDOFF, yy, ( int ) GPH_RNDOFF, zz );
# else
               fprintf( datfle, fpformat, yy, " " );
               fprintf( datfle, fpformat, zz, "\n" );
# endif
               xx += dx;
               jj++ ;
            };

            fclose( datfle );

            printf( CLEAR_LINE );
            printf( "\r gnuplot data - file '%s' ", datptr );
            printf( timefrm, tmestr );
         };
/*............................................................................*/
# endif /* EVL_TIMERSP == 1 */

         hh++ ;
      };/* next hh */

      break;

     case 3: /* >----------------- save spectral response ------------------> */
     case 5: /* >--- evaluate spectral response [ special frequencies ] ----> */

      if ( null < ( state->idx[null] ))
      {
         printf( "\n" );

         do
         {
            printf( "\n Please enter frequency interval "
               ">--> [ f1, f2 ] ( GHz ) >>\n" );

            if ( ZERO < xlower )
               ( csp->dfdbl ) = 1.e-9*xlower;

            strcpy(( csp->rqfrm ), "points" );
            strcpy(( csp->rqdbl ), ">----------------"
               "------------------------> [ enter f1 ]" );

            ( csp->lnlen ) = 79;
            ( csp->posit ) = 62;
/*............................................................................*/
            csp = txcnsl( csp );     /* [ text console ]                      */
/*.................................*/
            xlower = ( csp->indbl )*1.e+9;

            if ( ZERO < xupper )
               ( csp->dfdbl ) = 1.e-9*xupper;

            strcpy(( csp->rqfrm ), "points" );
            strcpy(( csp->rqdbl ), ">----------------"
               "------------------------> [ enter f2 ]" );

            ( csp->lnlen ) = 79;
            ( csp->posit ) = 62;
/*............................................................................*/
            csp = txcnsl( csp );     /* [ text console ]                      */
/*.................................*/
            xupper = ( csp->indbl )*1.e+9;

            printf( "\n frequency interval: [ %1.5e , %1.5e ] "
               "GHz\n ", ( 1.e-9*xlower ), ( 1.e-9*xupper ));
/*............................................................................*/
/* select frequencies: */

            if ( lbl == FIVE )
            {
               printf( "\n For special parameter evaluation:\n Select "
                  "any frequencies >--> [ escape: enter 0 ] ( GHz ) >>\n" );
   
               hh = null;
               frqlbl = null; do
               {
                  hh++;

                 reenter1:

                  strcpy(( csp->rqfrm ), "points" );

                  if( hh < 10 )
                     strcpy(( csp->rqdbl ), ">---------------> [ enter " );
                  else if ( hh < 100 )
                     strcpy(( csp->rqdbl ), ">--------------> [ enter " );
                  else
                     strcpy(( csp->rqdbl ), ">-------------> [ enter " );

                  strcat(( csp->rqdbl ), lotos( hh, null ));
                  strcat(( csp->rqdbl ), ". frequency; f1 <= f <= f2 ]" );

                  if( 1.0e-277 < frq[frqlbl] )
                     ( csp->dfdbl ) = 1.e-9*frq[frqlbl];

                  ( csp->lnlen ) = 79;
                  ( csp->posit ) = 62;
/*............................................................................*/
                  csp = txcnsl( csp );     /* [ text console ]                */
/*.......................................*/
                  frequency = ( csp->indbl )*1.e+9;

                  if (( xlower <= frequency )&&( frequency <= xupper )&&
                      ( frequency != ZERO ))
                  {
                     frq[frqlbl] = frequency;
                     frqlbl++ ;
                  }
                  else if ( frequency != ZERO )
                  {
                     printf( " Frequency out of domain - please reenter "
                        "%6ld. frequency [ GHz ]: ", hh );
                     goto reenter1;
                  };
               } while (( frequency != ZERO )&&( frqlbl < EVL_SPF ));
            }; /* end if lbl == FIVE */
/*............................................................................*/
/* confirm keyboard input: */

            printf( "\n" );

            strcpy(( csp->rqfrm ), "bracket" );
            strcpy(( csp->rqstr ), "Input correct ? >----------" );
            strcat(( csp->rqstr ), "---------------------> [y/n]" );
            strcpy(( csp->dfstr ), "y" );

            ( csp->lnlen ) = 79;
            ( csp->posit ) = 62;
/*............................................................................*/
            csp = txcnsl( csp );      /* [ text console ]                     */
/*..................................*/
            strcpy( ptr, ( csp->instr));

         } while (( *ptr == 'n' )||( *ptr == 'N' )); /* end do */
/*............................................................................*/
         strcpy( unit_ord, ( vpt->yunit ));
         strcat( unit_ord, "*" );
         strcat( unit_ord, ( vpt->xunit ));
         strcpy( unit_abs, "1/" );
         strcat( unit_abs, ( vpt->xunit ));
      };
/*............... evaluation mode / frequency domain specified ...............*/










/*............................................................................*/
/* Fourier transformations: */

      printf( "\n Please wait a moment !" );

      if ( lbl == THREE )
         printf( "\n [ Fourier transforms "
           "- writing data on files '%s' ... ]\n ", fleptr[null]);
      else if ( lbl == FIVE )
         printf( "\n [ Fourier transforms ]\n " );

      kk = ( fpt->ttlg[null] );
       
      nn = TWO;
      while (( nn < kk )
           &&( nn < FTR_SIZE ))
         nn *= TWO;

      mm = nn / TWO;
      pp = ( nn - kk ) / TWO;
      qq = pp + kk;

      hh = ONE;
      while( hh <= ( state->idx[null] ))
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

# if READ_REVERSE_ == 1
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
# endif
         ( fpt->tt[hh] ) = ( fpt->t[hh] ) + nn*( fpt->dt[null] );
         ( fpt->dt[hh] ) = ( fpt->dt[null] );
         ( fpt->ttlg[hh] ) = nn;
         ( fpt->mult[hh] ) = ONE;

         ( fpt->p ) = hh;
         ( fpt->q ) = hh;
         strcpy(( fpt->opt ), option);
/*............................................................................*/
         fpt = fftrf( fpt );                  /*  Fast Fourier transformation */
/*..........................................*/
         hh++ ;
      }; /* end while( hh <= state->idx[null] ) */
             
      dx = ( fpt->ds[null] );

# if WRITE_REVERSE == 1
      ii = mm;
      xx = ( fpt->s[null] );
# else
      ii = null;
      xx = ZERO;
# endif

      if( xx < xlower ) /* floor: the largest integer not geater than argument*/
         kk = ( long ) floor (( xlower - xx ) / dx );
      else
         kk = null;

      if( nn < kk )
         kk = nn;

      init = ii + kk;
      x1 = xx + kk*dx;

      if( xx < xupper ) /* ceil: the smallest integer not less than argument */
         kk = ( long ) ceil (( xupper - xx ) / dx );
      else
         kk = null;

      if( nn < kk )
         kk = nn;

      final = ii + kk + ONE;
      x2 = xx + kk*dx;

/* Fourier transformations terminated */
/*............................................................................*/










/*............................................................................*/
/* save spectral response: */

      if ( lbl == THREE )
      {
         hh = ONE;
         while ( hh <= state->idx[null] )
         {
            h_ = hh - ONE;

            strcpy( s_type[h_], vpt->name );
            strcpy( s_text[h_], vpt->text );

            specfle[h_] = fopen( fleptr[h_],"w+");

            fprintf( specfle[h_], spformat, s_type[h_] );
            fprintf( specfle[h_], "%s_%s%ld%s\n", s_text[h_],
               "[evaluated_index", ( long ) state->idx[hh], "]" );
            fprintf( specfle[h_], spformat, unit_abs );
            fprintf( specfle[h_], spformat, unit_ord );

            fprintf( specfle[h_], fpformat, x1, "\n" );
            fprintf( specfle[h_], fpformat, x2, "\n" );
            fprintf( specfle[h_], fpformat, dx, "\n" );
            fprintf( specfle[h_], i_format, ( long ) ( final - init ));

# if EVL_WRTFREQ == 1
            xx = x1;
# endif
            ii = init;
            jj = init;
            while ( ii < final )
            {
               if ( nn <= jj )
                  jj -= nn; 
# if EVL_WRTFREQ == 1
               fprintf( specfle[h_], fpformat, xx, "  " );
	       xx += dx;
# endif
               fprintf( specfle[h_], fpformat, fpt->r[hh][jj], "  " );
               fprintf( specfle[h_], fpformat, fpt->i[hh][jj], "\n" );

               jj++ ;
               ii++ ;
            };

/* close spectrum with date and time of creation: */

            nseconds = time( timer );
/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); ': */

            strncpy( ctmptr, ctime( &nseconds ), 24 );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

            TIMEFORM( tmestr, ctmptr );

            fprintf( specfle[h_], "\n%s%s%s\n%.24s", "DSC"
            " spectrum file '", fleptr[h_], "' terminated: ", tmestr );

            fprintf( specfle[h_], "\n" );

            fclose( specfle[h_] );

            printf( "\n\n %s%s%s%s", "DSC"
            " spectrum file '", fleptr[h_], "' terminated: \n ", tmestr );

# if EVL_FREQRSP == 1          
/*----------------------------------------------------------------------------*/
/* graphics, spectrum [ absolute value]: */

            xx = x1;
            ii = init;
            jj = init;
            kk = null;
            ll = null;
            while ( ii < final )
            {
               if ( nn <= jj )
                  jj -= nn; 

               yy = ( fpt->r[hh][jj] );
               zz = ( fpt->i[hh][jj] );

               spl.vct[kk][null] = xx;
               spl.vct[kk][ONE] = sqrt( yy*yy + zz*zz );

               mm = null; do
               {
                  spl.dmn[ll] = xx;
                  xx += ( dx / ( EVL_GPHINTP + ONE ));
                  ll++ ;
                  mm++ ;
               } while(( mm <= EVL_GPHINTP )&&( ii < ( final - ONE )));
               ii++ ;
               jj++ ;
               kk++ ;
            };

            if( GPH_POINTS <= ll )
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Graphics needs too many points !!!" );
               printf( "\n Number n = %ld exceeds macro GPH_POINTS = %ld "
               "in funtion %s. ", (long) ll, (long) GPH_POINTS, "graphp(*)" );
               printf( "\n [ Change macro in compliance "
                  "with memory resources.]\n" );

               exit( EXIT_FAILURE );
            };

            ( spt->mm ) = kk; /* number of support points */
            ( spt->nn ) = ll; /* number of interpolated points */
/*............................................................................*/
            spt = spline( spt );            /* spline interpolation:          */
/*.......................................*//*  spectrum [abs.value]           */
            if ( spt == null )                                      
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling spline function" );
               printf( "\n for computing spc.ref[%ld][] !!!", (long) h_ );
               printf( "\n [ program stopped.]\n" );

               exit( EXIT_FAILURE );
            }
            else
            {
               jj = null;
               while ( jj < ll )
               {
                  gph.vct[jj][null] = spl.dmn[jj];
                  gph.vct[jj][ONE]  = spl.fct[jj];
                  jj++ ;
               };
               strcpy( gph.file, fleptr[h_] );
               strcpy( gph.name, state->name );
               strcpy( gph.text, "spectrum_[abs.value]" );
               strcpy( gph.xunit, "Hz" );
               strcpy( gph.yunit, "volts" );
               ( gpt->nn ) = ll;
/*............................................................................*/
               ind = graphp( gpt );               /* graphics file:           */
/*..............................................*//* spcectrum [abs.value]    */
            };   
/*............................................................................*/
/* graphics, spectrum [real part]: */

            xx = x1;
            ii = init;
            jj = init;
            kk = null;
            ll = null;
            while ( ii < final )
            {
               if ( nn <= jj )
                  jj -= nn; 

               spl.vct[kk][null] = xx;
               spl.vct[kk][ONE]  = ( fpt->r[hh][jj] );
               mm = null;
               do  
               {
                  spl.dmn[ll] = xx;
                  xx += ( dx / ( EVL_GPHINTP + ONE ));
                  ll++ ;
                  mm++ ;
               } while(( mm <= EVL_GPHINTP )&&( ii < ( final - ONE )));
               ii++ ;
               jj++ ;
               kk++ ;
            };

            ( spt->mm ) = kk;
            ( spt->nn ) = ll;
/*............................................................................*/
            spt = spline( spt );            /* spline interpolation:          */
/*.......................................*//*  spectrum [abs.value]           */
            if ( spt == null )                                      
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling spline function" );
               printf( "\n for computing spc.ref[%ld][] !!!", (long) h_ );
               printf( "\n [ program stopped.]\n" );

               exit( EXIT_FAILURE );
            }
            else
            {
               jj = null;
               while ( jj < ll )
               {
                  gph.vct[jj][null] = spl.dmn[jj];
                  gph.vct[jj][ONE]  = spl.fct[jj];
                  jj++ ;
               };
               strcpy( gph.file, fleptr[h_] );
               strcat( gph.file, "r" );
               strcpy( gph.name, state->name );
               strcpy( gph.text, "spectrum_[real_part]" );
               strcpy( gph.xunit, "Hz" );
               strcpy( gph.yunit, "volts" );
               ( gpt->nn ) = ll;
/*............................................................................*/
               ind = graphp( gpt );              /* graphics file:            */
/*.............................................*//* spcectrum [abs.value]     */
            };   
/*............................................................................*/
/* graphics, spectrum [imaginary part]: */

            xx = x1;
            ii = init;
            jj = init;
            kk = null;
            ll = null;
            while ( ii < final )
            {
               if ( nn <= jj )
                  jj -= nn;

               spl.vct[kk][null] = xx;
               spl.vct[kk][ONE]  = ( fpt->i[hh][jj] );
               mm = null;
               do
               {
                  spl.dmn[ll] = xx;
                  xx += ( dx / ( EVL_GPHINTP + ONE ));
                  ll++ ;
                  mm++ ;
               } while(( mm <= EVL_GPHINTP )&&( ii < ( final - ONE )));
               ii++ ;
               jj++ ;
               kk++ ;
            };

            ( spt->mm ) = kk;
            ( spt->nn ) = ll;
/*............................................................................*/
            spt = spline( spt );            /* spline interpolation:          */
/*.......................................*//*  spectrum [abs.value]           */
            if ( spt == null )                                      
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling spline function" );
               printf( "\n for computing spc.ref[%ld][] !!!", (long) h_ );
               printf( "\n [ program stopped.]\n" );

               exit( EXIT_FAILURE );
            }
            else
            {
               jj = null;
               while ( jj < ll )
               {
                  gph.vct[jj][null] = spl.dmn[jj];
                  gph.vct[jj][ONE]  = spl.fct[jj];
                  jj++ ;
               };
               strcpy( gph.file, fleptr[h_] );
               strcat( gph.file, "i" );
               strcpy( gph.name, state->name );
               strcpy( gph.text, "spectrum_[imag.part]" );
               strcpy( gph.xunit, "Hz" );
               strcpy( gph.yunit, "volts" );
               ( gpt->nn ) = ll;
/*............................................................................*/
               ind = graphp( gpt );              /* graphics file:            */
/*.............................................*//* spcectrum [abs.value]     */
            };   
/*............................................................................*/
# endif
            hh++ ;
         };/* next hh */

         break;

      };/* end if lbl == THREE */
/*.................. end case 3 [ spectral response saved ] ..................*/










/*............................................................................*/
/* case 5, cont'd [ evaluate special parameters ]:                            */

      hh = ONE;
      while ( hh <= state->idx[null] )
      {
         h_ = hh - ONE;

         maxspec[h_] = ZERO;
         freqmax[h_] = ZERO;

         for ( ll = null ; ll < frqlbl ; ll++ )
            aspc[h_][ll] = ZERO;

         frequency = x1;

         ii = init;
         jj = init;
         while ( ii < final )
         {
            if ( nn <= jj ) 
               jj -= nn;

            real = fpt->r[hh][jj];
            imag = fpt->i[hh][jj];
            norm = sqrt( real*real + imag*imag);

            if ( maxspec[h_] < norm )
            {
               maxspec[h_] = norm;
               freqmax[h_] = frequency;
            };

            if ( ii < ( final - ONE ))
            {
               kk = jj + ONE;

               if ( nn <= kk )
                  kk -= nn;

               for( ll=null; ll<frqlbl; ll++ )
               {
                  if (( frequency <= frq[ll] )&&
                      ( frq[ll] <= ( frequency + ( fpt->ds[hh] ))))
/* interpolation: */
                  {                                      
                     xx = ( frq[ll] - frequency )/dx;
                     yy = ( frequency + dx - frq[ll] )/dx;

                                  /* dx = fpt->ds[hh], cf. above: */
                                  /* block 'if ( hh == ONE ){}'*/

                     rspc[h_][ll]  = ( xx*fpt->r[hh][kk] + yy*fpt->r[hh][jj] );
                     ispc[h_][ll]  = ( xx*fpt->i[hh][kk] + yy*fpt->i[hh][jj] );
                     aspc[h_][ll]  = ( rspc[h_][ll]*rspc[h_][ll] );
                     aspc[h_][ll] += ( ispc[h_][ll]*ispc[h_][ll] );
                     aspc[h_][ll]  = sqrt( aspc[h_][ll] );
                  }; /* end if freq...  */
               };/* next ll [ frq label]*/
            };/* end if ( ii < final - ONE ) */ 

            frequency += ( fpt->ds[hh] );
            jj++ ;
            ii++ ;
         };/* end while ( ii < final ) */
         freqmax[h_] /= 1.e+09; /* in GHz */

         hh++ ;

      };/* next hh */ 

/* tab5: */

      printf( "\n\n Special evaluated spectral parameters "
         "of file %s :", ( state->file ));
      printf( "\n\n  index|abs. maximum :at frequency " );

      for ( ll=null ; (( ll <frqlbl )&&( ll<THREE )) ; ll++ )
         printf( "|ampl.at GHz->" );

      printf( "\n       |[%.13s]       [%.3s]", unit_ord, "GHz" );

      for ( ll=null ; (( ll<frqlbl )&&( ll<THREE )) ; ll++ )
         printf( "|%.7e", ( 1.e-09*frq[ll] )); 
         
      printf( "\n --------------------------------------"
              "----------------------------------------" );

      for ( hh = ONE; hh <= state->idx[null]; hh++ )
      {
         h_ = hh - ONE;

         printf( "\n %6d|%.7e:%.7e", state->idx[hh],
                      maxspec[h_], freqmax[h_] );

         for ( ll=null; (( ll<frqlbl )&&( ll<THREE )); ll++ )
            printf( "|%.7e", aspc[h_][ll] );
      };

      jj = ( int ) (( frqlbl - THREE + ( FIVE - ONE )) / FIVE );

      ii = null;
      while ( ii < jj ) 
      {
         ll = THREE + ii*FIVE;
         if ( ll < frqlbl )
         {
            printf( "\n\n  index" );
            while (( ll < THREE + ( ii + ONE )*FIVE )&&( ll < frqlbl ))
            {
               printf( "|ampl.at GHz->" );
               ll++;
            };
         };
         ll = THREE + ii*FIVE;
         if ( ll < frqlbl )
         {
            printf( "\n       " );
            while (( ll < THREE + ( ii + ONE )*FIVE )&&( ll < frqlbl ))
            {
               printf( "|%.7e", ( 1.0e-09*frq[ll] ));
               ll++;
            };
         };
         printf( "\n --------------------------------------"
                 "----------------------------------------" );

         for ( hh=ONE; hh<=state->idx[null]; hh++ )
         {
            h_ = hh - ONE;
            ll = THREE + ii*FIVE; 
            if ( ll < frqlbl ) 
            {
               printf( "\n %6d", state->idx[hh] );  
               while (( ll < THREE + ( ii + ONE )*FIVE )&&( ll < frqlbl ))
               {
                  printf( "|%.7e", aspc[h_][ll] );
                  ll++;
               };
            };
         };/* next hh */
         ii++;
      };/* end while ii < jj */

/* menu5: */
      
      printf( "\n\n Save on '%s' file ... ", PLOT_FORMAT );
      printf( "\r\t\t\t\t\t\t          Select  | " );
      printf( "\n\t\t\t\t\t\t                  V " );

      printf( "\n * absolute maxima >-----------"
         "--------------------------------> [1]" ); 

      for ( ll=null; ll<frqlbl; ll++ )
      {
         printf( "\n * amplitudes at %.13e GHz >----------"
            "-----------> [%1ld]", ( 1.0e-09*frq[ll] ), ( ll+TWO ));
      }; 
      printf( "\n\n None / Return to main menu >-"
         "---------------------------------> [0]" );

      printf( "\n\n Select values to be saved "
         "[ Return: enter null ] >------------> \n" );

      flb[null] = null; do
      {
         printf( " >----> enter %6ld. label ..............."
            "....................: ", ( flb[null] + ONE ));
         scanf( "%s", ptr );

         ll = strtol( ptr, endp, DEC );
         if (( null < ll )&&( ll < frqlbl + TWO ))
         {
            printf( " >----> enter name of file to be created "
               "for label %-6d......: ", ( short ) ll );
            scanf( "%s", ptr );
              
            strcpy( fleptr[flb[null]], ptr );
            flb[null]++;
            flb[flb[null]] = ll;
         };
      } while (( null < ll )&&( flb[null] < MAX_PERIOD ));

      for ( ll=ONE; ll<=flb[null]; ll++ )
      {
         strcpy( gph.file, fleptr[ll-ONE] );
         strcpy( gph.name, vpt->name );

         if ( ONE == flb[ll] ) 
            strcpy( gph.text, "[maxima]" );

         if ( TWO <= flb[ll] ) 
            strcpy( gph.text, "---" );

         strcpy( gph.xunit, "index" );
         strcpy( gph.yunit, unit_ord );

         for ( hh=ONE; hh<=state->idx[null]; hh++ )
         {
            h_ = hh - ONE;
            gph.vct[h_][null] = state->idx[hh];

            if ( flb[ll] == ONE )  
               gph.vct[h_][ONE] = maxspec[h_];
            else if ( ONE < flb[ll] )
               gph.vct[h_][ONE] = aspc[h_][flb[ll]-TWO];
         };

         ( gpt->nn ) = state->idx[null];
/*............................................................................*/
         ind = graphp( gpt );       /* generate plot file(s) */
/*................................*/

         if ( ind == null )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on calling function graphp(*) !!!" );
            printf( "\n [ Please verify source code. ]" );
         };

      }; /* next ll [ file label ] */

      break;
/*................ end case 5 [ evaluate spectral response ] .................*/










/*.................. case 4 [ evaluate special frequencies ] .................*/
     case 4:

      jj = null;

      for ( hh=ONE; hh<=( state->idx[null] ); hh++ )
      {
         h_ = hh - ONE;

        absmaxm[h_] = ZERO;
        absmean[h_] = ZERO;
        meanenv[h_] = ZERO;
        timemax[h_] = ZERO;
           extr[h_] = null;

         evl1 = ZERO;
         dlt1 = ZERO;

         for ( ii=null; ii<( fpt->ttlg[null] ); ii++ )
         {
            if( (vpt->read) == ONE )
            {
               evl0 = ( fpt->r[hh][ii] );
            }
            else
            {
               xx = ( fpt->r[hh][ii] );
               yy = ( fpt->i[hh][ii] );

               evl0 = sqrt( xx*xx + yy*yy );
            };

            dlt0 = evl0 - evl1;
            norm = fabs( evl0 );

            absmean[h_] += norm;

            if ( absmaxm[h_] < norm )
            {
               absmaxm[h_] = norm;                    
               timemax[h_] = ( fpt->t[null] ) + ii*( fpt->dt[null] );
            };

            /* ( vpt->read ) = 1: rel. extremum */
            /* ( vpt->read ) = 2: rel. maximum of abs val */
            /* compute mean.enveloppe */

            if ((( (vpt->read) == ONE )&&( dlt0*dlt1 < ZERO ))||
                (( ZERO < dlt1 )&&( dlt0 < ZERO )&&( (vpt->read) == TWO )))
            {
               extr[h_]++;
               meanenv[h_] += fabs( evl1 );
            };

            evl1 = evl0;

            if ( null < ii ) 
               dlt1 = dlt0;

         }; /* next ii */

         absmean[h_] /= ( fpt->ttlg[null] );
    
         if ( null < extr[h_] )
         {
            meanenv[h_] /= extr[h_];
            modulation[h_] = 100.*(fabs( absmaxm[h_]-meanenv[h_] ))/meanenv[h_];
         }
         else
         {
            meanenv[h_] = ZERO;
            modulation[h_] = ZERO;
         };

/* oscillation amlitude */

         norm = absmean[h_];

         evl1 = ZERO;
         dlt1 = ZERO;

         jj = -ONE;

         for ( ii=null; ii<( fpt->ttlg[null] ); ii++ )
         {
            evl0 = ( fpt->r[hh][ii] );
            dlt0 = evl0 - evl1;

            if (( dlt0*dlt1 ) < - ( 1.e-13*norm )) /* rel.extremum ! */
            {
               pp = ( ii - TWO );
            /* xx = fpt->t[null] + pp*fpt->dt[null]; */
               xx = pp*( fpt->dt[null] );
               qq = null; do
               {
                  xx += ( fpt->dt[null] );
                  yy = ( fpt->r[hh][pp] ); 
                  gss.mr[qq][0] = 1.;        
                  gss.mi[qq][0] = ZERO;
                  gss.mr[qq][1] = xx;     
                  gss.mi[qq][1] = ZERO;
                  gss.mr[qq][2] = xx*xx;
                  gss.mi[qq][2] = ZERO;
                  gss.mr[qq][3] = yy;
                  gss.mi[qq][3] = ZERO;
                  pp++;
                  qq++;
               }  while ( qq < THREE );

               ( gjp->rank ) = qq;
               ( gjp->opt ) = 'e'; /* option: 'e'quation */
/*............................................................................*/
               gjp = gssjrd( gjp );                    /* compute extremum in */
/*..................................................*//* parabolic approximtn */

               xx = -(.5*gss.zr[1][0]/gss.zr[2][0] ); /* abscissa of extremum */
               yy = gss.zr[0][0] + gss.zr[1][0]*xx + gss.zr[2][0]*xx*xx;
                                                     /* valu of extremum */
               if ( jj >= null )
               {
                  fpt->r[hh][jj] = .5*( xx + ww );
                  fpt->i[hh][jj] = fabs( yy - zz ); /* maximum oscillation */
               };
               ww = xx;
               zz = yy;
               jj++;
            };

            evl1 = evl0;

            if ( null < ii )
               dlt1 = dlt0;

         }; /* next ii */
         extr[h_] = jj;
      }; /* next hh */

     tabl4:

      printf( "\n\n Special parameters evaluated in file '%s':",
         ( state->file ));
      printf( "\n\n index | abs. maximum : found at time|    abs. mean|  "
         "mean envlp.|   modulation" );
      printf( "\n       | [%10s] : [ %10s]|    [%7s]|  [%9s]| %12s",
        (vpt->yunit), (vpt->xunit), (vpt->yunit), (vpt->yunit), "[ percent]" ); 
      printf( "\n --------------------------------------"
         "----------------------------------------" );
      for ( hh=ONE; hh<=( state->idx[null] ); hh++ )
      {
         h_ = hh - ONE;
         printf("\n %6d| %.7e: %.7e|%.7e|%.7e|%.7e", ( state->idx[hh] ), 
         absmaxm[h_], timemax[h_], absmean[h_], meanenv[h_], modulation[h_] );
      }; /* next hh */
      printf( "\n --------------------------------------"
              "----------------------------------------" );
/* menu4: */

      printf( "\n\n Save on '%s' file ... ", PLOT_FORMAT );
      printf( "\r\t\t\t\t\t\t          Select  | " );
      printf( "\n\t\t\t\t\t\t                  V " );

      printf( "\n * absolute maxima >--------------"
                "-----------------------------> [1]" );
      printf( "\n * absolute means >---------------"
                "-----------------------------> [2]" );
      printf( "\n * mean enveloppe >---------------"
                "-----------------------------> [3]" );
      printf( "\n * modulation >-------------------"
                "-----------------------------> [4]" );
      printf( "\n * oscillation >------------------"
                "-----------------------------> [5]" );
      printf( "\n\n None / Return to main menu >---"
              "-------------------------------> [0]" );
      printf( "\n\n                                  "
              "                                 ] <- ?" );
      printf( "\r Please select item >---------------->"
                   " [ enter 1,2,... or 0 ] >> [" );

      scanf( "%s", ptr );
      ll = strtol( ptr, endp, DEC );

      if (( ONE <= ll )
        &&( ll <= FIVE ))  
      {
         printf( " >----> Please enter file name >----"
                 "---------------------------> " );
         scanf( "%s", ptr );

         strcpy( gph.name, vpt->name );
         strcpy( gph.text, vpt->text ); 

         if ( ll <= FOUR )
         {
            strcpy( gph.file, ptr );
            strcpy( gph.xunit, "index" );
            strcpy( gph.yunit, vpt->yunit );

            for ( hh=ONE; hh<=( state->idx[null] ); hh++ ) 
            {
               h_ = hh - ONE;

               gph.vct[h_][null] = ( state->idx[hh] ); 

               switch ( ll ) 
               {
                 case 1:
                  gph.vct[h_][ONE] = absmaxm[h_];
                  break;

                 case 2: 
                  gph.vct[h_][ONE] = absmean[h_];
                  break;

                 case 3:
                  gph.vct[h_][ONE] = meanenv[h_]; 
                  break;

                 case 4: 
                  gph.vct[h_][ONE] = modulation[h_]; 
                  break;

                 default: 
                  break;
               };
            }; /* next hh */

            ( gpt->nn ) = ( state->idx[null] );
/*............................................................................*/
            ind = graphp( gpt );           /*                                 */
/*.......................................*/
         }
         else if ( ll == FIVE )
         {
            for ( hh=ONE; hh<=( state->idx[null] ); hh++ )
            {
               h_ = hh - ONE;

               strcpy( gph.file, ptr );
               strcat( gph.file, lotos( hh, null ));
               strcpy( gph.xunit, "seconds" );
               strcpy( gph.yunit, ( vpt->yunit ));

               for ( ii=null; ii<extr[h_]; ii++ )
               {
                  gph.vct[ii][null] = ( fpt->r[hh][ii] );
                  gph.vct[ii][ONE]  = ( fpt->i[hh][ii] );
               };

               ( gpt->nn ) = extr[h_];
/*............................................................................*/
               ind = graphp( gpt );              /*                           */
/*.............................................*/
            };
         };

         printf( "\n --------------------------------------"
                  "----------------------------------------" );

         goto tabl4;
      };

      break;

/*............................................................................*/
    
     default:
      break;

   }; /* end switch ( lbl ) */
/*............................................................................*/
   
   ( csp->dfopt ) = 7; /* next default menu option: End of program  */
   goto menu;
}     
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/*----------------------------------------------------------------------------*/
# undef EVL_ITEMS
# undef EVL_DISP
# undef EVL_GPHINTP
# undef EVL_TIMERSP
# undef EVL_FREQRSP
/********************** end of function postdrv(*) ****************************/
