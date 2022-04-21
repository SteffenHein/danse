/* [ file: excite.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function excite(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   DSC excitation parameter input function                                    *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPT_LEVEL 2
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
# define EXC_DISP 1
/*----------------------------------------------------------------------------*/
/* structure excitation exc: cf. function scattr.h                            */
/*----------------------------------------------------------------------------*/
/* There are different ways to interlace excited with incident quantities:    */
     
# define EXC_IMPOSE 0 /* 0: additive superposition [ recommended ]            */
                      /* 1: direct imposition                                 */
                      /* 2: imposition with scattered quantities subtracted   */
                      /*    [ fixes total quantities ].                       */
# define EXC_MAXFLD 1 /* 1: write maximum excited field into exc.mx           */
/*----------------------------------------------------------------------------*/
/* constants, stabilizing ONE [ to be modified with care ]:                   */

# ifndef PRECISION
   # define PRECISION ( double )( 1.000e-15 )
# endif

# if EXC_IMPOSE != 0
   # define ONE_MINUS ( 1. - 33.*PRECISION ) 
# endif
/*----------------------------------------------------------------------------*/
static short exclbl[DSC_JOBS] = {null};
/*============================================================================*/

short excite( const short jj )
{
/* allusions: */
/*
   extern struct solverstat solver;
   extern struct excitation exc;

# if DSC_DOMAIN == 1
   extern DSC_FIELDS fld;
# else
   extern DSC_FIELDS fld[];
# endif
*/
/* declarations: */

   static struct solverstat
     *state = &solver;

   static short 
      ii = null,
      kk = null;

/* prototypes: */

   char 
     *lotos ( long lngint, char length );
/*----------------------------------------------------------------------------*/
/* allusions: */
/*
   extern short exclbl[];
*/
/* declarations: */

   static FILE 
     *excitfle = NULL;

   static const short
      lbnd =  0, 
      cbnd =  8, /* <~ cbnd = clns - 2 */
      clns = 10; /* columns */

   static const char
      slsh = 47,
      bslsh = 92;

   static short
      pp = null,
      ind = null,
      llns = null;

# if EXC_MAXFLD == 1
   static double
      uu = ZERO,
      vv = ZERO;
# endif

   static char 
      ptr[STS_SIZE] = {null},
      type[SHS_SIZE] = {null},
      unit[SHS_SIZE] = {null},
     *excptr = EXCITATION_FILE,     /* cf. main program solver.c */
    **endp = NULL;

   static const char 
     *scformat = "%80s";
/*----------------------------------------------------------------------------*/
/* memory allocations: - */
/*............................................................................*/
/* parameter reset: */

   for ( ii=null; ii<exc.ne; ii++ )
   {
      exc.me[ii] = null;
      exc.pe[ii] = null;
      exc.er[ii] = null;
      exc.ei[ii] = null;
   };
   exc.ne = null;

   for ( ii=null; ii<exc.nh; ii++ )
   {
      exc.mh[ii] = null;
      exc.ph[ii] = null;
      exc.hr[ii] = null;
      exc.hi[ii] = null;
   };
   exc.nh = null;

   for ( ii=null; ii<exc.nn; ii++ )
   {
      exc.fr[ii] = ZERO;
   };

   exc.nn = null;

   exc.dt = ZERO;
   exc.ht = ZERO;
   exc.rt = ZERO;

   exc.ne = null;
   exc.nh = null;

   exc.mx = ZERO;
   exc.sq = ZERO;

# if DSC_HCRMDE != 0
   exc.hcn = null;
   exc.hcdt = ZERO;
   exc.hcht = ZERO;
   exc.hcrt = ZERO;

   exc.nhc = null;
   exc.ntf = null;
   exc.ntn = null;
# endif
/*............................................................................*/
# if DSC_LNGNAMES == 1
   strcpy(( state->exc ), ( state->prfx ));
   strcat(( state->exc ), excptr );
# else
   strncpy(( state->exc ), ( state->prfx ), VSS_SIZE );
   strncat(( state->exc ), excptr, ( SHS_SIZE - VSS_SIZE - THREE ));
# endif

   strcat(( state->exc ), lotos( exclbl[jj], null ));

   excitfle = fopen(( state->exc ), "r" );

   if ( excitfle == null )
   {
      fprintf( display, "\n\n Error on opening excitation file %s\n ",
         ( state->exc ));
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":unable_to_open_file_" , 
         ( LGS_SIZE - SHS_SIZE));
      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( excitfle, scformat, ptr );
   strncpy( exc.name, ptr, SHS_SIZE - ONE );
/*............................................................................*/
/* file check: */

   if ( null != strncmp( exc.name, top.name, THREE ))
   {
      fclose( excitfle );
      fprintf( display, "\n\n File error: Inconsistent system identifier "
         "in file %s !!!", ( state->exc ));
      fprintf( display, "\n --- Please verify ! ---\n" );
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":inconsistent_DSC_model_identifier_", 
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };
/*............................................................................*/
/* enter comment string: */

   fscanf( excitfle, scformat, ptr );
   strncpy( exc.text, ptr, ( STS_SIZE - ONE ));

   fscanf( excitfle, scformat, ptr ); /* string "________________________..." */
/*............................................................................*/
   fscanf( excitfle, scformat, ptr ); /* string "electric_excitation" */ 
/*............................................................................*/
/* enter time/frequency domain identifier string: */

   fscanf( excitfle, scformat, ptr );

   if ( null == strncmp( ptr, "TIME_DOMAIN", TWO ))
      pp = ONE;
   else if ( null == strncmp( ptr, "FREQUENCY_DOMAIN", TWO ))
      pp = TWO;
   else
   {
      fprintf( display, "\n\n Error message from function %s:",
         __func__ );
      fprintf( display, "\n Can't time/frequency domain identifier "
         "%s !!!", ptr );
      fprintf( display, "\n [ must be 'TIME_DOMAIN' or "
	 "'FREQUENCY_DOMAIN']." );
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":illegal_domain_",
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };
/*............................................................................*/
/* EM field excitation type: */

   fscanf( excitfle, scformat, ptr );
   strncpy ( type, ptr, SHS_SIZE );
   fscanf( excitfle, scformat, ptr ); /* string "parameters" */

   if ( null == strncmp( type, "ZERO_EM_FIELD______", TWO ))
   {
      strncpy ( exc.type, "ZERO_EM_FIELD______", SHS_SIZE );
      exc.lbl = 0;

      goto cont1;
   }
   else if ( null == strncmp( type, "STEADY_STATE_______", TWO ))
   {
      strncpy ( exc.type, "STEADY_STATE_______", SHS_SIZE );
      exc.lbl = 1;

      fscanf( excitfle, scformat, ptr ); 
      exc.nn = strtol( ptr, endp, DEC ); /* smoothing order */

      fscanf( excitfle, scformat, ptr );
      exc.dt = strtod( ptr, endp ); /* smoothing time */

      goto cont1;
   }
   else if ( null == strncmp( type, "DIRAC_PULSE________", TWO ))
   {
      strncpy ( exc.type, "DIRAC_PULSE________", SHS_SIZE );
      exc.lbl = 2;

      goto cont1;
   }
   else if ( null == strncmp( type, "HARMONIC_SINUSOIDAL", TWO ))
   {                                                 
      strncpy ( exc.type, "HARMONIC_SINUSOIDAL", SHS_SIZE );
      exc.lbl = 3;

      fscanf( excitfle, scformat, ptr );
      exc.fr[null] = strtod( ptr, endp ); /* exciting frequency [Hz]*/

      goto cont1;
   } 
   else if ( null == strncmp( type, "SMOOTH_HARMONIC____", TWO ))
   {
      strncpy ( exc.type, "SMOOTH_HARMONIC____", SHS_SIZE );
      exc.lbl = 4;

      fscanf( excitfle, scformat, ptr );
      exc.fr[null] = strtod( ptr, endp ); /* exciting frequency [Hz] */
      fscanf( excitfle, scformat, ptr );
      exc.dt = strtod( ptr, endp ); /* smoothing time */

      goto cont1;
   } 
   else if ( null == strncmp( type, "MULTIPLE_HARMONIC__", TWO ))
   {                                              
      strncpy ( exc.type, "MULTIPLE_HARMONIC__", SHS_SIZE );
      exc.lbl = 5;

      fscanf( excitfle, scformat, ptr );
      exc.nn = strtol( ptr, endp, DEC ); /* number of exc. frequences */

      if ( EXCFR < exc.nn )
      {
         fprintf( display, "\n\n Error message from function %s:",
            __func__ );
         fprintf( display, "\n Too many exciting frequencies defined "
            "in DSC mesh %s !!!", exc.name );
         fprintf( display, "\n [ Maximum number is %ld = macro EXCFR "
            "in SOLVER.CONF or %s.", ( long ) EXCFR, __func__ );
         fprintf( display, "\n - Change macro only in compliance "
            "with memory resources !\n " );
         strncpy(( state->fcterr ), __func__, SHS_SIZE );
         strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
         strncat(( state->errmsg ), ":too_many_frequencies_excited_", 
            ( LGS_SIZE - SHS_SIZE ));
         return null; /* abnormal return */
      };
/*............................................................................*/

      for ( ii=null; ii<exc.nn; ii++ )
      {
         fscanf( excitfle, scformat, ptr );
         exc.fr[ii] = strtod( ptr, endp ); /* exc. frequencies [Hz] */
      };
      goto cont1;
   }
   else if ( null == strncmp( type, "GAUSS_PULSE________", TWO ))
   {                                              
      strncpy ( exc.type, "GAUSS_PULSE________", SHS_SIZE );
      exc.lbl = 6;

      fscanf( excitfle, scformat, ptr );
      exc.rt = strtod( ptr, endp ); /* rise time */
      fscanf( excitfle, scformat, ptr );
      exc.dt = strtod( ptr, endp ); /* delay time */

      goto cont1;
   }
   else if ( null == strncmp( type, "OSCILLATORY_GAUSS__", TWO ))
   {                                             
      strncpy ( exc.type, "OSCILLATORY_GAUSS__", SHS_SIZE );
      exc.lbl = 7;

      fscanf( excitfle, scformat, ptr );
      exc.fr[null] = strtod( ptr, endp ); /* exciting frequency [Hz] */
      fscanf( excitfle, scformat, ptr );
      exc.rt = strtod( ptr, endp );       /* sigma time */
      fscanf( excitfle, scformat, ptr );
      exc.dt = strtod( ptr, endp );       /* delay time */

      goto cont1;
   }
   else if ( null == strncmp( type, "RAMP_PULSE_________", TWO ))  
   {                                              
      strncpy ( exc.type, "RAMP_PULSE_________", SHS_SIZE );
      exc.lbl = 8;

      fscanf( excitfle, scformat, ptr );
      exc.rt = strtod( ptr, endp ); /* rise time */

      goto cont1;
   }
   else if ( null == strncmp( type, "RECTANGULAR_PERIODC", TWO )) 
   {                                              
      strncpy ( exc.type, "RECTANGULAR_PERIODC", SHS_SIZE );
      exc.lbl = 9;

      fscanf( excitfle, scformat, ptr );
      exc.fr[null] = strtod( ptr, endp ); /* exciting frequency [Hz] */

      goto cont1;
   }
   else if ( null == strncmp( type, "SAW_TOOTH_PERIODIC_", TWO )) 
   {                                             
      strncpy ( exc.type, "SAW_TOOTH_PERIODIC_", SHS_SIZE );
      exc.lbl = 10;

      fscanf( excitfle, scformat, ptr );
      exc.fr[null] = strtod( ptr, endp ); /* exciting frequency [Hz] */
      exc.rt = 1./exc.fr[null];

      goto cont1;
   } 
   else if ( null == strncmp( type, "HEAVISIDE_STEP_____", TWO ))
   {
      strncpy ( exc.type, "HEAVISIDE_STEP_____", SHS_SIZE );
      exc.lbl = 11;

      goto cont1;
   }
   else if ( null == strncmp( type, "DSC_TIMESTEP_PRDC__", TWO ))
   {                                           
      strncpy ( exc.type, "DSC_TIMESTEP_PRDC__", SHS_SIZE );
      exc.lbl = 12;

      goto cont1;
   }
   else if ( null == strncmp( type, "NOISE_AT_RANDOM____", TWO )) 
   {                                           
      strncpy ( exc.type, "NOISE_AT_RANDOM____", SHS_SIZE );
      exc.lbl = 13;

      goto cont1;
   }
   else if ( null == strncmp( type, "MORLET_WAVELET_____", TWO ))
   {                                             
      strncpy ( exc.type, "MORLET_WAVELET_____", SHS_SIZE ); 
      exc.lbl = 14;

      fscanf( excitfle, scformat, ptr );
      exc.fr[null] = strtod( ptr, endp ); /* exciting frequency [Hz] */
      fscanf( excitfle, scformat, ptr );
      exc.rt = strtod( ptr, endp ); /* rise time */
      fscanf( excitfle, scformat, ptr );
      exc.dt = strtod( ptr, endp ); /* delay time */

      goto cont1;
   }
   else if ( null == strncmp( type, "WAVE_PACKET_HARMNC_", TWO ))
   {
      strncpy ( exc.type, "WAVE_PACKET_HARMNC_", SHS_SIZE );
      exc.lbl = 15;

      fscanf( excitfle, scformat, ptr ); /* smoothing order */
      exc.nn = strtol( ptr, endp, DEC ); /* [ ignored if exc.nn < 2 ] */
                      
      fscanf( excitfle, scformat, ptr ); /* smoothing time */
      exc.dt = strtod( ptr, endp );      

      if ( exc.dt < 1.e-277 )
      {
         fprintf( display, "\n\n Error message from function %s:",
	    __func__ );
         fprintf( display, "\n Illegal: Vanishing transition time "
            "exc.dt = %e < 1.e-277 !!!", exc.dt ); 
         strncpy(( state->fcterr ), __func__, SHS_SIZE );
         strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
         strncat(( state->errmsg ), "vanishing_transition_time_exc.dt_",
            ( LGS_SIZE - SHS_SIZE ));
         return null; /* abnormal return */
      };
/*............................................................................*/

      fscanf( excitfle, scformat, ptr );
      exc.ht = strtod( ptr, endp );       /* plateau time */
      fscanf( excitfle, scformat, ptr );
      exc.fr[null] = strtod( ptr, endp ); /* exciting frequency [Hz] */
         
      goto cont1;
   }
   else
   {
      strncpy ( exc.type, "___________________", SHS_SIZE );
      exc.lbl = null;

      fprintf( display, "\n\n Error message from function %s:",
	 __func__ ); 
      fprintf( display, "\n Unknown excitation type %s!!!\n ", type );
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":illegal_excitation_type_", 
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
/*............................................................................*/
  cont1:

   fscanf( excitfle, scformat, ptr ); /* string "excited_ports" */
   fscanf( excitfle, scformat, ptr ); /* string "number" */

   fscanf( excitfle, scformat, ptr ); /* string "E-ports" */
   fscanf( excitfle, scformat, ptr ); /* long integer string */
   exc.ne = strtol( ptr, endp, DEC );
      
   if ( EXCEP < exc.ne )
   {
      fprintf( display, "\n\n Error message from function %s:",
	 __func__ );
      fprintf( display, "\n\n Too many E ports excited "
         "in DSC mesh %s !!!", exc.name );
      fprintf( display, "\n [ Maximum number is %ld = macro EXCEP "
         "in SOLVER.CONF or %s.", ( long ) EXCEP, __func__ );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_E_ports_excited_", 
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( excitfle, scformat, ptr ); /* string "H-ports" */
   fscanf( excitfle, scformat, ptr ); /* long integer string */
   exc.nh = strtol( ptr, endp, DEC );

   if ( EXCHP < exc.nh )
   {
      fprintf( display, "\n\n Error message from function %s:",
	 __func__ );
      fprintf( display, "\n\n Too many H ports excited "
         "in DSC mesh %s !!!", exc.name );
      fprintf( display, "\n [ Maximum number is %ld = macro EXCHP "
         "in SOLVER.CONF or %s.", ( long ) EXCHP, __func__ );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_H_ports_excited_", 
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
/*............................................................................*/
# if DSC_HCRMDE != 0

   fscanf( excitfle, scformat, ptr ); /* string "thermal_excitation" */
   fscanf( excitfle, scformat, ptr ); /* string "TIME_DOMAIN" */
   fscanf( excitfle, scformat, ptr ); /* string "PASSIVE" [ or else ] */
   strncpy( type, ptr, SHS_SIZE );
   fscanf( excitfle, scformat, ptr ); /* string "parameters:" */

/* reset: */
   exc.hcl = null;
   exc.nhc = null;
   exc.ntf = null;
   exc.ntn = null;

   if ( null == strncmp( type, "PASSIVE____________", TWO ))
   {
      strncpy ( exc.hctp, "PASSIVE____________", SHS_SIZE );
      exc.hcl = 0;

      goto cont2;
   }
   else if ( null == strncmp( type, "FLOATING___________", TWO ))
   {
      strncpy ( exc.hctp, "FLOATING___________", SHS_SIZE );
      exc.hcl = 1;

      goto cont2;
   }
   else if ( null == strncmp( type, "STEADY_STATE_______", TWO ))
   {
      strncpy ( exc.hctp, "STEADY_STATE_______", SHS_SIZE );
      exc.hcl = 2;

      fscanf( excitfle, scformat, ptr );
      exc.hcn = strtol( ptr, endp, DEC ); /* smoothing order */
      fscanf( excitfle, scformat, ptr );
      exc.hcdt = strtod( ptr, endp ); /* smoothing time */

      goto cont2;
   }
   else if ( null == strncmp( type, "DIRAC_PULSE________", TWO ))
   {
      strncpy ( exc.hctp, "DIRAC_PULSE________", SHS_SIZE );
      exc.hcl = 3;

      goto cont2;
   }
   else if ( null == strncmp( type, "SMOOTH_PULSE______", TWO ))
   {
      strncpy ( exc.hctp, "SMOOTH_PULSE_______", SHS_SIZE );
      exc.hcl = 4;

      fscanf( excitfle, scformat, ptr );
      exc.hcdt = strtod( ptr, endp ); /* smoothing time */

      goto cont2;
   } 
   else if ( null == strncmp( type, "GAUSS_PULSE________", TWO ))
   {
      strncpy ( exc.hctp, "GAUSS_PULSE________", SHS_SIZE );
      exc.hcl = 5;

      fscanf( excitfle, scformat, ptr );
      exc.hcrt = strtod( ptr, endp ); /* sigma time */
      fscanf( excitfle, scformat, ptr );
      exc.hcdt = strtod( ptr, endp ); /* delay time */

      goto cont2;
   } 
   else if ( null == strncmp( type, "RAMP_PULSE_________", TWO ))
   {
      strncpy ( exc.hctp, "RAMP_PULSE_________", SHS_SIZE );
      exc.hcl = 6;

      fscanf( excitfle, scformat, ptr );
      exc.hcrt = strtod( ptr, endp ); /* rise time */

      goto cont2;
   } 
   else if ( null == strncmp( type, "HEAVISIDE_STEP_____", TWO ))
   {
      strncpy ( exc.hctp, "HEAVISIDE_STEP_____", SHS_SIZE );
      exc.hcl = 7;

      goto cont2;
   } 
   else if ( null == strncmp( type, "PLATEAU_SMOOTHED___", TWO ))
   {
      strncpy ( exc.hctp, "PLATEAU_SMOOTHED___", SHS_SIZE );
      exc.hcl = 8;

      fscanf( excitfle, scformat, ptr ); 
      exc.hcn = strtol( ptr, endp, DEC ); /* smoothing order */
      fscanf( excitfle, scformat, ptr );
      exc.hcht = strtod( ptr, endp );      /* hold time */
      fscanf( excitfle, scformat, ptr );
      exc.hcdt = strtod( ptr, endp );      /* decline time */

      goto cont2;
   }
   else if ( null == strncmp( type, "PERIODIC_PULSE_____", TWO ))
   {
      strncpy ( exc.hctp, "PERIODIC_PULSE_____", SHS_SIZE );
      exc.hcl = 9;

      fscanf( excitfle, scformat, ptr );
      exc.hcht = strtod( ptr, endp );      /* hold time */
      fscanf( excitfle, scformat, ptr );
      exc.hcdt = strtod( ptr, endp );      /* decline time */

      goto cont2;
   }
   else if ( null == strncmp( type, "DOUBLE_PERIODIC____", TWO ))
   {
      strncpy ( exc.hctp, "DOUBLE_PERIODIC____", SHS_SIZE );
      exc.hcl = 10;

      fscanf( excitfle, scformat, ptr );
      exc.hcht = strtod( ptr, endp );      /* 1st hold time */
      fscanf( excitfle, scformat, ptr );
      exc.hcdt = strtod( ptr, endp );      /* 1st decline time */

      fscanf( excitfle, scformat, ptr );
      exc.hcht2 = strtod( ptr, endp );     /* 2nd hold time */
      fscanf( excitfle, scformat, ptr );
      exc.hcdt2 = strtod( ptr, endp );     /* 2nd decline time */

      goto cont2;
   }
   else if ( null == strncmp( type, "SWITCHED_SOURCES___", TWO ))
   {
      strncpy ( exc.hctp, "SWITCHED_SOURCES___", SHS_SIZE );
      exc.hcl = 11;

      fscanf( excitfle, scformat, ptr );
      exc.hcht = strtod( ptr, endp );      /* 1st hold time */
      fscanf( excitfle, scformat, ptr );
      exc.hcdt = strtod( ptr, endp );      /* 1st decline time */

      fscanf( excitfle, scformat, ptr );
      exc.hcht2 = strtod( ptr, endp );     /* 2nd hold time */
      fscanf( excitfle, scformat, ptr );
      exc.hcdt2 = strtod( ptr, endp );     /* 2nd decline time */

      goto cont2;
   }
   else
   {
      strncpy (( exc.hctp ), ptr, SHS_SIZE );
      exc.hcl = null;

      fprintf( display, "\n\n Error message from function %s:",
         __func__ );
      fprintf( display, "\n Unknown thermal excitation type %s!!!\n ",
         type );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":illegal_therm_exc_type_", 
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
/*............................................................................*/
  cont2:

   fscanf( excitfle, scformat, ptr ); /* string "excited_ports".*/
   fscanf( excitfle, scformat, ptr ); /* string "number..." */

   fscanf( excitfle, scformat, ptr ); /* string "heat_currents" */
   fscanf( excitfle, scformat, ptr ); /* number of excited node temps */
   exc.nhc = strtol( ptr, endp, DEC );

   if ( EXCHC < exc.nhc )
   {
      fprintf( display, "\n\n Error message from function %s:",
         __func__ );
      fprintf( display, "\n\n Too many temperature nodes to be excited "
         "in DSC mesh %s !!!", exc.name );
      fprintf( display, "\n [ Maximum number is %ld = macro EXCHC "
         "in SOLVER.CONF or %s.", ( long ) EXCHC, __func__ );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_temp_nodes_excited_", 
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };

   fscanf( excitfle, scformat, ptr ); /* string "face_temperature" */
   fscanf( excitfle, scformat, ptr ); /* number of excited face temps */
   exc.ntf = strtol( ptr, endp, DEC );

   if ( EXCTF < exc.ntf )
   {
      fprintf( display, "\n\n Error message from function %s:",
         __func__ );
      fprintf( display, "\n\n Too many temperature faces to be excited "
         "in DSC mesh %s !!!", exc.name );
      fprintf( display, "\n [ Maximum number is %ld = macro EXCTF "
         "in SOLVER.CONF or %s.", ( long ) EXCTF, __func__ );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_temp_nodes_excited_", 
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };

   fscanf( excitfle, scformat, ptr ); /* string "node_temperature" */
   fscanf( excitfle, scformat, ptr ); /* number of excited temp nodes */
   exc.ntn = strtol( ptr, endp, DEC );

   if ( EXCTN < exc.ntn )
   {
      fprintf( display, "\n\n Error message from function %s:",
         __func__ );
      fprintf( display, "\n\n Too many temperature nodes to be excited "
         "in DSC mesh %s !!!", exc.name );
      fprintf( display, "\n [ Maximum number is %ld = macro EXCTN "
         "in SOLVER.CONF or %s.", ( long ) EXCTN, __func__ );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), ( state->exc ), SHS_SIZE );
      strncat(( state->errmsg ), ":too_many_temp_nodes_excited_", 
         ( LGS_SIZE - SHS_SIZE ));

      return null; /* abnormal return */
   };
# endif
/*............................................................................*/
/* the cell, face, and port labels, and exciation parameters are read, now: */

   fscanf( excitfle, scformat, ptr ); /* string ">>CELL&FACE_LABELS_etc..." */

   for ( ii=null; ii<exc.ne; ii++ )
   {
      fscanf( excitfle, scformat, ptr );
      exc.me[ii] = strtol( ptr, endp, DEC ); 
      fscanf( excitfle, scformat, ptr );
      exc.pe[ii] = strtol( ptr, endp, DEC );
      fscanf( excitfle, scformat, ptr );
      exc.er[ii] = strtod( ptr, endp );

      if ( pp == TWO )/* frequency_domain */
      {
         fscanf( excitfle, scformat, ptr );
         exc.ei[ii] = strtod( ptr, endp );
      }
      else
         exc.ei[ii] = ZERO;

      uu = exc.er[ii];
      vv = exc.ei[ii];

      uu = uu*uu + vv*vv;
/*............................................................................*/
# if EXC_MAXFLD == 1
      if( exc.mx < uu )
      {
         exc.mx = uu;
         strcpy( unit, "V" );
      };
# endif
/*............................................................................*/
      if ( exc.pe[ii] < null )
      {
         exc.pe[ii] *= ( -ONE );
         exc.er[ii] *= ( -1. );
         exc.ei[ii] *= ( -1. );
      };
   };

   for ( ii=null; ii<exc.nh; ii++ )
   {
      fscanf( excitfle, scformat, ptr );
      exc.mh[ii] = strtol( ptr, endp, DEC );
      fscanf( excitfle, scformat, ptr );
      exc.ph[ii] = strtol( ptr, endp, DEC );
      fscanf( excitfle, scformat, ptr );
      exc.hr[ii]  = strtod( ptr, endp );

      if ( pp == TWO ) /* frequency_domain */
      {
         fscanf( excitfle, scformat, ptr );
         exc.hi[ii] = strtod( ptr, endp );
      }
      else
         exc.hi[ii] = ZERO;

      uu = exc.hr[ii];
      vv = exc.hi[ii];

      uu = uu*uu + vv*vv;
/*............................................................................*/
# if EXC_MAXFLD == 1

      if( exc.mx < uu )
      {
         exc.mx = uu;
         strcpy( unit, "A" );
      };
# endif
/*............................................................................*/
      if ( exc.ph[ii] < null )
      {
         exc.ph[ii] *= ( -ONE );
         exc.hr[ii] *= ( -1. );
         exc.hi[ii] *= ( -1. );
      };

      switch( exc.ph[ii] )
      {
        case  3:
        case  4:
        case  7:
        case  8:
        case 11:
        case 12:
         exc.hr[ii] *= ( -1. );
         exc.hi[ii] *= ( -1. ); 
         break;

        default:
         break;
      };
   };

# if EXC_MAXFLD == 1
   exc.mx = sqrt( exc.mx );
# endif
/*............................................................................*/
# if DSC_HCRMDE != 0
	      
   for ( ii=null; ii<exc.nhc; ii++ )
   {
      fscanf( excitfle, scformat, ptr );
      exc.mhc[ii] = strtol( ptr, endp, DEC );
      fscanf( excitfle, scformat, ptr );
      exc.fhc[ii] = strtol( ptr, endp, DEC );
      fscanf( excitfle, scformat, ptr );
      exc.hc[ii]  = strtod( ptr, endp );
   };
   for ( ii=null; ii<exc.ntf; ii++ )
   {
      fscanf( excitfle, scformat, ptr );
      exc.mtf[ii] = strtol( ptr, endp, DEC );
      fscanf( excitfle, scformat, ptr );
      exc.ftf[ii] = strtol( ptr, endp, DEC );
      fscanf( excitfle, scformat, ptr );
      exc.tf[ii]  = strtod( ptr, endp );
   };
   for ( ii=null; ii<exc.ntn; ii++ )
   {
      fscanf( excitfle, scformat, ptr );
      exc.mtn[ii] = strtol( ptr, endp, DEC );
      fscanf( excitfle, scformat, ptr );
      exc.tn[ii]  = strtod( ptr, endp );
   };
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/

   fclose( excitfle );

# if EXC_MAXFLD == 1
   exc.mx = sqrt( exc.mx );
# endif

   if (( state->rmfle ) != ONE )
      fprintf( display, "\n entered: excitation file %s", ( state->exc ));
   else if (( state->rmfle ) == ONE )
   {
      for ( kk=(jj+ONE); kk<( state->nbrjbs ); kk++ )
      {
         if ( exclbl[kk] == exclbl[jj] )
         {
            fprintf( display, "\n entered: excitation file %s", ( state->exc ));
            goto disp;
         };
      };
      ind = remove( state->exc );

      if ( ind == null ) 
         fprintf( display, "\n entered and removed: "
            "excitation file %s", ( state->exc ));
   };

  disp:   

# if EXC_DISP == 1

   fprintf( display, "\n\n text ..................................: " );

   ind = strlen( exc.text );
   if ( 35 < ind )  
      fprintf( display, "\n %s\n", exc.text );
   else
      fprintf( display, "%s", exc.text );

   fprintf( display, "\n EM field excitation type ..............: %s ",
      exc.type );

   switch ( exc.lbl )
   {
     case 0: /* ZERO_EM-FIELD______ */
      break;

     default: /* e.g. STEADY_STATE____ */

      if( ZERO < exc.dt )
      {
         fprintf( display, "\n smoothing order .......................: "
            "%ld ", ( long ) exc.nn );

         if (( state->dmn ) == 'f' )
         {
            fprintf( display, "\n number of smoothing cycles ............: "
               "%ld ", ( long ) exc.dt );
# if EXC_MAXFLD != 0                       
            fprintf( display, "\n excitation amplitude [ abs. maximum ]..: "
               "%e %s", exc.mx, unit );
# endif
         }
         else
            fprintf( display, "\n smoothing time ........................: "
               "%.12e s", exc.dt );
      }; 
      break;

     case 2: /* DIRAC_PULSE________ */
      break;

     case 3:
     case 9:
     case 10:
      fprintf( display, "\n exciting frequency ....................: "
	 "%.12e Hz", exc.fr[null] );
      break; 

     case 4:
      fprintf( display, "\n exciting frequency ....................: "
         "%.12e Hz", exc.fr[null] );
      fprintf( display, "\n transition time .......................: "
         "%.12e s", exc.dt );
      break;

     case 5:
      fprintf( display, "\n " );
      for ( ii=null; ii<exc.nn; ii++ )
      {
         fprintf( display, "\n %2d. exciting frequency ................: "
                           " %.15e Hz", ii+ONE, exc.fr[ii] );
      };
      break;

     case 6:
      fprintf( display, "\n sigma time ............................: "
         "%.12e s", exc.rt );
      fprintf( display, "\n delay .................................: "
         "%.12e s", exc.dt );
      break;

     case 7:
     case 14:
      fprintf( display, "\n exciting frequency ....................: "
         "%.12e Hz", exc.fr[null] );
      fprintf( display, "\n sigma time ............................: "
         "%.12e s", exc.rt );
      fprintf( display, "\n delay .................................: "
         "%.12e s", exc.dt );
      break;

     case 8:
      fprintf( display, "\n rise time .............................: "
	 "%.12e s", exc.rt );
      break;
   };
/*............................................................................*/
# if DSC_HCRMDE != 0

   fprintf( display, "\n thermal excitation type ...............: %s ",
      exc.hctp );

   switch ( exc.hcl )
   {
     case 0: /* PASSIVE____________ */
      break;

     case 1: /* FLOATING___________ */
      break;

     case 2: /* STEADY_STATE_______ */

      if( ZERO < exc.hcdt )
      {
         fprintf( display, "\n smoothing order .......................: "
            "%ld", ( long ) exc.hcn );

         fprintf( display, "\n smoothing time ........................: "
            "%.12e s", exc.hcdt );
      };
      break;

     case 3: /* DIRAC_PULSE________ */
      break;

     case 4: /* SMOOTH_PULSE_______ */
      fprintf( display, "\n smoothing order .......................: "
         "%d", exc.hcn );
      fprintf( display, "\n smoothing time ........................: "
         "%.12e s", exc.hcdt );
      break;

     case 5: /* GAUSS_PULSE________ */
      fprintf( display, "\n sigma time ............................: "
	 "%.12e s", exc.hcrt );
      fprintf( display, "\n delay time ............................: "
	 "%.12e s", exc.hcdt );
      break;

     case 6: /* RAMP_PULSE_________ */
      fprintf( display, "\n rise time .............................: "
         "%.12e s", exc.hcrt );
      break;

     case 7: /* HEAVISIDE_STEP_____ */
      break;

     case 8: /* PACKET_SMOOTHED____ */
      fprintf( display, "\n smoothing order .......................: "
         "%d", exc.hcn );
      fprintf( display, "\n hold time .............................: "
         "%.12e s", exc.hcht );
      fprintf( display, "\n decline time ..........................: "
         "%.12e s", exc.hcdt );
      break;

     case 9: /* PERIODIC_PULSE_____ */
      fprintf( display, "\n Hold time .............................: "
         "%.12e s", exc.hcht );
      fprintf( display, "\n Decline time ..........................: "
         "%.12e s", exc.hcdt );
      break;

     case 10: /* PERIODIC_PULSE_____ */
      fprintf( display, "\n 1st hold time .........................: "
         "%.12e s", exc.hcht );
      fprintf( display, "\n 1st decline time ......................: "
         "%.12e s", exc.hcdt );
      fprintf( display, "\n 2nd hold time .........................: "
         "%.12e s", exc.hcht2 );
      fprintf( display, "\n 2nd decline time ......................: "
         "%.12e s", exc.hcdt2 );
      break;

     case 11: /* SWITCHED_SOURCES___ */
      fprintf( display, "\n 1st hold time .........................: "
         "%.12e s", exc.hcht );
      fprintf( display, "\n 1st decline time ......................: "
         "%.12e s", exc.hcdt );
      fprintf( display, "\n 2nd hold time .........................: "
         "%.12e s", exc.hcht2 );
      fprintf( display, "\n 2nd decline time ......................: "
         "%.12e s", exc.hcdt2 );
      break;

     default:
      break;
   };
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* E_field: */

   if ( null < exc.ne )
   {  
      fprintf( display, "\n\n E-field excitation at %05ld ports: >>---"
         "------------------------------------->", ( long ) exc.ne );
      llns = exc.ne/clns + ONE;

/* E_cells: */

      for ( kk=null; kk<llns; kk++ )
      {
         fprintf( display, "\n -> cell:" ); 
         pp = kk*clns;

         for ( ii=ONE; ii<=clns; ii++ )
         {
            if (( lbnd < kk )&&( cbnd < ii ))
            {
               fprintf( display, "   ***%c", bslsh );
               fprintf( display, "%6ld%c", exc.me[exc.ne-ONE], bslsh );
               goto E_ports;
            };    

            if ( pp < exc.ne )   
               fprintf( display, "%6ld%c", exc.me[pp], bslsh );  

            pp++;

            if ( exc.ne <= pp )  
               goto E_ports;
         };

        E_ports:   

         fprintf( display, "\n -> port:" );
         pp = kk*clns;

         for ( ii=ONE; ii<=clns; ii++ ) 
         { 
            if (( lbnd < kk )&&( cbnd < ii ))
            {
               fprintf( display, "   ***%c", slsh );
               fprintf( display, "%6ld%c", (long) exc.pe[exc.ne-ONE], slsh );
               goto H_field;
            };
            if ( pp < exc.ne )   
               fprintf( display, "%6ld%c", (long) exc.pe[pp], slsh );  

            pp++;

            if ( exc.ne <= pp )  
              goto H_field; 
         };
      };            
   };

  H_field:

   if ( null < exc.nh )
   {
      fprintf( display, "\n\n H-field excitation at %05ld ports: >>"
         "---------------------------------------->", ( long ) exc.nh );
      llns = exc.nh/clns + ONE;

/* H_cells: */

      for ( kk=null; kk<llns; kk++ )
      {
         fprintf( display, "\n -> cell:" );
         pp = kk*clns;

         for ( ii=ONE; ii<=clns; ii++ )
         {
            if (( lbnd < kk )&&( cbnd < ii ))
            {
               fprintf( display, "   ***%c", bslsh );
               fprintf( display, "%6ld%c", exc.mh[exc.nh-ONE], bslsh );
               goto H_ports;
            };

            if ( pp < exc.nh )   
               fprintf( display, "%6ld%c", exc.mh[pp], bslsh ); 

            pp++;

            if ( exc.nh <= pp )  
               goto H_ports;  
         };

        H_ports: 

         fprintf( display, "\n -> port:" );
         pp = kk*clns;

         for ( ii=ONE; ii<=clns; ii++ )
         {
            if (( lbnd < kk )&&( cbnd < ii ))
            {
               fprintf( display, "   ***%c", slsh );
               fprintf( display, "%6ld%c", (long) exc.ph[exc.nh-ONE], slsh );
               goto listo;
            };
            if ( pp < exc.nh )   
               fprintf( display, "%6ld%c", (long) exc.ph[pp], slsh ); 

            pp++;

            if ( exc.nh <= pp )  
               goto listo;
         };
      };   
   }; 

  listo: ;

   fprintf( display, "\n --------------------------------------"
      "----------------------------------------");
# endif

   return ONE; /* normal return */
}
/*============================================================================*/
# undef EXC_DISP
/*********************** end of function excite(*) ****************************/



# ifdef _Include
/*----------------------------------------------------------------------------*/
/* [ file: excfld.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function excfld(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   DSC EMfield excitation function                                            *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# ifdef OPTIMIZE
   # pragma OPT_LEVEL 3
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
/* structure excitation exc: cf. function scattr.h                            */
/*----------------------------------------------------------------------------*/
/* There are different ways to interlace excited with incident quantities:    */
     
# ifndef EXC_IMPOSE
   # define EXC_IMPOSE 0 /* 0: additive superposition [ recommended ]         */
# endif                  /* 1: direct imposition                              */
                         /* 2: imposition with scattered quantities subtracted*/
                         /*    [ yields fixed total quantities ].             */
# ifndef EXC_MAXFLD
   # define EXC_MAXFLD 1 /* 1: write maximum excited field into ( ept->mx )   */
# endif
/*----------------------------------------------------------------------------*/
/* constants, stabilizing ONE [ to be modified with care ]:                   */

# ifndef PRECISION
   # define PRECISION ( 1.e-15 )
# endif

# if EXC_IMPOSE != 0
   # define ONE_MINUS ( 1. - 33.*PRECISION ) 
# endif
/*============================================================================*/

DSC_FIELDS *\
excfld( void )
{
/* allusions: */
/*
   extern struct solverstat solver;
   extern struct excitation exc;

# if DSC_DOMAIN == 1
   extern DSC_FIELDS fld;
# else
   extern DSC_FIELDS fld[];
# endif
*/
/* declarations: */

   static struct solverstat
     *state = &solver;

   static long
      mm = null;

   static short 
      ii = null;

   static signed char
      pp = null;

   static struct
      excitation *ept;

   static DSC_FIELDS
      *ret;

/* prototypes: */

   char 
     *lotos ( long lngint, char length );
/*----------------------------------------------------------------------------*/
   if(( state->dmn ) == 't' )
   {
/*............................................................................*/
# if DSC_DOMAIN == 2
      fprintf( display, "\n\n Error message from function %s:",
         __func__ );
      fprintf( display, "\n Solver is compiled for frequency "
         "domain operation 'DSC_DOMAIN = 2'" );
      fprintf( display, "\n - yet run in time domain mode 1 !!!" );
      fprintf( display, "\n [ Change time/frequency domain macro "
         "DSC_DOMAIN in SOLVER.CONF from 2 to" );
      fprintf( display, "\n   1: time domain, or          "
         "\n   0: time & frequency domain. " );
      fprintf( display, "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
         "in_mode_DSC_DOMAIN=2] ", LGS_SIZE );
      return null; /* abnormal return */
/*............................................................................*/
# else /* if DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
      ept = &exc;
/*............................................................................*/
# if DSC_DOMAIN == 0
      ret = &fld[null];
# else
      ret = &fld;
# endif
/*............................................................................*/
      ii = null;
      while( ii < ( ept->ne ))
      {
         mm = ( ept->me[ii] );
         pp = ( ept->pe[ii] ) - ONE;
/*............................................................................*/
# if EXC_IMPOSE == 1
         ( ret->i[mm][pp] ) = ZERO;
# elif EXC_IMPOSE == 2
         ( ret->i[mm][pp] ) = - ONE_MINUS*( ret->r[mm][pp] );
# endif
/*............................................................................*/
         ( ret->i[mm][pp] ) += ( state->swing )*( ept->er[ii] );

         ii++ ;
      };

      ii = null; 
      while( ii < ( ept->nh ))
      {
         mm = ( ept->mh[ii] );
         pp = ( ept->ph[ii] ) - ONE;
/*............................................................................*/
# if EXC_IMPOSE == 1
         ( ret->i[mm][pp] ) = ZERO;
# elif EXC_IMPOSE == 2
         ( ret->i[mm][pp] ) = ONE_MINUS*( ret->r[mm][pp] );
# endif
/*............................................................................*/
         ( ret->i[mm][pp] ) += ( state->swing )*( ept->hr[ii] );

         ii++ ; 
      };

      return ret; /* normal return */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
   }
   else /* if (( state->dmn ) != 't' ): frequency domain */
   {
/*............................................................................*/
# if DSC_DOMAIN == 1
      fprintf( display, "\n\n Error message from function %s:", __func__);
      fprintf( display, "\n Solver is compiled for time domain "
         "operation, 'DSC_DOMAIN = 1'" );
      fprintf( display, "\n - yet run in frequency domain mode 2 !!!" );
      fprintf( display, "\n [ Change time/frequency domain macro "
         "DSC_DOMAIN in SOLVER.CONF from 1 to" );
      fprintf( display, "\n   2: frequency domain, or     "
         "\n   0: time & frequency domain. " );
      fprintf( display, "\n   Then re-compile and restart program. ]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
         "in_mode_DSC_DOMAIN=1] ", LGS_SIZE );
      return NULL;
/*............................................................................*/
# else /* DSC_DOMAIN == 0 or 2 */

      ept = &exc;
      ret = &fld[null];

      ii = null;
      while( ii < ( ept->ne ))
      {
         mm = ( ept->me[ii] );
         pp = ( ept->pe[ii] ) - ONE;
/*............................................................................*/
# if EXC_IMPOSE == 1
         ( ret->r[mm][pp] ) = ZERO;
         ( ret->i[mm][pp] ) = ZERO;
# elif EXC_IMPOSE == 2
         ( ret->r[mm][pp] ) = - ONE_MINUS*( ret->r[mm][pp] );
         ( ret->i[mm][pp] ) = - ONE_MINUS*( ret->i[mm][pp] );
#endif
/*............................................................................*/
         ( ret->r[mm][pp] ) += ( state->swing )*( ept->er[ii] );
         ( ret->i[mm][pp] ) += ( state->swing )*( ept->ei[ii] );

         ii++ ;
      };

      ii = null;
      while( ii < ( ept->nh ))
      {
         mm = ( ept->mh[ii] );
         pp = ( ept->ph[ii] ) - ONE;
/*............................................................................*/
# if EXC_IMPOSE == 1
         ( ret->r[mm][pp] ) = ZERO;
         ( ret->i[mm][pp] ) = ZERO;
# elif EXC_IMPOSE == 2
         ( ret->r[mm][pp] ) = ONE_MINUS*( ret->r[mm][pp] );
         ( ret->i[mm][pp] ) = ONE_MINUS*( ret->i[mm][pp] );
# endif
/*............................................................................*/
         ( ret->r[mm][pp] ) += ( state->swing )*( ept->hr[ii] );
         ( ret->i[mm][pp] ) += ( state->swing )*( ept->hi[ii] );

         ii++ ;
      };

      return ret; /* normal return */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
   }; /* end if (( state->dmn ) != 't' ), frequency domain */
} 
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
# if EXC_IMPOSE != 0
   # undef ONE_MINUS
# endif
# undef EXC_IMPOSE
/************** end of DSC system excitation function excfld(*) ***************/



/*----------------------------------------------------------------------------*/
/* [ file: exchcr.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function exchcr(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   DSC temperature distribution and heat current excitation function          *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPT_LEVEL 3
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
/*
# include "../math/consts.h"
*/
/*----------------------------------------------------------------------------*/
/* structure excitation exc: cf. function scattr.h                            */
/*----------------------------------------------------------------------------*/
/* There are different ways to interlace excited with incident quantities:    */
/*----------------------------------------------------------------------------*/
/* constants, stabilizing ONE [ to be modified with care ]:                   */

# ifndef PRECISION
   # define PRECISION   1.e-15
# endif
/*============================================================================*/

DSC_HCRRTS * \
exchcr( const short cc )
{
/* allusions: */
/*
   extern struct solverstat solver;
   extern struct excitation exc;
*/
/* declarations: */

   static struct solverstat
     *state = &solver;

   static long
      mm = null;

   static short 
      fc = null,
      ii = null;

   static struct
     excitation *ept;

   static DSC_HCRRTS
      *hre;

/* prototypes: - */
/*----------------------------------------------------------------------------*/
   if(( state->dmn ) == 't' )
   {
/*............................................................................*/
# if DSC_DOMAIN == 2
      fprintf( display, "\n\n Error message from function %s:", __func__ );
      fprintf( display, "\n Solver is compiled for frequency "
         "domain operation 'DSC_DOMAIN = 2'" );
      fprintf( display, "\n - yet run in time domain mode 1 !!!" );
      fprintf( display, "\n [ Change time/frequency domain macro "
         "DSC_DOMAIN in SOLVER.CONF from 2 to" );
      fprintf( display, "\n   1: time domain, or          "
         "\n   0: time & frequency domain. " );
      fprintf( display, "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
                            "in_mode_DSC_DOMAIN=2] ", LGS_SIZE );
      return NULL; /* abnormal return */
/*............................................................................*/
# else /* if DSC_DOMAIN == 0 or 1 */

      ept = &exc;
      hre = &hcr[cc];

      ii = null;
      while( ii < ( ept->nhc ))
      {
         mm = ( ept->mhc[ii] );
         fc = ( ept->fhc[ii] );
         ( hre->ic[mm][fc] ) = ( state->hcswg )*( ept->hc[ii] );
         ii++ ;
      };

      ii = null;
      while( ii < ( ept->ntf ))
      {
         mm = ( ept->mtf[ii] );
         fc = ( ept->ftf[ii] );
         ( hre->tf[mm][fc] ) = ( state->hcswg )*( ept->tf[ii] );

         ii++ ;
      };

      ii = null;
      while( ii < ( ept->ntn ))
      {
         mm = ( ept->mtn[ii] );
         ( hre->tn[mm] ) = ( state->hcswg )*( ept->tn[ii] );

         ii++ ;
      };

      return hre; /* normal return */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0, 1 */
/*............................................................................*/
   }
   else /* if (( state->dmn ) != 't' ): frequency domain */
   {
/*............................................................................*/
# if DSC_DOMAIN == 1
      fprintf( display, "\n\n Error message from function %s:", __func__ );
      fprintf( display, "\n Solver is compiled for time domain "
         "operation, 'DSC_DOMAIN = 1'" );
      fprintf( display, "\n - yet run in frequency domain mode 2 !!!" );
      fprintf( display, "\n [ Change time/frequency domain macro "
         "DSC_DOMAIN in SOLVER.CONF from 1 to" );
      fprintf( display, "\n   2: frequency domain, or     "
         "\n   0: time & frequency domain." );
      fprintf( display, "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
         "in_mode_DSC_DOMAIN=1] ", LGS_SIZE );

      return NULL; /* abnormal return */
/*............................................................................*/
# else /* if DSC_DOMAIN == 0 or 2 */

      ept = &exc;
      hre = &hcr[cc];

      ii = null;
      while( ii < ( ept->nhc ))
      {
         mm = ( ept->mhc[ii] );
         fc = ( ept->fhc[ii] );
         ( hre->ic[mm][fc] ) = ( state->hcswg )*( ept->hc[ii] );
         ii++ ;
      };

      ii = null;
      while( ii < ( ept->ntf ))
      {
         mm = ( ept->mtf[ii] );
         fc = ( ept->ftf[ii] );
         ( hre->tf[mm][fc] ) = ( state->hcswg )*( ept->tf[ii] );

         ii++ ;
      };

      ii = null;
      while( ii < ( ept->ntn ))
      {
         mm = ( ept->mtn[ii] );
         ( hre->tn[mm] ) = ( state->hcswg )*( ept->tn[ii] );

         ii++ ;
      };

      return hre; /* normal return */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
   }; /* end if (( state->dmn ) != 't' ), frequency domain */
}
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/******** end of DSC heat current and temp distrib function exchcr(*) *********/
# endif /* DSC_HCRMDE != 0 */
# endif /* _Include defined */
/**************** end of DSC system excitation file excite.h ******************/
