/* [ file: trnsnt.h ] */
# define DO_TRNFLD "trnfld(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function trnfld(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Controls time dependence of the DSC system excitation, as specified        *
*   in excitation function 'excite.h' [ DIRAC pulse, HARMONIC, ... ] .         *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"
# include "../math/consts.h"
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPT_LEVEL 3
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
# define TRN_MARGIN  9.980e-01 
# define TRN_SMOOTH  9.995e-01
# define TRN_BOUND   3.000e+01
/*============================================================================*/

short trnfld( const double ttime )
{
/* allusions: */
/*
   extern struct excitation exc;
*/
/* declarations: */

   static double 
      xx = ZERO,
      yy = ZERO,
      tt = ZERO;

   static int 
      ii = null;

   double fabs( double x );
   double cos( double x );
   double sin( double x );
   double exp( double x );

   void srand( unsigned int seed ); 
   int rand( void );  

/*----------------------------------------------------------------------------*/

   switch ( exc.lbl )
   {
     case 0: /* type: ZERO_EM-FIELD______ */
      return null;
      break;

     case 1: /* type: STEADY_STATE_______ */

      if( exc.dt <= ttime )
      {
         solver.swing = 1.;
         return ONE;
      }
      else if ( exc.dt/2. < ttime )
      {
         tt = 2.* (( exc.dt - ttime ) / exc.dt );

         xx = 1.;
         ii = null; do   /* tt^exc.nn */
         {
            xx *= tt;
            ii++;
         }  while( ii < exc.nn );

         solver.swing = 1. - .5*xx;
         return ONE;
      }
      else 
      {
         tt = 2. * ttime / exc.dt;

         xx = 1.;
         ii = null; do   /* tt^exc.nn */
         {
            xx *= tt;
            ii++;
         }  while( ii < exc.nn );

         solver.swing = .5*xx;
         return ONE;
      };
      break;

     case 2: /* type: DIRAC_PULSE________ */
    
      if ( TRN_MARGIN*solver.dt < ttime )
      {
         return null;
      } 
      else 
      {
         solver.swing = 1./solver.dt;
         return ONE;
      };
      break;

     case 3: /* type: HARMONIC_SINUSOIDAL */

      solver.swing = sin( 2.*PI*ttime*exc.fr[null] );  
      return ONE;
      break;
   
     case 4: /* type: SMOOTH_HARMONIC____ */

      if ( exc.dt <= ttime )
      {
         tt = ttime - exc.dt;
         solver.swing = cos( 2.*PI*tt*exc.fr[null] );
         return ONE;
      }
      else if ( exc.dt/2. < ttime )
      {
         tt = ttime - exc.dt;
         xx = 1. - .5*exc.dt/( - tt );
         solver.swing = ( 1. - .5*exp( xx ))*cos( 2.*PI*tt*exc.fr[null] );
         return ONE;
      }
      else if ( ZERO < ttime )
      {
         tt = ttime - exc.dt;
         xx = 1. - .5*exc.dt/ttime ;
         solver.swing = .5*exp( xx )*cos( 2.*PI*tt*exc.fr[null] );
         return ONE;
      }
      else
      {
         return null;
      };
      break;

     case 5: /* type: MULTIPLE_HARMONIC__ */

      solver.swing = ZERO;
      for ( ii=null ; ii<exc.nn ; ii++ )
      {
         solver.swing += sin( 2.*PI*ttime*exc.fr[ii] );
      };
      return ONE;
      break;

     case 6: /* type: GAUSS_PULSE________ */

      tt = ttime - exc.dt;
      xx = tt/exc.rt;
      xx *= xx;
      if ( TRN_BOUND < xx )
      {
         return null;
      }
      else
      {
         solver.swing = exp( -.5*xx );
         return ONE;
      };
      break;

     case 7: /* type: OSCILLATORY_GAUSS__ */

      tt = ttime - exc.dt;
      xx = tt/exc.rt;
      xx *= xx;
      if ( TRN_BOUND < xx )
      {
         return null;
      }
      else
      {
         solver.swing = exp( -.5*xx );
         solver.swing *= cos( 2.*PI*tt*exc.fr[null] );
         return ONE;
      }; 
      break;

     case 8: /* type: RAMP_PULSE_________ */

      if ( exc.rt < ttime )
      {
         solver.swing = 1.;
         return ONE;
      }
      else
      {
         solver.swing = ttime/exc.rt ;
         return ONE;
      };
      break;

     case 9: /* type: RECTANGULAR_PERIODC */

      xx = cos( 2.*PI*ttime*exc.fr[null] );

      if ( ZERO <= xx )  
      { 
         solver.swing = 1.;
         return ONE;
      } 
      else              
      {
         solver.swing = -1.;
         return ONE;
      }; 
      break;

     case 10: /* type: SAW_TOOTH_PERIODIC_ */

      solver.swing = fmod( ttime, exc.rt )/exc.rt - .5;
      return ONE;
      break;
    
     case 11: /* type: HEAVISIDE_STEP_____ */
   
      solver.swing = 1.;
      return ONE;
      break;

     case 12: /* type: DSC_TIMESTEP_PRDC__ */

      if ( TRN_MARGIN*solver.dt < ttime )
      {
         solver.swing *= ( -1. );
         return ONE;
      }
      else
      {
         solver.swing = 1./solver.dt;
         return ONE;
      };
      break;

/*   case 13: */      /* type: NOISE_AT_RANDOM____ [ random times ] */
/*                                                 
      if ( TRN_MARGIN*solver.dt < ttime )
      {
         xx = ( double ) rand( )/RAND_MAX;

         if  ( TRN_SMOOTH < fabs(xx) ) 
            solver.swing *= ( -1. );

         return ONE;
      }
      else
      {
         srand( ONE );
         solver.swing = 1.;
         return ONE;
      };
      break;
*/  
      case 14: /* type: MORLET_WAVELET_____ */

      tt = ttime - exc.dt;
      xx = tt/exc.rt;
      xx *= xx;
      if ( TRN_BOUND < xx )
      {
         return null;
      }
      else
      {
         yy = exc.fr[null];
         solver.swing = cos( 2.*PI*yy*tt );
         solver.swing -= exp( -.5*yy*yy );
         solver.swing *= exp( -.5*xx );
         return ONE;
      };
      break;

     case 15: /* type: WAVE_PACKET_HARMNC_ */

      if ( exc.ht + 2.*exc.dt < ttime )
      {
         return null;
      }
      else if ( exc.ht + exc.dt < ttime )
      {
         tt = 2.*( ttime - exc.ht - exc.dt )/exc.dt;
         if ( 1. < tt )
         {
            xx = 2. - tt;
            yy = xx; 
            ii = ONE; do   /* xx^exc.nn */
            {
               xx *= yy;
               ii++;
            }  while( ii < exc.nn );
         }
         else if ( tt <= 1. ) 
         {
            yy = tt;
            ii = ONE; do   /* yy = tt^exc.nn */
            {
               yy *= tt; 
               ii++;
            }  while( ii < exc.nn );
            xx = 2. - yy;
         };
         xx /= 2.;
         solver.swing = xx * cos( 2.*PI*ttime*exc.fr[null] );
         return ONE;
      }
      else if ( exc.dt < ttime )
      {
         solver.swing = cos( 2.*PI*ttime*exc.fr[null] );
         return ONE;
      }
      else if ( 0. < ttime )
      {
         tt = 2.*ttime/exc.dt;

         if ( 1. < tt )
         {
            xx = 2. - tt;
            yy = xx;
            ii = ONE; do   /* yy = xx^exc.nn */
            {
               yy *= xx;
               ii++;
            }  while( ii < exc.nn );
            xx = 2. - yy;
         }
         else if ( tt <= 1. )
         {
            yy = tt;
            ii = ONE; do   /* yy = tt^exc.nn */
            {
               yy *= tt;
               ii++;
            }  while( ii < exc.nn );
            xx = yy;
         };
         xx /= 2.;
         solver.swing = xx * cos( 2.*PI*ttime*exc.fr[null] );
         return ONE;
      }
      else
         return null;
      break;

     default: /* Unknown excitation type */

      fprintf( display, "\n Message from function %s : ", __func__ );
      fprintf( display, "\n Unknown excitation type !!!\n " );

      exit( EXIT_FAILURE );

      break;
   };
}
/*============================================================================*/
/************************* end of function trnfld(*) **************************/

/*----------------------------------------------------------------------------*/
/* [ file: trnhcr.h ] */
# define DO_TRNCRR "trnhcr(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function trnhcr(*), DANSE release 1.0.                              *
*   [ DANSE program package ]                                                  *
*                                                                              *
*   Controls time dependence of DSC system thermic excitation, as specified    *
*   in function file excite.h.                                                 *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
/*
# include "../math/consts.h"
*/
# ifndef TRN_MARGIN
   # define TRN_MARGIN  9.980e-01
# endif
# ifndef TRN_SMOOTH
   # define TRN_SMOOTH  9.995e-01
# endif
# ifndef TRN_BOUND
   # define TRN_BOUND   3.000e+01
# endif
/*============================================================================*/

short trnhcr( const double hctme )
{
/* allusions: */
/*
   extern struct excitation exc;
*/
/* declarations: */

   static double 
      xx = ZERO,
      yy = ZERO,
      tt = ZERO,
      tt2 = ZERO;

   static int 
      ii = null;

   double fabs( double x );
   double cos( double x );
   double sin( double x );
   double exp( double x );

   void srand( unsigned int seed ); 
   int rand( void );  
/*----------------------------------------------------------------------------*/
   switch ( exc.hcl )
   {
     case 0: /* type: PASSIVE____________ */
     case 1: /* type: FLOATING___________ */

      solver.sws = ONE; /* enable sources [switch ON] */
      return null;
      break;

     case 2: /* type: STEADY_STATE_______ */

      solver.sws = ONE; /* enable sources [switch ON] */

      if( exc.hcdt <= hctme )
      {
         solver.hcswg = 1.;
         return ONE;
      }
      else if ( exc.hcdt/2. < hctme )
      {
         tt = 2.* (( exc.hcdt - hctme ) / exc.hcdt );

         xx = 1.;
         ii = null; do   /* tt^exc.hcn */
         {
            xx *= tt;
            ii++;
         }  while( ii < exc.hcn );

         solver.hcswg = 1. - .5*xx;
         return ONE;
      }
      else 
      {
         tt = 2. * hctme / exc.hcdt;

         xx = 1.;
         ii = null; do   /* tt^exc.hcn */
         {
            xx *= tt;
            ii++;
         }  while( ii < exc.hcn );

         solver.hcswg = .5*xx;
         return ONE;
      };
      break;

     case 3: /* type: DIRAC_PULSE________ */

      solver.sws = ONE; /* enable sources [switch ON] */

      if ( TRN_MARGIN*solver.hcdt < hctme )
      {
         return null;
      }
      else
      {
         solver.hcswg = 1.;
         return ONE;
      };
      break;

     case 4: /* type: SMOOTH_PULSE_______ */

      solver.sws = ONE; /* enable sources [switch ON] */

      if ( exc.hcdt <= hctme )
      {
         tt = hctme - exc.hcdt;
         solver.hcswg = 1.;
         return ONE;
      }
      else if ( exc.hcdt/2. < hctme )
      {
         tt = hctme - exc.hcdt;
         xx = 1. - .5*exc.hcdt/( - tt );
         solver.hcswg = ( 1. - .5*exp( xx ));
         return ONE;
      }
      else if ( ZERO < hctme )
      {
         tt = hctme - exc.hcdt;
         xx = 1. - .5*exc.hcdt/hctme;
         solver.hcswg = .5*exp( xx );
         return ONE;
      }
      else
      {
         return null;
      };
      break;

     case 5: /* type: GAUSS_PULSE________ */

      solver.sws = ONE; /* enable sources [switch ON] */

      tt = hctme - exc.hcdt;
      xx = tt/exc.hcrt;
      xx *= xx;
      if ( TRN_BOUND < xx )
      {
         return null;
      }
      else
      {
         solver.hcswg = exp( -.5*xx );
         return ONE;
      };
      break;

     case 6: /* type: RAMP_PULSE_________ */

      solver.sws = ONE; /* enable sources [switch ON] */

      if ( exc.hcrt < hctme )
      {
         solver.hcswg = 1.;
         return ONE;
      }
      else
      {
         solver.hcswg = hctme/exc.hcrt ;
         return ONE;
      };
      break;

     case 7: /* type: HEAVISIDE_STEP_____ */
   
      solver.sws = ONE; /* enable sources [switch ON] */

      solver.hcswg = 1.;
      return ONE;
      break;

     case 8: /* type: PLATEAU_SMOOTHED___ */

      solver.sws = ONE; /* enable sources [switch ON] */

      if ( exc.hcht + 2.*exc.hcdt < hctme )
      {
         return null;
      }
      else if ( exc.hcht + exc.hcdt < hctme )
      {
         tt = 2.*( hctme - exc.hcht - exc.hcdt )/exc.hcdt;
         if ( 1. < tt )
         {
            xx = 2. - tt;
            yy = xx; 
            ii = ONE; do   /* xx^exc.hcn */
            {
               xx *= yy;
               ii++;
            }  while( ii < exc.hcn );
         }
         else if ( tt <= 1. ) 
         {
            yy = tt;
            ii = ONE; do   /* yy = tt^exc.hcn */
            {
               yy *= tt; 
               ii++;
            }  while( ii < exc.hcn );
            xx = 2. - yy;
         };
         xx /= 2.;
         solver.hcswg = xx;
         return ONE;
      }
      else if ( exc.hcdt < hctme )
      {
         solver.hcswg = 1.;
         return ONE;
      }
      else if ( 0. < hctme )
      {
         tt = 2.*hctme/exc.hcdt;

         if ( 1. < tt )
         {
            xx = 2. - tt;
            yy = xx;
            ii = ONE; do   /* yy = xx^exc.hcn */
            {
               yy *= xx;
               ii++;
            }  while( ii < exc.hcn );
            xx = 2. - yy;
         }
         else if ( tt <= 1. )
         {
            yy = tt;
            ii = ONE; do   /* yy = tt^exc.hcn */
            {
               yy *= tt;
               ii++;
            }  while( ii < exc.hcn );
            xx = yy;
         };
         xx /= 2.;
         solver.hcswg = xx;
         return ONE;
      }
      else
      {
         return null;
      };
      break;

     case 9: /* type: PERIODIC_PULSE___ */

      solver.sws = ONE; /* enable sources [switch ON] */

      if (( exc.hcht + exc.hcdt ) <= ( hctme-tt ))
      {
         tt += ( exc.hcht + exc.hcdt ); /* one more period */
         solver.hcswg = 1.;
         return ONE;
      }
      else if ( exc.hcht <= ( hctme-tt ))
      {
         solver.hcswg = 0.;
         return ONE;
      }
      else if ( hctme < exc.hcht )
      {
         tt = ZERO; /* initialize period */
         solver.hcswg = 1.;
         return ONE;
      }
      else
      {
         solver.hcswg = 1.;
         return ONE;
      };
      break;

     case 10: /* type: DOUBLE_PERIODIC____ */

      if (( exc.hcht + exc.hcdt ) <= ( hctme-tt ))
      {
         tt += ( exc.hcht + exc.hcdt ); /* one more period */
         solver.hcswg = 1.;
         solver.sws = ONE; /* sources ON */
      }
      else if ( exc.hcht <= ( hctme-tt ))
      {
         solver.hcswg = 0.;
         solver.sws = null; /* sources OFF */
      }
      else if ( hctme < exc.hcht )
      {
         tt = ZERO; /* initialize period */
         solver.hcswg = 1.;
         solver.sws = ONE; /* sources ON */
      }
      else
      {
         solver.hcswg = 1.;
         solver.sws = ONE; /* sources ON */
      };

      if (( exc.hcht2 + exc.hcdt2 ) <= ( hctme-tt2 ))
      {
         tt2 += ( exc.hcht2 + exc.hcdt2 ); /* one more period */
         solver.hcswg *= 1.;
         solver.sws *= ONE;
         return ONE;
      }
      else if ( exc.hcht2 <= ( hctme-tt2 ))
      {
         solver.hcswg *= 0.;
         solver.sws *= null; /* sources OFF */
         return ONE;
      }
      else if ( hctme < exc.hcht2 )
      {
         tt2 = ZERO; /* initialize 2nd period */
         solver.hcswg *= 1.;
         solver.sws *= ONE;
         return ONE;
      }
      else
      {
         solver.hcswg *= 1.;
         solver.sws *= ONE;
         return ONE;
      };
      break;

     case 11: /* type: SWITCHED_SOURCES___ */

      if (( exc.hcht + exc.hcdt ) <= ( hctme-tt ))
      {
         tt += ( exc.hcht + exc.hcdt ); /* one more period */
         solver.sws = ONE; /* enable sources [switch ON] */
      }
      else if ( exc.hcht <= ( hctme-tt ))
      {
         solver.sws = null; /* sources OFF */
      }
      else if ( hctme < exc.hcht )
      {
         tt = ZERO; /* initialize period */
         solver.sws = ONE; /* sources ON */
      }
      else
      {
         solver.sws = ONE; /* sources ON */
      };

      if (( exc.hcht2 + exc.hcdt2 ) <= ( hctme-tt2 ))
      {
         tt2 += ( exc.hcht2 + exc.hcdt2 ); /* one more period */
         solver.sws *= ONE;
         return ONE;
      }
      else if ( exc.hcht2 <= ( hctme-tt2 ))
      {
         solver.sws *= null;
         return ONE;
      }
      else if ( hctme < exc.hcht2 )
      {
         tt2 = ZERO; /* initialize 2nd period */
         solver.sws *= ONE;
         return ONE;
      }
      else
      {
         solver.sws *= ONE;
         return ONE;
      };
      break;

     default: /* Unknown excitation type */

      fprintf( display, "\n Message from function %s : ", __func__ );
      fprintf( display, "\n Unknown excitation type !!!\n " );

      exit( EXIT_FAILURE );

      break;
   };
}
/*============================================================================*/
# endif /* DSC_HCRMDE != 0 */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/************************* end of function trnhcr(*) **************************/
