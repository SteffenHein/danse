/* [ file: wvepar.h ] */
# define DO_WVEPAR "wvepar(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function wvepar(*); DANSE release 1.0.                              *
*   This function returns wave propagation parameters of fundamental mode      *
*   for waveguide types *type = 'r'ectangular, 'c'ircular, or, 'e'lliptic.     *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: April 14, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# ifndef MY_VAC_
   # define MY_VAC_ ( 1.256637061e-06 )
# endif
# ifndef EPS_VAC
   # define EPS_VAC ( 8.854190000e-12 )
# endif
/*----------------------------------------------------------------------------*/
# include "ellwve.h"
/*----------------------------------------------------------------------------*/
typedef struct
{
   double 
      z0, z1, c0, c1, vp, vg, w0, w1, fc, wc, wg, zh, ze, q;
}  WVGDPAR;
/*============================================================================*/

WVGDPAR *\
wvepar( char *type, double aa, double bb,
        double eps, double my, double ff )
{
/* declarations: */

   static WVGDPAR wve = { null },
                 *str = &wve;

   static double ss = ZERO;

   static const double bound = 1.e-15;

/* allusions: */

   double fabs ( double x );
/*----------------------------------------------------------------------------*/
   wve.z0 = sqrt( MY_VAC_/EPS_VAC );
   wve.c0 = 1./sqrt( MY_VAC_*EPS_VAC );
   wve.w0 = wve.c0/ff;

   if ( bound < fabs( eps ) )
   {
      ss = sqrt( eps );

      wve.z1 = wve.z0/ss;
      wve.c1 = wve.c0/ss;
   }
   else
   {
      wve.z1 = wve.z0;
      wve.c1 = wve.c0;
   };

   if ( bound < fabs( my ) )
   {
      ss = sqrt( my );

      wve.z1 *= ss;
      wve.c1 /= ss;
   };
   wve.w1 = wve.c1/ff;

/* critical wavelength/frequency */

   if ( null == strncmp( type, "rectangular", THREE ))
   { /* rectangular waveguide: */
      wve.wc = 2.*aa;             /* TE10 mode cutoff [free] wave length */
   }
   else if ( null == strncmp( type, "circular", THREE ))
   { /* circular waveguide: */ 
      wve.wc = PI*aa/1.84118;     /* TE11 mode cutoff [free] wave length */
   }
   else if ( null == strncmp( type, "elliptical", THREE ))
   { /* elliptical waveguide: */
      wve.wc = ellwve( aa, bb );  /* TE11c mode cutoff [free] wave length */
   }
   else if ( null == strncmp( type, "tem", THREE ))
   { /* TEM [Transversal ElectroMagnetic] */
      wve.wc = HUGE_VALF;         /* critical wavelength */
      wve.fc = ZERO;              /* critical frequency */
      wve.wg = wve.w1;            /* waveguide wavelength */
      wve.ze = wve.z1;
      wve.zh = wve.z1;
      wve.vp = wve.c1;            /* phase velocity */
      wve.vg = wve.c1;            /* group velocity */
      wve.q  = ZERO;

      return str;
   }
   else if (( null == strncmp( type, "micro_strip" , THREE ))||
            ( null == strncmp( type, "micro_strip" , THREE )))
   { /* micro_strip [TEM]: */
      wve.wc = HUGE_VALF; /* critical wavelength */
      wve.fc = ZERO;              /* critical frequency */
      wve.wg = wve.w1;            /* waveguide wavelength */
      wve.ze = wve.z1;
      wve.zh = wve.z1;
      wve.vp = wve.c1;            /* phase velocity */
      wve.vg = wve.c1;            /* group velocity */
      wve.q  = ZERO;

      return str;
   }
   else if ( null == strncmp( type, "coax", THREE ))
   { /* coaxial line: */
      wve.wc = PI*( aa + bb )/2.; /* critical wavelength  */
      wve.fc = wve.c1/wve.wc;     /* critical frequency   */
      wve.wg = wve.w1;            /* waveguide wavelength */
      wve.ze = wve.z1;
      wve.zh = wve.z1;
      wve.vp = wve.c1;            /* phase velocity       */
      wve.vg = wve.c1;            /* group velocity       */
      wve.q  = ZERO;

      return str;
   };
      
   wve.fc = wve.c1/wve.wc;        /* critical frequency     */
   wve.q = wve.fc/ff; 

   if ( wve.q < 1. ) 
   {
      wve.q = sqrt( 1. - wve.q*wve.q );      

      if ( ZERO < fabs( wve.q ) )
      {
         wve.wg = wve.w1/wve.q;   /* waveguide wavelength    */
         wve.zh = wve.z1/wve.q;   /* TE mode field impedance */
      };
      wve.ze = wve.z1*wve.q;      /* TH mode field impedance */
      wve.vp = wve.c1/wve.q;      /* phase velocity          */
      wve.vg = wve.c1*wve.q;      /* group velocity          */
   }
   else
   {
      wve.ze = ZERO;
      wve.zh = ZERO;
      wve.vp = HUGE_VALF;
      wve.vg = ZERO;
   };

   return str;
}
/*============================================================================*/
/************************ End of function wvepar(*) ***************************/
