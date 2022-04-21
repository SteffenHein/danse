/* [ file: ellwve.h ] */
# define DO_ELLWVE "ellwve(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function ellwve(*), DANSE release 1.0.                              *
*   [ in DANSE program package ]                                               *
*                                                                              *
*   Returns the cutoff wavelength of the eH11 principle mode in an elliptic    *
*   waveguide, given the principle axes aa, bb of the [ transversal ] cross    *
*   section ellipse.                                                           *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define USE_LC2S 0 /* 0: use function lc2a(*) [ W.Krank ] */
                    /* 1: use function lc2s(e) [ Markuwitz ] */
/*----------------------------------------------------------------------------*/
/*
struct wvegde  

   double 
      z0,     
      c0,
      w0,
      wc;
};
static struct wvegde wve = {ZERO};
*/
# if USE_LC2S == 1 /* use lc2s(*) [ Markuwitz ] */
/*----------------------------------------------------------------------------*/
# define DO_LC2S "lc2s(*)"
/*******************************************************************************
*                                                                              *
*  ANSI C function lc2s(*)                                                     *
*                                                                              *
*  Returns the ratio of the critical wavelength of the eH11 principle mode     *
*  in an elliptic waveguide to the circumference                               *
*                            s = 2q*E(e^2)/e                                   *
*  of the transverse cross section ellipse, as a function of the normalized    *
*  excentricity e. q denotes the focal distance and E(*) the complete ellip-   *
*  tic integral of the second kind.                                            *
*  [ cf. Markuwitz, Waveguide Handbook, p.83 ]                                 *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

double lc2s( double xx )
{
   static const double
      c0  =   5.4312997e-01,
      c1  =   3.8099888e-02,
      c2  = - 5.1181479e-01,
      c3  =   2.5083461e+00,
      c4  = - 3.9046635e+00,
      c5  =   2.1572924e+00;

   static double 
      yy  = ZERO; 

/*----------------------------------------------------------------------------*/

   if (( 0. <= xx )&&( xx < 1. ))
   {
      yy = ( c0 + xx*( c1 + xx*( c2 + xx*( c3 + xx*( c4 + xx*c5 )))));

      return yy;
   }
   else
   {
      printf( "\n\n Error message from function %s:", DO_LC2S );
      printf( "\n Transferred argument out of domain: x = %.16e !!!", xx );
      printf( "\n Legal domain:  0. <= x < 1"
         "\n [ check calling program ].\n " ); 
      exit( EXIT_FAILURE );
   };
}
/*************************** end of function lc2s(*) **************************/
# include "../math/ellipt.h"

# else /* if USE_LC2S != 1: use lc2a(*) [ W. Krank ] */

# define DO_LC2A "lc2a(*)"
/*******************************************************************************
*                                                                              *
*  ANSI C function lc2a(*)                                                     *
*                                                                              *
*  Returns the ratio of the critical wavelength of the principle ( eH11 )      *
*  mode of elliptic waveguide to the principle axis a as a function of b/a.    *
*  [ cf. Wolfgang Krank: Ueber die Theorie und Technik des elliptischen        *
*  Wellrohrleiters, Dissertation, TH Aachen/Telefunken 1982 ].                 *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

double lc2a( double xx ) /* xx = b/a */
{
   static const double
      c0  =   1.6610403e+00,
      c1  =   5.9173598e-02,
      c2  = - 2.8889208e-01,
      c3  =   7.4688422e-01,
      c4  = - 6.3208825e-01,
      c5  =   1.1496640e-01,
      c6  =   4.4892501e-02;

   static double 
      yy  = ZERO; 

/*----------------------------------------------------------------------------*/

   if (( 0. <= xx )&&( xx <= 1. ))
   {
      yy = ( c0 + xx*( c1 + xx*( c2 + xx*( c3 + xx*( c4 + \
                  xx*( c5 + xx*c6 ))))));
   }
   else
   {
      printf( "\n\n Error message from function %s:", DO_LC2A );
      printf( "\n Transferred argument out of domain: x = %.16e !!!", xx );
      printf( "\n Legal domain:  0. <= x <= 1"
         "\n [ check calling program ].\n " ); 
      exit( EXIT_FAILURE );
   };

   return yy;
}
/*************************** end of function lc2a(*) **************************/

# endif /* USE_LC2S != 1 */

/*============================================================================*/

double ellwve ( double aa, double bb )
{
/* allusions: */
/*
   extern struct wvegde wve;
*/
/* declarations: */

# if USE_LC2S == 1
   static double
      ee = ZERO,
      e2 = ZERO,
      qq = ZERO,
      rr = ZERO;
# endif

   static double
      ss = ZERO,
      wc = ZERO;

/* allusions: */

   double sqrt( double xx );
# if USE_LC2S == 1
   double lc2s( double xx );
   double ellipt( double xx );
# else
   double lc2a( double xx );
# endif
/*----------------------------------------------------------------------------*/

# if USE_LC2S == 1
   rr = aa*aa;
   ss = bb*bb;

   qq = sqrt( rr - ss ); /* the focal distance [ q/2 the semifocal distance ] */

   e2 = 1. - ss/rr; /* ee^2 */
   ee = sqrt( e2 ); /* the excentricity */
/*............................................................................*/
   ss = 2.*qq*ellipt( e2 )/ee;   /*                                           */
/*.............................*/
                       /* lc2s(*): the ratio of the [ prinicple eH11-mode ] */
   wc = ss*lc2s( ee ); /* critical wavelength to the elliptic circumference */
                       /* s = 2qE(e^2)/e, as function the normalized excent-*/
                       /* ricity e, q denoting the focal distance. */ 
                       /* [ cf. Markuwitz: Waveguide Handbook, p.83 ] */
# else
   ss = bb/aa;         /* lc2a(*): the ratio of the [ principle eH11 mode ] */
   wc = aa*lc2a( ss ); /* critical wavelength to the principle axis of the */
# endif                /* ellipse as a function of a/b. */
                       /* [ cf. Wolfgang Krank, Thesis, RWTH Aachen 1982 ] */
   return wc;
}
/*============================================================================*/
/************************ end of function ellwve(*) ***************************/
