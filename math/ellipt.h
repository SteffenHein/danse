/* [ file: ellipt.h ] */
# define DO_ELLIPT "ellipt(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function ellip(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Returns the elliptic normal integral E( mm ) for any real mm in the        *
*   interval [ 0., 1. ] using the polynomial approximation (17.3.36) of        *
*   Abramowitz / Stegun:                                                       *
*   Handbook of Mathematical Functions, nineth printing, p. 592 .              *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*============================================================================*/

double ellipt( double mm )
{
   static const double
      a1 = .44325141463,
      a2 = .06260601220,
      a3 = .04757383546,
      a4 = .01736506451,
      
      b1 = .24998368310,
      b2 = .09200180037,
      b3 = .04069697526,
      b4 = .00526449639,

      lbound = .00000000000001,
      ubound = .99999999999999; 

   static double 
      mm1 = ZERO,
      yy = ZERO;
  
   double log ( double xx );

/*---------------------------------------------------------------------------*/

   if (( lbound <= mm )&&( mm <= ubound ))
   {
      mm1 = 1. - mm;

      yy = ( mm1*( b1 + mm1*( b2 + mm1*( b3 + mm1*b4 )))) * log( 1./ mm1 );
      yy += ( 1. + mm1*( a1 + mm1*( a2 + mm1*( a3 + mm1*a4 ))));
   }
   else if (( 0. <= mm )&&( mm < lbound ))
   {
      yy = PI/2.;
   }
   else if (( ubound < mm )&&( mm <= 1. ))
   {
      yy = 1.;
   }
   else
   {
      fprintf( stderr,
         "\n Error message from function %s :", DO_ELLIPT );
      fprintf( stderr,
         "\n\n Transferred argument m out of domain: m = %.16e !!!", mm );
      fprintf( stderr,
         "\n Observe domain of elliptic function E(m): 0. <= m <= 1\n "
         "[ check calling process.]\n" );
      exit( EXIT_FAILURE );
   };

   return yy;
}
/*============================================================================*/
/************************** end of function ellipt(*) *************************/
