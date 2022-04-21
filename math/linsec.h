/* [ file: linsec.h ] */
# define DO_LINSEC "linsec(*)"
/*******************************************************************************
*                                                                              *
*   ANSI-C function linsec(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Computes intersection point  ( px , py ) of straight lines characterized   *
*   by equations [ Hesse's normal form for lines ]:                            *
*                                                                              *
*                              < p , u^> = r                                   *
*                              < p , v^> = s     ,                             *
*                                                                              *
*   where  <.,.> denotes the scalar product and  u^, v^ are the unit vectors   *
*   pointing from the origin to the respective line  and perpendicular to it.  *
*   alfa  and  beta  are the angles of  u^ and  v^, respectively , with  the   *
*   positive x-axis.                                                           *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*---------------------------------------------------------------------------*/
struct linesect  
{
   double px,
	  py;
};
static struct linesect lin = {ZERO};
/*============================================================================*/

int linsec ( double alfa, double beta, double r, double s )
{
   static double det = ZERO,
		 dtx = ZERO,
		 dty = ZERO,
	     epsilon = ZERO;

   double sin ( double x );
   double cos ( double x );
   double fabs( double x );
/*---------------------------------------------------------------------------*/
   epsilon = 1.e-277;

   det = cos(alfa)*sin(beta)-sin(alfa)*cos(beta);

   if ( fabs( det ) < epsilon )
   {
      printf("\n\n Warning from function '%s' : ",DO_LINSEC);
      printf("\n\n sin(alfa-beta) = ZERO; parallel lines ");
      printf("\n - no point of intersection !!! ");
      return null;
   }
   else
   {
      dtx = r*sin(beta) - s*sin(alfa);
      dty = s*cos(alfa) - r*cos(beta);
      lin.px = dtx / det;
      lin.py = dty / det;
      return ONE;
   };
}
/*============================================================================*/
/************************* end of function linsec(*) **************************/
