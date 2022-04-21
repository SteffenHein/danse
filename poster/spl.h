/* [ file: spl.h ] */
/* Update: July 12, 2001 */
/*----------------------------------------------------------------------------*/
/* structure type definition of function spline(*) */

typedef struct
{
   signed char
      rtn;

   long
      mm, /* number of support points */
      nn; /* number of interpolated points */

   double
      intgr, fmin, fmax,
      vct[SPL_SPNTS][TWO], /* given points */
      dmn[SPL_INTPL],      /* interpolated points */
      fct[SPL_INTPL],
      drv[SPL_INTPL];

} SPLINES;
/*************************** end of file spl.h *******************************/
