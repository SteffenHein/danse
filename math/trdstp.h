/* [ file: trdstp.h ] */
/* Update: May 01, 2007 */
/*----------------------------------------------------------------------------*/
/* structure type definition header of function triads(*) */

typedef struct
{
   short
      rtn;

   char
      opt;

   double
      det,
      n[THREE],
      v[THREE][THREE],
      uv[THREE][THREE],
      vu[THREE][THREE];

} TRIADS;
/**************************** end of file trdstp.h ****************************/
