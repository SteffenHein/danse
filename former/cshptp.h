/* [ file: cshptp.h ] */
/* Update: 13 April 2007 */
/*-----------------------------------------------------------------------------
*//* geometrical cell parameter transfer structure
*//* All units are international units (mks), if not otherwise specified.
------------------------------------------------------------------------------*/
typedef struct
{
   signed char /* any return character [ may be used as an error flag, e.g. ] */
      rtn;

   signed char /* any option */
      opt;

   char
      skew; /* geometric skew indicator [ no skew:0 / skew:1] */

/*............................................................................*/
/* input parameters: */

   double
      c[CRNRS][DIMNS]; /* corner point coordinates */

   long
      cell;

   char
      face, port;

   TOPOLOGY *tpt;

   PARAMETERS *ppt;
/*............................................................................*/
/* computed geometric parameters: */

   double
      vol, /* cell volume [ m^3 ] */
      skv, /* cell skew volume [ m^3 ] */
      xn, yn, zn;

   double
      px[PORTS],
      py[PORTS],
      pz[PORTS],
      pm[PORTS],

      xp[PORTS],
      yp[PORTS],
      zp[PORTS],

      fx[FACES],
      fy[FACES],
      fz[FACES],
      fm[FACES],

      xf[FACES],
      yf[FACES],
      zf[FACES],
      vp[FACES]; /* face pyramide volume: vp[jj] = ( f[jj] | b[jj/2] )/6 */

   double
      e[PORTS][DIMNS], /* edge vectors */
     eu[PORTS][DIMNS], /* edge vectors relative to cell basis ub[] */

      p[PORTS][DIMNS], /* port vectors */
     pu[PORTS][DIMNS], /* port vectors relative to cell basis ub[] */

      f[FACES][DIMNS], /* face vectors; inner or outer - depending on */
                       /* face label !!! */
     nf[FACES][DIMNS], /* outer face normal vectors */
     fu[FACES][DIMNS], /* face vectors relative to cell basis ub[] */
     fa[DIMNS][DIMNS], /* opposite face vector arithmetic means */

      a[DIMNS][DIMNS], /* area vector matrix A = arv[i][j] */
     au[DIMNS][DIMNS], /* area vector matrix relative to cell basis ub[] */
     ai[DIMNS][DIMNS], /* Inverse are vector matrix A^-1 */

      b[DIMNS][DIMNS], /* node vector matrix B = ndv[i][j] */
     ub[DIMNS][DIMNS], /* Gram-Schmidt orthonormalzd node vectors [cell basis]*/
     bu[DIMNS][DIMNS], /* normalized Node vector matrix relative to ub[] */
     bi[DIMNS][DIMNS], /* inverse node vector basis B^-1 */

     cs[DIMNS][DIMNS], /* cell shape skew vectors */
     cu[DIMNS][DIMNS], /* cell shape skew vectors relative to cell basis */

     uv[DIMNS][DIMNS], /* any orthonormal basis */
     vu[DIMNS][DIMNS]; /* any coordinate vectors [ with respect to on ] */

} CSHAPE;
/*************************** end of file cshapt.h *****************************/
