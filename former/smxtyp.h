/* [ file: smxtyp.h ] */
/*-----------------------------------------------------------------------------
*//* s-matrix data transfer structure [ in header 'smxctp.h' ], DANSE 1.0
*//* All units are international units (mks), if not otherwise specified.
*//* [ Update: April 13, 2007 ]                                 Steffen Hein
------------------------------------------------------------------------------*/
typedef struct   
{
   signed char /* any return character [ may be used as an error flag, e.g. ] */
      rtn;

   char 
      opt; /* smatrx(*) options: 't' time_step, 'r': real s-parameters [ in */
           /* time-domain ], 'c': complex s-parameters [ in frequency-dmn ] */
   char
      loss, /* loss indicator [ lossless cell:0 / lossy cell:1 ] */
      skew, /* geometric skew indicator [ no skew:0 / skew:1] */
      isotrop; /* media isotropy indicator [ unisotropic: 0 / isotropic: 1 ] */ 

   short 
      med; /* media label */

   double
      adm, dt, dp, fr, omega, no, tr, tg, ld, vol, skv;

                       /* dt  = DSC time step                 [sec^-1       ] */
                       /* adm = y = field admittance          [Ohm^-1       ] */
                       /* no  = plasma electron density       [1/m^3        ] */
                       /* tr  = "      " relaxation time      [sec^-1       ] */
                       /* tg  = gyrmagn. relaxation time      [sec^-1       ] */
                       /* ld  = LANDE factor                  [dimensionless] */
   double
     mi[DIMNS],        /* mi = int.mag.fl.density,plasma[Tesla]               */
     ms[DIMNS],        /* ms = saturat. magnetization   [Tesla]               */
     hg[DIMNS],        /* hg = internal (static) magn.fld.[A/m]               */
     vp[FACES];        /* vp[jj] pyramide volume ( F[jj], b[jj/2] )/6         */

   double
    epr[DIMNS][DIMNS], /* ep = rel. permittivity tensor                       */
    epi[DIMNS][DIMNS],
    myr[DIMNS][DIMNS], /* my = rel. permeability "                            */
    myi[DIMNS][DIMNS],
     ke[DIMNS][DIMNS], /* ke = el.  conductivity "                            */
     km[DIMNS][DIMNS], /* km = magn conductivity "                            */

      c[CRNRS][DIMNS], /* corner point coordinates */

      e[PORTS][DIMNS], /* edge vectors */
     eu[PORTS][DIMNS], /* edge vectors relative to cell basis ub[]            */

      p[PORTS][DIMNS], /* port vectors */
     pu[PORTS][DIMNS], /* port vectors relative to cell basis ub[]            */

      f[FACES][DIMNS], /* face vectors */
     fu[FACES][DIMNS], /* face vectors relative to cell basis ub[]            */
     fa[DIMNS][DIMNS], /* opposite face vector means                          */

      a[DIMNS][DIMNS], /* area vector matrix A = arv[i][j]                    */
     au[DIMNS][DIMNS], /* area vector matrix relative to cell basis ub[]      */
     ai[DIMNS][DIMNS], /* Inverse are vector matrix A^-1                      */

      b[DIMNS][DIMNS], /* node vector matrix B = ndv[i][j]                    */
     bu[DIMNS][DIMNS], /* normalized Node vector matrix relative to ub[]      */
     ub[DIMNS][DIMNS], /* Gram-Schmidt orthonormalzd node vectors [cell basis]*/
     bi[DIMNS][DIMNS], /* inverse node vector basis B^-1                      */

     cs[DIMNS][DIMNS], /* cell geometric skew vectors                         */
     cu[DIMNS][DIMNS]; /* cell geometric skew vectors relative to cell basis  */

   double
     se[SORDR][SORDR], /* S-parameters [ electric type, time domain ]         */
    ser[DIMNS][DIMNS], /* S-parameters [ electric type, frq. dmn, real part ] */
    sei[DIMNS][DIMNS], /* S-parameters [ electric type, frq. dmn, imag.part ] */

     sm[SORDR][SORDR], /* S-parameters [ magnetic type ]                      */
    smr[DIMNS][DIMNS], /* S-parameters [ magnetic type, frq. dmn, real part ] */
    smi[DIMNS][DIMNS], /* S-parameters [ magnetic type, frq. dmn, imag.part ] */

     ge[DIMNS][DIMNS], /* gyroelectric bias                                   */
     gm[DIMNS][DIMNS], /* gyromagnetic bias                                   */

    te0[DIMNS][DIMNS], /* A*(ke/2+(EP+GE)/dt)*(B^-1)/4y                       */
    te1[DIMNS][DIMNS], /* A*(ke/2-(EP+GE)/dt)*(B^-1)/4y                       */
    tei[DIMNS][DIMNS], /* TE0^-1                                              */
    tet[DIMNS][DIMNS], /* (TE0^-1)*TE1/y                                      */
    tm0[DIMNS][DIMNS], /* y*A*(km/2+(MY+GM)/dt)*(B^-1)/4                      */
    tm1[DIMNS][DIMNS], /* y*A*(km/2-(MY+GM)/dt)*(B^-1)/4                      */
    tmi[DIMNS][DIMNS], /* TM0^-1                                              */
    tmt[DIMNS][DIMNS]; /* y*(TM0^-1)*TM1                                      */

   char 
     name[STS_SIZE],
     text[STS_SIZE], 
     etyp[SHS_SIZE], 
     mtyp[SHS_SIZE],
     getp[SHS_SIZE],
     gmtp[SHS_SIZE];

}  S_MATRIX;
/*************************** end of file smxctp.h *****************************/
