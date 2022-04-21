/* [ file: hcrstp.h ] */
/*-----------------------------------------------------------------------------
*//* heat current s-matrix data transfer structure
*//* DANSE relaease 1.0.
*//* All units are international units (mks), if not otherwise specified.
*//* [ Update: July 22, 2007 ]                           Steffen Hein
*//*---------------------------------------------------------------------------
*//* heat current s-matrix data transfer structure
------------------------------------------------------------------------------*/
typedef struct
{
   signed char /* any return character [ may be used as an error flag, e.g. ] */
      rtn;

   char 
      opt; /* hcrsmx(*) options: 't' time_step, 'r': real s-parameters [ in */
           /* time-domain ], 'c': complex s-parameters [ in frequency-dmn ] */
   char
      orth, /* orthogonal cell indicator [ non-orth:0 / orth:1 ] */
      skew, /* geometric skew indicator [ no skew:0 / skew:1 ] */
      loss, /* loss indicator [ lossless cell:0 / lossy cell:1 ] */
      isotrop; /* media isotropy indicator [ unisotropic: 0 / isotropic: 1 ] */ 

   short 
      med; /* media label */
/*............................................................................*/
/* input parameters: */

   double
      dt,  /* time step [sec] */
      dp,  /* phase shift [rad] */
      ke,  /* electric conductivity [A/(V*m) = S/m] */
      km,  /* magnetic conductivity [V/(A*m) = Ohm/m] */
      kh,  /* heat current conductivity [W/(K*m)] */
      cv;  /* heat capacity [J/(K*m^3)] */
/*............................................................................*/
# if DSC_FLDMDE != 0
   double
      rm, /* mean mass density [Kg/(m^3)] */
      tm, /* mean temperature [C] */
      bm, /* mean expansion coefficient [1/K] */
      cm, /* mean [ adiabatic ] compression coefficient [1/Pa] */
      ny, /* dynamic viscosity [Kg/(sec*m)] */
      q1, /* Cp/Cv - 1 [dimensionless] */
      td, /* dissipation time constant [sec] */
      LL; /* characteristic length [in Prandtl turbulence model, e.g.] */

   double
      gr[THREE], /* gravitational acceleration [m/(sec^2)]*/ 
      gp[THREE]; /* pressure gradient [Kg/((sec*m)^2)] */ 

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
 double
      c[CRNRS][THREE]; /* corner point coordinates */
/*............................................................................*/
/* computed geometric parameters: */

   double
      vol, /* cell volume [m^3] */
      skv; /* cell skew volume [m^3] */

   double
     fm[FACES],        /* size of face */
     vp[FACES];        /* vp[jj] pyramide volume ( F[jj], b[jj/2] )/6 */

   double
      e[PORTS][THREE], /* edge vectors */
     eu[PORTS][THREE], /* edge vectors relative to cell basis ub[] */

      p[PORTS][THREE], /* port vectors */
     pu[PORTS][THREE], /* port vectors relative to cell basis ub[] */

      f[FACES][THREE], /* face vectors */
     fu[FACES][THREE], /* face vectors relative to cell basis ub[] */
     fa[THREE][THREE], /* opposite face vector arithmetic means */

      a[THREE][THREE], /* area vector matrix A = arv[i][j] */
     au[THREE][THREE], /* area vector matrix relative to cell basis ub[] */
     ai[THREE][THREE], /* Inverse are vector matrix A^-1 */

      b[THREE][THREE], /* node vector matrix B = ndv[i][j] */
     bu[THREE][THREE], /* normalized Node vector matrix relative to ub[] */
     ub[THREE][THREE], /* Gram-Schmidt orthonormalzd node vectors [cell basis]*/
     bi[THREE][THREE], /* inverse node vector basis B^-1 */

     cs[THREE][THREE], /* cell shape skew vectors */
     cu[THREE][THREE]; /* cell shape skew vectors relative to cell basis */
/*............................................................................*/
/* computed [ temperature updating ] s-parameters [-> elf.smx<N>] */

   double
     hdt,  /* heat diffusion time step */
      ct;  /* heating coefficient [ dt/cv ] */

# if DSC_FLDMDE != 0
   double
     fdt,  /* fluid [ viscous ] diffusion time step */
      ft; /* fluid updating coefficient [ dt/rm ] */
# endif

/* face form vectors [ (F)*(B^-1) ] */
   double
     s[FACES][THREE];  /* face form vectors [ (F)*(B^-1) ] */

   char
     name[STS_SIZE],
     text[STS_SIZE], 
     ttyp[STS_SIZE]; 

} HCRSMX;
/*************************** end of file hcrstp.h *****************************/
