/* [ file: setmtp.h ] */
/* Update: October 14, 2007 */
/*----------------------------------------------------------------------------*/
typedef struct
{            
   signed char
      rtn, /* any returned character */
      opt; /* any returned option */

   short
      idx; /* media index */

   long
      ci, /* initial cell, to that media label <idx> shall be assigned */
      cf; /* final cell, to that ... */

/* E/H field relevant media parameters: */
   double
      eps, /* rel. permittivity tensor */
      myr, /* rel. permeability tensor */
      ke,  /* electric conductivity tns. [A/(V*m)] */
      km;  /* magnetic conductivity tns. [V/(A*m)] */
/*............................................................................*/
# if DSC_HCRMDE != 0
/* thermal [diffusion] media parameters: */

   double
      cv, /* specific heat [J/(K*m^3)] */
      kh; /* heat current conductivity [W/(K*m)] */
/*............................................................................*/
# if DSC_FLDMDE != 0
/* fluid media parameters: */

   short
      cnn; /* fluid connected component index; POSITIVE integer */

   double
      rm, /* mean mass density [Kg/m^3] */
      tm, /* mean temperature [C] */
      bm, /* mean thermal expansion coefficient [1/K] */
      cm, /* mean [ adiabatic ] compression coefficient [1/Pa] */
      ny, /* dynamic viscosity [Kg/(sec*m)] */
      q1, /* Cp/Cv - 1 [dimensionless] */
      td, /* dissipation time constant [sec] */
      LL; /* characteristic length [in Prandtl turbulence model, e.g.] */

   double
      gp[THREE], /* pressure gradient [Pa/m]*/
      gr[THREE]; /* gravitation acceleration [m/sec^2]*/

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* Gyrotropic media parameters: */

   double
      no,        /* plasma electron density  [1/m^3] */
      tr,        /* plasma curr.relax.time [seconds] */
      tg,        /* spin-spin  relax.time  [seconds] */
      ld,        /* LANDE factor */
      mi[DIMNS], /* int. magn. flux, plasma  [Tesla] */
      ms[DIMNS], /* saturation magnetization [Tesla] */
      hg[DIMNS]; /* internal magnetic field    [A/m] */

   char
      type[SHS_SIZE];
} MEDIUM;
/*************************** end of file setmtp.h *****************************/
