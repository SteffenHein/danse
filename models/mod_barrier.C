/* [ file: model.c ] */
# define DSC_MODEL "mod_barrier.C"
/*******************************************************************************
*                                                                              *
*   DANSE - Discrete Approximation of the Navier Stokes Equations              *
*                                                                              *
*   DSC model generation and evaluation functions model.c and modval.c         *
*   Prototype: 'flow against inclined board'                                   *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 15, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define MOD_DEFLT 3 /* default models */
/*----------------------------------------------------------------------------*/
/*
*//* THIS MODEL MAKES ONLY SENSE WITH CFD OPTIONS COMPILED \
*//* compile FORMER.C with DSC_HCRMDE 1, DSC_FLDMDE 1 in CONFIG.H
*/
/*____________________________________________________________________________*/


/*____________________________________________________________________________*/
# ifdef _PRE_MODEL
/* THIS SECTION IS COMPILED ONLY WITH OPTION -D_PRE_MODEL */
/*******************************************************************************
*                                                                              *
*   DANSE - Discrete Approximation of the Navier Stokes Equations              *
*                                                                              *
*   DSC model generation function model.c                                      *
*   [ Prototype: 'flow against inclined board' ]                               *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define _POSIX_SOURCE 1 /* some headers of the POSIX.1 standard will be used */
/*----------------------------------------------------------------------------*/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
# include <unistd.h> /* system specification header, cf. sysconf( ), etc. */
# include <time.h>   /* cf. time( ), ctime( ), asctime( ), localtime( ), etc. */
/*----------------------------------------------------------------------------*/
# include "./math/maths.h"  /* 'my' computation environment headers */
# include "./math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "./CONFIG.H" 
/*----------------------------------------------------------------------------*/
# if DSC_FLDMDE == 0
# error !!!! THIS MODEL MAKES ONLY SENSE WITH CFD OPTIONS COMPILED \
- compile FORMER.C with DSC_HCRMDE 1, DSC_FLDMDE 1 in CONFIG.H
# endif
/*____________________________________________________________________________*/
/* Edit and customize this header for FORMER.C configuration: */
# include "./former/FORMER.CONF" 
/*----------------------------------------------------------------------------*/
/* A macro that gives a conveniently readable form to the string returned */
/* with [ the pointer of ] function ctime(&nseconds) */
# include "./tools/TIMEFORM.M"
/*____________________________________________________________________________*/
# if USE_NCURSES == 1
   # include <ncurses.h>
# endif 
/*____________________________________________________________________________*/
# ifndef DSC_ADJCLL
   # define DSC_ADJCLL 1 
# endif
/*----------------------------------------------------------------------------*/
# define MOD_ENBLSKN 1 /* 1: enable [ skin effect ] heat source bnd faces */
# define MOD_ENBLESC 1 /* 1: enable envrmnt-surface heat conductivity */
/*----------------------------------------------------------------------------*/
/* excitation type: EXC_MAXW = "DIRAC_PULSE", "HARMONIC", "MULTIPLE_HARMONIC" */
/*
# define EXC_MAXW  "DIRAC_PULSE"
# define EXC_MAXW  "GAUSS_PULSE"
# define EXC_MAXW  "HARMONIC"
# define EXC_MAXW  "SMOOTH_HARMONIC"
# define EXC_MAXW  "OSCILLATING_GAUSS__" [e.g.]

# define EXC_HCRR  "PASSIVE____________"
# define EXC_HCRR  "PERIODIC_PULSE_____"
# define EXC_HCRR  "DOUBLE_PERIODIC____" [e.g.]
*/
# define EXC_MAXW  "OSCILLATING_GAUSS__"
# define EXC_HCRR  "PASSIVE____________"
/*----------------------------------------------------------------------------*/
# include "./former/formtp.h"
# include "./tools/txctyp.h"
/*----------------------------------------------------------------------------*/
# define BLC_P_BLIMTS 2 /* 1: display block limits , first period [k=0]       */
                        /* 2: "       "     "      , all periods              */
/*----------------------------------------------------------------------------*/
# define ELECTRIC_WALL -1
# define MAGNETIC_WALL -2

# ifndef EPS_VAC
   # define EPS_VAC ( 8.8541878170e-12 ) /* vacuum permittivity [A*sec/(V*m)] */
# endif
# ifndef MY_VAC_
   # define MY_VAC_ ( 1.2566370614e-06 ) /* vacuum permeability [V*sec/(A*m)] */
# endif
/*----------------------------------------------------------------------------*/
/*
Transfer structure for DSC model dependent parameters

The following structure  'trf' of type 'transfer' may be arbitrarily modified,
as wanted, by  introducing model dependent variables ( e.g. grid character -
istics, lengths, material constants etc. ). This  structure  has been created
for parameter transfer between the model defining functions, viz. systop(*),
syssmx(*), sysbnd(*), sysexc(*), and sysval(*).
*/
/*                      | maximum number M,                                   */
/*                      V                                                     */
# define MOD_BLOCKS    10 /* maximum M, BLOCK(0),...,BLOCK(M)                 */
# define MOD_DOMAINS    3 /* maximum M, blc[1],...,blc[M]                     */
# define MOD_LAYERS     5 /* layers [in z-direction]                          */
# define MOD_OPERATS   10 /* operations etc., trf.c[0],...,trf.c[M]           */
# define MOD_DIVISNS   20 /* divisions etc., trf.n[0],...,trf.n[M]            */
# define MOD_PARMTRS   40 /* parameters trf.s[0],...,trf.s[M]                 */
# define MOD_BOUNDRS  500 /* maximum M, trf.bdn[1],...,trf.bdn[M-1]           */
# define MOD_EXCITES  500 /* maximum M, trf.exn[1],...,trf.exn[M-1]           */
# define MOD_EVLUATE   50 /* maximum M, trf.val[0],...,trf.val[M-1]           */
# define MOD_ARCS       3 /* maximum M, trf.alfa[0],...,trf.alfa[M]           */
# define MOD_POINTS    20 /* maximum M, pnt[0],...,pnt[M-1]                   */
/*............................................................................*/
# if DSC_HCRMDE != 0
   # define MOD_HCBNDRS 4000
   # define MOD_HCEXCTS  200
   # define MOD_HCEVALS   50
# endif
/*............................................................................*/
struct transfer 
{
   double 
      r1, r2, dz,      /* special geometric and electric parameters */
      lgt, trv,      
      epsr, kppa, tgd,
      rsqr1, rsqr2,    /* skin effect surface resistance */
      uu, nrm;         /* special electric parameters */
/*............................................................................*/
/* this part is canonical: */

   double 
      fr, cwp, /* frequency, cw power */
      a, b, ra, ri; /* canonical waveguide parameters etc.*/

   double 
      s[MOD_PARMTRS+ONE], /* real model parameters [ lengths, arcs, e.g.] */
      bdd[MOD_BOUNDRS], /* boundary permittivity */
      euu[MOD_EXCITES], /* excitation voltage */
      alfa[MOD_ARCS+ONE];

   long
      pp,
      pm,
      bii[MOD_BOUNDRS], /* initial boundary cell label */
      eii[MOD_EXCITES], /* initial excited cell label */
      vii[MOD_EVLUATE]; /* initial evaluated cell label */

   short
      bnd, exc, val, /* boundary, excited and evlt'd port section counters */
      c[MOD_OPERATS+ONE], /* operation [='c'omputation ] modes */
      n[MOD_DIVISNS+ONE], /* integer model parameters [ divisions, e.g.] */
      bnn[MOD_BOUNDRS], /* number of periods in boundary port ... */
      bpp[MOD_BOUNDRS], /* boundary port period */
      enn[MOD_EXCITES], /* number of periods in excited port */
      epp[MOD_EXCITES], /* excited port period */
      vnn[MOD_EVLUATE], /* number of periods in evaluated port */
      vpp[MOD_EVLUATE]; /* excited port period */

   char 
      ref, /* reference line indicator [ 0 or 1 ] */
      bff[MOD_BOUNDRS], /* boundary face label [ 0,...,5 ] */
      eff[MOD_EXCITES], /* electrically excited face label [ 0,...,5 ] */
      vpt[MOD_EVLUATE], /* evaluated electric port index [ 1,...,12 ] */
      btp[MOD_BOUNDRS][STS_SIZE], /* waveguide boundary type */
      etp[MOD_EXCITES][STS_SIZE], /* excitation type */
      vtp[MOD_EVLUATE][STS_SIZE], /* type of evlt'd port ["e","h",...]*/
      vtx[MOD_EVLUATE][STS_SIZE], /* text string characterizing the port */

      ctx[MOD_OPERATS+ONE][STS_SIZE], /* computation mode[text strings */
      ntx[MOD_DIVISNS+ONE][STS_SIZE], /* xy-divisions text strings */
      mtx[MOD_DOMAINS+ONE][STS_SIZE], /*  z-divisions text strings */
      stx[MOD_PARMTRS+ONE][STS_SIZE], /* parameter text strings */

     *wgtype; /* waveguide type [ char string ] */
/*............................................................................*/
# if DSC_HCRMDE != 0  
/* model parameters for temperature propagation [ heat current mode ]: */  

   double
      temp, tempi,         /* thermal operation parameters */
      kh0, kh1, kh2, kh3,  /* thermal material parameters [ heat conduct.] */
      cv0, cv1, cv2, cv3,  /* thermal material parameters [ heat capacit.] */
      tm0, tev,            /* mean air and environment temperatures [Celsius] */
      tcool,               /* coolant temperature [Celsius] */
      shc,                 /* surface heat conductivity [W/(m^2*K)] */
      rm0,                 /* rm0: mean density [Kg/(m^3)] */
      bm0,                 /* bm0: thermal expans. coeff [1/K] */
      cm0,                 /* cm0: adiabatic compression coeff [1/Pa] */
      ny0,                 /* ny0: dynamic  viscosity [Kg/(m*sec)] */
      q1,                  /* q1:  Cp/Cv - 1 [dimensionless] */
      gr,                  /* gravitational acceleration [m/sec^2] */
      td, LL,              /* dissipation time [sec], charact length[m] */
      abs0, rad0, rad1,    /* radiation paramaters */
      hcrr,                /* imposed heat current density */
      hcht, hcdt,          /* heat current excitation hold & delay time */
      hcht2, hcdt2,        /* dito, for 2nd pulse [ if need is ] */
      nus,                 /* Nusselt number [ at no-slip boundary face ] */
      ti,                  /* inflow fluid temperature */
      uf,                  /* inflow fluid velocity */
      vf,                  /* outflow fluid velocity */

      crsdt;               /* coarsening time period [sec] */

   short
      bhc, ehc, vhc,    /* thermal boundary, excitation, and eval'n counter */
      bcn[MOD_HCBNDRS], /* number of periods */
      bcp[MOD_HCBNDRS], /* period */
      ecn[MOD_HCEXCTS], /* number of excited temperature nodes */
      ecp[MOD_HCEXCTS], /* period */
      vcn[MOD_HCEVALS], /* number of excited temperature nodes */
      vcp[MOD_HCEVALS]; /* period */

   char
      bcf[MOD_HCBNDRS], /* hcurr or temp boundary cell face index */
      ecf[MOD_HCEXCTS], /* excited hcurr or temp cell face index */
      vcf[MOD_HCEVALS], /* evaluated hcurr or temp cell face index */
      bct[MOD_HCBNDRS][STS_SIZE], /* thermal boundary type */
      ect[MOD_HCEXCTS][STS_SIZE], /* thermal excitation type */
      vct[MOD_HCEVALS][STS_SIZE], /* thermal and/or fluid evaluation type */
      vcx[MOD_HCEVALS][STS_SIZE]; /* any text characterizing evlt'd port */

   long
      bci[MOD_HCBNDRS], /* hcurr or temperature boundary; initial cell index */
      eci[MOD_HCEXCTS], /* hcurr or temperature excit'n; initial cell index */
      vci[MOD_HCEVALS]; /* hcurr or temperature evalt'n; initial cell index */

   double
      bcr[MOD_HCBNDRS], /* bcr, bcs: any real parameters, such as bndry heat */
      bcs[MOD_HCBNDRS], /* current density, temperature, or anything else */
      ecr[MOD_HCEXCTS]; /* excitet heat current density, temperature, or else */

# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
};
/*----------------------------------------------------------------------------*/
static struct transfer 
   trf = {ZERO};

static double 
   pnt[MOD_POINTS][THREE] = {{ZERO}};
/*----------------------------------------------------------------------------*/
static FORMSTATE *spt = NULL;
/*******************************************************************************
*                                                                              *
*   This function needs the quadrangular and triangular block definition       *
*   functions, on header file "puzzab.h", which are included here.             *
*                                                                              *
*   Munich, September 2000                                  Steffen Hein       *
*                                                                              *
*******************************************************************************/
# include "./tools/puzzab.h" /* block generation functions */
# include "./tools/wvepar.h"
# include "./tools/panels.h"
/*----------------------------------------------------------------------------*/
# include "./tools/inputs.h"
/*============================================================================*/

void deflt_operts( ) /* default operations */
{ 
/* allusions: - */
/*
   extern FORMSTATE *spt;

   extern struct transfer trf;
   extern struct blcstruc blc;
*/
   static short
      ll = LINLEN;

   void \
      cpylne( char lne[], const char *txt, const char *brc, short linlen );
/*----------------------------------------------------------------------------*/
/* write only CONNECTED (!) STRINGS into the string copy function strcpy(*)   */
/* [ 2nd. argument       |  ]                                                 */
/*                       V                                                    */
   strcpy( trf.ctx[0], "operations_[computation_modes]" );

   cpylne( trf.ctx[1], "reference_line/_structure", "void", ll );
   cpylne( trf.ctx[2], "in_channel:_no-slip/_free_slip_boundaries", "0/1", ll );

   trf.c[0] =  2; /* number of operation modes */
   trf.c[1] =  1; /* structure [ fixed ] */

/*.............................. defaults ....................................*/
# if MOD_DEFLT == 1  
   trf.c[2] = 1;
# elif MOD_DEFLT == 2
   trf.c[2] = 1;
# elif MOD_DEFLT == 3
   trf.c[2] = 1;
# elif MOD_DEFLT == 4
   trf.c[2] = 1;
# else /* [ cf. MOD_DEFLT == 1 ] */
   trf.c[2] = 1;
# endif
/*............................................................................*/

   return;
}
/*============================================================================*/

short rvise_operts( void )
{ 
/*
   extern FORMSTATE *spt;

   extern struct transfer trf;
   extern struct blcstruc blc;
*/
/* operation marks: */

   trf.c[ONE] = ONE;
/* hence, the following check is not needed */
/* [ it is left anyway as a canonical request that doesn't harm  ] */  

   if( null == trf.c[ONE] )
   {
      trf.ref = ONE;
      strcpy(( spt->tpt->text ), "reference_line" );
   }
   else
   {
      trf.c[ONE] = ONE;
      trf.ref = null;
      strcpy(( spt->tpt->text ), "flow_against_inclined_board" );
   };

   if ( trf.c[TWO] != null )
      trf.c[TWO] = ONE;
/*
   if ( trf.c[THREE] != null )
      trf.c[THREE] = ONE;

   if ( trf.c[FOUR] != null )
      trf.c[FOUR] = ONE;

   if (( trf.c[FIVE] != null )
     &&( trf.c[FIVE] != ONE ))
      trf.c[FIVE] = TWO;

   if ( trf.c[SIX] != null )
      trf.c[SIX] = ONE;
*/
   return null;
}
/*============================================================================*/

void deflt_divsns( void ) /* default configuration and divisions              */
{ 
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      ll = LINLEN;

/* prototypes: */

   void \
      cpylne( char lne[], const char *txt, const char *brc, short linlen );
/*----------------------------------------------------------------------------*/
/* _z_divisions:                                                              */
                                                                          
/* write only CONNECTED (!) STRINGS into the string copy function strcpy(*)   */
/* [ 2nd. argument       |  ]                                                 */
/*                       V                                                    */
   strcpy( trf.mtx[0], "vertical_[z]domains" );

   cpylne( trf.mtx[1], "vertical_divisions,_domain__null", "", ll );
/*
   cpylne( trf.mtx[2], "vertical_divisions,_domain__1", "", ll );
   cpylne( trf.mtx[3], "vertical_divisions,_domain__2", "", ll );
...
   cpylne( trf.mtx[8], "vertical_divisions,_domain__7", "=m6;_dependent", ll );
   cpylne( trf.mtx[9], "vertical_divisions,_domain__8", "=m5;_dependent", ll );
   cpylne( trf.mtx[10], "vertical_divisions,_domain__9", "=m4;_dependent", ll );
   cpylne( trf.mtx[11], "vertical_divisions,_domain_10", "=m3;_dependent", ll );
   cpylne( trf.mtx[12], "vertical_divisions,_domain_11", "=m2;_dependent", ll );
*/
/* z-divisions. */
   blc.m[null] = 1;

/*.............................. defaults ....................................*/
# if MOD_DEFLT == 1
   blc.m[1] = 1; /* number of divisions, domain 0 [ source line ] */
# elif MOD_DEFLT == 2
   blc.m[1] = 1;
# elif MOD_DEFLT == 3
   blc.m[1] = 1;
# elif MOD_DEFLT == 4
   blc.m[1] = 1;
# else /* [ cf. MOD_DEFLT == 1 ] */
   blc.m[1] = 1;
# endif
/*............................................................................*/
   
/* xy_divisions, identifiers etc.:                                            */

/* write only CONNECTED (!) STRINGS in the string copy function strcpy(*) !!! */
/* [ 2nd. argument       |  ]                                                 */
/*                       V                                                    */
   strcpy( trf.ntx[0] , "divisions,_identifiers,_operation_marks,_etc." );

   cpylne( trf.ntx[1], "bc_divisions,_block____1", "", ll );
   cpylne( trf.ntx[2], "bc_divisions,_block____2", "", ll );
   cpylne( trf.ntx[3], "bc_divisions,_blocks_3,5", "", ll );
   cpylne( trf.ntx[4], "bc_divisions,_block____4", "", ll );
   cpylne( trf.ntx[5], "ab_divisions,_blocks_3,4", "", ll );
   cpylne( trf.ntx[6], "bc_divisions,_blocks_1,2", "dependent", ll );
   cpylne( trf.ntx[7], "bc_divisions,_blocks_6,7", "dependent", ll );

   trf.n[0] = 7; /* number M of parameters trf.n[1],...,trf.n[M] */

/* xy_divisions: */
/*.............................. defaults ....................................*/
# if MOD_DEFLT == 1
   trf.n[1]  =  20;
   trf.n[2]  =  50;
   trf.n[3]  =   8;
   trf.n[4]  =   5;
   trf.n[5]  =   1;
# elif MOD_DEFLT == 2
   trf.n[1]  =  20;
   trf.n[2]  =  50;
   trf.n[3]  =   8;
   trf.n[4]  =   5;
   trf.n[5]  =   1;
# elif MOD_DEFLT == 3
   trf.n[1]  =  20;
   trf.n[2]  =  50;
   trf.n[3]  =   8;
   trf.n[4]  =   5;
   trf.n[5]  =   1;
# elif MOD_DEFLT == 4
   trf.n[1]  =  25;
   trf.n[2]  = 130;
   trf.n[3]  =   8;
   trf.n[4]  =   7;
   trf.n[5]  =   1;
# else /* [ cf. MOD_DEFLT == 1 ] */
   trf.n[1]  =  20;
   trf.n[2]  =  50;
   trf.n[3]  =   8;
   trf.n[4]  =   5;
   trf.n[5]  =   1;
# endif
/*............................................................................*/

   return;
}
/*============================================================================*/

short rvise_divsns( void )
{ 
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */
/*----------------------------------------------------------------------------*/
/* operation marks: */

   if ( trf.ref == ONE ) /* reference line */
      blc.m[0] = 1;
   else
      blc.m[0] = 1;

/* dependencies: */

   trf.n[6] = 2*trf.n[3] + trf.n[4];
   trf.n[7] = 2*trf.n[3];

   return ONE;
}
/*----------------------------------------------------------------------------*/
# ifndef CELSIUS_TO_KELVIN
   # define CELSIUS_TO_KELVIN ( 273.15 )
# endif
/*============================================================================*/

void deflt_params( void ) /* default geometrical and operating parameters     */
{                         /* DSC model: circular waveguide dielectric taper   */
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      ll = LINLEN;

/* prototypes: */

   void \
      cpylne( char lne[], const char *txt, const char *brc, short linlen );
/*----------------------------------------------------------------------------*/
/* lengths:                                                                   */
/*                                                                            */
/* write only CONNECTED (!) STRINGS in the string copy function strcpy(*) !!! */
/* [ 2nd. argument       |  ]                                                 */
/*                       V                                                    */
/* parameter strings:                                                         */
                                                                               
   strcpy( trf.stx[0], "Parameters" );

   cpylne( trf.stx[1], "channel_depth", "m", ll );
   cpylne( trf.stx[2], "channel_width", "m", ll );
   cpylne( trf.stx[3], "inflow_length", "m", ll );
   cpylne( trf.stx[4], "outflow_length", "m", ll );
   cpylne( trf.stx[5], "board_width", "m", ll );
   cpylne( trf.stx[6], "board_thickness", "m", ll );
   cpylne( trf.stx[7], "board_channel_angle", "DEG", ll );

# if DSC_FLDMDE != 0
   if ( TEMPGGE == 0 )
      cpylne( trf.stx[8], "temperature", "K", ll );
   else
      cpylne( trf.stx[8], "temperature", "deg_C", ll );
# else
   cpylne( trf.stx[8], "temperature", "inactive", ll );
# endif

   cpylne( trf.stx[9], "mean_density,_fluid","Kg/m^3", ll );
   cpylne( trf.stx[10], "thermal_expansion_coefficient,_fluid","1/K", ll );
   cpylne( trf.stx[11], "adiabatic_compression_coefficient,_fluid","1/Pa", ll );
   cpylne( trf.stx[12], "heat_conductivity,_fluid","W/(m*K)", ll );
   cpylne( trf.stx[13], "heat_capacity,_fluid","J/(Kg*K)", ll );
   cpylne( trf.stx[14], "dynamic_viscosity,_fluid","Kg/(m*s)", ll );
   cpylne( trf.stx[15], "gravitational_acceleration","m/sec^2;_fixed", ll );
   cpylne( trf.stx[16], "dissipation_time","sec", ll );
   cpylne( trf.stx[17], "characteristic_length","m", ll );
   cpylne( trf.stx[18], "start_fluid_dynamic_computations_with_that_delay",
      "sec", ll );
   cpylne( trf.stx[19], "inflow_fluid_velocity", "m/s", ll );
   cpylne( trf.stx[20], "coarsening_period_[LES_filter]", "sec", ll );
/*............................................................................*/
# if DSC_HCRMDE != 0
   cpylne( trf.stx[21], "PROPOSED_heat&fluid_time_step_[s]",
      "later_requested", ll );
# else
   cpylne( trf.stx[21], "proposed_heat&fluid_time_step_[s]",
      "inactive", ll );
# endif
/*............................................................................*/
/* number of parameters: */

   trf.s[0] = 21;

/*.............................. defaults ....................................*/
# if DSC_FLDMDE != 0
   if ( TEMPGGE == 0 )
      trf.s[8] = 2.93000e+02; /* inflow temperature [Kelvin] */
   else
      trf.s[8] = 2.00000e+01; /* inflow temperature [Celsius] */
# endif
/*............................................................................*/
# if MOD_DEFLT == 1
/* fluid: Water at 20 deg C */
/* stable with fluid flow time step = 2.0e-05 s */
/* and LES coarsening period 1.00e-04 s */

   trf.s[1] = 2.00000e-02; /* channel depth */
   trf.s[2] = 5.00000e-01; /* channel width */
   trf.s[3] = 5.00000e-01; /* inflow length [m] */
   trf.s[4] = 1.50000e+00; /* outflow length [m] */
   trf.s[5] = 6.00000e-02; /* board length [m] */
   trf.s[6] = 1.00000e-02; /* board width [m] */
   trf.s[7] = 3.50000e+01; /* angle alpha [DEG] */

   trf.s[9] = 9.98200e+02; /* fluid mean density [Kg/m^3]*/
   trf.s[10] = 1.00000e-07; /* thermal expans. coeff, fluid [1/K]*/
   trf.s[11] = 0.00000e+00; /* adiabatic compression coeff, fluid [1/K]*/
   trf.s[12] = 5.98000e-01; /* heat conductivity, fluid [W/(m*K)] */
   trf.s[13] = 4.18220e+03; /* heat capacity, fluid [J/(Kg*K)] */
   trf.s[14] = 1.00200e-03; /* dynamic viscosity, fluid [Kg/(m*sec)]*/
   trf.s[15] = 0.00000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[16] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[17] = 0.00000e+00; /* characteristic length [m] */
   trf.s[18] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[19] = 3.00000e+01; /* inflow fluid velocity [m/s] */          /* !!! */
   trf.s[20] = 1.00000e-04; /* coarsening time period [s], may work with this */
   trf.s[21] = 2.00000e-05; /* [proposed] heat and fluid time step [s] */

# elif MOD_DEFLT == 2
/* fluid: Water at 20 deg C */
/* stable with fluid flow time step = 1.0e-06 s */
/* and LES coarsening period 1.00e-05 s */

   trf.s[1] = 2.00000e-02; /* channel depth */
   trf.s[2] = 3.00000e-01; /* channel width */
   trf.s[3] = 5.00000e-01; /* inflow length [m] */
   trf.s[4] = 1.50000e+00; /* outflow length [m] */
   trf.s[5] = 1.00000e-01; /* board length [m] */
   trf.s[6] = 1.00000e-02; /* board width [m] */
   trf.s[7] = 3.50000e+01; /* angle alpha [DEG] */

   trf.s[9] = 9.98200e+02; /* fluid mean density [Kg/m^3]*/
   trf.s[10] = 1.00000e-07; /* thermal expans. coeff, fluid [1/K]*/
   trf.s[11] = 0.00000e+00; /* adiabatic compression coeff, fluid [1/K]*/
   trf.s[12] = 5.98000e-01; /* heat conductivity, fluid [W/(m*K)] */
   trf.s[13] = 4.18220e+03; /* heat capacity, fluid [J/(Kg*K)] */
   trf.s[14] = 1.00200e-03; /* dynamic viscosity, fluid [Kg/(m*sec)]*/
   trf.s[15] = 0.00000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[16] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[17] = 0.00000e+00; /* characteristic length [m] */
   trf.s[18] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[19] = 1.00000e+02; /* inflow fluid velocity [m/s] */          /* !!! */
   trf.s[20] = 1.00000e-05; /* coarsening time period [s], may work with this */
   trf.s[21] = 1.00000e-06; /* [proposed] heat and fluid time step [s] */

# elif MOD_DEFLT == 3
/* fluid: Water at 20 deg C */
/* stable with fluid flow time step = 5.00e-06 s */
/* and LES coarsening period 2.00e-04 s */

   trf.s[1] = 2.00000e-02; /* channel depth */
   trf.s[2] = 3.00000e-01; /* channel width */
   trf.s[3] = 5.00000e-01; /* inflow length [m] */
   trf.s[4] = 1.50000e+00; /* outflow length [m] */
   trf.s[5] = 1.00000e-01; /* board length [m] */
   trf.s[6] = 1.00000e-02; /* board width [m] */
   trf.s[7] = 5.50000e+01; /* angle alpha [DEG] */

   trf.s[9] = 9.98200e+02; /* fluid mean density [Kg/m^3]*/
   trf.s[10] = 1.00000e-07; /* thermal expans. coeff, fluid [1/K]*/
   trf.s[11] = 0.00000e+00; /* adiabatic compression coeff, fluid [1/K]*/
   trf.s[12] = 5.98000e-01; /* heat conductivity, fluid [W/(m*K)] */
   trf.s[13] = 4.18220e+03; /* heat capacity, fluid [J/(Kg*K)] */
   trf.s[14] = 1.00200e-03; /* dynamic viscosity, fluid [Kg/(m*sec)]*/
   trf.s[15] = 0.00000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[16] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[17] = 0.00000e+00; /* characteristic length [m] */
   trf.s[18] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[19] = 7.00000e+01; /* inflow fluid velocity [m/s] */          /* !!! */
   trf.s[20] = 2.00000e-04; /* coarsening time period [s], may work with this */
   trf.s[21] = 5.00000e-06; /* [proposed] heat and fluid time step [s] */

# elif MOD_DEFLT == 4
/* fluid: Water at 20 deg C */
/* stable with fluid flow time step = 1.0e-06 s */
/* and LES coarsening period 1.00e-05 s */

   trf.s[1] = 1.00000e-03; /* channel depth */
   trf.s[2] = 4.00000e-02; /* channel width */
   trf.s[3] = 5.00000e-02; /* inflow length [m] */
   trf.s[4] = 2.50000e-01; /* outflow length [m] */
   trf.s[5] = 1.50000e-02; /* board length [m] */
   trf.s[6] = 1.00000e-03; /* board width [m] */
   trf.s[7] = 4.50000e+01; /* angle alpha [DEG] */

   trf.s[9] = 9.98200e+02; /* fluid mean density [Kg/m^3]*/
   trf.s[10] = 1.00000e-07; /* thermal expans. coeff, fluid [1/K]*/
   trf.s[11] = 0.00000e+00; /* adiabatic compression coeff, fluid [1/K]*/
   trf.s[12] = 5.98000e-01; /* heat conductivity, fluid [W/(m*K)] */
   trf.s[13] = 4.18220e+03; /* heat capacity, fluid [J/(Kg*K)] */
   trf.s[14] = 1.00200e-03; /* dynamic viscosity, fluid [Kg/(m*sec)]*/
   trf.s[15] = 0.00000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[16] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[17] = 0.00000e+00; /* characteristic length [m] */
   trf.s[18] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[19] = 5.00000e+00; /* inflow fluid velocity [m/s] */          /* !!! */
   trf.s[20] = 1.00000e-04; /* coarsening time period [s], may work with this */
   trf.s[21] = 5.00000e-06; /* [proposed] heat and fluid time step [s] */

# else /* [ cf. MOD_DEFLT == 1 ] */
/* fluid: Water at 20 deg C */
/* stable with fluid flow time step = 2.0e-05 s */
/* and LES coarsening period 1.00e-04 s */

   trf.s[1] = 2.00000e-02; /* channel depth */
   trf.s[2] = 5.00000e-01; /* channel width */
   trf.s[3] = 5.00000e-01; /* inflow length [m] */
   trf.s[4] = 1.50000e+00; /* outflow length [m] */
   trf.s[5] = 6.00000e-02; /* board length [m] */
   trf.s[6] = 1.00000e-02; /* board width [m] */
   trf.s[7] = 3.50000e+01; /* angle alpha [DEG] */

   trf.s[9] = 9.98200e+02; /* fluid mean density [Kg/m^3]*/
   trf.s[10] = 1.00000e-07; /* thermal expans. coeff, fluid [1/K]*/
   trf.s[11] = 0.00000e+00; /* adiabatic compression coeff, fluid [1/K]*/
   trf.s[12] = 5.98000e-01; /* heat conductivity, fluid [W/(m*K)] */
   trf.s[13] = 4.18220e+03; /* heat capacity, fluid [J/(Kg*K)] */
   trf.s[14] = 1.00200e-03; /* dynamic viscosity, fluid [Kg/(m*sec)]*/
   trf.s[15] = 0.00000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[16] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[17] = 0.00000e+00; /* characteristic length [m] */
   trf.s[18] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[19] = 3.00000e+01; /* inflow fluid velocity [m/s] */          /* !!! */
   trf.s[20] = 1.00000e-04; /* coarsening time period [s], may work with this */
   trf.s[21] = 2.00000e-05; /* [proposed] heat and fluid time step [s] */
# endif
/*............................................................................*/

   return;
}
/*============================================================================*/

short rvise_params( void )
{ 
/* allusions:

   extern FORMSTATE *spt;
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: [void] */
/* prototypes: */

   void \
      cpylne( char lne[], const char *txt, const char *brc, short linlen );
/*----------------------------------------------------------------------------*/
/* waveguide type: "strip", "coax", "rectangular", "circular", "elliptic",etc.*/

   strcpy( trf.wgtype, "coaxial" );
/*............................................................................*/
   strcpy(( spt->ppt->text ), ( spt->tpt->text ));
/*............................................................................*/
/* time/frequency domain: */

   strcpy(( spt->ppt->domain ), "frequency_domain" );

   trf.fr = 1.e+08; /* not needed [yet to please the Maxwell field solver] */ 
   ( spt->ppt->fr ) = trf.fr;

   trf.b = trf.s[1];
   trf.a = trf.s[2];

   trf.alfa[1] = PI*trf.s[7]/180.;
   trf.alfa[2] = .5*PI - trf.alfa[1];
/*............................................................................*/
# if DSC_HCRMDE != 0
/* thermal parameters */ 
   ( spt->hcrstart ) = ZERO;

   if (( spt->hcdt ) < 1.e-277 )
      ( spt->hcdt ) = trf.s[21]; /* a proposed time step [ default ] */

/* thermal parameters: */ 
   trf.temp = trf.s[8]; 
/*
*//* trf.cv0: heat capacity [J/(K*m^3) !!!]
*//* trf.kh0: heat conductivity [W/(K*m)]
*//*
*/
   trf.kh0 = trf.s[12];
   trf.cv0 = trf.s[13]*trf.s[9]; /* [ heat capacity per volume ] */
/*............................................................................*/
# if DSC_FLDMDE != 0
   
/* coarsening period [ fluid component no 1 ]: */   

   trf.crsdt = trf.s[20];
   ( spt->ppt->fcp->crsdt[ONE] ) = trf.crsdt;

/* fluid parameters */ 
/*
*//* trf.tm0: mean temperature [Celsius] 
*//* trf.rm0: mean density [Kg/(m^3)]
*//* trf.bm0: thermal expansion coefficient [1/K]
*//* trf.cm0: [ adiabatic ] compression coefficient [1/Pa]
*//* trf.ny0: dynamic viscosity [Kg/(m*sec)]
*//* trf.q1: Cp/Cv - 1 [dimensionless]
*//* trf.td: dissipation time constant [sec]
*//* trf.gr: gravitation acceleration [m/(sec^2)]
*//* trf.LL: characteristic length [ for turbulence models - here not used ]
*/
/* A delay between the starts of heat diffusive and advective */
/* [ fluid dynamic ] operations may speed up the algorithm */

   ( spt->fldstart ) = trf.s[17];

   trf.tm0 = trf.s[8];
   trf.rm0 = trf.s[9];
   trf.bm0 = trf.s[10];
   trf.cm0 = trf.s[11]; 
   trf.ny0 = trf.s[14];

   trf.q1 = 0.4;

   trf.gr = trf.s[15];

/* dissipation time: */
   trf.td = trf.s[16];

/* characteristic length: */
   trf.LL = trf.s[17];

/* inflow fluid velocity: */
   trf.uf = trf.s[19];

/* outflow fluid velocity: */
   trf.vf = trf.uf;

/* Nusselt number at no-slip boundary faces */
   trf.nus = 1.;

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   return ONE;
}
/*============================================================================*/

void set_z_bases( )
{ 
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: void */
/* prototypes: void */
/*----------------------------------------------------------------------------*/

   blc.z[0] = ZERO;
   blc.z[1] = blc.z[0] + trf.s[1]; /* channel depth */

   return;
}
/*============================================================================*/

/*******************************************************************************
*                                                                              *
*   DANSE release 1.0r3 ]                                                      *
*                                                                              *
*   Point definition function cords(*)                                         *
*   [ DSC model prototype: 'flow aigainst inclined board' ]                  *
*                                                                              *
*   This function fixes the coordinates of the model [ geometry ] supporting   *
*   points pnt[n][]  ( n = 0,...,MOD_POINTS-ONE ) of the DSC model, viz.       *
*                                                                              *
*                pnt[n][0] = x, pnt[n][1] = y, pnt[n][2] = z,                  *
*                                                                              *
*   which are used to build up the block structure [ in function blocks(*) ].  *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define BLC_P_POINTS  0        /* BLC_P_POINTS = 1: print point coordinates  */
/*============================================================================*/

short cords( short layer )
{
/* allusions: */
/*   
   extern FORMSTATE *spt;
   extern BLOCSTR blc;
   extern struct transfer trf;
   extern double pnt[MOD_POINTS][THREE];
*/
/* declarations: */
/* canonical variables: */

   static short
      ii = null,
      jj = null, 
      points = null;

   static double
      uu = ZERO;
/*............................................................................*/
   static double 
      xx0 = ZERO,
      xx1 = ZERO,
      xx2 = ZERO,
      xx3 = ZERO,
      xx4 = ZERO,
      xx5 = ZERO,
      xx6 = ZERO,
      xx7 = ZERO,
      xx8 = ZERO,
      xx9 = ZERO,
      xx10 = ZERO,
      xx11 = ZERO,
      xx12 = ZERO,
      xx13 = ZERO,
      xx14 = ZERO,

      yy0 = ZERO,
      yy1 = ZERO,
      yy2 = ZERO,
      yy3 = ZERO,
      yy4 = ZERO,
      yy5 = ZERO,
      yy6 = ZERO,
      yy7 = ZERO,
      yy8 = ZERO,
      yy9 = ZERO,
      yy10 = ZERO,
      yy11 = ZERO,
      yy12 = ZERO,
      yy13 = ZERO,
      yy14 = ZERO,

      zz0 = ZERO;

   static double
      ss1 = ZERO,
      cc1 = ZERO,
      ss2 = ZERO,
      cc2 = ZERO;
/*............................................................................*/
/* prototypes: */

   double 
      sin( double ),
      cos( double ),
      sqrt( double );      
/*---------------------------------------------------------------------------*/
/* This part is canonical: take care with modifying it [ better you leave it */
/* untouched ] */

   uu = ZERO;

   jj = null;
   ii = null;
   while ( ii < layer )
   {
      while (( ii == blc.base[jj+ONE] )&&( jj < blc.m[null] ))
      {
         jj += ONE;
      };
      uu += blc.dz[jj];
      ii++;
   };
/*............................................................................*/
/* radii: --- */

/*............................................................................*/
/* bow [ arc] dependent quantities: --- */

   ss1 = sin( trf.alfa[1] );
   cc1 = cos( trf.alfa[1] );
   ss2 = sin( trf.alfa[2] );
   cc2 = cos( trf.alfa[2] );
/*............................................................................*/
/* xy-coordinates: */

   xx0 = ZERO;   
   yy0 = ZERO;

   xx1 = -.5*trf.s[5]*cc1 + .5*trf.s[6]*cc2;
   yy1 = +.5*trf.s[5]*ss1 + .5*trf.s[6]*ss2;

   xx2 = xx1 - trf.s[6]*cc2;
   yy2 = yy1 - trf.s[6]*ss2;

   xx3 = - xx1;
   yy3 = - yy1;

   xx4 = - xx2;
   yy4 = - yy2;

   xx5 = -.5*trf.s[2]*cc1 - .5*trf.s[6]/ss1;
   yy5 = +.5*trf.s[2];

   xx6 = - xx5;
   yy6 = - yy5;

   xx7 = xx6 - trf.s[6]/ss1;
   yy7 = yy6;

   xx8 = - xx7;
   yy8 = yy5;

   xx9 = xx5;
   yy9 = yy6;

   xx10 = xx6;
   yy10 = yy5;

   xx11 = xx5 - trf.s[3];
   yy11 = yy5;

   xx12 = xx11;
   yy12 = yy6;

   xx13 = xx6 + trf.s[4];
   yy13 = yy6;

   xx14 = xx13;
   yy14 = yy5;
/*............................................................................*/
/* z-coordinates: */

   zz0 = uu;
/*............................................................................*/
/* The total number of points: */

   points = 14;
/*............................................................................*/
/* pnt[][0]   =   x           pnt[][1]   =   y           pnt[][2]   =   z     */

   pnt[0][0]  = xx0;          pnt[0][1]  = yy0;          pnt[0][2]  = zz0;
   pnt[1][0]  = xx1;          pnt[1][1]  = yy1;          pnt[1][2]  = zz0;
   pnt[2][0]  = xx2;          pnt[2][1]  = yy2;          pnt[2][2]  = zz0;
   pnt[3][0]  = xx3;          pnt[3][1]  = yy3;          pnt[3][2]  = zz0;
   pnt[4][0]  = xx4;          pnt[4][1]  = yy4;          pnt[4][2]  = zz0;
   pnt[5][0]  = xx5;          pnt[5][1]  = yy5;          pnt[5][2]  = zz0;
   pnt[6][0]  = xx6;          pnt[6][1]  = yy6;          pnt[6][2]  = zz0;
   pnt[7][0]  = xx7;          pnt[7][1]  = yy7;          pnt[7][2]  = zz0;
   pnt[8][0]  = xx8;          pnt[8][1]  = yy8;          pnt[8][2]  = zz0;
   pnt[9][0]  = xx9;          pnt[9][1]  = yy9;          pnt[9][2]  = zz0;
   pnt[10][0] = xx10;         pnt[10][1] = yy10;         pnt[10][2] = zz0;
   pnt[11][0] = xx11;         pnt[11][1] = yy11;         pnt[11][2] = zz0;
   pnt[12][0] = xx12;         pnt[12][1] = yy12;         pnt[12][2] = zz0;
   pnt[13][0] = xx13;         pnt[13][1] = yy13;         pnt[13][2] = zz0;
   pnt[14][0] = xx14;         pnt[14][1] = yy14;         pnt[14][2] = zz0;

   return points; 
}
/************************ End of function cords(*)****************************/










/*******************************************************************************
*                                                                              *
*   DANSE release 1.0r3                                                        *
*                                                                              *
*   Block definition function blocks(*)                                        *
*   [ DSC model prototype: 'flow against inclined board' ]                   *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "./tools/blcutls.h" /* macros BLC_EXCITE(*), _BOUNDARY(*), */
                              /* BLC_TRIVIAL(*), and BLC_DECLARE(*)  */
# include "./tools/setmed.h"  /* media switching function setmed(*) */
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE != 0
# define SETMED(NN) \
{ \
   ( mdp->idx ) = (NN); \
 \
   if ((NN) == TWO ) \
   { \
      ( mdp->ci ) = lbl.ci; \
      ( mdp->cf ) = lbl.cf; \
 \
      ( mdp->cnn ) = ONE; \
 \
      ( mdp->eps ) = 1.; \
      ( mdp->myr ) = 1.; \
      ( mdp->ke ) = ZERO; \
      ( mdp->km ) = ZERO; \
 \
      ( mdp->kh ) = trf.kh0; \
      ( mdp->cv ) = trf.cv0; \
 \
      ( mdp->rm ) = trf.rm0; \
      ( mdp->tm ) = trf.tm0; \
      ( mdp->bm ) = trf.bm0; \
      ( mdp->cm ) = trf.cm0; \
      ( mdp->ny ) = trf.ny0; \
      ( mdp->q1 ) = trf.q1; \
      ( mdp->td ) = trf.td; \
      ( mdp->LL ) = trf.LL; \
      ( mdp->gr[0] ) = ZERO; \
      ( mdp->gr[0] ) = ZERO; \
      ( mdp->gr[1] ) = - trf.gr; \
      ( mdp->gr[2] ) = ZERO; \
      ( mdp->gp[0] ) = ZERO; \
      ( mdp->gp[1] ) = ZERO; \
      ( mdp->gp[2] ) = ZERO; \
   } \
   else if ((NN) == THREE ) \
   { \
      ( mdp->ci ) = lbl.ci; \
      ( mdp->cf ) = lbl.cf; \
      ( mdp->eps ) = ZERO; \
      ( mdp->myr ) = ZERO; \
      ( mdp->ke ) = ZERO; \
      ( mdp->km ) = ZERO; \
 \
      ( mdp->kh ) = trf.kh1; \
      ( mdp->cv ) = trf.cv1; \
 \
      ( mdp->rm ) = ZERO; \
      ( mdp->tm ) = ZERO; \
      ( mdp->bm ) = ZERO; \
      ( mdp->ny ) = ZERO; \
      ( mdp->q1 ) = ZERO; \
      ( mdp->td ) = ZERO; \
      ( mdp->LL ) = ZERO; \
      ( mdp->gr[0] ) = ZERO; \
      ( mdp->gr[1] ) = ZERO; \
      ( mdp->gr[2] ) = ZERO; \
      ( mdp->gp[0] ) = ZERO; \
      ( mdp->gp[1] ) = ZERO; \
      ( mdp->gp[2] ) = ZERO; \
   }; \
 \
   mdp = setmed( mdp ); \
}
# else /* if DSC_FLDMDE == 0 */
# define SETMED(NN) \
{ \
   ( mdp->idx ) = (NN); \
 \
   if ((NN) == TWO ) \
   { \
      ( mdp->ci ) = lbl.ci; \
      ( mdp->cf ) = lbl.cf; \
      ( mdp->eps ) = 1.; \
      ( mdp->myr ) = 1.; \
      ( mdp->ke ) = ZERO; \
      ( mdp->km ) = ZERO; \
 \
      ( mdp->kh ) = trf.kh0; \
      ( mdp->cv ) = trf.cv0; \
   } \
   else if ((NN) == THREE ) \
   { \
      ( mdp->ci ) = lbl.ci; \
      ( mdp->cf ) = lbl.cf; \
      ( mdp->eps ) = ZERO; \
      ( mdp->myr ) = ZERO; \
      ( mdp->ke ) = ZERO; \
      ( mdp->km ) = ZERO; \
 \
      ( mdp->kh ) = trf.kh1; \
      ( mdp->cv ) = trf.cv1; \
   }; \
 \
   mdp = setmed( mdp ); \
}
# endif /* DSC_FLDMDE == 0 */
# else /* if DSC_HCRMDE == 0 */
# define SETMED(NN) \
{ \
   ( mdp->idx ) = (NN); \
 \
   if ((NN) == TWO ) \
   { \
      ( mdp->ci ) = lbl.ci; \
      ( mdp->cf ) = lbl.cf; \
      ( mdp->eps ) = 1.; \
      ( mdp->myr ) = 1.; \
      ( mdp->ke ) = ZERO; \
      ( mdp->km ) = ZERO; \
   } \
   else if ((NN) == THREE ) \
   { \
      ( mdp->ci ) = lbl.ci; \
      ( mdp->cf ) = lbl.cf; \
      ( mdp->eps ) = ZERO; \
      ( mdp->myr ) = ZERO; \
      ( mdp->ke ) = ZERO; \
      ( mdp->km ) = ZERO; \
   }; \
 \
   mdp = setmed( mdp ); \
}
# endif /* DSC_HCRMDE == 0 */
/*===========================================================================*/

short blocks( short layer, char *option )
{
/*----------------------------------------------------------------------------*/
/* options != 'd'ivisions or 'p'oints        */
/* e.g. options 't'opology or 'c'oordinates: */
/*............................................................................*/
/* declarations: */

   static long
      ii = null,
      nn = null;

   static short 
      pp = null,
      jj = null;

   static const char
      e_wall = 'e';
/*
      m_wall = 'm';
*/
/* function prototypes: - */

   static MEDIUM
     *mdp = NULL;
     
/* function prototypes: - */

   MEDIUM
      *setmed( MEDIUM *mdp );

   BLC_DECLARE( ); /* part of declarations, cf.bldclr.h */

/* [ functions declared in blcdcl.h ]: */
/*
   void smedium( short medium_idx, long init_cell, long final_cell, \
                 double eps, double myr, double ke, double km, \
                 double kh, double cv, char *opt );

   BLC_EXCITE( char *type, long init_cell, short number_of_cells,
               short period, char face, double excitation_amplitude );
   BLC_BOUNDARY( char *type, long init_cell, short number_of_cells,
                 short period, char face, double permittivity );
   BLC_EVALUATE( char *type, long init_cell, short number_of_cells,
                 short period, char port, char *description );
   BLC_HCRBOUND( char *type, long init_cell, short number_of_cells,
                 short period, char face, double current_or_temperature );
   BLC_HCREVLTE( char *type, long init_cell, short number_of_cells,
                 short period, char face, char desription ["temperature_node",
                 e.g.] );
*/
/* function prototypes: - */
/* memory allocations: - */
/* ... */
/*----------------------------------------------------------------------------*/
   BLC_LIMITS( 'o' ); /* don't remove on modifying DSC model !!!              */
/*............................................................................*/
/* 
The following are dummy statements [ without computational impact ] that have
only been introduced to avoid subsidary warning messages from the C compiler
[option -Wall] recalling 'unused variables' [ a,b,...,r,d; e.g.]. The constants
referred to are usually declared in macro DECLARE( ) and are normally used in
the block connection macros Einlesen(*), CONNECT(*), etc. [ based on function
interface(*) ].
*/
   nn = forward; /* -> */
   nn = backward; /* <- */
   nn = a;
   nn = b;
   nn = c;
   nn = d;
   nn = ab;
   nn = bc;
   nn = ad;
   nn = dc;
   nn = ac;
   nn = ewall;
   nn = mwall;
   nn = e_wall;
/*............................................................................*/
   if( layer == null ) /* don't remove */
   {
      mdp = setmed( NULL );
      ( mdp->opt ) = *option;

      trf.exc = null;
      trf.bnd = null;
      trf.val = null;
# if DSC_HCRMDE != 0 
      trf.bhc = null;
      trf.ehc = null;
      trf.vhc = null;
# endif
   };
/*............................................................................*/
/*
    printf( "\n layer %d ", layer );
    scanf( "%s", ptr );
*/
/*.............................. waveguide blocks ............................*/
   BLOCK(1);

   blc.ab = trf.n[6];
   blc.bc = trf.n[1];

   CONNECT( 1, ab, 0, blc.ab,
            0, ewall, 0, 0 );

   CONNECT( 1, bc, 0, blc.bc,
            0, ewall, 0, 0 );

   CONNECT( 1, ad, 0, blc.bc,
            0, ewall, 0, 0 );
   
   if ( trf.c[2] == null )
   {
/* no-slip boundary conditions: */
      ii = lbl.m+ONE;
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "no_slip", ii, nn, pp, null, 1. );

      ii = lbl.m+blc.ab; 
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "no_slip", ii, nn, pp, ONE, 1. );
   }
   else
   {
/* free slip boundary conditions: */
      ii = lbl.m+ONE;
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "slip", ii, nn, pp, null, ZERO );

      ii = lbl.m+blc.ab; 
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "slip", ii, nn, pp, ONE, ZERO );
   };

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*blc.bc;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

/* inflow conditions: */

   ii = lbl.m + ONE;
   nn = blc.ab;
   pp = ONE;

   BLC_BOUNDARY( "inflow", ii, nn, pp, TWO, ( trf.s[19] ));

/* evaluate fluid velocity */

   ii = lbl.m + blc.ab/2;
   BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-inflow" );

   if ( *option == 'c' )
   {
      POINT( a, 11 );
      POINT( b, 12 );
      POINT( c, 9 );
      POINT( d, 5 );
   };

   QUDRL(1);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(2);

   blc.ab = trf.n[6];
   blc.bc = trf.n[2];

   CONNECT( 2, ab, 0, blc.ab,
            0, ewall, 0, 0 );

   CONNECT( 2, bc, 0, blc.bc,
            0, ewall, 0, 0 );

   CONNECT( 2, ad, 0, blc.bc,
            0, ewall, 0, 0 );
   
   if ( trf.c[2] == null )
   {
/* no-slip boundary conditions: */
      ii = lbl.m+ONE;
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "no_slip", ii, nn, pp, null, 1. );

      ii = lbl.m+blc.ab; 
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "no_slip", ii, nn, pp, ONE, 1. );
   }
   else
   {
/* free slip boundary conditions: */
      ii = lbl.m+ONE;
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "slip", ii, nn, pp, null, ZERO );

      ii = lbl.m+blc.ab; 
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "slip", ii, nn, pp, ONE, ZERO );
   };

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*blc.bc;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

/* outflow boundary conditions: */

   ii = lbl.m+ONE;
   nn = blc.ab;
   pp = ONE;
   BLC_BOUNDARY( "outflow", ii, nn, pp, TWO, trf.s[19] );

/* evaluate fluid velocity */

   ii = lbl.m+( blc.ab/2 );
   jj = null; do
   {
      BLC_EVALUATE( "un", ( ii+jj ), ONE, ONE, 'x', "x-outflow" );
      BLC_EVALUATE( "un", ( ii+jj ), ONE, ONE, 'y', "y-outflow" );
   } while (( ++jj ) < ONE );

   if ( *option == 'c' )
   {
      POINT( a, 13 );
      POINT( b, 14 );
      POINT( c, 10 );
      POINT( d, 6 );
   };

   QUDRL(2);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(3);

   blc.ab = trf.n[5];
   blc.bc = trf.n[3];

   CONNECT( 3, ab, 0, blc.ab,
            0, ewall, 0, 0 );

   CONNECT( 3, dc, 0, blc.ab,
            0, ewall, 0, 0 );

   CONNECT( 3, bc, 0, 0,
            1, dc, 0, forward );
   
/* no-slip boundary conditions at board: */

   ii = lbl.m+blc.ab*( blc.bc-ONE )+ONE;
   nn = blc.ab;
   pp = ONE;
   BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, 1. );

   if ( trf.c[2] == null )
   {
/* no-slip boundary conditions: */
      ii = lbl.m+ONE;
      nn = blc.ab;
      pp = ONE;
      BLC_BOUNDARY( "no_slip", ii, nn, pp, TWO, 1. );
   }
   else
   {
/* free slip boundary conditions: */
      ii = lbl.m+ONE;
      nn = blc.ab;
      pp = ONE;
      BLC_BOUNDARY( "slip", ii, nn, pp, TWO, ZERO );
   };

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*blc.bc;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

   if ( *option == 'c' )
   {
      POINT( a, 8 );
      POINT( c, 2 );
      POINT( d, 1 );
   };

   QUDRL(3);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(4);

   blc.ab = trf.n[5];
   blc.bc = trf.n[3];

   CONNECT( 4, ab, 0, blc.ab,
            0, ewall, 0, 0 );

   CONNECT( 4, dc, 0, blc.ab,
            0, ewall, 0, 0 );

   CONNECT( 4, bc, 0, 0,
            2, dc, 0, forward );
   
/* no-slip boundary conditions at board: */

   ii = lbl.m+blc.ab*( blc.bc-ONE )+ONE;
   nn = blc.ab;
   pp = ONE;
   BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, 1. );

   if ( trf.c[2] == null )
   {
/* no-slip boundary conditions: */
      ii = lbl.m+ONE;
      nn = blc.ab;
      pp = ONE;
      BLC_BOUNDARY( "no_slip", ii, nn, pp, TWO, 1. );
   }
   else
   {
/* free slip boundary conditions: */
      ii = lbl.m+ONE;
      nn = blc.ab;
      pp = ONE;
      BLC_BOUNDARY( "slip", ii, nn, pp, TWO, ZERO );
   };

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*blc.bc;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

   if ( *option == 'c' )
   {
      POINT( a, 7 );
      POINT( c, 4 );
      POINT( d, 3 );
   };

   QUDRL(4);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(5);

   blc.ab = trf.n[3];
   blc.bc = 2*blc.ab;

   CONNECT( 5, ab, null, blc.ab,
            1, dc, null, forward );

   CONNECT( 5, ac, 0, blc.ab,
            3, bc, 0, forward );

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*( blc.ab+ONE )/TWO;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

   TRNGL(5);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(6);

   blc.ab = trf.n[3];
   blc.bc = 2*blc.ab;

   CONNECT( 6, ab, null, blc.ab,
            2, dc, null, forward );

   CONNECT( 6, ac, 0, blc.ab,
            4, bc, 0, forward );

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*( blc.ab+ONE )/TWO;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

   TRNGL(6);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(7);

   blc.ab = trf.n[4];
   blc.bc = trf.n[7]; /* = 2*trf.n[3] */

   CONNECT( 7, dc, 0, blc.ab,
            0, ewall, 0, 0 );

   CONNECT( 7, ab, 0, blc.ab,
            1, dc, trf.n[3], forward );

   CONNECT( 7, ad, 0, blc.bc,
            5, bc, 0, forward );

   CONNECT( 7, dc, blc.ab, blc.ab,
            4, dc, 0, forward );

/* no-slip boundary conditions at board: */

   ii = lbl.m+blc.ab*( blc.bc-ONE )+ONE;
   nn = blc.ab;
   pp = ONE;
   BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, 1. );

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*blc.bc;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

   if ( *option == 'c' )
   {
      POINT( c, 3 );
   };

   QUDRL(7);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(8);

   blc.ab = trf.n[4];
   blc.bc = trf.n[7]; /* = 2*trf.n[3] */

   CONNECT( 8, dc, 0, blc.ab,
            0, ewall, 0, 0 );

   CONNECT( 8, ab, 0, blc.ab,
            2, dc, trf.n[3], forward );

   CONNECT( 8, ad, 0, blc.bc,
            6, bc, 0, forward );

   CONNECT( 8, dc, blc.ab, blc.ab,
            3, dc, 0, forward );

/* no-slip boundary conditions at board: */

   ii = lbl.m+blc.ab*( blc.bc-ONE )+ONE;
   nn = blc.ab;
   pp = ONE;
   BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, 1. );

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*blc.bc;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

   if ( *option == 'c' )
   {
      POINT( c, 1 );
   };

   QUDRL(8);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(9);

   blc.ab = trf.n[3];
   blc.bc = trf.n[7]; /* = 2*trf.n[3] */

   CONNECT( 9, bc, 0, blc.bc,
            0, ewall, 0, 0 );

   CONNECT( 9, ab, 0, blc.ab,
            1, dc, ( trf.n[3]+trf.n[4] ), forward );

   CONNECT( 9, ad, 0, blc.bc,
            7, bc, 0, forward );

   CONNECT( 9, dc, 0, blc.ab,
            4, ad, blc.ab, backward );

   if ( trf.c[2] == null )
   {
/* no-slip boundary conditions: */
      ii = lbl.m+blc.ab; 
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "no_slip", ii, nn, pp, ONE, 1. );
   }
   else
   {
/* free-slip boundary conditions: */
      ii = lbl.m+blc.ab; 
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "slip", ii, nn, pp, ONE, ZERO );
   };

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*blc.bc;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

   QUDRL(9);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   BLOCK(10);

   blc.ab = trf.n[3];
   blc.bc = trf.n[7]; /* = 2*trf.n[3] */

   CONNECT( 10, bc, 0, blc.bc,
             0, ewall, 0, 0 );

   CONNECT( 10, ab, 0, blc.ab,
             2, dc, ( trf.n[3]+trf.n[4] ), forward );

   CONNECT( 10, ad, 0, blc.bc,
             8, bc, 0, forward );

   CONNECT( 10, dc, 0, blc.ab,
             3, ad, blc.ab, backward );

   if ( trf.c[2] == null )
   {
/* no-slip boundary conditions: */
      ii = lbl.m+blc.ab;
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "no_slip", ii, nn, pp, ONE, 1. );
   }
   else
   {
/* free-slip boundary conditions: */
      ii = lbl.m+blc.ab;
      nn = blc.bc;
      pp = blc.ab;
      BLC_BOUNDARY( "slip", ii, nn, pp, ONE, ZERO );
   };

/* free slip boundary conditions at top and bottom: */

   ii = lbl.m+ONE;
   nn = blc.ab*blc.bc;
   pp = ONE;

   BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
   BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );

   QUDRL(10);

/* switch to fluid media parameter */
   SETMED(2);
/*............................................................................*/
   goto terminal; ;
/*............................................................................*/

  terminal:

   BLC_LIMITS( 'c' ); /* don't remove on modifying DSC model !!!              */

   return null;
} 
/***************** end of block definition function blocks(*) *****************/



/* [ function: systop ] */
/*******************************************************************************
*                                                                              *
*   DANSE - Discrete Approximation of the Navier Stokes Equations              *
*                                                                              *
*   Model dependent topological mesh generation function systop(*)             *
*                                                                              *
*   For any ( indexed ) family of points, labeled  ( one-to-one ) by           *
*   arbitrary  non-negative  integers , this  subroutine defines the           *
*   topological structure of a DSC mesh  by associating to each mesh           *
*   cell  ( indexed i = tpt->mi,...,tpt->mf )  the labels of its eight         *
*   vertex points , tpt->cm[i][k]  ( k=0,...,7 ) .                             *
*   The labeling scheme ( index k ) must be consistent with the mesh           *
*   point indexing conventions used in 'smatrx.h' .                            *
*   Also, electric and magnetic walls  may be defined in the present           *
*   subroutine  by  setting  the  neighbouring  mesh  identifiers of           *
*   adjacent ports , tpt->mn[i][j]  ( j=1,...,12 ) , to -1 or -2 ,             *
*   respectively.                                                              *
*                                                                              *
*   On the basis of the vertex points assignment fixed in 'systop(*)'          *
*   the complete DSC  network  ( connecting  neighbouring cells  and           *
*   ports ) is recovered in function 'linker.h'.                               *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define TPL_PRINT    0 /* TPL_PRINT = 1,2,...: printing support */
		        /* [ for model debugging purposes ] */

# define TPL_BOTWALLS 1 /* [don't:0] set [electric:1/magnetic:2] */
                        /* walls at bottom of structure */
# define TPL_TOPWALLS 1 /* [don't:0] set [electric:1/magnetic:2] */
                        /* walls at top of structure */

# define TPL_SANDWCHK 1 /* 1: perform sandwich check */
# define TPL_COORDCHK 0 /* 1: check coordinate periodicity */

# define ELECTRIC_WALL ( -1 )
# define MAGNETIC_WALL ( -2 )
/*----------------------------------------------------------------------------*/
# ifndef DSC_ADJCLL
   # define DSC_ADJCLL 1 
# endif
/*============================================================================*/

TOPOLOGY *
systop( FORMSTATE *state )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern TOPOLOGY *tpt;
   extern BLOCSTR blc;
   extern struct labels lbl;
*/
/* declarations: */

   static TOPOLOGY
     *tpt = null;

   static short
      ii = null,
      jj = null,
      kk = null,
      ind = null; 

   static long 
      hh = null,
      ll = null,
      mm = null,
      nn = null,
      pp = null,
      qq = null,
      rr = null;

   static char 
     *sysptr = DSC_MODEL;
/*
   static char 
     *ptr,
    **endp = null;
*/
/*............................................................................*/
# if TPL_COORDCHK == 1
   static double
      uu = ZERO,
      vv = ZERO,
      bound1 = ZERO;

  double 
      fabs( double x );
# endif
/*............................................................................*/
/* prototyping: */

   short 
      input( char *option ),
      cords( short layer ),
      blocks( short layer, char *option ),
      qudrl( short p, short q, char *option ),
      trngl( short p, char *option );
/*............................................................................*/
/* copy FORMSTATE [ state from formdrv ]: */

   tpt = ( state->tpt );
   spt = state;
/*............................................................................*/
/* copy identifiers_and_comments: */

   strncpy(( tpt->name ), sysptr, STS_SIZE );
   strncpy(( tpt->text ), "DSC_structure", STS_SIZE );
/*............................................................................*/
/* enter DSC operation parameters: */

/*........................................*/
   ind = input( "operations" );
/*........................................*/

   if ( ind == ONE )
   {
      printf( "\n\n Message from %s :", __func__ );
      printf( "\n Error on calling parameter input function input(*) "
         "[ option 'operations'] !!!\n" );
      exit( EXIT_FAILURE );
   }
   else if ( ind == -ONE ) 
   {
      ( tpt->rtn ) = -ONE;
      return tpt;
   };
/*............................................................................*/
/* enter DSC mesh divisions: */

/*........................................*/
   ind = input( "divisions" );
/*........................................*/

   if ( ind == ONE )
   {
      printf( "\n\n Message from %s :", __func__ );
      printf( "\n Error on calling parameter input function input(*) "
         "[ option 'divisions'] !!!\n" );
      exit( EXIT_FAILURE );
   }
   else if ( ind == -ONE ) 
   {
      ( tpt->rtn ) = -ONE;
      return tpt;
   };
/*............................................................................*/
/* base-layer indices: */

   blc.base[null] = null;
   ii = null; 
   while( ii < blc.m[null] )
   {  
      ii++;
      blc.base[ii] = blc.base[ii-ONE] + blc.m[ii];
   };
   while( ii < MOD_DOMAINS )
   {
      ii++;
      blc.base[ii] = blc.base[ii-ONE] + ONE;
   };
/*............................................................................*/
/* generating mesh topology: */

   ii = null;
   while ( ii <= MOD_BLOCKS )
   {
      blc.indnt[ii] = null;
      blc.leave[ii] = null;
      ii++;
   };
      
   lbl.mf = null;
   lbl.pf = null;

   blc.mm = null;
   blc.pp = null;

   ii = null;
   while ( ii < blc.m[null] ) /* blc.m[null] = number of domains */
   {
      blc.dmn = ii;   /* the actual domain index */
      blc.cov = null; /* cover flag [ set to ONE on a virtual cover layer ] */

      jj = blc.base[ii];
      while( jj <= blc.base[ii+ONE] )
      {
         blc.lbl = null; /* reset block counter to null */
         blc.lay = jj;   /* the actual layer index */
         lbl.m = blc.mm; /* the yet counted mesh cells */
         lbl.p = blc.pp; /* the yet counted vertex points */

         lbl.mi = null;  /* with the first call of qudrl(*), or trngl(*) */
         lbl.pi = null;  /* these numbers are set to ihe initial mesh cell */
                         /* index or vertex point, resp., of the actual layer */

         if ( jj == blc.base[ii+ONE] )
         {
            blc.cov = ONE;
/*............................................................................*/
# if BLC_P_BLIMTS == 2 
            printf( "\n\n base%03d, layer%03d "
               "[ cover - virtually removed ]: ", ( ii+ONE ), jj );
# endif
/*............................................................................*/
         };
/*............................................................................*/
         ind = clblcs( null, MOD_BLOCKS );    /* clear blc.border[][][]       */
                                             /*                               */
         ind = blocks( jj, "topology" );    /*                                */
/*........................................*/
         if ( ind == ONE )      /* Error */
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on block definition function blocks(*) "
               "[ option 'topology' ] !!!\n" );
            exit( EXIT_FAILURE );
         }
         else if ( ind == - ONE ) /* Escape */ 
         {
            ( tpt->rtn ) = -ONE;
            return tpt;
         };

         if (( jj == blc.base[ii] )&&( blc.cov == null ))
         {
            blc.dm[ii] = lbl.m  - lbl.mi + ONE;
            blc.dp[ii] = lbl.pf - lbl.pi + ONE;/* lbl.pf  may actually exceed */
         };                                    /* lbl.p  on  change of domain */
                                               /* h->h+1: jj = blc.base[h]    */
/*............................................................................*/
# if TPL_PRINT == 1
         if ( jj == blc.base[ii+ONE] )
         {
            printf( "\n\n domain:       %3d\n layer:        %3d", ii, jj );
            printf( "\n mesh period:  %3d\n point period: %3d\n", blc.dm[ii],
               blc.dp[ii] );
            printf( "\n Please acknowledge !"
               "\n [ enter any character ] >---> " );
            scanf( "%s", ptr );
         };
# endif
/*............................................................................*/
         hh = null ; 
         while(( ++hh ) <= blc.dm[ii] )
         {
            nn = blc.mm + hh;
            if ( blc.base[ii] < jj )          /* cover last layer by setting  */
            {                                 /* upper vertex points          */
               ll = blc.mm + hh - blc.dm[ii];
               kk = null; do
               {
                  ( tpt->cm[ll][kk+FOUR] ) = ( tpt->cm[nn][kk] );
               } while(( ++kk ) < FOUR );
            };
/*............................................................................*/
# if TPL_BOTWALLS == 1
            if ( jj == null )  
	    {
# if DSC_ADJCLL == 0 /* electric wall at face 4 */
	       ( tpt->mn[nn][4] ) = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1 /* electric wall at ports 3[+ONE] and 7[+ONE] */
	       ( tpt->mn[nn][3] ) = ELECTRIC_WALL;
	       ( tpt->mn[nn][6] ) = ELECTRIC_WALL;
# endif
	    };
# elif TPL_BOTWALLS == 2
            if ( jj == null )  
	    {
# if DSC_ADJCLL == 0 /* magnetic wall at face 4 */
	       ( tpt->mn[nn][4] ) = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1 /* magnetic wall at ports 3[+ONE] and 7[+ONE] */
	       ( tpt->mn[nn][3] ) = MAGNETIC_WALL;
	       ( tpt->mn[nn][6] ) = MAGNETIC_WALL;
# endif
	    };
# endif
/*............................................................................*/
# if TPL_TOPWALLS == 1
   	    if ( jj == blc.base[blc.m[null]] - ONE )
	    {
# if DSC_ADJCLL == 0 /* electric wall at face 5 */
	       ( tpt->mn[nn][5] ) = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1 /* electric wall at ports 1[+ONE] and 4[+ONE] */
	       ( tpt->mn[nn][1] ) = ELECTRIC_WALL;
	       ( tpt->mn[nn][4] ) = ELECTRIC_WALL;
# endif
	    };
# elif TPL_TOPWALLS == 2
   	    if ( jj == blc.base[blc.m[null]] - ONE )
	    {
# if DSC_ADJCLL == 0 /* magnetic wall at face 5 */
	       ( tpt->mn[nn][5] ) = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1 /* magnetic wall at ports 1[+ONE] and 4[+ONE] */
	       ( tpt->mn[nn][1] ) = MAGNETIC_WALL;
	       ( tpt->mn[nn][4] ) = MAGNETIC_WALL;
# endif
	    };
# endif
/*............................................................................*/
         }; /* while hh <= blc.dm[ii] */

         if ( jj < blc.base[ii+ONE] )
         {
            blc.mm += blc.dm[ii];
            blc.pp  = lbl.pi + blc.dp[ii] - ONE;
         }; /* end if jj < blc.base[ii+ONE] */
         jj++;
      }; /* while jj <= blc.base[ii+ONE] */
      ii++;
   }; /* while ii < blc.m[null] */

   ( tpt->mi ) = ONE;
   ( tpt->mf ) = blc.mm;
/*..........................................................................*/
# if TPL_SANDWCHK == 1
/* perform sandwich check: */
   printf( "\n\n %s: sandwich_check started ", __func__ );
   
   mm = null;
   ii = null;
   while ( ii < blc.m[null] )
   {
      jj = blc.base[ii];
      while( jj < blc.base[ii+ONE] )
      {
         hh = null;
         while (( ++hh ) <= blc.dm[ii] )
         {
            nn = mm + hh;

            kk = null; do
            {
               pp = ( tpt->cm[nn][kk] );
               qq = ( tpt->cm[nn][kk+FOUR] );
               rr = pp + blc.dp[ii]; /* periodic point */

               if ( qq != rr )
               {
                  printf( "\n\n Message from function %s :", __func__ ); 
                  printf( "\n\n DSC mesh sandwich check failed !!!" );   
                  printf( "\n [ Corner point periodicity condition hurt by "
                     "points no.%ld and %ld != %ld,\n   on cell no.%ld.]",
                     pp, qq, rr, nn );
                  ( tpt->rtn ) = ONE;
                  return tpt;
               };
/*............................................................................*/
# if TPL_COORDCHK == 1
/* periodic coordinates: */

               ll = null; do
               {
                  uu = ( spt->ppt->cpt->c[pp][ll] );
                  vv = ( spt->ppt->cpt->c[qq][ll] );
                  if ( bound1 < fabs( uu - vv ))
                  {
                     printf( "\n\n Message from function %s :", __func__ );
                     printf( "\n DSC mesh sandwich check failed !!!" );
                     printf( "\n [ Corner point periodicity condition hurt by "
                        "points no.%ld & %ld,\n   on cell no.%ld.]",
                        pp, qq, nn );

                     ( tpt->rtn ) = ONE;
                     return tpt;
                  };/* end if ... */
               } while(( ++ll ) < TWO );
# endif /* TPL_COORDCHK == 1 */
/*............................................................................*/
            } while(( ++kk ) < FOUR );
         };
         mm += blc.dm[ii];
         jj++;
      }; /* while ( jj < blc.base[ii+ONE] ) */
      ii++;
   }; /* while ( ii < blc.m[null] ) */
   
   printf( "\r %s: sandwich_check successfully terminated.", __func__ );
   
# endif /* TPL_SANDWCHK == 1 */
/*............................................................................*/
# if TPL_PRINT == 3
   printf( "\n\n tpt->mf = %d\n please acknowledge: ", ( tpt->mf ));
   scanf( "%s", ptr );
# endif
/*............................................................................*/
# if TPL_PRINT == 4
   for ( mm=( tpt->mi ); mm<=( tpt->mf ); mm++ )
   {
      printf( "\n mesh %d : (%d,%d,%d,%d,%d,%d,%d,%d)", mm,
         tpt->cm[mm][0], tpt->cm[mm][1], tpt->cm[mm][2], tpt->cm[mm][3],
         tpt->cm[mm][4], tpt->cm[mm][5], tpt->cm[mm][6], tpt->cm[mm][7] );
   };
   scanf( "%s", ptr );
# endif
/*............................................................................*/
   ( tpt->rtn ) = null;
   return tpt;
}  
/*==================== end of function body systop(*) ========================*/
# undef TPL_PRINT
# undef TPL_BOTWALLS
# undef TPL_TOPWALLS
# undef TPL_SANDWCHK
# undef TPL_COORDCHK
/*********************** end of function systop(*) ****************************/


/* [ function: syssmx ] */
/*******************************************************************************
*                                                                              *
*   DANSE - Discrete Approximation of the Navier Stokes Equations              *
*                                                                              *
*   Model dependent coordinates and media transfer function syssmx(*)          *
*                                                                              *
*   Given a topological DSC grid  by specifying ( in subroutine 'systop(*)')   *
*   vertex point  indices  tpt->cm[h][k]  ( k=0,...,7 )  for all  DSC meshes   *
*   (indexed h=tpt->mi,...,tpt->mf), this subroutine fixes the mesh dimensions *
*   by defining the vertex point coordinates ppt->cpt->c[i][j] ( j=0,1,2; in   *
*   x,y and z direction, respectively - index i runs  over  the  used vertex   *
*   point indices ).                                                           *
*                                                                              *
*   Also, the material mesh parameters, viz. the rel. permittivity and per-    *
*   meability tensors med.ep, med.my and the  electric and magnetic conduc-    *
*   tivity tensors, med.ke [Ohm^-1 m^-1], med.km [Ohm m^-1] , are specified    *
*   in this subroutine. Media index med.idx[nn] indicates a set of material    *
*   parameters ( med.ke[med.idx][] etc. ) characterizing cell number nn.       *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
# include "./math/linsec.h"
# include "./math/lotos.h"
*/
/*----------------------------------------------------------------------------*/
# define SMX_REPEATER 1 /* 1: enable s-matrix repeater function [ and take    */
                        /* care for not hurting the periodicity condition ! ] */
/* repeater function limited to domains D: SMX_BOTLIMIT < D < SMX_TOPLIMIT */
# define SMX_TOPLIMIT 2 
# define SMX_BOTLIMIT ( -1 )
/*----------------------------------------------------------------------------*/
# define SMX_SANDWCHK  1   /* 1: enable sandwich check                        */
# if SMX_SANDWCHK == 1
   # define SMX_COORDCHK 0 /* 1: enable coordinate periodicity check          */
# endif                    /* [ coordinate check requires sandwich check ]    */
/*----------------------------------------------------------------------------*/
# define SMX_POINTLOG  1   /* 1: enable point coordinate logfile, pnt[][]     */
# define SMX_POINTPLT  1   /* N: enable plot routine for points pnt[][] up to */
                           /* layer N-1 [ disable with SMX_POINTPLT 0 ]       */
/*----------------------------------------------------------------------------*/
# define SMX_SYSPLT    1   /* 1: enable sysplot function                      */
# if SMX_SYSPLT == 1
   # include "./former/sspltp.h"
# endif
/*============================================================================*/

PARAMETERS *
syssmx( FORMSTATE *state )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern PARAMETERS *ppt;
   extern BLOCSTR blc;
   extern struct labels lbl;
*/
/* declarations: */

   static PARAMETERS
     *ppt = null;

   static TOPOLOGY
     *tpt;

   static TXCNSL
     *csp;

   static char
      ptr[STS_SIZE] = {null},
      ctmptr[STS_SIZE] = {null},
      tmestr[STS_SIZE] = {null};

   static const char
     *timefrm = "created: %.24s ",
     *sysptr = DSC_MODEL;
/*
   static char **endp = null;
*/
   static long
     hh = null,
     mm = null,
     nn = null;

   static short
      ii = null,
      jj = null,
      kk = null,
     ind = null;
/*............................................................................*/
# if SMX_SANDWCHK == 1
   static long
      pp = null,
      qq = null,
      rr = null;

   # if SMX_COORDCHK == 1
      static signed char
         ll = null;

      static double
         uu = ZERO,
         vv = ZERO,
         bound1 = 1.e-12;
   # endif
# endif
/*............................................................................*/
# if SMX_POINTLOG == 1
   static time_t
      nseconds = null,
     *timer = null;

   static char
     *pntlog = "pnt.log",
      logptr[STS_SIZE] = {null};

   static FILE
     *pointlog;
# endif
/*............................................................................*/
# if SMX_POINTPLT != 0
   static double
      xx = ZERO,
      yy = ZERO,
      xmin = ZERO,
      xmax = ZERO,
      ymin = ZERO,
      ymax = ZERO;

   static char
     *pntplt = "gnu.pnt",
     *pntdat = ".dat",
      pltfle[STS_SIZE] = {null},
      datfle[STS_SIZE] = {null};

   static FILE
     *pointplt,
     *pointdat;
# endif
/*............................................................................*/
# if SMX_SYSPLT == 1
   static SYSPLT 
     *ssp;
# endif
/*............................................................................*/
/* prototypes: */

   time_t time( time_t *timer );
   char *ctime( const time_t *timer );

   double 
      cos( double x ),
      sin( double x ),
      tan( double x ),
      fabs( double x );

   int linsec( double a1 , double a2 , double r1 , double r2 );

   TOPOLOGY
     *systop( FORMSTATE *state );

   char
     *lotos( long n, char m );

   short 
      input( char *option ),
      cords( short layer ),
      blocks( short layer, char *option ),
      qudrl( short p, short q, char *option ),
      trngl( short p, char *option );

   TXCNSL
     *txcnsl( TXCNSL *csp );
/*............................................................................*/
# if SMX_SYSPLT == 1
   SYSPLT
     *sysplt( SYSPLT *ssp );
# endif
/*----------------------------------------------------------------------------*/
/* copy FORMSTATE [ state from formdrv: mesh topology information, e.g.]: */

   spt = state;
/*............................................................................*/
/*    systop( spt );     *//* usually done by calling process 'formdrv(*)' */
/*......................*/

   tpt = spt->tpt;
   ppt = spt->ppt;
/*............................................................................*/
/* copy identifiers_and_comments: */

   strncpy(( ppt->name ), sysptr, STS_SIZE );
   strncpy(( ppt->text ), ( tpt->text ), STS_SIZE );
/*............................................................................*/
/* initialize text console */

   csp = txcnsl( null );
/*............................................................................*/
# if SMX_SYSPLT == 1
/* initialize plot function sysplt(*): */

   ssp = sysplt( null ); 
/*............................................................................*/
/* configure sysplt(*): */

   ( ssp->tpt ) = tpt;
   ( ssp->ppt ) = ppt;
    strcpy(( ssp->flbl ), ( spt->flbl ));

   ( ssp->media ) = 1; /* N: that number of media shall be plotted */

/*   For N = number of media, define media indices mdidx[0],...,mdidx[N-1] 
*//*
*//* ( ssp->mdidx[0] ) = 0; K: media index,
*//* for example:
*//* K=1 for standard medium [ vacuum ],
*//* K=0 for trivial cells,
*//* [trivial (K=0) requires macro SPL_PLOTTRIVL in file FORMER.CONF set to 1]
*//*
*//* ( ssp->mdidx[0] ) = 2; K: media index [ K=2 for any Epsr1 != 1., e.g.]
*//* ( ssp->mdidx[1] ) = 3; L: media index [ L=3 for any Epsr2 != 1., etc.]
*//* ( ssp->mdidx[2] ) = 4; M: media index [ M=4 ... ]
*/
   ( ssp->mdidx[0] ) = 2; /* fluid dielectric */

# endif /* SMX_SYSPLT == 1 */
/*............................................................................*/
/* call parameters: */
/*............................................................................*/
   ind = input( "parameters" );   /* call parameters                          */
/*..............................*/

   if ( ind == ONE )
   {
      printf( "\n\n Message from %s :", DO_SYSSMX );
      printf( "\n Error on input function input(*) "
         "[ option 'parameters'] !!!\n " );
      exit( EXIT_FAILURE );
   }
   else if ( ind == -ONE ) 
   {
      ppt->rtn = -ONE;
      return ppt;
   };
/*............................................................................*/
/* z - increments: */

   ii = null;
   while ( ii < blc.m[null] )
   {
      jj = ii + ONE;
      if( null < blc.m[jj] )
         blc.dz[ii] =  ( blc.z[jj] - blc.z[ii] ) / blc.m[jj];
      else
         blc.dz[ii] = ZERO;
      ii++;
   };
/*............................................................................*/
   # if SMX_SYSPLT == 1
/* sysplot: */
   printf( "\n\n Create plot files for selected layers:\n" );
   printf( "\n Please enter point layers to be plotted "
      " [ Escape/end: enter -ONE ]\n" );

   ( csp->dflng ) = -ONE; /* default */

   jj = blc.base[blc.m[null]];
   ii = null;
   do
   {
      ii++ ;
      strcpy( ptr, lotos( ii, null ));
      strcpy( csp->rqlng, "Enter " );
      strcat( csp->rqlng, ptr );
      strcat( csp->rqlng, ". layer" );
/*............................................................................*/
      csp = txcnsl( csp );        /* input on text console                    */
/*..............................*/
      ( csp->dflng ) = ( csp->inlng ) + ONE;; /* next default */
      ( ssp->lay[ii] ) = ( csp->inlng );
   }  while (( null <= ( csp->inlng ))&&(( ssp->lay[ii] ) < jj )&&( ii <= jj ));
   if ( csp->inlng == -ONE )
     ii--;
   ( ssp->lay[null] ) = ii;
# endif /* SMX_SYSPLT == 1 */
/*............................................................................*/
# if SMX_POINTLOG == 1
/* supporting point logfile: */

   strcpy ( logptr, pntlog );
   pointlog = fopen( logptr, "w" );
/*
   printf( "\n opened: Supporting points logfile '%s'", logptr );
*/
   fprintf( pointlog, "%s\n\n", "SYSSMX_log_file: "
      "support point coordinates");
   fprintf( pointlog, "%s\n", ( ppt->name ));
   fprintf( pointlog, "%s\n", ( ppt->text ));
   fprintf( pointlog, "point_layers________: %3d \n",
      blc.base[blc.m[null]]+ONE );
# endif
/*............................................................................*/
/* generating mesh geometry: */

   ii = null;
   while ( ii <= MOD_BLOCKS )
   {
      blc.indnt[ii] = null;
      blc.leave[ii] = null;
      ii++;
   };
/*............................................................................*/
/* 'trivial' medium [ indexed null ]: */

   ( ppt->mpt->idx[null] ) = null; /* med.idx[null]: maximum index */
                                   /* of predefined media */
   ( ppt->mpt->tg[null] ) = ZERO; /* ppt in struct FORMSTATE points to struct */
/*............................................................................*/
# if DSC_HCRMDE != 0
   ( ppt->mpt->cv[null] ) = ZERO; /* PARAMETERS, and therein mpt points to */
   ( ppt->mpt->kh[null] ) = ZERO; /* struct media */
# endif
/*............................................................................*/

   kk = null;
   do
   {
      ppt->mpt->ep[null][kk] = ZERO;
      ppt->mpt->my[null][kk] = ZERO;
      ppt->mpt->ke[null][kk] = ZERO;
      ppt->mpt->km[null][kk] = ZERO;
      ppt->mpt->ms[null][kk] = ZERO;
      ppt->mpt->hg[null][kk] = ZERO;
      kk++;
   }  while ( kk < DIMNS );
   do
   {
      ppt->mpt->ep[null][kk] = ZERO;
      ppt->mpt->my[null][kk] = ZERO;
      ppt->mpt->ke[null][kk] = ZERO;
      ppt->mpt->km[null][kk] = ZERO;
      kk++;
   }  while ( kk < SIX );

/* canonical medium [ vacuum, indexed ONE ]: */

   ( ppt->mpt->idx[null] ) = ONE;
   ( ppt->mpt->tg[ONE] ) = ZERO;

# if DSC_HCRMDE != 0
   ( ppt->mpt->cv[ONE] ) = ZERO;
   ( ppt->mpt->kh[ONE] ) = ZERO;
# endif

   kk = null;
   do
   {
      ppt->mpt->ep[ONE][kk] = ONE;
      ppt->mpt->my[ONE][kk] = ONE;
      ppt->mpt->ke[ONE][kk] = ZERO;
      ppt->mpt->km[ONE][kk] = ZERO;
      ppt->mpt->ms[ONE][kk] = ZERO;
      ppt->mpt->hg[ONE][kk] = ZERO;
      kk++;
   }  while ( kk < DIMNS );
   do
   {
      ppt->mpt->ep[ONE][kk] = ZERO;
      ppt->mpt->my[ONE][kk] = ZERO;
      ppt->mpt->ke[ONE][kk] = ZERO;
      ppt->mpt->km[ONE][kk] = ZERO;
      kk++;
   }  while ( kk < SIX );
/*............................................................................*/
/* reset s-matrix repeater flags: */

   hh = null; /* tpt: in struct FORMSTATE pointer to struct TOPOLOGY */
   do
   {
      ( ppt->rep->smx[hh] ) = null;
/*............................................................................*/
# if DSC_HCRMDE != null
      ( ppt->rhp->smx[hh] ) = null;
# endif
/*............................................................................*/
   } while (( hh++ ) < ( tpt->mf ));
/*............................................................................*/
/* reset cell and vertex point counters */

   lbl.mf = null;
   lbl.pf = null;

   blc.mm = null;
   blc.pp = null;

   blc.dm[null] = ( tpt->mf ) - ( tpt->mi ) + ONE;
/*............................................................................*/
/* blc.dm[N]: eventually set to the number of mesh cells ['cell period'] in   */
/* any layer of domain N [ analogous to blc.dp[N] for vertex points ] */

   ii = null;
   while ( ii < blc.m[null] ) /* blc.m[null] = number of domains */
   {
      blc.dmn = ii;   /* the actual domain index */
      blc.cov = null; /* cover flag [ set to ONE on a virtual cover layer ] */

      jj = blc.base[ii];
      while( jj <= blc.base[ii+ONE] )
      {
         blc.lay = jj;   /* the actual layer index */
         lbl.m = blc.mm; /* the yet counted mesh cells */
         lbl.p = blc.pp; /* the yet counted vertex points */

         lbl.mi = null;  /* with the first call of qudrl(*), or trngl(*) */
         lbl.pi = null;  /* these numbers are set to ihe initial mesh cell */
                         /* index or vertex point, resp., of the actual layer */

         if ( jj < blc.base[ii+ONE] ) /* initialize material tensors */
         {
            nn = blc.mm + ONE;
            mm = blc.mm + blc.dm[ii]; 
            while( nn <= mm )
            {
               ( ppt->mpt->idx[nn] ) = ONE; /* default: canonical medium */
               nn++;                        /* [vacuum] */
            }; /* while nn <= ... */
         }
         else if ( jj == blc.base[blc.m[null]] )
         {
            nn = blc.mm + ONE;
            mm = blc.mm + blc.dm[ii];
            while( nn <= mm )
            {
               ( ppt->mpt->idx[nn] ) = null; /* set trivial medium on */
               nn++;                         /* last cover [ without */
            }; /* while nn <= ... */         /* computational effect ] */
         };                                  /* - Only for sysplt(*) function */

         if ( jj == blc.base[ii+ONE] )
         {
            blc.cov = ONE;
/*............................................................................*/
# if SMX_PRINT == 1
            printf( "\n\n base%03d, layer%03d "
                    "[ cover - virtually removed ]: ", ii+ONE, jj );
# endif
/*............................................................................*/
         };

/*............................................................................*/
         ind = cords( jj );      /* define supporting points, increments etc. */
/*.............................*/

         if ( MOD_POINTS < ind )
         {
            printf( "\n\n Message from function %s :", DO_SYSSMX );
            printf( "\n\n Too many support points !!! " );
            printf( "\n [ The maximum number is %d = macro MOD_POINTS "
                    " ", MOD_POINTS );
            printf( "\n   in function '%s'. ] ", DO_SYSSMX );
            printf( "\n   - Change macro in compliance with memory"
                    " resources. \n" );
            exit( EXIT_FAILURE );
         };
/*............................................................................*/
# if SMX_POINTLOG == 1
/* supporting point logfile:    */

         if (( jj < blc.base[ii+ONE] )||( jj == blc.base[blc.m[null]] ))
         {
            fprintf( pointlog, "\nlayer_______________: %3d", jj );
            fprintf( pointlog, "\nnumber_of_points____: %3d\n\n", ind );

            for ( kk = null; kk<ind; kk++ )
            {
               fprintf( pointlog, "point%03d: %+.15e %+.15e %+.15e\n", kk,
                  pnt[kk][0], pnt[kk][1], pnt[kk][2] );
            };
         };
# endif
/*............................................................................*/
# if SMX_POINTPLT != 0
         if ( jj < SMX_POINTPLT )
	 {
            xmin = 1.0e+77;
            xmax = -xmin;
            ymin = xmin;
            ymax = xmax;
	    
            for ( kk = null; kk<ind; kk++ )
            {
               if( pnt[kk][0] < xmin )
	          xmin = pnt[kk][0];
               if( xmax < pnt[kk][0] )
	          xmax = pnt[kk][0];
               if( pnt[kk][1] < ymin )
	          ymin = pnt[kk][1];
               if( ymax < pnt[kk][1] )
	          ymax = pnt[kk][1];
            };
            xx = xmax - xmin;
            yy = ymax - ymin;

	    if ( xx < yy )
	       xx = yy;
	    else if ( yy < xx )
	       yy = xx;

            xx = .5*( xmin + xmax );
            xmin = xx - .55*yy;
            xmax = xx + .55*yy;
            xx = .5*( ymin + ymax );
            ymin = xx - .55*yy;
            ymax = xx + .55*yy;

            strcpy ( pltfle, pntplt );
            strcat ( pltfle, ( lotos( jj, null )));
            pointplt = fopen( pltfle, "w" );
            strcpy ( datfle, pltfle );
            strcat ( datfle, pntdat );
            pointdat = fopen( datfle, "w" );

            fprintf( pointplt, "set title '%s %s'\n",
               ( ppt->name ), ( ppt->text ));
            fprintf( pointplt, "set xrange [%+.15e:%+.15e]\n", xmin, xmax );
            fprintf( pointplt, "set yrange [%+.15e:%+.15e]\n", ymin, ymax );
            fprintf( pointplt, "set xlabel '%s'\n", "m" );
            fprintf( pointplt, "set ylabel '%s'\n", "m" );
            fprintf( pointplt, "set size square\n" );
            fprintf( pointplt, "set grid\n" );
            fprintf( pointplt, "set border\n" );
            fprintf( pointplt, "%s %c\n", "plot", 92 );
            fprintf( pointplt, "'%s'\n", datfle );

            for ( kk = null; kk<ind; kk++ )
            {
               fprintf( pointdat, "%+.15e, %+.15e\n",
	          pnt[kk][null], pnt[kk][ONE] );
            };

            nseconds = time( timer );
            strncpy( ctmptr, ( ctime( &nseconds )), 24 );

            TIMEFORM( tmestr, ctmptr );

            fprintf( pointplt, "\n#SYSSMX point plot file %s", pntplt );
            fprintf( pointplt, "\n#created: %24s\n", tmestr );

	    fclose( pointplt );
	    fclose( pointdat );
	 };
# endif /* SMX_POINTPLT != 0 */
/*............................................................................*/
         ind = clblcs( null, MOD_BLOCKS );           /* clear blc.border[]... */
                                                    /*  [ c.f. blocks.h ]     */
         ind = blocks( jj , "coordinates" );       /* generate blocks         */
/*...............................................*/
         if ( ind == ONE )             /* Error */
         {
            printf( "\n\n Message from function %s :", DO_SYSSMX );
            printf( "\n Error on block definition function"
               " blocks(*) [ option 'c'oordinates ] !!!\n " );
            exit( EXIT_FAILURE );
         }
         else if ( ind == -ONE ) /* Escape */
         {
            ( ppt->rtn ) = -ONE;
            return ppt;
         };
/*............................................................................*/
# if SMX_REPEATER == 1  /* repeat s-matrix labels   */

         if (( blc.base[ii] < jj )
           &&( jj < blc.base[ii+ONE] )
           &&( SMX_BOTLIMIT < ii )
           &&( ii < SMX_TOPLIMIT ))
         {
            nn = lbl.mi;
            mm = lbl.mf;
            while( nn <= mm )
            {
               hh = nn - blc.dm[ii]; /* periodic cell !!! */
               ( ppt->rep->smx[nn] ) = hh;
/*............................................................................*/
# if DSC_HCRMDE != 0
               ( ppt->rhp->smx[nn] ) = hh;
# endif
/*............................................................................*/
               nn++ ;                        /*[ function formdrv(*) will     */
            };                               /*  recursively set rep.smx[nn]  */
         };                                  /*  to the s-matrix label of the */
                                             /*  FIRST cell in the here just  */
                                             /*  generated backward pointing  */
# endif /* SMX_REPEATER == 1 */              /*  chain of periodic cells.]    */
/*............................................................................*/
/*
         printf( "\n please acknowlege...: " );
         scanf( "%s", ptr );
*/
/*............................................................................*/
         if (( jj == blc.base[ii] )&&( blc.cov == null ))
         {
            blc.dm[ii] = lbl.m  - lbl.mi + ONE;
            blc.dp[ii] = lbl.pf - lbl.pi + ONE; /* lbl.pf may actually exceed */
         };                                     /* lbl.p  on a domain change, */
                                                /* ii->ii+1: jj = blc.base[ii]*/
     
         if (( jj < blc.base[ii+ONE] )||( jj == blc.base[blc.m[null]] ))
         {
/*............................................................................*/
# if SMX_SYSPLT == 1
            kk = null;
            while ( kk < ( ssp->lay[null] ))
            {
               kk++;
               if ( jj == ( ssp->lay[kk] ))
               {
                  ( ssp->mi ) = lbl.mi;
                  ( ssp->mf ) = lbl.mf;
                  ( ssp->pm ) = blc.dm[ii];
                  ( ssp->layer ) = jj;

/*............................................................................*/
                  ssp = sysplt( ssp );               /* 'gnuplot' of layer jj */
/*.................................................*/

               };
            };
# endif /* SMX_SYSPLT == 1 */
/*............................................................................*/
            blc.mm += blc.dm[ii];
            blc.pp = lbl.pi + blc.dp[ii] - ONE;

         };/* end if ( jj < blc.base[ii+ONE] ) */
         jj++;
      }; /* end while ( jj <= blc.base[ii+ONE] ) */
      ii++;
   };/* end while ( ii < blc.m[null] ) */
/*............................................................................*/
/* stop time: */

   nseconds = time( timer );
   strncpy( ctmptr, ( ctime( &nseconds )), 24 );

   TIMEFORM( tmestr, ctmptr );
/*............................................................................*/
# if SMX_POINTLOG ==1
/*............................................................................*/
# if SMX_SYSPLT == 1
   if ( null < ( ssp->lay[null] ))
      printf( "\n\n" );
# endif
/*............................................................................*/
   fprintf( pointlog, "\nSYSSMX logfile %s ", logptr );
   fprintf( pointlog, "created:\n%24s\n", tmestr );

   fclose ( pointlog );

   printf( CLEAR_LINE ); 
   printf( "\r Supporting point logfile %s ", logptr );
   printf( timefrm, tmestr );

# else /* if SMX_POINTLOG != 1 */

   printf( CLEAR_LINE ); 
   printf( "\r Coordinates written with function %s ", SYSSMX_DO );
   printf( timefrm, tmestr );

# endif /* SMX_POINTLOG ... */
/*............................................................................*/
# if SMX_SYSPLT == 1

   if ( null == ( ssp->lay[null] ))
   {
      printf( "\n\n Create plot file [ global mesh ]:\n" );
      ( ssp->mi ) = ( tpt->mi );
      ( ssp->mf ) = ( tpt->mf );
      ( ssp->pm ) = blc.dm[0];
      ( ssp->layer ) = -ONE;
      ( ssp->toplr ) = blc.base[blc.m[null]]; /* the top layer */
/*............................................................................*/
      ssp = sysplt( ssp );      /* 'gnuplot' of entire DSC mesh system        */
/*............................*/
   };
# endif /* SMX_SYSPLT == 1 */
/*............................................................................*/
# if SMX_SANDWCHK == 1
/* perform sandwich check: */
   
   mm = null;
   ii = null;
   while ( ii < blc.m[null] )
   {
      jj = blc.base[ii];
      while( jj < blc.base[ii+ONE] )
      {
         hh = ONE;
         while ( hh <= blc.dm[ii] )
         {
            nn = mm + hh;

            kk = null; 
            while( kk < FOUR )
            {
               pp = tpt->cm[nn][kk];
               qq = tpt->cm[nn][kk+FOUR];
               rr = pp + blc.dp[ii]; /* periodic point */

               if ( qq != rr )
               {
                  printf( "\n\n Message from function %s :", DO_SYSSMX ); 
                  printf( "\n\n DSC mesh sandwich check failed !!!");   
                  printf( "\n [ Corner point periodicity condition hurt by "
                     "points no.%ld and %ld != %ld,\n   on cell no.%ld. ]",
                     pp, qq, rr, nn );

                  ppt->rtn = ONE;
                  return ppt;
               };
/*............................................................................*/
# if SMX_COORDCHK == 1
/* periodicity of coordinates: */

               ll = null;
               while( ll < TWO )
               {
                  uu = ppt->cpt->c[pp][ll];
                  vv = ppt->cpt->c[qq][ll];
                  if ( bound1 < fabs( uu - vv ) )
                  {
                     printf( "\n\n Message from function '%s' : ", DO_SYSSMX );
                     printf( "\n\n DSC mesh sandwich check failed !!! " );
                     printf( "\n [ Corner point periodicity condition hurt by "
                             "points no.%ld and %ld,\n   on cell no.%ld. ]",
                              pp, qq, nn );

                     ppt->rtn = ONE;
                     return ppt;
                  }; 
                  ll++;
               };/* while ll < TWO */

# endif /* SMX_COORDCHK == 1 */
/*............................................................................*/
               kk++;
            };/* while kk < FOUR  */
            hh++;
         };
         mm += blc.dm[ii];
         jj++;
      };/* while jj < blc.base[ii+ONE]  */
      ii++;
   };/* while ii < blc.m[null] */
# endif /* SMX_SANDWCHK == 1 */
/*............................................................................*/

   ( ppt->rtn ) = null;
   return ppt;
}
/*==================== end of function body syssmx(*) ========================*/
# undef SMX_SANDWCHK
# undef SMX_COORDCHK
# undef SMX_REPEATER 
# undef SMX_BOTLIMIT 
# undef SMX_TOPLIMIT 
# undef SMX_POINTLOG
# undef SMX_POINTPLT
# undef SMX_SYSPLT
/*********************** end of function syssmx(*) ****************************/


/* [ function: sysbnd ] */
/*******************************************************************************
*                                                                              *
*   ANSI-C function sysbnd(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   This subroutine configurates DSC system boundary file 'bndry...'.          *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "./tools/bndrfl.h"
# include "./tools/cshptp.h"
/*----------------------------------------------------------------------------*/
# ifndef MOD_ENBLRAD
   # define MOD_ENBLRAD 1 /* 1: enable heat radiative boundary conditions */
# endif
# ifndef MOD_ENBLESC
   # define MOD_ENBLESC 1 /* 1: enable envrnmt-surface heat conductivity */
# endif
# ifndef MOD_ENBLSKN
   # define MOD_ENBLSKN 1 /* 1: enable [ skin effect ] heat source bnd faces */
# endif
/*----------------------------------------------------------------------------*/
# if MOD_ENBLRAD == 1
/* The Stefan - Boltzmann constant [ W/(K^4*m^2) ]: */
# ifndef STEFAN_BOLTZ
   # define STEFAN_BOLTZ ( 5.6705100000e-08 )
# endif
/* Temperature shift [ degree Celsius to Kelvin ]: */
# ifndef CELSIUS_TO_KELVIN
   # define CELSIUS_TO_KELVIN ( 273.15 )
# endif
# endif
/*============================================================================*/

BOUNDARIES *
sysbnd( FORMSTATE *state )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static TOPOLOGY *tpt;
   static PARAMETERS *ppt;
   static BOUNDARIES *bpt;

   static double
      sigma = ZERO; /* equival. surface conductivity 1/ZL [Siemens]  */

   static short
      pp = null,
      prt = null,
      ports = null;

   static char
      face = null;

   static long
      ii = null,
      mm = null,
      nn = null;

   static short
      jj = null,
      kk = null;

   static char 
      tpp[SHS_SIZE] = {null},
     *sysptr = DSC_MODEL;

   WVGDPAR 
     *wve = NULL;

   CSHAPE
     *crd = NULL;
/*
   static char 
      ptr[SHS_SIZE] = {null},
    **endp = null;
*/
/*............................................................................*/
# if DSC_HCRMDE != 0
   static double
      ss = ZERO,
      tt = ZERO,
      hcv = ZERO;
/*............................................................................*/
# if MOD_ENBLSKN == 1
   static short
      ll = null;

   static double
      admce = ZERO;
# endif
/*............................................................................*/
# if MOD_ENBLRAD == 1
   static double
      evt = ZERO,
      t4c = ZERO;
# elif MOD_ENBLESC == 1
   static double
      evt = ZERO;
# endif
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* prototypes: */

   WVGDPAR
     *wvepar( char *type, double a,  double b,
            double eps, double my, double f );

   CSHAPE
     *cshape( CSHAPE *crd );

# if DSC_HCRMDE != 0
   double
      pow( double x, double y ); /* function x^y */
# endif
/*----------------------------------------------------------------------------*/
/* copy FORMSTATE [ state from formdrv: mesh topology information, e.g.]: */

   tpt = ( state->tpt );
   ppt = ( state->ppt );
   bpt = ( state->bpt );
   spt = state;
/*............................................................................*/
/* initialize struct CSHAPE *crd: */

   crd = cshape( NULL );

  ( crd->tpt ) = tpt;
  ( crd->ppt ) = ppt;
/*............................................................................*/
# if DSC_HCRMDE != 0
/* the characteristic field admittance */
   admce = sqrt( EPS_VAC / MY_VAC_ );
# endif
/*............................................................................*/
/* copy identifiers_and_comments: */

   strcpy(( bpt->name ), sysptr );
   strcpy(( bpt->text ), ( tpt->text ));
/*...........................................................................*/
   printf( "\n ======================================="
              "=======================================" );
   PRBLDCLR( "\n" );
   printf( "\r %*s", 78, "SYSBND" );
   PRNORMAL( "" );

   if( MOD_BOUNDRS < trf.bnd )
   {
      fprintf( stderr, "\n\n Error message from function %s :", __func__ );
      fprintf( stderr, "\n Too many boundary faces defined !!!" );
      fprintf( stderr, "\n [ Number %ld exceeds maximum number %ld",
         ( long ) trf.bnd, ( long ) MOD_BOUNDRS );
      fprintf( stderr, "\n   = macro MOD_BOUNDRS in transfer structure trf" );
      fprintf( stderr, "\n   - Change macro in compliance with memory "
         "resources.]\n" );
      exit( EXIT_FAILURE );
   };
/*............................................................................*/
   ( bpt->n ) = null;

   ports = trf.bnd; /* number of electric boundary ports [ or port sections ] */

   prt = null;
   while ( prt < ports ) 
   {
/*............................................................................*/
/* boundary type; conversion to lowercase: */

      strncpy( tpp, trf.btp[prt], SHS_SIZE );
      kk = strlen( tpp );

      if ( SHS_SIZE < kk )
         kk = SHS_SIZE;

      jj = null; 
      do
      {
         tpp[jj] = tolower( tpp[jj] );
      } while (( jj++ ) < kk );

      strncpy( trf.btp[prt], tpp, SHS_SIZE );
/*............................................................................*/
      if ( null == strncmp( tpp, "coaxial", THREE ))
      {
         wve = wvepar( "coax", 2*trf.ra, 2*trf.ri, trf.bdd[prt], 1., trf.fr );
         sigma = 1./( wve->z1 );
      }
      else if ( null == strncmp( tpp, "rectangular", THREE ))
      {
         wve = wvepar( "rectangular", trf.a, trf.b, trf.bdd[prt], 1., trf.fr );
         sigma = 1./( wve->zh );
      }
      else if ( null == strncmp( tpp, "circular", THREE ))
      {
         wve = wvepar( "circular", trf.a, trf.b, trf.bdd[prt], 1., trf.fr );
         sigma = 1./( wve->zh );
      }
      else if ( null == strncmp( tpp, "elliptic", THREE ))
      {
         wve = wvepar( "elliptic", trf.a, trf.b, trf.bdd[prt], 1., trf.fr );
         sigma = 1./( wve->zh );
      }
      else if (( null == strncmp( tpp, "skin_effect", TWO ))
             ||( null == strcmp( tpp, "sr" )))
      {
         /* resistive [ e.g skin effect ] losses: */

         if ( fabs( trf.bdd[prt] ) < 1.e-77 )
            sigma = 1.000e+77;
         else
            sigma = ( 1./trf.bdd[prt] );
      };
/*............................................................................*/
      nn = trf.bnn[prt]; /* the number of periods */

      mm = ( bpt->n ) + nn;
      if ( BNDAP <= mm )
      {
         fprintf( stderr, "\n\n ERROR message from function %s :", __func__ );
         fprintf( stderr, "\n Too many aperiodic boundary faces defined !!!"  );
         fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
            ( long ) BNDAP  );
         fprintf( stderr, "\n   = macro BNDAP in file FORMER.CONF" );
         fprintf( stderr, "\n   - Change macro in compliance with "
            "memory ressources.]\n\n" );

         exit( EXIT_FAILURE );
      };

      ii = trf.bii[prt]; /* the initial cell */
      pp = trf.bpp[prt]; /* the period */
      face = trf.bff[prt];

      mm = null;
      while ( mm < nn )
      {
/*............................................................................*/
         bndrfl( ii, face, sigma );                  /*                       */
/*.................................................*/
         ii += pp; /* add the period */
         mm++;
      }; /* while ( mm < trf.bnn[prt] ) */
      prt++ ;
   }; /* while ( prt < ports ) */
/*............................................................................*/
# if DSC_HCRMDE != 0
   ss = ZERO;
   tt = ZERO;

   ( bpt->ntf ) = null; /* temperature boundary face counter */
   ( bpt->ntn ) = null; /* temperature node counter */
   ( bpt->nhc ) = null; /* heat current boundary face counter */
   ( bpt->nsc ) = null; /* envrmnt-surface heat cond face counter */
   ( bpt->nsk ) = null; /* heat source boundary face counter */
/*............................................................................*/
# if DSC_FLDMDE != 0

   ( bpt->nns ) = null; /* no-slip boundary face counter */
   ( bpt->nsl ) = null; /* free slip boundary face counter */
   ( bpt->nif ) = null; /* inflow boundary face counter */
   ( bpt->nof ) = null; /* outflow boundary face counter */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
   ports = trf.bhc; /* the number of boundary temperature/ heat current ports */

   prt = null;
   while ( prt < ports )
   {
/*............................................................................*/
/* boundary type; conversion to lowercase: */

      strncpy( tpp, trf.bct[prt], SHS_SIZE );
      kk = strlen( tpp );

      if ( SHS_SIZE < kk )
         kk = SHS_SIZE;

      jj = null;
      do
      {
         tpp[jj] = tolower( tpp[jj] );
      } while (( jj++ ) < kk );

      strncpy( trf.bct[prt], tpp, SHS_SIZE );
/*............................................................................*/
/* fixed face temperature: */
	   
      if ( null == strcmp( tpp, "tf" ))
      {
         nn = trf.bcn[prt]; /* the number of [subsequent] bound temp cells */

         mm = ( bpt->ntf ) + nn;
         if ( BNDTF <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
	       __func__ );
            fprintf( stderr, "\n Too many fixed temperature faces "
	       "defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDTF  );
            fprintf( stderr, "\n   = macro BNDTF in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial boundary temperature cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */
         hcv = trf.bcr[prt]; /* the imposed face temperature [degree Celsius] */

         mm = null;
         while ( mm < nn )
         {
            jj = ( bpt->ntf );
            ( bpt->mtf[jj] ) = ii;
            ( bpt->ftf[jj] ) = face;
            ( bpt->tf[jj] ) = hcv;
            ( bpt->ntf )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end: fixed face temperatures */
/*............................................................................*/
/* fixed incident heat current densities: */

      if (( null == strncmp( tpp, "heat_current", THREE  ))
        ||( null == strcmp( tpp, "hc" )))
      {
         nn = trf.bcn[prt]; /* the number of [subsq] heat curr boundary cells */

         mm = ( bpt->nhc ) + nn;
         if ( BNDHC <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
               __func__ );
            fprintf( stderr, "\n Too many incident heat current faces "
               "defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDHC  );
            fprintf( stderr, "\n   = macro BNDHC in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial heat current boundary cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */

         if ( 1.e-277 < trf.abs0 )
            hcv = trf.abs0*trf.bcr[prt]; /* + light absorbed [ W/m^2 ] */
         else
            hcv = trf.bcr[prt];          /* + light absorbed [ W/m^2 ] */

         mm = null;
         while ( mm < nn )
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            jj = ( bpt->nhc );
            ( bpt->mhc[jj] ) = ii;
            ( bpt->fhc[jj] ) = face;
            ( bpt->hc[jj] ) = hcv*( crd->fm[(short)face] );
            ( bpt->rd[jj] ) = ZERO;
            ( bpt->nhc )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end incident heat current boundary conds */ 
/*............................................................................*/
# if MOD_ENBLRAD == 1
/* Stefan-Boltzman heat radiation [ exchanged at boundary faces ]: */

      if (( null == strncmp( tpp, "stefan-boltzmann", TWO )) 
        ||( null == strncmp( tpp, "radiation", THREE ))
        ||( null == strcmp( tpp, "rd" ))
        ||( null == strcmp( tpp, "sb" )))
      {
         nn = trf.bcn[prt]; /* the number of [subsq] heat rad boundary cells */

         mm = ( bpt->nhc ) + nn;
         if ( BNDHC <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
               __func__ );
            fprintf( stderr, "\n Too many Stefan-Boltzmann radiation "
               "- incident"  );
            fprintf( stderr, "\n heat current - boundary faces defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDHC  );
            fprintf( stderr, "\n   = macro BNDHC in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial heat rad boundary cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */
         evt = trf.bcr[prt]; /* environment temperature */
         t4c = STEFAN_BOLTZ*trf.rad0*trf.rad1;
         hcv = ( t4c * pow(( evt + CELSIUS_TO_KELVIN ), 4. ));
                      /* incident [ environmental ] heat radiation density */
         mm = null;
         while ( mm < nn )
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            jj = ( bpt->nhc );
            ( bpt->mhc[jj] ) = ii;
            ( bpt->fhc[jj] ) = face;
            ( bpt->hc[jj] ) = hcv*( crd->fm[(short)face] );
            ( bpt->rd[jj] ) = - t4c*( crd->fm[(short)face] );
            ( bpt->nhc )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end heat radiation boundary conditions */ 
# endif /* MOD_ENBLRAD == 1 */
/*............................................................................*/
# if MOD_ENBLESC == 1
/* surface-environment conductivity [approximate convection, e.g.]:*/

      if ( null == strcmp( tpp, "sf" ))
      {
         nn = trf.bcn[prt]; /* the number of [subsq] heat rad boundary cells */

         mm = ( bpt->nsc ) + nn;
         if ( BNDSC <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :", 
               __func__ );
            fprintf( stderr, "\n Too many environment thermal contact "
               "faces defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDSC  );
            fprintf( stderr, "\n   = macro BNDSC in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial heat rad boundary cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */
         evt = trf.bcr[prt]; /* environment temperature [Celsius] */
         hcv = trf.shc;      /* envrmt-surface heat conductivity [W/(K*m^2)] */

         mm = null;
         while ( mm < nn )
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            jj = ( bpt->nsc );
            ( bpt->msc[jj] ) = ii;
            ( bpt->fsc[jj] ) = face;
/* surface heat conductance [W/K]: */
            ( bpt->sc[jj] ) = hcv*( crd->fm[(short)face] );
/* reference temperature [Celsius] */
            ( bpt->tr[jj] ) = evt;
            ( bpt->nsc )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end: surface heat conductivities */ 
# endif /* MOD_ENBLESC == 1 */
/*............................................................................*/
# if MOD_ENBLSKN == 1
/* skin effect lossy heat sources: */

      if (( null == strncmp( tpp, "sources", THREE  ))
        ||( null == strcmp( tpp, "sc" )))
      {
         nn = trf.bcn[prt]; /* the number of [subsq] heat source cells */

         mm = ( bpt->nsk ) + nn;
         if ( BNDSK <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :", 
               __func__ );
            fprintf( stderr, "\n Too many skin effect lossy faces "
               "defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDSK  );
            fprintf( stderr, "\n   = macro BNDSK in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial heat source cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */
         hcv = trf.bcr[prt]; /* the surface resistance [ R_square ] */

         mm = null;
         while ( mm < nn )
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            jj = ( bpt->nsk );
            ( bpt->msk[jj] ) = ii;
            ( bpt->fsk[jj] ) = face;
            ( bpt->rs[jj] ) = hcv;

/* ss = ( 1/Z_field )*sqrt( R_square * the size of the face ): */

            ss = admce*sqrt( hcv*( crd->fm[(short)face] ));

/* [essentially] compute ss*( crd->vu ) with ( crd->vu ) := adj(P^-1), */
/* the adjoint port vector matrix [ in a local ON coordinate sytem ]:  */

            kk = null; do
            {
               ll = null; do
               {
                  ( bpt->sk[jj][kk][ll] ) = ss * \
                     ( ONE - TWO*(( face+ll )%TWO )) * \
                        ( crd->vu[kk][(( ll+ONE )%TWO )] );
/*
same as:
                  ( bpt->sk[jj][kk][(( ll+ONE )%TWO )] ) = ss * \
                     ( TWO*(( face+ll )%TWO ) - ONE ) * \
                        ( crd->vu[kk][ll] );
*/
               } while(( ++ll )< TWO );
            } while(( ++kk )< TWO );

            ( bpt->nsk )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end: electric/magnetic heat sources [on faces] */
# endif /* MOD_ENBLSKN == 1 */
/*............................................................................*/
# if DSC_FLDMDE != 0
/* no-slip boundary conditions: */

      if (( null == strncmp( tpp, "no_slip", TWO ))
        ||( null == strcmp( tpp, "ns" )))
      {
         nn = trf.bcn[prt]; /* the number of [subsequent] bound temp cells */

         mm = ( bpt->nns ) + nn;
         if ( BNDNS <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
               __func__ );
            fprintf( stderr, "\n Too many no-slip boundary faces "
               "defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDNS  );
            fprintf( stderr, "\n   = macro BNDNS in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial boundary temperature cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */
	 ss = trf.bcr[prt]; /* Nusselt number */

         mm = null;
         while ( mm < nn )
         {
            jj = ( bpt->nns );
            ( bpt->mns[jj] ) = ii;
            ( bpt->fns[jj] ) = face;
/*............................................................................*/
# if NUSSELT != 0
            ( bpt->nus[jj] ) = ss;
# endif
/*............................................................................*/
            ( bpt->nns )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end, no-slip boundary conds */
/*............................................................................*/
/* free slip boundary conditions: */

       if (( null == strncmp( tpp, "free_slip", THREE ))
         ||( null == strncmp( tpp, "sl", TWO  )))
      {
         nn = trf.bcn[prt]; /* the number of [subsequent] bound temp cells */

         mm = ( bpt->nsl ) + nn;
         if ( BNDSL <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :", 
               __func__ );
            fprintf( stderr, "\n Too many free slip boundary faces "
               "defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDSL  );
            fprintf( stderr, "\n   = macro BNDSL in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial boundary temperature cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */

         mm = null;
         while ( mm < nn )
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            jj = ( bpt->nsl );
            ( bpt->msl[jj] ) = ii;
            ( bpt->fsl[jj] ) = face;

            kk = null; do
            {
               if ( 1.e-277 < ( crd->fm[( short )face] ))
               {
                  ( bpt->nf[jj][kk] ) = ( crd->nf[( short )face][kk] );
               }
               else
                  ( bpt->nf[jj][kk] ) = ZERO;
            } while(( ++kk ) < THREE );

            ( bpt->nsl )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end, tangential slip boundary conds */
/*............................................................................*/
/* inflow boundary conditions: */

      if (( null == strncmp( tpp, "inflow", TWO ))
        ||( null == strncmp( tpp, "if", TWO ))
        ||( null == strncmp( tpp, "uf", TWO )))
      {
         nn = trf.bcn[prt]; /* the number of [subsequent] bound temp cells */

         mm = ( bpt->nif ) + nn;
         if ( BNDIF <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
               __func__ );
            fprintf( stderr, "\n Too many inflow boundary faces "
               "defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDIF  );
            fprintf( stderr, "\n   = macro BNDIF in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial boundary temperature cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */
	 ss = trf.bcr[prt]; /* inflow velocity [m/s] */
	 tt = trf.bcs[prt]; /* inflow temperature */

         mm = null;
         while ( mm < nn )
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            jj = ( bpt->nif );

            ( bpt->mif[jj] ) = ii;
            ( bpt->fif[jj] ) = face;

            ( bpt->ti[jj] ) = tt; /* inflow temperature */

            kk = null; do
            {
               if ( 1.e-277 < ( crd->fm[(short)face] ))
                  ( bpt->uf[jj][kk] ) = - ( ss*( crd->nf[(short)face][kk] ));
               else
                  ( bpt->uf[jj][kk] ) = ZERO;
            } while(( ++kk ) < THREE );
            ( bpt->nif )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end, inflow boundary conds */
/*............................................................................*/
/* outflow boundary conditions: */

      if (( null == strncmp( tpp, "outflow", TWO ))
        ||( null == strncmp( tpp, "of", TWO )))
      {
         nn = trf.bcn[prt]; /* the number of [subsequent] bound temp cells */

         mm = ( bpt->nof ) + nn;
         if ( BNDOF <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
               __func__ );
            fprintf( stderr, "\n Too many outflow boundary faces "
               "defined !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld", 
               ( long ) BNDOF  );
            fprintf( stderr, "\n   = macro BNDOF in file FORMER.CONF" );
            fprintf( stderr, "\n   - Change macro in compliance with "
               "memory ressources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         ii = trf.bci[prt]; /* the initial boundary temperature cell [index] */
         pp = trf.bcp[prt]; /* the period */
         face = trf.bcf[prt]; /* the [common] face index */
         ss = trf.bcr[prt]; /* outflow velocity [m/s] */

         mm = null;
         while ( mm < nn )
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            jj = ( bpt->nof );
            ( bpt->mof[jj] ) = ii;
            ( bpt->fof[jj] ) = face;

            kk = null; do
            {
               if ( 1.e-277 < ( crd->fm[(short)face] ))
               {
                  ( bpt->vf[jj][kk] ) = ss*( crd->nf[(short)face][kk] );
               }
               else
                  ( bpt->vf[jj][kk] ) = ZERO;
            } while(( ++kk ) < THREE );
            ( bpt->nof )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.bcn[prt] ) */
      }; /* end, outflow boundary conds */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      prt++ ;
   }; /* while ( prt < ports ) */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/

   PRBLDCLR( "\n" );
   printf( "\r %*s", 78, "SYSBND" );
   PRNORMAL( "" );
   ( bpt->rtn ) = null;
   return bpt;
}
/*==================== end of function body sysbnd(*) ========================*/


/* [ function: sysexc ] */
/*******************************************************************************
*                                                                              *
*   ANSI-C function sysexc(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0.]                                                            *
*                                                                              *
*   DSC model dependent excitation parameter function:                         *
*   This subroutine defines DSC field excitation parameters,                   *
*   which are written into static struct excitation 'exc' .                    *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/* excitation type: EXC_MAXW = "DIRAC_PULSE", "HARMONIC", "MULTIPLE_HARMONIC" */
/*
# define EXC_MAXW  "DIRAC_PULSE"
# define EXC_MAXW  "GAUSS_PULSE"  
# define EXC_MAXW  "HARMONIC"
# define EXC_HCRR  "PASSIVE"
# define EXC_HCRR  "PERIODIC_PULSE"

# include "./tools/cshptp.h"
*/
# include "./tools/excfld.h"
/*============================================================================*/

EXCITATION *
sysexc( FORMSTATE *state )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static TOPOLOGY *tpt;
   static PARAMETERS *ppt;
   static EXCITATION *ept;

   static short
      pp = null,
      prt = null,
      ports = null;

   static long 
      ii = null,
      mm = null,
      nn = null;

   static short
      jj = null,
      kk = null;

   static char 
      face = null,
      tpp[SHS_SIZE] = {null},
     *sysptr = DSC_MODEL,
     *excptr = EXC_MAXW; 

   static CSHAPE
      *crd = NULL; 
/*
   static char 
     *ptr,
    **endp =  null;
*/
   static double
      xx = ZERO,
      yy = ZERO,
      ux = ZERO,
      uy = ZERO,
      uz = ZERO,
      uu = ZERO,
      rr = ZERO,
      phase = ZERO;

   TOPOLOGY *
      systop( FORMSTATE *state );

   PARAMETERS *
      syssmx( FORMSTATE *state );

   int excfld( long c, char f, char ftype, 
               double x, double y, double z, double phase );

   WVGDPAR *
      wvepar( char *type, double a,  double b,
              double eps, double my, double f ); 

   CSHAPE *
      cshape( CSHAPE* crd );
/*----------------------------------------------------------------------------*/
/* copy FORMSTATE [ state from formdrv: mesh topology information, e.g.]: */

   tpt = ( state->tpt );
   ppt = ( state->ppt );
   ept = ( state->ept );
   spt = state;
/*............................................................................*/
/* initialize struct CSHAPE *crd: */

   crd = cshape( NULL );

  ( crd->tpt ) = tpt;
  ( crd->ppt ) = ppt;
/*............................................................................*/
/* copy identifiers_and_comments: */

   strcpy(( ept->name ), sysptr );
   strcpy(( ept->text ), ( tpt->text));
   strcpy(( ept->type ), excptr );
/*............................................................................*/
/*                                           [ usually called in 'formdrv.h' ]
   syssmx( state );
*/ 
/*............................................................................*/

   printf( "\n ======================================="
              "=======================================" );
   PRBLDCLR( "\n" );
   printf( "\r %*s", 78, "SYSEXC" );
   PRNORMAL( "" );

   if( MOD_EXCITES < trf.exc )
   {
      fprintf( stderr, "\n\n Error message from function %s :", __func__ );
      fprintf( stderr, "\n Too many excitation faces defined !!!" );
      fprintf( stderr, "\n [ Number %ld exceeds maximum number %ld",
         ( long ) trf.exc, ( long ) MOD_EXCITES );
      fprintf( stderr, "\n   = macro MOD_EXCITES in transfer structure trf" );
      fprintf( stderr, "\n   - Change macro in compliance with memory "
         "resources.]\n" );

      exit( EXIT_FAILURE );
   };

   if( *( ppt->domain ) == 'f' )
   {
      strcpy(( ept->type ), "STEADY_STATE_______" );
      ( ept->dt ) = 1000;
      ( ept->nn ) = 2 ; /* smoothing order */
   } /* end if *( ppt->domain ) == 'f'requency_domain */
   else /* if *( ppt->domain ) == 't'ime_domain       */
   {
      strcpy(( ept->type ), excptr );

      if ( null != strncmp( excptr, "DIRAC", TWO ))
      {
         ( ept->fr[null] ) = trf.fr;
         ( ept->rt ) = 1./( 2.*PI*trf.fr );

         if ( null == strncmp( excptr, "WAVE_PACKET", TWO ))
         {
            ( ept->dt ) = ( ept->rt );
            ( ept->ht ) = 50./trf.fr; /* plateau time for exc.type */
                                      /* "WAVE_PACKET_" */
            ( ept->nn ) = 2;          /* smoothing order [ 2 <= n ] */
         }
         else
            ( ept->dt ) = 3.*( ept->rt );
      };
   }; /* end if *ppt->domain != 'f'requency_domain */
/*............................................................................*/

   ( ept->ne ) = null;
   ( ept->nh ) = null;

   ports = trf.exc;

   prt = null; 
   while ( prt < ports ) 
   {
      ii = trf.eii[prt];
      nn = trf.enn[prt];
      pp = trf.epp[prt];
      uu = trf.euu[prt];
      face = trf.eff[prt];
/*............................................................................*/
/* excitation type; conversion to lowercase: */

      strncpy( tpp, trf.etp[prt], SHS_SIZE );
      kk = strlen( tpp );

      if ( SHS_SIZE < kk )
         kk = SHS_SIZE;

      jj = null;
      do
      {
         tpp[jj] = tolower( tpp[jj] );
      } while (( jj++ ) < kk );

      strncpy( trf.etp[prt], tpp, SHS_SIZE );
/*............................................................................*/
      mm = null;
      while ( mm < nn )
      {
         if ( null == strcmp( tpp, "e"  ))
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );       /*                                     */
/*...................................*/  
            xx = ( crd->xf[(short)face] );
            yy = ( crd->yf[(short)face] );
         /* zz = ( crd->zf[(short)face] ); not used */

/* TEM mode in coaxial waveguide centered at (0,0) */

            rr = xx*xx + yy*yy;
            ux = uu*xx/rr;
            uy = uu*yy/rr;
            uz = ZERO;
            phase = ZERO;

/* TE10 mode in rectangular waveguide */
/*
            ux = ZERO;
            uy = uu*cos( PI*xx/trf.a );
            uz = ZERO;
            phase = ZERO;
*/
/*............................................................................*/
            excfld( ii, face, 'e', ux , uy, uz, phase );       /*             */
/*...........................................................*/
         };/* end, E-field excitations */
         ii += pp; /* jump one period */
         mm++ ;
      }; /* while ( mm < trf.enn[prt] )  */
      prt++;
   }; /* while ( prt < ports ) */
/*............................................................................*/
# if DSC_HCRMDE != 0

   strcpy (( ept->hctp ), EXC_HCRR );

   ( ept->hcht ) = trf.hcht;
   ( ept->hcdt ) = trf.hcdt;

   ( ept->hcht2 ) = trf.hcht2;
   ( ept->hcdt2 ) = trf.hcdt2;

   ( ept->nhc ) = null; /* heat current excitation face counter */
   ( ept->ntf ) = null; /* temperature face excitation counter */
   ( ept->ntn ) = null; /* temperature node excitation counter */

   ports = trf.ehc; /* the number of boundary temperature/ heat current ports */

   prt = null;
   while ( prt < ports )
   {
/*............................................................................*/
/* excitation type; conversion to lowercase: */

      strncpy( tpp, trf.ect[prt], SHS_SIZE );
      kk = strlen( tpp );

      if ( SHS_SIZE < kk )
         kk = SHS_SIZE;

      jj = null;
      do
      {
         tpp[jj] = tolower( tpp[jj] );
      } while (( jj++ ) < kk );

      strncpy( trf.ect[prt], tpp, SHS_SIZE );
/*............................................................................*/
/* imposed heat current density: */

      if (( null == strncmp( tpp, "heat_current", THREE ))
        ||( null == strcmp( tpp, "hc" )))
      {
         nn = trf.ecn[prt]; /* the number of [subsq] heat curr boundary cells */

         mm = ( ept->nhc ) + nn;
         if ( EXCHC <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
               __func__ );
            fprintf( stderr, "\n Too many heat currents excited !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld",
               ( long ) EXCHC  );
            fprintf( stderr, "\n   = macro EXCHC in file FORMER.CONF ]" );
            fprintf( stderr, "\n Change macro in compliance with "
               "memory ressources !\n\n " );

            exit( EXIT_FAILURE );
         };

         ii = trf.eci[prt]; /* the initial heat current boundary cell [index] */
         pp = trf.ecp[prt]; /* the period */
         face = trf.ecf[prt]; /* the [common] face index */
   
         uu = trf.ecr[prt]; /* imposed heat current density [ W/m^2 ] */

         mm = null;
         while ( mm < nn )
         {
            ( crd->cell ) = ii;
            ( crd->face ) = face;
            ( crd->opt ) = 'f';
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            jj = ( ept->nhc );
            ( ept->mhc[jj] ) = ii;
            ( ept->fhc[jj] ) = face;
            ( ept->hc[jj] ) = uu * ( crd->fm[(short)face] );
            ( ept->nhc )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.ecn[prt] ) */
      } /* end, heat current excitation */
/*............................................................................*/
/* temperature excited on boundary: */
      else if ( null == strcmp( tpp, "tf"  ))
      {
         nn = trf.ecn[prt]; /* the number of [subsequent] bound temp cells */

         mm = ( ept->ntf ) + nn;
         if ( EXCTF <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
               __func__ );
            fprintf( stderr, "\n Too many temperature faces excited !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld",
               ( long ) EXCTF  );
            fprintf( stderr, "\n   = macro EXCTF in file FORMER.CONF ]" );
            fprintf( stderr, "\n Change macro in compliance with "
               "memory ressources !\n\n " );

            exit( EXIT_FAILURE );
         };

         ii = trf.eci[prt]; /* the initial boundary temperature cell [index] */
         pp = trf.ecp[prt]; /* the period */
         face = trf.ecf[prt]; /* the [common] face index */
         uu = trf.ecr[prt]; /* the imposed face temperature [degree Celsius] */

         mm = null;
         while ( mm < nn )
         {
            jj = ( ept->ntf );
            ( ept->mtf[jj] ) = ii;
            ( ept->ftf[jj] ) = face;
            ( ept->tf[jj] ) = uu;
            ( ept->ntf )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.ecn[prt] ) */
      } /* end, face temperatue excitation */
/*............................................................................*/
/* temperature excited in node: */
      else if ( null == strcmp( tpp, "tn"  ))
      {
         nn = trf.ecn[prt]; /* the number of [subsequent] bound temp cells */

         mm = ( ept->ntn ) + nn;
         if ( EXCTN <= mm )
         {
            fprintf( stderr, "\n\n ERROR message from function %s :",
               __func__ );
            fprintf( stderr, "\n Too many temperature nodes excited !!!"  );
            fprintf( stderr, "\n [ Number exceeds maximum number %ld",
               ( long ) EXCTN  );
            fprintf( stderr, "\n   = macro EXCTN in file FORMER.CONF ]" );
            fprintf( stderr, "\n Change macro in compliance with "
               "memory ressources !\n\n " );

            exit( EXIT_FAILURE );
         };

         ii = trf.eci[prt]; /* the initial boundary temperature cell [index] */
         pp = trf.ecp[prt]; /* the period */
         uu = trf.ecr[prt]; /* the imposed node temperature [degree Celsius] */

         mm = null;
         while ( mm < nn )
         {
            jj = ( ept->ntn );
            ( ept->mtn[jj] ) = ii;
            ( ept->tn[jj] ) = uu;
            ( ept->ntn )++; /* count this cell */
            ii += pp; /* jump one period */
            mm++;
         }; /* while ( mm < trf.ecn[prt] ) */
      }; /* end, node temperatue excitation */
/*............................................................................*/
      prt++;
   }; /* while ( prt < ports ) */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/

   PRBLDCLR( "\n" );
   printf( "\r %*s", 78, "SYSEXC" );
   PRNORMAL( "" );
   ( ept->rtn ) = null;
   return ept;
}
/*==================== end of function body sysexc(*) ========================*/


/* [ function: sysval ] */
/*******************************************************************************
*                                                                              *
*   ANSI-C function sysval(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0.]                                                            *
*                                                                              *
*   DSC model dependent evaluation parameter function:                         *
*   This subroutine defines DSC field evaluation parameters,                   *
*   which are written into static struct evaluation 'val' .                    *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define EVL_TDITER 32768 /* number of iterations [ cycles ] - time domain */
# define EVL_TDINIT     1 /* first evaluated cycle - time domain */
# define EVL_TDSTOP 32768 /* last computed Maxwfld cycle - time domain */
# define EVL_TDRPTR     2 /* cycle length [ repetition rate ] - time domain */
# define EVL_FRITER 30000 /* number of interations [ cycles ] - freq domain */
# define EVL_FRINIT    21 /* first evaluated cycle - frequency domain */
# define EVL_FRSTOP   200 /* last computed Maxwfld cycle - frequency domain */
# define EVL_FRRPTR    20 /* cycle length [ repetition rate ] - freq domain */
# define EVL_HTINIT     1 /* first evaluated cycle - heat_&_fluids */
# define EVL_HTSTOP 30000 /* last computed cycle - heat_&_fluids */
# define EVL_HTRPTR     2 /* cycle length [ repetition rate ] - heat_&_fluids */
/*----------------------------------------------------------------------------*/
/*
# include "./tools/cshapt.h"
*/
/*============================================================================*/

EVALUATION *
sysval( FORMSTATE *state )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static TOPOLOGY *tpt;
   static PARAMETERS *ppt;
   static EVALUATION *vpt;

   static FILE *evlprts;

   static signed char
      sgn = null,
      prt = null,
      ports = null;

   static const signed char
      face[PORTS] = { 3, 5, 2, 4, 5, 1, 4, 0, 1, 3, 0, 2 },
      hport[PORTS] = { 10, -5, -12, 7, 2, -9, -4, 11, 6, -1, -8, 3 };

   static short 
      pp = null,
      qq = null,
      jj = null,
      kk = null;

   static long
      ii = null,
      mm = null,
      nn = null;

   static char 
      tpp[SHS_SIZE] = {null},
      ctmptr[STS_SIZE] = {null},
      tmestr[STS_SIZE] = {null};

   static char 
     *sysptr = DSC_MODEL,
     *mode   = "individual",
     *prtlog,
     *ptr;

   time_t 
      nseconds = null,
     *timer = null;
/*
   static char 
    **endp = null;
*/
   CSHAPE
     *crd = NULL;
/*............................................................................*/
# if DSC_HCRMDE != 0
   static signed char
      fce = null;
# endif
/*............................................................................*/
/* prototyping: */

   TOPOLOGY
     *systop( FORMSTATE *state );

   PARAMETERS
     *syssmx( FORMSTATE *state );

   CSHAPE
     *cshape( CSHAPE *crd );

   char
     *lotos( long n, char m );

   time_t time( time_t *timer);
   char *ctime( const time_t *timer );
/*----------------------------------------------------------------------------*/
/* memory allocations: */ 

   ptr = ( char * ) calloc( LGS_SIZE, ONE );
   prtlog = ( char * ) calloc( STS_SIZE, ONE );
/*............................................................................*/
/* copy FORMSTATE [ state from formdrv: mesh topology information, e.g.]: */

   tpt = ( state->tpt );
   ppt = ( state->ppt );
   vpt = ( state->vpt );
   spt = state;
/*............................................................................*/
/* initialize struct CSHAPE *crd: */

   crd = cshape( NULL );

  ( crd->tpt ) = tpt;
  ( crd->ppt ) = ppt;
/*............................................................................*/
/* copy identifiers_and_comments: */

   strcpy(( vpt->name ), sysptr );
   strcpy(( vpt->text ), ( tpt->text ));
/*............................................................................*/
/*                                           [ usually called in formdrv.h ]
   syssmx( state );
*/
/*............................................................................*/

   printf( "\n ======================================="
              "=======================================" );
   PRBLDCLR( "\n" );
   printf( "\r %*s", 78, "SYSVAL" );
   PRNORMAL( "" );

   strcpy(( vpt->mode_ep ), mode );
   strcpy(( vpt->mode_hp ), mode );
/*............................................................................*/
/* store evaluated ports: */

   strcpy( prtlog, "evl.log" );
   strcat( prtlog, ( state->flbl ));

   evlprts = fopen( prtlog, "w" );
  
   fprintf( evlprts, "Evaluated_ports_log" );
   fprintf( evlprts, "\n%s ", ( vpt->name ));
   fprintf( evlprts, "\n%s ", ( vpt->text ));
   fprintf( evlprts, "\n%s ", "[ Types, positions, directions, dimensions "
      "of evaluated quantities, etc. ]" );
/*............................................................................*/

   ( vpt->nep ) = null;
   ( vpt->nhp ) = null;

   ports = trf.val;

   prt = null;
   while ( prt < ports ) 
   {
      qq = trf.vpt[prt];
      ii = trf.vii[prt]; /* the initial cell index of a port [in this section]*/
      nn = trf.vnn[prt]; /* the number of ports [in this section] */
      pp = trf.vpp[prt]; /* the counting period of this section */
/*............................................................................*/
/* evaluation type; conversion to lowercase: */

      strncpy( tpp, trf.vtp[prt], SHS_SIZE );
      kk = strlen( tpp );

      if ( SHS_SIZE < kk )
         kk = SHS_SIZE;

      jj = null;
      do
      {
         tpp[jj] = tolower( tpp[jj] );
      } while (( jj++ ) < kk );

      strncpy( trf.vtp[prt], tpp, SHS_SIZE );
/*............................................................................*/
      strcpy( ptr, "[ " );
      strcat( ptr, trf.vtx[prt] );

      if ( null == strcmp( tpp, "e"  ))
         strcat( ptr, " - unit: Volts ]" );
      else if ( null == strcmp( tpp, "h"  ))
         strcat( ptr, " - unit: Amperes * sqrt( MY_VAC_/EPS_VAC ) ]" );

      mm = null;
      while ( mm < nn )
      {
         if ( null == strcmp( tpp, "e"  )) /* E-field port */
         {
            sgn = ONE;

            ( vpt->mep[( vpt->nep )] ) = ii;
            ( vpt->pep[( vpt->nep )] ) = qq;

            ( vpt->nep )++ ;

            if ( EVLEP <= ( vpt->nep ))
            {
               fprintf( stderr, "\n\n ERROR message from function %s :",
                  __func__ );
               fprintf( stderr, "\n Too many E-ports to be evaluated !!!"  );
               fprintf( stderr, "\n [ Number exceeds maximum number %ld",
                  ( long ) EVLEP  );
               fprintf( stderr, "\n   = macro EVLEP in file FORMER.CONF ]" );
               fprintf( stderr, "\n Change macro in compliance with "
                  "memory ressources !\n\n " );

               exit( EXIT_FAILURE );
            };

            ( crd->cell ) = ii;
/*............................................................................*/
            crd = cshape( crd );          /*                                  */
/*......................................*/
            printf( "\n E-port%d ( %ld | %d ) %s,",
               ( vpt->nep), ii, qq, ptr );
            printf( " length/m: %.7e", ( crd->pm[(qq-ONE)] ));

            fprintf( evlprts, "\n\nE-port%d: ( %ld | %d ) %s",
               ( vpt->nep ), ii, qq, ptr );
            fprintf( evlprts, "\nposition_/m: (% .7e, % .7e, % .7e )", \
               ( crd->xp[(qq-ONE)] ), ( crd->yp[(qq-ONE)] ), \
               ( crd->zp[(qq-ONE)] ));
            fprintf( evlprts, "\ndirection/m: (% .7e, % .7e, % .7e )", \
               ( crd->px[(qq-ONE)] )*sgn, ( crd->py[(qq-ONE)] )*sgn, \
               ( crd->pz[(qq-ONE)] )*sgn );
            fprintf( evlprts, "\nlength___/m:  % .7e", ( crd->pm[(qq-ONE)] ));
            fprintf( evlprts, "\nface%d__/m^2:  % .7e",
               face[(qq-ONE)], ( crd->fm[face[(qq-ONE)]] ));
         } /* end, port voltage evaluation */
         else if (( null == strcmp( trf.vtp[prt], "h"  )) /* H-field port */
                ||( null == strcmp( trf.vtp[prt], "H" ))) /* evaluation */
         {
            ( vpt->mhp[( vpt->nhp )] ) = ii;
            ( vpt->php[( vpt->nhp )] ) = qq;

            ( vpt->nhp )++ ;

            if ( EVLHP <= ( vpt->nhp ))
            {
               fprintf( stderr, "\n\n ERROR message from function %s :",
                  __func__ );
               fprintf( stderr, "\n Too many H-ports to be evaluated !!!"  );
               fprintf( stderr, "\n [ Number exceeds maximum number %ld",
                  ( long ) EVLHP  );
               fprintf( stderr, "\n   = macro EVLHP in file FORMER.CONF ]" );
               fprintf( stderr, "\n Change macro in compliance with "
                  "memory ressources !\n\n " );

               exit( EXIT_FAILURE );
            };

            kk = hport[( qq-ONE )];

            if ( null < kk )
               sgn = ONE;
            else
               sgn = -ONE;
            kk = abs(kk);

            ( crd->cell ) = ii;
/*............................................................................*/
            crd = cshape( crd );           /*                                 */
/*.......................................*/
            printf( "\n H-port%d ( %ld | %d ) %s,",
               ( vpt->nhp), ii, qq, ptr );
            printf( " length/m: %.7e", ( crd->pm[(kk-ONE)] ));

            fprintf( evlprts, "\n\nH-port%d: ( %ld | %d ) %s",
               ( vpt->nhp ), ii, qq, ptr );
            fprintf( evlprts, "\nposition_/m: (% .7e, % .7e, % .7e )", \
               ( crd->xp[(kk-ONE)] ), ( crd->yp[(kk-ONE)] ), \
               ( crd->zp[(kk-ONE)] )); 
            fprintf( evlprts, "\ndirection/m: (% .7e, % .7e, % .7e )",
               ( crd->px[(kk-ONE)] )*sgn, ( crd->py[(kk-ONE)] )*sgn,
               ( crd->pz[(kk-ONE)] )*sgn );
            fprintf( evlprts, "\nlength___/m:  % .7e", ( crd->pm[(kk-ONE)] ));
            fprintf( evlprts, "\nface%d__/m^2:  % .7e", face[(kk-ONE)],
               ( crd->fm[face[(kk-ONE)]] ));
         }; /* end, port current evaluation */
         ii += pp; /* add the period */
         mm++ ;
      }; /* while ( mm < trf.vnn[prt] ) */
      prt++;
   };/* while ( prt < ports ) */
/*............................................................................*/
# if DSC_HCRMDE != 0
/* evaluated thermal ports: */

   ( vpt->nhc ) = null;
   ( vpt->ntf ) = null;
   ( vpt->ntn ) = null;
/*............................................................................*/
# if DSC_FLDMDE != 0
   ( vpt->nun ) = null;
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
   ports = trf.vhc;

   prt = null; 
   while ( prt < ports ) 
   {
      qq = trf.vcf[prt];
      ii = trf.vci[prt];
      nn = trf.vcn[prt];
      pp = trf.vcp[prt];
      fce = trf.vcf[prt]; /* the [common] face index */
/*............................................................................*/
/* evaluation type; conversion to lowercase: */

      strncpy( tpp, trf.vct[prt], SHS_SIZE );
      kk = strlen( tpp );

      if ( SHS_SIZE < kk )
         kk = SHS_SIZE;

      jj = null;
      do
      {
         tpp[jj] = tolower( tpp[jj] );
      } while (( jj++ ) < kk );

      strncpy( trf.vct[prt], tpp, SHS_SIZE );
/*............................................................................*/
      strcpy( ptr, "[ " );
      strcat( ptr, trf.vcx[prt] );

      if (( null == strcmp( tpp, "tn"  ))
        ||( null == strcmp( tpp, "tf" )))
         strcat( ptr, " - unit: DEG Celsius ]" );
      else if ( null == strcmp( tpp, "un"  ))
         strcat( ptr, " - unit: meters/seconds ]" );

      mm = null;
      while ( mm < nn )
      {
         if ( null == strcmp( tpp, "tn"  ))
         {
/* node temperature evaluation: */

            ( vpt->mtn[( vpt->ntn )] ) = ii;
            ( vpt->ntn )++ ;

            if ( EVLTN <= ( vpt->ntn ))
            {
               fprintf( stderr, "\n\n ERROR message from function %s :",
                  __func__ );
               fprintf( stderr, "\n Too many temperature nodes to be "
                  "evaluated !!!"  );
               fprintf( stderr, "\n [ Number exceeds maximum number %ld",
                  ( long ) EVLTN  );
               fprintf( stderr, "\n   = macro EVLTN in file FORMER.CONF ]" );
               fprintf( stderr, "\n Change macro in compliance with "
                  "memory ressources !\n\n " );

               exit( EXIT_FAILURE );
            };

            ( crd->cell ) = ii;
/*............................................................................*/
            crd = cshape( crd );           /*                                 */
/*.......................................*/
            printf( "\n TN-node%d ( %ld | %c ) %s,",
               ( vpt->ntn ), ii, '*', ptr );

            fprintf( evlprts, "\n\nTN-node%d: ( %ld | %c ) %s",
               ( vpt->ntn ), ii, '-', ptr );
            fprintf( evlprts, "\nposition_/m: (% .7e, % .7e, % .7e )",
               ( crd->xn ), ( crd->yn ), ( crd->zn ));
            fprintf( evlprts, "\nvolume_/m^3:   %.7e", ( crd->vol ));
         } /* end, node temperature evaluation */
         else if ( null == strcmp( tpp, "tf"  ))
         {
/* face temperature evaluation: */
            ( vpt->mtf[( vpt->ntf )] ) = ii;
            ( vpt->ftf[( vpt->ntf )] ) = fce; ;
            ( vpt->ntf )++ ;

            if ( EVLTF <= ( vpt->ntf ))
            {
               fprintf( stderr, "\n\n ERROR message from function %s :",
                  __func__ );
               fprintf( stderr, "\n Too many temperature faces to be "
                  "evaluated !!!"  );
               fprintf( stderr, "\n [ Number exceeds maximum number %ld",
                  ( long ) EVLTF  );
               fprintf( stderr, "\n   = macro EVLTF in file FORMER.CONF ]" );
               fprintf( stderr, "\n Change macro in compliance with "
                  "memory ressources !\n\n " );

               exit( EXIT_FAILURE );
            };

            ( crd->cell ) = ii;
/*............................................................................*/
            crd = cshape( crd );           /*                                 */
/*.......................................*/
            printf( "\n TF-node%d ( %ld | %d ) %s,",
               ( vpt->ntf ), ii, fce, ptr );

            fprintf( evlprts, "\n\nTF-node%d: ( %ld | %d ) %s",
               ( vpt->ntf ), ii, fce, ptr );
            fprintf( evlprts, "\nposition_/m: (% .7e, % .7e, % .7e )",
               ( crd->xf[fce] ), ( crd->yf[fce] ), ( crd->zf[fce] ));
            fprintf( evlprts, "\ndirection/m: (% .7e, % .7e, % .7e )",
               ( crd->fx[fce] ), ( crd->fy[fce] ), ( crd->fz[fce] ));
            fprintf( evlprts, "\nface%d__/m^2:  % .7e", fce,
               ( crd->fm[fce] ));
         }; /* end, face temperature evaluation */
/*............................................................................*/
# if DSC_FLDMDE != 0
         if (( null == strncmp( tpp, "velocity", THREE ))
           ||( null == strcmp( tpp, "ve" ))
           ||( null == strcmp( tpp, "un" )))
         {
/* nodal velocity evaluation: */
            ( vpt->mun[( vpt->nun )] ) = ii;
            ( vpt->cun[( vpt->nun )] ) = ( signed char )( fce-120 ); /* ASCII */
            ( vpt->nun )++ ;             /* ( x, y, z ) to integer ( 0, 1, 2 )*/

            if ( EVLUN <= ( vpt->nun ))
            {
               fprintf( stderr, "\n\n ERROR message from function %s :",
                  __func__ );
               fprintf( stderr, "\n Too many nodal flows to be "
                  "evaluated !!!"  );
               fprintf( stderr, "\n [ Number exceeds maximum number %ld",
                  ( long ) EVLUN  );
               fprintf( stderr, "\n   = macro EVLUN in file FORMER.CONF ]" );
               fprintf( stderr, "\n Change macro in compliance with "
                  "memory ressources !\n\n " );

               exit( EXIT_FAILURE );
            };

            ( crd->cell ) = ii;
/*............................................................................*/
            crd = cshape( crd );           /*                                 */
/*.......................................*/
            printf( "\n UN-node%d ( %ld | %c ) %s,",
               ( vpt->nun ), ii, fce, ptr );

            fprintf( evlprts, "\n\nUN-node%d: ( %ld | %c ) %s",
               ( vpt->nun ), ii, fce, ptr );
            fprintf( evlprts, "\nposition_/m: (% .7e, % .7e, % .7e )",
               ( crd->xn ), ( crd->yn ), ( crd->zn ));
            fprintf( evlprts, "\nvolume_/m^3:   %.7e", ( crd->vol ));
/*
*//*        fprintf( evlprts, "\nU-component:   %d[%c]", ( fce-120 ), fce );
*//*        fprintf( evlprts, "\nU-component:   %c", fce );
*/
         }; /* end, nodal fluid velocity evaluation */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
         ii += pp; /* add the period */
         mm++ ;
      }; /* while ( mm < trf.vcp[prt] ) */
      prt++;
   }; /* while ( prt < ports ) */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/*
   printf( "\n\n Additional standing wave evaluation "
      "on source line: >>------------>" );
   printf( "\n Please enter number of ports [ Escape: "
      "enter null ] >>------------> \n" );
   printf( " >>---> ........................"
      "...................................: " );
   scanf( "%s", ptr );
   kk = strtol( ptr, endp, DEC );

   evl = THREE;
   ii = ONE;
   while (( null < evl )&&( evl < kk ))
   {
      printf( "\r >>---> enter %2d. section ...."
         ".....................................: ", ii );
      scanf( "%s", ptr );
      evl = strtol( ptr, endp, DEC );
      ( vpt->mep[( vpt->nep )] ) = trf.vli[null] + evl*blc.dm[null];
      ( vpt->pep[( vpt->nep )] ) = trf.vpt[null];
      ( vpt->nep )++ ;

      if ( EVLEP <= ( vpt->nep ))
      {
         fprintf( stderr, "\n\n ERROR message from function %s :",
            __func__ );
         fprintf( stderr, "\n Too many E-ports to be evaluated !!!"  );
         fprintf( stderr, "\n [ Number exceeds maximum number %ld",
            ( long ) EVLEP  );
         fprintf( stderr, "\n   = macro EVLEP in file FORMER.CONF ]" );
         fprintf( stderr, "\n Change macro in compliance with "
            "memory ressources !\n\n " );

         exit( EXIT_FAILURE );
      };

      ii++;
      evl += TWO;
   };
*/
/*............................................................................*/
/* stop time: */

   nseconds = time( timer );
   strncpy( ctmptr, ( ctime( &nseconds )), 24 );

   TIMEFORM( tmestr, ctmptr );
/*............................................................................*/
   fprintf( evlprts, "\n\nEvaluated ports logfile %s ", prtlog );
   fprintf( evlprts, "created:\n%24s\n", tmestr );

   fclose( evlprts );

   if( null == strncmp(( spt->ppt->domain ), "time_domain", ONE ))
   {
# if EVL_TDITER != 0
      ( vpt->n ) = EVL_TDITER;
      ( vpt->ni ) = EVL_TDINIT;
      ( vpt->r ) = EVL_TDRPTR;
# endif
      ;
# ifdef EVL_TDSTOP 
# if EVL_TDSTOP != 0
      ( vpt->nf ) = EVL_TDSTOP;
# else
      ( vpt->nf ) = ( vpt->n );
# endif
# else /* EVL_TDSTOP not defined */
      ( vpt->nf ) = ( vpt->n );
# endif /* EVL_TDSTOP not defined */
   }
   else if( null == strncmp(( spt->ppt->domain ), "frequency_domain", ONE ))
   {
# if EVL_FRITER != 0
      ( vpt->n ) = EVL_FRITER;
      ( vpt->ni ) = EVL_FRINIT;
      ( vpt->r ) = EVL_FRRPTR;
# endif
        ;
# ifdef EVL_FRSTOP
# if EVL_FRSTOP != 0
      ( vpt->nf ) = EVL_FRSTOP;
# else
      ( vpt->nf ) = ( vpt->n );
# endif
# else /* EVL_FRSTOP not defined */
      ( vpt->nf ) = ( vpt->n );
# endif /* EVL_FRSTOP not defined */
   };
/*............................................................................*/
# if DSC_HCRMDE != 0

   # if DSC_INTLCE == 0 /* separate internal loops [ cycles ] */
/*............................................................................*/
      # ifdef EVL_HTINIT /* [ the first evaluated cycle, heat and fluids ] */
         # if EVL_HTINIT != 0
/* the first evaluated cycle, heat and fluids */
            ( vpt->nj ) = EVL_HTINIT;
         # else
            ( vpt->nj ) = ONE;
         # endif
      # else /* EVL_HTINIT not defined */
         ( vpt->nj ) = ONE;
      # endif /* EVL_HTINIT not defined */
/*............................................................................*/
      # ifdef EVL_HTSTOP /* [ the last computed cycle, heat and fluids ] */
         # if EVL_HTSTOP != 0
            ( vpt->nt ) = EVL_HTSTOP;
         # else
            ( vpt->nt ) = ( vpt->n );
         # endif
      # else /* EVL_HTSTOP not defined */
         ( vpt->nt ) = ( vpt->n );
      # endif /* EVL_HTSTOP not defined */
/*............................................................................*/
      if ( trf.ref == null )
      {
      # ifdef EVL_HTRPTR /* [ the internal repetiton rate, heat and fluids ] */
         # if EVL_HTRPTR != 0
            ( vpt->rc ) = EVL_HTRPTR;
         # else
            ( vpt->rc ) = -ONE;
         # endif
      # else /* EVL_HTRPTR not defined */
         ( vpt->rc ) = ( vpt->r );
      # endif /* EVL_HTRPTR not defined */
   };
   # elif DSC_INTLCD == 1 /* interlaced internal loops */
   ( vpt->nj ) = ( vpt->ni );
   ( vpt->nt ) = ( vpt->nf );
   ( vpt->rc ) = ( vpt->r );
   # endif /* DSC_INTLCD == 1 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   PRBLDCLR( "\n" );
   printf( "\r %*s", 78, "SYSVAL" );
   PRNORMAL( "" );

   ( vpt->rtn ) = null;
   return vpt;
}
/*==================== end of function body sysval(*) ========================*/
# undef EVL_TDITER
# undef EVL_TDINIT
# undef EVL_TDSTOP
# undef EVL_TDRPTR
# undef EVL_FRITER
# undef EVL_FRINIT
# undef EVL_FRSTOP
# undef EVL_FRRPTR
# undef EVL_HTINIT
# undef EVL_HTSTOP
# undef EVL_HTRPTR
/*********************** end of function sysval(*) ****************************/
# endif /* [ end of section compiled with option -D_PRE_MODEL ] */
/*____________________________________________________________________________*/



/*____________________________________________________________________________*/
# ifdef _POST_MODEL
/* SECTION COMPILED IN OPTION -D_POST_MODEL */
/*******************************************************************************
*                                                                              *
*   ANSI-C function modval(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   release 1.0.                                                               *
*                                                                              *
*   This function determines special [ i.e. DSC model dependent ] evaluation   *
*   modes for files dsc.val<n>                                                 *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define _POSIX_SOURCE 1 /* some headers of the POSIX.1 standard will be used */
/*----------------------------------------------------------------------------*/
# define MDV_PLOT_SWR 1 /*1: create standing waves plots at frequencies frq[j]*/
# define MDV_SPCLPRTS 3 /* number of evalt'd special ports: 0<= n <MAX_PERIOD */
# define MDV_AUTOEVAL 0 /*1: automatic port labeling for standing wave eval.  */
# define MDV_SPLINTPL 1 /*>0: spline interpolate spc.format [references etc.] */
# define MDV_GPHINTPL 5 /*>0: spline interpl. gph.format [gnu.graphics,e.g.]  */
# define MDV_NORMALZE 0 /*1: normalize S-parmts to abs.value 1 [ cautious ! ] */
/*----------------------------------------------------------------------------*/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
# include <unistd.h> /* system specification header, cf. sysconf( ), etc. */
# include <time.h>   /* cf. time( ), ctime( ), asctime( ), localtime( ), etc. */
/*----------------------------------------------------------------------------*/
# include "./math/maths.h"  /* 'my' computation environment headers */
# include "./math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "./CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize this header for POSTER configuration: */
# include "./poster/POSTER.CONF" 
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
   # include <ncurses.h>
# endif 
/*----------------------------------------------------------------------------*/
# include "./tools/txctyp.h"  
# include "./poster/posttp.h"  /* typedefs: POSTSTATE, EVALUATE, SPLINES, etc.*/
/*----------------------------------------------------------------------------*/
# define LARGE_LOG_VAL ( 1.e+3 )
/*----------------------------------------------------------------------------*/
static SPLINES spl = {null};
static GRAPHICS gph = {null};
/*----------------------------------------------------------------------------*/
/* 'my_terminal' configuration: */

# if USE_NCURSES == 1
   # include <termcap.h>     /* terminal type header */
   static char *term;        /* terminal type string */ 

   # define CLSCREEN {\
     printf( "%s", tgetstr( "cl", null )); \
   }

   # define PRBLDCLR(a) {\
     printf( "%s%s", tgetstr( "md", null ), (a)); /* bold clear output */ \
   }

   # define PRINVERS(a) {\
     printf( "%s%s", tgetstr( "mr", null ), (a)); /* inverse */ \
   }

   # define PRNORMAL(a) {\
     printf( "%s%s", tgetstr( "me", null ), (a)); /* back to normal output */ \
   }
# else
   # define CLSCREEN {\
     printf( "\f" );\
   }

   # define PRBLDCLR(a) {\
     printf( "%s", (a));\
   }

   # define PRINVERS(a) {\
     printf( "%s", (a));\
   }

   # define PRNORMAL(a) {\
     printf( "%s", (a));\
   }
# endif
/*----------------------------------------------------------------------------*/
/* configure function poster(*) [ scheme 30D ]: */
/*...........................................................................*/
# define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X' ] */
/*...........................................................................*/
# define PST_MESH  0      /* 1: plot mesh cell system */
# define PST_WALLS 1      /* 1: plot E-walls, 2: plot E- and M-walls */
# define PST_ISOSCLE 1    /* 1: equally scaled axes */
/*...........................................................................*/
# define PST_TMPPLT  0    /* Disable/enable temperature plot mode [0/1/...] */
# define PST_PRSPLT  4    /* Disable/enable pressure plot mode [0/1/...] */
# define PST_FLWPLT  4    /* Disable/enable velocity plot mode [0/1/...] */
# define PST_PRSINIT 1
/*...........................................................................*/
# if CMPRSSBL == 0
/* Disable density plot mode */
   # define PST_DNSPLT  0
   # define PST_DNSINIT 0
# else /* if CMPRSSBL != 0 */
/* Disable/enable density plot mode [0/1/...] */
   # define PST_DNSPLT  4
   # define PST_DNSINIT 1  
# endif /* CMPRSSBL != 0 */
/*...........................................................................*/
/* scales: */
/*
# define PST_SCALE1 ( 3.000e-02 )
# define PST_SCALE2 ( 2.000e-03 )
*/
# define PST_SCALE1 ( 1.000e-01 )
# define PST_SCALE2 ( 2.000e-03 )
/*...........................................................................*/
# ifndef PST_SCLFLW /* scale of fluid flow [ maximum fluid velocity , e.g. ] */
# if MOD_DEFLT == 1
   # define PST_SCLFLW ( 3.000e+01 )
# elif MOD_DEFLT == 2
   # define PST_SCLFLW ( 1.000e+02 )
# elif MOD_DEFLT == 3
   # define PST_SCLFLW ( 1.000e+02 )
# elif MOD_DEFLT == 4
   # define PST_SCLFLW ( 5.000e+00 )
# else
   # define PST_SCLFLW ZERO
# endif /* MOD_DEFLT == ... */
# endif /* not defined PST_SCLFLW */
/*...........................................................................*/
# include "./tools/pstprc.fld" /* some postprocessing tools */
/*============================================================================*/

POSTSTATE *
modval( POSTSTATE *state )
{
/* allusions: */
/*
   extern SPLINES spl;
   extern GRAPHICS gph;
*/
/* declarations: */

   static FILE
      *evalfle = NULL,
      *timefle[MAX_PERIOD],
      *specfle[MAX_PERIOD];

   static FFT *fpt;
   static TXCNSL *csp;
   static SPLINES *spt = &spl;
   static GRAPHICS *gpt = &gph;
   static EVALUATE *vpt;

   static char     /* operation marks                                         */
      opm1 = null, /* opm1 = [0] 1: reference file [not] charged              */
      opm2 = null; /* opm2 = [0] 1: reference is [not] mean inc. wave         */

   static short
      j_ = null,
      k_ = null,
      ll = null,
      ind = null,
      lbl = null,
      frqlbl = null,

   envl[MAX_PERIOD] = {null};

   static long 
      hh = null,
      h_ = null,
      ii = null,
      i_ = null,
      jj = null,
      kk = null,
      mm = null,
      nn = null,
      pp = null,
      qq = null,
      init = null,
      final = null;

   static short 
     *idx = NULL;

   static char 
      ptr[STS_SIZE] = {null},
      tmeptr[STS_SIZE] = {null},
      spcfle[STS_SIZE] = {null},
      spcpfx[VSS_SIZE] = "spc.",
      option[SHS_SIZE] = "forward", /* Fourier transformation option */
      d_format[STS_SIZE] = {null},
     *s_format = "%s\n",
     *i_format = "%ld\n",
    **endp = null;

   static char 
      absc_unit[SHS_SIZE] = {null},
      ordn_unit[SHS_SIZE] = {null},
      t_type[MAX_PERIOD][SHS_SIZE] = {{null}},
      t_text[MAX_PERIOD][STS_SIZE] = {{null}},
      s_type[MAX_PERIOD][SHS_SIZE] = {{null}},
      s_text[MAX_PERIOD][STS_SIZE] = {{null}},
      fleptr[MAX_PERIOD][STS_SIZE] = {{null}};

   static double
      xx = ZERO,
      yy = ZERO,
      zz = ZERO,
      x1 = ZERO,
      x2 = ZERO,
      dx = ZERO,
      df = ZERO,
      norm = ZERO,
      real = ZERO,
      imag = ZERO,
      evl0 = ZERO,
      evl1 = ZERO,
      dlt0 = ZERO,
      dlt1 = ZERO,
      minm = ZERO,
      maxm = ZERO,
      xlower = ZERO,
      xupper = ZERO,
      frequency = ZERO;

   static double 
      frq[EVL_SPF] = {ZERO},
      vswr[EVL_SPF] = {ZERO},
      refl[EVL_SPF] = {ZERO},
      mean[EVL_SPF] = {ZERO},
      rref[EVL_SPF] = {ZERO},
      iref[EVL_SPF] = {ZERO},
      absmaxm[MAX_PERIOD] = {ZERO},
      timemax[MAX_PERIOD] = {ZERO},
      absmean[MAX_PERIOD] = {ZERO},
      meanenv[MAX_PERIOD] = {ZERO},
      maxspec[MAX_PERIOD] = {ZERO},
      freqmax[MAX_PERIOD] = {ZERO},
      modulation[MAX_PERIOD] = {ZERO},
      rspc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      ispc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      aspc[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      rcpl[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      icpl[MAX_PERIOD][EVL_SPF] = {{ZERO}},
      acpl[MAX_PERIOD][EVL_SPF] = {{ZERO}};

   static time_t
      nseconds = null,
     *timer = null;

/* system parameters and function prototypes */

   time_t
      time( time_t *timer );

   char
     *ctime( const time_t *timer );

# ifndef _CCBUG
   char
     *strcpy( char *ptr1, const char *ptr2 ),
     *strcat( char *ptr1, const char *ptr2 ),
     *strncat(char *ptr1, const char *ptr2, size_t n );
# endif

/* mathematical function prototypes */

   double 
      sqrt( double x ),
      log10( double x ),
      ceil( double x ),
      floor( double x );

/* user defined function prototypes: */

   FFT *
      fftrf( FFT *fpt );

   int 
      dspval( void );

   SPLINES
     *spline( SPLINES *spt );

   int 
      graphp( GRAPHICS *gpt );

   POSTSTATE 
     *pstprc( POSTSTATE *stp );

   char
     *lotos( long, char );

   TXCNSL 
     *txcnsl( TXCNSL *csp );

   EVALUATE
     *readval( POSTSTATE *state, char option );
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* get the terminal info: */

   term = ( char *) getenv( "TERM" ); /* get the terminal */

   ind = tgetent( null, term );

   if( ONE != ind )
   {
      fprintf( stderr, "Error on getting the termcap info\n" ); 
      exit( EXIT_FAILURE );
   };
# endif
/*............................................................................*/
/* string assignments [ initializations ], memory allocation etc. */
   
/* double format: */

   if ( null == strncmp( PLOT_FORMAT, "SPLINE" ,THREE ))
      strncpy( d_format ,"%+.15E%s", TEN );
   else
      strncpy( d_format , "%+.15e%s", TEN );

   strcpy( gph.format, PLOT_FORMAT );
/*............................................................................*/
/* initialize text console: */

   csp = txcnsl( null );
/*............................................................................*/
/* Error messages: */

# if LINUX_C_SNTX != 1
   hh = strlen ( LNGSTR );
   if ( hh < MAX_PERIOD * ( DEC + ONE ) ) /* <- length of long string LNGSTR  */
   {                                      /*    initializer in 'consts.h'     */
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD > %d !!!", 250/DEC );
      printf( "\n [ Filename pointer array 'fleptr[]' may perturb memory.]" );
   };
# endif

   mm = MAX_PERIOD;
   if ( FTR_NMBR < mm )                /* <- FTR_NMBR: macro in POSTER.CONF */
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD = %d !!!", MAX_PERIOD );
      printf( "\n [ Required is MAX_PERIOD <= %d = FTR_NMBR "
        "<- macro defined in FORMER.CONF.]", FTR_NMBR );
   };

   if ( FTR_NMBR < mm )                /* <- FTR_NMBR: macro in POSTER.CONF */
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Oversized macro MAX_PERIOD = %d !!!", MAX_PERIOD );
      printf( "\n [ Required is MAX_PERIOD <= %d = FTR_NMBR "
              "<- macro defined in 'dmnsnd.h'.]", FTR_NMBR );
   };
/*............................................................................*/
/* initialize pointers to structures [ of type EVALUATE and FFT ]: */

   vpt = ( state->vpt );
   fpt = ( vpt->fpt );
   idx = ( state->idx );
/*............................................................................*/
/* evaluation mode lbl = 1,2,... */

   if ( null == strncmp( vpt->exctyp, "HARMONIC", FIVE ))
      lbl = 4;
   else if ( null == strncmp( vpt->exctyp, "SMOOTH_HARMONIC", FIVE ))
      lbl = 4;
   else if ( null == strncmp( vpt->exctyp, "STEADY_STATE", FIVE ))
      lbl = 5;
   else
      lbl = 6;
/*............................................................................*/
   PRBLDCLR( "\n" );
   printf( " %*s", 78, "MODVAL" );
   PRNORMAL( "" );

   if (( ONE < lbl )
     &&( lbl < SEVEN ))
   {
/* enter special ports 1,2, ... */

      ii = null;
      while(( ii < ( state->fldprd ))
          &&( ii < MDV_SPCLPRTS )
          &&( ii < MAX_PERIOD ))
      {
         ii++;
         idx[ii] = ii;
      };
      idx[null] = ii;

      while (( ii < ( state->fldprd ))
           &&( ii < MAX_PERIOD ))
      {
         ii++;
         hh = ii - MDV_SPCLPRTS;
         h_ = hh - ONE;
/*............................................................................*/
# if MDV_AUTOEVAL == 1
         idx[ii] = ii;
# else
/*............................................................................*/
         if ( hh == ONE )
            printf( "\n VSWR evaluation: Please enter indices "
               "[ Escape: enter null ] >---->\n" );
        next_index:

         printf( " >----> enter %6ld. index ..................."
            "....................: ", hh );
         scanf( "%s", ptr );
         idx[ii] = ( short ) strtol( ptr, endp, DEC );
/*............................................................................*/
/* check if index is in domain: */

         if (( state->fldprd ) < idx[ii] )
         {
            printf( "\n This is not a Maxwell field port "
               "in selected file %s !!!", ( state->file ));
            if ( idx[ii] <= ( state->period ))
               printf( "\n [ It is a heat or fluid port that can be "
                  "separately evaluated.]" );
            printf( "\n\n >----> Enter new index >------"
               "------------------------------------>\n");
            goto next_index;
         };
# endif /* if MDV_AUTOEVAL != 1 */
/*............................................................................*/
         if ( null < idx[ii] )
         {
            if (( ONE < lbl )
              &&( lbl < FOUR ))
            {
               printf( " >----> enter name of file to be created"
                       " for index %-6d..........: ", idx[ii] );
               scanf( "%s", ptr );
               strcpy( fleptr[h_], ptr );
            };
            idx[null] = ii;

            if ( MAX_PERIOD <= idx[null] )
            {
               printf( "\n The last evaluation index %d has been "
                  "accepted.", idx[idx[null]] );
               printf( "\n [ Maximum number is %ld ",
                  ( long ) MAX_PERIOD );
               printf( "= macro MAX_PERIOD "
                  "in header '%s'.]\n", "POSTER.CONF" );

              goto read_values;
            };
         }
         else /* if ( idx[ii] <= null ) */
            goto read_values;
      }; 

     read_values:
/*............................................................................*/
      vpt = readval( state, 'r' );       /*                                   */
/*.....................................*/
   }; /* end if (( ONE < lbl )&&( lbl < SEVEN )) */
/*........... values dsc.val<n> read & copied for selected ports ...........*/









/*............................................................................*/
/* evaluation [ options labelled lbl ]: */

   switch ( lbl )
   {
     case null:
      return state;

     case 1:

      printf( "\n" );
      break;

     case 2: /* >------ save time response ---------------------------------> */

      printf( "\n\n Please wait a moment !" );
      printf( "\n [ Writing data into files %s ... ]\n",
         fleptr[null] );

      for ( hh=ONE; hh<=idx[null]; hh++ )
      {
         h_ = hh - ONE;

         strncpy( t_type[h_], ( vpt->name ), SHS_SIZE );
         strncpy( t_text[h_], ( vpt->text ), STS_SIZE );

         timefle[h_] = fopen( fleptr[h_], "w+" );

         fprintf( timefle[h_], s_format, t_type[h_] );
         fprintf( timefle[h_], "%s%s%d%s\n", t_text[h_], "_[evaluated"
            "_index", idx[hh], "]" );
         fprintf( timefle[h_], s_format, absc_unit );
         fprintf( timefle[h_], s_format, ordn_unit );

         fprintf( timefle[h_], d_format, ( fpt->t[null] ), "\n" );
         fprintf( timefle[h_], d_format, ( fpt->tt[null] ), "\n" );
         fprintf( timefle[h_], d_format, ( fpt->dt[null] ), "\n" );
         fprintf( timefle[h_], i_format, ( fpt->ttlg[null] ));

         for ( i_=null ; i_<( fpt->ttlg[null] ); i_++ )
         {
            fprintf( timefle[h_], d_format, ( fpt->r[hh][i_] ), "  " );
            fprintf( timefle[h_], d_format, ( fpt->i[hh][i_] ), "\n" );

# if LIST_FILE == 1
            printf( "  %ld:   %1.16e  + ( %1.16e )*j \n",
                    i_, ( fpt->r[hh][i_] ), ( fpt->i[hh][i_] ));
# endif
         };

         nseconds = time( timer );
         strcpy( tmeptr, ctime( &nseconds ));
         fprintf( timefle[h_], "\n%s%s%s%s\n", "DSC time response "
            "file '", fleptr[h_], "' terminated: ", tmeptr );
         printf( " %s%s%s%s", "DSC time response "
            "file '", fleptr[h_], "' terminated: ", tmeptr );

         fclose( timefle[h_] );
      };/* next hh */

      break;
/*...................... end case2 [ save time response ] ....................*/









     case 3: /* >------------ save spectral response -----------------------> */
     case 6: /* >---------- evaluate spectral response ---------------------> */

      if ( null < idx[null]  )
      {
         if (( lbl == SIX ) 
           &&( null != strncmp(( vpt->text ), "reference", NINE )))
         {                    /* enter external reference spectrum  'spc.ref' */
            strcpy( spcfle, spcpfx );
            strcat( spcfle, "ref" );

            printf( "\n" );

           open_reference1:

            evalfle = fopen( spcfle, "r+" );

            if ( evalfle == null )
            {
               printf( "\n Reference file %s not found "
                  "in present directory:", spcfle );
               printf( "\n Please re-enter filename [ Escape: "
                  "enter null ] >----> " );
               scanf( "%s", spcfle );

               if ( *spcfle == '0' )
                  goto frq_interval;
               else
                  goto open_reference1;
            }
            else if ( evalfle != null )
            {
               printf( "\n opened: reference spectrum file %s", spcfle ); 
               fscanf( evalfle, "%s", ptr );

               if ( null != strncmp( ptr, ( vpt->name ), THREE ))
               {
                  printf( "\n\n Error message from function %s :", __func__ );
                  printf( "\n\n Incompatible system identifier '%s'", ptr );
                  printf( "\n in reference spectrum, file %s !!!", spcfle );
                  printf( " [ overriding.]\n" );

                  fclose( evalfle );

                  goto frq_interval;
               };
       
               fscanf( evalfle, "%s", ptr );
               fscanf( evalfle, "%s", absc_unit );
               fscanf( evalfle, "%s", ordn_unit );

/* the domain lower bound: */
               fscanf( evalfle, "%s", ptr );
               xlower = strtod( ptr, endp ); /* the lower frequency bound */
	        
/* the domain upper bound: */
               fscanf( evalfle, "%s", ptr );
               xupper = strtod( ptr, endp ); /* the upper frequency bound */

/* the increment: */
               fscanf( evalfle, "%s", ptr );
/*............................................................................*/
# if MDV_SPLINTPL == 0
               df = strtod( ptr, endp );
# else
               dx = strtod( ptr, endp );
# endif
/*............................................................................*/
/* the number of sample points: */

               fscanf( evalfle, "%s", ptr );
               frqlbl = strtol( ptr, endp, DEC );

               if ( EVL_SPF < frqlbl )
                  frqlbl = EVL_SPF;

               if ( SPL_INTPL < frqlbl )
                  frqlbl = SPL_INTPL;
               
               jj = null;
               while( jj < frqlbl )
               {
/*............................................................................*/
# if EVL_WRTFREQ == 1
                  fscanf( evalfle, "%s", ptr ); /* read the frequency */
# endif
/*............................................................................*/
                  fscanf( evalfle, "%s", ptr ); /* the real part */
                  rref[jj] = strtod( ptr, endp );
                  fscanf( evalfle, "%s", ptr ); /* the imaginary part */
                  iref[jj] = strtod( ptr, endp );

                  jj++ ;   
               };
               printf( "\r entered: reference spectrum %s .      \n",
                  spcfle );
               fclose( evalfle );
               opm1 = ONE;
            }
            else /* case: no reference spectrum found or lbl != SIX */
               goto frq_interval;
         }
         else /* case: create reference spectrum */
              /* specify frequency domain */
         {
           frq_interval:

            opm1 = null;

            strncpy( ordn_unit, ( vpt->yunit ), SHS_SIZE );
            strcat( ordn_unit, "*" );
            strncat( ordn_unit, ( vpt->xunit ), SHS_SIZE );
            strcpy( absc_unit, "1/" );
            strcat( absc_unit, ( vpt->xunit ));

            printf( "\n\n Please enter frequency interval "
               "[ f1, f2 ] / GHz\n\n" );

            strcpy(( csp->rqfrm ), "points" );
            strcpy(( csp->rqdbl ), "Enter lower frequeny " );
            strcat(( csp->rqdbl ), ">>---------------------------> f1 / GHz" );
/*............................................................................*/
            csp = txcnsl( csp );        /* text console                       */
/*....................................*/
            xlower = ( csp->indbl );
            xlower *= ( 1.0e+09 );

            strcpy(( csp->rqfrm ), "points" );
            strcpy(( csp->rqdbl ), "Enter upper frequeny " );
            strcat(( csp->rqdbl ), ">>---------------------------> f2 / GHz" );
/*............................................................................*/
            csp = txcnsl( csp );        /* text console                       */
/*....................................*/
            xupper = ( csp->indbl );
            xupper *= ( 1.0e+09 );

            printf( "\n frequency interval: [ %1.5e , %1.5e ] "
               "GHz\n\n", ( 1.0e-09*xlower ), ( 1.0e-09*xupper ));

            strcpy(( csp->rqfrm ), "bracket" );
            strcpy(( csp->rqstr ), "Input correct >>---------" );
            strcat(( csp->rqstr ), "------------------------> [ Y/n ] ?" );
            strcpy(( csp->dfstr ), "y" );
/*............................................................................*/
            csp = txcnsl( csp );        /* text console                    */
/*....................................*/
            strcpy( ptr, ( csp->instr ));

            if (( *ptr == 'n' )||( *ptr == 'N' )) 
               goto frq_interval;
         };
      };
/* evaluation parameters [ reference spectrum, frequency domain ] specified.  */
/*............................................................................*/










/*............................................................................*/
/* Fourier transformations: */

      printf( "\n Please wait a moment !" );

      if ( lbl == THREE )
         printf( "\n [ Fourier transforms "
            "- writing data on files '%s' ... ]\n", fleptr[null] );
      if ( lbl == SIX )
         printf( "\n [ Fourier transforms ]\n" );

      kk = ( fpt->ttlg[null] );

      nn = TWO;
      while (( nn < kk )&&( nn < FTR_SIZE ))
         nn *= TWO;

      mm = nn / TWO;
      pp = ( nn - kk ) / TWO;
      qq = pp + kk;

      hh = ONE; 
      while( hh <= idx[null] )
      {
         h_ = hh - ONE;

         ii = nn - ONE;
         while ( qq <= ii )
         {
            ( fpt->r[hh][ii] ) = ZERO;
            ( fpt->i[hh][ii] ) = ZERO;
            ii-- ;
         };
         while ( pp <= ii )
         {
            ( fpt->r[hh][ii] ) = ( fpt->r[hh][ii-pp] );
            ( fpt->i[hh][ii] ) = ( fpt->i[hh][ii-pp] );
            ii-- ;
         };
         while ( null <= ii )
         {
            ( fpt->r[hh][ii] ) = ZERO;
            ( fpt->i[hh][ii] ) = ZERO;
            ii-- ;
         };
/*............................................................................*/
# if READ_REVERSE_ == 1 /* interchange upper [negative] with lower [positive] */
                        /* frequencies: */
         ii = pp; do
         {
            jj = ii + mm;

            xx = ( fpt->r[hh][jj] );
            yy = ( fpt->i[hh][jj] );
	    
            ( fpt->r[hh][jj] ) = ( fpt->r[hh][ii] );
            ( fpt->i[hh][jj] ) = ( fpt->i[hh][ii] );

            ( fpt->r[hh][ii] ) = xx;
            ( fpt->i[hh][ii] ) = yy;

            ii++ ;
         }  while ( ii < mm );

         ( fpt->t[hh] ) = - mm*( fpt->dt[null] );
# else
         ( fpt->t[hh] ) = ZERO;
# endif /* READ_REVERSE */
/*............................................................................*/
         ( fpt->tt[hh] ) = ( fpt->t[hh] ) + nn*( fpt->dt[null] );
         ( fpt->dt[hh] ) = ( fpt->dt[null] );
         ( fpt->ttlg[hh] ) = nn;
         ( fpt->mult[hh] ) = ONE;
         ( fpt->p ) = hh;
         ( fpt->q ) = hh;
         strcpy( fpt->opt, option);
/*............................................................................*/
         fpt = fftrf( fpt );      /* Fast Fourier transform                   */
/*..............................*/
         hh++ ;
      }; /* next hh */

      dx = ( fpt->ds[null] );
/*............................................................................*/
# if WRITE_REVERSE == 1
      ii = mm;
      xx = ( fpt->s[null] );
# else
      ii = null;
      xx = ZERO;
# endif
/*............................................................................*/
      if (( MDV_SPLINTPL == 0 )
        &&( opm1 == ONE ))
      {
         x1 = xlower + dx/2.;
         x2 = xupper - dx/2.;
      }
      else
      {
         x1 = xlower;
         x2 = xupper;
      };

      if( xx < x1 )
         kk = ( long ) floor (( x1 - xx ) / dx ); /* kk = largest integer not */
      else                                        /* greater than argument */
         kk = null;

      if( nn < kk )
         kk = nn;
      
      init = ii + kk;
      x1 = xx + kk*dx; /* [1]; cf. remark [2], below */

      if( xx < x2 )
         kk = ( long ) ceil (( x2 - xx ) / dx ); /* kk = smallest integer not */
      else                                       /* less than argument */
         kk = null;

      if( nn < kk )
         kk = nn;

      final = ii + kk + ONE;
      x2 = xx + kk*dx;

      kk = final - init;

/* Fourier transformations terminated */
/*............................................................................*/










/*............................................................................*/
      if ( lbl == THREE ) /* save spectral response */ 
      {
         hh = ONE;
         while( hh <= idx[null] )
         {
            h_ = hh - ONE;

            strncpy( s_type[h_], ( vpt->name ), SHS_SIZE );
            strncpy( s_text[h_], ( vpt->text ), STS_SIZE );

            specfle[h_] = fopen( fleptr[h_], "w+" );

            fprintf( specfle[h_], s_format, s_type[h_] );
            fprintf( specfle[h_], "%s%s%d%s\n", s_text[h_],"_[evaluated"
               "_index", idx[hh], "]" );
            fprintf( specfle[h_], s_format, absc_unit );
            fprintf( specfle[h_], s_format, ordn_unit );
            fprintf( specfle[h_], d_format, x1, "\n" );
            fprintf( specfle[h_], d_format, x2, "\n" );
            fprintf( specfle[h_], d_format, dx, "\n" );
            fprintf( specfle[h_], i_format, kk );

# if EVL_WRTFREQ == 1
            frequency = x1;
# endif
            jj = init;
            ii = init;
            while ( ii < final )
            {
               if ( nn <= jj )
                  jj -= nn;
 
# if EVL_WRTFREQ == 1
               fprintf( specfle[h_], d_format, frequency, "  " );
	       frequency += dx;
# endif
               fprintf( specfle[h_], d_format, ( fpt->r[hh][jj] ), "  " );
               fprintf( specfle[h_], d_format, ( fpt->i[hh][jj] ), "\n" );

# if LIST_FILE == 1
               printf( "  %d:   %1.16le  + ( %1.16le )*j \n",
                  jj, ( fpt->r[hh][jj] ), ( fpt->i[hh][jj] ));
# endif
               jj++ ;
               ii++ ;
            };

            nseconds = time( timer );
            strcpy( tmeptr, ctime( &nseconds ));
            fprintf( specfle[h_], "\n%s%s%s%s\n", "DSC "
               "spectrum file '", fleptr[h_], "' terminated: ", tmeptr );

            printf( "\n\n %s%s%s%s", "DSC "
               "spectrum file '", fleptr[h_], "' terminated: ", tmeptr );

            fclose( specfle[h_] );

            hh++ ;
         };
         break;
      };/* end if lbl == THREE */
/*............................................................................*/
/* spectral response saved */
/*............................................................................*/










/*............................................................................*/
/* case 6, cont'd [ spectral response, spline interpolated ]:                 */
/*............................................................................*/
# if MDV_SPLINTPL == 0
/*
      printf( "\n frqlbl = %.7d, kk = %.7d", frqlbl, kk );
      printf( "\n xlower = %.15le , x1 = %.15le", xlower, x1 );
      printf( "\n xupper = %.15le , x2 = %.15le", xupper, x2 );
      printf( "\n df     = %.15le , dx = %.15le", df, dx );
      scanf( "%s", ptr );
*/
      if( opm1 == ONE )
      {
         if(( frqlbl != kk )||
            ( 1.e-7 < fabs( ( df - dx ) / dx ))||
            ( 1.e-7 < fabs( ( xlower - x1 ) / dx ))||
            ( 1.e-7 < fabs( ( xupper - x2 ) / dx )))
         {
            printf( "\n\n Error message from function %s :", __func__ );
            printf( "\n\n Evaluated file does not match reference file" );
            printf( " !!!\n " );
            break;
         };
      };

      frqlbl = kk;
      xlower = x1;
      xupper = x2;

# elif MDV_SPLINTPL < 10

      if( opm1 == null ) /* no reference file entered [ yet generated, e.g. ] */
      { 
         frqlbl = kk*MDV_SPLINTPL;
         if ( EVL_SPF < frqlbl )
            frqlbl = EVL_SPF;
      };

      dx = ( xupper - xlower ) / ( frqlbl - ONE );

# elif MDV_SPLINTPL >= 30

      if( opm1 == null ) /* no reference file entered [ yet generated, e.g. ] */
      {
         frqlbl = MDV_SPLINTPL;
         if ( EVL_SPF < frqlbl )
            frqlbl = EVL_SPF;
      };

      dx = ( xupper - xlower ) / ( frqlbl - ONE );
# endif
/*............................................................................*/
      if ( SPL_SPNTS < kk )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Too many supp. points in frequency domain !!!" );
         printf( "\n [ Number %ld exceeds maximum number of spline "
            "interpolation points %ld", kk, (long) SPL_SPNTS );
         printf( "\n   = macro SPL_INTPL in function %s.", "spline(*)" );
         printf( "\n   - Change macro in compliance with memory "
            "resources.]\n" );
         exit( EXIT_FAILURE );
      }
      else if ( SPL_INTPL < frqlbl )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Too many number of frequency points !!!" );
         printf( "\n [ Number %ld exceeds maximum number of spline "
            "interpolation points %ld", ( long ) frqlbl, (long) SPL_INTPL );
         printf( "\n   = macro SPL_INTPL in function %s.", "spline(*)" );
         printf( "\n   - Change macro in compliance with memory "
            "resources.]\n" );
         exit( EXIT_FAILURE );
      }
      else if ( EVL_SPF < frqlbl )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Too many number of frequency points !!!" );
         printf( "\n [ Number %ld exceeds maximum number %ld",
            ( long ) frqlbl, (long) EVL_SPF );
         printf( "\n   = macro EVL_SPF in function %s.", __func__ );
         printf( "\n   - Change macro in compliance with memory "
            "resources.]\n" );
         exit( EXIT_FAILURE );
      };

      if( xlower < x1 ) /* [2]: this can only happen due to roundoff error in */
         xlower = x1;   /* formula [1].                                       */
      if( x2 < xupper ) /* [ analogous statement ]                            */
         xupper = x2;

      df = ( xupper - xlower )/( frqlbl - ONE );
      frequency = xlower;

      ii = ONE;
      jj = frqlbl-ONE;
      frq[null] = xlower;
      spl.dmn[null] = xlower;
      while ( ii < jj ) 
      {
         frequency += df;
         frq[ii] = frequency;
         spl.dmn[ii] = frequency;
         ii++;
      };                    /* [3]; spl.dmn[] must lie strictly within the    */
      frq[ii] = xupper;     /* closed interval                                */
      spl.dmn[ii] = xupper; /* [ spl.vct[0][0], spl.vct[frqlbl-ONE][0] ],     */
                            /* which ist granted by virtue of [2] above       */
                            /* and condition [4], below.                      */
/*............................................................................*/
      hh = ONE;
      while( hh <= idx[null] )
      {
         h_ = hh - ONE;
         
/* spectrum, port hh [real part]: */

         frequency = x1;
         jj = init;
         ii = init;
         while ( ii < final )
         {
            if ( nn <= jj )
               jj -= nn;

            spl.vct[ii-init][null] = frequency;
            spl.vct[ii-init][ONE] = ( fpt->r[hh][jj] );

            frequency += ( fpt->ds[hh] );

            jj++ ;
            ii++ ;
         };
         spl.vct[kk-ONE][null] = x2; /* [4]: cf. comment [3] above */
/*............................................................................*/
# if MDV_SPLINTPL == 0

         j_ = null; 
         while( j_ < frqlbl )
         {
            rspc[h_][j_] = spl.vct[j_][ONE];
            j_++ ;
         };
# else
         spl.nn = frqlbl;
         spl.mm = kk;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation:                 */
/*.................................*//* real part                             */
         if (( spt->rtn ) == ONE )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on calling spline function" );
            printf( "\n for computing rspc[%ld][] !!!", h_ );
            printf( "\n [ program stopped.]\n" );

            exit( EXIT_FAILURE );
         };

         j_ = null; 
         while( j_ < frqlbl )
         {
            rspc[h_][j_] = spl.fct[j_];
            j_++ ;
         };
# endif
/*............................................................................*/
/* spectrum, port hh [maginary part]: */

         jj = init;
         ii = init;
         while ( ii < final )
         {
            if ( nn <= jj )
               jj -= nn;
               
            spl.vct[ii-init][ONE] = ( fpt->i[hh][jj] );

            jj++ ;
            ii++ ;
         };
/*............................................................................*/
# if MDV_SPLINTPL == 0

         j_ = null;
         while( j_ < frqlbl )
         {
            ispc[h_][j_] = spl.vct[j_][ONE];
            j_++ ;
         };
# else
         spl.nn = frqlbl;
         spl.mm = kk;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation:                 */
/*.................................*//* imaginary part                        */
         if (( spt->rtn ) == ONE )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on calling spline function" );
            printf( "\n for computing ispc[%ld][] !!!", h_ );
            printf( "\n [ program stopped.]\n" );

            exit( EXIT_FAILURE );
         };

         j_ = null; 
         while ( j_ < frqlbl )
         {
            ispc[h_][j_] = spl.fct[j_];
            j_++ ;
         };
# endif /* MDV_SPLINTPL != 0 */
/*............................................................................*/
/* spectrum, port hh [absolute value]: */

         jj = init;
         ii = init;
         while ( ii < final )
         {
            if ( nn <= jj )
               jj -= nn;

            real = ( fpt->r[hh][jj] );
            imag = ( fpt->i[hh][jj] );
            spl.vct[ii-init][ONE] = sqrt( real*real + imag*imag );

            jj++ ;
            ii++ ;
         };

         maxspec[h_] = ZERO;
         freqmax[h_] = ZERO;
/*............................................................................*/
# if MDV_SPLINTPL == 0

         j_ = null;
         while ( j_ < frqlbl )
         {
            norm = fabs( spl.vct[j_][ONE] );
            aspc[h_][j_] = norm;

            if ( maxspec[h_] < norm )
            {
               maxspec[h_] = norm;
               freqmax[h_] = spl.vct[j_][null];
            };
            j_++ ;
         };
# else
         spl.nn = frqlbl;
         spl.mm = kk;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation:                 */
/*.................................*//* absolute value                        */
         if (( spt->rtn ) == ONE )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on calling spline function" );
            printf( "\n for computing aspc[%ld][] !!!", h_ );
            printf( "\n [ program stopped.]\n" );

            exit( EXIT_FAILURE );
         };

         j_ = null; 
         while ( j_ < frqlbl )
         {
            norm = fabs( spl.fct[j_] );
            aspc[h_][j_] = norm;

            if ( maxspec[h_] < norm )
            {
               maxspec[h_] = norm;
               freqmax[h_] = spl.dmn[j_];
            };
            j_++ ;
         };
# endif /* MDV_SPLINTPL != 0 */
/*............................................................................*/
         freqmax[h_] /= 1.e+09; /* in GHz */

         hh++ ;

      };/* next hh */
/*......................... spectral response ready ..........................*/










/*............................................................................*/
/* generate reference spectrum: */
/*
      printf( "\n\n                                   "
              "                                  ) <- ?" );
      printf( "\r Generate reference spectrum ? >----"
              "----------------> [ y/n ] >--> (" );
      scanf( "%s", ptr );

      if (( *ptr == 'y' )||( *ptr == 'Y' ))
      {
*/ /* } */
      printf( "\n" );

      if ( null == strncmp(( vpt->text ), "reference", NINE ))
      {
         strcpy( spcfle, spcpfx );
         strcat( spcfle, "ref" );
         strcat( spcfle, ( state->flbl ));

         evalfle = fopen( spcfle, "w+" );

         if ( evalfle == null )
         {
            printf( "\n Error on opening reference spectrum file"
               " %s !", spcfle );
            printf( "\n [overriding: spectrum not saved.]\n " );
            break;
         };

         fprintf( evalfle, s_format, ( vpt->name ));
         fprintf( evalfle, s_format, "reference_spectrum" );
         fprintf( evalfle, s_format, absc_unit );
         fprintf( evalfle, s_format, ordn_unit );
         fprintf( evalfle, d_format, xlower, "\n" );
         fprintf( evalfle, d_format, xupper, "\n" );
         fprintf( evalfle, d_format, dx, "\n" );
         fprintf( evalfle, i_format, frqlbl );

# if EVL_WRTFREQ == 1
	 frequency = xlower;
# endif
         j_ = null;
         while ( j_ < frqlbl )
         {
# if EVL_WRTFREQ == 1
            fprintf( evalfle, d_format, frequency, "  " );
	    frequency += dx;
# endif
            fprintf( evalfle, d_format, rspc[null][j_], "  " );
            fprintf( evalfle, d_format, ispc[null][j_], "\n" );
            j_++ ;
         };

         nseconds = time(timer);
         strcpy( tmeptr, ctime( &nseconds ));

         fprintf( evalfle, "\n# SPECTRUM file %s\n", spcfle );
         fprintf( evalfle, "# created:%s", tmeptr );

         fclose( evalfle );

         printf( "\n" );
         printf( CLEAR_LINE );
         printf( "\r REFSPC.file %s ", spcfle );
         printf( "created: %.24s ", tmeptr );

         if ( idx[null] <= MDV_SPCLPRTS )
            break;
      }; /* end if ( null == strncmp( vpt->text, "reference", NINE )) */
/*............................................................................*/
/* reference spectrum ready [ and stored ] */
/*............................................................................*/











/*............................................................................*/
/* display on screen [ console ]: */

      printf( "\n\n Special evaluated spectral parameters "
         "of file %s :", ( state->file ));
      printf( "\n\n  index|abs. maximum :at frequency " );

      for ( j_= null; (( j_< frqlbl )
        &&( j_< THREE )); j_++ )
         printf( "|ampl.at GHz->" );

      printf( "\n       |[%.13s]       [%.3s]", ordn_unit, "GHz" );

      for ( j_= null; (( j_< frqlbl )
        &&( j_< THREE )); j_++ )
         printf( "|%.7e", frq[j_]/1.e+09 );

      printf( "\n --------------------------------------"
              "----------------------------------------" );

      hh = ONE;
      while( hh <= idx[null] )
      {
         h_ = hh - ONE;
         printf( "\n %6d|%.7e:%.7e", idx[hh], maxspec[h_], freqmax[h_] );

         for ( j_=null ; (( j_< frqlbl )&&( j_< THREE )) ; j_++ )
            printf( "|%.7e", aspc[h_][j_] );

         hh++ ;
      };/* next hh */

      jj = ( int ) (( frqlbl - THREE + ( FIVE - ONE )) / FIVE );

      ii = null;
      while ( ii < jj )
      {
         j_ = THREE + ii*FIVE;
         if ( j_ < frqlbl )
         {
            printf( "\n\n  index" );
            while (( j_ < THREE + ( ii + ONE )*FIVE )&&( j_ < frqlbl ))
            {
               printf( "|ampl.at GHz->" );
               j_++;
            };
         };
         j_ = THREE + ii*FIVE;
         if ( j_ < frqlbl )
         {
            printf( "\n       " );
            while (( j_ < THREE + ( ii + ONE )*FIVE )&&( j_ < frqlbl ))
            {
               printf( "|%.7e", ( 1.e-09*frq[j_] ));
               j_++;
            };
         };
         printf( "\n -----------------------------------"
            "-------------------------------------------" );

         hh = ONE;
         while( hh <= idx[null] )
         {
            h_ = hh - ONE;
            j_ = THREE + ii*FIVE;
            if ( j_ < frqlbl )
            {
               printf( "\n %6d", idx[hh] );
               while (( j_ < THREE + ( ii + ONE )*FIVE )&&( j_ < frqlbl ))
               {
                  printf( "|%.7e", aspc[h_][j_] );
                  j_++;
               };
            };
            hh++ ;
         };
         ii++;
      };/* end while ii < jj */

      printf( "\n -----------------------------------"
         "-------------------------------------------" );
/*............................ display ready ................................*/










/*............................................................................*/
      opm2 = ONE;

      if (( opm1 == null )
        &&( MDV_SPCLPRTS < idx[null] ))  
      {
         printf( "\n\n normalizing S-parameters to absolute value"
            "\n of incident wave amplitude.\n " );
      }
      else if (( opm1 == null )
             &&( idx[null] <= MDV_SPCLPRTS ))
      {
         opm2 = null;
         printf( "\n\n normalizing S-parameters relative to port1.\n " ); 
      };
         
/* = number of vswr evaluation points: */

      kk = idx[null] - MDV_SPCLPRTS;

/* number of spline interpolataion points between support points */
/* spl.vct[][null]: */

      if ( null < kk )
      {
         i_ = MDV_GPHINTPL;
   
         if ( GPH_POINTS <= kk * ( i_ + ONE ) ) 
            i_ = ( short ) (( double)( GPH_POINTS / kk ) - ONE );
      };

      j_ = null; 
      while ( j_ < frqlbl )
      {
         if ( null < kk ) /* i.e. if ( MDV_SPCLPRTS < idx[null] ) */ 
         {
            ii = null;
            hh = ONE;
            h_ = null;
            while( hh < kk )
            {
               mm = hh + MDV_SPCLPRTS;
               nn = hh + MDV_SPCLPRTS + ONE;
               spl.vct[h_][null] = idx[mm];
               spl.vct[h_][ONE]  = aspc[mm-ONE][j_];
               spl.dmn[ii] = idx[mm];
               ii++;
               for ( k_= ONE; k_<= i_ ; k_++ )
               {
                  x2 = ( double ) k_ / ( i_ + ONE ) ;
                  x1 = 1. - x2;
                  spl.dmn[ii]  = x1*idx[mm] + x2*idx[nn];
                  ii++;
               };
               h_ = hh;
               hh++ ;
            };/* next hh */
            spl.vct[h_][null] = idx[idx[null]];
            spl.vct[h_][ONE]  = aspc[idx[null]-ONE][j_];
            spl.dmn[ii] = idx[idx[null]];
            ii++ ;

            spl.nn = ii;
            spl.mm = kk;
/*............................................................................*/
            spt = spline( spt );        /* spline interpolation: stand. wave  */
/*....................................*//* over pnts.idx, at frequency frq[j_]*/
            if (( spt->rtn ) == ONE )
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling spline function !!!" );
               printf( "\n [ program stopped.]\n" );
	       
               exit( EXIT_FAILURE );
            };

            minm = +1.e+277;
            maxm = - minm;

            jj = null;
            while( jj < ii )
            {
               if ( fabs( spl.fct[jj] ) < minm )
                  minm = fabs( spl.fct[jj] );
               if ( maxm < fabs( spl.fct[jj] ) )
                  maxm = fabs( spl.fct[jj] );
               jj++ ;
            };

            mean[j_] = .5*( minm + maxm ); /* mean incident wave amplitude */
/*............................................................................*/
# if MDV_PLOT_SWR == 1

            if ( j_ == null )
            {
               jj = null;
               while ( jj < ii )
               {
                  gph.vct[jj][null] = spl.dmn[jj];
                  gph.vct[jj][ONE]  = spl.fct[jj];
                  jj++ ;
               };
               strcpy( gph.file, "stw" );
               strcat( gph.file, lotos( j_, null ));
               strcat( gph.file, state->flbl );
               strcpy( gph.name, state->name );
               strcpy( gph.text, "standing_wave" );
               strcpy( gph.xunit, "index" );
               strcpy( gph.yunit, "volts" );
               gph.nn = ii;
/*............................................................................*/
               ind = graphp( gpt );               /* graphics file: std.wave  */
/*..............................................*//* at frequency frq[j_]     */
            };
# endif /* MDV_PLOT_SWR == 1 */
/*............................................................................*/
         }
         else /* if ( idx[null] <= MDV_SPCLPRTS ) *//* normalization to port1 */
            mean[j_] = aspc[null][j_];

         if ( opm1 == null ) /* absent reference spectrum:  */
                             /* normalize to incident wave  */
                             /* [ mean absolute amplitude ] */
         {
            real = mean[j_];
            imag = ZERO;
            norm = real*real + imag*imag; 

            if ( null < kk ) /* i.e standing wave ports to be evaluated       */
            {                
               if ( 1.e-277 < minm )
               {
                  yy = maxm/minm;
                  vswr[j_] = yy;
                  zz = ( yy - 1. )/( yy + 1. );
                  refl[j_] = 100.*zz;
               }
               else
               {
                   vswr[j_] = HUGE_VALF;
                   refl[j_] = 100.; 
               };
/*  
               rspc[null][j_] -= real;             ???  - inspection !!!    
               ispc[null][j_] -= imag;             
*/
            };
         }
         else if ( opm1 == ONE )  /* normalization to reference spectrum      */
         {
            real = rref[j_];
            imag = iref[j_];
            norm = real*real + imag*imag;

            rspc[null][j_] -= real;
            ispc[null][j_] -= imag;
            xx = ( real*rspc[null][j_] + imag*ispc[null][j_] ) / norm;
            yy = ( real*ispc[null][j_] - imag*rspc[null][j_] ) / norm;
            zz = sqrt( xx*xx + yy*yy );

            if ( 1.e-277 < fabs( 1. - zz ) )
            {
               vswr[j_] = ( 1. + zz )/( 1. - zz );
               refl[j_] = 100.*zz;
            } 
            else
            {
               vswr[j_] = HUGE_VALF;
               refl[j_] = 100.;
            };
         };

         ii = MDV_SPCLPRTS;

         if ( idx[null] < ii )
            ii = idx[null]; 

         h_= null;
         while( h_< ii )
         {
            xx = ( real*rspc[h_][j_] + imag*ispc[h_][j_] ) / norm;
            yy = ( real*ispc[h_][j_] - imag*rspc[h_][j_] ) / norm;

            rcpl[h_][j_] = xx;
            icpl[h_][j_] = yy;

            zz = xx*xx + yy*yy; 

            if ( 1.e-277 < zz )
            {
/*............................................................................*/
# if MDV_NORMALZE == 1 /* normalize S-parameters to abs.value one */

               zz = sqrt( zz ); 
               rcpl[h_][j_] /= zz;
               icpl[h_][j_] /= zz;

               acpl[h_][j_] = 0.;
# else
               acpl[h_][j_] = 10.*log10( zz ); 
# endif
/*............................................................................*/
            }
            else
               acpl[h_][j_] = - LARGE_LOG_VAL;
            h_++ ;
         };
         j_++ ;
      };/* while ( j_ < frqlbl ) */
      
/* s-parameters ready */
/*............................................................................*/
/* save S-parameters: */

      ii = MDV_SPCLPRTS;

      if ( idx[null] < ii )
         ii = idx[null];
       
      h_ = null;
      while( h_ < ii )
      {
         strcpy( ptr, "s" );
         strcat( ptr, lotos ( h_+ONE, null ) );
         strcat( ptr, "1_" );
         strcat( ptr, state->flbl );
         
         strcpy( spcfle, spcpfx );
         strcat( spcfle, ptr );

         evalfle = fopen( spcfle, "w+" );

         jj = strlen( ptr );
         strcat( ptr, "_<<" );
         strncat( ptr, vpt->text, ( STS_SIZE - SIX - jj ));
         strcat( ptr, ">>" );
              
         fprintf( evalfle, s_format, ( vpt->name ));
         fprintf( evalfle, s_format, ptr );
         fprintf( evalfle, s_format, absc_unit );
         fprintf( evalfle, s_format, "---" );
         fprintf( evalfle, d_format, xlower, "\n" );
         fprintf( evalfle, d_format, xupper, "\n" );
         fprintf( evalfle, d_format, dx, "\n" );
         fprintf( evalfle, i_format, frqlbl );
/*............................................................................*/
# if EVL_WRTFREQ == 1
	 frequency = xlower;
# endif
/*............................................................................*/
         j_ = null;
         while ( j_ < frqlbl )
         {
/*............................................................................*/
# if EVL_WRTFREQ == 1
            fprintf( evalfle, d_format, frequency, "  " );
	    frequency += dx;
# endif
/*............................................................................*/
            fprintf( evalfle, d_format, rcpl[h_][j_], "  " );
            fprintf( evalfle, d_format, icpl[h_][j_], "\n" );
            j_++ ;
         };

         nseconds = time(timer);
         strcpy( tmeptr, ctime( &nseconds ));

         fprintf( evalfle, "\n# SPECTRUM file %s\n", spcfle );
         fprintf( evalfle, "# created:%s", tmeptr );

         fclose( evalfle );

         printf( "\n" );
         printf( CLEAR_LINE );
         printf( "\r SPECTR.file %s", spcfle );
         printf( " created: %.24s ", tmeptr );

         h_++ ;
      };
/*............................................................................*/
      i_ = MDV_GPHINTPL;         /* number of interpoltation points           */
                                 /* between frequencies frq                   */
      if ( GPH_POINTS <= frqlbl * ( i_ + ONE ) )
             i_ = ( short )(( double )( GPH_POINTS / frqlbl ) - ONE );

      ii = null;
      j_ = null;
      while( j_< frqlbl-ONE )
      {
         spl.vct[j_][null] = frq[j_];

         if ( opm2 == ONE )
            spl.vct[j_][ONE] = vswr[j_];

         spl.dmn[ii] = frq[j_];
         ii++;
         for ( k_= ONE ; k_<= i_ ; k_++ )
         {
            x2 = ( double ) k_ / ( i_ + ONE );
            x1 = 1. - x2;
            spl.dmn[ii]  = x1*frq[j_] + x2*frq[j_+ONE];
            ii++ ;
         };
         j_++ ;
      }; /* until j_ = frqlbl - ONE */
      spl.vct[j_][null] = frq[j_];

      if ( opm2 == ONE )
         spl.vct[j_][ONE] = vswr[j_];

      spl.dmn[ii] = frq[j_];
      ii++;

      if ( opm2 == ONE )
      {
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation: VSWR            */
/*.................................*/
         ii = null;
         j_ = null;
         while( j_ < frqlbl-ONE )
         {
            spl.vct[j_][ONE] = refl[j_];
            gph.vct[ii][null] = spl.dmn[ii];
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][null] = spl.dmn[ii];
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++;
         }; /* until j_ = frqlbl - ONE */
         spl.vct[j_][ONE] = refl[j_];
         gph.vct[ii][null] = spl.dmn[ii];
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++ ;

         strcpy( gph.file, "swr" );
         strcat( gph.file, state->flbl );
         strcpy( gph.name, state->name );
         strcpy( gph.text, "voltage_standing_wave_ratio " );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "---" );
         gph.nn = ii;

         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );              /* graphics file: VSWR             */
                                          /*                                  */
         spt = spline( spt );            /* spline interpolation: reflexion   */
/*.....................................*/
         ii = null;
         j_ = null; 
         while( j_ < frqlbl-ONE )
         {
            zz = fabs( refl[j_]/100. );

            if ( 0. < zz )
               spl.vct[j_][ONE] = 20.*log10( zz );
            else                                  /* Replacing LARGE_LOG_VAL  */
               spl.vct[j_][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                                  /* made spline(*) function  */
                                                  /* inoperative              */
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++;
         }; /* until j_ = frqlbl - ONE */

         zz = fabs( refl[j_]/100. );

         if ( 0. < zz )
            spl.vct[j_][ONE] = 20.*log10( zz );
         else                                  /* Replacing LARGE_LOG_VAL  */
            spl.vct[j_][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                               /* made spline(*) function  */
                                               /* inoperative              */
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++ ;

         strcpy( gph.file, "ref" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "reflexion " );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "0/0" );

         gph.nn = ii;
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );              /* graphics file: reflex.[%], port1*/
                                          /*                                  */
         spt = spline( spt );            /* spline interpolation: refl.[dB],p1*/
/*.....................................*/
         ii = null;
         j_ = null;
         while( j_ < frqlbl-ONE )
         {
            zz = fabs( refl[j_]/100. );

            if( zz < 1. )
               spl.vct[j_][ONE] = 10.*log10( 1. - zz*zz );
            else                                  /* Replacing LARGE_LOG_VAL  */
               spl.vct[j_][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                                  /* made spline(*) function  */
                                                  /* inoperative              */
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++;
         }; /* until j_ = frqlbl - ONE */

         zz = fabs( refl[j_]/100. );

         if ( zz < 1. )
            spl.vct[j_][ONE] = 10.*log10( 1. - zz*zz );
         else                                  /* Replacing LARGE_LOG_VAL  */
            spl.vct[j_][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                               /* makes spline(*) function */
                                               /* inoperative              */
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++ ;

         strcpy( gph.file, "rtl" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "return_loss" );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "dB" );

         gph.nn = ii;
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );              /* graphics file: return loss [dB] */
                                          /*                                  */
         spt = spline( spt );            /* spline intpl.: insert.loss [dB],p1*/
/*.....................................*/
      }; /* end if ( opm2 == ONE ) */

      ii = null;
      j_ = null; 
      while( j_ < frqlbl-ONE )
      {
         if ( null < MDV_SPCLPRTS )
            spl.vct[j_][ONE] = acpl[null][j_];

         if ( opm2 == ONE ) 
         {
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ; 
            };
         };
         j_++ ;
      }; /* until j_ = frqlbl - ONE */

      if ( null < MDV_SPCLPRTS )
         spl.vct[j_][ONE] = acpl[null][j_];

      if ( opm2 == ONE ) 
      {
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++ ;

         strcpy( gph.file, "isl" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "insertion_loss " );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "dB" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );              /* graphics file: insert.loss [dB] */
/*.......................................*/
      }; /* end if ( opm2 == ONE ) */
/*....... reflexion, VSWR,..., insertion loss etc. now have been saved .......*/










/*............................................................................*/
/* plot s-parameters, special ports: [ in PLOT_FORMAT, e.g. "SPLINE" ] */

      ll = null;
      while (( ll < MDV_SPCLPRTS )
           &&( ll < idx[null] ))
      {
         strcpy( ptr, lotos( ll+ONE, null ));

         spl.mm = frqlbl;
         spl.nn = ii;
/*............................................................................*/
         spt = spline( spt );        /* spline interpolation: |sn1|           */
/*.................................*/ 
         ii = null;
         j_ = null;
         while ( j_ < frqlbl-ONE )
         {
            spl.vct[j_][ONE] = rcpl[ll][j_];
            gph.vct[ii][null] = spl.dmn[ii];
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++ ;
         }; /* until j_ = frqlbl-ONE */

         spl.vct[j_][ONE] = rcpl[ll][j_];
         gph.vct[ii][null] = spl.dmn[ii];
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++;

         strcpy( gph.file, "s" );
         strcat( gph.file, ptr ); 
         strcat( gph.file, "1_" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "|S" );
         strcat( gph.text, ptr );
         strcat( gph.text, "1|" );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "dB" );

         gph.nn = ii;
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );                /* graphics file: |sn1|          */
                                            /*                                */
         spt = spline( spt );              /* spline interpolation: real(sn1) */
/*.......................................*/
         ii = null;
         j_ = null;
         while ( j_ < frqlbl-ONE )
         {
            spl.vct[j_][ONE] = icpl[ll][j_];
            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_ ++;
         }; /* until j_ = frqlbl-ONE */
         spl.vct[j_][ONE] = icpl[ll][j_];
         gph.vct[ii][ONE] = spl.fct[ii];
         ii++;

         strcpy( gph.file, "r" );
         strcat( gph.file, ptr );
         strcat( gph.file, "1_" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "real(s" );
         strcat( gph.text, ptr );
         strcat( gph.text, "1)" );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "--" );

         gph.nn = ii;
         spl.nn = ii;
         spl.mm = frqlbl;
/*............................................................................*/
         ind = graphp( gpt );                /* graphics file: real(sn1)      */
                                            /*                                */
         spt = spline( spt );              /* spline interpolation: imag(sn1) */
/*.......................................*/
         ll++;

         ii = null;
         j_ = null;
         while ( j_ < frqlbl-ONE )
         {
            if (( ll < MDV_SPCLPRTS  )
              &&( ll < idx[null] ))
               spl.vct[j_][ONE] = acpl[ll][j_];

            gph.vct[ii][ONE] = spl.fct[ii];
            ii++;
            k_ = ONE;
            while(  k_ <= i_ )
            {
               gph.vct[ii][ONE] = spl.fct[ii];
               ii++ ;
               k_++ ;
            };
            j_++ ;
         }; /* until j_ = frqlbl-ONE */

         if (( ll < MDV_SPCLPRTS  )
           &&( ll < idx[null] ))
            spl.vct[j_][ONE] = acpl[ll][j_];

         gph.vct[ii][ONE] = spl.fct[ii];
         ii++;

         strcpy( gph.file, "i" );
         strcat( gph.file, ptr );
         strcat( gph.file, "1_" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "imag(s" );
         strcat( gph.text, ptr );
         strcat( gph.text, "1)" );
         strcpy( gph.xunit, "Hz" );
         strcpy( gph.yunit, "--" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );                /* graphics file: imag(sn1)      */
/*.........................................*/

      };/* end while ( ll < SPECIAL_ ... ) */
      break;
/*.................... s-parameter computation finished ......................*/










/*............................................................................*/
     case 4:

      hh = ONE;
      while ( hh <= idx[null] )
      {
         h_ = hh - ONE;

         evl1 = ZERO;
         dlt1 = ZERO;

         absmaxm[h_] = ZERO;
         absmean[h_] = ZERO;
         meanenv[h_] = ZERO;
         timemax[h_] = ZERO;
            envl[h_] = null;

         i_ = null;
         while ( i_ < ( fpt->ttlg[null] ))
         {
            evl0 = ( fpt->r[hh][i_] );
            norm = fabs( evl0 );
            dlt0 = evl0 - evl1;

            absmean[h_] += norm;

            if ( absmaxm[h_] < norm )
            {
               absmaxm[h_] = norm;
               timemax[h_] = ( fpt->t[null] ) + i_*( fpt->dt[null] );
            };

            if ( dlt0*dlt1 < ZERO )/* rel.extremum: compute mean.envlp.*/
            {
               envl[h_]++;
               meanenv[h_] += fabs( evl1 );
            };

            if ( null < i_ )
               dlt1 = dlt0;

            evl1 = evl0;

            i_++ ;
         }; /* next i_ */

         absmean[h_] /= ( fpt->ttlg[null] );

         if ( null < envl[h_] )
         {
            meanenv[h_] /= envl[h_];
            modulation[h_] = 100.*(fabs( absmaxm[h_] - meanenv[h_] ));

            if ( 1.e-277 < fabs(meanenv[h_] ) )
               modulation[h_] /= meanenv[h_];
            else
               modulation[h_] = HUGE_VALF;
         }
         else
         {
            meanenv[h_] = ZERO;
            modulation[h_] = ZERO;
         };
         hh++ ;
      }; /* next hh */

/* tabl4: */

      printf( "\n\n Special parameters evaluated in file %s :",
         ( state->file ));
      printf( "\n\n index | abs. maximum : found at time|    abs. mean|  "
         "mean envlp.|   modulation" );
      printf( "\n       | [%10s] : [ %10s]|    [%7s]|  [%9s]| %12s",
         ( vpt->yunit ), ( vpt->xunit ), ( vpt->yunit ), ( vpt->yunit ),
         "[ percent]" );
      printf( "\n -----------------------------------"
         "-------------------------------------------" );

      hh = ONE;
      while ( hh <= idx[null] )
      {
         h_ = hh - ONE;
         printf( "\n %6d| %.7e: %.7e|%.7e|%.7e|%.7e", idx[hh], absmaxm[h_],
         timemax[h_], absmean[h_], meanenv[h_], modulation[h_] );
         hh++ ;
      };
      printf( "\n -----------------------------------"
         "-------------------------------------------" );

/* menu4: */

      kk = idx[null] - MDV_SPCLPRTS;

      i_ = MDV_GPHINTPL; /* number of interpolation points */

      ii = null;
      h_ = null;
      hh = ONE;
      while ( hh < kk )
      {
         mm = hh + MDV_SPCLPRTS;
         nn = hh + MDV_SPCLPRTS + ONE;

         spl.vct[h_][null] = idx[mm];
         spl.vct[h_][ONE]  = meanenv[mm-ONE];
         spl.dmn[ii] = spl.vct[h_][null];

         k_ = ONE;
         while ( k_ <= i_ )
         {
            x2 = ( double ) k_ / ( i_ + ONE );
            x1 = 1. - x2;
            spl.dmn[ii]  = x1*idx[mm] + x2*idx[nn];
            ii++ ;
            k_++ ;
         };
         h_ = hh;
         hh++ ;
      }; /* until hh = kk; h_ = kk - ONE */
      spl.vct[h_][null] = idx[idx[null]];
      spl.vct[h_][ONE]  = meanenv[idx[null]-ONE];
      spl.dmn[ii] = idx[idx[null]];
      ii++ ;

      printf( "\n %s: VSWR spline interpolation started.", __func__ );

      spl.mm = kk;
      spl.nn = ii;
/*............................................................................*/
      spt = spline( spt );      /* spline interpolation: VSWR                 */
/*............................*/

      printf( "\r %s: VSWR spline interpolation terminated.", __func__ );

      spl.fmax = fabs( spl.fmax );
      spl.fmin = fabs( spl.fmin );

      yy = spl.fmax/spl.fmin;               /* voltage standing wave ratio    */
      vswr[null] = yy;

      xx =  ( yy - 1. )/( yy + 1. );        /* reflexion factor               */
      refl[null] = xx;

      if ( fabs( xx ) < .75 )
         printf( "\n\n Source line matching:" );
      else
         printf( "\n\n Source line matching [ approximate values ]:" );

      printf( "\n\n VSWRatio  = % .12e ", yy );

      zz = 20.*log10( xx );                 /* reflexion [dB]                 */
      printf( "\n\n REFLEXION = % .12e dB [ %.12e percent ] ", zz, 100.*xx );

      yy = sqrt( 1. - xx*xx );              /* insertion factor               */
      zz = 20.*log10( yy );                 /* insertion loss [dB]            */
      printf( "\n INSERTION = % .12e dB [ %.12e percent ] ", zz, 100.*yy );

      printf( "\n\n S-parameters: \n" );

      zz = 20.*log10( xx );                 /* S11 = reflexion [dB]           */
      printf( "\n S[1,1] = % .12e dB ", zz );

      xx = .5*( spl.fmax + spl.fmin );      /* incident wave amplitude        */
      xx *= xx;                              

      hh = ONE;
      while ( hh < MDV_SPCLPRTS )
      {
         yy = meanenv[ hh ];
         zz = 10.*log10( yy*yy / xx );
         printf( "\n S[%ld,1] = % .12e dB ", hh+ONE, zz );
         hh++ ;
      };

      printf( "\n\n Please acknowledge [ enter any character ]: " );
      scanf( "%s", ptr );

      strcpy( gph.file, "stw" );
      strcat( gph.file, ( state->flbl ));
      strcpy( gph.name, ( state->name ));
      strcpy( gph.text, "standing_wave" );
      strcpy( gph.xunit, "index" );
      strcpy( gph.yunit, "volts" );

      ii = null;
      hh = ONE;
      while ( hh < kk )
      {
         gph.vct[ii][null] = spl.dmn[ii];
         gph.vct[ii][ONE]  = spl.fct[ii];

         k_ = ONE;
         while ( k_ <= i_ )
         {
            gph.vct[ii][null] = spl.dmn[ii];
            gph.vct[ii][ONE]  = spl.fct[ii];
            ii++ ;
            k_++ ;
         };
         hh++ ;
      }; /* next hh */
      gph.vct[ii][null] = spl.dmn[ii];
      gph.vct[ii][ONE]  = spl.fct[ii];
      ii++;

      gph.nn = ii;
/*............................................................................*/
      ind = graphp( gpt );      /* graphics file: standing wave               */
/*............................*/

      break;

     case 5:

      if ( null == strncmp( vpt->text, "reference", NINE ))
      {
         strcpy( spcfle, spcpfx );
         strcat( spcfle, "ref" );
         strcat( spcfle, state->flbl );

         evalfle = fopen( spcfle, "w+" );

         if ( evalfle == null )
         {
            printf( "\n Error on opening reference file "
               "%s !", spcfle );
            printf( "\n [overriding: values not saved.]\n " );
            break;
         };

         fprintf( evalfle, "%s\n", ( vpt->name ));
         fprintf( evalfle, "%s\n", "reference_spectrum" );
         fprintf( evalfle, "%s\n", ( vpt->xunit ));
         fprintf( evalfle, "%s\n", ( vpt->yunit ));
	 
         fprintf( evalfle, d_format, ( double )( fpt->t[null] ), "\n" );
         fprintf( evalfle, d_format, ( double )( fpt->tt[null] ), "\n" );
         fprintf( evalfle, d_format, ( double )( fpt->dt[null] ), "\n" );
         fprintf( evalfle, "%ld \n", ( fpt->ttlg[null] ));
/*............................................................................*/
# if EVL_WRTFREQ == 1
	 frequency = ( fpt->t[null] );
# endif
/*............................................................................*/
         jj = null;
         while ( jj < ( fpt->ttlg[null] ))
         {
/*............................................................................*/
# if EVL_WRTFREQ == 1
            fprintf( evalfle, d_format, frequency, "  " );
	    frequency += ( fpt->dt[null] );
# endif
/*............................................................................*/
            fprintf( evalfle, d_format, ( fpt->r[ONE][jj] ), "  " );
            fprintf( evalfle, d_format, ( fpt->i[ONE][jj] ), "\n" );
            jj++ ;
         };

         nseconds = time( timer );
         strcpy( tmeptr, ctime( &nseconds ));

         fprintf( evalfle, "\n# DSC reference %s\n", spcfle );
         fprintf( evalfle, "# created:%s", tmeptr );

         fclose( evalfle );

	 printf( "\n" );
         printf( CLEAR_LINE );
         printf( "\n Reference file %s ", spcfle );
         printf( "created:\n %.24s\n ", tmeptr );

      } /* end if ( null == strncmp( vpt->text, "reference", NINE )) */
      else if ( null < idx[null]  )
      {
/* enter external reference <spcpfx>.ref: */

         strcpy( spcfle, spcpfx );
         strcat( spcfle, "ref" );

         printf( "\n" );

        open_reference2:

         evalfle = fopen( spcfle, "r+" );

         if ( evalfle == null )
         {
            printf( "\n Reference file %s not found "
               "in present directory:", spcfle );
            printf( "\n Please re-enter filename [ Escape: "
               "enter null ] >----> " );
            scanf( "%s", spcfle );

            if ( *spcfle == '0' )
               break;
            else
               goto open_reference2;
         }
         else if ( evalfle != null )
         {
            printf( "\n opened: reference file %s ", spcfle );

            fscanf( evalfle, "%s", ptr );

            if ( null != strncmp( ptr, vpt->name, THREE ))
            {
               printf( "\n\n Error message from function %s :", __func__ );
               printf( "\n\n Incompatible system identifier '%s'", ptr );
               printf( "\n on reference spectrum, file %s !!!", spcfle );
               printf( " [ overriding. ]\n" );

               fclose( evalfle );

               break;
            };

            fscanf( evalfle, "%s", ptr );
            fscanf( evalfle, "%s", absc_unit );
            fscanf( evalfle, "%s", ordn_unit );
            fscanf( evalfle, "%s", ptr );/* xlower */
            xlower = strtod( ptr, endp );

            fscanf( evalfle, "%s", ptr );/* xupper */
            xupper = strtod( ptr, endp );

            fscanf( evalfle, "%s", ptr ); /* dx */
            dx = strtod( ptr, endp );

            fscanf( evalfle, "%s", ptr );
            ii = strtol( ptr, endp, DEC );

            xx = xlower;
            jj = null;
            while( jj < ii )
            {
/*............................................................................*/
# if EVL_WRTFREQ == 1
               fscanf( evalfle, "%s", ptr );
# endif
/*............................................................................*/
               fscanf( evalfle, "%s", ptr );
               rref[jj] = strtod( ptr, endp );
               fscanf( evalfle, "%s", ptr );
               iref[jj] = strtod( ptr, endp );
               xx += dx;

               jj++ ;
            };
            xupper = xx - dx;
            printf( "\r entered: reference %s .      \n",
               spcfle );
            fclose( evalfle );
            opm1 = ONE;
         };

         jj = null; 
         while( jj < ii )
         {
            real = rref[jj];
            imag = iref[jj];
            norm = real*real + imag*imag;

            ( fpt->r[ONE][jj] ) -= real;
            ( fpt->i[ONE][jj] ) -= imag;

            rspc[null][jj] = ( fpt->r[ONE][jj] );
            ispc[null][jj] = ( fpt->i[ONE][jj] );

            if ( 1.e-277 < norm )
            {
               xx = ( real*rspc[null][jj] + imag*ispc[null][jj] ) / norm;
               yy = ( real*ispc[null][jj] - imag*rspc[null][jj] ) / norm;
               zz = sqrt( xx*xx + yy*yy );

               if ( 1.e-277 < fabs( 1. - zz ) )
               {
                  vswr[jj] = ( 1. + zz )/( 1. - zz );
                  refl[jj] = 100.*zz;
               } 
               else
               {
                  vswr[jj] = HUGE_VALF;
                  refl[jj] = 100.;
               };
            }
            else
            {
               xx = HUGE_VALF;
               yy = HUGE_VALF;
               zz = HUGE_VALF;
               vswr[jj] = HUGE_VALF;
               refl[jj] = HUGE_VALF;
            };

            kk = MDV_SPCLPRTS;

            if ( idx[null] < kk )
               kk = idx[null]; 

            hh = ONE;
            while( hh <= kk )
            {
               h_ = hh - ONE;

               real = rref[jj];
               imag = iref[jj];
               norm = real*real + imag*imag;

               rspc[h_][jj] = ( fpt->r[hh][jj] );
               ispc[h_][jj] = ( fpt->i[hh][jj] );

               if ( 1.e-277 < norm )
               {
                  xx = ( real*rspc[h_][jj] + imag*ispc[h_][jj] ) / norm;
                  yy = ( real*ispc[h_][jj] - imag*rspc[h_][jj] ) / norm;
                  zz = xx*xx + yy*yy;

                  if ( 1.e-277 < zz )
                  {
/*............................................................................*/
# if MDV_NORMALZE == 1

                     zz = sqrt( zz );
                     xx /= zz;
                     yy /= zz;
                  
                     acpl[h_][jj] = 0.;
# else
                     acpl[h_][jj] = 10.*log10( zz );
# endif
/*............................................................................*/
                  }
                  else
                     acpl[h_][jj] = - LARGE_LOG_VAL;
               }
               else
               {
                  xx = HUGE_VALF;
                  yy = HUGE_VALF;
                  zz = HUGE_VALF; 

                  acpl[h_][jj] = LARGE_LOG_VAL;
               };

               rcpl[h_][jj] = xx;
               icpl[h_][jj] = yy;

               hh++ ;
            };
            jj++ ;

         };/* next jj */

         jj = null; xx = xlower; 
         while( jj < ii )
         {
            gph.vct[jj][ONE] = vswr[jj];
            gph.vct[jj][null] = xx;

            xx += dx;
            jj++ ;
         };/* next jj */

         strcpy( gph.file, "swr" );
         strcat( gph.file, state->flbl );
         strcpy( gph.name, state->name );
         strcpy( gph.text, "voltage_standing_wave_ratio" );
         strcpy( gph.xunit, absc_unit );
         strcpy( gph.yunit, "---" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );        /* graphics file: VSWR                   */
/*.................................*/
         jj = null; xx = xlower;
         while( jj < ii )
         {
            gph.vct[jj][ONE] = refl[jj];
            gph.vct[jj][null] = xx;

            xx += dx ;
            jj++ ;
         };/* next jj */

         strcpy( gph.file, "ref" );
         strcat( gph.file, state->flbl );
         strcpy( gph.name, state->name );
         strcpy( gph.text, "reflexion" );
         strcpy( gph.xunit, absc_unit );
         strcpy( gph.yunit, "0/0" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );        /* graphics file: reflexion [0/0]        */
/*.................................*/
         jj = null; xx = xlower;
         while( jj < ii )
         {
            zz = refl[jj];

            if( 1.e-277 < zz  )
               gph.vct[jj][ONE] = 20.*log10( zz/100. );
            else                                  /* Replacing LARGE_LOG_VAL  */
               gph.vct[jj][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                                  /* made graphp(*) function  */
                                                  /* inoperative              */
            gph.vct[jj][null] = xx;

            xx += dx;
            jj++ ;
         };/* next jj */

         strcpy( gph.file, "rtl" );
         strcat( gph.file, ( state->flbl ));
         strcpy( gph.text, "return_loss" );
         strcpy( gph.xunit, absc_unit );
         strcpy( gph.yunit, "dB" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );        /* graphics file: return loss [dB]       */
/*.................................*/
         jj = null; xx = xlower;
         while( jj < ii )
         {
            zz = refl[jj]/100.;
       
            if( zz < 1. )
               gph.vct[jj][ONE] = 10.*log10( 1. - zz*zz );
            else                                  /* Replacing LARGE_LOG_VAL  */
               gph.vct[jj][ONE] = -LARGE_LOG_VAL; /* by HUGE_VALF at times */
                                                  /* made graphp(*) function  */
                                                  /* inoperative              */
            gph.vct[jj][null] = xx;

            xx += dx;
            jj++ ;
         };/* next jj */

         strcpy( gph.file, "isl" );
         strcat( gph.file, state->flbl );
         strcpy( gph.text, "insertion_loss" );
         strcpy( gph.xunit, absc_unit );
         strcpy( gph.yunit, "dB" );
         gph.nn = ii;
/*............................................................................*/
         ind = graphp( gpt );        /* graphics file: insert.loss [dB]       */
/*.................................*/
/* still as above:
         kk = MDV_SPCLPRTS;
	 
         if ( idx[null] < kk )
            kk = idx[null]; 
*/
         hh = ONE;
         while( hh <= kk )
         {
            h_ = hh - ONE;
            jj = null; xx = xlower;
            while( jj < ii )
            {
               gph.vct[jj][ONE] = acpl[h_][jj];
               gph.vct[jj][null] = xx;

               xx += dx;
               jj++ ;
            };/* next jj */
            
            strcpy( ptr, "s" );
            strcat( ptr, lotos( hh, null ) );
            strcat( ptr , "1_" );
            strcat( ptr, state->flbl );
            strcpy( gph.file, ptr );
            strcpy( gph.text, ptr );
/*
            strcat( gph.text, "[steady_state]" );
*/
            strcat( gph.text, "_" );
            strcat( gph.text, vpt->text );
            strcpy( gph.xunit, absc_unit );
            strcpy( gph.yunit, "dB" );
            gph.nn = ii;
/*............................................................................*/
            ind = graphp( gpt );           /* graphics file: s[k,1] [dB]      */
/*.......................................*/
/* complex s-parameters : */

            strcpy( spcfle, spcpfx );
            strcat( spcfle, ptr );

            evalfle = fopen( spcfle, "w+" );

            if ( evalfle == null )
            {
               printf( "\n Error on opening file "
                  "%s !", spcfle );
               printf( "\n [overriding: values not saved.]\n " );
            }
            else
            {
               fprintf( evalfle, "%s\n", ( vpt->name ));
               fprintf( evalfle, "%s_%s\n", ptr, ( vpt->text ));
               fprintf( evalfle, "%s\n", ( vpt->xunit ));
               fprintf( evalfle, "%s\n", "---" );
               fprintf( evalfle, d_format, xlower, "\n" );
               fprintf( evalfle, d_format, xupper, "\n" );
               fprintf( evalfle, d_format, dx, "\n" );
               fprintf( evalfle, "%ld\n", ii );
/*............................................................................*/
# if EVL_WRTFREQ == 1
	       frequency = xlower;
# endif
/*............................................................................*/
               jj = null;
               while ( jj < ii )
               {
/*............................................................................*/
# if EVL_WRTFREQ == 1
                  fprintf( evalfle, d_format, frequency, "  " );
	          frequency += dx;
# endif
/*............................................................................*/
                  fprintf( evalfle, d_format, rcpl[h_][jj], "  " );
                  fprintf( evalfle, d_format, icpl[h_][jj], "\n" );
                  jj++ ;
               };

               nseconds = time( timer );
               strcpy( tmeptr, ctime( &nseconds ));

               fprintf( evalfle, "\n# DSC file %s\n", spcfle );
               fprintf( evalfle, "# created:%s", tmeptr );

               fclose( evalfle );

            }; /* end if evalfle != null */  

            hh++ ;
         }; /* end while ( hh <= kk [ kk <= MDV_SPCLPRTS ) ] */
      }; /* end if ( null < idx[null]  ) */

      break;
/*............................................................................*/

     default:
      break;
   }; /* end switch ( lbl ) */

   PRBLDCLR( "\r" );
   printf( " %*s", 78, "MODVAL" );
   PRNORMAL( "" );
   ( state->rtn ) = null;

   return state;
}
/*============================================================================*/
# undef LARGE_LOG_VAL
/************************ end of function modval(*) ***************************/


/*----------------------------------------------------------------------------*/
# endif /* [ end of section compiled with option -D_POST_MODEL ] */
/*----------------------------------------------------------------------------*/
# undef EXC_MAXW
# undef EXC_HCRR
/*************************** end of file model.c ******************************/
