/* [ function: model.c ] */
# define DSC_MODEL "mod_ccoax.G"
/*******************************************************************************
*                                                                              *
*   DSC model generation and evaluation functions model.c and modval.c         *
*   Prototype: 'convection in coaxial line'                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 27, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define MOD_DEFLT 7 /* default models */
/*----------------------------------------------------------------------------*/
# define MOD_FLOWS 1 /* evaluate flows - none:0, only vertical [y]:1, all:2 */
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
/*____________________________________________________________________________*/



/*____________________________________________________________________________*/
# ifdef _PRE_MODEL
/* THIS SECTION IS COMPILED ONLY WITH OPTION -D_PRE_MODEL */
/*******************************************************************************
*                                                                              *
*   DSC model generation function model.c                                      *
*   Prototype: 'convection in coaxial line'                                    *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
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
/* Edit and customize this header for FORMER.C configuration: */
# include "./former/FORMER.CONF" 
/*----------------------------------------------------------------------------*/
/* A macro that gives a conveniently readable form to the string returned */
/* with [ the pointer of ] function ctime(&nseconds) */
# include "./tools/TIMEFORM.M"
/*----------------------------------------------------------------------------*/
# ifndef DSC_ADJCLL
   # define DSC_ADJCLL 1 
# endif
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
   # include <ncurses.h>
# endif 
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
/*******************************************************************************
*                                                                              *
*   TRANSFER STRUCTURE for DSC model mod_ccoax                                 *
*                                                                              *
*   This structure defines any set of model specific parameters [ shared by    *
*   the DSC model generating functions systop(*),...,sysval(*) ]. - It is,     *
*   thus, in general highly model dependent (!), but can likewise be nearly    *
*   identical for different DSC models.                                        *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 20, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
The following structure 'trf' of type 'transfer' may be arbitrarily modified as
suitable, by introducing model dependent parameters [ e.g. grid characterist -
ics, lengths, material constants etc.]. This structure is, hence, especially
designed to define [ any shared ] parameters that have to be exchanged between
the DSC model defining functions systop(*) syssmx(*), sysbnd(*), sysexc(*), and
sysval(*).
*/
/*                      | maximum number M,                                   */
/*                      V                                                     */
# define MOD_BLOCKS    10 /* maximum M, BLOCK(0),...,BLOCK(M)                 */
# define MOD_DOMAINS    5 /* maximum M, blc[1],...,blc[M]                     */
# define MOD_LAYERS   100 /* layers [in z-direction]                          */
# define MOD_OPERATS   10 /* operations etc., trf.c[0],...,trf.c[M]           */
# define MOD_DIVISNS   20 /* divisions etc., trf.n[0],...,trf.n[M]            */
# define MOD_PARMTRS   40 /* parameters trf.s[0],...,trf.s[M]                 */
# define MOD_BOUNDRS 1000 /* maximum M, trf.bdn[1],...,trf.bdn[M-1]           */
# define MOD_EXCITES  100 /* maximum M, trf.exn[1],...,trf.exn[M-1]           */
# define MOD_EVLUATE   50 /* maximum M, trf.val[0],...,trf.val[M-1]           */
# define MOD_ARCS      10 /* maximum M, trf.alfa[0],...,trf.alfa[M]           */
# define MOD_POINTS    20 /* maximum M, pnt[0],...,pnt[M-1]                   */
/*............................................................................*/
# if DSC_HCRMDE != 0
   # define MOD_HCBNDRS 1000
   # define MOD_HCEXCTS  100
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
*   functions, on header file "puzzab.h", which are included below.            *
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

   cpylne( trf.ctx[1], "reference_line/_waveguide", "0/1", ll );
   cpylne( trf.ctx[2], "time/_frequency_domain", "0/1" ,ll );
   cpylne( trf.ctx[3], "disable/_enable_internal_convection", "0,1", ll );
   cpylne( trf.ctx[4], "apply_skin_effect_heat_sources_to_metal/_gas",
      "0,1", ll );
   cpylne( trf.ctx[5], "apply_sources_to_inner/_outer/_all_conductors",
      "0/1/2", ll );
   cpylne( trf.ctx[6], "no_slip/free_slip_boundary_cds_at_conductors",
      "0/1", ll );
   cpylne( trf.ctx[7], "no_slip/free_slip/_outflow_bd_conds,_axial",
      "0/1/2", ll );

   trf.c[0] =  7; /* number M of operation modes */

   trf.c[1] = 1; /* reference line / structure [0/1] */
   trf.c[2] = 1; /* time / frequency domain [0/1] */
   trf.c[3] = 1; /* disable / enable internal convection [0/1] */

/*.............................. defaults ....................................*/
# if MOD_DEFLT == 1  
   trf.c[3] = 1;
   trf.c[4] = 0; /* apply skin effect heat sources to metal */
   trf.c[5] = 2; /* skin effect losses on all conductors */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 2
   trf.c[3] = 1;
   trf.c[4] = 1; /* apply skin effect heat sources to gas */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 3
   trf.c[3] = 0;
   trf.c[4] = 0; /* apply skin effect heat sources to metal */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 4
   trf.c[3] = 1;
   trf.c[4] = 1; /* apply skin effect heat sources to gas */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 5
   trf.c[3] = 1;
   trf.c[4] = 1; /* apply skin effect heat sources to gas */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 6
   trf.c[3] = 1;
   trf.c[4] = 1; /* apply skin effect heat sources to gas */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 7
   trf.c[3] = 1;
   trf.c[4] = 1; /* apply skin effect heat sources to gas */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 8
   trf.c[3] = 1;
   trf.c[4] = 1; /* apply skin effect heat sources to gas */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 9
   trf.c[3] = 1;
   trf.c[4] = 1; /* apply skin effect heat sources to gas */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# elif MOD_DEFLT == 10
   trf.c[3] = 1;
   trf.c[4] = 1; /* apply skin effect heat sources to gas */
   trf.c[5] = 0; /* only heat sources on inner conductor */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# else /* [ cf. MOD_DEFLT == 1 ] */
   trf.c[3] = 1;
   trf.c[4] = 0; /* apply skin effect heat sources to metal */
   trf.c[5] = 2; /* skin effect losses on all conductors */
   trf.c[6] = 0; /* no slip boundary conditions at conductors */
   trf.c[7] = 1; /* free slip boundary conditions, axial */
# endif
/*............................................................................*/

   return;
}
/*=================== end of function deflt_operts(*) ========================*/

short rvise_operts( void )
{ 
/*
   extern FORMSTATE *spt;

   extern struct transfer trf;
   extern struct blcstruc blc;
*/

/* operation marks: */

   if( null == trf.c[ONE] )
   {
      trf.ref = ONE;
      strcpy(( spt->tpt->text ), "reference_line" );
   }
   else
   {
      trf.c[ONE] = ONE;
      trf.ref = null;
      strcpy(( spt->tpt->text ), "coaxial_line__" );
   };

   if ( trf.c[TWO] != null )
      trf.c[TWO] = ONE;

   if ( trf.c[THREE] != null )
      trf.c[THREE] = ONE;

   if ( trf.c[FOUR] != null )
      trf.c[FOUR] = ONE;

   if (( trf.c[FIVE] != null )
     &&( trf.c[FIVE] != ONE ))
      trf.c[FIVE] = TWO;

   if ( trf.c[SIX] != null )
      trf.c[SIX] = ONE;

   return null;
}
/*=================== end of function rvise_operts(*) ========================*/

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
   cpylne( trf.mtx[2], "vertical_divisions,_domain__1", "", ll );
   cpylne( trf.mtx[3], "vertical_divisions,_domain__2", "", ll );
   cpylne( trf.mtx[4], "vertical_divisions,_domain__3", "", ll );
   cpylne( trf.mtx[5], "vertical_divisions,_domain__4", "", ll );
   cpylne( trf.mtx[6], "vertical_divisions,_domain__5", "", ll );
   cpylne( trf.mtx[7], "vertical_divisions,_domain__6", "", ll );
   cpylne( trf.mtx[8], "vertical_divisions,_domain__7", "=m6;_dependent", ll );
   cpylne( trf.mtx[9], "vertical_divisions,_domain__8", "=m5;_dependent", ll );
   cpylne( trf.mtx[10], "vertical_divisions,_domain__9", "=m4;_dependent", ll );
   cpylne( trf.mtx[11], "vertical_divisions,_domain_10", "=m3;_dependent", ll );
   cpylne( trf.mtx[12], "vertical_divisions,_domain_11", "=m2;_dependent", ll );

/* z-divisions. */
   blc.m[null] = 4;
/*.............................. defaults ....................................*/
# if MOD_DEFLT == 1
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  10; /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 2
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  10; /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 3
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  10; /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 4
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  10; /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 5
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  3;  /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 6
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  1;  /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 7
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  1;  /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 8
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  10;  /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 9
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  1;  /* .......*/
   blc.m[4] =  10; /* termination line */
# elif MOD_DEFLT == 10
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  10;  /* .......*/
   blc.m[4] =  10; /* termination line */
# else /* [ cf. MOD_DEFLT == 1 ] */
   blc.m[1] =  20; /* number of divisions, domain 0 [ source line ] */
   blc.m[2] =  10; /* "      "  divisions, domain 1 [ buffer ] */
   blc.m[3] =  10; /* .......*/
   blc.m[4] =  10; /* termination line */
# endif
/*............................................................................*/
/* xy_divisions, identifiers etc.:                                            */

/* write only CONNECTED (!) STRINGS in the string copy function strcpy(*) !!! */
/* [ 2nd. argument       |  ]                                                 */
/*                       V                                                    */
   strcpy( trf.ntx[0] , "divisions,_identifiers,_operation_marks,_etc." );

   cpylne( trf.ntx[1], "ab_divisions,_blocks__1,5", "", ll );
   cpylne( trf.ntx[2], "bc_divisions,_block___1", "", ll );
   cpylne( trf.ntx[3], "ab_divisions,_block___3", "", ll );
   cpylne( trf.ntx[4], "bc_divisions,_block___6", "", ll );
   cpylne( trf.ntx[5], "sum,_inner_divisions", "dependent", ll );
   cpylne( trf.ntx[6], "sum,_outer__divisions", "dependent", ll );

   trf.n[0] = 6; /* number M of parameters trf.n[1],...,trf.n[M] */

/* xy_divisions: */
/*.............................. defaults ....................................*/
# if MOD_DEFLT == 1
   trf.n[1]  =  4;
   trf.n[2]  =  8;
   trf.n[3]  =  8;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 2
   trf.n[1]  =  4;
   trf.n[2]  =  8;
   trf.n[3]  =  8;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 3
   trf.n[1]  =  4;
   trf.n[2]  =  8;
   trf.n[3]  =  8;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 4
   trf.n[1]  =  6;
   trf.n[2]  = 10;
   trf.n[3]  = 14;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 5
   trf.n[1]  =  6;
   trf.n[2]  = 10;
   trf.n[3]  = 14;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 6
   trf.n[1]  =  6;
   trf.n[2]  = 10;
   trf.n[3]  = 14;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 7
   trf.n[1]  =  6;
   trf.n[2]  = 12;
   trf.n[3]  = 14;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 8
   trf.n[1]  =  6;
   trf.n[2]  = 12;
   trf.n[3]  = 14;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 9
   trf.n[1]  =  6;
   trf.n[2]  = 12;
   trf.n[3]  = 14;
   trf.n[4]  =  1;
# elif MOD_DEFLT == 10
   trf.n[1]  =  6;
   trf.n[2]  = 12;
   trf.n[3]  = 14;
   trf.n[4]  =  1;
# else /* [ cf. MOD_DEFLT == 1 ] */
   trf.n[1]  =  4;
   trf.n[2]  =  8;
   trf.n[3]  =  8;
   trf.n[4]  =  1;
# endif
/*............................................................................*/

   return;
}
/*=================== end of function deflt_divsns(*) ========================*/

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
      blc.m[0] = 2;
   else
      blc.m[0] = 4;

/* dependencies: */

   blc.m[4] = blc.m[2];
   trf.n[5] = 2*trf.n[1] + trf.n[3];
   trf.n[6] = 4*trf.n[2] + trf.n[5];

   return ONE;
}
/*=================== end of function rvise_divsns(*) ========================*/
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

   cpylne( trf.stx[1], "frequency", "Hz", ll );
   cpylne( trf.stx[2], "outer_diameter,_coaxial_line", "", ll );
   cpylne( trf.stx[3], "inner_diameter,_coaxial_line", "", ll );
   cpylne( trf.stx[4], "wall_thickness,_mantle", "", ll );
   cpylne( trf.stx[5], "wall_thickness,_inner_tube", "", ll );
   cpylne( trf.stx[6], "length,_source_line", "", ll );
   cpylne( trf.stx[7], "lengths,_buffer_and_termination", "", ll );
   cpylne( trf.stx[8], "length,_cooled_section", "", ll );
   cpylne( trf.stx[9], "electric_conductivity,_waveguide", "S/m", ll );
   cpylne( trf.stx[10], "equivalent_square_resistance", \
      "Ohms;_dependent" , ll );
   cpylne( trf.stx[11], "penetration_depth", "m;_dependent", ll );
   cpylne( trf.stx[12], "heat_conductivity,_mantle", "W/(K*m)", ll );
   cpylne( trf.stx[13], "heat_capacity,_mantle", "J/(K*m^3)", ll );

# if ELDE_HCRMDE != 0   
   if ( TEMPGGE == 0 )
      cpylne( trf.stx[14], "cooling_temperature", "K", ll );
   else
      cpylne( trf.stx[14], "cooling_temperature", "deg_C", ll );
# else
   cpylne( trf.stx[14], "cooling_temperature",
      "inactive", ll );
# endif

   cpylne( trf.stx[15], "environment-surface_heat_conductivity",
      "W/(K*m^2)", ll );

# if ELDE_HCRMDE != 0   
   if ( TEMPGGE == 0 )
      cpylne( trf.stx[16], "environment_temperature", "K", ll );
   else
      cpylne( trf.stx[16], "environment_temperature", "deg_C", ll );
# else
   cpylne( trf.stx[16], "environment_temperature",
      "inactive", ll );
# endif

   cpylne( trf.stx[17], "equiv._CW_power_input", "Watts", ll );
   cpylne( trf.stx[18], "excitation_gauge", "change_with_care_!", ll );
   cpylne( trf.stx[19], "cooled_outer_surface",
      "m^2;_dependent", ll );

# if ELDE_HCRMDE != 0   
   if ( TEMPGGE == 0 )
      cpylne( trf.stx[20], "initial_mean_temperature,_gas","K", ll );
   else
      cpylne( trf.stx[20], "initial_mean_temperature,_gas","deg_C", ll );
# else
   cpylne( trf.stx[20], "initial_mean_temperature,_gas",
      "inactive", ll );
# endif

   cpylne( trf.stx[21], "mean_density,_gas","Kg/m^3", ll );
   cpylne( trf.stx[22], "thermal_expansion_coefficient,_gas","1/K", ll );
   cpylne( trf.stx[23], "adiabatic_compression_coefficient,_gas","1/Pa", ll );
   cpylne( trf.stx[24], "heat_conductivity,_gas","W/(m*K)", ll );
   cpylne( trf.stx[25], "heat_capacity,_gas","J/(Kg*K)", ll );
   cpylne( trf.stx[26], "dynamic_viscosity,_gas","Kg/(m*s)", ll );
   cpylne( trf.stx[27], "gravitational_acceleration","m/sec^2;_fixed", ll );
   cpylne( trf.stx[28], "dissipation_time","sec", ll );
   cpylne( trf.stx[29], "characteristic_length","m", ll );
   cpylne( trf.stx[30], "start_fluid_dynamic_computations_with_that_delay",
      "sec", ll );
   cpylne( trf.stx[31], "coarsening_period_[LES_filter]", "sec", ll );
/*............................................................................*/
# if DSC_HCRMDE != 0
   cpylne( trf.stx[32], "PROPOSED_heat&fluid_time_step_[s]",
      "later_requested", ll );
# else
   cpylne( trf.stx[32], "PROPOSED_heat&fluid_time_step_[s]",
      "inactive", ll );
# endif
/*............................................................................*/
   cpylne( trf.stx[33], "Nusselt_number_at_no-slip_boundary_faces", "", ll );

/* default values:                                                           */
/*............................................................................*/
   trf.s[0]  = 32;               /* parameters [ international units ] */

/*.............................. defaults ....................................*/
# if DSC_HCRMDE != 0
   if ( TEMPGGE == 0 )
   {
      trf.s[14] = 3.13000e+02; /* cooling temperature [Kelvin] */
      trf.s[16] = 2.93000e+02; /* environment temperature [Kelvin] */
      trf.s[20] = 2.73000e+02; /* initial mean temperature, gas [Kelvin] */
   }
   else
   {
      trf.s[14] = 4.00000e+01; /* cooling temperature [Celsius] */
      trf.s[16] = 2.00000e+01; /* environment temperature [Celsius] */
      trf.s[20] = 0.00000e+00; /* initial mean temperature, gas [Celsius] */
   };
# endif
/*............................................................................*/
# if MOD_DEFLT == 1
/* [ stable with heat and fluid time step 5.00e-03 s */
/*   and LES coarsening period 3.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-01; /* length, cooled section */
   trf.s[9] = 4.18000e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 2.00000e+05; /* effective CW power input */
   trf.s[18] = 1.46190e+02; /* excitation gauge */

/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.43000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 1.00000e+02; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 3.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 5.00000e-03; /* [proposed] heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 2
/* [ stable with heat and fluid time step 5.00e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-01; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 1.60620e+05; /* effective CW power input */
   trf.s[18] = 1.59959e+02; /* excitation gauge */

/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 1.50000e+01; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 5.00000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 3
/* [ stable with heat and fluid time step 5.00e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-01; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 1.60620e+05; /* effective CW power input */
   trf.s[18] = 1.59959e+02; /* excitation gauge */
   
/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 5.00000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 4
/* [ stable with heat and fluid time step 5.00e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-01; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 1.60620e+05; /* effective CW power input */
   trf.s[18] = 1.59959e+02; /* excitation gauge */
   
/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 5.00000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 5
/* [ stable with heat and fluid time step 2.50e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 6.00000e-02; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 1.60620e+05; /* effective CW power input */
   trf.s[18] = 1.25286e+02; /* excitation gauge */
   
/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 2.50000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 6
/* [ stable with heat and fluid time step 2.50e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-02; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 1.60620e+05; /* effective CW power input */
   trf.s[18] = 1.40527e+02; /* excitation gauge */
   
/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 2.50000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 7
/* [ stable with heat and fluid time step 2.50e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-02; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 1.60620e+05; /* effective CW power input */
   trf.s[18] = 1.25286e+02; /* excitation gauge */
   
/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 2.50000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 8
/* [ stable with heat and fluid time step 2.00e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 1.00000e-01; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 1.60620e+05; /* effective CW power input */
   trf.s[18] = 1.25286e+02; /* excitation gauge */
   
/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 2.00000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 9
/* [ stable with heat and fluid time step 2.50e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 5.50000e+07; /* frequency [ Hz ] */
   trf.s[2] = 3.45000e-01; /* outer diameter */
   trf.s[3] = 2.09000e-01; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 6.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-02; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 5.00000e+06; /* effective CW power input */
   trf.s[18] = 2.50000e+01; /* excitation gauge */
   
/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 2.00000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# elif MOD_DEFLT == 10
/* [ stable with heat and fluid time step 2.50e-03 s */
/*   and LES coarsening period 2.00e-02 s ] */

   trf.s[1] = 5.50000e+07; /* frequency [ Hz ] */
   trf.s[2] = 3.45000e-01; /* outer diameter */
   trf.s[3] = 2.09000e-01; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 6.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-01; /* length, cooled section */
   trf.s[9] = 2.86500e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 5.00000e+06; /* effective CW power input */
   trf.s[18] = 2.50000e+01; /* excitation gauge */
   
/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.99000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 0.00000e+00; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 2.00000e-03; /* [proposed ]heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */

# else /* [ cf. MOD_DEFLT == 1 ] */
/* [ stable with heat and fluid time step 5.00e-03 s */
/*   and LES coarsening period 3.00e-02 s ] */

   trf.s[1] = 1.00000e+08; /* frequency [ Hz ] */
   trf.s[2] = 2.30000e-01; /* outer diameter */
   trf.s[3] = 9.99000e-02; /* inner diameter */
   trf.s[4] = 5.00000e-03; /* wall thickness */
   trf.s[5] = 2.00000e-03; /* wall thickness, inner tube */
   
   trf.s[6] = 4.00*trf.s[2]; /* length, source line */
   trf.s[7] = 0.50*trf.s[2]; /* length, buffer and termination */
   trf.s[8] = 2.00000e-01; /* length, cooled section */
   trf.s[9] = 4.18000e+07; /* electric conductivity, waveguide */
                           /* [ value for Cu, 125 deg_C: 4.18e+07 S/m ] */
   trf.s[12] = 3.93000e-02; /* heat conductivity, mantle */
                            /* value for Cu, 125 deg_C: 393 W/(m*K) */
   trf.s[13] = 3.51000e+02; /* heat capacity, mantle */
                            /* value for Cu, 125 deg_C: 3.51e+06 J/(Kg*K) */
                            /* [ the last two values are drastically */
                            /*   scaled in order to avoid stiffness ] */

   trf.s[15] = 4.06900e+00; /* environment surface conductance [W/(K*m^2)]*/
   trf.s[17] = 3.00000e+05; /* effective CW power input */
   trf.s[18] = 1.46190e+02; /* excitation gauge */

/* trf.s[19], cooled outer surface [ dependent ] */

   trf.s[21] = 9.99800e-01; /* mean density, gas [Kg/m^3]*/
   trf.s[22] = 2.83800e-03; /* thermal expans. coeff, gas [1/K]*/
   trf.s[23] = 1.00000e-05; /* adiabatic compression coeff, gas [1/K]*/
   trf.s[24] = 2.43000e-02; /* heat conductivity, gas [W/(m*K)] */
   trf.s[25] = 1.00720e+03; /* heat capacity, gas [J/(Kg*K)] */
   trf.s[26] = 2.10000e-05; /* dynamic viscosity, gas [Kg/(m*sec)]*/
   trf.s[27] = 9.81000e+00; /* gravitational acceleration [m/s^2]*/
   trf.s[28] = 1.00000e+02; /* dissipation time [sec] */
   trf.s[29] = 2.30000e-01; /* characteristic length [m] */
   trf.s[30] = 0.00000e+00; /* start fluid operations with that delay */
   trf.s[31] = 2.00000e-02; /* coarsening time [s] - may work with that */
   trf.s[32] = 5.00000e-03; /* [proposed] heat and fluid time step [s] */
   trf.s[33] = 1.00000e+00; /* Nusselt number at no-slip boundary */
# endif
/*............................................................................*/

   return;
}
/*=================== end of function deflt_params(*) ========================*/

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

   if( null == trf.c[TWO] )
      strcpy(( spt->ppt->domain ), "time_domain" );
   else
      strcpy(( spt->ppt->domain ), "frequency_domain" );

   trf.fr = trf.s[1];

   ( spt->ppt->fr ) = trf.fr;

   trf.ra = .5*trf.s[2];
   trf.ri = .5*trf.s[3];

   trf.cwp = trf.s[17];  /* effective CW input power */

/* excitation amplitude */ 

   trf.uu = sqrt( trf.s[18]*trf.cwp );
   trf.s[19] = PI*( trf.s[2] + 2.*trf.s[4] )*trf.s[8];

   trf.s[10] = sqrt( PI*MY_VAC_*trf.fr/trf.s[9] ); /* R_square [ ohms ] */
   trf.s[11] = 1./sqrt( PI*MY_VAC_*trf.fr*trf.s[9] ); /* penetration depth */

   trf.rsqr1 = trf.s[10];
/*............................................................................*/
# if DSC_HCRMDE != 0
/* A delay between the electric and heat conductive computations */
/* may be useful [ to speed up the algorithm, for instance ] */

   ( spt->hcrstart ) = ZERO;

/* default time step: */
   if (( spt->hcdt ) < 1.e-277 )
      ( spt->hcdt ) = trf.s[32]; /* proposed time step [ default ] */

/* heat conductivity, mantle: */
   trf.kh1 = trf.s[12]; 
/* heat capacity per volume (!), mantle: */
   trf.cv1 = trf.s[13];
/* environment-surface heat conductivity [W/(K*m^2)] */
   trf.shc = trf.s[15]; 
/* environment temperature [Celsius]: */
   trf.tev = trf.s[16]; 
/* contact temperature, cooling channel: */
   trf.temp = trf.s[14];
/* gas parameters: */ 
/*
*//* trf.kh0: heat conductivity [W/(K*m)]
*//* trf.cv0: heat capacity [J/(K*m^3) !!!]
*/
   trf.kh0 = trf.s[24];
   trf.cv0 = trf.s[21]*trf.s[25];
/*............................................................................*/
# if DSC_FLDMDE != 0
/* fluid parameters: */ 
/*
*//* trf.tm0: mean temperature [Celsius] 
*//* trf.rm0: mean density [Kg/(m^3)]
*//* trf.bm0: thermal expansion coefficient [1/K]
*//* trf.cm0: [ adiabatic ] compression coefficient [1/Pa]
*//* trf.ny0: dynamic viscosity [Kg/(m*sec)]
*//* trf.q1: Cp/Cv - 1 [dimensionless]
*//* trf.td: dissipation time constant [sec]
*//* trf.gr: gravitation acceleration [m/(sec^2)]
*//* trf.LL: a reference length [ for turbulence models - presently not used ]
*/
/* gas parameters: */ 

   if ( trf.c[THREE] == null )
   {
      trf.s[22] = ZERO;
      trf.s[26] = ZERO;
      trf.s[27] = ZERO;
   };

   trf.tm0 = trf.s[20];
   trf.rm0 = trf.s[21];
   trf.bm0 = trf.s[22];
   trf.cm0 = trf.s[23]; 
   trf.ny0 = trf.s[26];

   trf.q1 = 0.4;

   trf.gr = trf.s[27];
   trf.td = trf.s[28];
   trf.LL = trf.s[29];

/* coarsening time period */
   trf.crsdt = trf.s[31]; 
   ( spt->ppt->fcp->crsdt[ONE] ) = trf.crsdt;

/* Nusselt number at no-slip boundary [ optional - here not used ] */
   trf.nus = trf.s[33]; 

/* A delay between the electric and advective [ fluid dynamic ] */
/* computations may be useful [ to speed up the algorithm, e.g.]: */
   ( spt->fldstart ) = trf.s[30];

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   return ONE;
}
/*=================== end of function rvise_params(*) ========================*/

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
   blc.z[1] = blc.z[0] + trf.s[6]; /* source line */

   if ( trf.ref == null )
   {
      blc.z[2] = blc.z[1] + trf.s[7]; /* + buffer */
      blc.z[3] = blc.z[2] + trf.s[8]; /* + length, cooled section */
      blc.z[4] = blc.z[3] + trf.s[7]; /* + lenght, termination */
   };

   return;
}
/*=================== end of function set_z_bases(*) ========================*/

/*******************************************************************************
*                                                                              *
*   Function cords(*)                                                          *
*   [ DSC model prototype: 'convection in coaxial line'; DANSE release 1.0.]   *
*                                                                              *
*   This function defines, in a DSC model dependent way, coordinates of        *
*   support points pnt[n][]  ( n = 0,...,MOD_POINTS-ONE ), viz.:               *
*                                                                              *
*                pnt[n][0] = x, pnt[n][1] = y, pnt[n][2] = z,                  *
*                                                                              *
*   used to build up the block structure in function 'blocks(*)'               *
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

      zz   = ZERO,
      zz0 = ZERO,
      alfa = ZERO,
      ss = ZERO;

   static short
      ii = null,
      jj = null, 
      points = null;

/* prototypes: */

   double 
      sin( double ),
      cos( double ),
      sqrt( double );      
/*---------------------------------------------------------------------------*/
/* This part is canonical: take care with modifying it [ better you leave it */
/* untouched ] */

   ss = ZERO;
   zz = ZERO;

   jj = null;
   ii = null;
   while ( ii < layer )
   {
      while (( ii == blc.base[jj+ONE] )&&( jj < blc.m[null] ))
      {
         jj += ONE;
      };
      zz += blc.dz[jj];
      ii++;
   };
/*............................................................................*/
/* radii: --- */

/*............................................................................*/
/* bows [ arcs ]: --- */

   trf.alfa[null] = PI;
   ss = (( double ) trf.n[1] )/trf.n[5]; 
   trf.alfa[1] = ss*trf.alfa[null];
   ss = (( double ) trf.n[3] )/trf.n[5]; 
   trf.alfa[2] = ss*trf.alfa[null];
   trf.alfa[3] = trf.alfa[1];
   ss = (( double ) trf.n[1] )/trf.n[6]; 
   trf.alfa[4] = ss*trf.alfa[null];
   ss = 2.*((( double ) trf.n[2] )/trf.n[6] ); 
   trf.alfa[5] = ss*trf.alfa[null];
   ss = (( double ) trf.n[3] )/trf.n[6]; 
   trf.alfa[6] = ss*trf.alfa[null];
   trf.alfa[7] = trf.alfa[5];
   trf.alfa[8] = trf.alfa[4];
/*............................................................................*/
/* xy-coordinates: */

   xx0 = ZERO;   
   yy0 = ZERO;

   alfa = - trf.alfa[null]/2.;
   xx1 = trf.ri*cos( alfa );
   yy1 = trf.ri*sin( alfa );

   alfa += trf.alfa[1];
   xx2 = trf.ri*cos( alfa );
   yy2 = trf.ri*sin( alfa );

   alfa += trf.alfa[2];
   xx3 = trf.ri*cos( alfa );
   yy3 = trf.ri*sin( alfa );

   alfa += trf.alfa[3];
   xx4 = trf.ri*cos( alfa );
   yy4 = trf.ri*sin( alfa );

   alfa = - trf.alfa[null]/2.;
   xx5 = trf.ra*cos( alfa );
   yy5 = trf.ra*sin( alfa );

   alfa = - trf.alfa[null]/2.;
   xx13 = ( trf.ri - trf.s[5] )*cos( alfa );
   yy13 = ( trf.ri - trf.s[5] )*sin( alfa );

   alfa = trf.alfa[null]/2.;
   xx14 = ( trf.ri - trf.s[5] )*cos( alfa );
   yy14 = ( trf.ri - trf.s[5] )*sin( alfa );

   alfa = - trf.alfa[null]/2 + trf.alfa[4];
   xx6 = trf.ra*cos( alfa );
   yy6 = trf.ra*sin( alfa );

   alfa += trf.alfa[5];
   xx7 = trf.ra*cos( alfa );
   yy7 = trf.ra*sin( alfa );

   alfa += trf.alfa[6];
   xx8 = trf.ra*cos( alfa );
   yy8 = trf.ra*sin( alfa );

   alfa += trf.alfa[7];
   xx9 = trf.ra*cos( alfa );
   yy9 = trf.ra*sin( alfa );

   alfa += trf.alfa[8];
   xx10 = trf.ra*cos( alfa );
   yy10 = trf.ra*sin( alfa );

   alfa = - trf.alfa[null]/2.;
   xx11 = ( trf.ra + trf.s[4] )*cos( alfa );
   yy11 = ( trf.ra + trf.s[4] )*sin( alfa );

   alfa = trf.alfa[null]/2.;
   xx12 = ( trf.ra + trf.s[4] )*cos( alfa );
   yy12 = ( trf.ra + trf.s[4] )*sin( alfa );
/*............................................................................*/
/* z-coordinates: */

   zz0 = zz;
/*............................................................................*/
/* The total number of points: */

   points = 15;
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
/*====================== end of function cords(*) ===========================*/



/*******************************************************************************
*                                                                              *
*   Block definition function                                                  *
*   [ DSC model prototype: 'convection in coaxial line'; DANSE release 1.0.]   *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 29, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "./tools/blcutls.h" /* macros BLC_EXCITE(*), _BOUNDARY(*), */
                              /* BLC_TRIVIAL(*), and BLC_DECLARE(*)  */
# include "./tools/setmed.h"  /* media switching function setmed(*) */
/**************************** macro SETMED(*) *********************************/
# if DSC_HCRMDE != 0
/*............................................................................*/
# if DSC_FLDMDE != 0
# define SETMED(NN) \
{ \
   ( mdp->idx ) = (NN); \
 \
   if ((NN) == TWO ) \
   { \
      ( mdp->ci ) = lbl.ci; /* initial cell index */\
      ( mdp->cf ) = lbl.cf; /* final cell index */\
 \
      ( mdp->cnn ) = ONE;   /* index of fluid [ connected ] component */\
 \
      ( mdp->eps ) = 1.;    /* relative permittivity */\
      ( mdp->myr ) = 1.;    /* relative permeability */\
      ( mdp->ke ) = ZERO;   /* electric conductivity [S/m] */\
      ( mdp->km ) = ZERO;   /* magnetic conductivity [O/m] */\
 \
      ( mdp->kh ) = trf.kh0; /* thermal conductivity [W/(m*K)] */\
      ( mdp->cv ) = trf.cv0; /* specific heat per volume [J/(m^3)] */\
 \
      ( mdp->rm ) = trf.rm0; /* mean fluid density [Kg/(m^3)]*/\
      ( mdp->tm ) = trf.tm0; /* mean fluid temperature [DEG Celsius] */\
      ( mdp->bm ) = trf.bm0; /* thermal expansion coefficient [1/K] */\
      ( mdp->cm ) = trf.cm0; /* adiabatic compression coefficient [1/Pa] */\
      ( mdp->ny ) = trf.ny0; /* dynamic viscosity [Kg/(s*m)] */\
      ( mdp->q1 ) = trf.q1;  /* Cp/Cv - 1 [ dimensionless ] */\
      ( mdp->td ) = trf.td; /* relaxation time [s]; usually ZERO */\
      ( mdp->LL ) = trf.LL; /* characteristic length [m]; usulally ZERO */\
      ( mdp->gr[0] ) = ZERO; /* gravit acceleration [m/(s^2)]; x-component */\
      ( mdp->gr[1] ) = - trf.gr; /* dito; y-component */\
      ( mdp->gr[2] ) = ZERO; /* dito; z-component */\
      ( mdp->gp[0] ) = ZERO; /* pressure gradient [N/(m^2)]; x-component */\
      ( mdp->gp[1] ) = ZERO; /* dito; y-component */\
      ( mdp->gp[2] ) = ZERO; /* dito; z-component */\
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
/*............................................................................*/
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
/*............................................................................*/
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
/************************ end of macro SETMED(*) ******************************/

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
/*
   ptr[null] = null; 
*/
/*............................................................................*/
   if( layer == null ) /* don't remove */
   {
      mdp = setmed( NULL );
      ( mdp->opt ) = *option;

      trf.exc = null;
      trf.bnd = null;
      trf.val = null;
/*............................................................................*/
# if DSC_HCRMDE != 0 
      trf.bhc = null;
      trf.ehc = null;
      trf.vhc = null;
# endif
/*............................................................................*/
   };
/*............................................................................*/
/*
    printf( "\n layer %d ", layer );
    scanf( "%s", ptr );
*/
/*.............................. waveguide blocks ............................*/
   BLOCK(1);

   blc.ab = trf.n[1];
   blc.bc = trf.n[2];

   if (( layer < blc.base[2] )
     ||( blc.base[3] <= layer ))
   {
      CONNECT( 1, ab, 0, blc.ab,
               0, ewall, 0, 0 );

      CONNECT( 1, dc, 0, blc.ab,
               0, ewall, 0, 0 );
   };

   CONNECT( 1, ad, 0, blc.bc,
            0, mwall, 0, 0 );
   
   if ( layer == null )
   {
/* the port excitation [ coaxial port ]: */

      ii = lbl.m+ONE; /* the initial mesh cell [ index ] */
      nn = blc.ab*blc.bc; /* the number of periods */
      pp = ONE; /* the period */

      BLC_EXCITE( "e", ii, nn, pp, FOUR, trf.uu ); /* "e": E field excitation */

/* boundary conditions at that port: */

      BLC_BOUNDARY( "coax", ii, nn, pp, FOUR, 1. );
   }
   else if ( layer == ( blc.base[blc.m[null]] - ONE ))
   {
/* boundary conditions at termination: */

      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE;

      BLC_BOUNDARY( "coax", ii, nn, pp, FIVE, 1. );
   };

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
/* fluid flow boundary conditions */
/* at outer conductor */

      ii = lbl.m+ONE;
      nn = blc.ab;
      pp = ONE;

      if ( trf.c[6] == null )
      {
         BLC_BOUNDARY( "no_slip", ii, nn, pp, TWO, trf.s[33] );
      }
      else
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, TWO, ZERO );
      };

/* at inner conductor */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
      nn = blc.ab;
      pp = ONE;

      if ( trf.c[6] == null )
      {
         BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, trf.s[33] );
      }
      else
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, THREE, ZERO );
      };

/* free slip fluid conditions: */

      ii = lbl.m+ONE;
      nn = blc.bc;
      pp = blc.ab;

      BLC_BOUNDARY( "slip", ii, nn, pp, null, ZERO );

/* outer skin effect lossy boundary: */

      ii = lbl.m+ONE; /* the initial mesh cell [ index ] */
      nn = blc.ab; /* the number of periods */
      pp = ONE; /* the period */

      BLC_BOUNDARY( "skin", ii, nn, pp, TWO, trf.rsqr1 );

/* the pertinent thermal sources: */

      if (( trf.c[4] == ONE ) /* apply skin effect thermal losses to gas */
        &&(( trf.c[5] == ONE )
         ||( trf.c[5] == TWO )))
      {
         BLC_BOUNDARY( "sc", ii, nn, pp, TWO, trf.rsqr1 );
      };

/* inner skin effect lossy boundary: */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
      nn = blc.ab;
      pp = ONE;

      BLC_BOUNDARY( "skin", ii, nn, pp, THREE, trf.rsqr1 );

/* the pertinent thermal sources: */

      if (( trf.c[4] == ONE ) /* apply skin effect thermal losses to gas */
        &&(( trf.c[5] == null )
         ||( trf.c[5] == TWO )))
      {
         BLC_BOUNDARY( "sc", ii, nn, pp, THREE, trf.rsqr1 );
      };
   };

   if ( layer == blc.base[2] )
   {
/* fluid flow boundary conditions: */

      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE; 

      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "inflow", ii, nn, pp, FOUR, ZERO );
      };
   };

   if ( layer == ( blc.base[3]-ONE ))
   { 
/* fluid flow boundary conditions: */

      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE;

      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "outflow", ii, nn, pp, FIVE, ZERO );
      };
   };

   if ( layer == ( blc.base[ONE] - ONE ))
   {
/* the E field evaluation: */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
      BLC_EVALUATE( "e", ii, ONE, ONE, 8, "reference" ); /* 'e': E field */
   };
/*............................................................................*/
# if MOD_DEFLT == 3
   ii = ONE;
# else
   ii = blc.m[3]/2;
# endif
/*............................................................................*/
   if ( layer == ( blc.base[2]+ii ))
   {
/* temperature evaluation: */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
      BLC_EVALUATE( "tf", ii, ONE, ONE, THREE, "inner_1" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_inner_1" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_inner_1" );
# endif
# if MOD_FLOWS != 0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_inner_1" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_inner_1" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_inner_1" );
# endif
      ii = lbl.m+ONE;
      BLC_EVALUATE( "tf", ii, ONE, ONE, TWO, "outer_1" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_outer_1" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_outer_1" );
# endif
# if MOD_FLOWS != 0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_outer_1" );
# endif
      ii = lbl.m+blc.ab*blc.bc/2;
      BLC_EVALUATE( "tn", ii, ONE, ONE, null, "bulk_1" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_bulk_1" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_bulk_1" );
# endif
# if MOD_FLOWS != 0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_bulk_1" );
# endif
   };

   if ( *option == 'c' )
   {
      qdl.aab = trf.alfa[4];
      qdl.adc = trf.alfa[1];

      POINT( a, 5 );
      POINT( b, 6 );
      POINT( c, 2 );
      POINT( d, 1 );
   };

   QUDRL(1);

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
      SETMED(2);
   };
/*............................................................................*/
   BLOCK(2);

   blc.ab = trf.n[2];
   blc.bc = 2*trf.n[2];

   CONNECT( 2, ab, null, blc.ab,
            1, bc, blc.ab, backward );

   if (( layer < blc.base[2] )
     ||( blc.base[3] <= layer ))
   {
      CONNECT( 2, bc, 0, blc.bc,
               0, ewall, 0, 0 );
   };

   if ( layer == null )
   {
/* the excitation: */

      ii = lbl.m+ONE; /* the initial mesh cell [ index ] */
      nn = ( blc.ab+ONE )*blc.ab/2; /* the number of periods */
      pp = ONE; /* the period */

      BLC_EXCITE( "e", ii, nn, pp, FOUR, trf.uu ); /* E field excitation */

/* the boundary: */

      BLC_BOUNDARY( "coax", ii, nn, pp, FOUR, 1. );
   }
   else if ( layer == ( blc.base[blc.m[null]] - ONE ))
   {
      ii = lbl.m+ONE;
      nn = ( blc.ab+ONE )*blc.ab/2; /* the number of periods */
      pp = ONE;

      BLC_BOUNDARY( "coax", ii, nn, pp, FIVE, 1. );
   };

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
      ii = lbl.m;
      jj = blc.ab; do
      {
         ii += jj; /* the initial mesh cell [ index ] */
         nn = ONE; /* the number of periods */
         pp = ONE; /* the period */

/* fluid flow boundary conditions */
/* at outer conductor */

         if ( trf.c[6] == null ) /* no-slip */
         {
            BLC_BOUNDARY( "no_slip", ii, nn, pp, ONE, trf.s[33] );
            BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, trf.s[33] );
         }
         else /* free slip */
         {
            BLC_BOUNDARY( "slip", ii, nn, pp, ONE, ZERO );
            BLC_BOUNDARY( "slip", ii, nn, pp, THREE, ZERO );
         };

/* outer skin effect lossy boundaries: */

         BLC_BOUNDARY( "skin", ii, nn, pp, ONE, trf.rsqr1 );
         BLC_BOUNDARY( "skin", ii, nn, pp, THREE, trf.rsqr1 );

         if (( trf.c[4] == ONE ) /* apply skin effect thermal losses to gas */
           &&(( trf.c[5] == ONE )
            ||( trf.c[5] == TWO )))
         {    
            BLC_BOUNDARY( "sc", ii, nn, pp, ONE, trf.rsqr1 );
            BLC_BOUNDARY( "sc", ii, nn, pp, THREE, trf.rsqr1 );
         };
      } while ( null < ( --jj ));
   }; /* end if ... */

   if ( layer == blc.base[2] )
   { 
/* fluid flow boundary conditions: */

      ii = lbl.m+ONE;
      nn = ( blc.ab+ONE )*blc.ab/2;
      pp = ONE; 

      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "inflow", ii, nn, pp, FOUR, ZERO );
      };
   };

   if ( layer == ( blc.base[3]-ONE ))
   { 
/* fluid flow boundary conditions: */
	   
      ii = lbl.m+ONE;
      nn = ( blc.ab+ONE)*blc.ab/2;
      pp = ONE; 

      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "outflow", ii, nn, pp, FIVE, ZERO );
      };
   };

   if ( *option == 'c' )
   {
      tgl.abc = trf.alfa[5];

      POINT( c, 7 );
   };

   TRNGL(2);

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
      SETMED(2);
   };
/*............................................................................*/
   BLOCK(3);

   blc.ab = trf.n[3];
   blc.bc = trf.n[2];

   CONNECT( 3, ad, 0, blc.bc,
            2, ac, blc.bc, backward );

   if (( layer < blc.base[2] )
     ||( blc.base[3] <= layer ))
   {
      CONNECT( 3, ab, 0, blc.ab,
               0, ewall, 0, 0 );
      
      CONNECT( 3, dc, 0, blc.ab,
               0, ewall, 0, 0 );
   };

   if ( layer == null )
   {
/* the excitation: */

      ii = lbl.m+ONE; /* the initial mesh cell [ index ] */
      nn = blc.ab*blc.bc; /* the number of periods */
      pp = ONE; /* the period */

      BLC_EXCITE( "e", ii, nn, pp, FOUR, trf.uu );

/* the boundary: */

      BLC_BOUNDARY( "coax", ii, nn, pp, FOUR, 1. );
   }
   else if ( layer == ( blc.base[blc.m[null]] - ONE ))
   {
      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE;

      BLC_BOUNDARY( "coax", ii, nn, pp, FIVE, 1. );
   };

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
/* fluid flow boundary conditions */
/* at outer conductor */

      ii = lbl.m+ONE;
      nn = blc.ab;
      pp = ONE;

      if ( trf.c[6] == null ) /* no-slip */
      {
         BLC_BOUNDARY( "no_slip", ii, nn, pp, TWO, trf.s[33] );
      }
      else /* free slip */
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, TWO, ZERO );
      };

/* at inner conductor */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
      nn = blc.ab;
      pp = ONE;

      if ( trf.c[6] == null )
      {
         BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, trf.s[33] );
      }
      else
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, THREE, ZERO );
      };

/* outer skin effect lossy boundary: */

      ii = lbl.m+ONE; /* the initial mesh cell [ index ] */
      nn = blc.ab; /* the number of periods */
      pp = ONE; /* the period */

      BLC_BOUNDARY( "skin", ii, nn, pp, TWO, trf.rsqr1 );
      
/* the pertinent thermal sources: */

      if (( trf.c[4] == ONE ) /* apply skin effect losses to gas */
        &&(( trf.c[5] == ONE )
         ||( trf.c[5] == TWO )))
      {
         BLC_BOUNDARY( "sc", ii, nn, pp, TWO, trf.rsqr1 );
      };

/* inner skin effect lossy boundary: */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
      nn = blc.ab;
      pp = ONE;

      BLC_BOUNDARY( "skin", ii, nn, pp, THREE, trf.rsqr1 );

/* the pertinent thermal sources: */

      if (( trf.c[4] == ONE ) /* apply skin effect losses to gas */
        &&(( trf.c[5] == null )
         ||( trf.c[5] == TWO )))
      {
         BLC_BOUNDARY( "sc", ii, nn, pp, THREE, trf.rsqr1 ); 
      };
   };
  
   if ( layer == blc.base[2] )
   { 
/* fluid flow boundary conditions: */

      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE; 

      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "inflow", ii, nn, pp, FOUR, ZERO );
      };
   };

   if ( layer == ( blc.base[3]-ONE ))
   { 
/* fluid flow boundary conditions: */
	   
      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE; 
      
      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "outflow", ii, nn, pp, FIVE, ZERO );
      };
   };
/*............................................................................*/
# if MOD_DEFLT == 3
   ii = ONE;
# else
   ii = blc.m[3]/2;
# endif
/*............................................................................*/
   if ( layer == ( blc.base[2]+ii ))
   {
/* temperature evaluation: */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + blc.ab/2;
      BLC_EVALUATE( "tf", ii, ONE, ONE, THREE, "inner_2" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_inner_2" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_inner_2" );
# endif
# if MOD_FLOWS !=0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_inner_2" );
# endif

      ii = lbl.m+blc.ab/2;
      BLC_EVALUATE( "tf", ii, ONE, ONE, TWO, "outer_2" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_outer_2" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_outer_2" );
# endif
# if MOD_FLOWS != 0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_outer_2" );
# endif

      ii = lbl.m+blc.ab*blc.bc/2;
      BLC_EVALUATE( "tn", ii, ONE, ONE, null, "bulk_2" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_bulk_2" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_bulk_2" );
# endif
# if MOD_FLOWS != 0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_bulk_2" );
# endif
   };

   if ( *option == 'c' )
   {
      qdl.adc = trf.alfa[2];
      qdl.aab = trf.alfa[6];

      POINT( b, 8 );
      POINT( c, 3 );
   };

   QUDRL(3);

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
      SETMED(2);
   };
/*............................................................................*/
   BLOCK(4);

   blc.ab = trf.n[2];
   blc.bc = 2*trf.n[2];

   CONNECT( 4, ab, null, blc.ab,
            3, bc, blc.ab, backward );

   if (( layer < blc.base[2] )
     ||( blc.base[3] <= layer ))
   {
      CONNECT( 4, bc, 0, blc.bc,
               0, ewall, 0, 0 );
   };

   if ( layer == null )
   {
/* the excitation: */

      ii = lbl.m+ONE; /* the initial mesh cell [ index ] */
      nn = ( blc.ab+ONE )*blc.ab/2; /* the number of periods */
      pp = ONE; /* the period */

      BLC_EXCITE( "e", ii, nn, pp, FOUR, trf.uu ); /* E-field excitation */

/* the boundary: */

      BLC_BOUNDARY( "coax", ii, nn, pp, FOUR, 1. );
   }
   else if ( layer == ( blc.base[blc.m[null]] - ONE ))
   {
      ii = lbl.m+ONE;
      nn = ( blc.ab+ONE )*blc.ab/2; /* the number of periods */
      pp = ONE;

      BLC_BOUNDARY( "coax", ii, nn, pp, FIVE, 1. );
   };

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
      ii = lbl.m;
      jj = blc.ab; do
      {
         ii += jj; /* the initial mesh cell [ index ] */
         nn = ONE; /* the number of periods */
         pp = ONE; /* the period */

/* fluid flow boundary conditions */
/* at outer conductor */

         if ( trf.c[6] == null ) /* no-slip */
         {
            BLC_BOUNDARY( "no_slip", ii, nn, pp, ONE, trf.s[33] );
            BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, trf.s[33] );
         }
         else /* free slip */
         {
            BLC_BOUNDARY( "slip", ii, nn, pp, ONE, ZERO );
            BLC_BOUNDARY( "slip", ii, nn, pp, THREE, ZERO );
         };

/* outer skin effect lossy boundaries: */

         BLC_BOUNDARY( "skin", ii, nn, pp, ONE, trf.rsqr1 );
         BLC_BOUNDARY( "skin", ii, nn, pp, THREE, trf.rsqr1 );

         if (( trf.c[4] == ONE ) /* apply skin effect thermal losses to gas */
           &&(( trf.c[5] == ONE )
            ||( trf.c[5] == TWO )))
         {
            BLC_BOUNDARY( "sc", ii, nn, pp, ONE, trf.rsqr1 );
            BLC_BOUNDARY( "sc", ii, nn, pp, THREE, trf.rsqr1 );
         };
      } while ( null < ( --jj ));
   }; /* end if ... */

   if ( layer == blc.base[2] )
   { 
/* fluid flow boundary conditions: */

      ii = lbl.m+ONE;
      nn = ( blc.ab+ONE)*blc.ab/2;
      pp = ONE; 

      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "inflow", ii, nn, pp, FOUR, ZERO );
      };
   };

   if ( layer == ( blc.base[3]-ONE ))
   { 
/* fluid flow boundary conditions: */
	   
      ii = lbl.m+ONE;
      nn = ( blc.ab+ONE)*blc.ab/2;
      pp = ONE; 

      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "outflow", ii, nn, pp, FIVE, ZERO );
      };
   };

   if ( *option == 'c' )
   {
      tgl.abc = trf.alfa[7];
      
      POINT( c, 9 );
   };

   TRNGL(4);

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
      SETMED(2);
   };
/*............................................................................*/
   BLOCK(5);

   blc.ab = trf.n[1];
   blc.bc = trf.n[2];

   CONNECT( 5, ad, 0, blc.bc,
            4, ac, blc.bc, backward );

   if (( layer < blc.base[2] )
     ||( blc.base[3] <= layer ))
   {
      CONNECT( 5, ab, 0, blc.ab,
               0, ewall, 0, 0 );

      CONNECT( 5, dc, 0, blc.ab,
               0, ewall, 0, 0 );
   };

   CONNECT( 5, bc, 0, blc.bc,
            0, mwall, 0, 0 );

   if ( layer == null )
   {
/* the excitation: */

      ii = lbl.m+ONE; /* the initial mesh cell [ index ] */
      nn = blc.ab*blc.bc; /* the number of periods */
      pp = ONE; /* the period */

      BLC_EXCITE( "e", ii, nn, pp, FOUR, trf.uu ); /* E-field excitation */

/* the boundary: */

      BLC_BOUNDARY( "coax", ii, nn, pp, FOUR, 1. );
   }
   else if ( layer == ( blc.base[blc.m[null]] - ONE ))
   {
      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE;

      BLC_BOUNDARY( "coax", ii, nn, pp, FIVE, 1. );
   };

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
/* fluid flow boundary conditions */
/* at outer conductor */

      ii = lbl.m+ONE;
      nn = blc.ab;
      pp = ONE;

      if ( trf.c[6] == null )
      {
         BLC_BOUNDARY( "no_slip", ii, nn, pp, TWO, trf.s[33] );
      }
      else
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, TWO, ZERO );
      };

/* at inner conductor */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
      nn = blc.ab;
      pp = ONE;

      if ( trf.c[6] == null )
      {
         BLC_BOUNDARY( "no_slip", ii, nn, pp, THREE, trf.s[33] );
      }
      else
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, THREE, ZERO );
      };

/* free slip conditions at symmetry plane: */

      ii = lbl.m+blc.ab;
      nn = blc.bc;
      pp = blc.ab;

      BLC_BOUNDARY( "slip", ii, nn, pp, ONE, ZERO );

/* outer skin effect lossy boundary: */

      ii = lbl.m+ONE; /* the initial mesh cell [ index ] */
      nn = blc.ab; /* the number of periods */
      pp = ONE; /* the period */

      BLC_BOUNDARY( "skin", ii, nn, pp, TWO, trf.rsqr1 );

/* the pertinent thermal sources: */

      if (( trf.c[4] == ONE ) /* apply skin effect losses to gas */
        &&(( trf.c[5] == ONE )
         ||( trf.c[5] == TWO )))
      {
         BLC_BOUNDARY( "sc", ii, nn, pp, TWO, trf.rsqr1 );
      };

/* inner skin effect lossy boundary: */

      ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
      nn = blc.ab;
      pp = ONE;

      BLC_BOUNDARY( "skin", ii, nn, pp, THREE, trf.rsqr1 );

/* the pertinent thermal sources: */

      if (( trf.c[4] == ONE ) /* apply skin effect losses to gas */
        &&(( trf.c[5] == null )
         ||( trf.c[5] == TWO )))
      {
         BLC_BOUNDARY( "sc", ii, nn, pp, THREE, trf.rsqr1 ); 
      };
   };
  
   if ( layer == blc.base[2] )
   { 
/* fluid flow boundary conditions: */

      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE; 

      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FOUR, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "inflow", ii, nn, pp, FOUR, ZERO );
      };
   };

   if ( layer == ( blc.base[3]-ONE ))
   { 
/* fluid flow boundary conditions: */
	   
      ii = lbl.m+ONE;
      nn = blc.ab*blc.bc;
      pp = ONE; 
      
      if ( trf.c[7] == ONE )
      {
         BLC_BOUNDARY( "slip", ii, nn, pp, FIVE, ZERO );
      }
      else if ( trf.c[7] == TWO )
      {
         BLC_BOUNDARY( "outflow", ii, nn, pp, FIVE, ZERO );
      };
   };
/*............................................................................*/
# if MOD_DEFLT == 3
   ii = ONE;
# else
   ii = blc.m[3]/2;
# endif
/*............................................................................*/
   if ( layer == ( blc.base[2]+ii ))
   {
/* temperature evaluation: */

      ii = lbl.m+blc.ab*blc.bc;
      BLC_EVALUATE( "tf", ii, ONE, ONE, THREE, "inner_3" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_inner_3" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_inner_3" );
# endif
# if MOD_FLOWS != 0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_inner_3" );
# endif
      
      ii = lbl.m+blc.ab;
      BLC_EVALUATE( "tf", ii, ONE, ONE, TWO, "outer_3" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_outer_3" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_outer_3" );
# endif
# if MOD_FLOWS != 0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_outer_3" );
# endif

      ii = lbl.m+blc.ab*blc.bc/2;
      BLC_EVALUATE( "tn", ii, ONE, ONE, null, "bulk_3" );
# if MOD_FLOWS == 2
      BLC_EVALUATE( "un", ii, ONE, ONE, 'x', "x-flow_bulk_3" );
      BLC_EVALUATE( "un", ii, ONE, ONE, 'z', "z-flow_bulk_3" );
# endif
# if MOD_FLOWS != 0
      BLC_EVALUATE( "un", ii, ONE, ONE, 'y', "y-flow_bulk_3" );
# endif
   };

   if ( *option == 'c' )
   {
      qdl.adc = trf.alfa[3];
      qdl.aab = trf.alfa[8];

      POINT( b, 10 );
      POINT( c, 4 );
   };

   QUDRL(5);

   if (( blc.base[2] <= layer )
     &&( layer < blc.base[3] ))
   {
      SETMED(2);
   };
/*............................................................................*/
   if (( blc.base[2] <= layer )
     &&( layer <= blc.base[3] ))
   {
      ii = leave( 5, blc.base[3], blc.base[blc.m[null]], layer );

      if( ii == ONE )
         goto terminal;
/*............................................................................*/
      BLOCK(6);

      blc.ab = trf.n[6];
      blc.bc = trf.n[4];

      jj = null;
      CONNECT( 6, dc, 0, trf.n[1],
               1, ab, 0, forward );

      jj += trf.n[1];
      CONNECT( 6, dc, jj, jj+2*trf.n[2],
               2, bc, 0, forward );

      jj += 2*trf.n[2];
      CONNECT( 6, dc, jj, jj+trf.n[3],
               3, ab, 0, forward );

      jj += trf.n[3];
      CONNECT( 6, dc, jj, jj+2*trf.n[2],
               4, bc, 0, forward );

      jj += 2*trf.n[2];
      CONNECT( 6, dc, jj, jj+trf.n[1],
               5, ab, 0, forward );

      CONNECT( 6, ab, 0, blc.ab,
               0, ewall, 0, 0 );

      CONNECT( 6, bc, 0, blc.bc,
               0, mwall, 0, 0 );

      CONNECT( 6, ad, 0, blc.bc,
               0, mwall, 0, 0 );

      if ( layer == blc.base[2] )
         qdl.bot = e_wall;

      if ( layer == ( blc.base[3] - ONE ))
         qdl.top = e_wall;
      
/* outer thermal sources pertinent to skin effect lossy boundary */

      if (( trf.c[4] == null ) /* apply skin effect losses to metal */
        &&(( trf.c[5] == ONE )
         ||( trf.c[5] == TWO )))
      {
         ii = lbl.m + ( blc.bc-ONE )*blc.ab + ONE;
         nn = blc.ab;
         pp = ONE;
         BLC_BOUNDARY( "sc", ii, nn, pp, THREE, trf.rsqr1 );
      };
/*............................................................................*/
# if MOD_ENBLESC == 1
/* convection [ environment-surface heat conductance ]: */

      ii = lbl.m+ONE;
      nn = blc.ab;
      pp = ONE;
/*
      BLC_BOUNDARY( "sf", ii, nn, pp, TWO, trf.tev );
*/
      BLC_BOUNDARY( "tf", ii, nn, pp, TWO, trf.s[14] );
# endif
/*............................................................................*/
/* temperature evaluation: */

      if ( layer == ( blc.base[2] + blc.m[3]/2 ))
      {
/* evaluate mantle face temperatures: */

         BLC_EVALUATE( "tf", ( lbl.m + ONE ), ONE, ONE, TWO, "mantle_1" );
         BLC_EVALUATE( "tf", ( lbl.m + blc.ab/2 ), ONE, ONE, TWO, "mantle_2" );
         BLC_EVALUATE( "tf", ( lbl.m + blc.ab ), ONE, ONE, TWO, "mantle_3" );
      };

      if ( *option == 'c' )
      {
         qdl.aab = trf.alfa[null];
         qdl.adc = trf.alfa[null];

         POINT( a, 11 );
         POINT( b, 12 );
      };

      QUDRL(6);

      if (( blc.base[2] <= layer )
        &&( layer < blc.base[3] ))
      {
         SETMED(3);
      };
/*............................................................................*/
      BLOCK(7);

      blc.ab = trf.n[5];
      blc.bc = trf.n[4];

      jj = null;
      CONNECT( 7, ab, 0, trf.n[1],
               1, dc, 0, forward );

      jj += trf.n[1];
      CONNECT( 7, ab, jj, jj+trf.n[3],
               3, dc, 0, forward );

      jj += trf.n[3];
      CONNECT( 7, ab, jj, jj+trf.n[1],
               5, dc, 0, forward );

      CONNECT( 7, dc, 0, blc.ab,
               0, ewall, 0, 0 );

      CONNECT( 7, bc, 0, blc.bc,
               0, mwall, 0, 0 );

      CONNECT( 7, ad, 0, blc.bc,
               0, mwall, 0, 0 );

      if ( layer == blc.base[2] )
         qdl.bot = e_wall;

      if ( layer == ( blc.base[3] - ONE ))
         qdl.top = e_wall;
      
/* inner thermal sources pertinent to skin effect lossy boundary */

      if (( trf.c[4] == null ) /* apply skin effect losses to metal */
        &&(( trf.c[5] == null )
         ||( trf.c[5] == TWO )))
      {
         ii = lbl.m+ONE;
         nn = blc.ab;
         pp = ONE;
         BLC_BOUNDARY( "sc", ii, nn, pp, TWO, trf.rsqr1 );
      };
/*............................................................................*/
/* temperature evaluation: */

      if ( layer == ( blc.base[2] + blc.m[3]/2 ))
      {
/* evaluate mantle face temperatures: */
/*
         BLC_EVALUATE( "tf", ( lbl.m + ONE ), ONE, ONE, TWO, "tube_1" );
         BLC_EVALUATE( "tf", ( lbl.m + blc.ab/2 ), ONE, ONE, TWO, "tube_2" );
         BLC_EVALUATE( "tf", ( lbl.m + blc.ab ), ONE, ONE, TWO, "tube_3" );
*/
      };

      if ( *option == 'c' )
      {
         qdl.aab = trf.alfa[null];
         qdl.adc = trf.alfa[null];

         POINT( d, 13 );
         POINT( c, 14 );
      };

      QUDRL(7);

      if (( blc.base[2] <= layer )
        &&( layer < blc.base[3] ))
      {
         SETMED(3);
      };
/*............................................................................*/
   }; /* end if layers in cooled waveguide section */
/*
   goto terminal;
*/
/*............................................................................*/

  terminal:

   BLC_LIMITS( 'c' ); /* don't remove on modifying DSC model !!!              */

   return null;
} 
/*====================== end of function blocks(*) ===========================*/


/* [ function: systop ] */
/*******************************************************************************
*                                                                              *
*   ANSI-C function systop(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0.]                                                            *
*                                                                              *
*   Model dependent DSC mesh generation subroutine                             *
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
/*----------------------------------------------------------------------------*/
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
*   ANSI-C function syssmx(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0.]                                                            *
*                                                                              *
*   DSC model dependent coordinates and media transfer function                *
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

 /* ptr = ( char * ) calloc( STS_SIZE , ONE ); */
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
        /*  zz = ( crd->zf[(short)face] ); [ not required ] */

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
/* approximately TEM mode [ in stripline ] */
/*
            ux = ZERO;
            uy = uu;
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
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 24, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define EVL_TDITER 32768 /* number of iterations [ cycles ] - time domain */
# define EVL_TDINIT     1 /* first evaluated cycle - time domain */
# define EVL_TDSTOP 32768 /* last computed Maxwfld cycle - time domain */
# define EVL_TDRPTR     2 /* cycle length [ repetition rate ] - time domain */
# define EVL_FRITER 12000 /* number of interations [ cycles ] - freq domain */
# define EVL_FRINIT    21 /* first evaluated cycle - frequency domain */
# define EVL_FRSTOP   200 /* last computed Maxwfld cycle - frequency domain */
# define EVL_FRRPTR    20 /* cycle length [ repetition rate ] - freq domain */
# define EVL_HTINIT     1 /* first evaluated cycle - heat_&_fluids */
# define EVL_HTSTOP 12000 /* last computed cycle - heat_&_fluids */
# define EVL_HTRPTR    10 /* cycle length [ repetition rate ] - heat_&_fluids */
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
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0.]                                                            *
*                                                                              *
*   This function determines special [ i.e. DSC model dependent ] evaluation   *
*   modes for files dsc.val<n>                                                 *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 22, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define MDV_PLOT_SWR 1 /*1: create standing waves plots at frequencies frq[j]*/
# define MDV_SPCLPRTS 3 /* number of evalt'd special ports: 0<= n <MAX_PERIOD */
# define MDV_AUTOEVAL 0 /*1: automatic port labeling for standing wave eval.  */
# define MDV_SPLINTPL 1 /*>0: spline interpolate spc.format [references etc.] */
# define MDV_GPHINTPL 5 /*>0: spline interpl. gph.format [gnu.graphics,e.g.]  */
# define MDV_NORMALZE 0 /*1: normalize S-parmts to abs.value 1 [ cautious ! ] */
/*----------------------------------------------------------------------------*/
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
/* configure function pstprc(*): */
/*...........................................................................*/
# define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X' ] */
/*...........................................................................*/
# define PST_MESH  0      /* 1: plot mesh cell system */
# define PST_WALLS 2      /* 1: plot E-walls, 2: plot E- and M-walls */
# define PST_ISOSCLE 1    /* 1: equally scaled axes */
/*...........................................................................*/
# define PST_TMPPLT  4    /* Disable/enable temperature plot mode [0/1/...] */
# define PST_PRSPLT  3    /* Disable/enable pressure plot mode [0/1/...] */
# define PST_FLWPLT  4    /* Disable/enable velocity plot mode [0/1/...] */
# define PST_PRSINIT 1
/*...........................................................................*/
# if CMPRSSBL == 0
/* Disable density plot mode */
   # define PST_DNSPLT  0
   # define PST_DNSINIT 0
# else /* if CMPRSSBL != 0 */
/* Disable/enable density plot mode [0/1/...] */
   # define PST_DNSPLT  3
   # define PST_DNSINIT 1  
# endif /* CMPRSSBL != 0 */
/*...........................................................................*/
/* scale, vectors: */
/* PST_SCALE1: ratio of PST_SCLFLW or, if this is ZERO, of max | flow | to */
/* the geometric mesh extension */
# define PST_SCALE1 ( 1.000e-01 )

/* scale, vector pointers: */ 
/* PST_SCALE2: the ratio of the vector pointer size to the mesh extension */
# define PST_SCALE2 ( 1.000e-02 )
/*...........................................................................*/
# ifndef PST_SCLFLW /* scale of fluid flow [ maximum fluid velocity , e.g. ] */
# if MOD_DEFLT == 1
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 2
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 3
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 4
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 5
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 6
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 7
   # undef PST_PAUSE
   # define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X'] */
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 8
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 9
   # undef PST_PAUSE
   # define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X'] */
   # define PST_SCLFLW ( 1.000e-01 )
# elif MOD_DEFLT == 10
   # undef PST_PAUSE
   # define PST_PAUSE ( -1 ) /* pause/seconds [ in gnuplot command 'pause X'] */
   # define PST_SCLFLW ( 1.000e-01 )
# else
   # define PST_SCLFLW ZERO
# endif /* MOD_DEFLT == ... */
# endif /* not defined PST_SCLFLW */
/*...........................................................................*/
# include "./tools/pstprc.fld"
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

   static const char
     *timefrm = "created: %.24s ";

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

   static const char 
     *spcpfx = "spc.",
     *option = "forward", /* Fourier transformation option */
     *spformat = "%s\n",
     *ipformat = "%ld\n",
     *scformat = "%80s";

   static char 
      ptr[STS_SIZE] = {null},
      tmeptr[STS_SIZE] = {null},
      spcfle[STS_SIZE] = {null},
      dpformat[STS_SIZE] = {null},
      absc_unit[STS_SIZE] = {null},
      ordn_unit[STS_SIZE] = {null},
      t_type[MAX_PERIOD][STS_SIZE] = {{null}},
      t_text[MAX_PERIOD][STS_SIZE] = {{null}},
      s_type[MAX_PERIOD][STS_SIZE] = {{null}},
      s_text[MAX_PERIOD][STS_SIZE] = {{null}},
      fleptr[MAX_PERIOD][STS_SIZE] = {{null}},
    **endp = null;

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
      strncpy( dpformat ,"%+.15E%s", TEN );
   else
      strncpy( dpformat , "%+.15e%s", TEN );

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
               printf( " >----> enter name of file to be created "
                  "for index %-6d..........: ", idx[ii] );
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

         strncpy( t_type[h_], ( vpt->name ), STS_SIZE );
         strncpy( t_text[h_], ( vpt->text ), STS_SIZE );

         timefle[h_] = fopen( fleptr[h_], "w+" );

         fprintf( timefle[h_], spformat, t_type[h_] );
         fprintf( timefle[h_], "%s%s%d%s\n", t_text[h_], "_[evaluated"
            "_index", idx[hh], "]" );
         fprintf( timefle[h_], spformat, absc_unit );
         fprintf( timefle[h_], spformat, ordn_unit );

         fprintf( timefle[h_], dpformat, ( fpt->t[null] ), "\n" );
         fprintf( timefle[h_], dpformat, ( fpt->tt[null] ), "\n" );
         fprintf( timefle[h_], dpformat, ( fpt->dt[null] ), "\n" );
         fprintf( timefle[h_], ipformat, ( fpt->ttlg[null] ));

         for ( i_=null ; i_<( fpt->ttlg[null] ); i_++ )
         {
            fprintf( timefle[h_], dpformat, ( fpt->r[hh][i_] ), "  " );
            fprintf( timefle[h_], dpformat, ( fpt->i[hh][i_] ), "\n" );

# if LIST_FILE == 1
            printf( "  %ld:   %1.16e  + ( %1.16e )*j \n",
                    i_, ( fpt->r[hh][i_] ), ( fpt->i[hh][i_] ));
# endif
         };

         nseconds = time( timer );
         strcpy( tmeptr, ctime( &nseconds ));
         fprintf( timefle[h_], "\n%s %s %s %.24s\n", "DSC time response file",
            fleptr[h_], "terminated:", tmeptr );
         printf( " %s %s %s %.24s", "DSC time response file",
            fleptr[h_], "terminated:", tmeptr );

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
               fscanf( evalfle, scformat, ptr );

               if ( null != strncmp( ptr, ( vpt->name ), THREE ))
               {
                  printf( "\n\n Error message from function %s :", __func__ );
                  printf( "\n\n Incompatible system identifier '%s'", ptr );
                  printf( "\n in reference spectrum, file %s !!!", spcfle );
                  printf( " [ overriding.]\n" );

                  fclose( evalfle );

                  goto frq_interval;
               };
       
               fscanf( evalfle, scformat, ptr );
               fscanf( evalfle, scformat, absc_unit );
               fscanf( evalfle, scformat, ordn_unit );

/* the domain lower bound: */
               fscanf( evalfle, scformat, ptr );
               xlower = strtod( ptr, endp ); /* the lower frequency bound */
	        
/* the domain upper bound: */
               fscanf( evalfle, scformat, ptr );
               xupper = strtod( ptr, endp ); /* the upper frequency bound */

/* the increment: */
               fscanf( evalfle, scformat, ptr );
/*............................................................................*/
# if MDV_SPLINTPL == 0
               df = strtod( ptr, endp );
# else
               dx = strtod( ptr, endp );
# endif
/*............................................................................*/
/* the number of sample points: */

               fscanf( evalfle, scformat, ptr );
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
                  fscanf( evalfle, scformat, ptr ); /* read the frequency */
# endif
/*............................................................................*/
                  fscanf( evalfle, scformat, ptr ); /* the real part */
                  rref[jj] = strtod( ptr, endp );
                  fscanf( evalfle, scformat, ptr ); /* the imaginary part */
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

            strncpy( ordn_unit, ( vpt->yunit ), STS_SIZE );
            strcat( ordn_unit, "*" );
            strncat( ordn_unit, ( vpt->xunit ), STS_SIZE );
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

            strncpy( s_type[h_], ( vpt->name ), STS_SIZE );
            strncpy( s_text[h_], ( vpt->text ), STS_SIZE );

            specfle[h_] = fopen( fleptr[h_], "w+" );

            fprintf( specfle[h_], spformat, s_type[h_] );
            fprintf( specfle[h_], "%s%s%d%s\n", s_text[h_],"_[evaluated"
               "_index", idx[hh], "]" );
            fprintf( specfle[h_], spformat, absc_unit );
            fprintf( specfle[h_], spformat, ordn_unit );
            fprintf( specfle[h_], dpformat, x1, "\n" );
            fprintf( specfle[h_], dpformat, x2, "\n" );
            fprintf( specfle[h_], dpformat, dx, "\n" );
            fprintf( specfle[h_], ipformat, kk );

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
               fprintf( specfle[h_], dpformat, frequency, "  " );
	       frequency += dx;
# endif
               fprintf( specfle[h_], dpformat, ( fpt->r[hh][jj] ), "  " );
               fprintf( specfle[h_], dpformat, ( fpt->i[hh][jj] ), "\n" );

# if LIST_FILE == 1
               printf( "  %d:   %1.16le  + ( %1.16le )*j \n",
                  jj, ( fpt->r[hh][jj] ), ( fpt->i[hh][jj] ));
# endif
               jj++ ;
               ii++ ;
            };

            nseconds = time( timer );
            strcpy( tmeptr, ctime( &nseconds ));
            fprintf( specfle[h_], "\n%s %s %s %.24s\n", "DSC spectrum file",
               fleptr[h_], "terminated:", tmeptr );

            printf( "\n\n %s %s %s %.24s", "DSC spectrum file",
               fleptr[h_], "terminated:", tmeptr );

            fclose( specfle[h_] );

            hh++ ;
         };
         break;
      };/* end if lbl == THREE */
/*............................................................................*/
/* spectral response saved */
/*............................................................................*/
/* still case 6 [ cont'd ] - spectral response; spline interpolation:         */
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
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many supp. points in frequency domain !!!" );
         fprintf( stderr, "\n [ Number %ld exceeds maximum number of spline "
            "interpolation points %ld", kk, (long) SPL_SPNTS );
         fprintf( stderr, "\n   = macro SPL_INTPL in function %s.", "spline(*)" );
         fprintf( stderr, "\n   - Change macro in compliance with memory "
            "resources.]\n" );

         exit( EXIT_FAILURE );
      }
      else if ( SPL_INTPL < frqlbl )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many number of frequency points !!!" );
         fprintf( stderr, "\n [ Number %ld exceeds maximum number of spline "
            "interpolation points %ld", ( long ) frqlbl, (long) SPL_INTPL );
         fprintf( stderr, "\n   = macro SPL_INTPL in function %s.", "spline(*)" );
         fprintf( stderr, "\n   - Change macro in compliance with memory "
            "resources.]\n" );

         exit( EXIT_FAILURE );
      }
      else if ( EVL_SPF < frqlbl )
      {
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Too many number of frequency points !!!" );
         fprintf( stderr, "\n [ Number %ld exceeds maximum number %ld",
            ( long ) frqlbl, (long) EVL_SPF );
         fprintf( stderr, "\n   = macro EVL_SPF in function %s.", __func__ );
         fprintf( stderr, "\n   - Change macro in compliance with memory "
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

         fprintf( evalfle, spformat, ( vpt->name ));
         fprintf( evalfle, spformat, "reference_spectrum" );
         fprintf( evalfle, spformat, absc_unit );
         fprintf( evalfle, spformat, ordn_unit );
         fprintf( evalfle, dpformat, xlower, "\n" );
         fprintf( evalfle, dpformat, xupper, "\n" );
         fprintf( evalfle, dpformat, dx, "\n" );
         fprintf( evalfle, ipformat, frqlbl );

# if EVL_WRTFREQ == 1
	 frequency = xlower;
# endif
         j_ = null;
         while ( j_ < frqlbl )
         {
# if EVL_WRTFREQ == 1
            fprintf( evalfle, dpformat, frequency, "  " );
	    frequency += dx;
# endif
            fprintf( evalfle, dpformat, rspc[null][j_], "  " );
            fprintf( evalfle, dpformat, ispc[null][j_], "\n" );
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
         printf( timefrm, tmeptr );

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
              
         fprintf( evalfle, spformat, ( vpt->name ));
         fprintf( evalfle, spformat, ptr );
         fprintf( evalfle, spformat, absc_unit );
         fprintf( evalfle, spformat, "---" );
         fprintf( evalfle, dpformat, xlower, "\n" );
         fprintf( evalfle, dpformat, xupper, "\n" );
         fprintf( evalfle, dpformat, dx, "\n" );
         fprintf( evalfle, ipformat, frqlbl );
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
            fprintf( evalfle, dpformat, frequency, "  " );
	    frequency += dx;
# endif
/*............................................................................*/
            fprintf( evalfle, dpformat, rcpl[h_][j_], "  " );
            fprintf( evalfle, dpformat, icpl[h_][j_], "\n" );
            j_++ ;
         };

         nseconds = time(timer);
         strcpy( tmeptr, ctime( &nseconds ));

         fprintf( evalfle, "\n# SPECTRUM file %s\n", spcfle );
         fprintf( evalfle, "# created:%s", tmeptr );

         fclose( evalfle );

         printf( "\n" );
         printf( CLEAR_LINE );
         printf( "\r SPECTR.file %s ", spcfle );
         printf( timefrm, tmeptr );

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
/* reflexion, VSWR,..., insertion loss etc. have been saved */
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
/*............................................................................*/

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
	 
         fprintf( evalfle, dpformat, ( double )( fpt->t[null] ), "\n" );
         fprintf( evalfle, dpformat, ( double )( fpt->tt[null] ), "\n" );
         fprintf( evalfle, dpformat, ( double )( fpt->dt[null] ), "\n" );
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
            fprintf( evalfle, dpformat, frequency, "  " );
	    frequency += ( fpt->dt[null] );
# endif
/*............................................................................*/
            fprintf( evalfle, dpformat, ( fpt->r[ONE][jj] ), "  " );
            fprintf( evalfle, dpformat, ( fpt->i[ONE][jj] ), "\n" );
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

            fscanf( evalfle, scformat, ptr );

            if ( null != strncmp( ptr, vpt->name, THREE ))
            {
               printf( "\n\n Error message from function %s :", __func__ );
               printf( "\n\n Incompatible system identifier '%s'", ptr );
               printf( "\n on reference spectrum, file %s !!!", spcfle );
               printf( " [ overriding. ]\n" );

               fclose( evalfle );

               break;
            };

            fscanf( evalfle, scformat, ptr );
            fscanf( evalfle, scformat, absc_unit );
            fscanf( evalfle, scformat, ordn_unit );
            fscanf( evalfle, scformat, ptr );/* xlower */
            xlower = strtod( ptr, endp );

            fscanf( evalfle, scformat, ptr );/* xupper */
            xupper = strtod( ptr, endp );

            fscanf( evalfle, scformat, ptr ); /* dx */
            dx = strtod( ptr, endp );

            fscanf( evalfle, scformat, ptr );
            ii = strtol( ptr, endp, DEC );

            xx = xlower;
            jj = null;
            while( jj < ii )
            {
/*............................................................................*/
# if EVL_WRTFREQ == 1
               fscanf( evalfle, scformat, ptr );
# endif
/*............................................................................*/
               fscanf( evalfle, scformat, ptr );
               rref[jj] = strtod( ptr, endp );
               fscanf( evalfle, scformat, ptr );
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
               fprintf( evalfle, dpformat, xlower, "\n" );
               fprintf( evalfle, dpformat, xupper, "\n" );
               fprintf( evalfle, dpformat, dx, "\n" );
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
                  fprintf( evalfle, dpformat, frequency, "  " );
	          frequency += dx;
# endif
/*............................................................................*/
                  fprintf( evalfle, dpformat, rcpl[h_][jj], "  " );
                  fprintf( evalfle, dpformat, icpl[h_][jj], "\n" );
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
/*=================== end of function body modval(*) =========================*/
# undef LARGE_LOG_VAL
/************************ end of function modval(*) ***************************/


/*----------------------------------------------------------------------------*/
# endif /* [ end of section compiled with option -D_POST_MODEL ] */
/*----------------------------------------------------------------------------*/
# undef EXC_MAXW
# undef EXC_HCRR
/*************************** end of file model.c ******************************/
