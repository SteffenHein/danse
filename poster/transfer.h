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
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
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
