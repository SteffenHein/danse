/* [ file: formtp.h ] */
/*----------------------------------------------------------------------------*/
/* typedef header of program former.c */
/* [ Update: October 23, 2007 ] */
/*----------------------------------------------------------------------------*/
# ifndef DSC_ADJCLL
   # define DSC_ADJCLL 1
# endif
/*----------------------------------------------------------------------------*/
/* mesh topology defining structure [ cf. header 'ecrrtp.h'] */

typedef struct
{ 
   char
      rtn;

   char 
      name[STS_SIZE],
      text[STS_SIZE];

   long
      mi, mf, /* initial and final mesh cell indices */
      ci, cf, /* initial and final vertex point indices     */
/*............................................................................*/
# if NFCES != 0
      fi, ff, /* initial and final face indices */
      nf[NFCES+ONE][TWO], /* nf[I][0/1]: cells adj to face I */
# endif
/*............................................................................*/
# if DSC_ADJCLL == 0
      mn[NODES+ONE][FACES], /* cn[C][F]: cell adj to cell C, face F  */
# elif DSC_ADJCLL == 1
      mn[NODES+ONE][PORTS], /* cn[C][P]: cell adj to cell C, port P  */
# endif
/*............................................................................*/
      cm[NODES+ONE][CRNRS]; /* vertex point identifier               */
                   
   signed char
      fn[NODES+ONE][FACES], /* fn[C][F]: face adjacent to cell C, face F */
      pn[NODES+ONE][PORTS]; /* connected [Maxwell field] port identifier */

} TOPOLOGY;
/*----------------------------------------------------------------------------*/
/* the corner point coordinates structure type: */
struct coordinates
{
   signed char
      rtn;
/*                                          
   char 
      name[STS_SIZE],
      text[STS_SIZE];
*/
   double
      c[CPNTS+ONE][DIMNS]; /* c[][0]=x, c[][1]=y, c[][2]=z */
};
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
      type[STS_SIZE];
} MEDIUM;
/*----------------------------------------------------------------------------*/
/* the material constants structure type: */
struct media
{            
   signed char
      rtn;

   short
      idx[NODES+ONE];     /* media identifier [ label associated to cell ] */

/* E/H field relevant media parameters: */ 
   double
      ep[MEDIA+ONE][SIX], /* rel. permittivity tensor */
      my[MEDIA+ONE][SIX], /* rel. permeability tensor */
      ke[MEDIA+ONE][SIX], /* electric conductivity tns. [A/(V*m)] */
      km[MEDIA+ONE][SIX]; /* magnetic conductivity tns. [V/(A*m)] */
/*............................................................................*/
# if DSC_HCRMDE != 0
/* thermal media parameters: */

   double
      cv[MEDIA+ONE], /* specific heat [J/(K*m^3)] */
      kh[MEDIA+ONE]; /* heat current conductivity [W/(K*m)] */
/*............................................................................*/
# if DSC_FLDMDE != 0
/* fluid media parameters: */

   double
      rm[MEDIA+ONE], /* mean mass density [Kg/m^3] */
      tm[MEDIA+ONE], /* mean temperature [C] */
      bm[MEDIA+ONE], /* mean thermal expansion coefficient [1/K] */
      cm[MEDIA+ONE], /* mean [ adiabatic ] compression coefficient [1/Pa] */
      ny[MEDIA+ONE], /* dynamic viscosity [Kg/(sec*m)] */
      q1[MEDIA+ONE], /* Cp/Cv - 1 [dimensionless] */
      td[MEDIA+ONE], /* [dissipation] time constant [sec] */
      LL[MEDIA+ONE]; /* characteristic length [in Prandtl turb model, e.g.]*/

   double
      gr[MEDIA+ONE][THREE], /* gravitation acceleration [ m/(sec^2) ]*/
      gp[MEDIA+ONE][THREE]; /* pressure gradient [ Pa/m ]*/

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* Gyrotropic media parameters: */

   double
      mi[MEDIA+ONE][DIMNS], /* int. magn. flux, plasma  [Tesla] */
      ms[MEDIA+ONE][DIMNS], /* saturation magnetization [Tesla] */
      hg[MEDIA+ONE][DIMNS], /* internal magnetic field    [A/m] */
      no[MEDIA+ONE],        /* plasma electron density  [1/m^3] */
      tr[MEDIA+ONE],        /* plasma curr.relax.time [seconds] */
      tg[MEDIA+ONE],        /* spin-spin  relax.time  [seconds] */
      ld[MEDIA+ONE];        /* LANDE factor */

   char
      type[MEDIA+ONE][STS_SIZE];
};
/*----------------------------------------------------------------------------*/
/* repeater function structure type */
/* [ stores repeatedly used identical s-parameters ]: */
struct repeater
{
   long
      smx[NODES+ONE];
};
/*----------------------------------------------------------------------------*/
# if ( DSC_HCRMDE != 0 )\
   &&( DSC_FLDMDE != 0 )

struct fconnect
{
   short
      cnn[NODES+ONE]; /* 0<N: cnn[N] = fluid connected component of node N */
                      /* cnn[null] = counter of fluid connected components */
   double
      crsdt[NFCNN+ONE]; /* coarsening time period */
};
# endif
/*----------------------------------------------------------------------------*/
/* parameter reference structure type: */
typedef struct
{
   signed char
      rtn;

   char
      name[STS_SIZE],
      text[STS_SIZE];

   char
      domain[STS_SIZE];

   double
      fr;

   struct coordinates
     *cpt;

   struct media
     *mpt;

   struct repeater
     *rep;
/*............................................................................*/
# if DSC_HCRMDE != 0
   struct repeater 
     *rhp; /* -> thermal & fluid s-parameter repeater */
/*............................................................................*/
# if DSC_FLDMDE != 0
   struct fconnect
     *fcp; /* -> fluid connected components structure */
# endif
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
} PARAMETERS;
/*----------------------------------------------------------------------------*/
/* the boundary structure type: */
typedef struct 
{
   signed char
      rtn;

   char
      name[STS_SIZE],
      text[STS_SIZE];

   long 
      n, p, 
      m[BNDAP], cb[BNDPR], cp[BNDPR];

   unsigned char
      f[BNDAP], p0[BNDPR], p1[BNDPR];

   signed char 
      pp0[BNDPR],pp1[BNDPR];

   double
      phi,

      r00[BNDAP], r01[BNDAP],
      r10[BNDAP], r11[BNDAP],

      i00[BNDAP], i01[BNDAP],
      i10[BNDAP], i11[BNDAP];
/*............................................................................*/
# if DSC_HCRMDE != 0
   long
      ntf, /* number of temperature boundary faces */
      ntn, /* number of temperature nodes */
      nhc, /* the number of heat current boundary faces */
      nsc, /* number of surface heat conducting faces */
      nsk; /* number of [ skin effect ] heat source bndry nodes */

   long           /* pertinent mesh cell indices */            
      mtf[BNDTF], /* temperature boundary cell index */
      mtn[BNDTN], /* temperature node [i.e. cell] index */
      mhc[BNDHC], /* heat current boundary cell index */
      msc[BNDSC], /* surface heat conducting cell index */
      msk[BNDSK]; /* [ skin effect ] heat source cell index */

   signed char    /* pertinent face indices */
      ftf[BNDTF], /* temperature boundary face index */
      fhc[BNDHC], /* heat current boundary face index */
      fsc[BNDSC], /* surface heat conducting face index */
      fsk[BNDSK]; /* [ skin effect ] heat source face index */

   double        /* pertinent real parameters */
      tf[BNDTF], /* value of face temperature */
      tn[BNDTN], /* value of node temperature */
      hc[BNDHC], /* imposed incoming heat current [per face] */
      sc[BNDSC], /* surface heat conductance [W/(m^2*K)]*/
      tr[BNDSC], /* reference [e.g.envrnmnt] temperature */
      rd[BNDHC], /* outgoing heat radiation coeff [hc -= t4*T^4] */
      rs[BNDSK], /* skin effect surface resistivity [ R_square x */
      sk[BNDSK][TWO][TWO]; /* x penetration depth, Ohm*m ] */
                           /* skin effect heat source coefficients */
/*............................................................................*/
# if DSC_FLDMDE != 0
   long
      nns, /* number of no-slip boundary faces */
      nsl, /* number of free slip boundary faces */
      nif, /* number of inflow boundary faces */
      nof; /* number of outflow boundary faces */

   long /* pertinent cells [ mesh cell index ] */                      
      mns[BNDNS],
      msl[BNDSL],
      mif[BNDIF],
      mof[BNDOF];

   signed char /* pertinent faces */
      fns[BNDNS],
      fsl[BNDSL],
      fif[BNDIF],
      fof[BNDOF];

   double
      ti[BNDIF], /* inflow temperature */
      uf[BNDIF][THREE], /* inflow fluid velocity vector [m/s] */
      vf[BNDOF][THREE], /* outflow fluid velocity vector [m/s] */
      nf[BNDSL][THREE]; /* face normal vector */
/*............................................................................*/
# if NUSSELT != 0
   double
      nus[BNDNS];       /* Boundary layer Nusselt number */
# endif /* NUSSELT != 0 */
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
} BOUNDARIES;
/*----------------------------------------------------------------------------*/
/* the excitation mode structure type */
typedef struct  
{  
   signed char
      rtn;

   char
      name[STS_SIZE],
      text[STS_SIZE];

   char
      type[STS_SIZE];

   short       
      ne, nh, nn;

   long       
      me[EXCEP], mh[EXCHP];

   signed char 
      pe[EXCEP], ph[EXCHP];

   double     
      rt, ht, dt,            /* rise_time, hold_time, delay_time */
      fr[EXCFR],             /* max. number of exciting frequencies */
      er[EXCEP], ei[EXCEP],
      hr[EXCHP], hi[EXCHP];
/*............................................................................*/
# if DSC_HCRMDE != 0 /* excited thremal nodes and faces */
   char
      hctp[STS_SIZE];

   short
      hcn; /* any purpose integer [ such as smoothing order, e.g.] */

   double
      hcrt, /* rise time */
      hcht, /* hold time */
      hcdt, /* delay time */
      hcrt2, /* 2nd rise time [ option DOUBLE_PERIODIC____ ] */
      hcht2, /* 2nd hold time */
      hcdt2; /* 2nd delay time */
    
   short
      nhc, /* number of excited heat currents */
      ntf, /* number of excited temperature faces ] */
      ntn; /* number of excited temperature nodes */

   long
      mhc[EXCHC], /* excited cell index [ heat currents ] */
      mtf[EXCTF], /* excited cell index [ face temperatures ] */
      mtn[EXCTN]; /* excited cell index [ node temperatures ] */

   short 
      fhc[EXCHC], /* excited face index [ heat currents ] */
      ftf[EXCTF]; /* excited face index [ face temperatures ] */

   double 
      hs[EXCHC], /* nodal heat source density */
      hc[EXCHC], /* excited face current */
      tf[EXCTF], /* excited face temperature */
      tn[EXCTN]; /* excited node temperature */

# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
} EXCITATION;
/*----------------------------------------------------------------------------*/
/* the evaluation mode structure type: */
typedef struct  
{
   signed char
      rtn;

   char
      name[STS_SIZE],
      text[STS_SIZE];

   short       
      r, /* internal repetition rate */
      nep, /* number of evaluated E ports */
      nhp, /* number of evaluated H ports */
      nen, /* number of evaluated E nodes */
      nhn; /* number of evaluated H nodes */

   long   
      n,          /* number of connection-reflection cycles */
      ni,         /* 1st evaluated cycle */
      nf;         /* last computed cycle [ Maxwell field ] */

   long   
      mep[EVLEP], /* evaluated E port cell indices */
      mhp[EVLHP], /* evaluated H port cell indices */
      men[EVLEN], /* evaluated E node cell indices */
      mhn[EVLHN]; /* evalueted H node cell indices */

   signed char 
      pep[EVLEP], /* evaluated E port [port] indices */
      php[EVLHP], /* evaluated H port indices */
      cen[EVLEN], /* evaluated E node indices */
      chn[EVLHN]; /* evaluated H node indices */

   char 
      mode_ep[SHS_SIZE],
      mode_hp[SHS_SIZE],
      mode_en[SHS_SIZE],
      mode_hn[SHS_SIZE];
/*............................................................................*/
# if DSC_HCRMDE != 0 /* evaluated thermal nodes and faces */

   short  
      rc;         /* internal repetition rate [ heat and fluids ]*/

   short  
      nhc,        /* number of evaluated heat currents */
      ntf,        /* number of evaluated face temperatures */
      ntn;        /* number of evaluated node temperatures */

   long
      nj,         /* 1st evaluated cycle [ heat and fluids ] */
      nt;         /* last computed cycle [ heat and fluids ] */

   long
      mhc[EVLHC], /* evaluated cell index [ heat currents ] */
      mtf[EVLTF], /* evaluated cell index [ face temperatures ] */
      mtn[EVLTN]; /* evaluated cell index [ node temperatures ] */

   signed char
      fhc[EVLHC], /* evaluated face index [ heat currents ] */
      ftf[EVLTF]; /* evaluated face index [ face temperatures ] */
/*............................................................................*/
# if DSC_FLDMDE != 0 /* evaluated fluids */
   short  
      nun;        /* number of evaluated nodal velocities */

   long
      mun[EVLTN]; /* mesh cell index [ nodal velocities ] */

   signed char
      cun[EVLTF]; /* component index [ nodal velocities ] */

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
} EVALUATION;
/*----------------------------------------------------------------------------*/
/* driver function formdrv(*) state structure */
/* [ reflects actually charged topology, parameter, boundary conditions, ..., */
/*   file names, file labels etc.]: */
typedef struct
{
   signed char
      rtn;

   char
      domain[STS_SIZE], /* domain indicator */
      file[STS_SIZE],
      flbl[FOUR]; /* index n of files elf.top<n>,..., elf.val<n> */

   short
      item,
      stat;

   double
      dt, dp, fr; /* dt: Maxwell field DSC time step [ time domain ] */
/*............................................................................*/
# if DSC_HCRMDE != 0
   double
      hcrstart; /* start thermal computations with that delay [ seconds ] */

   double
      hcdt, hcdp; /* heat current time step/phase shift */
/*............................................................................*/
# if DSC_FLDMDE != 0
   double
      fldstart; /* start fluid dynamic computations with that delay [ sec ]*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/

   TOPOLOGY
     *tpt;

   PARAMETERS
     *ppt;

   BOUNDARIES
     *bpt;

   EXCITATION
     *ept;

   EVALUATION
     *vpt;

} FORMSTATE;
/*************************** end of file formtp.h *****************************/
