/* [ file: solvtp.h ] */
/* Update: 04 November 2007 */
/*----------------------------------------------------------------------------*/
# ifndef DSC_ADJCLL
   # define DSC_ADJCLL 0
# endif
/*----------------------------------------------------------------------------*/
/* The mesh topology structure [ essentially assigning adjacent cells, faces, */
/* and ports to cells faces and ports ] */
/*............................................................................*/
struct topology
{
   signed char
      rtn;

   char
      name[SHS_SIZE], /* model identifier */
      text[STS_SIZE]; /* any comment */

   long
      n; /* the total number of cells [i.e. nodes] */
/*............................................................................*/
# if DSC_ADJCLL == 0 /* <- default; recommended [more economical] usage */
   long
      m[NODES+ONE][FACES];    /* the cell adjacent to face [2nd argument] */
# else /* DSC_ADJCLL == 1 */ /* of cell [1st argument] */
   long
      m[NODES+ONE][PORTS];    /* the cell adjacent to port [2nd argument] */
# endif                       /* of cell [1st argument] */
/*............................................................................*/
   signed char
      f[NODES+ONE][FACES],    /* the face adjacent to face [2nd argument] */
      p[NODES+ONE][PORTS];    /* the port adjacent to port [2nd argument] */
                              /* of cell [1st argument] */
}; /* [ end of struct topology ] */
/*----------------------------------------------------------------------------*/
/* Maxwell field scattering parameters structure */
/*............................................................................*/
struct tlmsmx
{ 
   signed char
      rtn;

   char
      name[SHS_SIZE],
      text[STS_SIZE];

   char
      etyp[NRSMX+ONE], /* type of s-matrix related to Ampere's law */
      mtyp[NRSMX+ONE]; /* type of s-matrix Faraday's law */
                       /* [diagonal, symmetric, asymmetric, or trivial, e.g.] */
   double
      adm, /* Maxwell field admittance [ = sqrt( EPS_VAC/MY_VAC_) ] */
      shr,
      shi;

/* S-parameters: */

   double
      se[NRSMX+ONE][SENTR][SENTR], /* E-field s-parameters [to Ampere's law] */
      sh[NRSMX+ONE][SENTR][SENTR], /* H-field s-parameters [to Faraday's law] */

      sge[GESMX][DIMNS][DIMNS], /* gyroelectric s-parameters */
      sgm[GMSMX][DIMNS][DIMNS], /* gyromagnetic s-parameters */

/* deflection parameters: */

      de[GESMX][DIMNS][DIMNS], /* E-field */
      dm[GMSMX][DIMNS][DIMNS], /* H-field */

/* gyroelectric&magnetic exponential matrices eue&m = exp[X*Ue&m]: */

      eue[GESMX][DIMNS][DIMNS],
      eum[GMSMX][DIMNS][DIMNS];
   
/* indices */

   long
      hh[NODES+ONE], /* smx.hh[m] = s-matrix pertinent to cell indexed m */
      me[NRSMX+ONE],
      mm[NRSMX+ONE];

   long
      ge, gm; /* gyroelectric and -magnetic cell counter */

}; /* [ end of struct tlmsmx ] */
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
/* thermal scattering parameters structure */
/*............................................................................*/
struct hcrsmx
{ 
   signed char
      rtn;

/* type of heat current cell [ trivial, source, e.g.] */

   char
      ttyp[HCSMX+ONE], /* type of thermal node [ trivial, hcurr, fluid, e.g.] */
      scs[HCSMX+ONE];  /* electromagnetic heat source indicator */

   long
      hh[HCNDS+ONE]; /* hcs.hh[m] = thermal updating parameter [set] */
                     /* pertinent to cell indexed m */
   double
      adm, /* field admittance */
   /* hcdt, *//* thermal time step [ seconds ] */
   /* dt, *//* time step [ seconds ] */
      fr, /* frequency [ Hz ] */
      ph, /* phase shift [ radians ] */
      shr,
      shi;

/* other [ temporal, spatial, geometric,... ] parameters: */

   double
      ct[HCSMX+ONE], /* normalized temperature updating coefficient ( dt/cv ) */
     vol[HCSMX+ONE], /* cell volume [m^3]*/
      sf[HCSMX+ONE], /* cell surface [m^2] */
      fm[HCSMX+ONE][FACES], /* fm[][n] area of face n */
      f[HCSMX+ONE][FACES][THREE]; /* f[][n][] exterior face vector of face n */

/* material parameters: */

   double
      kh[HCSMX+ONE], /* heat conductivity [W/(K*m)]*/
      cv[HCSMX+ONE], /* heat capacity per VOLUME [J/(K*m^3) !!!] */
      ke[HCSMX+ONE], /* electric conductivity */
      km[HCSMX+ONE]; /* magnetic conductivity */
/*............................................................................*/
# if DSC_FLDMDE != 0
/* fluid updating parameters structure */

   short
      cnn[HCNDS+ONE]; /* hcs.cnn[N] = fluid connected component, cell N [>0] */
                      /* hcs.cnn[null] counts the fluid connected components */
   signed char
      fctype[HCNDS+ONE][FACES]; /* fluid [ boundary ] face type identifier */
                                /* hcs.fctype[m][f]: type of face f in cell m */
   double
      ft[HCSMX+ONE], /* normalized flow updating coefficient (dt/rm) */
      rm[HCSMX+ONE], /* mean density [Kg/(m^3)] */
      tm[HCSMX+ONE], /* mean temperature [Celsius] */
      bm[HCSMX+ONE], /* mean expansion coefficient [1/K] */
      cm[HCSMX+ONE], /* mean [ adiabatic ] compression coefficient [1/Pa] */
      ny[HCSMX+ONE], /* dynamic viscosity [Kg/(sec*m)] */
      q1[HCSMX+ONE], /* Cp/Cv - 1 [dimensionless] */
      td[HCSMX+ONE], /* dissipation time [sec] */
      dc[HCSMX+ONE]; /* dissipation factor [1-hcdt/td; dimensionless] */
/*............................................................................*/
/* Prandtl turbulence models */
# if (( TURBMOD == 1 ) \
    ||( TURBMOD == 2 ))
   double
      rl2[HCSMX+ONE]; /* square of characteristic length; rl2 = rm*LL^2; */
# endif                                         /* LL: cf. struct media */
/*............................................................................*/
   double
      gr[HCSMX+ONE][THREE], /* gravitational acceleration [ m/sec^2 ] */
      gp[HCSMX+ONE][THREE]; /* pressure gradient [ Pa/m ] */

   double
      tmean[NFCNN+ONE],  /* mean temperature of fluid connected components */
      crsdt[NFCNN+ONE],  /* time interval at that coarsening takes place */
      eptlt[NFCNN+ONE],  /* total potential energy in fluid cnn components */
      volume[NFCNN+ONE]; /* total volume of fluid connected components */
/*............................................................................*/
# if CMPRSSBL != 0
   double
      mass[NFCNN+ONE], /* total mass in fluid connected components */
      msdr[NFCNN+ONE]; /* mass [ defect ] ratio in fluid connected components */
# endif /* CMPRSSBL != 0 */
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
   double
       s[HCSMX+ONE][FACES][DIMNS], /* face vectors ( F^t ) * adj( B^-1 ) */
      bi[HCSMX+ONE][DIMNS][DIMNS]; /* matrix adj ( B^-1 ) */

/* double
       b[HCSMX+ONE][DIMNS][DIMNS], *//* matrix adj( B ) */
    /* pb[HCSMX+ONE][DIMNS][DIMNS][DIMNS]; */
                                /* pb[ii][jj][kk] = adj(B^-1)*<b[ii]|...>[jj] */
/* deflection parameters: */
/*
      dc[HCSMX+ONE][DIMNS][DIMNS];
*/
}; /* [ end of struct hcrsmx ] */
/*----------------------------------------------------------------------------*/
# endif /* DSC_HCRMDE != 0 */
/*----------------------------------------------------------------------------*/
struct boundary
{
   signed char
      rtn;

   char
      name[SHS_SIZE],
      text[STS_SIZE];

   double
      r, i,
      r00[BNDAP], r01[BNDAP], /* boundary face reflection matrix, real part */
      r10[BNDAP], r11[BNDAP], /* [ time or frequency domain ] */
      i00[BNDAP], i01[BNDAP], /* boundary face reflection matrix, imag.part */
      i10[BNDAP], i11[BNDAP]; /* [ frequency domain ] */
                   
   long
      n, /* number of aperiodic boundary faces */
      p; /* number of periodic boundary faces */

   long
      m[BNDAP], cb[BNDPR], cp[BNDPR]; /* mesh cell indices for aperiodic [m] */
                                      /* and periodic boundary faces [cb, cp] */
   signed char
      f[BNDAP], /* aperiodic boundary face indices */
     p0[BNDPR],  p1[BNDPR], /* periodic  boundary port indices */
    pp0[BNDPR], pp1[BNDPR]; /* pertinent periodic port indices */
/*............................................................................*/
# if DSC_HCRMDE != 0
/* numbers: */

   long
      ntf, /* number of fixed temperature boundary faces */
      ntn, /* number of fixed temperature nodes */
      nhc, /* number of faces with imposed incident heat current */ 
      nsc, /* number of surface heat conducting faces */ 
      nsk; /* number of [skin effect] heat source faces */

/* pertinent mesh cell indices: */

   long
      mtf[BNDTF], /* fixed face temperature cell indices */
      mtn[BNDTN], /* fixed node temperature cell indices */
      mhc[BNDHC], /* imposed [face-] heat current cell indices */
      msc[BNDSC], /* surface heat conducting cell indices */
      msk[BNDSK]; /* [ skin effect ] heat src bndry cell inices */

/* face labels: */

   signed char
      ftf[BNDTF], /* fixed temperature face indices */
      fhc[BNDHC], /* heat current boundary face indices */
      fsc[BNDSC], /* surface heat conducting cell faces */
      fsk[BNDSK]; /* [ skin effekt ] heat src face indices */

/* real parameters: */

   double
      tf[BNDTF],
      tn[BNDTN],
      hc[BNDHC], /* imposed incoming radiation current */
      rd[BNDHC], /* outgoing heat radiation coeff [hc += rd*Tn^4]*/
      sc[BNDSC], /* envrmnt-surface heat conductance */
      tr[BNDSC], /* reference temperature */
      hs[BNDSK], /* smoothed heat sources at skin effect faces */
      rs[BNDSK], /* skin effect surface resistivity [ R_square x */
      sk[BNDSK][TWO][TWO]; /* x penetration depth, Ohm*m ] */
                           /* skin effect heat source coefficients */
/*............................................................................*/
# if DSC_FLDMDE != 0
   signed char
      fcetype[NODES][FACES];

   long
      nns, /* number of no-slip boundary faces */
      nsl, /* number of free slip boundary faces */
      nif, /* number of inflow boundary faces */
      nof; /* number of outflow boundary faces */

   long /* pertinent cells [ mesh cell indices ] */
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
      ti[BNDIF],        /* inflow fluid temperature */
      uf[BNDIF][THREE], /* inflow fluid velocity vector [m/s] */
      vf[BNDOF][THREE], /* outflow fluid velocity vector [m/s] */
      nf[BNDSL][THREE]; /* face normal vector */
/*............................................................................*/
# if NUSSELT != 0
   double
      nus[BNDNS]; /* Nusselt number at no-slip boundary face */
# endif
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
}; /* [ end of struct boundary ] */
/*----------------------------------------------------------------------------*/
/* The excitation parameters structure */
/*............................................................................*/
struct excitation
{  
   signed char
      rtn;

   char
      name[SHS_SIZE],
      text[STS_SIZE],
      type[SHS_SIZE];

   double
      rt, ht, dt, /* rise time, hold time, delay time */
      mx, sq,     /* maximum, squared sum of exc. fields */
      er[EXCEP],  /* excited voltages, real part */
      ei[EXCEP],  /* excited voltages, imaginary part */
      hr[EXCHP],  /* excited currents, real part */
      hi[EXCHP],  /* excited currents, imaginary part */
      fr[EXCFR];  /* frequencies [ multiple harmonic excit., e.g.] */

   long
      ne, /* number of excited E ports */
      nh, /* number of excited h ports */
      nn; /* any label, e.g. number of exc.freqs. */

   long
      me[EXCEP], /* electrically exc. cell indices */
      mh[EXCHP]; /* magnetically "    "    "       */

   signed char
      pe[EXCEP], /* electrically exc. port indices */
      ph[EXCHP]; /* magnetically "    "    "       */

   signed char
      lbl; /* excitation type label [ index ] */
/*............................................................................*/
# if DSC_HCRMDE != 0
/* additional heat current nodes to be excited: */

   signed char
      hcl; /* excitation type label [index] */

   char
      hctp[SHS_SIZE]; /* thermal excitation type string */

   short
      hcn; /* any purpose integer [ such as smoothing order, e.g.] */

   double
      hcrt, /* temperature rise time */
      hcht, /* temperature hold time */
      hcdt, /* temperature delay time */
      hcrt2, /* 2nd temperature rise time [ option DOUBLE_PERIODIC____, e.g.] */
      hcht2, /* '' hold time */
      hcdt2; /* '' delay time */

   long
      nhc, /* number of excited heat currents [ incident at faces ] */
      ntf, /* number of excited face temperatures */
      ntn; /* number of excited node temperatures */

   long
      mhc[EXCHC], /* excited heat current cell index */
      mtf[EXCTF], /* excited face temperature cell index */
      mtn[EXCTN]; /* excited node temperature cell index */

   signed char 
      fhc[EXCHC], /* excited heat current face index */
      ftf[EXCTF]; /* excited temperature face index */

   double
      hs[EXCHC], /* nodal source intensity */
      hc[EXCHC], /* excited node temperature */
      tf[EXCTF], /* excited face temperature */
      tn[EXCTN]; /* excited face current */

# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
}; /* [ end of struct excitation ] */
/*----------------------------------------------------------------------------*/
/* The evaluation parameters structure */
/*............................................................................*/
struct evaluation
{
   signed char
      rtn;

   char
      name[SHS_SIZE],
      text[STS_SIZE];

   double
      dt, fr, ph; /* printing time step, operating frequency, phase shift */

   double
      epr[EVLEP], epi[EVLEP],
      hpr[EVLHP], hpi[EVLHP],
      enr[EVLEN], eni[EVLEN],
      hnr[EVLHN], hni[EVLHN];

   short
      r;  /* repetition [averaging] rate for Maxwell field computation */

   short
      nep, /* number of evaluated E-ports */
      nen, /* number of evaluated E-nodes */
      nhp, /* number of evaluated H-ports */
      nhn; /* numper of evaluated H-nodes */

   long
      n,  /* number of iteration cycles */
      ni, /* the first evaluated cycle [ Maxwell field ] */
      nf; /* the last computed cycle [ Maxwell field ] */
      
   long
      mep[EVLEP],
      mhp[EVLHP],
      men[EVLEN],
      mhn[EVLHN];

   signed char
      pep[EVLEP],
      php[EVLHP],
      cen[EVLEN],
      chn[EVLHN];

   char
      read;

   char
      mode_ep[SHS_SIZE],
      mode_hp[SHS_SIZE],
      mode_en[SHS_SIZE],
      mode_hn[SHS_SIZE];
/*............................................................................*/
# if DSC_HCRMDE != 0
/* thermal nodes to be evaluated: */

   short
      rc; /* internal repetition [ averaging ] rate for heat and fluid comput.*/

   short
      ntf[DSC_HCRMDE], /* number of evaluated temperature faces */
      ntn[DSC_HCRMDE], /* number of evaluated temperature nodes */
      nhc[DSC_HCRMDE]; /* number of evaluated current faces */

   long
      nj, /* the first evaluated cycle [ heats and fluids ] */
      nt; /* the last computed cycle [ heats and fluids ] */
/*............................................................................*/
   long
      mtf[DSC_HCRMDE][EVLTF], /* evaluated cell index [ face evaluation ] */
      mtn[DSC_HCRMDE][EVLTN], /* evaluated cell index [ node evaluation ] */
      mhc[DSC_HCRMDE][EVLHC]; /* evaluated cell index [ node evaluation ] */

   signed char
      ftf[DSC_HCRMDE][EVLTF], /* evaluated face index [ temperatures ] */
      fhc[DSC_HCRMDE][EVLHC]; /* evaluated face index [ heat currents ] */

   double
      tf[DSC_HCRMDE][EVLTF], /* evaluated face temperature */
      tn[DSC_HCRMDE][EVLTN], /* evaluated node temperature */
      hc[DSC_HCRMDE][EVLHC]; /* evaluated face current */
/*............................................................................*/
# if DSC_FLDMDE != 0
/* fluid nodes to be evaluated: */

   short
      nun[DSC_HCRMDE]; /* number of evaluated fluid velocity nodes */

   long
      mun[DSC_HCRMDE][EVLUN]; /* evaluated fluid velocity mesh cell index, */

   signed char
      cun[DSC_HCRMDE][EVLUN]; /* component index [0:x, 1:y, 2:z], */

   double
      un[DSC_HCRMDE][EVLUN]; /* evaluated fluid velocity component */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
}; /* [ end of struct evaluation ] */
/*----------------------------------------------------------------------------*/
/* The Maxwell field structure */
/*............................................................................*/
typedef struct
{
   signed char
      rtn;

   double                   /* time domain:        | frequency domain:        */
      r[NODES+ONE][PORTS],  /* reflected voltages  | voltages, real part      */
      i[NODES+ONE][PORTS];  /* incident  "         | "       , imaginary part */

} DSC_FIELDS;
/*----------------------------------------------------------------------------*/
/* The 'gauge' [ deflection ] field structure */
/*............................................................................*/
typedef struct
{
   signed char
      rtn;

   double
      r[NODES+ONE][THREE],
      i[NODES+ONE][THREE];

} DSC_GGEFLD;
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
/* The thermal&fluid fields structure */
/*............................................................................*/
typedef struct
{
   signed char
      rtn;

   double
      qn[HCNDS+ONE],        /* heat source density */
      tn[HCNDS+ONE],        /* nodal temperatures */
      tf[HCNDS+ONE][FACES]; /* face temperatures */

   double
      ic[HCNDS+ONE][FACES], /* heat currents over the faces */
      tt[HCNDS+ONE][FACES]; /* truncated and normalized temperature gradients */
                            /* on the faces */
/*............................................................................*/
# if DSC_FLDMDE != 0

   double
      pn[HCNDS+ONE],        /* nodal pressure */
      pf[HCNDS+ONE][FACES]; /* face pressure */

   double /* pressure exchange parameters */
      pg[HCNDS+ONE][FACES], /* < fc | grad( p ) > */
      pt[HCNDS+ONE][FACES]; /* truncated pressure gradients on the faces */

   double
      un[HCNDS+ONE][THREE], /* nodal fluid velocity */
      uf[HCNDS+ONE][FACES][THREE]; /* face fluid velocity */

   double
      iu[HCNDS+ONE][FACES][THREE], /* < fc | grad( uf[*] ) > */
      ut[HCNDS+ONE][FACES][THREE]; /* truncated and normalized flow gradients */
                                   /* on the faces */
/*............................................................................*/
# if CMPRSSBL != 0
/* compressible flow */

   double
      rn[HCNDS+ONE],        /* nodal fluid density */
      rf[HCNDS+ONE][FACES]; /* face fluid density */

   double
      rg[HCNDS+ONE][FACES], /* < fc | grad( r ) > */
      rt[HCNDS+ONE][FACES]; /* truncated density gradients on the faces */

# endif /* CMPRSSBL != 0 */
/*............................................................................*/
# if NUSSELT != 0

   double
      nus[HCNDS+ONE];       /* Nusselt number */

# endif /* NUSSELT != 0 */
/*............................................................................*/
# ifdef TURBMOD
/* Prandtl turbulence models */
# if (( TURBMOD == 1 ) \
    ||( TURBMOD == 2 ))
   double
      myt[HCNDS+ONE]; /* Prandtl turbulent viscosity */
# endif /* TURBMOD != 0 */
# endif /* defined TURBMOD */
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
} DSC_HCRRTS;
/*----------------------------------------------------------------------------*/
# endif /* DSC_HCRMDE != 0 */
/*----------------------------------------------------------------------------*/
/* The gyrotropic deflection field structure: */
typedef struct
{
   signed char
      rtn;

   double
      de[GENDS][DIMNS], /* gyrel. deflect.volt. */
      ge[GENDS][DIMNS], /* gyrel. current */
      dm[GMNDS][DIMNS], /* gyrmg. deflect.volt. */
      gm[GMNDS][DIMNS], /* gyrmg. current */
      emx,              /* max. gyroelel. current density */
      mmx;              /* max. gyromagn. current density */

} DSC_DEFLCT;
/*----------------------------------------------------------------------------*/
/* The cell state field structure */
/*
*//* typedef struct
*//* {
*//*    signed char
*//*       rtn;
*//*
*//*    double
*//*       b[DIMNS][DIMNS],
*//*       f[FACES][DIMNS];
*//*
*//*   double
*//*      r[PORTS],
*//*      i[PORTS];
*//*
*//*   DSC_FIELDS *fpt;
*//*   DSC_GGEFLD *gfp;
*//*
*//* # if DSC_HCRMDE != 0
*//*   double 
*//*      rc[DSC_HCRMDE][FACES],
*//*      ic[DSC_HCRMDE][FACES];
*//*
*//*   double
*//*      q[DSC_HCRMDE],
*//*      m[DSC_HCRMDE],
*//*      e[DSC_HCRMDE][THREE],
*//*      h[DSC_HCRMDE][THREE],
*//*      v[DSC_HCRMDE][THREE];
*//*
*//*   struct hcrsmx *hsp[DSC_HCRMDE];
*//*   DSC_HCRRTS *hci[DSC_HCRMDE];
*//*   DSC_HCRRTS *hre[DSC_HCRMDE];
*//*
*//* # endif
*//*
*//* } CLLSTS;
*/
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
/*............................................................................*/
/* The cell cluster structure [ embracing a set of spatial finite differences */
/* and canonical values derived from the ( non-orthogonal ) cell geometry ]  */
/* - A structure used in some frequently re-iterarated standard operations */
/*............................................................................*/
typedef struct
{
   double
      node, /* actual node value [ so long left fixed ] */
      nupd, /* new node value [ to be actually updated ] */
      dvgr, /* mean divergence of gradient [ surface integral: SUM over */
            /* the cell faces fc ( grdf[fc] | face vector[fc] ) ] */
      face[FACES], /* actual face value */
      grad[THREE], /* nodal gradient */
      db[THREE], /* nodal gradient in node vector representation */
                 /* db[j] = < b[j] | grad(*) > */
      grdf[FACES][THREE]; /* actual gradient on faces */

   long
      hh; /* the [ actual ] pertinent s-parameter index */

   unsigned char
      par; /* parameter type ['t'emperature, 'p'ressure, 'v'elocity, 'f'low, */
           /* or 'd'ensity */

   struct hcrsmx
     *hsp;
} CLUSTER;
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*----------------------------------------------------------------------------*/
struct solverstat
{
   signed char
      rtn;

   double
      swing, /* swing [ weighting factor during Maxwell field excitation ] */
      ttime, /* Maxwell field actual process time */
      dt,    /* Maxwell field DSC time step [ seconds ] */
      fr,    /* Maxwell field frequency [ Hz ] */
      ph;    /* Maxwell field phase shift [ radians ] */

   long
      nn,    /* number of [ external ] iteration loops */
      mxwnn; /* Maxwell field iterations counter */

   long
      logps, /* log file position indicator */
      stofs,
      endfle;

   short
      rstart,
      nbrjbs; /* number of jobs */

   char
      rmfle, /* 'remove files' [ operation mark ] */
      dmn,   /* time/frequency domain identifier  */
     *opt,   /* [program startup] command line options string */
     *top,
     *smx,
     *bnd,
     *exc,
     *val,
     *prfx,
     *logfle,
     *fcterr,
     *errmsg;

   struct topology 
     *tpt;

   struct tlmsmx
     *spt;

   struct boundary
     *bpt;

   struct excitation
     *ept;

   struct evaluation
     *vpt;

   DSC_FIELDS
     *inc,
     *out;
/*............................................................................*/
# if DSC_DOMAIN == 0
   DSC_FIELDS
     *gfp;

   DSC_DEFLCT
     *dfp;
# elif DSC_DOMAIN == 1
   DSC_GGEFLD
     *gfp;

   DSC_DEFLCT
     *dfp;
# endif
/*............................................................................*/
# if DSC_HCRMDE != 0

   long
      hcrnn; /* thermal iterations counter */ 

   double
      hcrstart; /* start thermal operations with that delay [ seconds ] */

   double
      hctme, /* heat current process time */
      hcdt,  /* heat current DSC time step */
      hcswg, /* heat current swing [weighting factor during excitation phase] */
      hsk;   /* heat sources [sum over skin effect lossy boundary faces, e.g.]*/

   char
      sws,   /* switched electric sources: OFF=0, ON=1 */
      hclbl; /* heat current field label [ index ] */
      
   struct hcrsmx
     *hsp;

   DSC_HCRRTS
     *hci[DSC_HCRMDE],
     *hre[DSC_HCRMDE];

   CLUSTER
     *tmp;
/*............................................................................*/
# if DSC_FLDMDE != 0

   FILE
     *dscstr, /* flow parameters log */
     *prsstr, /* pressure log */
     *tmpfst, /* temperature field */
     *prsfst, /* pressure field */
# if CMPRSSBL != 0
     *dnsfst, /* density field */
# endif
     *velfst; /* velocity field */

   char
      lesfltr;

   char
      dsclog[STS_SIZE], /* flow parameters log file name */
      prslog[STS_SIZE], /* pressure log file name */
      tmpdsp[STS_SIZE], /* temperature field file name */
      prsdsp[STS_SIZE], /* pressure field file name */
# if CMPRSSBL != 0
      dnsdsp[STS_SIZE], /* density field file name */
# endif
      veldsp[STS_SIZE], /* velocity log file name */
      dnslog[STS_SIZE]; /* density log file name */

   long
      lwfcl, /* lower fluid cell label */
      upfcl, /* upper fluic cell label */
      fldnn, /* fluid dynamic iterations counter */
      tmpfp, /* temperature field position indicator [ returned by ftell(*)] */
      tmpef, /* temperature field eof position indicator [ dito ] */
      prsfp, /* pressure field position indicator [ dito ] */
      prsef, /* pressure field eof position indicator [ dito ] */
# if CMPRSSBL != 0
      dnsfp, /* density field position indicator [ dito ] */
      dnsef, /* density field eof position indicator [ dito ] */
# endif
      velfp, /* fluid velocity field position indicator [ dito ] */
      velef; /* fluid velocity field eof position indicator [ dito ] */

   double
      fldstart; /* start fluid dynamic operations with that delay [s] */

   CLUSTER
     *prs,
# if CMPRSSBL != 0
     *dns,
     *flw[THREE],
# endif
     *vel[THREE];
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
};
/************************** end of file solvtp.h ******************************/
