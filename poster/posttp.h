/* [ file: evaltp.h ] */
/* Update: February 22, 2005 */
/*----------------------------------------------------------------------------*/
# define TYPE_FFT 1
typedef struct
{                     
   signed char
      rtn;

   char
      opt[SHS_SIZE];

   short 
      p, q,
      mult[FTR_NMBR+ONE];

   long
      ttlg[FTR_NMBR+ONE], stlg[FTR_NMBR+ONE];

/* distributions, real and imaginary parts: */
   double
      r[FTR_NMBR+ONE][FTR_SIZE+ONE],
      i[FTR_NMBR+ONE][FTR_SIZE+ONE],

      t[FTR_NMBR+ONE], tt[FTR_NMBR+ONE], dt[FTR_NMBR+ONE],
      s[FTR_NMBR+ONE], ss[FTR_NMBR+ONE], ds[FTR_NMBR+ONE];

/* normalization constant: */
   double
      nor;

} FFT;
/*----------------------------------------------------------------------------*/
# define TYPE_EVALUATE 1
typedef struct
{
   signed char
      rtn;

   char
      name[SHS_SIZE],
      text[STS_SIZE];

   short       
      r, /* repetition rate [ number of averaging iterations per cycle ] */
      nep,
      nhp,
      nen, 
      nhn;

   long        
      n,  /* number of iteration cyles */
      ni, /* first evaluated cycle */
      nf; /* last computed Maxwell field cycle */

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
      mode_ep[SHS_SIZE],
      mode_hp[SHS_SIZE],
      mode_en[SHS_SIZE],
      mode_hn[SHS_SIZE];

   double      
      dt, tc, ti, tf,   /* time_step, cycle -, initial-, final time */
      fr;               /* frequency */

   double
      ep[EVLEP],       hp[EVLHP],
      en[EVLEN],       hn[EVLHN];

   char
      read,
      file[SHS_SIZE],
      domain[SHS_SIZE],
      xunit[SHS_SIZE],
      yunit[SHS_SIZE],
      bndtyp[SHS_SIZE],
      exctyp[SHS_SIZE];
/*...........................................................................*/
# if DSC_HCRMDE != 0 /* thermal evaluation */
   char
      excctp[SHS_SIZE];

   double
      cdt, ctc, ctj, ctf; /* time_step, cylcle_time, initial_time, final_time */

   short  
      rc,               /* [ thermal ] repetition rate */
      ntf[DSC_HCRMDE], /* number of evaluated temperature faces */
      ntn[DSC_HCRMDE], /* number of evaluated temperature nodes */
      nhc[DSC_HCRMDE]; /* number of evaluated heat currrents faces */

   long
      nj, /* initial thermal cycle */
      nt; /* last computed thermal cycle */
   
   long
      mtf[DSC_HCRMDE][EVLTF], /* evaluated face temperature cell indices */
      mtn[DSC_HCRMDE][EVLTN], /* evaluated node temperature cell indices */
      mhc[DSC_HCRMDE][EVLHC]; /* evaluated heat current cell indices */

   signed char
      ftf[DSC_HCRMDE][EVLTF], /* evaluated face temperature face indices */
      fhc[DSC_HCRMDE][EVLHC]; /* evaluated heat current face indices */
/*...........................................................................*/
# if DSC_FLDMDE != 0 /* fluid flows */
   short  
      nun[DSC_HCRMDE];        /* number of evaluated nodal velocities */

   long
      mun[DSC_HCRMDE][EVLUN]; /* mesh cell indices */

   signed char
      cun[DSC_HCRMDE][EVLUN]; /* velocity components */

# endif /* DSC_FLDMDE != 0 */
/*...........................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*...........................................................................*/
   FFT *fpt;

} EVALUATE;
/*----------------------------------------------------------------------------*/
# define TYPE_SPLINES 1
typedef struct
{
   signed char
      rtn;

   long
      mm, /* number of supporting points */
      nn; /* number of interpolated points */

   double
      intgr, fmin, fmax,
      vct[SPL_SPNTS][TWO], /* given points */
      dmn[SPL_INTPL],      /* interpolated points */
      fct[SPL_INTPL],
      drv[SPL_INTPL];

} SPLINES;
/*----------------------------------------------------------------------------*/
# define TYPE_GRAPHICS 1
typedef struct
{
   signed char
      rtn, /* return operation mark: 0: returm with error */
      dsp; /* display operation mark: 1 display some file saving messages */

   char
      name[STS_SIZE],
      text[STS_SIZE];

   char
      file[STS_SIZE],
      curve[GPH_CURVES][STS_SIZE],
      format[SHS_SIZE];

   char
      xunit[SHS_SIZE],
      yunit[SHS_SIZE];

   short
      nc; /* nc = number of graphics [ 'curves' ] */

   long
      nn,
      np[GPH_CURVES]; /* np = number of sample points */

   double
      xmin,
      xmax,
      ymin,
      ymax,
      vct[GPH_POINTS][GPH_CURVES+ONE];

} GRAPHICS;
/*----------------------------------------------------------------------------*/
# define TYPE_POSTSTATE 1
typedef struct
{
   signed char
      rtn;

   char
      name[SHS_SIZE],
      text[STS_SIZE],
      file[SHS_SIZE],
      flbl[FOUR];

   short
      period,
      fldprd,
      hcrprd,
      idx[MAX_PERIOD+ONE];

   long
      fleofs;

   EVALUATE
     *vpt;

} POSTSTATE;
/*************************** end of file posttp.h *****************************/
