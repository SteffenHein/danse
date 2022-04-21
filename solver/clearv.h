/* [ file: clearv.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C subroutine clearv(*)                                                *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   DSC system parameter reset function                                        *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* the following macro should be defined in "../CONFIG.H"              */
# ifndef DSC_ADJCLL      /* assign neighbouring cell index top.mn[][k] thus:  */
   # define DSC_ADJCLL 0 /* 0: to faces [ k is a face index; 0 <= k < 6 ]     */
# endif                  /* 1: to ports [ k is a port index; 0 <= k < 12 ]    */
/*============================================================================*/
short
clearv( long nlm, long slm, long gel, long ges, long gml, long gms, long xel,
        long xhl, long bnp, long bnn, long vep, long ven, long vhp, long vhn )
{ 
/* allusions: */
/*
   extern struct topology top;
   extern struct tlmsmx smx;
   extern struct boundary bnd;
   extern struct excitation exc;
   extern struct evaluation val;

# if DSC_DOMAIN == 0
   extern DSC_FIELDS fld[];
   extern DSC_DEFLCT dfl;

# elif DSC_DOMAIN == 1
   extern DSC_FIELDS fld;
   extern DSC_GGEFLD ggf;
   extern DSC_DEFLCT dfl;

# else
   extern DSC_FIELDS fld[];
# endif

# if DSC_HCRMDE != 0
   extern DSC_HCRRTS hcr[];

# if DSC_HCGAUGE == 2
   extern DSC_GGECRR ggc[];
# endif
# endif
*/
/* declarations: */

   static struct solverstat
     *state = &solver;
     
   static long 
      hh = null, 
      ii = null;

   static short 
      jj = null,
      kk = null,
      pp = null;

   static signed char
      cc = null,
      prd = null;
/*----------------------------------------------------------------------------*/
   ( state->dmn ) = 't'; /* default: time domain computation */

   ( state->swing ) = ZERO; /* [ time dependent ] EMfield excitation factor */
   ( state->ttime ) = ZERO; /* actual time [ counter ] */
   ( state->dt ) = ZERO; /* micro time step */
   ( state->fr ) = ZERO; /* frequency */
   ( state->ph ) = ZERO; /* phase shift */

   ( state->nn ) = null;
/*............................................................................*/
# if DSC_HCRMDE != 0 /* heat propagation process parameters */
   ( state->hctme ) = ZERO; /* actual time [ counter ] */
   ( state->hcdt ) = ZERO;  /* micro time step */
   ( state->hcswg ) = ZERO; /* [ time dependent ] thermal excitation factor */
   ( state->hsk ) = ZERO; /* heat sources [ sum over skin effect lossy */ 
                          /* boundaries ] */

   ( state->hcrnn ) = -ONE; /* reset thermal interations counter */
   ( state->hcrstart ) = ZERO;  /* start delay for thermal computations */

   ( state->sws ) = null; /* 'switched sources' operation mark */
   ( state->hclbl ) = null;
/*............................................................................*/
# if DSC_FLDMDE != 0
   ( state->fldnn ) = -ONE; /* reset fluid interations counter */
   ( state->fldstart ) = ZERO;  /* start delay for fluid dynamic computations */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   if ( nlm != null )
   {
      if ( NODES < nlm )
         nlm = NODES;

      if ( NRSMX < slm )
         slm = NRSMX;

      prd = 2*( null < bnd.p );

      hh = ONE; 
      while( hh <= nlm )
      {
/*............................................................................*/
         pp = null;
         do
         {
            top.m[hh][pp] = null;
            top.f[hh][pp] = null;
            top.p[hh][pp] = null;
         } while (( ++pp ) < FACES );
         do
         {
/*............................................................................*/
# if DSC_ADJCLL == 1
            top.m[hh][pp] = null;
# endif
/*............................................................................*/
            top.p[hh][pp] = null;
         } while (( ++pp ) < PORTS );
/*............................................................................*/
         cc = null; do
         {
            ii = hh + cc*nlm;
/*............................................................................*/
# if DSC_DOMAIN == 1

            pp = null; do
            {
               fld.r[ii][pp] = ZERO;
               fld.i[ii][pp] = ZERO;
            } while (( ++pp ) < PORTS );

            kk = null; do
            {
               ggf.r[ii][kk] = ZERO;
               ggf.i[ii][kk] = ZERO;
            } while (( ++kk ) < DIMNS );
/*............................................................................*/
# elif DSC_DOMAIN == 2 
            jj = null; do
            {
               pp = null; do
               {
                  fld[jj].r[ii][pp] = ZERO;
                  fld[jj].i[ii][pp] = ZERO;
               } while (( ++pp ) < PORTS );
            } while (( ++jj ) < TWO );
/*............................................................................*/
# elif DSC_DOMAIN == 3 
            jj = null; do
            {
               pp = null; do
               {
                  fld[jj].r[ii][pp] = ZERO;
                  fld[jj].i[ii][pp] = ZERO;
               } while (( ++pp ) < PORTS );
            } while (( ++jj ) < SIX );
/*............................................................................*/
# else /* default: DSC_DOMAIN == 0 */
            jj = null; do
            {
               pp = null; do
               {
                  fld[jj].r[ii][pp] = ZERO;
                  fld[jj].i[ii][pp] = ZERO;
               } while (( ++pp ) < PORTS );
            } while (( ++jj ) < TWO );
# endif /* DSC_DOMAIN == ... */
/*............................................................................*/
         } while (( ++cc ) < prd );
         hh++ ;
      }; /* next hh */
/*............................................................................*/
# if DSC_HCRMDE != 0

      hh = null; 
      while( hh <= HCNDS )
      {
/*............................................................................*/
# if DSC_FLDMDE != 0

         hcs.cnn[hh] = null; /* hcs.cnn[hh]: fluid cnn. comp., cell hh [>0] */
                             /* hcd.cnn[null]: counter of fluid cnn.compts. */
         pp = null; do
         {
            hcs.fctype[hh][pp] = null; /* fluid boundary type identifier */
         } while (( ++pp ) < FACES );
# endif
/*............................................................................*/
         cc = null; do
         {
            hcr[cc].tn[hh] = ZERO;
            hcr[cc].qn[hh] = ZERO;
/*............................................................................*/
# if DSC_FLDMDE != 0
            hcr[cc].pn[hh] = ZERO;
/*............................................................................*/
# if CMPRSSBL != 0
            hcr[cc].rn[hh] = ZERO; /* nodal fluid density */
# endif
/*............................................................................*/
            jj = null; do
            {
               hcr[cc].un[hh][jj] = ZERO;
            } while (( ++jj ) < DIMNS );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
            pp = null; do
            {
               hcr[cc].tf[hh][pp] = ZERO;
               hcr[cc].tt[hh][pp] = ZERO;
               hcr[cc].ic[hh][pp] = ZERO;
/*............................................................................*/
# if DSC_FLDMDE != 0
                                              /* fluid cell face */
               hcr[cc].pf[hh][pp] = ZERO;
               hcr[cc].pg[hh][pp] = ZERO; /* <pressure gradiend | face vector>*/
               hcr[cc].pt[hh][pp] = ZERO; /* returned correction parameters */
/*............................................................................*/
# if CMPRSSBL != 0
               hcr[cc].rf[hh][pp] = ZERO; /* face density */
               hcr[cc].rg[hh][pp] = ZERO; /* <density gradient | face vector> */
               hcr[cc].rt[hh][pp] = ZERO; /* returned correction parameters */
# endif
/*............................................................................*/
               jj = null; do
               {
                  hcr[cc].uf[hh][pp][jj] = ZERO;
                  hcr[cc].ut[hh][pp][jj] = ZERO;
                  hcr[cc].iu[hh][pp][jj] = ZERO;
               } while (( ++jj ) < DIMNS );
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
            } while (( ++pp ) < FACES );
/*............................................................................*/
# if DSC_HCGAUGE == 2
            kk = null; do
            {
               ggc[cc].gc[hh][kk] = ZERO;
            } while (( ++kk ) < FACES );
# endif
/*............................................................................*/
         } while (( ++cc ) < DSC_HCRMDE );
         hh++ ;
      }; /* next hh, while ( hh <= HCNDS ) */
# endif /* DSC_HCRMDE == 0 */
/*............................................................................*/
      top.n = null;

      hh = ONE; 
      while( hh <= slm )
      {
         jj = null; do
         {
            kk = null; do
            {
               smx.se[hh][jj][kk] = ZERO;
               smx.sh[hh][jj][kk] = ZERO;
               kk++ ;
            } while ( kk < SENTR );
            jj++ ;
         } while ( jj < SENTR );
         hh++ ;
      }; /* next hh */

      smx.hh[null] = null;
      smx.shr = ZERO;
      smx.shi = ZERO;
   };

   if ( gml != null )
   {
      if ( GMNDS < gml )
         gml = GMNDS;

     if ( GMSMX < gms ) 
         gms = GMSMX;

      prd = 2*( null < bnd.p ); 
         
      hh = null; 
      while( hh < gml )
      { 
         cc = null;
         while( cc < prd )
         {
            ii = hh + cc*gml;
/*............................................................................*/
# if DSC_DOMAIN == 0
            jj = null; do
            {
               dfl.dm[ii][jj] = ZERO;
               dfl.gm[ii][jj] = ZERO;
               jj++;
            } while (( ++jj ) < THREE );
# elif DSC_DOMAIN == 1
            jj = null; do
            {
               dfl.dm[ii][jj] = ZERO;
               dfl.gm[ii][jj] = ZERO;
            } while (( ++jj ) < THREE );
# endif
/*............................................................................*/
            cc++ ;
         };
         hh++ ;
      };

      hh = null;
      while( hh < gms )
      {
         jj = null; do
         {
            kk = null; do
            {
               smx.dm[hh][jj][kk] = ZERO;
               smx.sgm[hh][jj][kk] = ZERO;
               smx.eum[hh][jj][kk] = ZERO;
            } while (( ++kk ) < THREE );
         } while (( ++jj ) < THREE );
         hh++ ;
      };
/*............................................................................*/
# if DSC_DOMAIN == 0
      dfl.mmx = ZERO;
# elif DSC_DOMAIN == 1
      dfl.mmx = ZERO;
# endif
/*............................................................................*/
      smx.gm   = null;
   };

   if ( gel != null )
   {
      if ( GENDS < gel )
         gel = GENDS;

      if ( GESMX < ges )
         ges = GESMX;

      prd = 2*( null < bnd.p );

      hh = null; 
      while( hh < gel )
      {
         cc = null;
         while( cc < prd )
         {
            ii = hh + cc*gel;
/*............................................................................*/
# if DSC_DOMAIN == 0
            jj = null; do
            {
               dfl.de[ii][jj] = ZERO;
               dfl.ge[ii][jj] = ZERO;
            } while (( ++jj ) < THREE );
# elif DSC_DOMAIN == 1
            jj = null; do
            {
               dfl.de[ii][jj] = ZERO;
               dfl.ge[ii][jj] = ZERO;
            } while (( ++jj ) < THREE );
# endif
/*............................................................................*/
            cc++ ;
         }; 
         hh++ ;
      };

      hh = null;
      while( hh < ges )
      {
         jj = null; do
         {
            kk = null; do
            {
               smx.de[hh][jj][kk] = ZERO;
               smx.sge[hh][jj][kk] = ZERO;
               smx.eue[hh][jj][kk] = ZERO;
            } while (( ++kk ) < THREE );
         } while (( ++jj ) < THREE );
         hh++ ;
      };
/*............................................................................*/
# if DSC_DOMAIN == 0
      dfl.emx = ZERO;
# elif DSC_DOMAIN == 1
      dfl.emx = ZERO;
# endif
/*............................................................................*/
      smx.ge  = null;
   };

   exc.lbl = -ONE;

   if ( xel != null )
   {
      if ( EXCEP < xel )
         xel = EXCEP;

      pp = null;
      while( pp < xel )
      {
         exc.me[pp] = null;
         exc.pe[pp] = null;
         exc.er[pp] = ZERO;
/*............................................................................*/
# if DSC_DOMAIN != 1
         exc.ei[pp] = ZERO;
# endif
/*............................................................................*/
         pp++ ;
      };
      exc.ne = null;
   };

   if ( xhl != null )
   {
      if ( EXCHP < xhl )
         xhl = EXCHP;

      pp = null;
      while( pp < xhl )
      {
         exc.mh[pp] = null;
         exc.ph[pp] = null;
         exc.hr[pp] = ZERO;
         exc.hi[pp] = ZERO;
         pp++ ;
      };
      exc.nh = null; 
   };

   if (( exc.ne == null)&&( exc.nh == null ))
   {
      if ( EXCFR < exc.nn )
         exc.nn = EXCFR;

      jj = null;
      while( jj < exc.nn )
      {
         exc.fr[jj] = ZERO;
         jj++ ;
      };
      exc.nn = null;
      exc.rt = ZERO;
      exc.dt = ZERO;
      exc.ht = ZERO;
      exc.mx = ZERO;
   };

   if ( bnp != null )
   {
      if ( BNDPR < bnp )
         bnp = BNDPR;

      pp = null; 
      while( pp < bnp )
      {
         bnd.cb[pp]   = null;
         bnd.cp[pp]   = null;
         bnd.p0[pp]   = null;
         bnd.p1[pp]   = null;
         bnd.pp0[pp]  = null;
         bnd.pp1[pp]  = null;
         pp++ ;
      };
      bnd.r = ZERO;
      bnd.i = ZERO;
      bnd.p = null;
   };
   
   if ( bnn != null )
   {
      if (  BNDAP < bnn )
         bnn = BNDAP;

      pp = null;
      while( pp < bnn )
      {
         bnd.m[pp]    = null;
         bnd.f[pp]    = null;
         bnd.r00[pp]  = ZERO;
         bnd.r01[pp]  = ZERO;
         bnd.r10[pp]  = ZERO;
         bnd.r11[pp]  = ZERO;

         bnd.i00[pp]  = ZERO;
         bnd.i01[pp]  = ZERO;
         bnd.i10[pp]  = ZERO;
         bnd.i11[pp]  = ZERO;
         pp++ ;
      };
      bnd.n = null;
   };

   if ( vep != null )
   { 
      if ( EVLEP < vep )
         vep = EVLEP;

      pp = null; 
      while( pp < vep )
      {
         val.mep[pp] = null;
         val.pep[pp] = null;
         val.epr[pp] = ZERO;
         val.epi[pp] = ZERO;
         pp++ ;
      };
      val.nep = null;
   };

   if ( ven != null )
   {
      if ( EVLEN < ven )
         ven = EVLEN;

      pp = null; 
      while( pp < ven )
      {
         val.men[pp] = null;
         val.cen[pp] = null;
         val.enr[pp] = ZERO;
         val.eni[pp] = ZERO;
         pp++ ;
      };
      val.nen = null;
   };

   if ( vhp != null )
   {
      if ( EVLHP < vhp )
         vhp = EVLHP;

      pp = null; 
      while( pp < vhp )
      {
         val.mhp[pp] = null;
         val.php[pp] = null;
         val.hpr[pp] = ZERO;
         val.hpi[pp] = ZERO;
         pp++ ;
      };
      val.nhp = null; 
   };

   if ( vhn != null )
   {
      if ( EVLHN < vhn )
         vhn = EVLHN;

      pp = null;
      while( pp < vhn )
      {
         val.mhn[pp] = null;
         val.chn[pp] = null;
         val.hnr[pp] = ZERO;
         val.hni[pp] = ZERO;
         pp++ ;
      };
      val.nhn = null;
   };

   if (( val.nep == null)&&( val.nen == null )&&
       ( val.nhp == null)&&( val.nhn == null))
   {
      val.read = null;
      val.n  = null;
      val.ni = null;
      val.r  = null;
      val.dt = ZERO; 
   };
/*............................................................................*/
# if DSC_HCRMDE != 0
   ii = null;
   while( ii < bnd.nsk )
   {
      bnd.msk[ii] = null;
      bnd.fsk[ii] = null;
      bnd.hs[ii] = ZERO;
      bnd.rs[ii] = ZERO;
      bnd.sk[ii][null][null] = ZERO;
      bnd.sk[ii][null][ONE] = ZERO;
      bnd.sk[ii][ONE][null] = ZERO;
      bnd.sk[ii][ONE][ONE] = ZERO;
      ii++;
   };
   bnd.nsk = null;

   ii = null;
   while( ii < bnd.nhc )
   {
      bnd.mhc[ii] = null;
      bnd.fhc[ii] = null;
      bnd.hc[ii] = ZERO;
      ii++;
   };
   bnd.nhc = null;

   ii = null;
   while( ii < bnd.ntf )
   {
      bnd.mtf[ii] = null;
      bnd.ftf[ii] = null;
      bnd.tf[ii] = ZERO;
      ii++;
   };
   bnd.ntf = null;

   ii = null;
   while( ii < bnd.ntn )
   {
      bnd.mtn[ii] = null;
      bnd.tn[ii] = ZERO;
      ii++;
   };
   bnd.ntn = null;

   ii = null;
   while( ii < bnd.nsc )
   {
      bnd.msc[ii] = null;
      bnd.fsc[ii] = null;
      bnd.sc[ii] = ZERO;
      bnd.tr[ii] = ZERO;
      ii++;
   };
   bnd.nsc = null;

   ii = null;
   while( ii < exc.nhc )
   {
      exc.mhc[ii] = null;
      exc.fhc[ii] = null;
      exc.hc[ii] = ZERO;
      ii++;
   };
   exc.nhc = null;

   ii = null;
   while( ii < exc.ntf )
   {
      exc.mtf[ii] = null;
      exc.ftf[ii] = null;
      exc.tf[ii] = ZERO;
      ii++;
   };
   exc.ntf = null;

   ii = null;
   while( ii < exc.ntn )
   {
      exc.mtn[ii] = null;
      exc.tn[ii] = ZERO;
      ii++;
   };
   exc.ntn = null;

   kk = null;
   while ( kk < DSC_HCRMDE )
   {  
      ii = null;
      while( ii < val.nhc[kk] )
      {
         val.mhc[kk][ii] = null;
         val.fhc[kk][ii] = null;
         val.hc[kk][ii] = ZERO;
         ii++;
      };
      val.nhc[kk] = null;

      ii = null;
      while( ii < val.ntf[kk] )
      {
         val.mtf[kk][ii] = null;
         val.ftf[kk][ii] = null;
         val.tf[kk][ii] = ZERO;
         ii++;
      };
      val.ntf[kk] = null;

      ii = null;
      while( ii < val.ntn[kk] )
      {
         val.mtn[kk][ii] = null;
         val.tn[kk][ii] = ZERO;
         ii++;
      };
      val.ntn[kk] = null;
      kk++ ;
   }  while ( kk < DSC_HCRMDE );
# endif  
/*............................................................................*/
   return ONE;
}
/*============================================================================*/
/********** end of DSC system parameter clearing function clearv(*) ***********/
