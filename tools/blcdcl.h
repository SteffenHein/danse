/* [ file: blcdcl.h ] */
/* Update: April 09, 2022 */
/*----------------------------------------------------------------------------*/
/* variable declaration, macro and auxiliary function definition header of */
/* function blocks(*) */

# define BLC_DECLARE( ) \
 \
/* allusions: */ \
/* \
   extern FORMSTATE *spt; \
   extern BLOCSTR blc; \
   extern struct transfer trf; \
   extern struct labels lbl; \
   extern struct quadrl qdl; \
   extern struct triangle tgl; \
*/ \
 \
/* declarations: */ \
 \
   time_t \
      nseconds = null, \
     *timer = null; \
 \
   static short \
      cntii = null; \
\
   static char \
     *tmeptr, \
     *blclog = "blc.log"; \
 \
   static char \
      ptr[STS_SIZE] = {null}; \
 \
   static const signed char \
      a = null, \
      b = ONE, \
      c = TWO, \
      d = THREE, \
      ab = null, \
      bc = ONE, \
      ad = TWO, \
      dc = THREE, \
      ac = TWO, \
      forward = ONE, \
      backward = -ONE, \
      ewall = -ONE, \
      mwall = -TWO; \
 \
/* prototyping: */ \
 \
   short qudrl( short, short, char *option ); \
   short trngl( short, char *option ); \
 \
   void quadrangle ( short Block, short Layer, char *option ); \
   void triangle( short Block, short Layer, char *option ); \
 \
   void interface ( short ActualBlock, signed char ActualSide, \
                               short StartActual, short StopActual, \
                    short FormerBlock, signed char FormerSide, \
                               short StartFormer, signed char Sense ); \
   int leave( short block, \
           short layer1, short layer2, short layer ); \
 \
   void readout( short, char ); \
   void vertex( signed char, short ); \
 \
   void smedium( short medium_idx, long init_cell, long final_cell, \
                 double eps, double myr, double ke, double km, \
                 double kh, double cv, char *option ); 
/*----------------------------------------------------------------------------*/
/* configure Maxwell field boundaries: */
/*============================================================================*/

void blc_boundary( char tp, long ii, short nn, short pp, char ff, double dd )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      jj = null;

/* prototypes: */

# if DSC_HCRMDE != 0
   void blc_hcrbound( char tp, long ii, short nn, short pp,
      char ff, double hh );
# endif
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.bnd;

      if( MOD_BOUNDRS <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many boundary blocks defined !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_BOUNDRS.",
            ( long ) jj, ( long ) MOD_BOUNDRS );
         printf( "\n - Change macro MOD_BOUNDRS in model.c"
            "\n   in compliance with memory ressources.]\n" );
         exit( EXIT_FAILURE );
      };

      tp = tolower(tp);

      switch(tp) /* check for legal boundary types: */
      {
        case 'c': /* coaxial [TEM mode] */
        case 'r': /* rectangular waveguide [TEmn modes] */
        case 's': /* circular waveguide [TEmn modes] */
        case 'e': /* elliptic waveguide [c_TEmn modes] */
        case 'o': /* resistive surface losses [skin effect, e.g.] */
         break;

        case 'f': /* fixed face temperatures */
        case 'n': /* fixed node temperatures */
        case 'h': /* heat current density imposed on boundary face */
        case 'b': /* Stefan-Boltzmann [heat radiative] environmental exchange */
        case 'x': /* dito */
        case 'k': /* fixed surface heat conductance [ apprx. conv., e.g.] */
        case 'l': /* heat sources, on [or adjacent to] skin effect lossy bnd */
        case 'u': /* fluid slip boundary condition */
/*............................................................................*/
# if DSC_HCRMDE == 0

         printf( "\n\n Message from function blc_boundary(*):" );
         printf( "\n Thermal boundary type '%c' remains inoperative", (tp));
         printf( "\n due to uncompiled heat propagation mode !!!" );
         printf( "\n [ To implement heat propagation " );
         printf( "\n   set macro DSC_HCRMDE in file CONFIG.H to 1" );
         printf( "\n   and recompile program package." );
         printf( "\n Program continues [yet ignoring heat "
            "propagation modes]." );
         return;
# else /* DSC_HCRMDE != 0 */
/*............................................................................*/
# if DSC_FLDMDE == 0 

         if ( tp == 'u' )
         {
            printf( "\n\n Message from function blc_boundary(*):" );
            printf( "\n Fluid flow boundary type '%c' remains "
              "inoperative", (tp));
            printf( "\n due to uncompiled fluid flow computation mode !!!" );
            printf( "\n [ To implement fluid flow computation " );
            printf( "\n   set macro DSC_FLDMDE in file CONFIG.H to 1" );
            printf( "\n   and recompile program package." );
            printf( "\n Program continues [yet ignoring fluid "
               "flow computation modes]." );
            return;
         };
# endif /* DSC_FLDMDE == 0 */
/*............................................................................*/
         blc_hcrbound( tp, ii, nn, pp, ff, dd );
         return;

# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         break;

        default:
         printf( "\n Message from function blc_boundary(*):" );
         printf( "\n Unknown boundary type '%c' !!!", (tp));
         printf( "\n Legal electric boundary types are:" );
         printf( "\n C,c - coaxial [TEM mode]" ); 
         printf( "\n R,r - rectangular waveguide [TEmn modes]" ); 
         printf( "\n S,s - circular waveguide [TEmn modes]" ); 
         printf( "\n E,e - elliptic waveguide [c_TEmn modes]" ); 
         printf( "\n O,e - resitive [e.g. skin effect lossy] boundary" ); 
/*............................................................................*/
# if DSC_HCRMDE != 0
         printf( "\n Legal thermal boundary types are:" );
         printf( "\n F,f - fixed face temperature" ); 
         printf( "\n N,n - fixed node temperature" ); 
         printf( "\n H,h - heat current density imposed on face" ); 
         printf( "\n B,b - Stefan-Boltzmann [radiative] environmental "
            "exchange" ); 
         printf( "\n X,x - [dito - same as preceding]" ); 
         printf( "\n K,k - fixed surface heat resistivity "
            "[apprx. convection]" ); 
         printf( "\n L,l - heat sources on [or adjacent to] skin effect "
            "lossy face" );
/*............................................................................*/
# if DSC_FLDMDE != 0 
         printf( "\n U,u - fluid flow slip boundary face" ); 
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
         break;
      };

      trf.btp[jj] = (tp); /* c/r etc.: coaxial/rectangl.waveguide etc.[type] */
      trf.bii[jj] = (ii); /* initial boundary cell label */
      trf.bnn[jj] = (nn); /* number of periods */
      trf.bpp[jj] = (pp); /* period */
      trf.bff[jj] = (ff); /* cell face index */
      trf.bdd[jj] = (dd); /* any real parameter [ permittivity, e.g.] */
      trf.bnd++ ; /* count this boundary section */
   }; /* end if blc.cov == null */
   return;
}
/*----------------------------------------------------------------------------*/
/* configure Maxwell field excitation: */
/*============================================================================*/

void blc_excite( char tp, long ii, short nn, short pp, char ff, double uu )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short 
      jj = null;

/* prototypes: */

# if DSC_HCRMDE != 0
   void blc_hcrexcit( char tp, long ii, short nn, short pp, 
     char ff, double hh );
# endif
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.exc;

      if( MOD_EXCITES <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many blocks excited !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_EXCITES.",
            ( long ) jj, ( long ) MOD_EXCITES );
         printf( "\n - Change macro MOD_EXCITES in model.c"
            "\n   in compliance with memory resources.]\n" );
         exit( EXIT_FAILURE );
      };

      tp = tolower(tp);

      switch(tp) /* check for legal excitation options: */
      {
        case 'e': /* electric field excitation */
        case 'h': /* magnetic field excitation */
         break;

        case 'c': /* heat current density excited [imposed] on face */
        case 'f': /* face temperatures excitation */
        case 'n': /* node temperatures excitation */

# if DSC_HCRMDE == 0

         printf( "\n\n Message from function blc_excite(*):" );
         printf( "\n Thermal excitation type '%c' remains inoperative", (tp));
         printf( "\n due to uncompiled heat propagation mode !!!" );
         printf( "\n [ To implement heat propagation " );
         printf( "\n   set macro DSC_HCRMDE in file CONFIG.H to 1" );
         printf( "\n   and recompile program package." );
         printf( "\n Program continues [yet ignoring heat "
            "propagation modes]." );
         return;

# else /* DSC_HCRMDE != 0 */

         blc_hcrexcit( tp, ii, nn, pp, ff, uu );
         return;

# endif /* DSC_HCRMDE != 0 */

         break;

        default:
         printf( "\n Message from function blc_excite(*):" );
         printf( "\n Unknown excitation type '%c' !!!", (tp));
         printf( "\n Legal field exciation types are:" );
         printf( "\n E,e - electric field" );
         printf( "\n H,h - magnetic field" );

# if DSC_HCRMDE != 0

         printf( "\n Legal thermal excitation types are:" );
         printf( "\n F,f - face temperature excitation" );
         printf( "\n N,n - node temperature excitation" );

# endif /* DSC_HCRMDE != 0 */

         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
         break;
      }; /* end switch(tp) */

      trf.etp[jj] = (tp); /* type of excitation port */
      trf.eii[jj] = (ii); /* initial excited face label */
      trf.enn[jj] = (nn); /* number of periods */
      trf.epp[jj] = (pp); /* period */
      trf.eff[jj] = (ff); /* cell face index */
      trf.euu[jj] = (uu); /* voltage normalization factor */
      trf.exc++ ; /* count this excitation section */
   }; /* end if blc.cov == null */
   return;
}
/*----------------------------------------------------------------------------*/
/* configure Maxwell field evaluation: */
/*============================================================================*/

void blc_evaluate( char tp, long ii, short nn, short pp, char pt, char *txt )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short 
      jj = null;

/* prototypes: */

# if DSC_HCRMDE != 0 /* heat current modes */
   void blc_hcrevlte( char tp, long ii, short nn, short pp, 
      char ff, char *txt );
# endif
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.val;

      if( MOD_EVLUATE <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many blocks to be evaluated !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_EVLUATE.",
            ( long ) jj, ( long ) MOD_EVLUATE );
         printf( "\n - Change macro MOD_EVLUATE in model.c"
            "\n   in compliance with memory resources.]\n" );
         exit( EXIT_FAILURE );
      };

      tp = tolower(tp);

      switch(tp) /* check for legal evaluation options: */
      {
        case 'e': /* electric field evaluation */
        case 'h': /* magnetic field evaluation */
         break;

        case 'f': /* face temperatures evaluation */
        case 'n': /* node temperatures evaluation */
        case 'u': /* nodal velocity evaluation */
/*............................................................................*/
# if DSC_HCRMDE == 0
         printf( "\n\n Message from function blc_evaluate(*):" );
         printf( "\n Thermal evaluation type '%c' remains inoperative", (tp));
         printf( "\n due to uncompiled heat propagation mode !!!" );
         printf( "\n [ To implement heat propagation " );
         printf( "\n   set macro DSC_HCRMDE in file CONFIG.H to 1" );
         printf( "\n   and recompile program package." );
         printf( "\n Program continues [yet ignoring heat "
            "propagation modes]." );
         return;
# else /* DSC_HCRMDE != 0 */
/*............................................................................*/
# if DSC_FLDMDE == 0
         if ( jj == 'u' )
         {
            printf( "\n\n Message from function blc_evaluate(*):" );
            printf( "\n Fluid velocity evaluation, type '%c', remains "
                    "inoperative", (tp));
            printf( "\n due to uncompiled fluid flow computation mode !!!" );
            printf( "\n [ To implement fluid flow computation " );
            printf( "\n   set macro DSC_FLDMDE in file CONFIG.H to 1" );
            printf( "\n   and recompile program package." );
            printf( "\n Program continues [yet ignoring fluid "
               "flow computation modes]." );
            return;
         };
# endif /* DSC_FLDMDE == 0 */
/*............................................................................*/
         blc_hcrevlte( tp, ii, nn, pp, pt, txt );
         return;

# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         break;

        default:
         printf( "\n Message from function blc_evaluate(*):" );
         printf( "\n Unknown evaluation type '%c' !!!", (tp));
         printf( "\n Legal field evaluation modes are:" );
         printf( "\n E,e - electric field" );
         printf( "\n H,h - magnetic field" );

# if DSC_HCRMDE != 0
         printf( "\n Legal thermal evaluation modes are:" );
         printf( "\n F,f - face temperature evaluation" );
         printf( "\n N,n - node temperature evaluation" );
# endif /* DSC_HCRMDE != 0 */

         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
         break;
      }; /* end switch(tp) */

      if(( null < jj )&&( null == strncmp( "reference", txt, THREE )))
      { /* set reference port at first place: */
         while( null < jj )
         {
            jj-- ;
            trf.vtp[jj+ONE] = trf.vtp[jj];
            trf.vii[jj+ONE] = trf.vii[jj];
            trf.vnn[jj+ONE] = trf.vnn[jj];
            trf.vpp[jj+ONE] = trf.vpp[jj];
            trf.vpt[jj+ONE] = trf.vpt[jj];
            strncpy( trf.vtx[jj+ONE], trf.vtx[jj], STS_SIZE );
         };
      };
         
      trf.vtp[jj] = (tp); /* type of evaluation port */
      trf.vii[jj] = (ii); /* initial evaluated face label */
      trf.vnn[jj] = (nn); /* number of periods */
      trf.vpp[jj] = (pp); /* period */
      trf.vpt[jj] = (pt); /* cell port index */
      strncpy( trf.vtx[jj], txt, STS_SIZE );
      trf.val++ ; /* count this evaluation section */
   }; /* end if blc.cov == null */
   return;
}
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0 /* heat current modes */
/*----------------------------------------------------------------------------*/
/* configure heat current and temperature boundaries: */
/*============================================================================*/

void blc_hcrbound( char tp, long ii, short nn, short pp, char ff, double hh )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      jj = null;
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.bhc; /* counter */

      if( MOD_HCBNDRS <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many heat current boundary blocks defined !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_HCBNDRS.",
            ( long ) jj, ( long ) MOD_HCBNDRS );
         printf( "\n - Change macro MOD_HCBNDRS in model.c"
            "\n   in compliance with memory ressources.]\n" );
         exit( EXIT_FAILURE );
      };

      tp = tolower(tp);

      switch(tp) /* check for legal boundary types: */
      {
        case 'f': /* fixed face temperatures */
        case 'n': /* fixed node temperatures */
        case 'h': /* heat current density imposed on face */
        case 'b': /* Stefan-Boltzmann [heat radiative] environmental exchange */
        case 'x': /* dito */
        case 'k': /* fixed surface heat resistivity [apprx. convection] */
        case 'l': /* heat sources, on [or adjacent to] skin effect lossy bnd */
        case 'u': /* fluid slip boundary conditions */
         break;

        default:
         printf( "\n Message from function blc_hcrbound(*):" );
         printf( "\n Unknown boundary type '%c' !!!", (tp));
         printf( "\n Legal thermal boundary types are:" );
         printf( "\n F,f - fixed face temperature" ); 
         printf( "\n N,n - fixed node temperature" ); 
         printf( "\n H,h - heat current density imposed on face" ); 
         printf( "\n B,b - Stefan-Boltzmann heat radiative face" );
         printf( "\n X,x - same as preceding" );
         printf( "\n L,l - heat sources on [or adjacent to] skin effect "
            "lossy face" );
# if DSC_FLDMDE != 0
         printf( "\n U,u - fluid slip boundary conditions on face" ); 
# endif
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
         break;
      }; /* end switch(tp) */

      trf.bct[jj] = (tp); /* heat current boundary type */
      trf.bci[jj] = (ii); /* initial boundary cell label */
      trf.bcn[jj] = (nn); /* number of periods */
      trf.bcp[jj] = (pp); /* period */
      trf.bcf[jj] = (ff); /* cell face index */
      trf.bch[jj] = (hh); /* heat current density, temperature, or suchlike */
      trf.bhc++ ;
   }; /* end if blc.cov == null */
   return;
}
/*----------------------------------------------------------------------------*/
/* configure heat current and temperature excitation: */
/*============================================================================*/

void blc_hcrexcit( char tp, long ii, short nn, short pp, char ff, double hh )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short 
      jj = null;
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.ehc;

      if( MOD_HCEXCTS <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many current blocks excited !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_HCEXCTS.",
            ( long ) jj, ( long ) MOD_HCEXCTS );
         printf( "\n - Change macro MOD_HCEXCTS in model.c"
            "\n   in compliance with memory resources.]\n" );
         exit( EXIT_FAILURE );
      };

      tp = tolower(tp);

      switch(tp) /* check for legal thermal excitation types */
      {
        case 'c': /* heat current density excited [imposed] on face */
        case 'f': /* face temperature excitation */
        case 'n': /* node temperature excitation */
         break;

        default:
         printf( "\n Message from function blc_hcrexcit(*):" );
         printf( "\n Unknown excitation type '%c' !!!", (tp));
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
         break;
      }; /* end switch(tp) */

      trf.ect[jj] = (tp); /* type of excitation port */
      trf.eci[jj] = (ii); /* initial excited cell label */
      trf.ecn[jj] = (nn); /* number of periods */
      trf.ecp[jj] = (pp); /* period */
      trf.ecf[jj] = (ff); /* cell face index */
      trf.ech[jj] = (hh); /* heat current density, temperature, or suchlike */
      trf.ehc++ ;
   }; /* end if blc.cov == null */
   return;
}
/*----------------------------------------------------------------------------*/
/* configure heat current and temperature evaluation: */
/*============================================================================*/

void blc_hcrevlte( char tp, long ii, short nn, short pp, char ff, char *txt )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short 
      jj = null;
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.vhc;

      if( MOD_HCEVALS <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many heat current blocks to be evaluated !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_HCEVALS.",
            ( long ) jj, ( long ) MOD_HCEVALS );
         printf( "\n - Change macro MOD_HCEXCTS in model.c"
            "\n   in compliance with memory resources.]\n" );
         exit( EXIT_FAILURE );
      };

      tp = tolower(tp);

      switch(tp) /* check for legal evaluation options: */
      {
        case 'f': /* face temperatures evaluation */
        case 'n': /* node temperatures evaluation */
        case 'u': /* nodal velocities evaluation */
         break;

        default:
         printf( "\n Message from function blc_evaluate(*):" );
         printf( "\n Unknown evaluation type '%c' !!!", (tp));
         printf( "\n Legal thermal evaluation modes are:" );
         printf( "\n F,f - face temperature evaluation" );
         printf( "\n N,n - node temperature evaluation" );
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
         break;
      }; /* end switch(tp) */

      trf.vct[jj] = (tp); /* type of evaluation port */
      trf.vci[jj] = (ii); /* initial evaluated cell label */
      trf.vcn[jj] = (nn); /* number of periods */
      trf.vcp[jj] = (pp); /* period */
      trf.vcf[jj] = (ff); /* cell face index, or component */
      strncpy( trf.vcx[jj], txt, STS_SIZE );
      trf.vhc++ ;
   }; /* end if blc.cov == null */
   return;
}
/*============================================================================*/
# endif /* DSC_HCRMDE != 0 */
/*----------------------------------------------------------------------------*/
/* macros: */
/*----------------------------------------------------------------------------*/
# define BLC_TRIVIAL(AA,BB,TP) \
{ \
   if (( blc.cov == null ) \
     &&((AA) <= blc.dmn ) \
     &&( blc.dmn <= (BB))) \
   { \
      if(((TP) == 'T' )||((TP) == 't' )) \
         tgl.trv = e_wall; \
      else \
         qdl.trv = e_wall; \
   }; \
}
/*----------------------------------------------------------------------------*/
# define BLC_BOUNDARY(TP,II,NN,PP,FF,DD) \
{ \
   if ( *option == 'c' ) \
      blc_boundary((TP),(II),(NN),(PP),(FF),(DD)); \
}
/*............................................................................*/
# define BOUNDARY(TP,II,NN,PP,FF,DD) \
   BLC_BOUNDARY((TP),(II),(NN),(PP),(FF),(DD))
/*----------------------------------------------------------------------------*/
# define BLC_EXCITE(TP,II,NN,PP,FF,UU) \
{ \
   if ( *option == 'c' ) \
      blc_excite((TP),(II),(NN),(PP),(FF),(UU)); \
}
/*............................................................................*/
# define EXCITE(TP,II,NN,PP,FF,UU) \
   BLC_EXCITE((TP),(II),(NN),(PP),(FF),(UU)) \
/*----------------------------------------------------------------------------*/
# define BLC_EVALUATE(TP,II,NN,PP,PT,TXT) \
{ \
   if ( *option == 'c' ) \
      blc_evaluate((TP),(II),(NN),(PP),(PT),(TXT)); \
}
/*............................................................................*/
# define EVALUATE(TP,II,NN,PP,PT,TXT) \
   BLC_EVALUATE((TP),(II),(NN),(PP),(PT),(TXT))
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
/*............................................................................*/
# define BLC_HCRBOUND(TP,II,NN,PP,FF,HH) \
{ \
   if ( *option == 'c' ) \
      blc_hcrbound((TP),(II),(NN),(PP),(FF),(HH)); \
}
/*............................................................................*/
# define HCRBOUND(TP,II,NN,PP,FF,HH) \
{ \
   if ( *option == 'c' ) \
      blc_hcrbound((TP),(II),(NN),(PP),(FF),(HH)); \
}
/*----------------------------------------------------------------------------*/
# define BLC_HCREXCIT(TP,II,NN,PP,FF,HH) \
{ \
   if ( *option == 'c' ) \
      blc_hcrexcit((TP),(II),(NN),(PP),(FF),(HH)); \
}
/*............................................................................*/
# define HCREXCIT(TP,II,NN,PP,FF,HH) \
{ \
   if ( *option == 'c' ) \
      blc_hcrexcit((TP),(II),(NN),(PP),(FF),(HH)); \
}
/*----------------------------------------------------------------------------*/
# define BLC_HCREVLTE(TP,II,NN,PP,FF,TXT) \
{ \
   if ( *option == 'c' ) \
      blc_hcrevlte((TP),(II),(NN),(PP),(FF),(TXT)); \
}
# define HCREVLTE(TP,II,NN,PP,FF,TXT) \
{ \
   if ( *option == 'c' ) \
      blc_hcrevlte((TP),(II),(NN),(PP),(FF),(TXT)); \
}
# endif
/*************************** end of file bldclr.h *****************************/
