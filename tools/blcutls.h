/* [ file: blcutls.h ] */
/* Update: March 19, 2022 */
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

void blc_boundary( char *tp, long ii, short nn, short pp, char ff, double dd )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      jj = null,
      kk = null,
      ll = null;

   static char
      tpp[SHS_SIZE] = {null};

/* prototypes: */

# if DSC_HCRMDE != 0
   void blc_hcrbound( char *tp, long ii, short nn, short pp,
      char ff, double rr );
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
/*............................................................................*/
/* string conversion to lower case: */

      strncpy( tpp, tp, SHS_SIZE );
      ll = strlen( tpp );

      kk = null;
      do
      {
         tpp[kk] = tolower( tpp[kk] );
      } while(( kk++ ) < ll );
/*............................................................................*/
/* electric boundary conditions */
/* [equivalent surface impedances]: */
      
      if (( null == strncmp( tpp, "coaxial", THREE )) /* TEM mode */
        ||( null == strncmp( tpp, "rectangular", THREE )) /* TE10 mode */
        ||( null == strncmp( tpp, "circular", THREE )) /* circ wg, TE11 mode */
        ||( null == strncmp( tpp, "elliptic", THREE )) /* ell wg, TE11 mode */
        ||( null == strncmp( tpp, "surface", THREE )) /* R_square, surf impdce*/
        ||( null == strncmp( tpp, "skin_effect", TWO )) /* dito */
        ||( null == strcmp( tpp, "sr" ))) /* dito */
      {
         goto copy1;
      }                              /* thermal and fluid boundary conditions */
      else if (( null == strncmp( tpp, "heat_current", THREE )) /* fixed heat */
                                     /* Stefan-Boltzmann heat radiation */
             ||( null == strncmp( tpp, "stefan_boltzmann", THREE )) 
             ||( null == strncmp( tpp, "rad", THREE )) /* dito */
             ||( null == strncmp( tpp, "sources", THREE )) /* heat sources */
             ||( null == strncmp( tpp, "no_slip", TWO )) /* no-slip bnd cnds */
             ||( null == strncmp( tpp, "ns", TWO )) /* dito */
             ||( null == strncmp( tpp, "slip", TWO )) /* free slip bdry conds */
             ||( null == strncmp( tpp, "free_slip", TWO )) /* dito */
             ||( null == strncmp( tpp, "inflow", TWO )) /* inflow bndry cds */
             ||( null == strncmp( tpp, "outflow", TWO )) /* outflow bdry cds */
             ||( null == strcmp( tpp, "if" )) /* inflow */
             ||( null == strcmp( tpp, "uf" )) /* dito */
             ||( null == strcmp( tpp, "of" )) /* outflow */
	     ||( null == strcmp( tpp, "tf" )) /* fixed face temperature [C] */
             ||( null == strcmp( tpp, "tn" )) /* fixed node temperature [C] */
             ||( null == strcmp( tpp, "hc" )) /* fixed heat current */
             ||( null == strcmp( tpp, "rd" )) /* stefan-boltzmann rad */
             ||( null == strcmp( tpp, "sb" )) /* dito */
             ||( null == strcmp( tpp, "sf" )) /* surfce heat res [W/(K*m^2)] */
             ||( null == strcmp( tpp, "sc" ))) /* heat sources */
      {
/*............................................................................*/
# if DSC_HCRMDE == 0
         printf( "\n\n Message from function blc_boundary(*):" );
         printf( "\n Thermal boundary type '%s' remains inoperative", tp );
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
         if (( null == strncmp( tpp, "no_slip", TWO ))
           ||( null == strncmp( tpp, "free_slip", TWO ))
           ||( null == strncmp( tpp, "inflow", TWO ))
           ||( null == strncmp( tpp, "outflow", TWO ))
           ||( null == strcmp( tpp, "ns" ))
           ||( null == strcmp( tpp, "sl" ))
           ||( null == strcmp( tpp, "if" ))
           ||( null == strcmp( tpp, "uf" ))
           ||( null == strcmp( tpp, "of" )))
         {
            printf( "\n\n Message from function blc_boundary(*):" );
            printf( "\n Fluid flow boundary type '%s' remains "
              "inoperative", tp );
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
         blc_hcrbound( tpp, ii, nn, pp, ff, dd );
         return;
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      }
      else
      {
         printf( "\n Message from function blc_boundary(*):" );
         printf( "\n Unknown boundary type '%s' !!!", tp );
         printf( "\n Legal electric boundary types are:" );
         printf( "\n COAX, coax - coaxial [TEM mode]" ); 
         printf( "\n RECT, rect - rectangular waveguide [TEmn modes]" ); 
         printf( "\n CIRC, circ - circular waveguide [TEmn modes]" ); 
         printf( "\n ELL, ell - elliptic waveguide [c_TEmn modes]" ); 
         printf( "\n SR, sr - resitive [e.g. skin effect lossy] boundary" ); 
         printf( "\n SK, sk - [ dito -- same as preceeding ]" ); 
/*............................................................................*/
# if DSC_HCRMDE != 0
         printf( "\n Legal thermal boundary types are:" );
         printf( "\n TF, tf - fixed face temperature" ); 
         printf( "\n TN, tn - fixed node temperature" ); 
         printf( "\n HC, hc - heat current density imposed on face" ); 
         printf( "\n SB, sb - Stefan-Boltzmann [radiative] environmental "\
            "exchange" ); 
         printf( "\n RD, rd - [dito -- same as preceeding]" ); 
         printf( "\n SF, sf - fixed surface heat resistivity "\
            "[approximating convection, e.g.]" ); 
         printf( "\n SC, sc - thermal heat sources on [ or adjacent to ] "\
                 "skin effect\n          "\
                 "lossy [ electric ] boundary face [ of type sk ]" );
/*............................................................................*/
# if DSC_FLDMDE != 0 
         printf( "\n NS, ns - no-slip boundary conditions" ); 
         printf( "\n SL, sl - free slip boundary conditions" ); 
         printf( "\n IF, if - inflow bndry conditions ..." ); 
         printf( "\n OF, of - outflow bndry conditions ..." ); 
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
      };

     copy1:
      strncpy( trf.btp[jj], tpp, SHS_SIZE );

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

void blc_excite( char *tp, long ii, short nn, short pp, char ff, double uu )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      jj = null,
      kk = null,
      ll = null;

   static char
      tpp[SHS_SIZE] = {null};

/* prototypes: */

# if DSC_HCRMDE != 0
   void blc_hcrexcit( char *tp, long ii, short nn, short pp, 
      char ff, double rr );
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
/*............................................................................*/
/* string conversion to lower case: */

      strncpy( tpp, tp, SHS_SIZE );
      ll = strlen( tpp );

      kk = null;
      do
      {
         tpp[kk] = tolower( tpp[kk] );
      } while(( kk++ ) < ll );
/*............................................................................*/
/* excitation type */

      if (( null == strcmp( tpp, "e" )) /* electic field [port voltage; V] */
        ||( null == strcmp( tpp, "h" ))) /* magnetic field [port current; A*Z0]*/
      {
         goto copy1;
      }                                 
      else if (( null == strcmp( tpp, "tf" )) /* face temperature [Celsius]*/
             ||( null == strcmp( tpp, "tn" )) /* node temperature [Celsius]*/
                                         /* incident heat current density */
             ||( null == strncmp( tpp, "heat_current", THREE ))
             ||( null == strcmp( tpp, "hc" ))) /* dito */
      {
/*............................................................................*/
# if DSC_HCRMDE == 0
         printf( "\n\n Message from function blc_excite(*):" );
         printf( "\n Thermal excitation type '%s' remains inoperative", tp );
         printf( "\n due to uncompiled heat propagation mode !!!" );
         printf( "\n [ To implement heat propagation " );
         printf( "\n   set macro DSC_HCRMDE in file CONFIG.H to 1" );
         printf( "\n   and recompile program package." );
         printf( "\n Program continues [yet ignoring heat "
            "propagation modes]." );
         return;
# else /* DSC_HCRMDE != 0 */
         blc_hcrexcit( tpp, ii, nn, pp, ff, uu );
         return;
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      }
      else
      {
         printf( "\n Message from function blc_excite(*):" );
         printf( "\n Unknown excitation type '%s' !!!", tp );
         printf( "\n Legal field exciation types are:" );
         printf( "\n E, e - electric field" );
         printf( "\n H, h - magnetic field" );
/*............................................................................*/
# if DSC_HCRMDE != 0
         printf( "\n Legal thermal excitation types are:" );
         printf( "\n TF, tf - face temperature excitation" );
         printf( "\n TN, tn - node temperature excitation" );
         printf( "\n HC, hc - heat current density excitation" );
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
      };

     copy1:
      strncpy( trf.etp[jj], tpp, SHS_SIZE );

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

void blc_evaluate( char *tp, long ii, short nn, short pp, char pt, char *txt )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      jj = null,
      kk = null,
      ll = null;

   static char
      tpp[SHS_SIZE] = {null};

/* prototypes: */
/*............................................................................*/
# if DSC_HCRMDE != 0 /* heat current modes */
   void blc_hcrevlte( char *tp, long ii, short nn, short pp, 
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
/*............................................................................*/
/* string conversion to lower case: */

      strncpy( tpp, tp, SHS_SIZE );
      ll = strlen( tpp );

      kk = null;
      do
      {
         tpp[kk] = tolower( tpp[kk] );
      } while(( kk++ ) < ll );
/*............................................................................*/
/* evaluation type: */

      if (( null == strcmp( tpp, "e" )) /* electric field: port voltage [V] */
        ||( null == strcmp( tpp, "h" ))) /* magnetic field: port current [A*Z0] */
      {
         goto copy1;
      }                                      /* evaluate: */
      else if (( null == strcmp( tpp, "tf" )) /* face temperature [C] */
             ||( null == strcmp( tpp, "tn" )) /* node temperature [C] */
             ||( null == strcmp( tpp, "un" ))) /* nodal fluid velocity [m/sec] */
      {
/*............................................................................*/
# if DSC_HCRMDE == 0
         printf( "\n\n Message from function blc_evaluate(*):" );
         printf( "\n Thermal evaluation type '%s' remains inoperative", tp );
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
         if ( null == strcmp( tpp, "un" )) /* nodal fluid velocity [m/sec] */
         {
            printf( "\n\n Message from function blc_evaluate(*):" );
            printf( "\n Fluid velocity evaluation, type '%s', remains "
                    "inoperative", tp );
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
         blc_hcrevlte( tpp, ii, nn, pp, pt, txt );
         return;
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      }
      else
      {
         printf( "\n Message from function blc_evaluate(*):" );
         printf( "\n Unknown evaluation type '%s' !!!", tp );
         printf( "\n Legal field evaluation modes are:" );
         printf( "\n E, e - electric field" );
         printf( "\n H, h - magnetic field" );
/*............................................................................*/
# if DSC_HCRMDE != 0
         printf( "\n Legal thermal evaluation modes are:" );
         printf( "\n TF, tf - face temperature evaluation" );
         printf( "\n TN, tn - node temperature evaluation" );
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
      };

     copy1:

      if (( null < jj )\
        &&( null == strncmp( "reference", txt, THREE )))
      { /* set reference port at first place: */
         while( null < jj )
         {
            jj-- ;
            strncpy( trf.vtp[jj+ONE], trf.vtp[jj], SHS_SIZE );
            strncpy( trf.vtx[jj+ONE], trf.vtx[jj], STS_SIZE );
	    
            trf.vii[jj+ONE] = trf.vii[jj];
            trf.vnn[jj+ONE] = trf.vnn[jj];
            trf.vpp[jj+ONE] = trf.vpp[jj];
            trf.vpt[jj+ONE] = trf.vpt[jj];
         };
      };
         
      strncpy( trf.vtp[jj], tpp, SHS_SIZE ); /* type of evaluated port */
      strncpy( trf.vtx[jj], txt, STS_SIZE );

      trf.vii[jj] = (ii); /* initial evaluated face label */
      trf.vnn[jj] = (nn); /* number of periods */
      trf.vpp[jj] = (pp); /* period */
      trf.vpt[jj] = (pt); /* cell port index */

      trf.val++ ; /* count this evaluation section */
   }; /* end if blc.cov == null */
   return;
}
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0 /* heat current modes */
/*----------------------------------------------------------------------------*/
/* configure heat current and temperature boundaries: */
/*============================================================================*/

void blc_hcrbound( char *tp, long ii, short nn, short pp, char ff, double rr )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      jj = null,
      kk = null,
      ll = null;

   static char
      tpp[SHS_SIZE] = {null};
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.bhc; /* block counter */

      if( MOD_HCBNDRS <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many thermal boundary blocks defined !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_HCBNDRS.",
            ( long ) jj, ( long ) MOD_HCBNDRS );
         printf( "\n - Change macro MOD_HCBNDRS in model.c"
            "\n   in compliance with memory ressources.]\n" );
         exit( EXIT_FAILURE );
      };
/*............................................................................*/
/* string conversion to lower case: */

      strncpy( tpp, tp, SHS_SIZE );
      ll = strlen( tpp );

      kk = null;
      do
      {
         tpp[kk] = tolower( tpp[kk] );
      } while(( kk++ ) < ll );
/*............................................................................*/
/* thermal boundary type: */

      if (( null == strncmp( tpp, "heat_current", THREE )) /* fixed heat crr */
        ||( null == strncmp( tpp, "stefan_boltzmann", THREE )) 
                                          /* Stefan-Boltzmann heat radiation */
        ||( null == strncmp( tpp, "rad", THREE )) /* dito */
        ||( null == strncmp( tpp, "sources", THREE )) /* heat sources */
        ||( null == strncmp( tpp, "no_slip", TWO )) /* no-slip */
        ||( null == strncmp( tpp, "slip", TWO )) /* free slip */
        ||( null == strncmp( tpp, "free_slip", TWO )) /* dito */
        ||( null == strncmp( tpp, "inflow", TWO )) /* inflow */
        ||( null == strncmp( tpp, "outflow", TWO )) /* outflow */
        ||( null == strcmp( tpp, "ns" )) /* no-slip */
        ||( null == strcmp( tpp, "sl" )) /* free slip */
        ||( null == strcmp( tpp, "if" )) /* inflow */
        ||( null == strcmp( tpp, "uf" )) /* dito, as preceeding  */
        ||( null == strcmp( tpp, "of" )) /* outflow */
        ||( null == strcmp( tpp, "tf" )) /* fixed face temperature [C] */
        ||( null == strcmp( tpp, "tn" )) /* fixed node temperature [C] */
        ||( null == strcmp( tpp, "hc" )) /* fixed incident heat current */
        ||( null == strcmp( tpp, "rd" )) /* stefan-boltzmann radiation */
        ||( null == strcmp( tpp, "sb" )) /* dito, as preceeding */
        ||( null == strcmp( tpp, "sf" )) /* surface heat resistce [W/(K*m^2)] */
        ||( null == strcmp( tpp, "sc" ))) /* surface heat sources [ from skin */
                                          /* effect, e.g. ] */
      {
/*............................................................................*/
# if DSC_FLDMDE == 0 
         if (( null == strncmp( tpp, "no_slip", TWO ))
           ||( null == strncmp( tpp, "free_slip", TWO ))
           ||( null == strncmp( tpp, "inflow", TWO ))
           ||( null == strncmp( tpp, "outflow", TWO ))
           ||( null == strcmp( tpp, "ns" ))
           ||( null == strcmp( tpp, "sl" ))
           ||( null == strcmp( tpp, "if" ))
           ||( null == strcmp( tpp, "uf" ))
           ||( null == strcmp( tpp, "of" )))
         {
            printf( "\n\n Message from function blc_boundary(*):" );
            printf( "\n Fluid flow boundary type '%s' remains "
              "inoperative", tp );
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
         goto copy1;
      }
      else
      {
         printf( "\n Message from function blc_hcrbound(*):" );
         printf( "\n Unknown boundary type '%s' !!!", tp );
         printf( "\n Legal thermal boundary types are:" );
         printf( "\n TF, tf - fixed face temperature" ); 
         printf( "\n TN, tn - fixed node temperature" ); 
         printf( "\n HC, hc - heat current density imposed on face" ); 
         printf( "\n SB, sb - Stefan-Boltzmann heat radiative face" );
         printf( "\n RD, rd - same as preceeding" );
         printf( "\n SC, sc - heat sources on [or adjacent to] skin effect "
            "lossy face" );
/*............................................................................*/
# if DSC_FLDMDE != 0
         printf( "\n SL,sl - free slip boundary conditions on face" ); 
         printf( "\n IF,if - fluid inflow conditions ..." ); 
         printf( "\n OF,of - fluid outflow conditions ..." ); 
# endif
/*............................................................................*/
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
      };

     copy1:
      strncpy( trf.bct[jj], tpp, SHS_SIZE );

      trf.bci[jj] = (ii); /* initial boundary cell label */
      trf.bcn[jj] = (nn); /* number of periods */
      trf.bcp[jj] = (pp); /* period */
      trf.bcf[jj] = (ff); /* cell face index */
      trf.bcr[jj] = (rr); /* any real parameter, such as a heat current */
                          /* density, temperature, or anything else */
      trf.bhc++ ;
   }; /* end if blc.cov == null */
   return;
}
/*----------------------------------------------------------------------------*/
/* configure thermal excitation: */
/*============================================================================*/

void blc_hcrexcit( char *tp, long ii, short nn, short pp, char ff, double rr )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      jj = null,
      kk = null,
      ll = null;

   static char
      tpp[SHS_SIZE] = {null};
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.ehc;

      if( MOD_HCEXCTS <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many thermal blocks excited !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_HCEXCTS.",
            ( long ) jj, ( long ) MOD_HCEXCTS );
         printf( "\n - Change macro MOD_HCEXCTS in model.c"
            "\n   in compliance with memory resources.]\n" );
         exit( EXIT_FAILURE );
      };
/*............................................................................*/
/* string conversion to lower case: */

      strncpy( tpp, tp, SHS_SIZE );
      ll = strlen( tpp );

      kk = null;
      do
      {
         tpp[kk] = tolower( tpp[kk] );
      } while(( kk++ ) < ll );
/*............................................................................*/
/* thermal excitation type: */

      if (( null == strncmp( tpp, "heat_current", THREE ))
        ||( null == strcmp( tpp, "tf" )) /* face temperature [Celsius]*/
        ||( null == strcmp( tpp, "tn" )) /* node temperature [Celsius]*/
        ||( null == strcmp( tpp, "hc" ))) /* incident heat current */
      {
         goto copy1;
      }
      else
      {
         printf( "\n Message from function blc_hcrexcit(*):" );
         printf( "\n Unknown excitation type '%s' !!!", tp );
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
      }; 

     copy1:
      strncpy( trf.ect[jj], tpp, SHS_SIZE ); /* type of excitation port */

      trf.eci[jj] = (ii); /* initial excited cell label */
      trf.ecn[jj] = (nn); /* number of periods */
      trf.ecp[jj] = (pp); /* period */
      trf.ecf[jj] = (ff); /* cell face index */
      trf.ecr[jj] = (rr); /* any real parameter, such as a heat current */
                          /* density, temperature, or anything else */
      trf.ehc++ ;
   }; /* end if blc.cov == null */
   return;
}
/*----------------------------------------------------------------------------*/
/* configure thermal evaluation: */
/*============================================================================*/

void blc_hcrevlte( char *tp, long ii, short nn, short pp, char ff, char *txt )
{
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct transfer trf;
*/
/* declarations: */

   static short
      jj = null,
      kk = null,
      ll = null;

   static char
      tpp[SHS_SIZE] = {null};
/*----------------------------------------------------------------------------*/
   if( blc.cov == null )
   {
      jj = trf.vhc;

      if( MOD_HCEVALS <= jj )
      {
         printf( "\n\n Error mesage from function blocks(*) in model.c:" );
         printf( "\n Too many thermal blocks to be evaluated !!!" );
         printf( "\n [ The number %ld attains upper limit %ld"
            "\n = macro MOD_HCEVALS.",
            ( long ) jj, ( long ) MOD_HCEVALS );
         printf( "\n - Change macro MOD_HCEXCTS in model.c"
            "\n   in compliance with memory resources.]\n" );
         exit( EXIT_FAILURE );
      };
/*............................................................................*/
/* string conversion to lower case: */

      strncpy( tpp, tp, SHS_SIZE );
      ll = strlen( tpp );

      kk = null;
      do
      {
         tpp[kk] = tolower( tpp[kk] );
      } while(( kk++ ) < ll );
/*............................................................................*/
/* thermal evaluation type: */

      if (( null == strcmp( tpp, "tf" )) /* face temperature [C] */
        ||( null == strcmp( tpp, "tn" )) /* node temperature [C] */
        ||( null == strcmp( tpp, "un" ))) /* nodal fluid velocity [m/sec] */
      {
         goto copy1; 
      }
      else
      {
         printf( "\n Message from function blc_evaluate(*):" );
         printf( "\n Unknown evaluation type '%s' !!!", tp );
         printf( "\n Legal thermal evaluation modes are:" );
         printf( "\n TF, tf - face temperature evaluation" );
         printf( "\n TN, tn - node temperature evaluation" );
         printf( "\n [ Please correct DSC model "
            "and recompile source file model.c.]\n" );
         exit( EXIT_FAILURE );
      };

     copy1: 
      strncpy( trf.vct[jj], tpp, SHS_SIZE); /* type of evaluation port */
      strncpy( trf.vcx[jj], txt, STS_SIZE );

      trf.vci[jj] = (ii); /* initial evaluated cell label */
      trf.vcn[jj] = (nn); /* number of periods */
      trf.vcp[jj] = (pp); /* period */
      trf.vcf[jj] = (ff); /* cell face index, or component */

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
{ \
   if ( *option == 'c' ) \
      blc_boundary((TP),(II),(NN),(PP),(FF),(DD)); \
}
/*----------------------------------------------------------------------------*/
# define BLC_EXCITE(TP,II,NN,PP,FF,UU) \
{ \
   if ( *option == 'c' ) \
      blc_excite((TP),(II),(NN),(PP),(FF),(UU)); \
}
/*............................................................................*/
# define EXCITE(TP,II,NN,PP,FF,UU) \
{ \
   if ( *option == 'c' ) \
      blc_excite((TP),(II),(NN),(PP),(FF),(UU)); \
}
/*----------------------------------------------------------------------------*/
# define BLC_EVALUATE(TP,II,NN,PP,PT,TXT) \
{ \
   if ( *option == 'c' ) \
      blc_evaluate((TP),(II),(NN),(PP),(PT),(TXT)); \
}
/*............................................................................*/
# define EVALUATE(TP,II,NN,PP,PT,TXT) \
{ \
   if ( *option == 'c' ) \
      blc_evaluate((TP),(II),(NN),(PP),(PT),(TXT)); \
}
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
