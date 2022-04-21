/* [ File: panels.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function header panels.h                                            *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define BLC_DCLARE  0
# define BLC_INCLN   0
# define BLC_CENTER  0
# define BLC_SMEDIA  1
# define BLC_ZCONE   1
/*----------------------------------------------------------------------------*/
# ifndef DSC_ADJCLL
   # define DSC_ADJCLL 1
# endif
/*----------------------------------------------------------------------------*/
# ifdef MOD_DOMAINS
  # define BLC_DOMAINS MOD_DOMAINS
# else
  # define BLC_DOMAINS 50
# endif
# ifdef MOD_POINTS
  # define BLC_SUPPNTS MOD_POINTS
# else
  # define BLC_SUPPNTS 200
# endif
# ifdef MOD_BLOCKS
   # define BLC_BLOCKS MOD_BLOCKS
# else
   # define BLC_BLOCKS 100
# endif
# ifdef MOD_LAYERS
   # define BLC_LAYERS MOD_LAYERS
# else
   # define BLC_LAYERS 200
# endif
# ifndef BLC_P_BLIMTS
   # define BLC_P_BLIMTS 2
# endif
# ifdef PUZZ_POINTS
   # define BLC_PZZPNTS PUZZ_POINTS
# else
   # define BLC_PZZPNTS 501
# endif
/*----------------------------------------------------------------------------*/
typedef struct
{
   char 
     cov;

   double
     zz;

   double
      z[BLC_DOMAINS+ONE],
     dz[BLC_DOMAINS+ONE];

   long
      pm, pp, mm, mi, mf, pi, pf, ci, cf,
      dm[BLC_DOMAINS+ONE], /* the number of mesh cells per layer in domain[*] */
      dp[BLC_DOMAINS+ONE], /* the number of vertex points per  "  "         " */
      border[BLC_BLOCKS+ONE][FOUR][BLC_PZZPNTS];

   short 
      dmn, /* actual domain index on function calls cords(*) and blocks(*)    */
      lay, /* actual layer  "     "  "        "     "        "   "            */
      lbl, /* cell counter actualized by block generating functions qudrl(*)  */
      ab,  /* and trngl(*) - ab, bc: divisions transferred to these functions */
      bc,  /* [ blc.ab = ab_Teilung, blc.bc = bc_Teilung ]                    */
      m[BLC_DOMAINS+ONE],
   base[BLC_DOMAINS+ONE],
    spl[BLC_LAYERS+ONE],
  indnt[BLC_BLOCKS+ONE],
  leave[BLC_BLOCKS+ONE];

} BLOCSTR;
/*----------------------------------------------------------------------------*/
static BLOCSTR blc = {null};
/*----------------------------------------------------------------------------*/
static long 
   Grenzen_tmp[BLC_BLOCKS+ONE][FOUR][BLC_PZZPNTS] = {{{null}}};

short
   ab_Teilung = null,
   bc_Teilung = null;

FILE *blclmts; \

# define INTERFACE(AA,BB,CC,DD,EE,FF,GG,HH) \
  interface((AA),(BB),(CC),(DD),(EE),(FF),(GG),(HH))

# define CONNECT(AA,BB,CC,DD,EE,FF,GG,HH) \
  interface((AA),(BB),(CC),(DD),(EE),(FF),(GG),(HH))

# define connect(AA,BB,CC,DD,EE,FF,GG,HH) \
  interface((AA),(BB),(CC),(DD),(EE),(FF),(GG),(HH))

# define EINLESEN( BB,CC,DD,EE,FF,GG,HH) \
{ \
   if (( blc.ab == null ) \
     ||( blc.bc == null )) \
   { \
      if (( ab_Teilung != null ) \
        &&( bc_Teilung != null )) \
      { \
         blc.ab = ab_Teilung; \
         blc.ab = bc_Teilung; \
      } \
      else \
      { \
         fprintf( stderr, "\n\n Message from macro %s " \
                          "in function blocks(*):", "EINLESEN" ); \
         if ( blc.ab == null ) \
            fprintf( stderr, "\n illeagal blc.ab = 0 !!!" ); \
          else if ( blc.bc == null ) \
            fprintf( stderr, "\n illeagal blc.bc = 0 !!!" ); \
          fprintf( stderr, "\n [ Check parameters in block no. %d " \
            "of 'model.c' ]\n", blc.lbl ); \
          exit( EXIT_FAILURE ); \
      } \
   }; \
   interface(blc.lbl,(BB),(CC),(DD),(EE),(FF),(GG),(HH)); \
}

# define Einlesen( BB,CC,DD,EE,FF,GG,HH) \
   EINLESEN( BB,CC,DD,EE,FF,GG,HH) \

# define einlesen( BB,CC,DD,EE,FF,GG,HH) \
   EINLESEN(BB,CC,DD,EE,FF,GG,HH) \

# define DO_INTERFACE "interface(*)"
/*============================================================================*/

void interface ( short ActualBlock, signed char ActualSide,
                               short StartActual, short StopActual,
                 short FormerBlock, signed char FormerSide,
                               short StartFormer, signed char Sense )
{ 
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct quadrl qdl;
   extern struct triangle tgl;
*/
/* declarations: */

   static short 
      ii = null,
      jj = null;

/*----------------------------------------------------------------------------*/
   if( ActualBlock != blc.lbl )
   {
      printf( "\n\n Message from block generation "
         "function %s :", DO_INTERFACE );
      printf( "\n Block label error !!!" );
      printf( "\n Transferred block label %d differs "
         "from actual block label %d.", ActualBlock, blc.lbl );
      printf( "\n [ Check parameters in block no. %d "
         "of 'model.c' ]\n", blc.lbl );
      exit( EXIT_FAILURE );
   };

   if ( BLC_BLOCKS < ActualBlock )
   {
      printf( "\n\n Message from block generation "
         "function %s :", DO_INTERFACE );
      printf( "\n Too many blocks defined !!!");
      printf( "\n [ Blocklabel exceedes maximum %d "
         "= macro BLC_BLOCKS. ", BLC_BLOCKS );
      printf( "\n - Change macro in compliance with memory resources.]\n");
      exit( EXIT_FAILURE );
   };

   switch ( ActualSide )
   {  
     case 0: /* triangle or quadrangle: side ab */

      switch ( FormerSide )    
      { 
        case -1: /* e_wall */ 
         for ( ii=StartActual; ii<=StopActual; ii++ )
         {
            qdl.ab[ii] = 'e';
            tgl.ab[ii] = 'e';  
         };
         break;
                   
        case -2: /* m_wall */
         for ( ii=StartActual; ii<=StopActual; ii++ )
         {
            qdl.ab[ii] = 'm';
            tgl.ab[ii] = 'm';
         };
         break;
                
        default :  

         if ( blc.ab < StopActual )
         {
            printf( "\n\n Error in block generation "
               "function %s :", DO_INTERFACE );
            printf( "\n In block %d, side ab: Excessive point "
               "label ( > blc.ab )", blc.lbl );
            printf( "\n copied from neighboring block %d !!!", FormerBlock );
            printf( "\n [ Please correct block structure in calling "
               "function blocks(*).]\n" );
            exit( EXIT_FAILURE );
         }; 

         for ( ii=StartActual; ii<=StopActual; ii++ )
         {
            jj = Sense*(ii-StartActual)+StartFormer;

            if ( jj < null )
            {
               printf( "\n\n Error in block generation "
                  "function %s :", DO_INTERFACE );
               printf( "\n In block %d, side ab: Negative point "
                  "label", blc.lbl );
               printf( "\n copied from neighboring block %d !!!", FormerBlock );
               printf( "\n [ Please correct block structure in calling "
                  "function blocks(*).]\n" );
               exit( EXIT_FAILURE );
            };

            qdl.cab[ii] = blc.border[FormerBlock][FormerSide][jj];
            tgl.cab[ii] = blc.border[FormerBlock][FormerSide][jj];
         };
         break;  
      };
      break;

     case 1 : /* triangle or quadrangle: side bc*/

      switch ( FormerSide )
      { 
        case -1: /* e_wall */  
         for ( ii=StartActual; ii<=StopActual; ii++ )
         {
            qdl.bc[ii] = 'e';
            tgl.bc[ii] = 'e';
         }; 
         break;
                   
        case -2:  /* m_wall */
         for ( ii=StartActual; ii<=StopActual; ii++ )
         {
            qdl.bc[ii] = 'm';
            tgl.bc[ii] = 'm';
         };
         break;
                   
        default:  

         if (( blc.bc < StopActual )&&( 2*blc.ab < StopActual ))
         {
            printf( "\n\n Error in block generation "
               "function %s :", DO_INTERFACE );
            printf( "\n In block %d, side bc: Excessive point "
               "label ( > blc.bc ) ", blc.lbl );
            printf( "\n copied from neighboring block %d !!!", FormerBlock );
            printf( "\n [ Please correct block structure in calling "
               "function blocks(*).]\n" );
            exit( EXIT_FAILURE );
         };

         for ( ii=StartActual; ii<=StopActual; ii++ )
         {  
            jj = Sense*(ii-StartActual)+StartFormer;

            if ( jj < null )
            {
               printf( "\n\n Error in block generation "
                  "function %s :", DO_INTERFACE );
               printf( "\n In block %d, side bc: Negative point "
                  "label", blc.lbl );
               printf( "\n copied from neighboring block %d !!!", FormerBlock );
               printf( "\n [ Please correct block structure in calling "
                  "function blocks(*).]\n" ); 
               exit( EXIT_FAILURE );
            };

            qdl.cbc[ii] = blc.border[FormerBlock][FormerSide][jj];
            tgl.cbc[ii] = blc.border[FormerBlock][FormerSide][jj];
         };
         break;
      };   
      break;
 
     case 2: /* triange: side ac, or quadrangle: side ad */ 

      switch ( FormerSide )
      { 
        case -1: /* e_wall */  
         for ( ii=StartActual; ii<=StopActual; ii++ )
         {
            qdl.ad[ii] = 'e';
            tgl.ac[ii] = 'e';
         };
         break; 
                   
        case -2: /* m_wall */  
         for ( ii=StartActual; ii<=StopActual; ii++ )
         {
            qdl.ad[ii] = 'm';
            tgl.ac[ii] = 'm';
         };
         break;
                   
        default :  

         if (( blc.bc < StopActual )&&( blc.ab < StopActual ))
         {
            printf( "\n\n Error in block generation "
               "function %s :", DO_INTERFACE );
            printf( "\n In block %d, side ac ( triangle ) / "
               "side ad ( quadrangle ):", blc.lbl );
            printf( "\n Excessive point label ( blc.ab < ii or "
               "blc.bc < ii ) !!! " );
            printf( "\n copied from neighboring block %d.", FormerBlock );
            printf( "\n [ Please correct block structure in calling "
               "function blocks(*).]\n" );
            exit( EXIT_FAILURE );
         };

         for ( ii=StartActual; ii<=StopActual; ii++ )
         { 
            jj = Sense*(ii-StartActual)+StartFormer;

            if ( jj < null )
            {
               printf( "\n\n Error in block generation "
                  "function %s :", DO_INTERFACE );
               printf( "\n In block %d, side ac ( triangle ) / "
                  "side ad ( quadrangle ):", blc.lbl );
               printf( "\n Negative point label copied from neighboring "
                  "block %d !!!", FormerBlock );
               printf( "\n [ Please correct block structure in calling "
                  "function blocks(*).]\n" );
               exit( EXIT_FAILURE );
            };

            qdl.cad[ii] = blc.border[FormerBlock][FormerSide][jj];
            tgl.cac[ii] = blc.border[FormerBlock][FormerSide][jj];
         };
         break;
      };    
      break;

     case 3: /* quadrangle: side dc */

      switch ( FormerSide )
      { 
        case -1: /* e_wall */  

         for ( ii=StartActual; ii<=StopActual; ii++ )
            qdl.dc[ii] = 'e';

         break;
                   
        case -2: /* m_wall */
         for ( ii=StartActual; ii<=StopActual; ii++ )
            qdl.dc[ii] = 'm';
         break;
                   
        default :  

         if ( blc.ab < StopActual )
         {
            printf( "\n\n Error in block generation "
               "function %s :", DO_INTERFACE );
            printf( "\n In block %d, side dc: Excessive point label "
               "( > blc.ab )", blc.lbl );
            printf( "\n copied from neighboring block %d !!!", FormerBlock );
            printf( "\n [ Please correct block structure in calling "
               "function blocks(*).]\n" );
            exit( EXIT_FAILURE );
         }; 

         for ( ii=StartActual; ii<=StopActual; ii++ )
         { 
            jj = Sense*(ii-StartActual)+StartFormer;

            if ( jj < null )
            {
               printf( "\n\n Error in block generation "
                  "function %s :", DO_INTERFACE );
               printf( "\n In block %d, side dc: Negative point "
                  "label ", blc.lbl );
               printf( "\n copied from neighboring block %d !!!", FormerBlock );
               printf( "\n [ Please correct block structure in calling "
                  "function blocks(*).]\n" );
               exit( EXIT_FAILURE );
            };

            qdl.cdc[ii] = blc.border[FormerBlock][FormerSide][jj];
         };
         break;
      };    
      break;
   };
   return;
}
/*====================== end of function interface(*) ========================*/
# define QUADRANGLE(NN) \
  quadrangle((NN), blc.lay, option )

# define QUDRL(NN) \
  quadrangle((NN), blc.lay, option )

# define VIERECK(NN) \
  quadrangle((NN), blc.lay, option )
/*============================================================================*/

void quadrangle( short nn, short layer, char *option )
{ 
   static long
      ii = null;

   void readout( short Block, char Typ );

   if ( BLC_BLOCKS < nn )
   {
      printf( "\n\n Too many blocks defined in function blocks(*) !!!" );
      printf( "\n [ The maximum block number is %d "
         "= macro BLC_BLOCKS ", BLC_BLOCKS );
      printf( "\n   in program section 'model.c' ]" );
      printf( "\n - Change macro in compliance with "
         "memory resources.\n" );
      exit( EXIT_FAILURE );
   };
 
   if ( nn != blc.lbl )
   {
      printf( "\n\n Block label error: " );
      printf( "\n Label '%d' in BLOCK( %d ) ", blc.lbl, blc.lbl );
      printf( " differs from argument '%d' in QUDRL(*).", nn );
      printf( "\n [ Please correct in block structure function "
         "blocks(*).]\n" );
      exit( EXIT_FAILURE );
   };

   if ( blc.ab == null )
   {
      printf( "\n\n Warning from function 'quadrangle(*)':" );
      printf( "\n\n In quadranglular block %d:  blc.ab = null !!!", blc.lbl );

      if ( ab_Teilung != null )
      {
         printf( "\n  [ overrides, taking value ab_Teilung ]\n\n " );
         blc.ab = ab_Teilung;
      }
      else
         printf( "\n  [ overrides ]\n\n " );
   };

   if ( blc.bc == null )
   {
      printf( "\n\n Warning from macro 'QUDRL(*)' : " );
      printf( "\n\n In quadranglr.block %d:  blc.bc = null !!!", blc.lbl );

      if ( bc_Teilung != null )
      {
         printf( "\n  [ overrides, taking value bc_Teilung ] \n\n " );
         blc.bc = bc_Teilung;
      }
      else
         printf( "\n  [ overrides ]\n\n " );
   };

   if (( BLC_P_BLIMTS == ONE )&&( layer == null )&&( *option == 't' ))
   {
      printf( "\n quadranglr.block %d\t--->\tfirst cell: %ld",
         nn, ( lbl.m + ONE ));

      fprintf( blclmts, "\nquadranglr.block %d", nn );
      fprintf( blclmts, "\t--->\tfirst cell: %ld", ( lbl.m + ONE ));
   }
   else if (( BLC_P_BLIMTS == TWO )&&( *option == 't' ))
   {
      printf( "\n quadranglr.block %d\t--->\tfirst cell: %ld",
         nn, ( lbl.m + ONE ));

      if ( layer < blc.base[ blc.m[null] ] )
      {
         fprintf( blclmts, "\nquadranglr.block %d", nn );
         fprintf( blclmts, "\t--->\tfirst cell: %ld", ( lbl.m + ONE ));
      };
   };
 
   for( ii=null; ii<BLC_PZZPNTS; ii++ ) 
   { 
      tgl.cab[ii]    = null;
      tgl.cac[ii]    = null;
      tgl.cbc[2*ii]  = null;
      tgl.cbc[2*ii+1]= null;
      tgl.ab[ii]     = null;
      tgl.ac[ii]     = null;
      tgl.bc[2*ii]   = null;
      tgl.bc[2*ii+1] = null;
   };

   blc.ci = lbl.m + ONE; /* the first cell index of block nn */

/*............................................................................*/
   qudrl( blc.ab, blc.bc, option );   /*                                      */
/*..................................*/
 
   blc.cf = lbl.m; /* the last cell index of block nn */
 
   readout( nn, null );

   return;
}
/*===================== end of function quadrangle(*) ========================*/
# define TRIANGLE(NN) \
  triangle((NN), blc.lay, option )

# define TRNGL(NN) \
  triangle((NN), blc.lay, option )

# define DREIECK(NN) \
{ \
  blc.ab = ab_Teilung; \
  blc.bc = bc_Teilung; \
  triangle((NN), blc.lay, option ); \
}
/*============================================================================*/

void triangle( short nn, short layer, char *option )
{ 
   static long
      ii = null;

   void readout( short Block, char Typ );

   if ( BLC_BLOCKS < nn )
   {
      printf( "\n\n Too many blocks defined in function blocks(*) !!!" );
      printf( "\n [ The maximum block number is %d "
         "= macro BLC_BLOCKS ", BLC_BLOCKS );
      printf( "\n   in program section 'model.c' ]" );
      printf( "\n - Change macro in compliance with "
         "memory resources.\n" );
      exit( EXIT_FAILURE );
   };

   if ( nn != blc.lbl )
   {
      printf( "\n\n Block labeling error:" );
      printf( "\n Label '%d' in BLOCK( %d ) ", blc.lbl, blc.lbl );
      printf( "differs from argument '%d' in TRNGL(*) !!!", nn );
      printf( "\n [ Please correct in block structure function "
         "blocks(*).]\n" );
      exit( EXIT_FAILURE );
   };

   if ( blc.ab == null )
   {
      printf( "\n\n Warning from macro 'TRNGL(*)':" );
      printf( "\n\n In triangular block %d:  blc.ab = null !!!", blc.lbl );

      if ( ab_Teilung != null )
      {
         printf( "\n  [ overrides, taking value ab_Teilung ]\n\n " );
         blc.ab = ab_Teilung;
      }
      else
         printf( "\n  [ overrides ]\n\n " );
   };

   if (( BLC_P_BLIMTS == ONE )&&( layer == null )&&( *option == 't' ))
   {
      printf( "\n triangular block %d\t--->\tfirst cell: %ld",
         nn, ( lbl.m + ONE ));

      fprintf( blclmts, "\ntriangular block %d", nn );
      fprintf( blclmts, "\t--->\tfirst cell: %ld", ( lbl.m + ONE ));
   }
   else if (( BLC_P_BLIMTS == TWO )&&( *option == 't' ))
   {
      printf( "\n triangular block %d\t--->\tfirst cell: %ld",
         nn, ( lbl.m + ONE ));

      if ( layer < blc.base[ blc.m[null] ] )
      {
         fprintf( blclmts, "\ntriangular block %d", nn );
         fprintf( blclmts, "\t--->\tfirst cell: %ld", ( lbl.m + ONE ));
      };
   };

   for( ii=null; ii<BLC_PZZPNTS; ii++ )
   {
      qdl.cab[ii] = null;
      qdl.cad[ii] = null;
      qdl.cbc[ii] = null;
      qdl.cdc[ii] = null;
      qdl.ab[ii]  = null;
      qdl.ad[ii]  = null;
      qdl.bc[ii]  = null;
      qdl.dc[ii]  = null;
   };
 
   blc.ci = lbl.m + ONE; /* the first cell index of block nn */

/*............................................................................*/
   trngl( blc.ab, option );   /*                                              */
/*..........................*/

   blc.cf = lbl.m; /* the last cell index of block nn */

   readout( nn, ONE );

   return;
}
/*====================== end of function triangle(*) =========================*/

void readout ( short Block, char Typ )
{ 
/* allusions: */
/*
   extern BLOCSTR blc;
   extern struct quadrl qdl;
   extern struct triangle tgl;
*/
/* declarations: */

   static short
      ii = null;

   static signed char
      ab = null,
      bc = ONE, 
      ad = TWO,
      dc = THREE,
      ac = TWO;

/*----------------------------------------------------------------------------*/

   if ( BLC_BLOCKS < Block )
   {
      printf( "\n\n Message from block boundary transfer function "
         "'readout(*)':" );
      printf( "\n Too many blocks defined !!!" );
      printf( "\n [ Blocklabel exceedes maximum %d "
         "= macro BLC_BLOCKS. ", BLC_BLOCKS );
      printf( "\n - Change macro in compliance with memory resources.]\n" );
      exit( EXIT_FAILURE );
   };

   switch ( Typ )
   {  
     case 0: /* Viereck */
      for ( ii=null ; ii<=blc.ab ; ii++ )
      {
         blc.border[Block][ab][ii] = qdl.sab[ii];
         blc.border[Block][dc][ii] = qdl.sdc[ii];
      };
                   
      for ( ii=null; ii<=blc.bc; ii++ )
      {
         blc.border[Block][bc][ii] = qdl.sbc[ii];
         blc.border[Block][ad][ii] = qdl.sad[ii];
      };
      break;
      
     case 1: /* Dreieck */
      for ( ii=null; ii<=blc.ab; ii++ )
      {
         blc.border[Block][ab][ii] = tgl.sab[ii];
         blc.border[Block][ac][ii] = tgl.sac[ii];
      };
      for ( ii=null; ii<=2*blc.ab; ii++ )
      {
         blc.border[Block][bc][ii] = tgl.sbc[ii];
      };
      break;
   };

   return;
}
/*======================= end of function readout(*) =========================*/
# define POINT(AA,BB) \
  point((AA),(BB))

# define PUNKT(AA,BB) \
  point((AA),(BB))

# define Punkt(AA,BB) \
  point((AA),(BB))

# define punkt(AA,BB) \
  point((AA),(BB))
/*============================================================================*/

void point( signed char vertex, short pp )
{
/* allusions: */
/*
   extern struct quadrl qdl;
   extern struct triangle tgl;
   extern double pnt[BLC_SUPPNTS][THREE];
*/
/*----------------------------------------------------------------------------*/

   if ( BLC_SUPPNTS <= pp )
   {
      printf( "\n\n Message from point transfer function 'vertex(*)':" );
      printf( "\n Too many points defined !!!" );
      printf( "\n [ Point index 'pp' exceedes maximum %d "
         "= macro BLC_SUPPNTS . ", BLC_SUPPNTS );
      printf( "\n - Change macro in compliance with memory resources.]\n" );
      exit( EXIT_FAILURE );
   };

   switch ( vertex )
   { 
     case 0: /* point A */
      qdl.ax=pnt[pp][0];
      qdl.ay=pnt[pp][1];
      qdl.az=pnt[pp][2];
      tgl.ax=qdl.ax;
      tgl.ay=qdl.ay;
      tgl.az=qdl.az;
      break;
     
     case 1: /* point B */ 
      qdl.bx=pnt[pp][0];
      qdl.by=pnt[pp][1];
      qdl.bz=pnt[pp][2];
      tgl.bx=qdl.bx;
      tgl.by=qdl.by;
      tgl.bz=qdl.bz;
      break;

     case 2: /* point C */ 
      qdl.cx=pnt[pp][0];
      qdl.cy=pnt[pp][1];
      qdl.cz=pnt[pp][2];
      tgl.cx=qdl.cx;
      tgl.cy=qdl.cy;
      tgl.cz=qdl.cz;
      break;

     case 3 : /* point D */ 
      qdl.dx=pnt[pp][0];
      qdl.dy=pnt[pp][1];
      qdl.dz=pnt[pp][2];
      break;
   };
   return;
}
/*======================== end of function point(*) ==========================*/

/*******************************************************************************
*                                                                              *
*   Function indnt(*)                                                          *
*                                                                              *
*   On calling this function, cell counting restarts at block bloc             *
*   if blc.base[base] <= layer                                                 *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

int indnt( short bloc, short base, short layer, char *option )
{ 
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern BLOCSTR blc;
   extern struct labels lbl;
*/
/* declarations: */

   static short
      cntii = null,
      cntjj = null;

   static long
      cntll = null;

   cntii = bloc;
   if ( cntii != blc.lbl )
   {
      printf( "\n\n Error in block %d !!!", blc.lbl );
      printf( "\n [ block label differs from 1st. argument in indnt(%d,...)",
                                                                    cntii );
      printf( "\n   - please check function blocks(*).]\n" );
      exit( EXIT_FAILURE );
   };

   cntjj = blc.base[base];
   if ( cntjj <= layer )
   {
      if( null < blc.indnt[bloc] )
      {
         if ( lbl.mi == null )
            lbl.mi = blc.mm + ONE;

         cntll = lbl.mi;
         while( cntll <= lbl.mf )
         {
            if ( *option == 'c' ) /* option "coordinates" */
               ( spt->ppt->mpt->idx[cntll] ) = ONE;
            else /* if ( *option == 't'opology ) */
	    {
               cntii = null; do
/*............................................................................*/
# if DSC_ADJCLL == 0
               {
                  ( spt->tpt->mn[cntll][cntii] ) = null;
               } while (( ++cntii ) < FACES );
# elif DSC_ADJCLL == 1
               {
                  ( spt->tpt->mn[cntll][cntii] ) = null;
               } while (( ++cntii ) < PORTS );
# endif
/*............................................................................*/
            };
            cntll++;
         }; /* while cntll <= lbl.mf */

         lbl.m  = lbl.mi - ONE;
         lbl.mi = null;
         lbl.mf = lbl.m;
/*
*//* some constellations require this:
*//*
*//*     lbl.pi = null;
*/
      };
      if ( layer == null )
         blc.indnt[bloc] = TWO;
      else
         blc.indnt[bloc] = ONE;
   }
   else
      blc.indnt[bloc] = null;
   
   return blc.indnt[bloc];
}
/*======================== end of function indnt(*) ==========================*/
# define INDNT( BLC, BAS ) \
{ \
   indnt((BLC), (BAS), layer, option ); \
}
/******************************************************************************
*                                                                             *
*   Function leave(*)                                                         *
*                                                                             *
*   Returns ONE at a second passage at base layer1 and then at every passage  *
*   if layer1 <= layer <= layer2.                                             *
*   In all other cases the function returns null;                             *
*   Note that layer1, layer2 must be be base layers !                         *
*   The function can thus be used for leaving the block generation module     *
*   past any block, for layer1 <= layer <= layer2.                            *
*   To this end, the function is called after a quadrangle(*) or triangle(*)  *
*   function, and a 'goto terminal' statement  activated if the returned      * 
*   value is ONE.                                                             *
*                                                                             *
******************************************************************************/

/*============================================================================*/

int leave( short block,
           short layer1, short layer2, short layer )
{ 
   if ( block != blc.lbl )
   {
      printf( "\n\n Error in block %d !!!", blc.lbl );
      printf( "\n [ block label differs from 1st. argument "
         "in function leave(%d,...)  ", block );
      printf( "\n   - please check function blocks(*).]\n" );
      exit( EXIT_FAILURE );
   };

   if ((layer1 <= layer )&&( layer <= layer2 ))
   {
      switch( blc.leave[block] )
      {
        case 1:
         return ONE; /* this value should activate a 'goto terminal' */
                     /* statement in blocks(*) after the function call */
        default:
         blc.leave[block] = ONE;
         break;
      };
   }
   else
      blc.leave[block] = null;

   return null;
}
/*======================== end of function leave(*) ==========================*/
# define LEAVE(AA,BB,CC,DD,EE) { \
   (EE) = leave((AA),(BB),(CC),(DD)); }

/*******************************************************************************
*                                                                              *
*   Blocks clearing function  [ resets Grenzen_tmp[][][] to null ]             *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

short clblcs( short m , short n )
{
/* allusiobns: */
/*
   extern long 
      Grenzen_tmp[BLC_BLOCKS+ONE][FOUR][BLC_PZZPNTS];
*/
/* declarartions: */

   static short 
      ii, jj, kk;

   if ( m < null )
      m = null;
   
   if ( BLC_BLOCKS < n )
      n = BLC_BLOCKS;

   for ( ii = m; ii <= n; ii++ ) 
   {
      for ( jj = null; jj < FOUR; jj++ )
      {
         for ( kk = null; kk < BLC_PZZPNTS; kk++ )
         {
            Grenzen_tmp[ii][jj][kk] = null;
         };
      };
   };

   blc.ab = null;
   blc.bc = null;
   ab_Teilung = null;
   bc_Teilung = null;

   return n;
}
/*======================= end of function clblcs(*) ==========================*/
# if BLC_SMEDIA == 1

# define medium( LBL, EPS, MY, KE, KM ) \
    smedium((LBL), blc.ci, blc.cf, (EPS), (MY), (KE), (KM), 0, 0, option )

# define MEDIUM( LBL, EPS, MY, KE, KM ) \
    smedium((LBL), blc.ci, blc.cf, (EPS), (MY), (KE), (KM), 0, 0, option )

# define swtchmd( LBL, CI, CF, EPS, MY, KE, KM ) \
    smedium((LBL), (CI), (CF), (EPS), (MY), (KE), (KM), 0, 0, option )

# define SWTCHMD( LBL, CI, CF, EPS, MY, KE, KM ) \
    smedium((LBL), (CI), (CF), (EPS), (MY), (KE), (KM), 0, 0, option )

# define SMEDIUM( LBL, CI, CF, EPS, MY, KE, KM, KH, CV ) \
    smedium((LBL), (CI), (CF), (EPS), (MY), (KE), (KM), (KH), (CV), option )
/*============================================================================*/

void smedium( short medium_idx, long init_cell, long final_cell,
              double eps, double myr, double ke, double km,
              double kh, double cv, char *option )
{ 
/* allusions: */
/*
   extern FORMSTATE *spt;
*/
/* declarations: */

   static long
      ll = null;

   static short
      ii = null;

   static char
      pp = null;
/*----------------------------------------------------------------------------*/
   if( blc.cov != null )
      return;

   if (( *option == 'c' )
     &&(( spt->ppt->mpt->idx[null] ) < ( medium_idx - ONE )))
   {
      printf( "\n\n Error message from media switching function swtchmd(*):" );
      printf( "\n Illegal medium index %d !!!", medium_idx );
      printf( "\n Medium (indexed) %d has not yet been defined.", 
         ( medium_idx - ONE ));
      printf( "\n [ Media must be labelled in natural order by subsequent"
         " integers" );
      printf( "\n   - i.e. without omitting a number.]\n" );
      exit( EXIT_FAILURE );
   };
/*............................................................................*/
   if( null < init_cell )
   {
      for( ll=init_cell; ll<=final_cell; ll++ )
      {
         if( *option == 'c' ) /* option: 'c'oordinates */
            ( spt->ppt->mpt->idx[ll] ) = medium_idx;

         if( medium_idx == null ) /* trivial cell */
         { /* encapsulate cell with electric walls */
# if DSC_ADJCLL == 0
            for ( pp=null; pp<FACES; pp++ )
               ( spt->tpt->mn[ll][(int)pp] ) = -ONE;
# elif DSC_ADJCLL == 1
            for ( pp=null; pp<PORTS; pp++ )
               ( spt->tpt->mn[ll][(int)pp] ) = -ONE;
# endif
         };
      };
   };

/* [ mpt->idx[null]: the number ( = maximum index ) of yet defined media ] */

   if (( *option == 'c' )
     &&(( spt->ppt->mpt->idx[null] ) < medium_idx ))
   {
      ( spt->ppt->mpt->idx[null] ) = medium_idx;
      ( spt->ppt->mpt->tg[medium_idx] ) = ZERO;

      if ( fabs( eps ) < 1.e-277 )
         strcpy(( spt->ppt->mpt->type[medium_idx] ), "trivial_E" );
      else if ( fabs( myr ) < 1.e-277 )
         strcpy(( spt->ppt->mpt->type[medium_idx] ), "trivial_M" );
      else
         strcpy(( spt->ppt->mpt->type[medium_idx] ), "non-trv_EM" );

# if DSC_HCRMDE != 0
      ( spt->ppt->mpt->kh[medium_idx] ) = kh;
      ( spt->ppt->mpt->cv[medium_idx] ) = cv;

      if ( fabs( kh ) < 1.e-277 )
         strcat(( spt->ppt->mpt->type[medium_idx] ), "_trv_hcr" );
      else if ( fabs( cv ) < 1.e-277 )
         strcat(( spt->ppt->mpt->type[medium_idx] ), "_trv_hcr" );
# endif /* if DSC_HCRMDE != 0 */

      ii = null;
      do
      {
         ( spt->ppt->mpt->ep[medium_idx][ii] ) = eps;
         ( spt->ppt->mpt->my[medium_idx][ii] ) = myr;
         ( spt->ppt->mpt->ke[medium_idx][ii] ) = ke;
         ( spt->ppt->mpt->km[medium_idx][ii] ) = km;
         ( spt->ppt->mpt->ms[medium_idx][ii] ) = ZERO;
         ( spt->ppt->mpt->hg[medium_idx][ii] ) = ZERO;
      }  while(( ++ii ) < THREE );
      do
      {
         ( spt->ppt->mpt->ep[medium_idx][ii] ) = ZERO;
         ( spt->ppt->mpt->my[medium_idx][ii] ) = ZERO;
         ( spt->ppt->mpt->ke[medium_idx][ii] ) = ZERO;
         ( spt->ppt->mpt->km[medium_idx][ii] ) = ZERO;
      }  while(( ++ii ) < SIX );
   };
 
   return;
}
/*====================== end of function smedium(*) ==========================*/
# endif /* BLC_SMEDIA == 1 */

# if BLC_ZCONE == 1

# define DO_ZCONE "zcone(*)"
/*******************************************************************************
*                                                                              *
*   Function zcone(*)                                                          *
*                                                                              *
*   This function adjusts the z coordinates of all mesh cell points labelled   *
*   init to final in placing these points on the rotation surface ( cone )     *
*   around the vertical axis (i.e. in z direction ) that passes through given  *
*   points (z1,r1), (z2,r2).                                                   *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

void zcone( long init, long final,
            double z1, double r1, double z2, double r2 )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
*/
/* declarations: */

   static double 
      xx = ZERO,
      yy = ZERO,
      rr = ZERO,
      ss = ZERO,
      tt = ZERO;

   static long 
      mm = null,
      pp = null;

   static char 
      cc = null;

   double sqrt( double xx );

/*----------------------------------------------------------------------------*/

   if( fabs(r1-r2) < 1.e-77 )
   {
      printf( "\n\n Error message from function %s :", DO_ZCONE );
      printf( "\n Coincident radii r1 = r2 = %.12e !!!", r1 );
      printf( "\n The cone is degenerated to a cylinder." );
      printf( "\n The z-coordinates cannot be readjusted." );
      printf( "\n [ Check the calling function.]\n" );
      exit( EXIT_FAILURE );
   };

   mm = init;
   while( mm <= final )
   {
      cc = null; do
      {
         pp = ( spt->tpt->cm[mm][(int)cc] );
         xx = ( spt->ppt->cpt->c[pp][null] );
         yy = ( spt->ppt->cpt->c[pp][ONE] );
   
         rr = sqrt( xx*xx + yy*yy );
   
         ss = ( r1 - rr ) / ( r1 - r2 );
         tt = 1. - ss;

         ss *= z2;
         ss += ( tt*z1 );

         ( spt->ppt->cpt->c[pp][TWO] ) = ss;
      } while(( ++cc ) < FOUR ); 
      mm++ ;
   };
}
/*======================== end of function zcone(*) ==========================*/
# endif

# if BLC_INCLN == 1

# define DO_INCLN "incln(*)"
/*******************************************************************************
*                                                                              *
*   Function incln(*)                                                          *
*   [ in DANSE program package ]                                               *
*                                                                              *
*   In option  "v(u)"  ( u, v = x, y, z; u!=v ) this function redefines the    *
*   v-coordinate cor.c[jj][.] of all points jj, which occur as basis corner    *
*   points in cells mm, mm+1,...,nn  ( i.e.  jj = top.cm[mm][0,...,3],...      *
*   top.cm[nn][0,.,3] ) as                                                     *
*                                                                              *
*               v_mean + dv*( u[jj] - u_mean )/( u_max - u_min) ,              *
*                                                                              *
*   where u[jj] , u_mean and v_mean respectively denote the u-coordinate of    *
*   point jj, and the interval midpoints of the u- and v- coordinate ranges    *
*                                                                              *
*                   [ u_min, u_max ] and [ v_min, v_max ]                      *
*                                                                              *
*   of all points ii = top.cm[mm][0,...,3],...,top.cm[nn][0,...,3] .           *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

void incln( char *option, long mm , long nn , double dv )  
{
/* allusions: */
/*
   extern FORMSTATE *spt;
*/
/* declarations: */   

   static long 
      ii = null,
      jj = null;

   static double  
      u_min = ZERO,
      u_max = ZERO,
      u_mean = ZERO,
      v_min = ZERO,
      v_max = ZERO,
      v_mean = ZERO;

   static signed char
      idx1 = null,
      idx2 = null,
      kk = null;

/*----------------------------------------------------------------------------*/
   if ( mm == null )
   {
      printf("\n\n Warning from function '%s' : ", DO_INCLN );
      printf("\n\n Initial cell index null transfered. " );
      printf("\n [ overriding. ]\n");
      return;    
   };

   if ( nn < mm )
   {
      printf( "\n\n Error message from function %s :", DO_INCLN );
      printf( "\n Final cell index n=%ld precedes initial "
         "cell index m=%ld !!!", nn, mm );
      printf("\n [ Please correct function call:" );
      printf("\n   Transfer legal 2nd. and 3rd. "
         "arguments, m <= n, to '%s'.]\n", DO_INCLN );
      exit( EXIT_FAILURE );
   };

   if ( null == strncmp ( option, "x(y)", FOUR ))
   {
      idx1 = 1; 
      idx2 = 0;
   }
   else if ( null == strncmp ( option, "y(z)", FOUR ))
   {
      idx1 = 2;
      idx2 = 1;
   }
   else if ( null == strncmp ( option, "z(x)", FOUR ))
   {
      idx1 = 0;
      idx2 = 2;
   }
   else if ( null == strncmp ( option, "x(z)", FOUR ))
   {
      idx1 = 2;
      idx2 = 0;
   }
   else if ( null == strncmp ( option, "y(x)", FOUR ))
   {
      idx1 = 0;
      idx2 = 1;
   }
   else if ( null == strncmp ( option, "z(y)", FOUR ))
   {
      idx1 = 1;
      idx2 = 2;
   }
   else if ( null == strncmp ( option, "xy", TWO ))
   {
      idx1 = 1;
      idx2 = 0;
   }
   else if ( null == strncmp ( option, "yz", TWO ))
   {
      idx1 = 2;
      idx2 = 1;
   }
   else if ( null == strncmp ( option, "zx", TWO ))
   {
      idx1 = 0;
      idx2 = 2;
   }
   else if ( null == strncmp ( option, "xz", TWO ))
   {
      idx1 = 2;
      idx2 = 0;
   }
   else if ( null == strncmp ( option, "yx", TWO ))
   {
      idx1 = 0;
      idx2 = 1;
   }
   else if ( null == strncmp ( option, "zy", TWO ))
   {
      idx1 = 1;
      idx2 = 2;
   }
   else
   {
      printf( "\n\n Error message from function %s :", DO_INCLN );
      printf( "\n Unknown option %c%s%c !!!", 34, option, 34 );
      printf( "\n [ The legal options are: " ); 
      printf( "%cx(y)%c, %cy(z)%c, %cz(x)%c, %cz(y)%c, %cy(x)%c, %cx(z)%c.",
         34,34,34,34,34,34,34,34,34,34,34,34,34 ); 
      printf( "\n - Please correct function call:" );
      printf( "\n   Transfer legal option as 1st. "
         "argument to '%s'.]\n", DO_INCLN );
      exit( EXIT_FAILURE );
   };

   u_min = 1.e+277;
   v_min =   u_min;
   u_max = - u_min;
   v_max = - v_min; 

   ii = mm;
   while ( ii <= nn )
   {
      kk = null; do
      {
         jj = ( spt->tpt->cm[ii][kk] );
         if ( spt->ppt->cpt->c[jj][idx1] < u_min )
            u_min = ( spt->ppt->cpt->c[jj][idx1] ); 

         if ( spt->ppt->cpt->c[jj][idx1] > u_max )
            u_max = ( spt->ppt->cpt->c[jj][idx1] );

         if ( spt->ppt->cpt->c[jj][idx2] < v_min )
            v_min = ( spt->ppt->cpt->c[jj][idx2] );

         if ( spt->ppt->cpt->c[jj][idx2] > v_max )
            v_max = ( spt->ppt->cpt->c[jj][idx2] );

      } while (( ++kk ) < FOUR );
      ii++;
   };

   u_mean = .5*( u_min + u_max );
   v_mean = .5*( v_min + v_max ); 

   ii = mm;
   while ( ii <= nn )
   {
      kk = null; do
      {
         jj = ( spt->tpt->cm[ii][kk] );

         ( spt->ppt->cpt->c[jj][idx2] ) = v_mean +
            dv*(( spt->ppt->cpt->c[jj][idx1] ) - u_mean)/(u_max - u_min);

      } while (( ++kk ) < FOUR );
      ii++;
   };

   return;
}  
/*======================== end of function incln(*) ==========================*/
# endif /* BLC_INCLN == 1 */

# if BLC_CENTER == 1

/*============================================================================*/

double center( long cll, char fce, char crd )
{
   long 
      cp0 = null, 
      cp1 = null, 
      cp2 = null,
      cp3 = null;
   double 
      ss = ZERO;

   switch( fce )
   {
     case 0: 
      cp0=( spt->tpt->cm[cll][2] ); 
      cp1=( spt->tpt->cm[cll][6] );
      cp2=( spt->tpt->cm[cll][0] );
      cp3=( spt->tpt->cm[cll][4] );
      break;

     case 1:    
      cp0=( spt->tpt->cm[cll][5] );
      cp1=( spt->tpt->cm[cll][7] );
      cp2=( spt->tpt->cm[cll][1] );
      cp3=( spt->tpt->cm[cll][3] );
      break;

     case 2:
      cp0=( spt->tpt->cm[cll][4] );
      cp1=( spt->tpt->cm[cll][5] );
      cp2=( spt->tpt->cm[cll][0] );
      cp3=( spt->tpt->cm[cll][1] );
      break;

     case 3:
      cp0=( spt->tpt->cm[cll][3] );
      cp1=( spt->tpt->cm[cll][7] );
      cp2=( spt->tpt->cm[cll][2] );
      cp3=( spt->tpt->cm[cll][6] );
      break;

     case 4:
      cp0=( spt->tpt->cm[cll][1] );
      cp1=( spt->tpt->cm[cll][3] );
      cp2=( spt->tpt->cm[cll][0] );
      cp3=( spt->tpt->cm[cll][2] );
      break;

     case 5:
      cp0=( spt->tpt->cm[cll][7] );
      cp1=( spt->tpt->cm[cll][6] );
      cp2=( spt->tpt->cm[cll][4] );
      cp3=( spt->tpt->cm[cll][5] );
      break;
   };

   ss = .4*(( spt->ppt->cpt->c[cp0][crd-120] ) + \
            ( spt->ppt->cpt->c[cp1][crd-120] ) + \
            ( spt->ppt->cpt->c[cp2][crd-120] ) + \
            ( spt->ppt->cpt->c[cp3][crd-120] ));
   return ss;
}
/*======================== end of function center(*) =========================*/
# endif /* BLC_CENTER == 1 */
/*----------------------------------------------------------------------------*/
/* some useful macros called from blocks(*) at times: */
/*----------------------------------------------------------------------------*/
# define ZK_( BAS, TOP ) \
{ \
   ss = ( double )( blc.base[(TOP)] - layer ) / \
        ( blc.base[(TOP)] - blc.base[(BAS)] ); \
   tt =  1. - ss; \
}
/*----------------------------------------------------------------------------*/
# define BLOCK(NN) \
{ \
  blc.lbl = (NN); \
 \
  if ( BLC_BLOCKS < blc.lbl ) \
  { \
      printf( "\n\n Too many blocks defined in function blocks(*) !!!" ); \
      printf( "\n [ The maximum block number is %d" \
         " = macro BLC_BLOCKS ", BLC_BLOCKS ); \
      printf( "\n   in section 'syssmx.h' ]" ); \
      printf( "\n - Change macro in compliance wihth" \
         " memory resources.\n" ); \
      exit( EXIT_FAILURE ); \
  }; \
 \
  if (( null < BLC_P_BLIMTS )&&( blc.lay == null ) \
        &&( *option == 't')) \
  { \
   /*  printf( "\n block %d\t--->\tfirst cell: %ld", (NN), lbl.m + ONE ); */; \
  }; \
 \
  ab_Teilung = null; \
  bc_Teilung = null; \
}
/*----------------------------------------------------------------------------*/
# define BLC_LIMITS(NN) \
{ \
   cntii = ( short ) (NN); \
   if ( *option == 't' ) \
   { \
      if ( cntii == ( short ) 'o' ) \
      { \
         blclmts = fopen( blclog, "a" ); \
 \
         if (( blc.lay == null )&&( blc.cov == null )) \
         { \
            fclose( blclmts ); \
            blclmts = fopen( blclog, "w" ); \
 \
            fprintf( blclmts, "Block_limits" ); \
            fprintf( blclmts, "\n%s", ( spt->tpt->name )); \
            fprintf( blclmts, "\n%s\n\n", ( spt->tpt->text )); \
         }; \
 \
         if ((( BLC_P_BLIMTS == ONE )&&( blc.lay == null ))|| \
             ((BLC_P_BLIMTS == TWO )&&( blc.lay < blc.base[ blc.m[null] ] ))) \
         { \
            fprintf( blclmts, "layer %d", blc.lay ); \
            for ( cntii = null; cntii<blc.m[null]; cntii++ ) \
            { \
               if (( blc.lay<blc.base[cntii+ONE] )|| \
	           ( blc.lay == blc.base[blc.m[null]] )) \
               { \
                  if ( blc.lay == blc.base[cntii] ) \
                  { \
                     if ( blc.cov == ONE ) \
                     { \
                        fprintf( blclmts, ", base %d [ cover " \
                                          "- virtually removed ]", cntii ); \
                     } \
                     else if ( blc.cov == null ) \
                     { \
                        fprintf( blclmts, ", base %d", cntii ); \
                     }; \
                  };/* end if blc.lay == blc.base ... */ \
               };/* end if blc.lay < blc.base ... */ \
            }; \
            fprintf( blclmts, ":\n" ); \
         }; \
      } \
      else if ( cntii == ( short ) 'c' ) \
      { \
         if ((( BLC_P_BLIMTS == ONE )&&( blc.lay == null ))||\
             (( BLC_P_BLIMTS == TWO )&&( blc.lay < blc.base[ blc.m[null] ] ))) \
         { \
            fprintf( blclmts, "\n\n" ); \
         } \
         else if ( blc.lay == blc.base[ blc.m[null] ] ) \
         { \
            tmeptr = ( char * ) calloc( STS_SIZE, ONE ); \
 \
            nseconds = time( timer ); \
            tmeptr   = ctime( &nseconds ); \
 \
            blclmts = fopen( blclog, "a" ); \
 \
            fprintf( blclmts, "%s%s%s\n%s", "Block_limits_logfile '", \
                               blclog, "' created:", tmeptr ); \
            fprintf( blclmts, "\n\n" ); \
 \
            printf( "\n\n %s%s%s\n %.24s ", "Block limits logfile '", \
                    blclog, "' created: ", tmeptr ); \
         }; \
         fclose( blclmts ); \
      }; \
   }; \
}
/*----------------------------------------------------------------------------*/
# if BLC_DCLARE == 1

# define DECLARE( ) \
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
 \
   void readout( short, char ); \
   void vertex( signed char, short );

# endif /* end if BLC_DCLARE == 1 */
/*
# undef BLC_BLOCKS
*/
# undef BLC_DCLARE
# undef BLC_LAYERS
# undef BLC_DOMAINS
# undef BLC_SUPPNTS
# undef BLC_PZZPNTS
# undef BLC_INCLN
# undef BLC_ZCONE
# undef BLC_CENTER
# undef BLC_SMEDIA
# undef BLC_PZZPNTS
/************ End of DSC mesh ( block ) generation tools blocks.h *************/
