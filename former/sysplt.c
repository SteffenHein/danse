/* [ File: sysplt.c ] */
# define DO_SYSPLT "sysplt(*)"
/*******************************************************************************
*                                                                              *
*   GNUPLOT files creation function sysplot(*)                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   This function has originally been created by Ralf Gehring.                 *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 08, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
# include <time.h>           /* cf. time( ),ctime( ),asctime( ),localtime( )  */
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H"
/*----------------------------------------------------------------------------*/
# include "../former/FORMER.CONF" /* The FORMER.C configuration header */
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
   # include <ncurses.h>
# endif 
/*----------------------------------------------------------------------------*/
# ifndef DSC_ADJCLL
   # define DSC_ADJCLL 1 /* 0/1: define adjacent cells by faces [ports] */
# endif
/*----------------------------------------------------------------------------*/
# include "../former/formtp.h"
# include "../tools/txctyp.h"
/*----------------------------------------------------------------------------*/
# define SPL_INCLUDE 1
/*............................................................................*/
# if SPL_INCLUDE != 0
/* Include configuration header */
   # include "./sspltp.h"
# else /* SPL_INCLUDE == 0 */
/* Configure the following macros */
/*............................................................................*/
# ifdef USER_PATH
   # define SPL_PATH USER_PATH
# else
   # define SPL_PATH ""
# endif
/*............................................................................*/
# define SPL_PREFIX "gnu."
# define SPL_DISPLAY    1 /* 0<N: display plot sytem file creation messages */ 
                          /* 2:   display plot data file creation messages */
/*............................................................................*/
/* SPL_PLOTMEDIA - that maximum number of media */
# ifndef SPL_PLOTMEDIA
   # define SPL_PLOTMEDIA 10 /* K: plot that maximum number of media */
# endif /* not defined SPL_PLOTMEDIA */
/*............................................................................*/
/* SPL_LAYERS - that maximum number of plotted layers */
# ifndef SPL_LAYERS
   # define SPL_LAYERS 300
# endif /* not defined SPL_LAYERS */
/*............................................................................*/
/* SPL_PLOTBOTTM = [0] 1: [don't] plot walls on cell bottom */
# ifndef SPL_PLOTBOTTM 
   # define SPL_PLOTBOTTM 0
# endif /* not defined SPL_PLOTBOTTM */
/*............................................................................*/
/* SPL_PLOTTOPS = [0] 1: [don't] plot walls on cell tops */
# ifndef SPL_PLOTTOPS
   # define SPL_PLOTTOPS 0
# endif /* not defined SPL_PLOTTOPS */
/*............................................................................*/
/* SPL_PLOTTRIVL = [0] 1: [don't] plot trivial cells */
# ifndef SPL_PLOTTRIVL
   # define SPL_PLOTTRIVL 0
# endif /* not defined SPL_PLOTTRIVL */
/*............................................................................*/
/* SPL_PLOTEWLLS = 0/1/2: don't/plot/plot_dominant electric walls */
# ifndef SPL_PLOTEWLLS
   # define SPL_PLOTEWLLS 1
# endif /* not defined SPL_PLOTEWLLS */
/*............................................................................*/
/* SPL_PLOTMWLLS = 0/1/2: don't/plot/plot_dominant magnetic walls */
# ifndef SPL_PLOTMWLLS
   # define SPL_PLOTMWLLS 1
# endif /* not defined SPL_PLOTMWLLS */
/*............................................................................*/
/* SPL_RESCLE_XY_ = [0] 1: [don't] set uniform scales on xy axes */
# ifndef SPL_RESCLE_XY_
   # define SPL_RESCLE_XY_ 1 /* [0] 1: [don't] set uniform scales on xy axes */
# endif /* not defined SPL_RESCLE_XY_ */
/*............................................................................*/
/* SPL_RESCLE_XYZ = [0] 1: [don't] set uniform scales on all axes */
# ifndef SPL_RESCLE_XYZ
   # define SPL_RESCLE_XYZ 1
# endif /* not defined SPL_RESCLE_XYZ */
/*............................................................................*/
# define PLOT2_FILE "2p"
# define PLOT3_FILE "3p"
# define CELL2_FILE "2mesh"
# define CELL3_FILE "3mesh"

# define MED12_FILE "2trvl"
# define MED13_FILE "3trvl"
# define MED22_FILE "2epsr"
# define MED23_FILE "3epsr"

# define EWLL2_FILE "2ewll"
# define EWLL3_FILE "3ewll"
# define MWLL2_FILE "2mwll"
# define MWLL3_FILE "3mwll"
/*............................................................................*/
# ifndef ELECTRIC_WALL
   # define ELECTRIC_WALL -1
# endif
# ifndef MAGNETIC_WALL
   # define MAGNETIC_WALL -2
# endif
/*............................................................................*/
/* The structure [ type definition ] of DSC system plot function sysplt(*) */
typedef struct
{
   short
      rtn;

   short
      media,
      layer,
      toplr, /* the top layer, in general equals blc.base[blc.m[null]] */
      lay[SPL_LAYERS+ONE];

   long
      mi, mf, pm;

   char
      flbl[THREE],
      p2[SHS_SIZE],
      p3[SHS_SIZE],
      c2[SHS_SIZE],
      c3[SHS_SIZE];
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
   char
      ew2[SHS_SIZE],
      ew3[SHS_SIZE];
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
   char
      mw2[SHS_SIZE],
      mw3[SHS_SIZE];
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
   char
      md2[SPL_PLOTMEDIA][SHS_SIZE],
      md3[SPL_PLOTMEDIA][SHS_SIZE];

   short
      mdidx[SPL_PLOTMEDIA]; /* media index identifier: labels */
                            /* indices of media to be plotted */
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/

   TOPOLOGY
     *tpt;

   PARAMETERS
     *ppt;

} SYSPLT;
/*............................................................................*/
# endif /* SPL_INCLUDE == 0 */
/*----------------------------------------------------------------------------*/
static SYSPLT spl = {null};
/*----------------------------------------------------------------------------*/
# if SPL_DISPLAY != 0
   # include "dsptyp.h"
# endif
/*----------------------------------------------------------------------------*/
# define SPL_FPRNT2D(X,Y,Z) \
{ \
   fprintf((X), "%+.15e %+.15e%s", \
      ( ppt->cpt->c[tpt->cm[ii][(Y)]][0] ), \
      ( ppt->cpt->c[tpt->cm[ii][(Y)]][1] ), (Z)); \
}
/*............................................................................*/
# define SPL_FPRNT3D(X,Y,Z) \
{ \
   fprintf((X), "%+.15e %+.15e %+.15e%s", \
      ( ppt->cpt->c[tpt->cm[ii][(Y)]][0] ), \
      ( ppt->cpt->c[tpt->cm[ii][(Y)]][1] ), \
      ( ppt->cpt->c[tpt->cm[ii][(Y)]][2] ), (Z)); \
}
/*============================================================================*/

SYSPLT *\
sysplt ( SYSPLT *ssp )
{
/* allusions: */
/*
   extern FORMSTATE *spp;
*/
/* declarations: */

   static TOPOLOGY
     *tpt;

   static PARAMETERS
     *ppt;

   static SYSPLT
     *spp = &spl;

   static TXCNSL
     *csp;

   static double
      x_min   = ZERO,
      y_min   = ZERO,
      z_min   = ZERO,
      x_max   = ZERO,
      y_max   = ZERO,
      z_max   = ZERO,
      x_mean  = ZERO, 
      y_mean  = ZERO,
      z_mean  = ZERO,
      x_shift = ZERO,
      y_shift = ZERO,
      z_shift = ZERO,
      xx      = ZERO,
      yy      = ZERO,
      zz      = ZERO,
      rscale  =  .95;  /* range scaling factor */

   static long 
      ii = null,
      mm = null,
      nn = null;

   static short
      jj = null;

   static signed char
      kk  = null,
      ll  = null,
      opt = null,
      lbl[TEN] = {null};

   static char
      ptr[STS_SIZE] = {null},
      pptr[STS_SIZE] = {null},
      cptr[STS_SIZE] = {null},
/*............................................................................*/
# if DSC_ADJCLL == 0
      vtx[FACES][FOUR] = {{null}},
# elif DSC_ADJCLL == 1
      vtx[PORTS][FOUR] = {{null}},
# endif
/*............................................................................*/
     *index,
     *prefix = SPL_PREFIX,
     *timeptr;
/*............................................................................*/
# if SPL_DISPLAY != 0
   static DSPLAY
     *dsp;

   static char
     *timefrm = " created: %.24s";

   static long
      range = null;
# endif
/*............................................................................*/
   static char 
      **endp = null; 

   time_t nseconds = null;
   time_t   *timer = null;

/* streams: */
       
   static FILE 
      *pltfle2d,
      *pltfle3d,
      *cllfle2d,
      *cllfle3d;
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
   static char
      pp = null,
      mdpt[SPL_PLOTMEDIA][STS_SIZE]; 

   static FILE 
      *medfle2d[SPL_PLOTMEDIA],
      *medfle3d[SPL_PLOTMEDIA];
# endif
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
   static char
      emk = null,
      eptr[STS_SIZE] = {null};

   static FILE 
     *elwfle2d,
     *elwfle3d;
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
   static char
      mmk = null,
      mptr[STS_SIZE] = {null};

   static FILE 
      *mwlfle2d,
      *mwlfle3d;
# endif
/*............................................................................*/
/* prototypes: */

   time_t time( time_t *timer );

   char
     *lotos( long n, char length );

   TXCNSL
     *txcnsl( TXCNSL *csp );
/*............................................................................*/
# if SPL_DISPLAY != 0
   DSPLAY
     *dsplay( DSPLAY *dsp );
# endif
/*----------------------------- start of job ---------------------------------*/
/* initialize struct DSPLAY *dsp [ if used ]: */

# if SPL_DISPLAY != 0
      dsp = dsplay( null );
# endif
/*............................................................................*/
/* initialize struct SYSPLT spl: */

   if ( ssp == null )
   { 
      ( spp->media ) = SPL_PLOTMEDIA;

      ( spp->mi ) = null;
      ( spp->mf ) = null;
      ( spp->pm ) = null;
      ( spp->layer ) = null;
      ( spp->toplr ) = null;

      jj = null;
      while( jj < SHS_SIZE )
      {
         ( spp->p2[jj] ) = null;
         ( spp->p3[jj] ) = null;
         ( spp->c2[jj] ) = null;
         ( spp->c3[jj] ) = null;
         jj++ ;
      };
      strncpy(( spp->p2 ), PLOT2_FILE, SHS_SIZE );
      strncpy(( spp->p3 ), PLOT3_FILE, SHS_SIZE );
      strncpy(( spp->c2 ), CELL2_FILE, SHS_SIZE );
      strncpy(( spp->c3 ), CELL3_FILE, SHS_SIZE );

      jj = null;
      while( jj <= SPL_LAYERS )
      {
        ( spp->lay[jj] ) = null;
        jj++ ;
      };
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
      jj = null;
      while( jj < SHS_SIZE )
      {
         ( spp->ew2[jj] ) = null;
         ( spp->ew3[jj] ) = null;
         jj++ ;
      };
      strncpy(( spp->ew2 ), EWLL2_FILE, SHS_SIZE );
      strncpy(( spp->ew3 ), EWLL3_FILE, SHS_SIZE );
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
      jj = null;
      while( jj < SHS_SIZE )
      {
         ( spp->mw2[jj] ) = null;
         ( spp->mw3[jj] ) = null;
         jj++ ;
      };
      strncpy(( spp->mw2 ), MWLL2_FILE, SHS_SIZE );
      strncpy(( spp->mw3 ), MWLL3_FILE, SHS_SIZE );
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
      jj = null;
      while( jj < SPL_PLOTMEDIA )
      {
         kk = null;
	 while( kk < SHS_SIZE )
	 {
            ( spp->md2[jj][kk] ) = null;
            ( spp->md3[jj][kk] ) = null;
	    kk++ ;
         };
	 jj++ ;
      };
      jj = null;
/*............................................................................*/
# if SPL_PLOTTRIVL != 0
      strcpy(( spp->md2[jj] ), MED12_FILE );
      strcpy(( spp->md3[jj] ), MED13_FILE );
      ( spp->mdidx[jj] ) = 0; /* medium 0: trivial cells */
      ++jj;
# endif /* SPL_PLOTTRIVL != 0 */
/*............................................................................*/
      kk = TWO;
      while( jj < SPL_PLOTMEDIA )
      {
/* default: medium with permittivity Epsr[jj+kk-1] */
/* [ eventually transferred from calling program ] */ 

         ( spp->mdidx[jj] ) = ( kk++ );

         strcpy( ptr, MED22_FILE );
         strcat( ptr, lotos( jj, null ));
         strcpy(( spp->md2[jj] ), ptr );
         strcpy( ptr, MED23_FILE );
         strcat( ptr, lotos( jj, null ));
         strcpy(( spp->md3[jj] ), ptr );
         jj++;
      };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      ( spp->rtn ) = null;
      return spp;
   }
   else
   {
      spp = ssp;

      tpt = ( ssp->tpt );
      ppt = ( ssp->ppt );

      if ( SPL_LAYERS < ( ssp->toplr ))
         ( spp->toplr ) = SPL_LAYERS;

      if ( SPL_PLOTMEDIA < ( ssp->media ))
      {
         printf( "\n\n Too many media defined in DSC model !!!" );
         printf( "\n [ the number, media = %d, exceeds the maximum number %d",
            ( ssp->media ), SPL_PLOTMEDIA );
         printf( "\n   = macro SPL_PLOTMEDIA, defined in file sspltp.h ]" );
         printf( "\n - Program continued with displaying only allowed "
            "maximum." );
         
         ( spp->media ) = SPL_PLOTMEDIA;
      };
   };
/*............................................................................*/
/* set buffer length = null: */

   ii = setvbuf( stdin, null, _IONBF, null );
   ii = setvbuf( stdout, null, _IONBF, null ); 
/*............................................................................*/
/* memory allocations: */

   index   = ( char *) calloc( VSS_SIZE, ONE );
   timeptr = ( char *) calloc( STS_SIZE, ONE );
/*............................................................................*/
   lbl[0] = 0;  lbl[1] = 1;  lbl[2] = 3 ; lbl[3] = 2;  lbl[4] = 0;
   lbl[5] = 4;  lbl[6] = 5;  lbl[7] = 7 ; lbl[8] = 6;  lbl[9] = 4;
/*............................................................................*/
# if DSC_ADJCLL == 0
   vtx[0][0]  = 0; vtx[0][1]  = 2; vtx[0][2]  = 6; vtx[0][3]  = 4;
   vtx[1][0]  = 1; vtx[1][1]  = 3; vtx[1][2]  = 7; vtx[1][3]  = 5;
   vtx[2][0]  = 1; vtx[2][1]  = 0; vtx[2][2]  = 4; vtx[2][3]  = 5;
   vtx[3][0]  = 3; vtx[3][1]  = 2; vtx[3][2]  = 6; vtx[3][3]  = 7;
   vtx[4][0]  = 0; vtx[4][1]  = 1; vtx[4][2]  = 3; vtx[4][3]  = 2;
   vtx[5][0]  = 4; vtx[5][1]  = 5; vtx[5][2]  = 7; vtx[5][3]  = 6;
# elif DSC_ADJCLL == 1
   vtx[0][0]  = 3; vtx[0][1]  = 2; vtx[0][2]  = 6; vtx[0][3]  = 7;
   vtx[1][0]  = 4; vtx[1][1]  = 5; vtx[1][2]  = 7; vtx[1][3]  = 6;
   vtx[2][0]  = 1; vtx[2][1]  = 0; vtx[2][2]  = 4; vtx[2][3]  = 5;
   vtx[3][0]  = 0; vtx[3][1]  = 1; vtx[3][2]  = 3; vtx[3][3]  = 2;
   vtx[4][0]  = 4; vtx[4][1]  = 5; vtx[4][2]  = 7; vtx[4][3]  = 6;
   vtx[5][0]  = 1; vtx[5][1]  = 3; vtx[5][2]  = 7; vtx[5][3]  = 5;
   vtx[6][0]  = 0; vtx[6][1]  = 1; vtx[6][2]  = 3; vtx[6][3]  = 2;
   vtx[7][0]  = 0; vtx[7][1]  = 2; vtx[7][2]  = 6; vtx[7][3]  = 4;
   vtx[8][0]  = 1; vtx[8][1]  = 3; vtx[8][2]  = 7; vtx[8][3]  = 5;
   vtx[9][0]  = 3; vtx[9][1]  = 2; vtx[9][2]  = 6; vtx[9][3]  = 7;
   vtx[10][0] = 0; vtx[10][1] = 2; vtx[10][2] = 6; vtx[10][3] = 4;
   vtx[11][0] = 1; vtx[11][1] = 0; vtx[11][2] = 4; vtx[11][3] = 5;
# endif
/*............................................................................*/
   x_shift = ZERO;
   y_shift = ZERO;
   z_shift = ZERO;
   rscale  = .95;
/*............................................................................*/

   if (( ssp->layer ) < null )
   {
/*............................................................................*/
      csp = txcnsl( null );    /* initialize text console                     */
/*...........................*/
      ( csp->items ) = 3;
      ( csp->dfopt ) = 0;
      strcpy(( csp->envmt ), "SYSPLOT" );
      strcpy(( csp->cmmnt ), "DSC system plot" );
      strcpy(( csp->tasks ), "Create:" );
      strcpy(( csp->mline[1] ), "2D - plot" );
      strcpy(( csp->mline[2] ), "3D - plot" );
      strcpy(( csp->mline[3] ), "2&3D - plot" );
      strcpy(( csp->mline[4] ), "Plot single layers" );
      strcpy(( csp->escpe ), "Continue/escape: enter 0" );
/*............................................................................*/
      csp = txcnsl( csp );     /* input on text console */
/*...........................*/
      opt = ( csp->option );

      if (( null < opt )&&( opt < FOUR ))
      {
         printf( "\n" );
         strcpy(( csp->rqstr ), "Rescale coordinate system ?" );
         strcpy(( csp->dfstr ), "n" );
/*............................................................................*/
         csp = txcnsl( csp );        /* input on text console */
/*.................................*/
         strcpy( ptr, ( csp->instr )); 

         if ( *ptr == 'y' )
         {
           scaling:

            strcpy(( csp->rqdbl ), "Please enter scaling factor"
                                " [ ZERO < factor ]" );
            ( csp->dfdbl ) = 1.;
/*............................................................................*/
            csp = txcnsl( csp );        /* input on text console              */
/*....................................*/
            rscale = ( csp->indbl );

            if ( rscale < 1.e-277 )
               goto scaling;

            strcpy(( csp->rqdbl ), "Enter x-shift [-1 <=shift <= 1.]" );
            ( csp->dfdbl ) = 0.;
/*............................................................................*/
            csp = txcnsl( csp );        /* input on text console              */
/*....................................*/
            x_shift = ( csp->indbl );

            strcpy(( csp->rqdbl ), "Enter y-shift [-1 <=shift <= 1.]" );
            ( csp->dfdbl ) = 0.;
/*............................................................................*/
            csp = txcnsl( csp );        /* input on text console              */
/*....................................*/
            y_shift = ( csp->indbl );

            if (( ONE < opt )&&( opt < FOUR ))
            {
               strcpy(( csp->rqdbl ), "Enter z-shift [-1 <=shift <= 1.]" );
               ( csp->dfdbl ) = 0.;
/*............................................................................*/
               csp = txcnsl( csp );     /* input on text console              */
/*....................................*/
               z_shift = ( csp->indbl );
            };
         };
      }
      else if ( opt == FOUR )
      {
         printf( "\n Please enter indices of layers to be plotted"
                 " [ Escape/end: enter -ONE ]:\n" );

         jj = ( ssp->toplr );
         ii = null;
         do
         {
            ii++ ;
            printf( "\t\t\t\t\t\t\t ....................? " );
            printf( "\r >----> enter %.3ld. layer >---"
               "----------------------------------> :", ii );
            scanf( "%s", ptr );
            ( spp->lay[ii] ) = strtol( ptr, endp, DEC );
         } while (( *ptr != '-' )&&(( spp->lay[ii] ) < jj )&&( ii <= jj ));
         if ( *ptr == '-' )
           ii--;
         ( spp->lay[null] ) = ii;
      };
   }
   else /* if ( null <= ( ssp->layer )) */
   {
/*............................................................................*/
      index = lotos(( ssp->layer), null );  /* converts ( int ) layer */
/*........................................*//* into ASCII string */ 
      opt = FOUR;
   };
  
   if (( opt == TWO )
     ||( opt == THREE )
     ||( opt == FOUR ))
   {
/*............................................................................*/
/* 3D-scale: */

      if ( opt == FOUR )
      {
         ll = EIGHT;
         mm = ( ssp->mi );
         nn = mm + ( ssp->pm ) - ONE;
      }
      else
      {
         ll = FOUR;

         mm = ( ssp->mi );
         nn = ( ssp->mf );
      };

      x_min = 1.e277;
      y_min = 1.e277;
      z_min = 1.e277;

      x_max = -x_min;
      y_max = -y_min;
      z_max = -z_min;
     
      for ( ii=mm; ii<=nn; ii++ )
      {
         for ( jj=null; jj<ll; jj++ )
         {   
            xx = ( ppt->cpt->c[tpt->cm[ii][jj]][0] );
            yy = ( ppt->cpt->c[tpt->cm[ii][jj]][1] );
            zz = ( ppt->cpt->c[tpt->cm[ii][jj]][2] );                 
               
            if( x_max < xx ) x_max = xx;
            if( y_max < yy ) y_max = yy;
            if( z_max < zz ) z_max = zz;
                                 
            if( xx < x_min ) x_min = xx;
            if( yy < y_min ) y_min = yy;
            if( zz < z_min ) z_min = zz;
         };                                             
      };

      xx = ( x_max - x_min )/2.;
      yy = ( y_max - y_min )/2.;
      zz = ( z_max - z_min )/2.;
/*............................................................................*/
# if SPL_RESCLE_XYZ == 1
      if( xx < yy ) xx = yy;
      if( xx < zz ) xx = zz;
      yy = xx;
      zz = xx;
# endif
/*............................................................................*/
      x_mean = ( x_min + x_max )/2.;
      x_min = x_mean + ( x_shift - 1./rscale )*xx;
      x_max = x_mean + ( x_shift + 1./rscale )*xx;

      y_mean = ( y_min + y_max )/2.;
      y_min = y_mean + ( y_shift - 1./rscale )*yy;
      y_max = y_mean + ( y_shift + 1./rscale )*yy;
      
      z_mean = ( z_min + z_max )/2.;
      z_min = z_mean + ( z_shift - 1./rscale )*zz;
      z_max = z_mean + ( z_shift + 1./rscale )*zz;
/*............................................................................*/
/* 3D-plot system file */

      strcpy( pptr, prefix );
      strcat( pptr, ssp->p3 );

      if ( opt == FOUR )
         strcat( pptr, index );
      else
         strncat( pptr, ( ssp->flbl ), sizeof( &( ssp->flbl )));

      pltfle3d = fopen( pptr, "w" );

      if ( pltfle3d == null )
      {
         printf( "\n\n Message from function %s :", DO_SYSPLT );
         printf( "\n Unknown error on opening file '%s' !!!", pptr );
         printf( "\n [ overrides and returns to calling program ].\n" );
         ( spp->rtn ) = ONE;
         return spp;
      };
/*
      printf( "\n opened: 3D-plot system file '%s'", pptr );
*/
/*............................................................................*/
      strcpy( cptr, prefix );
      strcat( cptr, ( ssp->c3 ));
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
      strcpy( eptr, prefix );
      strcat( eptr, ( ssp->ew3 ));
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
      strcpy( mptr, prefix );
      strcat( mptr, ( ssp->mw3 ));
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
      jj = null;
      while( jj < ( ssp->media ))
      {
         strcpy( mdpt[jj], prefix );
         strcat( mdpt[jj], ( ssp->md3[jj] ));
         jj++;
      };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      if ( opt == FOUR )
      {
         ll = FIVE;

         strcat( cptr, index );
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         strcat( eptr, index );
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         strcat( mptr, index );
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
         jj = null;
         while( jj < ( ssp->media ))
         {
            strcat( mdpt[jj], index );
            jj++;
         };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      }
      else
      {
         ll = TEN;

         strncat( cptr, ( ssp->flbl ), sizeof( &( ssp->flbl )));
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         strncat( eptr, ( ssp->flbl ), sizeof( &( ssp->flbl )));
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         strncat( mptr, ( ssp->flbl ), sizeof( &( ssp->flbl )));
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
         jj = null;
         while( jj < ( ssp->media ))
         {
            strncat( mdpt[jj], ( ssp->flbl ), sizeof( &( ssp->flbl )));
            jj++;
         };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      };

      fprintf( pltfle3d, "set parametric\n" );
/*
      if ( opt != FOUR )
         fprintf( pltfle3d, "set hidden3d\n" );
*/
      fprintf( pltfle3d, "set nokey\n" );
      fprintf( pltfle3d, "set style data lines\n" );
      fprintf( pltfle3d, "set title '3D-plot of %s %s",
         ( ppt->name ), ( ppt->text ));

      if ( opt == FOUR )
         fprintf( pltfle3d, " [ layer %s ]'\n", index ); 
      else
         fprintf( pltfle3d, "'\n" );

      fprintf( pltfle3d, "set xrange [%.15e:%.15e]\n", x_min, x_max );
      fprintf( pltfle3d, "set yrange [%.15e:%.15e]\n", y_min, y_max );
      fprintf( pltfle3d, "set zrange [%.15e:%.15e]\n", z_min, z_max );
      fprintf( pltfle3d, "set xlabel 'x/m'\n" );
      fprintf( pltfle3d, "set ylabel 'y/m'\n" );
      fprintf( pltfle3d, "set zlabel 'z/m'\n" );       
      fprintf( pltfle3d, "set nogrid\n" );
      fprintf( pltfle3d, "set border\n" );
      fprintf( pltfle3d, "xrot=60 \nzrot=0\n" );
      fprintf( pltfle3d, "splot %c\n", 92 ); /* dec 92 = ASCII '\' */

      fprintf( pltfle3d, "'%s%s' with lines", SPL_PATH, cptr );
/*............................................................................*/
# if SPL_PLOTEWLLS == 1
      fprintf( pltfle3d, ",%c\n'%s%s' with lines", 92, SPL_PATH, eptr );
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS == 1
      fprintf( pltfle3d, ",%c\n'%s%s' with lines", 92, SPL_PATH, mptr );
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
      jj = null;
      while( jj < ( ssp->media ))
      {
         fprintf( pltfle3d, ",%c\n'%s%s' with lines", 92, SPL_PATH, mdpt[jj] );
         jj++;
      };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
# if SPL_PLOTMWLLS == 2
      fprintf( pltfle3d, ",%c\n'%s%s' with lines", 92, SPL_PATH, mptr );
# endif
/*............................................................................*/
# if SPL_PLOTEWLLS == 2
      fprintf( pltfle3d, ",%c\n'%s%s' with lines", 92, SPL_PATH, eptr );
# endif
/*............................................................................*/
      fprintf( pltfle3d, "\npause -1 '[ hit return to continue ]'\n" );

      nseconds = time( timer );
      timeptr = ctime( &nseconds );

      fprintf( pltfle3d, "\n# 3D-plot system file '%s'", pptr );
      fprintf( pltfle3d, "\n# created: %.24s", timeptr );
      fprintf( pltfle3d, "\n# EOF" );

      fclose( pltfle3d ); 

# if SPL_DISPLAY != 0
      if ( opt == FOUR )
      {
         printf( "\n 3D-plot system file %s", pptr );
         printf( timefrm, timeptr );
      }
      else
      {
         strcpy( ptr, "\n 3D-plot system file '" );
         strcat( ptr, pptr );
         strcat( ptr, "' created: " );
         strncat( ptr, timeptr, 24 );
         range = nn  - mm + ONE;
/*............................................................................*/
/* start dsplay(*) function: */

         ( dsp->option ) = 's'; /* display 's'tart message */
         strcpy(( dsp->messge ), ptr );
         dsplay( dsp );

         ( dsp->option ) = 'm'; /* 'm'essage under runnig cursor */
         strcpy(( dsp->messge ), "[Generating plot files;"
            " please wait a moment]" );
         dsplay( dsp );

         ( dsp->option ) = 'c'; /* option: running 'c'ursor */
         ( dsp->range ) = range;
/*............................................................................*/
      };
# endif
/*............................................................................*/
/* 3D-plot points file: */ 

      cllfle3d = fopen( cptr,"w");

      if ( cllfle3d == null )
      {
         printf( "\n\n Message from function %s :", DO_SYSPLT );
         printf( "\n Unknown error on opening file '%s' !!!", cptr );
         printf( "\n [ overrides and returns to calling program ].\n " );
         spp->rtn = ONE;
         return spp;
      };
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
      elwfle3d = fopen( eptr, "w");

      if ( elwfle3d == null )
      {
         printf( "\n\n Message from function %s :", DO_SYSPLT );
         printf( "\n Unknown error on opening file '%s' !!!", eptr );
         printf( "\n [ overrides and returns to calling program ].\n " );
         spp->rtn = ONE;
         return spp;
      };
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
      mwlfle3d = fopen( mptr, "w");

      if ( mwlfle3d == null )
      {
         printf( "\n\n Message from function %s :", DO_SYSPLT );
         printf( "\n Unknown error on opening file '%s' !!!", mptr );
         printf( "\n [ overrides and returns to calling program ].\n " );
         spp->rtn = ONE;
         return spp;
      };
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
      jj = null;
      while( jj < ( ssp->media ))
      {
         medfle3d[jj] = fopen( mdpt[jj], "w" );

         if ( medfle3d[jj] == null )
         {
            printf( "\n\n Message from function %s :", DO_SYSPLT );
            printf( "\n Unknown error on opening file '%s' !!!", mdpt[jj] );
            printf( "\n [ overrides and returns to calling program ].\n " );
            spp->rtn = ONE;
            return spp;
         };
         jj++;
      };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      for ( ii=mm; ii<=nn; ii++ )
      {
         fprintf( cllfle3d, "# cell %ld\n", ii );
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         emk = null;
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         mmk = null;
# endif
/*............................................................................*/
         jj = null;
         while ( jj < ll )
         {
            kk = lbl[jj];
            fprintf( cllfle3d, "%+.15e %+.15e %+.15e\n",
               ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
               ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
               ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
            jj++;
         };
         fprintf( cllfle3d, "\n" );
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
         jj = null;
         while( jj < ( ssp->media ))
         {
            if (( ppt->mpt->idx[ii] ) == ( ssp->mdidx[jj] ))
            {
               fprintf( medfle3d[jj], "# cell %ld\n", ii );
               pp = null;
               while ( pp < ll )
               {
                  kk = lbl[(int)pp];
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  pp++;
               };

               if ( opt == FOUR )
               {
                  kk = THREE;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  kk = ONE;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  kk = TWO;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
               }
               else /* opt != FOUR */
               {
                  kk = THREE;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  kk = FIVE;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  kk = TWO;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  kk = SIX;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  kk = ONE;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  kk = SEVEN;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  kk = null;
                  fprintf( medfle3d[jj], "%+.15e %+.15e %+.15e\n\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
               };
            }; /* end if ... = mediaindex */
            jj++ ;
         }; /* end while */
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
# if DSC_ADJCLL == 0
         jj = FACES;
# elif DSC_ADJCLL == 1
         jj = PORTS;
# endif
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         jj = null;
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         jj = null;
# endif
/*............................................................................*/
# if DSC_ADJCLL == 0
         while ( jj < FACES )
         {
            switch( jj )
            {
/*............................................................................*/
# if SPL_PLOTBOTTM == 0
              case 4:
               break;
# endif
/*............................................................................*/
# if SPL_PLOTTOPS == 0
              case 5:
               break;
# endif
/*............................................................................*/
# elif DSC_ADJCLL == 1
         while ( jj < PORTS )
         {
            switch( jj )
            {
/*............................................................................*/
# if SPL_PLOTBOTTM == 0
              case 3:
              case 6:
               break;
# endif
/*............................................................................*/
# if SPL_PLOTTOPS == 0
              case 1:
              case 4:
               break;
# endif
/*............................................................................*/
# endif /* end if DSC_ADJCLL == 1 */

              default:

/*............................................................................*/
# if SPL_PLOTEWLLS != 0
               if (( tpt->mn[ii][jj] ) == ELECTRIC_WALL )
               {
                  if ( emk == null )
                  {
                     fprintf( elwfle3d, "# cell %ld\n", ii );
                     emk = ONE;
                  };

                  kk = vtx[jj][0];
                  fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                  kk = vtx[jj][1];
                  fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                  if ( opt != FOUR )
                  {
                     kk = vtx[jj][2];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                     kk = vtx[jj][3];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
               
                     kk = vtx[jj][0];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  }
/*............................................................................*/
# if SPL_PLOTBOTTM == 1
/*............................................................................*/
# if DSC_ADJCLL == 0 
                  else if ( jj == FOUR )
                  {
# elif DSC_ADJCLL == 1 
                  else if (( jj == 3 )
                         ||( jj == 6 ))
                  {
# endif
/*............................................................................*/
                     kk = vtx[jj][2];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                     kk = vtx[jj][3];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
               
                     kk = vtx[jj][0];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                     kk = vtx[jj][2];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                     kk = vtx[jj][1];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
               
                     kk = vtx[jj][3];
                     fprintf( elwfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  };
# endif /* SPL_PLOTBOTTM == 1 */
/*............................................................................*/
                  fprintf( elwfle3d, "\n" );
               };
# endif /* SPL_PLOTEWLLS != 0 */
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
               if (( tpt->mn[ii][jj] ) == MAGNETIC_WALL )
               {
                  if ( mmk == null )
                  {
                     fprintf( mwlfle3d, "# cell %ld\n", ii );
                     mmk = ONE;
                  };

                  kk = vtx[jj][0];
                  fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                  kk = vtx[jj][1];
                  fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                  if ( opt != FOUR )
                  {
                     kk = vtx[jj][2];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                     kk = vtx[jj][3];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
               
                     kk = vtx[jj][0];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  } 
# if SPL_PLOTBOTTM == 1
/*............................................................................*/
# if DSC_ADJCLL == 0 
                  else if ( jj == FOUR )
                  {
# elif DSC_ADJCLL == 1 
                  else if (( jj == 3 )
                         ||( jj == 6 ))
                  {
# endif
/*............................................................................*/
                     kk = vtx[jj][2];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                     kk = vtx[jj][3];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
               
                     kk = vtx[jj][0];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                     kk = vtx[jj][2];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

                     kk = vtx[jj][1];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
               
                     kk = vtx[jj][3];
                     fprintf( mwlfle3d, "%+.15e %+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));
                  };
# endif /* SPL_PLOTBOTTM == 1 */
/*............................................................................*/
                  fprintf( mwlfle3d, "\n" );
               };
# endif /* SPL_PLOTMWLLS != 0 */
/*............................................................................*/

               break;
            }; /* end switch ... */
            jj++ ;
         }; /* end while ... */
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         if( emk == ONE )
            fprintf( elwfle3d, "\n" );
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         if( mmk == ONE )
            fprintf( mwlfle3d, "\n" );
# endif
/*............................................................................*/
         if ( opt != FOUR )
         {
            jj = null;
            while( jj < FOUR ) /* plot vertical cell edges */
            {
               fprintf( cllfle3d, "%+.15e %+.15e %+.15e\n", 
                  ( ppt->cpt->c[tpt->cm[ii][jj]][0] ),
                  ( ppt->cpt->c[tpt->cm[ii][jj]][1] ),
                  ( ppt->cpt->c[tpt->cm[ii][jj]][2] ));
              
               kk = jj+FOUR;
               fprintf( cllfle3d, "%+.15e %+.15e %+.15e\n\n",
                  ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                  ( ppt->cpt->c[tpt->cm[ii][kk]][1] ),
                  ( ppt->cpt->c[tpt->cm[ii][kk]][2] ));

               jj++;
            }; /* end while ... */
         }; /* end if ( opt != FOUR ) */
/*............................................................................*/
# if SPL_DISPLAY != 0
         if ( opt != FOUR )
         {
            ( dsp->state ) = ( ii  - mm );
            dsplay( dsp );
         };
# endif
/*............................................................................*/
      }; /* next cell index ii */

      fclose( cllfle3d );
/*............................................................................*/
# if SPL_DISPLAY != 0
      if ( opt != FOUR )
      {
/* clear display: */

         printf( CLEAR_LINE );
      };
# endif
/*............................................................................*/
      nseconds = time( timer );
      timeptr = ctime( &nseconds );

/*............................................................................*/
# if SPL_DISPLAY == 2
      printf( "\n 3D-plot data - file %s", cptr );
      printf( timefrm, timeptr );
# endif
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
      fclose( elwfle3d );
# if SPL_DISPLAY == 2
      printf( "\n 3D-plot data - file %s", eptr );
      printf( timefrm, timeptr );
# endif
# endif /* SPL_PLOTEWLLS != 0 */
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
      fclose( mwlfle3d );
# if SPL_DISPLAY == 2
      printf( "\n 3D-plot data - file %s", mptr );
      printf( timefrm, timeptr );
# endif
# endif /* SPL_PLOTMWLLS != 0 */
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
      jj = null;
      while( jj < ( ssp->media ))
      {
         fclose( medfle3d[jj] );
# if SPL_DISPLAY == 2
         printf( "\n 3D-plot data - file %s", mdpt[jj] );
         printf( timefrm, timeptr );
# endif
         jj++ ;
      };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
   };
/*............................................................................*/
   if (( opt == ONE )
     ||( opt == THREE )
     ||( opt == FOUR ))
   {
/* 2D-scale: */

      if ( opt == FOUR )
      {
         mm = ( ssp->mi );
         nn = mm + ( ssp->pm ) - ONE;
      }
      else
      {
         mm = ( ssp->mi );
         nn = ( ssp->mf );
      };

         mm = ( ssp->mi );
         nn = mm + ( ssp->pm ) - ONE;

      x_min = 1.e277;
      y_min = 1.e277;

      x_max = -x_min;
      y_max = -y_min;

      for ( ii=mm; ii<=nn; ii++ )
      {
         for ( jj=null; jj<FOUR; jj++ )
         {
            kk = lbl[jj];

            xx = ( ppt->cpt->c[tpt->cm[ii][kk]][0] );
            yy = ( ppt->cpt->c[tpt->cm[ii][kk]][1] );

            if( x_max < xx )
               x_max = xx;

            if( y_max < yy ) 
               y_max = yy;

            if( xx < x_min ) 
               x_min = xx;

            if( yy < y_min ) 
               y_min = yy;
         };
      };

      xx = ( x_max - x_min )/2.;
      yy = ( y_max - y_min )/2.;

# if SPL_RESCLE_XY_ == 1
      if( xx < yy ) xx = yy;
      yy = xx;
# endif

      x_mean = ( x_min + x_max )/2.;
      x_min = x_mean + ( x_shift - 1./rscale )*xx;
      x_max = x_mean + ( x_shift + 1./rscale )*xx;

      y_mean = ( y_min + y_max )/2.;
      y_min = y_mean + ( y_shift - 1./rscale )*yy;
      y_max = y_mean + ( y_shift + 1./rscale )*yy;
/*............................................................................*/
/* 2D-plot system file: */

      strcpy( pptr, prefix );
      strcat( pptr, ( ssp->p2 ));

      if ( opt == FOUR )
         strcat( pptr, index );
      else
         strncat( pptr, ( ssp->flbl ), sizeof( &( ssp->flbl )));

      pltfle2d = fopen( pptr, "w" );

      if ( pltfle2d == null )
      {
         printf( "\n\n Message from function %s :", DO_SYSPLT );
         printf( "\n Unknown error on opening file '%s' !!!", pptr );
         printf( "\n [ overrides and returns to calling program ].\n" );
         ( spp->rtn ) = ONE;
         return spp;
      };
/*
      printf( "\n opened: 2D-plot system file '%s'", pptr );
*/
/*............................................................................*/

      strcpy( cptr, prefix );

/*............................................................................*/
# if SPL_PLOTEWLLS != 0
      strcpy( eptr, prefix );
# endif
# if SPL_PLOTMWLLS != 0
      strcpy( mptr, prefix );
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
      jj = null;
      while( jj < ( ssp->media ))
      {
         strcpy( mdpt[jj], prefix );
         jj++ ;
      };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      if ( opt == FOUR )
      {
         strcat( cptr, ( ssp->c3 ));
         strcat( cptr, index );
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         strcat( eptr, ( ssp->ew3 ));
         strcat( eptr, index );
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         strcat( mptr, ( ssp->mw3 ));
         strcat( mptr, index );
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
         jj = null;
         while( jj < ( ssp->media ))
         {
            strcat( mdpt[jj], ( ssp->md3[jj] ));
            strcat( mdpt[jj], index );
            jj++;
         };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      }
      else
      {
         strcat( cptr, ( ssp->c2 ));
         strncat( cptr, ( ssp->flbl ), sizeof( &( ssp->flbl )));
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         strcat( eptr, ( ssp->ew2 ));
         strncat( eptr, ( ssp->flbl ), sizeof( &(ssp->flbl )));
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         strcat( mptr, ( ssp->mw2 ));
         strncat( mptr, ( ssp->flbl ), sizeof( &(ssp->flbl )));
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
         jj = null;
         while( jj < ( ssp->media ))
         {
            strcat( mdpt[jj], ( spp->md2[jj] ));
            strncat( mdpt[jj], ( ssp->flbl ), sizeof( &( ssp->flbl )));
            jj++;
         };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      };

      fprintf( pltfle2d, "set title '2D-plot of %s %s",
         ( ppt->name ), ( ppt->text ));

      if ( opt == FOUR )
         fprintf( pltfle2d, " [ layer %s ]'\n", index );
      else 
         fprintf( pltfle2d, "'\n" );

      fprintf( pltfle2d, "set xrange [%.15e:%.15e]\n", x_min, x_max );
      fprintf( pltfle2d, "set yrange [%.15e:%.15e]\n", y_min, y_max );
      fprintf( pltfle2d, "set xlabel 'x/m'\n" );
      fprintf( pltfle2d, "set ylabel 'y/m'\n" );
      fprintf( pltfle2d, "set size square\n" );
      fprintf( pltfle2d, "set nogrid\n" );
      fprintf( pltfle2d, "set border\n" );
      fprintf( pltfle2d, "plot %c\n", 92 );

      fprintf( pltfle2d, "'%s%s' with lines", SPL_PATH, cptr );
/*............................................................................*/
# if SPL_PLOTEWLLS == 1
      fprintf( pltfle2d, ",%c\n'%s%s' with lines", 92, SPL_PATH, eptr );
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS == 1
      fprintf( pltfle2d, ",%c\n'%s%s' with lines", 92, SPL_PATH, mptr );
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
      jj = null;
      while( jj < ( ssp->media ))
      {
         fprintf( pltfle2d, ",%c\n'%s%s' with lines", 92, SPL_PATH, mdpt[jj] );
         jj++;
      };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
# if SPL_PLOTMWLLS == 2
      fprintf( pltfle2d, ",%c\n'%s%s' with lines", 92, SPL_PATH, mptr );
# endif
/*............................................................................*/
# if SPL_PLOTEWLLS == 2
      fprintf( pltfle2d, ",%c\n'%s%s' with lines", 92, SPL_PATH, eptr );
# endif
/*............................................................................*/
      fprintf( pltfle2d, "\npause -1 '[ hit return to continue ]'\n" );

      nseconds = time( timer );
      timeptr = ctime( &nseconds );

      fprintf( pltfle2d, "\n# 2D-plot system file '%s'", pptr );
      fprintf( pltfle2d, "\n# created:%.24s", timeptr );
      fprintf( pltfle2d, "\n# EOF" );

      fclose( pltfle2d );

# if SPL_DISPLAY != 0
      if( opt == THREE )
         printf( "\r 2D-plot system file %s", pptr );
      else
         printf( "\n 2D-plot system file %s", pptr );

      printf( timefrm, timeptr );
# endif
/*............................................................................*/
/* 2D-plot points file: */

      if ( opt != FOUR )
      {
         cllfle2d = fopen( cptr, "w" );

         if ( cllfle2d == null )
         {
            printf( "\n\n Message from function %s :", DO_SYSPLT );
            printf( "\n Unknown error on opening file '%s' !!!", cptr );
            printf( "\n [ overrides and returns to calling program ].\n" );
            ( spp->rtn ) = ONE;
            return spp;
         };
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         elwfle2d = fopen( eptr, "w" );

         if ( elwfle2d == null )
         {
            printf( "\n\n Message from function %s :", DO_SYSPLT );
            printf( "\n Unknown error on opening file '%s' !!!", eptr );
            printf( "\n [ overrides and returns to calling program ].\n" );
            ( spp->rtn ) = ONE;
            return spp;
         };
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         mwlfle2d = fopen( mptr, "w" );

         if ( mwlfle2d == null )
         {
            printf( "\n\n Message from function %s :", DO_SYSPLT );
            printf( "\n Unknown error on opening file '%s' !!!", mptr );
            printf( "\n [ overrides and returns to calling program ].\n" );
            ( spp->rtn ) = ONE;
            return spp;
         };
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
         jj = null;
         while( jj < ( ssp->media ))
         {
            medfle2d[jj] = fopen( mdpt[jj], "w" );

            if ( medfle2d[jj] == null )
            {
               printf( "\n\n Message from function %s :", DO_SYSPLT );
               printf( "\n Unknown error on opening file '%s' !!!", mdpt[jj] );
               printf( "\n [ overrides and returns to calling program ].\n" );
               ( spp->rtn ) = ONE;
               return spp;
            };
            jj++;
         }; 
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
         for ( ii=mm; ii<=nn; ii++ )
         {
            fprintf( cllfle2d, "# cell %ld\n", ii );
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
            emk = null;
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
            mmk = null;
# endif
/*............................................................................*/
            jj = null;
            while ( jj < FIVE )
            {
               kk = lbl[jj];
               fprintf( cllfle2d, "%+.15e %+.15e\n",
                  ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                  ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
               jj++;
            };
            fprintf( cllfle2d, "\n" );
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
            jj = null;
            while( jj < ( ssp->media ))
            {
               if (( ppt->mpt->idx[ii] ) == ( ssp->mdidx[jj] ))
               {
                  fprintf( medfle2d[jj], "# cell %ld\n", ii );
                  pp = null;
                  while ( pp < FIVE )
                  {
                     kk = lbl[(int)pp];
                     fprintf( medfle2d[jj], "%+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
                     pp++;
                  };

                  kk = THREE;
                  fprintf( medfle2d[jj], "%+.15e %+.15e\n\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
                  kk = ONE;
                  fprintf( medfle2d[jj], "%+.15e %+.15e\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
                  kk = TWO;
                  fprintf( medfle2d[jj], "%+.15e %+.15e\n\n",
                     ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                     ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
               };
               jj++ ;
            }; /* end while ( jj < ( ssp->media )) */
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
# if DSC_ADJCLL == 0
            jj = FACES;
# elif DSC_ADJCLL == 1
            jj = PORTS;
# endif
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
            jj = null;
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
            jj = null;
# endif
/*............................................................................*/
# if DSC_ADJCLL == 0
            while ( jj < FACES )
            {
               switch( jj )
               {
                 case 5:
                  break;
/*............................................................................*/
# if SPL_PLOTBOTTM == 0
                 case 4:
                  break;
# endif
/*............................................................................*/
# elif DSC_ADJCLL == 1
            while( jj < PORTS )
            {
               switch( jj )
               {
                 case 1:
                 case 4:
                  break;

/*............................................................................*/
# if SPL_PLOTBOTTM == 0
                 case 3:
                 case 6:
                  break;
# endif
/*............................................................................*/
# endif /* end if DSC_ADJCLL == 1 */
/*............................................................................*/

                 default: 
                  
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
                  if (( tpt->mn[ii][jj] ) == ELECTRIC_WALL )
                  {
                     if ( emk == null )
                     {
                        fprintf( elwfle2d, "# cell %ld\n", ii );
                        emk = ONE;
                     };

                     kk = vtx[jj][0];
                     fprintf( elwfle2d, "%+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                     kk = vtx[jj][1];
                     fprintf( elwfle2d, "%+.15e %+.15e\n\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
/*............................................................................*/
# if SPL_PLOTBOTTM == 1
/*............................................................................*/
# if DSC_ADJCLL == 0 
                     if ( jj == FOUR )
                     {
# elif DSC_ADJCLL == 1 
                     if (( jj == 3 )
                       ||( jj == 6 ))
                     {
# endif
/*............................................................................*/
                        kk = vtx[jj][2];
                        fprintf( elwfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][3];
                        fprintf( elwfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][0];
                        fprintf( elwfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][2];
                        fprintf( elwfle2d, "%+.15e %+.15e\n\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][1];
                        fprintf( elwfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][3];
                        fprintf( elwfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
                     };
# endif /* SPL_PLOTBOTTM == 1 */
/*............................................................................*/
                     fprintf( elwfle2d, "\n" );
                  };
# endif /* SPL_PLOTEWLLS != 0 */
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
                  if (( tpt->mn[ii][jj] ) == MAGNETIC_WALL )
                  {
                     if ( mmk == null )
                     {
                        fprintf( mwlfle2d, "# cell %ld\n", ii );
                        mmk = ONE;
                     };

                     kk = vtx[jj][0];
                     fprintf( mwlfle2d, "%+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                     kk = vtx[jj][1];
                     fprintf( mwlfle2d, "%+.15e %+.15e\n",
                        ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                        ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
/*............................................................................*/
# if SPL_PLOTBOTTM == 1
/*............................................................................*/
# if DSC_ADJCLL == 0
                     if ( jj == FOUR )
                     {
# elif DSC_ADJCLL == 1
                     if (( jj == 3 )
                       ||( jj == 6 ))
                     {
# endif
/*............................................................................*/
                        kk = vtx[jj][2];
                        fprintf( mwlfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][3];
                        fprintf( mwlfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][0];
                        fprintf( mwlfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
			
                        kk = vtx[jj][2];
                        fprintf( mwlfle2d, "%+.15e %+.15e\n\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][1];
                        fprintf( mwlfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));

                        kk = vtx[jj][3];
                        fprintf( mwlfle2d, "%+.15e %+.15e\n",
                           ( ppt->cpt->c[tpt->cm[ii][kk]][0] ),
                           ( ppt->cpt->c[tpt->cm[ii][kk]][1] ));
                     };

# endif /* SPL_PLOTBOTTM == 1 */
/*............................................................................*/
                     fprintf( mwlfle2d, "\n" );
                  };
# endif /* SPL_PLOTMWLLS != 0 */
/*............................................................................*/

                  break;
               }; /* end switch ... */
               jj++ ;
            }; /* end while ... */
/*............................................................................*/
         }; /* next cell index ii */

         fclose( cllfle2d );

         nseconds = time( timer );
         timeptr = ctime( &nseconds );

/*............................................................................*/
# if SPL_DISPLAY == 2
         printf( "\n 2D-plot data - file %s", cptr );
         printf( timefrm, timeptr );
# endif
/*............................................................................*/
# if SPL_PLOTEWLLS != 0
         fclose( elwfle2d );
   # if SPL_DISPLAY == 2
         printf( "\n 2D-plot data - file %s", eptr );
         printf( timefrm, timeptr );
   # endif
# endif
/*............................................................................*/
# if SPL_PLOTMWLLS != 0
         fclose( mwlfle2d );
   # if SPL_DISPLAY == 2
         printf( "\n 2D-plot data - file %s", mptr );
         printf( timefrm, timeptr );
   # endif
# endif
/*............................................................................*/
# if SPL_PLOTMEDIA != 0
         jj = null;
         while( jj < ( ssp->media ))
         {
            fclose( medfle2d[jj] );

# if SPL_DISPLAY == 2
            printf( "\n 2D-plot data - file %s", mdpt[jj] );
            printf( timefrm, timeptr );
# endif
            jj++;
         };
# endif /* SPL_PLOTMEDIA != 0 */
/*............................................................................*/
      }; /* end if opt != FOUR */
   };

   spp->rtn = null;
   return spp;
}
/*============================================================================*/
# undef SPL_INCLUDE
# undef SPL_PATH
# undef SPL_PREFIX
# undef SPL_DISPLAY
# undef SPL_LAYERS
# undef SPL_PLOTBOTTM
# undef SPL_PLOTEWLLS
# undef SPL_PLOTMWLLS
# undef SPL_PLOTMEDIA
# undef SPL_RESCLE_XY_
# undef SPL_RESCLE_XYZ
/***************** end of plot file configuration sysplt(*) *******************/
