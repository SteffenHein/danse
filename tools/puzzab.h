/* [ file: puzzab.h ] */
# define DO_PUZZLE "puzzab.h"
/*******************************************************************************
*                                                                              *
*   ANSI C function set puzzle.h                                               *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   A collection of C functions that permit the design of triangular           *
*   and quadrangular plane DSC submeshes                                       *
*                                                                              *
*   (C) SHEIN; Bad Aibling, February 2007                     Steffen Hein     *
*   [ Update: March 22, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# ifndef PUZZ_POINTS
  # define PUZZ_POINTS 501
# endif
/*----------------------------------------------------------------------------*/
# ifndef DSC_ADJCLL
   # define DSC_ADJCLL 0 
# endif
/*----------------------------------------------------------------------------*/
# include "../math/arcsct.h"




# define DO_QUDRL "qudrl(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function qudrl(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   DSC mesh generation subroutine for quadrangular domains [ type'AB']        *
*                                                                              *
*   In option[0] = 't' this function associates vertex points, indexed         *
*   top.cm[n][0,..,7], with mesh cells, labeled n=lbl.m+i (i=1,..,p*q).        *
*   ( label lbl.m - of type extern static - has to be transferred, if          *
*   necessary, before function execution ).                                    *
*   Boundary point indices for boundary AB are copied from integer             *
*   array qdl.cab[i] (i=1,...,p), if the number is greater than                *
*   null. Otherwise ( if qdl.cab[i] <= null ) a new vertex point               *
*   (index) is defined. After function execution all points                    *
*   of boundary AB are labeled in qdl.sab[i] ( i=0,...,p ).                    *
*   The same happens with boundary AD and the pertinent arrays                 *
*   qdl.cad[],qdl.sad[] and with BC, DC and qdl.cbc[], qdl.sbc[],              *
*   ..., qdl.sdc[] - in the latter two cases with p replaced by q.             *
*   Also, if qdl.ab[i] ( i=1,...,p ) is set to 'e' or 'm', respec-             *
*   tively, an electric or magnetic wall is defined at the i-th                *
*   boundary cell in the usual way ( viz. by associating -1, or -2             *
*   with the pertinent neighbouring port indicator top.mn ).                   *
*   - Analogous proceedings with AD, BC and DC.                                *
*   ( After function execution all wall pointers are reset to null. )          *
*                                                                              *
*   In *option = 'c' vertex point coordinates are written into                 *
*   array cor.c[n][].                                                          *
*   To this end, coordinates must be formerly defined as points                *
*   A = (qdl.ax,.ay,.az ),..., D = (qdl.dx,.dy,.dz), and arcs                  *
*   as qdl.aab, .aad, .abc, .aad in radians, i.e. in ]-2PI,2P[ .               *
*   - The sign of the arc indicates the sense of curvature .                   *
*   ( After function execution all arcs are reset to ZERO. )                   *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
called functions: argz(*), subarc(*), arcsct(*)
in header "../math/arcsct.h"
*/
/*----------------------------------------------------------------------------*/
# define QDL_PRINT 0       /* QDL_PRINT  = n: print vertex point coordinates  */
                           /*                 of mesh cell no. n              */
/*----------------------------------------------------------------------------*/
# ifndef QDL_BNDSPC
  # define QDL_BNDSPC 0    /* QDL_BNDSPC = 1: export boundary point spacing   */
# endif                    /*                 continuously to block interior  */
/*----------------------------------------------------------------------------*/
# ifndef ELECTRIC_WALL
   # define ELECTRIC_WALL ( -1 )
# endif
/*----------------------------------------------------------------------------*/
# ifndef MAGNETIC_WALL
   # define MAGNETIC_WALL ( -2 )
# endif
/*----------------------------------------------------------------------------*/
# ifndef EPSILON
   # define EPSILON ( 1.000e-277 )
# endif
/*----------------------------------------------------------------------------*/
# ifndef PRECISION
   # define PRECISION ( 1.000e-15 )
# endif
/*----------------------------------------------------------------------------*/
# define AB_POINTS PUZZ_POINTS
# define BC_POINTS PUZZ_POINTS
/*----------------------------------------------------------------------------*/
struct labels   
{
   long 
      m, /* the last mesh cell index [updated on each qudrl(*), trngl(*) call]*/
      p, /* the last vertex point index [updated ...""] */
      ci, /* the first mesh cell index of a generated block [...""] */
      cf, /* the last mesh cell index of a generated block [...""] */
      mi, /* the first mesh cell index in a blocks layer */  
      mf, /* the last mesh cell index in a blocks layer */
      pi, /* the first vertex point index in a blocks layer */
      pf; /* the last vertex point index in a blocks layer */
}; 
static struct labels lbl = {null};
/*----------------------------------------------------------------------------*/
struct quadrl   
{
   double 
    ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,
    aab,abc,aad,adc,saab,sabc,saad,sadc;

   long 
    cab[AB_POINTS],
    cbc[BC_POINTS],
    cad[BC_POINTS],
    cdc[AB_POINTS],
    sab[AB_POINTS],
    sbc[BC_POINTS],
    sad[BC_POINTS],
    sdc[AB_POINTS];

   char
     ab[AB_POINTS],
     bc[BC_POINTS],
     ad[BC_POINTS],
     dc[AB_POINTS],
    wab[AB_POINTS],
    wbc[BC_POINTS], 
    wad[BC_POINTS],
    wdc[AB_POINTS]; 

   char
      top, bot, trv;                 /* top, bot: boundary wall indicators    */
                                     /* trv: trivial node indicator           */
}; 
static struct quadrl qdl  = {ZERO};
/*----------------------------------------------------------------------------*/
# define SET_CORNER(i) \
{ \
   lbl.p++; \
 \
   if ( CPNTS < lbl.p ) \
   { \
      printf( "\n\n Error message from function %s :", __func__ ); \
      printf( "\n Corner point limit exceeded !!!" ); \
      printf( "\n [ Maximum number is %ld = macro CPNTS " \
         "in %s.", (long) CPNTS, "FORMER.CONF" ); \
      printf( "\n - Change macro only in compliance " \
         "with memory resources.]\n" ); \
      exit( EXIT_FAILURE ); \
   }; \
 \
   if ( lbl.pi <= null ) \
      lbl.pi = lbl.p; \
 \
   (i) = lbl.p; \
}
/*----------------------------------------------------------------------------*/
# define SET_NODE( ) \
{ \
   lbl.m++; \
 \
   if ( NODES < lbl.m ) \
   { \
      printf( "\n\n Error message from function %s :", __func__ ); \
      printf( "\n Mesh cell limit exceeded !!!" ); \
      printf( "\n [ Maximum number is %ld = macro NODES " \
         "in %s.", (long) NODES, "FORMER.CONF" ); \
      printf( "\n - Change macro only in compliance with memory " \
         "resources\n   and same macro in function scattr(*) " \
         "of program SOLVER.C.]\n" ); \
      exit( EXIT_FAILURE ); \
   }; \
 \
   if ( lbl.mi <= null ) \
      lbl.mi = lbl.m; \
}
/*============================================================================*/

short qudrl( short p, short q, char *option )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern struct labels lbl;
   extern struct quadrl qdl;
   extern struct subarc arc;
*/
/* declarartions: */

   static double
      xx = ZERO,
      yy = ZERO,
      rr = ZERO,
      ss = ZERO,

      px0 = ZERO,
      py0 = ZERO,
      pz0 = ZERO,
      px1 = ZERO,
      py1 = ZERO,
      pz1 = ZERO,

      qx0 = ZERO,
      qy0 = ZERO,
      qz0 = ZERO,
      qx1 = ZERO,
      qy1 = ZERO,
      qz1 = ZERO,

      dad = ZERO,
      dbc = ZERO,
      dab = ZERO,
      ddc = ZERO,

      alfa0 = ZERO,
      alfa1 = ZERO;

   static const double
      bound1 = 2.*PI*PRECISION;

   static long rp[AB_POINTS] = {null};

   static long 
      j = null,
      k = null,
      l = null,
      m = null;

   static short 
      h = null,
      i = null,
      ind = null;

   static char 
      o = null;

# if QDL_BNDSPC == 1
   static char 
      arcab = null,
      arcad = null,
      arcbc = null,
      arcdc = null;
# endif
    
# if (( QDL_PRINT != null )\
    ||( TGL_PRINT != null ))
   static char *ptr;
# endif

# if QDL_BNDSPC == 1

   static double
      xx0 = ZERO,
      xx1 = ZERO,
      xx2 = ZERO,
      xx3 = ZERO,
      xx4 = ZERO,
      xx5 = ZERO,
      xx6 = ZERO,
      xx7 = ZERO,

      yy0 = ZERO,
      yy1 = ZERO,
      yy2 = ZERO,
      yy3 = ZERO,
      yy4 = ZERO,
      yy5 = ZERO,
      yy6 = ZERO,
      yy7 = ZERO,

      zz0 = ZERO,
      zz1 = ZERO,
      zz2 = ZERO,
      zz3 = ZERO,

# if QDL_PRINT != null
      zz5 = ZERO,
      zz7 = ZERO,
# endif

      rr0 = ZERO,
      rr1 = ZERO,
      rr2 = ZERO,
      rr3 = ZERO,
      rr4 = ZERO,
      rr5 = ZERO,
      rr6 = ZERO,
      rr7 = ZERO,

      rx0 = ZERO,
      ry0 = ZERO,

      rx1 = ZERO,
      ry1 = ZERO,

      sx0 = ZERO,
      sy0 = ZERO,

      sx1 = ZERO,
      sy1 = ZERO,

      dd0 = ZERO,
      dd1 = ZERO,
      dd2 = ZERO,
      dd3 = ZERO,

      da0 = ZERO,
      da1 = ZERO,
      da2 = ZERO,
      da3 = ZERO;

   static double
      beta0 = ZERO,
      beta1 = ZERO;

   static short
      hh0 = null,
      hh1 = null,
      hl0 = null,
      hl1 = null,
      dh0 = null,
      dh1 = null,
      ii0 = null,
      ii1 = null,
      il0 = null,
      il1 = null,
      di0 = null,
      di1 = null;

# endif  /* end if QDL_BNDSPC == 1 */
/*............................................................................*/
/* prototypes: */

   double argz( double x, double y, short n );
   double subarc( double x, double y, double z, double w,
                                      double alfa, double r );
# if QDL_BNDSPC == 1
   double arcsct( double ax, double ay, double bx, double by, double alfa,
                  double cx, double cy, double dx, double dy, double beta );
# endif

# if (( QDL_PRINT != null )\
    ||( TGL_PRINT != null ))
   ptr = ( char *) calloc( SHS_SIZE , ONE );
# endif
/*----------------------------------------------------------------------------*/

   if ( AB_POINTS - ONE < p )
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Too many horizontal [u-] divisions "
         "in block generation routine !!!" );
      printf( "\n [ Maximum number is %ld = macro AB_POINTS - 1 "
         "in %s.", (long) AB_POINTS-ONE, __func__ );
      printf( "\n - Change macro in compliance with memory resources.]\n" );
      exit( EXIT_FAILURE );
   };

   if ( BC_POINTS - ONE < q )
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Too many vertical [v-] divisions "
         "in block generation routine !!!" );
      printf( "\n [ Maximum number is %ld = macro BC_POINTS - 1 "
         "in %s.", (long) BC_POINTS-ONE, __func__ );
      printf( "\n - Change macro in compliance with memory resources.]\n" );
      exit( EXIT_FAILURE );
   }; 

/*............................................................................*/
/* boundary point transfer for null divisions in A-B and D-C directions: */

   if ( p == null )
   {
      for ( h=null; h<=q; h++ )
      {
         if (( null < qdl.cad[h] )&&( null < qdl.cbc[h] )) 
         {
            if ( qdl.cad[h] != qdl.cbc[h] )
            {  
               printf( "\n\n Error message from function %s :", __func__ );
               printf( "\n Null divisions in A-B direction conflicting" );
               printf( "\n with boundary point transfer instructions on" );
               printf( "\n boundaries AD and BC,\n points %ld "
                  "and %ld !!!", qdl.cad[h], qdl.cbc[h] );
               printf( "\n\n [ Program execution stopped.]\n" );
               exit( EXIT_FAILURE );
            };
         };

         qdl.sbc[h] = qdl.cad[h];
         qdl.sad[h] = qdl.cbc[h]; 

         if ( qdl.ad[h] != null ) qdl.wbc[h]=qdl.ad[h];
         if ( qdl.bc[h] != null ) qdl.wad[h]=qdl.bc[h];

         qdl.ad[h] = null;
         qdl.bc[h] = null;

         if (( null < h )&&( h < q ))
         {
            qdl.cad[h] = null;
            qdl.cbc[h] = null;
         };
      };

      if ( null < qdl.cad[0] ) 
         qdl.sab[0] = qdl.cad[0];
      if ( null < qdl.cbc[0] ) 
         qdl.sab[0] = qdl.cbc[0];
      if ( null < qdl.cad[q] ) 
         qdl.sdc[0] = qdl.cad[q];
      if ( null < qdl.cbc[q] ) 
         qdl.sdc[0] = qdl.cbc[q];

      qdl.sabc = qdl.aad;
      qdl.saad = qdl.abc;

   };/* end if p == null */

/*............................................................................*/
/* boundary point transfer for null divisions in A-D and B-C directions: */

   if ( q == null )
   {
      for ( i=null; i<=p; i++ )
      {
         if (( null < qdl.cab[i] )&&( null < qdl.cdc[i] ))
         {
            if ( qdl.cab[i] != qdl.cdc[i] )
            {
               printf( "\n\n Error message from function %s :", __func__ );
               printf( "\n Null divisions in A-D direction conflicting" );
               printf( "\n with boundary point transfer instructions on" );
               printf( "\n boundaries AB and DC,\n points %ld "
                  "and %ld !!!", qdl.cab[i], qdl.cdc[i] );
               printf( "\n\n [ Program execution stopped.]\n" );
               exit( EXIT_FAILURE );
            };
         };

         qdl.sdc[i] = qdl.cab[i];
         qdl.sab[i] = qdl.cdc[i];
           
         if ( qdl.ab[i] != null ) 
            qdl.wdc[i]=qdl.ab[i];
         if ( qdl.dc[i] != null ) 
            qdl.wab[i]=qdl.dc[i];

         qdl.ab[i] = null;
         qdl.dc[i] = null;

         if (( null < i )&&( i < p )) 
         {
            qdl.cab[i] = null;
            qdl.cdc[i] = null;
         };   
      };

      if ( null < qdl.cab[0] ) 
         qdl.sad[0] = qdl.cab[0];
      if ( null < qdl.cdc[0] ) 
         qdl.sad[0] = qdl.cdc[0];
      if ( null < qdl.cab[p] ) 
         qdl.sbc[0] = qdl.cab[p];
      if ( null < qdl.cdc[p] ) 
         qdl.sbc[0] = qdl.cdc[p];

      qdl.sadc = qdl.aab;
      qdl.saab = qdl.adc;

   };/* end if q == null */

   if (( q == null )||( p == null ))
   {
      qdl.cad[0] = null;
      qdl.cbc[0] = null;
      qdl.cad[q] = null;
      qdl.cbc[q] = null;

      qdl.cab[0] = null;
      qdl.cdc[0] = null;
      qdl.cab[p] = null;
      qdl.cdc[p] = null;

      qdl.aad = ZERO;
      qdl.abc = ZERO;

      qdl.aab = ZERO;
      qdl.adc = ZERO;

      return null;
   }
   else
      lbl.ci = lbl.m + ONE;
/*............................................................................*/

   if ( *option == 'c' )         /* copy boundary point coordinates if        */
   {                             /* boundary points are defined [ indexed>0]  */
      if ( null < qdl.cab[null] )
      {
         k = qdl.cab[null];
         qdl.ax = spt->ppt->cpt->c[k][0];
         qdl.ay = spt->ppt->cpt->c[k][1];
         qdl.az = spt->ppt->cpt->c[k][2];
      };
      if ( null < qdl.cab[p] )
      {
         k = qdl.cab[p];
         qdl.bx = spt->ppt->cpt->c[k][0];
         qdl.by = spt->ppt->cpt->c[k][1];
         qdl.bz = spt->ppt->cpt->c[k][2];
      };
      if ( null < qdl.cad[null] )
      {
         k = qdl.cad[null];
         qdl.ax = spt->ppt->cpt->c[k][0];
         qdl.ay = spt->ppt->cpt->c[k][1];
         qdl.az = spt->ppt->cpt->c[k][2];
      };
      if ( null < qdl.cad[q] )
      {
         k = qdl.cad[q];
         qdl.dx = spt->ppt->cpt->c[k][0];
         qdl.dy = spt->ppt->cpt->c[k][1];
         qdl.dz = spt->ppt->cpt->c[k][2];
      };
      if ( null < qdl.cbc[null] )
      {
         k = qdl.cbc[null];
         qdl.bx = spt->ppt->cpt->c[k][0];
         qdl.by = spt->ppt->cpt->c[k][1];
         qdl.bz = spt->ppt->cpt->c[k][2];
      };
      if ( null < qdl.cbc[q] )
      {
         k = qdl.cbc[q];
         qdl.cx = spt->ppt->cpt->c[k][0];
         qdl.cy = spt->ppt->cpt->c[k][1];
         qdl.cz = spt->ppt->cpt->c[k][2];
      };
      if ( null < qdl.cdc[null] )
      {
         k = qdl.cdc[null];
         qdl.dx = spt->ppt->cpt->c[k][0];
         qdl.dy = spt->ppt->cpt->c[k][1];
         qdl.dz = spt->ppt->cpt->c[k][2];
      }; 
      if ( null < qdl.cdc[p] )
      {
         k = qdl.cdc[p];
         qdl.cx = spt->ppt->cpt->c[k][0];
         qdl.cy = spt->ppt->cpt->c[k][1];
         qdl.cz = spt->ppt->cpt->c[k][2];
      };

      xx = qdl.ax - qdl.bx;
      yy = qdl.ay - qdl.by;
      dab = sqrt( xx*xx + yy*yy );

      if ( dab < EPSILON )
         dab = ZERO; 

      ss = fabs( qdl.aab );

      if (( p <= ONE )||( ss < bound1 )||( dab < EPSILON ))
      {
         qdl.aab = ZERO;

# if QDL_BNDSPC == 1 
         arcab = null;  
# endif
      }
# if QDL_BNDSPC == 1
      else if (( PI + EPSILON ) < ss )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Illegal arc transferred on side AB !" ); 
         printf( "\n [ Arc AB absolutely exceeds PI.]\n" );  
         exit( EXIT_FAILURE );
      }
      else
         arcab = ONE;
# endif

      xx = qdl.ax - qdl.dx;
      yy = qdl.ay - qdl.dy;
      dad = sqrt( xx*xx + yy*yy );

      if ( dad < EPSILON )
         dad = ZERO;

      ss = fabs( qdl.aad );

      if (( q <= ONE )||( ss < bound1 )||( dad < EPSILON ))
      {
         qdl.aad = ZERO;
# if QDL_BNDSPC == 1
         arcad = null;
# endif
      }
# if QDL_BNDSPC == 1
      else if ( ( PI + EPSILON ) < ss )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Illegal arc transferred on side AD !" );
         printf( "\n [ Arc AD absolutely exceeds PI.]\n" );
         exit( EXIT_FAILURE );
      }
      else
         arcad = ONE;
# endif

      xx = qdl.bx - qdl.cx;
      yy = qdl.by - qdl.cy;
      dbc = sqrt( xx*xx + yy*yy );

      if ( dbc < EPSILON )
         dbc = ZERO;

      ss = fabs( qdl.abc );

      if (( q <= ONE )||( ss < bound1 )||( dbc < EPSILON ))
      {
         qdl.abc = ZERO;

# if QDL_BNDSPC == 1
         arcbc = null;
# endif
      }
# if QDL_BNDSPC == 1
      else if ( ( PI + EPSILON ) < ss )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Illegal arc transferred on side BC !" );
         printf( "\n [ Arc BC absolutely exceeds PI.]\n" );
         exit( EXIT_FAILURE );
      }
      else
         arcbc = ONE;
# endif

      xx = qdl.dx - qdl.cx;
      yy = qdl.dy - qdl.cy;
      ddc = sqrt( xx*xx + yy*yy );

      if ( ddc < EPSILON )
         ddc = ZERO;

      ss = fabs( qdl.adc );

      if (( p <= ONE )||( ss < bound1 )||( ddc < EPSILON )) 
      {
         qdl.adc = ZERO;

# if QDL_BNDSPC == 1
         arcdc = null;
# endif
      }
# if QDL_BNDSPC == 1
      else if ( ( PI + EPSILON ) < ss )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Illegal arc transferred on side DC !" );
         printf( "\n [ Arc DC absolutely exceeds PI.]\n" );
         exit( EXIT_FAILURE );
      }
      else
         arcdc = ONE;

      hh0 = ONE;
      rr1 = ZERO;
      xx1 = qdl.ax;
      yy1 = qdl.ay;
      zz1 = qdl.az;

      hh1 = ONE;
      rr3 = ZERO;
      xx3 = qdl.bx;
      yy3 = qdl.by;
      zz3 = qdl.bz;

# endif 

   };
     
   for ( i=null ; i<=p ; i++ )
   {
      rp[i] = null;
   };
        
/*............................................................................*/
/* here starts mesh generation part */

   for ( h=ONE ; h<=q ; h++ ) 
   {
      if ( *option == 'c' )
      {
# if QDL_BNDSPC == 1
         rr0 = rr1;

         if ( h == hh0 )
         {
            xx0 = xx1;
            yy0 = yy1;
            zz0 = zz1;

            while(( qdl.cad[hh0] == null )&&( hh0 < q )) 
               hh0++;
   
            if ( hh0 < q )            
            {
               xx1 = spt->ppt->cpt->c[qdl.cad[hh0]][0];
               yy1 = spt->ppt->cpt->c[qdl.cad[hh0]][1];
               zz1 = spt->ppt->cpt->c[qdl.cad[hh0]][2];
            }
            else if ( q <= hh0 )
            {
               xx1 = qdl.dx;
               yy1 = qdl.dy;
               zz1 = qdl.dz;
            };
            hh0++;
            hl0 = h;
            dh0 = hh0 - h;

            if ( EPSILON < dad )
            {
               xx = xx1 - xx0;
               yy = yy1 - yy0;
               dd0 = sqrt( xx*xx + yy*yy );
               
               if ( arcad == ONE )
               {
                  ss = dd0*sin( qdl.aad/2. )/dad;

                  if ( ( 1. - EPSILON ) < fabs( ss ) )
                  {
                     da0 =  PI;
                  }
                  else if ( fabs( ss ) < -1. + EPSILON )
                  {
                     da0 = -PI;
                  }
                  else
                    da0 = 2.*asin( dd0*sin( qdl.aad/2. )/dad );
               }
               else
                  da0 = ZERO;
            }
            else /* dad ~ ZERO */
            {
               dd0 = ZERO;
               da0 = ZERO;
            };
         };/* end if h == hh0 */

         if ( arcad == ONE )
         {
            rr1 += ( da0/( qdl.aad*dh0 ) );
         }
         else if ( EPSILON < dad )
         {
            rr1 += ( dd0/( dad*dh0 ) );
         }
         else
         {
            rr1 += ( ( double ) 1./q );
         };

         rr = ( double )( h-hl0 )/dh0;
/*............................................................................*/
         ss = subarc( xx0, yy0, xx1, yy1, da0, rr );       /*                 */
/*.......................................................*/
         px0 = arc.rx;
         py0 = arc.ry;
         pz0 = zz0 + rr*( zz1 - zz0 );

         rr = ( double ) ( h-hl0+ONE )/dh0;
/*............................................................................*/
         ss = subarc( xx0, yy0, xx1, yy1, da0, rr );       /*                 */
/*.......................................................*/
         px1 = arc.rx;
         py1 = arc.ry;
         pz1 = zz0 + rr*( zz1 - zz0 );

/* points q0, q1: */ 

         rr2 = rr3;

         if ( h == hh1 )
         {
            xx2 = xx3;
            yy2 = yy3;
            zz2 = zz3;

            while(( qdl.cbc[hh1] == null )&&( hh1 < q ))
               hh1++ ;

            if ( hh1 < q )
            {
               xx3 = spt->ppt->cpt->c[qdl.cbc[hh1]][0];
               yy3 = spt->ppt->cpt->c[qdl.cbc[hh1]][1];
               zz3 = spt->ppt->cpt->c[qdl.cbc[hh1]][2];
            }
            else if ( q <= hh1 )
            {
               xx3 = qdl.cx;
               yy3 = qdl.cy;
               zz3 = qdl.cz;
            };
            hh1++;
            hl1 = h;
            dh1 = hh1 - h;

            if ( EPSILON < dbc )
            {
               xx = xx3 - xx2;
               yy = yy3 - yy2;
               dd1 = sqrt( xx*xx + yy*yy );

               if ( arcbc == ONE ) 
               {
                  ss = dd1*sin( qdl.abc/2. )/dbc;

                  if ( fabs( ss ) > 1.-EPSILON )
                  {
                     da1 =  PI;
                  }
                  else if ( fabs( ss ) < ( -1. + EPSILON ))
                  {
                     da1 = -PI;
                  }
                  else
                    da1 = 2.*asin( dd1*sin( qdl.abc/2. )/dbc );
               } 
               else
                  da1 = ZERO;
            }
            else
            {
               dd1 = ZERO;
               da1 = ZERO;
            };
         }; /* end if h == hh1 */

         if ( arcbc == ONE )
         {
            rr3 += ( da1/( qdl.abc*dh1 ));
         }
         else if ( EPSILON < dbc )
         {
            rr3 += ( dd1/( dbc*dh1 ));
         }
         else
         {
            rr3 += (( double ) 1./q );
         };

         rr = ( double )( h-hl1 )/dh1;
/*............................................................................*/
         ss = subarc( xx2, yy2, xx3, yy3, da1, rr );       /*                 */
/*.......................................................*/
         qx0 = arc.rx;
         qy0 = arc.ry;
         qz0 = zz2 + rr*( zz3 - zz2 );

         rr = ( double ) ( h-hl1+ONE )/dh1;
/*............................................................................*/
         ss = subarc( xx2, yy2, xx3, yy3, da1, rr );       /*                 */
/*.......................................................*/
         qx1 = arc.rx;
         qy1 = arc.ry;
         qz1 = zz2 + rr*( zz3 - zz2 );

         rr = .5*( rr0 + rr2 );
         alfa0 = qdl.aab + rr*( qdl.adc - qdl.aab );

         if ( fabs( alfa0 ) < bound1 )
            alfa0 = ZERO;

         rr = .5*( rr1 + rr3 );
         alfa1 = qdl.aab + rr*( qdl.adc - qdl.aab );

         if ( fabs( alfa1 ) < bound1 )
            alfa1 = ZERO;

         ii0 = ONE;
         rr5 = ZERO;
         xx5 = qdl.ax;
         yy5 = qdl.ay;

# if QDL_PRINT != null
         zz5 = qdl.az;
# endif
         ii1 = ONE;
         rr7 = ZERO;
         xx7 = qdl.dx;
         yy7 = qdl.dy;

# if QDL_PRINT != null
         zz7 = qdl.dz;
# endif

# else /* if QDL_BNDSPC == null */

         rr = ( double ) ( h - ONE )/q;
/*............................................................................*/
         ss = subarc( qdl.ax, qdl.ay, qdl.dx, qdl.dy, qdl.aad, rr );    /*    */
/*....................................................................*/
         px0 = arc.rx;
         py0 = arc.ry;
         pz0 = qdl.az + rr*( qdl.dz - qdl.az );

/*............................................................................*/
         ss = subarc( qdl.bx, qdl.by, qdl.cx, qdl.cy, qdl.abc, rr );    /*    */
/*....................................................................*/
         qx0 = arc.rx;
         qy0 = arc.ry;
         qz0 = qdl.bz + rr*( qdl.cz - qdl.bz );

         alfa0 = qdl.aab + rr*( qdl.adc - qdl.aab );

         rr = ( double ) h/q;
/*............................................................................*/
         ss = subarc( qdl.ax, qdl.ay, qdl.dx, qdl.dy, qdl.aad, rr );    /*    */
/*....................................................................*/
         px1 = arc.rx;
         py1 = arc.ry;
         pz1 = qdl.az + rr*( qdl.dz - qdl.az );

/*............................................................................*/
         ss = subarc( qdl.bx, qdl.by, qdl.cx, qdl.cy, qdl.abc, rr );    /*    */
/*....................................................................*/
         qx1 = arc.rx;
         qy1 = arc.ry;
         qz1 = qdl.bz + rr*( qdl.cz - qdl.bz );

         alfa1 = qdl.aab + rr*( qdl.adc - qdl.aab );

# endif /* end if QDL_BNDSPC != 1 */

/*............................................................................*/
/* print border points: */

# if QDL_PRINT != null

         if ( QDL_PRINT == lbl.m + ONE )
         {
            printf( "\n\n Intermediate values from function"
                    " '%s' ", __func__ );
            printf( "\n [ generating cell no. %ld ]: " , lbl.m + ONE );
            printf( "\n\n p0 = ( %.5e , %.5e , %.5e )", px0, py0, pz0 );
            printf( "\n q0 = ( %.5e , %.5e , %.5e )", qx0, qy0, qz0 );
            printf( "\n p1 = ( %.5e , %.5e , %.5e )", px1, py1, pz1 );
            printf( "\n q1 = ( %.5e , %.5e , %.5e )", qx1, qy1, qz1 );
         };
# endif
      };

      j = null;
      k = null;
      l = null;
      m = null;
        
      for ( i=ONE ; i<=p ; i++ )
      {

# if QDL_BNDSPC == 1
         if ( *option == 'c' )
         {
            rr4 = rr5;

            if ( i == ii0 )
            {
               xx4 = xx5;
               yy4 = yy5;

               while(( qdl.cab[ii0] == null )&&( ii0 < p ))
                  ii0++;

               if ( ii0 < p )
               {
                  xx5 = spt->ppt->cpt->c[qdl.cab[ii0]][0];
                  yy5 = spt->ppt->cpt->c[qdl.cab[ii0]][1];
# if QDL_PRINT != null
                  zz5 = spt->ppt->cpt->c[qdl.cab[ii0]][2];
# endif
               }
               else if ( p <= ii0 )
               {
                  xx5 = qdl.bx;
                  yy5 = qdl.by;
# if QDL_PRINT != null
                  zz5 = qdl.bz;
# endif
               };
               ii0++;
               il0 = i;
               di0 = ii0 - i;

               if ( EPSILON < dab )
               {
                  xx = xx5 - xx4;
                  yy = yy5 - yy4;
                  dd2 = sqrt( xx*xx + yy*yy );

                  if ( arcab == ONE )
                  {
                     ss = dd2*sin( qdl.aab/2. )/dab;

                     if ( ( 1. - EPSILON ) < fabs( ss ) )
                     {
                        da2 =  PI;
                     }
                     else if ( fabs( ss ) < -1. + EPSILON )
                     {
                        da2 = -PI;
                     }
                     else
                       da2 = 2.*asin( dd2*sin( qdl.aab/2. )/dab );
                  }
                  else
                     da2 = ZERO; 
               }
               else
               {
                  dd2 = ZERO;
                  da2 = ZERO;
               };
            };

            if ( arcab == ONE )
            {
               rr5 += ( da2/( qdl.aab*di0 ));
            }
            else if ( EPSILON < dab )
            {
               rr5 += ( dd2/( dab*di0 ));
            }
            else
            {
               rr5 += (( double ) 1./p );
            };

            rr = ( double ) ( i-il0 )/di0;
/*............................................................................*/
            ss = subarc( xx4, yy4, xx5, yy5, da2, rr );       /*              */
/*..........................................................*/
            rx0 = arc.rx;
            ry0 = arc.ry;

            rr = ( double ) ( i-il0+ONE )/di0;
/*............................................................................*/
            ss = subarc( xx4, yy4, xx5, yy5, da2, rr );       /*              */
/*..........................................................*/
            rx1 = arc.rx;
            ry1 = arc.ry;

            rr6 = rr7;

            if ( i == ii1 )
            {
               xx6 = xx7;
               yy6 = yy7;

               while(( qdl.cdc[ii1] == null )&&( ii1 < p ))
                  ii1++;

               if ( ii1 < p )
               {
                  xx7 = spt->ppt->cpt->c[qdl.cdc[ii1]][0];
                  yy7 = spt->ppt->cpt->c[qdl.cdc[ii1]][1];
# if QDL_PRINT != null
                  zz7 = spt->ppt->cpt->c[qdl.cdc[ii1]][2];
# endif
               }
               else if ( p <= ii1 )
               {
                  xx7 = qdl.cx;
                  yy7 = qdl.cy;
# if QDL_PRINT != null
                  zz7 = qdl.cz;
# endif
               };
               ii1++;
               il1 = i;
               di1 = ii1 - i;

               if ( EPSILON < ddc )
               {
                  xx = xx7 - xx6;
                  yy = yy7 - yy6;
                  dd3 = sqrt( xx*xx + yy*yy );

                  if ( arcdc == ONE )
                  {
                     ss = dd3*sin( qdl.adc/2. )/ddc;

                     if ( ( 1. - EPSILON ) < fabs( ss ) )
                     {
                        da3 =  PI;
                     }
                     else if ( fabs( ss ) < ( -1. + EPSILON ))
                     {
                        da3 = -PI;
                     }
                     else
                       da3 = 2.*asin( dd3*sin( qdl.adc/2. )/ddc );
                  }
                  else
                     da3 = ZERO;
               }
               else
               {
                  dd3 = ZERO;
                  da3 = ZERO;
               };
            };

            if ( arcdc == ONE )
            {
               rr7 += ( da3/( qdl.adc*di1 ));
            }
            else if ( EPSILON < ddc )
            {
               rr7 += ( dd3/( ddc*di1 ));
            }
            else
            {
               rr7 += (( double ) 1./p );
            };

            rr = ( double ) ( i-il1 )/di1;
/*............................................................................*/
            ss = subarc( xx6, yy6, xx7, yy7, da3, rr );       /*              */
/*..........................................................*/
            sx0 = arc.rx;
            sy0 = arc.ry;
            
            rr = ( double ) ( i-il1+ONE )/di1;
/*............................................................................*/
            ss = subarc( xx6, yy6, xx7, yy7, da3, rr );       /*              */
/*..........................................................*/
            sx1 = arc.rx;
            sy1 = arc.ry;

# if QDL_PRINT != null
            if ( QDL_PRINT == lbl.m + ONE )
            {
               printf( "\n\n r0 = ( %.5e, %.5e, %.5e )", rx0, ry0, zz5 );
               printf( "\n s0 = ( %.5e, %.5e, %.5e )", sx0, sy0, zz7 );
               printf( "\n r1 = ( %.5e, %.5e, %.5e )", rx1, ry1, zz5 );
               printf( "\n s1 = ( %.5e, %.5e, %.5e )", sx1, sy1, zz7 );
               printf( "\n\n [ please acknowledge ] " );
               scanf( "%s", ptr );
            };
# endif
            rr = .5*( rr4 + rr6 );
            beta0 = qdl.aad + rr*( qdl.abc - qdl.aad );

            if ( fabs( beta0 ) < bound1 )
               beta0 = ZERO;

            rr = .5*( rr5 + rr7 );
            beta1 = qdl.aad + rr*( qdl.abc - qdl.aad );

            if ( fabs( beta1 ) < bound1 )
               beta1 = ZERO;
         };

# endif /* QDL_BNDSPC == 1 */

         SET_NODE( );

     /* set_j: */

         j = k;
         if ( null < j ) 
            goto set_k;
         if ( h == ONE ) 
            j = qdl.cab[i-ONE];     /* transfer point */
         if ( null < j ) 
            goto set_k;             /* from adjacent domains, */
         if ( i == ONE ) 
            j = qdl.cad[h-ONE];     /* interfaces */
         if ( null < j ) 
            goto set_k;

         j = rp[i-ONE];

         if ( null < j ) 
            goto set_k;

         SET_CORNER(j);

         if ( *option == 'c' )
         {

# if QDL_BNDSPC == 1 
/*............................................................................*/
            ss = arcsct( px0, py0, qx0, qy0, alfa0, rx0, ry0, sx0, sy0, beta0 );
/*............................................................................*/
            if ( ss < -PRECISION )
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling function arcsct(*) !!!" );
               printf( "\n [ Unable to find intersection point of arcs "
                  "AB, CD ; ab_idx: %d, ad_idx: %d,", i, h );
               printf( "\n returned value:\t    ss = %+.7e", ss ); 
               printf( "\n A = ( %+.5e, %+.5e ) ;  B = ( %+.5e, %+.5e ),"
                  "\n\t\t\t alfa0 = %+.7e rad,", px0, py0, qx0, qy0, alfa0 );
               printf( "\n C = ( %+.5e, %+.5e ) ;  D = ( %+.5e, %+.5e ),"
                  "\n\t\t\t beta0 = %+.7e rad.", rx0, ry0, sx0, sy0, beta0 );
               printf( "\n The actual point index is %ld on cell %ld.]\n",
                  j, lbl.m );
               exit( EXIT_FAILURE );
            };                        
            spt->ppt->cpt->c[j][0] = arc.xx;
            spt->ppt->cpt->c[j][1] = arc.yy;
            spt->ppt->cpt->c[j][2] = pz0 + ss*( qz0 - pz0 );

# else /* if QDL_BNDSPC != 1 */
            rr = ( double ) ( i - ONE )/p;
/*............................................................................*/
            ss = subarc( px0, py0, qx0, qy0, alfa0, rr );       /*            */
/*............................................................*/
            spt->ppt->cpt->c[j][0] = arc.rx;
            spt->ppt->cpt->c[j][1] = arc.ry;
            spt->ppt->cpt->c[j][2] = pz0 + rr*( qz0 - pz0 );
# endif
         };

        set_k:

         k = null;

         if ( null < k ) 
            goto set_l;
         if ( i == p ) 
            k = qdl.cbc[h-ONE];
         if ( null < k ) 
            goto set_l;
         if ( h == ONE ) 
            k = qdl.cab[i];
         if ( null < k ) 
            goto set_l;

         k = rp[i];

         if ( null < k ) 
            goto set_l;

         SET_CORNER(k); 

         if ( *option == 'c' )
         {

# if QDL_BNDSPC == 1
/*............................................................................*/
            ss = arcsct( px0, py0, qx0, qy0, alfa0, rx1, ry1, sx1, sy1, beta1 );
/*............................................................................*/
            if ( ss < -PRECISION )
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling function arcsct(*) !!!" );
               printf( "\n [ Unable to find intersection point of arcs "
                  "AB, CD ; ab_idx: %d, ad_idx: %d,", i, h );
               printf( "\n returned value:\t    ss = %+.7e", ss );
               printf( "\n A = ( %+.5e, %+.5e ) ;  B = ( %+.5e, %+.5e ),"
                  "\n\t\t\t alfa0 = %+.7e rad,", px0, py0, qx0, qy0, alfa0 );
               printf( "\n C = ( %+.5e, %+.5e ) ;  D = ( %+.5e, %+.5e ),"
                  "\n\t\t\t beta1 = %+.7e rad.", rx1, ry1, sx1, sy1, beta1 );
               printf( "\n The actual point index is %ld on cell %ld.]\n",
                  k, lbl.m );
               exit( EXIT_FAILURE );
            }; 

            spt->ppt->cpt->c[k][0] = arc.xx;
            spt->ppt->cpt->c[k][1] = arc.yy;
            spt->ppt->cpt->c[k][2] = pz0 + ss*( qz0 - pz0 );

# else /* if QDL_BNDSPC != 1 */
            rr = ( double ) i/p;
/*............................................................................*/
            ss = subarc( px0, py0, qx0, qy0, alfa0, rr );       /*            */
/*............................................................*/
            spt->ppt->cpt->c[k][0] = arc.rx;
            spt->ppt->cpt->c[k][1] = arc.ry;
            spt->ppt->cpt->c[k][2] = pz0 + rr*( qz0 - pz0 );
# endif 
         };

        set_l:

         l = m;

         if ( i == ONE ) 
            l = qdl.cad[h];
         if ( null < l ) 
            goto set_m;
         if ( h == q ) 
            l = qdl.cdc[i-ONE];
         if ( null < l ) 
            goto set_m;

         SET_CORNER(l);

         if ( *option == 'c' )
         {

# if QDL_BNDSPC == 1
/*............................................................................*/
            ss = arcsct( px1, py1, qx1, qy1, alfa1, rx0, ry0, sx0, sy0, beta0 );
/*............................................................................*/
            if ( ss < -PRECISION )
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling function arcsct(*) !!!" );
               printf( "\n [ Unable to find intersection point of arcs "
                  "AB, CD ; ab_idx: %d, ad_idx: %d,", i, h );
               printf( "\n returned value:\t    ss = %+.7e", ss );
               printf( "\n A = ( %+.5e, %+.5e ) ;  B = ( %+.5e, %+.5e ),"
                  "\n\t\t\t alfa1 = %+.7e rad,", px1, py1, qx1, qy1, alfa1 );
               printf( "\n C = ( %+.5e, %+.5e ) ;  D = ( %+.5e, %+.5e ),"
                  "\n\t\t\t beta0 = %+.7e rad.", rx0, ry0, sx0, sy0, beta0 );
               printf( "\n The actual point index is %ld on cell %ld.]\n",
                  l, lbl.m );
               exit( EXIT_FAILURE );
            }; 

            spt->ppt->cpt->c[l][0] = arc.xx;
            spt->ppt->cpt->c[l][1] = arc.yy;
            spt->ppt->cpt->c[l][2] = pz1 + ss* ( qz1 - pz1 );

# else /* if QDL_BNDSPC != 1 */
            rr = ( double ) ( i-ONE )/p;
/*............................................................................*/
            ss = subarc( px1, py1, qx1, qy1, alfa1, rr );       /*            */
/*............................................................*/
            spt->ppt->cpt->c[l][0] = arc.rx;
            spt->ppt->cpt->c[l][1] = arc.ry;
            spt->ppt->cpt->c[l][2] = pz1 + rr* ( qz1 - pz1 );
# endif
         };

        set_m:

         m = null;

         if ( h == q ) 
            m = qdl.cdc[i];
         if ( null < m ) 
            goto terminal;
         if ( i == p ) 
            m = qdl.cbc[h];
         if ( null < m ) 
            goto terminal;

         SET_CORNER(m);

         if ( *option == 'c' )
         {

# if QDL_BNDSPC == 1
/*............................................................................*/
            ss = arcsct( px1, py1, qx1, qy1, alfa1, rx1, ry1, sx1, sy1, beta1 );
/*............................................................................*/
            if ( ss < -PRECISION )
            {
               printf( "\n\n Message from function %s :", __func__ );
               printf( "\n Error on calling function arcsct(*) !!!" );
               printf( "\n [ Unable to find intersection point of arcs "
                  "AB, CD ; ab_idx: %d, ad_idx: %d,", i, h );
               printf( "\n returned value:\t    ss = %+.7e", ss );
               printf( "\n A = ( %+.5e, %+.5e ) ;  B = ( %+.5e, %+.5e ),"
                  "\n\t\t\t alfa1 = %+.7e rad,", px1, py1, qx1, qy1, alfa1 );
               printf( "\n C = ( %+.5e, %+.5e ) ;  D = ( %+.5e, %+.5e ),"
                  "\n\t\t\t beta1 = %+.7e rad.", rx1, ry1, sx1, sy1, beta1 );
               printf( "\n The actual point index is %ld on cell %ld.]\n",
                  m, lbl.m );
               exit( EXIT_FAILURE );
            }; 

            spt->ppt->cpt->c[m][0] = arc.xx;
            spt->ppt->cpt->c[m][1] = arc.yy;
            spt->ppt->cpt->c[m][2] = pz1 + ss*( qz1 - pz1 );

# else /* if QDL_BNDSPC != 1 */
            rr = ( double ) i/p;
/*............................................................................*/
            ss = subarc( px1, py1, qx1, qy1, alfa1, rr );       /*            */
/*............................................................*/
            spt->ppt->cpt->c[m][0] = arc.rx;
            spt->ppt->cpt->c[m][1] = arc.ry;
            spt->ppt->cpt->c[m][2] = pz1 + rr*( qz1 - pz1 );
# endif
         };

        terminal:

         if ( *option == 't' )
         {
            spt->tpt->cm[lbl.m][0] = j; 
            spt->tpt->cm[lbl.m][1] = k;
            spt->tpt->cm[lbl.m][2] = l;
            spt->tpt->cm[lbl.m][3] = m;

            ind = null - ( qdl.trv == 'e' ) - 2*( qdl.trv == 'm' ); 

# if DSC_ADJCLL == 0
            for ( o=null; o<FACES; o++ )   /* clear neighbouring face identif.*/
               spt->tpt->mn[lbl.m][(int)o] = ind;
# elif DSC_ADJCLL == 1
            for ( o=null ; o<PORTS ; o++ ) /* clear neighbouring port identif.*/
               spt->tpt->mn[lbl.m][(int)o] = ind;
# endif

/* electric and magnetic walls: */ 

            if ( ind == null )
            {
               if (( i == ONE )&&( qdl.ad[h] == 'e' ))  /* face 0 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][0] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][7] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][10] = ELECTRIC_WALL;
# endif
               };

               if (( i == ONE )&&( qdl.ad[h] == 'm' ))
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][0] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][7] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][10] = MAGNETIC_WALL;
# endif
               }; 

               if (( i == p )&&( qdl.bc[h] == 'e' ))    /* face 1 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][1] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][5] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][8] = ELECTRIC_WALL;
# endif
               };

               if (( i == p )&&( qdl.bc[h] == 'm' ))
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][1] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][5] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][8] = MAGNETIC_WALL;
# endif
               };

               if (( h == ONE )&&( qdl.ab[i] == 'e' ))  /* face 2 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][2] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][2] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][11] = ELECTRIC_WALL;
# endif
               };

               if (( h == ONE )&&( qdl.ab[i] == 'm' ))
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][2] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][2] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][11] = MAGNETIC_WALL;
# endif
               };

               if (( h == q )&&( qdl.dc[i] == 'e' ))    /* face 3 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][3] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][0] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][9] = ELECTRIC_WALL;
# endif
               };

               if (( h == q )&&( qdl.dc[i] == 'm' ))
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][3] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][0] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][9] = MAGNETIC_WALL;
# endif
               };

/* top_bottom: */

               if ( qdl.bot == 'e' )                    /* face 4 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][3] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][6] = ELECTRIC_WALL;
# endif
               };

               if ( qdl.bot == 'm' )
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][3] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][6] = MAGNETIC_WALL;
# endif
               };

               if ( qdl.top == 'e' )                    /* face 5 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][5] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][1] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# endif
               };

               if ( qdl.top == 'm' )
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][5] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][1] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# endif
               };
            };/* end if (( qdl.trv != 'e' )&&( qdl.trv != 'm' )) */
         };/* end if ( ...*option == 't'opology ) */  

/* set trivial node material parameters to ZERO: */

         if ( *option == 'c' ) 
         { 
            if (( qdl.trv == 'e' )||( qdl.trv == 'm' ))
               spt->ppt->mpt->idx[lbl.m] = null;

# if QDL_PRINT != null 
            if ( QDL_PRINT == lbl.m )
            {
               printf( "\n\n Intermediate values from function"
                       " '%s' ", __func__ );
               printf( "\n [ generated cell no. %ld ]: " , lbl.m );
               printf( "\n\n j = %ld : ( %.5e, %.5e, %.5e ) ", j,
                  spt->ppt->cpt->c[j][0],
                  spt->ppt->cpt->c[j][1],
                  spt->ppt->cpt->c[j][2] );

               printf( "\n k = %ld : ( %.5e, %.5e, %.5e ) ", k,
                  spt->ppt->cpt->c[k][0],
                  spt->ppt->cpt->c[k][1],
                  spt->ppt->cpt->c[k][2] );

               printf( "\n l = %ld : ( %.5e, %.5e, %.5e ) ", l,
                  spt->ppt->cpt->c[l][0],
                  spt->ppt->cpt->c[l][1],
                  spt->ppt->cpt->c[l][2] );

               printf( "\n m = %ld : ( %.5e, %.5e, %.5e ) ", m,
                  spt->ppt->cpt->c[m][0],
                  spt->ppt->cpt->c[m][1],
                  spt->ppt->cpt->c[m][2] );

               printf( "\n\n [ please acknowledge ] " );
               scanf( "%s", ptr );
            };
# endif
         };/* end if ( *option == 'c'oordinates ) */ 

/* boundaries: */

         if ( i == ONE )
         {
            qdl.sad[h-ONE] = j;
            qdl.sad[h]     = l;
         };

         if ( h == ONE )
         {
            qdl.sab[i-ONE] = j;
            qdl.sab[i]     = k;
         };
           
         if ( i == p )
         {
            qdl.sbc[h-ONE] = k;
            qdl.sbc[h]     = m;
         };
           
         if ( h == q )
         {
            qdl.sdc[i-ONE] = l;
            qdl.sdc[i]     = m;
         };

         rp[i-ONE] = l; /* overtake vertex point[2] for next line h+ONE */

      };/* next i<=p */

      rp[p] = m; /* overtake final vertex point[3] for next line h+ONE */

   };/* next h<=q */

/* reset wall and point indicators, boundaries: */

   for ( h=null ; h<=q ; h++ )
   {
      qdl.wad[h] = qdl.ad[h];
      qdl.wbc[h] = qdl.bc[h];

       qdl.ad[h] = null;
       qdl.bc[h] = null;

      qdl.cad[h] = null;
      qdl.cbc[h] = null;
   };
   for ( i=null ; i<=p ; i++ )  
   {                        
      qdl.wab[i] = qdl.ab[i];
      qdl.wdc[i] = qdl.dc[i];

       qdl.ab[i] = null; 
       qdl.dc[i] = null;

      qdl.cab[i] = null;
      qdl.cdc[i] = null;
   };

   qdl.saab = qdl.aab;
   qdl.sabc = qdl.abc;
   qdl.saad = qdl.aad;
   qdl.sadc = qdl.adc;

   qdl.aab = ZERO; /* reset angles */
   qdl.abc = ZERO;
   qdl.aad = ZERO;
   qdl.adc = ZERO;

   qdl.bot = null; /* reset wall indicators, */
   qdl.top = null; /* top and bottom */
   qdl.trv = null; /* trivial node indicator */

   if ( lbl.pf < lbl.p )
      lbl.pf = lbl.p; 

   if ( lbl.mf < lbl.m )
      lbl.mf = lbl.m;

   lbl.cf = lbl.m;

   return ONE; 
}
/*============================================================================*/
# undef EPSILON
# undef PRECISION
# undef ELECTRIC_WALL
# undef MAGNETIC_WALL
# undef AB_POINTS
# undef BC_POINTS
# undef SET_CORNER
# undef SET_NODE
# undef QDL_PRINT
# undef QDL_BNDSPC
/************************* end of function qudrl-ab ***************************/










# define DO_TRNGL "trngl(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function trngl(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   DSC mesh generation function for triangular domains                        *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
called functions: argz(*), subarc(*), arcsct(*)
in header "../math/arcsct.h"
*/
/*----------------------------------------------------------------------------*/
# define TGL_PRINT 0
/*----------------------------------------------------------------------------*/
# define SET_TRIANGLES 0  /* SET_TRIANGLES 1: generate only triangular cells  */
                          /* [ the number of cells will be n = p*p ].         */
                          /* SET_TRIANGLES 0: generate quadrangular cells     */
                          /* within block, and triangular cells along side bc */
                          /* [ the number of cells will be n = p*(p+1)/2 ].   */
/*----------------------------------------------------------------------------*/
# define AB_POINTS PUZZ_POINTS
/*----------------------------------------------------------------------------*/
# define ELECTRIC_WALL -1
# define MAGNETIC_WALL -2
/*----------------------------------------------------------------------------*/
struct triangle 
{
   double 
    ax,ay,az,bx,by,bz,cx,cy,cz, aab,abc,aac;

   long 
    cab[AB_POINTS],
    cbc[2*AB_POINTS],
    cac[AB_POINTS],
    sab[AB_POINTS],
    sbc[2*AB_POINTS],
    sac[AB_POINTS];

   char  
    ab[AB_POINTS],
    bc[2*AB_POINTS], 
    ac[AB_POINTS];
                      
   char top, bot, trv; /* top, bot: boundary wall indicators  */
                       /* trv: trivial node indicator         */
}; 
struct triangle tgl = {ZERO};
/*----------------------------------------------------------------------------*/
# define SET_CORNER(i) \
{ \
   lbl.p++; \
 \
   if ( CPNTS < lbl.p ) \
   { \
      printf( "\n\n Error message from function %s :", __func__ ); \
      printf( "\n Corner point limit exceeded !!!" ); \
      printf( "\n [ Maximum number is %ld = macro CPNTS " \
         "in %s.", (long) CPNTS, "FORMER.CONF" ); \
      printf( "\n - Change macro only in compliance " \
         "with memory resources.]\n" ); \
      exit( EXIT_FAILURE ); \
   }; \
 \
   if ( lbl.pi <= null ) \
      lbl.pi = lbl.p; \
 \
   (i) = lbl.p; \
}
/*----------------------------------------------------------------------------*/
# define SET_NODE( ) \
{ \
   lbl.m++; \
 \
   if ( NODES < lbl.m ) \
   { \
      printf( "\n\n Error message from function %s :", __func__ ); \
      printf( "\n Mesh cell limit exceeded !!!" ); \
      printf( "\n [ Maximum number is %ld = macro NODES " \
         "in %s.", (long) NODES, "FORMER.CONF" ); \
      printf( "\n - Change macro only in compliance with memory " \
         "resources\n   and same macro in function scattr(*) " \
         "of program SOLVER.C.]\n" ); \
      exit( EXIT_FAILURE ); \
   }; \
 \
   if ( lbl.mi <= null ) \
      lbl.mi = lbl.m; \
}
/*============================================================================*/

short trngl( short p, char *option )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern struct labels lbl;
   extern struct triangle tgl; 
   extern struct subarc arc; 
*/
/* declarartions: */

   static double 
      px0 = ZERO,
      py0 = ZERO,
      pz0 = ZERO,
      px1 = ZERO,
      py1 = ZERO,
      pz1 = ZERO,
      qx0 = ZERO,
      qy0 = ZERO,
      qz0 = ZERO,
      qx1 = ZERO,
      qy1 = ZERO,
      qz1 = ZERO,
      rr = ZERO,
      dbc = ZERO;
     /*  ss = ZERO, */

   static long 
      rp[AB_POINTS] = {null};

   static long
      j = null,
      k = null,
      l = null,
      m = null, 
      n = null;

   static short
      h = null,
      i = null,
      q = null,
      ind = null;

   static char
      o = null;

   double 
      argz( double x, double y, short n ),
      subarc( double x, double y, double z, double w,
                                  double alfa, double r );
/*----------------------------------------------------------------------------*/

   if ( AB_POINTS - ONE < p )
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Too many divisions in triangular block "
         "generation routine !!!" );
      printf( "\n [ Maximum number on side AB is %d = macro AB_POINTS - 1 "
         "in %s.", AB_POINTS - ONE, __func__ );
      printf( "\n - Change macro in compliance with memory resources.]\n" );
      exit( EXIT_FAILURE );
   };

   lbl.ci = lbl.m + ONE;
/*............................................................................*/

   if ( *option == 'c' )         /* copy boundary point coordinates if        */
   {                             /* boundary points are defined [ indexed>0]  */
      if ( null < tgl.cab[0] )
      {
         k = tgl.cab[0];
         tgl.ax = spt->ppt->cpt->c[k][0];
         tgl.ay = spt->ppt->cpt->c[k][1];
         tgl.az = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cab[p] )
      {
         k = tgl.cab[p];
         tgl.bx = spt->ppt->cpt->c[k][0];
         tgl.by = spt->ppt->cpt->c[k][1];
         tgl.bz = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cac[0] )
      {
         k = tgl.cac[0];
         tgl.ax = spt->ppt->cpt->c[k][0];
         tgl.ay = spt->ppt->cpt->c[k][1];
         tgl.az = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cac[p] )
      {
         k = tgl.cac[p];
         tgl.cx = spt->ppt->cpt->c[k][0];
         tgl.cy = spt->ppt->cpt->c[k][1];
         tgl.cz = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cbc[0] )
      {
         k = tgl.cbc[0];
         tgl.bx = spt->ppt->cpt->c[k][0];
         tgl.by = spt->ppt->cpt->c[k][1];
         tgl.bz = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cbc[2*p] )
      {
         k = tgl.cbc[2*p];
         tgl.cx = spt->ppt->cpt->c[k][0];
         tgl.cy = spt->ppt->cpt->c[k][1];
         tgl.cz = spt->ppt->cpt->c[k][2];
      };
   };
/*............................................................................*/
     
   for ( h=null ; h<=p ; h++ )
   {
      rp[h] = null;
   };

   for ( h=ONE ; h<=p ; h++ )
   {
      if ( *option == 'c' )
      {
         rr = ( double ) ( h-ONE )/p;
/*............................................................................*/
         subarc( tgl.ax, tgl.ay, tgl.cx, tgl.cy, tgl.aac, rr );     /*   */
/*.....................................................................*/
         px0 = arc.rx;
         py0 = arc.ry;
         pz0 = tgl.az + rr*( tgl.cz - tgl.az );

/*............................................................................*/
         subarc( tgl.bx, tgl.by, tgl.cx, tgl.cy, tgl.abc, rr );     /*   */
/*.....................................................................*/
         qx0 = arc.rx;
         qy0 = arc.ry;
         qz0 = tgl.bz + rr*( tgl.cz - tgl.bz );

         dbc = rr*tgl.abc;

         rr = ( double ) h/p;
/*............................................................................*/
         subarc( tgl.ax, tgl.ay, tgl.cx, tgl.cy, tgl.aac, rr );     /*   */
/*.....................................................................*/
         px1 = arc.rx;
         py1 = arc.ry;
         pz1 = tgl.az + rr*( tgl.cz- tgl.az );

/*............................................................................*/
         subarc( tgl.bx, tgl.by, tgl.cx, tgl.cy, tgl.abc, rr );     /*   */
/*.....................................................................*/
         qx1 = arc.rx;
         qy1 = arc.ry;
         qz1 = tgl.bz + rr*( tgl.cz - tgl.bz );

         dbc = rr*tgl.abc - dbc; 

# if TGL_PRINT != null  
         if ( TGL_PRINT == lbl.m + ONE )
         {
            printf( "\n\n Intermediate values from function"
                    " '%s' ", __func__ );
            printf( "\n [ generating cell no. %ld ]: " , lbl.m + ONE );
            printf( "\n\n p0 = ( %.5e , %.5e , %.5e )", px0, py0, pz0);
            printf( "\n q0 = ( %.5e , %.5e , %.5e )", qx0, qy0, qz0);
            printf( "\n p1 = ( %.5e , %.5e , %.5e )", px1, py1, pz1);
            printf( "\n q1 = ( %.5e , %.5e , %.5e )", qx1, qy1, qz1);
            printf( "\n\n [ please acknowledge ] " );
            scanf( "%s", ptr);
         };
# endif
      };/* end if *option == 'c'oordinates */

      q = p-h+ONE;
       
      j = null;
      k = null;
      l = null;    
      m = null;
      n = null;
        
      for ( i=ONE; i<=q; i++ )
      {
         SET_NODE( );

     /* set_j: */

         j = k;
         if ( null < j ) 
            goto set_k;
         if ( h == ONE ) 
            j = tgl.cab[i-ONE];
         if ( null < j ) 
            goto set_k;
         if ( i == ONE ) 
            j = tgl.cac[h-ONE];
         if ( null < j ) 
            goto set_k;

         j = rp[i-ONE]; 

         if ( null < j ) 
            goto set_k;

         SET_CORNER(j);

         if ( *option == 'c' )
         {
            rr = ( double ) ( i-ONE )/q;
/*............................................................................*/
            subarc( px0, py0, qx0, qy0, tgl.aab, rr );       /*          */
/*..............................................................*/
            spt->ppt->cpt->c[j][0] = arc.rx;
            spt->ppt->cpt->c[j][1] = arc.ry;
            spt->ppt->cpt->c[j][2] = pz0 + rr*( qz0 - pz0 );
         };

        set_k:
        
         k = null;
         if ( h == ONE ) 
            k = tgl.cab[i];
         if ( null < k ) 
            goto set_l;  
         if ( i == q ) 
            k = tgl.cbc[2*(h-ONE)];
         if ( null < k ) 
            goto set_l;

         k = rp[i];

         if ( k  > null ) 
            goto set_l; 

         SET_CORNER(k);

         if ( *option == 'c' )
         {
            rr = ( double ) i /q;
/*............................................................................*/
            subarc( px0, py0, qx0, qy0, tgl.aab, rr );       /*         */
/*..............................................................*/
            spt->ppt->cpt->c[k][0] = arc.rx;
            spt->ppt->cpt->c[k][1] = arc.ry;
            spt->ppt->cpt->c[k][2] = pz0 + rr*( qz0 - pz0 );
         };
         
        set_l:

         l = m;
         if ( null < l ) 
            goto set_m;
         if ( i == ONE ) 
            l = tgl.cac[h];
         if ( null < l ) 
            goto set_m;
         if ( i == q ) 
            l = tgl.cbc[2*h];
         if ( null < l ) 
            goto set_m;

         SET_CORNER(l);

         if ( *option == 'c' )
         {
            if ( ONE < q )
            {
               rr = ( double ) ( i-ONE )/( q-ONE );
/*............................................................................*/
               subarc( px1, py1, qx1, qy1, tgl.aab, rr );       /*      */
/*.................................................................*/
               spt->ppt->cpt->c[l][0] = arc.rx;
               spt->ppt->cpt->c[l][1] = arc.ry;
               spt->ppt->cpt->c[l][2] = pz1 + rr*( qz1 - pz1 );
            }
            else
            {
               spt->ppt->cpt->c[l][0] = px1;
               spt->ppt->cpt->c[l][1] = py1;
               spt->ppt->cpt->c[l][2] = pz1;
            };
         };

        set_m: 

         m = null;

         if ( i < q )
	 {
	    if ( i == q-ONE ) 
               m = tgl.cbc[2*h];

	    if ( null < m ) 
               goto set_n;

            SET_CORNER(m);

            if( *option == 'c' )
            { 
               rr = ( double ) i /( q-ONE );
/*............................................................................*/
               subarc( px1, py1, qx1, qy1, tgl.aab, rr );       /*       */
/*.................................................................*/
               spt->ppt->cpt->c[m][0] = arc.rx;
               spt->ppt->cpt->c[m][1] = arc.ry;
               spt->ppt->cpt->c[m][2] = pz1 + rr*( qz1 - pz1 );
            };
         };

        set_n:

         if ( i == q )
         {
            n = tgl.cbc[2*h-ONE];

            if ( null < n ) 
               goto terminal;

            SET_CORNER(n);

            if ( *option == 'c' )
            {
               rr = .5;
/*............................................................................*/
               subarc( qx0, qy0 , qx1, qy1, dbc , rr );       /*         */
/*...............................................................*/
               spt->ppt->cpt->c[n][0] = arc.rx;
               spt->ppt->cpt->c[n][1] = arc.ry;
               spt->ppt->cpt->c[n][2] = qz0 + rr*( qz1 - qz0 );
            };
         };

# if SET_TRIANGLES == 1
         if ( i < q )
         {
            SET_CORNER(n);
 
            if ( *option == 'c' )
            {
               spt->ppt->cpt->c[n][0] = .5*( spt->ppt->cpt->c[k][0] +\
                                             spt->ppt->cpt->c[l][0] );
               spt->ppt->cpt->c[n][1] = .5*( spt->ppt->cpt->c[k][1] +\
                                             spt->ppt->cpt->c[l][1] );
               spt->ppt->cpt->c[n][2] = .5*( spt->ppt->cpt->c[k][2] +\
                                             spt->ppt->cpt->c[l][2] );
            };
         };
# endif
  
        terminal:

         if ( *option == 't' )
         {

# if SET_TRIANGLES == 1

            spt->tpt->cm[lbl.m][0] = j;
            spt->tpt->cm[lbl.m][1] = k;
            spt->tpt->cm[lbl.m][2] = l;
            spt->tpt->cm[lbl.m][3] = n;
# else
            if ( i < q ) 
            {
	       spt->tpt->cm[lbl.m][0] = j;
	       spt->tpt->cm[lbl.m][1] = k;
	       spt->tpt->cm[lbl.m][2] = l;
	       spt->tpt->cm[lbl.m][3] = m;
            }
            else if ( i == q )
	    {
	       spt->tpt->cm[lbl.m][0] = j;
	       spt->tpt->cm[lbl.m][1] = k;
	       spt->tpt->cm[lbl.m][2] = l;
	       spt->tpt->cm[lbl.m][3] = n;
            };
# endif
            ind = null - ( tgl.trv == 'e' ) - 2*( tgl.trv == 'm' );

# if DSC_ADJCLL == 0
            for ( o=null; o<FACES; o++ )  /* clear neighbouring face identif.*/
               spt->tpt->mn[lbl.m][(int)o] = ind;
# elif DSC_ADJCLL == 1
            for ( o=null; o<PORTS; o++ ) /* clear neighbouring port identif.*/
               spt->tpt->mn[lbl.m][(int)o] = ind;
# endif

/* electric and magnetic walls: */

            if ( ind == null ) /* non-trivial cell */
            {
               if (( i == ONE )&&( tgl.ac[h] == 'e' ))      /* face 0 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][0] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][7] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][10] = ELECTRIC_WALL;
# endif
               };

               if (( i == ONE )&&( tgl.ac[h] == 'm' ))   
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][0] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][7] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][10] = MAGNETIC_WALL;
# endif
               };
       
               if (( i == q )&&( tgl.bc[2*h-ONE] == 'e' ))  /* face 1 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][1] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][5] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][8] = ELECTRIC_WALL;
# endif
               };

               if (( i == q )&&( tgl.bc[2*h-ONE] == 'm' ))
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][1] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][5] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][8] = MAGNETIC_WALL;
# endif
               };

               if (( h == ONE )&&( tgl.ab[i] == 'e' ))      /* face 2 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][2] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][2] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][11] = ELECTRIC_WALL;
# endif
               };

               if (( h == ONE )&&( tgl.ab[i] == 'm' )) 
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][2] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][2] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][11] = MAGNETIC_WALL;
# endif
               };

               if (( i == q )&&( tgl.bc[2*h] == 'e' ))      /* face 3 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][3] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][0] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][9] = ELECTRIC_WALL;
# endif
               };

               if (( i == q )&&( tgl.bc[2*h] == 'm' ))  
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][3] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][0] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][9] = MAGNETIC_WALL;
# endif
               };

/* top, bottom1: */

               if ( tgl.bot == 'e' )                        /* face 4 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][3] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][6] = ELECTRIC_WALL;
# endif
               };

               if ( tgl.bot == 'm' )                    
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][3] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][6] = MAGNETIC_WALL;
# endif
               };

               if ( tgl.top == 'e' )                       /* face 5 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][5] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][1] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# endif
               };

               if ( tgl.top == 'm' )
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][5] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][1] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# endif
               };
            };/* end if (( tgl.trv != 'e' )&&( tgl.trv != 'm' )) */ 
         };/* end if ( *option == 't'opology ) */

/* set trivial node material parameters to ZERO: */

         if ( *option == 'c' ) 
         {
            if (( tgl.trv == 'e' )||( tgl.trv == 'm' ))
               spt->ppt->mpt->idx[lbl.m] = null;

# if TGL_PRINT != null 
            if ( TGL_PRINT == lbl.m )
            {
               printf( "\n\n Intermediate values from function"
                       " '%s' ", __func__ );
               printf( "\n [ generated cell no. %ld ]: " , lbl.m );
               printf( "\n\n j = %ld : ( %.5e , %.5e , %.5e ) ", j,
                  spt->ppt->cpt->c[j][0],
                  spt->ppt->cpt->c[j][1],
                  spt->ppt->cpt->c[j][2] );

               printf( "\n k = %ld : ( %.5e , %.5e , %.5e ) ", k,
                  spt->ppt->cpt->c[k][0],
                  spt->ppt->cpt->c[k][1],
                  spt->ppt->cpt->c[k][2] );

               printf( "\n l = %ld : ( %.5e , %.5e , %.5e ) ", l,
                  spt->ppt->cpt->c[l][0],
                  spt->ppt->cpt->c[l][1],
                  spt->ppt->cpt->c[l][2] );

# if TGL_TRIANGLES == 1

               printf( "\n n = %ld : ( %.5e , %.5e , %.5e ) ", n,
                  spt->ppt->cpt->c[n][0],
                  spt->ppt->cpt->c[n][1],
                  spt->ppt->cpt->c[n][2] );
# else
               if( i < q )
               {
                  printf( "\n m = %ld : ( %.5e , %.5e , %.5e ) ", m,
                     spt->ppt->cpt->c[m][0],
                     spt->ppt->cpt->c[m][1],
                     spt->ppt->cpt->c[m][2] );
               }
               else if( i == q )
               {
                  printf( "\n n = %ld : ( %.5e , %.5e , %.5e ) ", n,
                     spt->ppt->cpt->c[n][0],
                     spt->ppt->cpt->c[n][1],
                     spt->ppt->cpt->c[n][2] );
               };
# endif
            };
# endif
         };/* end if *option == 'c'oordinates */ 

# if SET_TRIANGLES == 1

         if ( i < q )
         {
            SET_NODE( );

            if ( *option == 't' )
            { 
               spt->tpt->cm[lbl.m][0] = l;
               spt->tpt->cm[lbl.m][1] = n;
               spt->tpt->cm[lbl.m][2] = m;
               spt->tpt->cm[lbl.m][3] = k;

               ind = null - ( tgl.trv == 'e' ) - 2*( tgl.trv == 'm' );

# if DSC_ADJCLL == 1
               for ( o=null; o<=FACES ; o++ )  /* clear neighb. face identf. */
                  spt->tpt->mn[lbl.m][o] = ind;
# elif DSC_ADJCLL == 1
               for ( o=null ; o<=PORTS ; o++ ) /* clear neighb. port identf. */
                  spt->tpt->mn[lbl.m][o] = ind;
# endif
               if ( ind == null )
               {
                /* spt->tpt->bottom2: */

                  if ( tgl.bot == 'e' )                    /* face 4 */
                  {
# if DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][3] = ELECTRIC_WALL;
                     spt->tpt->mn[lbl.m][6] = ELECTRIC_WALL;
# endif
                  };

                  if ( tgl.bot == 'm' )
                  {
# if DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][3] = MAGNETIC_WALL;
                     spt->tpt->mn[lbl.m][6] = MAGNETIC_WALL;
# endif
                  };

                  if ( tgl.top == 'e' )                    /* face 5 */
                  {
# if DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][5] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][1] = ELECTRIC_WALL;
                     spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# endif
                  };

                  if ( tgl.top == 'm' )
                  {
# if DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][5] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][1] = MAGNETIC_WALL;
                     spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# endif
                  };
               };/* end if (( tgl.trv != 'e' )&&( tgl.trv != 'm' )) */
            };/* end if *option == 't'opology */ 

            if ( *option == 'c' )
            {
               if (( tgl.trv == 'e' )||( tgl.trv == 'm' ))
                  spt->ppt->mpt->idx[lbl.m] = null;

# if TGL_PRINT != null 
               if ( TGL_PRINT == lbl.m )
               {
                  printf( "\n\n Intermediate values from function"
                          " '%s' ", __func__ );
                  printf( "\n [ generated cell no. %ld ]: " , lbl.m );
                  printf( "\n\n k = %ld : ( %.5e , %.5e , %.5e ) ", k,
                     spt->ppt->cpt->c[k][0],
                     spt->ppt->cpt->c[k][1],
                     spt->ppt->cpt->c[k][2] );

                  printf( "\n l = %ld : ( %.5e , %.5e , %.5e ) ", l,
                     spt->ppt->cpt->c[l][0],
                     spt->ppt->cpt->c[l][1],
                     spt->ppt->cpt->c[l][2] );

                  printf( "\n m = %ld : ( %.5e , %.5e , %.5e ) ", m,
                     spt->ppt->cpt->c[m][0],
                     spt->ppt->cpt->c[m][1],
                     spt->ppt->cpt->c[m][2] );

                  printf( "\n n = %ld : ( %.5e , %.5e , %.5e ) ", n,
                     spt->ppt->cpt->c[n][0],
                     spt->ppt->cpt->c[n][1],
                     spt->ppt->cpt->c[n][2] );

                  printf( "\n\n [ please acknowledge ] " );
                  scanf( "%s",ptr);
              };
# endif
            };/* end if ( *option == 'c'oordinates */
         };/* end if ( i < q ) */

# endif  /* # endif SET_TRIANGLES == 1 */
      
         if ( h == ONE )
         {
            tgl.sab[i-ONE] = j;
            tgl.sab[i]     = k;
         };
         if ( i == q )
         {
            tgl.sbc[2*(h-ONE)] = k;
            tgl.sbc[2*h-ONE]   = n;
            tgl.sbc[2*h]       = l;
         };
         if ( i == ONE )
         {
            tgl.sac[h-ONE] = j;
            tgl.sac[h]     = l;
         };

         rp[i-ONE] = l;/* overtake vertex point [label 2] for next column h+1 */

      };/* next i */
   };/* next h */ 

/* reset wall and point indicators, boundaries: */

   for ( h=null ; h<=p ; h++ )
   {
      tgl.ab[h] = null; /* boundary ( wall ) indicators */
      tgl.ac[h] = null;
     tgl.cab[h] = null; /* boundary point transfer index */
     tgl.cac[h] = null;
   };
   for ( h=null ; h<=2*p ; h++ )
   {
      tgl.bc[h] = null; /* boundary (- wall ) indicators */
     tgl.cbc[h] = null; /* boundary point transfer index  */
   };

   tgl.aab = ZERO; /* reset angles */
   tgl.aac = ZERO;
   tgl.abc = ZERO;

   tgl.bot = null; /* reset wall indicators, */
   tgl.top = null; /* top and bottom         */
   tgl.trv = null; /* trivial node indicator */

   if ( lbl.pf < lbl.p )
      lbl.pf = lbl.p;

   if ( lbl.mf < lbl.m )
      lbl.mf = lbl.m;
           
   lbl.cf = lbl.m;

   return ONE;
}
/*============================================================================*/
# undef ELECTRIC_WALL
# undef MAGNETIC_WALL
# undef AB_POINTS
# undef SET_TRIANGLES
# undef SET_CORNER
# undef SET_NODE
# undef TGL_PRINT
/************************ end of function trngl(*) ****************************/
