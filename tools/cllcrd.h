/* [ file: cllcrd.h ] */
# define DO_CLLCRD "cllcrd(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function cllcrd(*), DANSE release 1.0.                              *
*   This function determines the center coordinates of a DSC cell [ labelled   *
*   by that variable ], and, optionally, the center coordinates of a face or   *
*   port of that cell, the domain of indices being [ for macros defined in     *
*   function linker(*) ]:                                                      *
*                                                                              *
*      1 <= cell <= NODES;  0 <= face < FACES;  1 <= port <= PORTS             *
*                                                                              *
*   In addition, the norm [Euclidean length] of the port vector is computed.   *
*   All values are written into a structure CLLCRD, a pointer to which is      *
*   returned.                                                                  *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "cllctp.h"
/*----------------------------------------------------------------------------*/
/*
typedef struct
{
   double xn, yn, zn, xf, yf, zf, xp, yp, zp, px, py, pz, pm, fx, fy, fz, fm;
} CLLCRD;
*/
/*----------------------------------------------------------------------------*/
typedef struct
{
   short
      rtn;

   char
      opt;

   double
      det,
      n[THREE],
      v[THREE][THREE],
      uv[THREE][THREE],
      vu[THREE][THREE];

} TRIADS;
static TRIADS trd = {null};
/*----------------------------------------------------------------------------*/
# define CLL_FACE( ) \
{ \
   ii = null; do \
   { \
      fce[ii] = ZERO; \
 \
      kk = null; do \
      { \
         nn = ( spt->tpt->cm[cell][vtx[kk]] ); \
         edge[kk][ii] = ( spt->ppt->cpt->c[nn][ii] ); \
         kk++ ; \
      } while ( kk < FOUR ); \
      edge[FOUR][ii] = edge[null][ii]; \
      do \
      { \
         edge[kk][ii] -= edge[kk-ONE][ii]; \
         kk--; \
      } while( null < kk ); \
      edge[null][ii] = edge[FOUR][ii]; \
 \
      ii++ ; \
   } while ( ii < THREE ); \
 \
   kk = null; do \
   { \
      ii = null; do \
      { \
         jj = (( kk+ONE ) % FOUR ); \
         ( trp->v[null][ii] ) = edge[kk][ii]; \
         ( trp->v[ONE][ii] ) = edge[jj][ii]; \
         ii++ ; \
      } while( ii < THREE ); \
 \
      ( trp->opt ) = 'c'; /* option: 'c'ross product */\
/*..............................*/ \
      trp = triads( trp );      /* call triads(*) in option 'c'ross_product */ \
/*............................*/ \
      ii = null; do \
      { \
         fce[ii] += ( trp->v[TWO][ii] ); \
         ii++ ; \
      } while( ii < THREE ); \
      kk++ ; \
   } while( kk < FOUR ); \
 \
   ii = null; do \
   { \
      fce[ii] /= FOUR; \
      ii++ ; \
   } while( ii < THREE ); \
}
/*============================================================================*/

CLLCRD *\
cllcrd( long cell, char face, char port )
{
/* allusions; */
/*
   extern FORMSTATE *spt;
*/
/* declarations: */

   static CLLCRD 
      crd = { null },
     *str = &crd;

   static TRIADS
     *trp = &trd;

   static long
      nn = null;

   static short
      ii = null,
      jj = null,
      kk = null,
      vtx[FOUR] = {null};

   static double
      xx = ZERO,
      yy = ZERO,
      zz = ZERO;

   static double
      fce[THREE] = {ZERO},
      edge[FIVE][THREE] = {{ZERO}};

/* prototypes: */

   TRIADS
      *triads( TRIADS *trp );

   double sqrt( double x ); 
   
/*----------------------------------------------------------------------------*/
/* center of cell [ 'node' center coordinates ]: */

   xx = ZERO;
   yy = ZERO;
   zz = ZERO;

   ii = null;
   do 
   {
      xx += spt->ppt->cpt->c[spt->tpt->cm[cell][ii]][null];
      yy += spt->ppt->cpt->c[spt->tpt->cm[cell][ii]][ONE];
      zz += spt->ppt->cpt->c[spt->tpt->cm[cell][ii]][TWO];
      ii++;
   } while ( ii < CRNRS );

/* the node [center] coordinates: */

   crd.xn = xx/CRNRS;
   crd.yn = yy/CRNRS;
   crd.zn = zz/CRNRS;
   
/*----------------------------------------------------------------------------*/
/* determine face vertices: */

   switch( face )
   {
     case 0: 
      vtx[0] = 0;
      vtx[1] = 2;
      vtx[2] = 6;
      vtx[3] = 4;
      break;
     
     case 1:
      vtx[0] = 1;
      vtx[1] = 3;
      vtx[2] = 7;
      vtx[3] = 5;
      break;

     case 2:
      vtx[0] = 0;
      vtx[1] = 4;
      vtx[2] = 5;
      vtx[3] = 1;
      break;

     case 3:
      vtx[0] = 2;
      vtx[1] = 6;
      vtx[2] = 7;
      vtx[3] = 3;
      break;

     case 4:
      vtx[0] = 0;
      vtx[1] = 1;
      vtx[2] = 3;
      vtx[3] = 2;
      break;

     case 5:
      vtx[0] = 4;
      vtx[1] = 5;
      vtx[2] = 7;
      vtx[3] = 6;
      break;
   };
/*............................................................................*/
/* center of face: */

   xx = ZERO;
   yy = ZERO;
   zz = ZERO;
   
   ii = null; do
   {
      xx += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][null] );
      yy += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][ONE] );
      zz += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][TWO] );
      ii++;
   } while ( ii < FOUR );
/*............................................................................*/
/* the face [center] coordinates: */

   crd.xf = xx/FOUR;
   crd.yf = yy/FOUR;
   crd.zf = zz/FOUR;
/*............................................................................*/
/* the face [surface] vector: */

   CLL_FACE( );

   crd.fx = fce[null];
   crd.fy = fce[ONE];
   crd.fz = fce[TWO];
   
   crd.fm = sqrt( crd.fx*crd.fx + crd.fy*crd.fy + crd.fz*crd.fz );

   if( port < ONE )
      return str;

/*----------------------------------------------------------------------------*/
/* determine port vectors: */

   switch( port )
   {
     case 1:
      vtx[0] = 3;
      vtx[1] = 7;
      vtx[2] = 2;
      vtx[3] = 6;
      break;

     case 2:
      vtx[0] = 5;
      vtx[1] = 7;
      vtx[2] = 4;
      vtx[3] = 6;
      break;

    case 3:
      vtx[0] = 1;
      vtx[1] = 5;
      vtx[2] = 0;
      vtx[3] = 4;
      break;

     case 4:
      vtx[0] = 1;
      vtx[1] = 3;
      vtx[2] = 0;
      vtx[3] = 2;
      break;

     case 5:
      vtx[0] = 6;
      vtx[1] = 7;
      vtx[2] = 4;
      vtx[3] = 5;
      break;

     case 6:
      vtx[0] = 3;
      vtx[1] = 7;
      vtx[2] = 1;
      vtx[3] = 5;
      break;
     
     case 7:
      vtx[0] = 2;
      vtx[1] = 3;
      vtx[2] = 0;
      vtx[3] = 1;
      break;

     case 8:
      vtx[0] = 2;
      vtx[1] = 6;
      vtx[2] = 0;
      vtx[3] = 4;
      break;

    case 9:
      vtx[0] = 5;
      vtx[1] = 7;
      vtx[2] = 1;
      vtx[3] = 3;
      break;

     case 10:
      vtx[0] = 6;
      vtx[1] = 7;
      vtx[2] = 2;
      vtx[3] = 3;
      break;

     case 11:
      vtx[0] = 4;
      vtx[1] = 6;
      vtx[2] = 0;
      vtx[3] = 2;
      break;

     case 12:
      vtx[0] = 4;
      vtx[1] = 5;
      vtx[2] = 0;
      vtx[3] = 1;
      break;
   };
/*............................................................................*/
/* center of port: */

   xx = ZERO;
   yy = ZERO;
   zz = ZERO;

   ii = null; do
   {
      xx += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][null] );
      yy += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][ONE] );
      zz += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][TWO] );
      ii++;
   } while ( ii < FOUR );

/* the port [center] coordinates: */

   crd.xp = xx/FOUR;
   crd.yp = yy/FOUR;
   crd.zp = zz/FOUR;

/*............................................................................*/
/* port vector: */

   xx = ZERO;
   yy = ZERO;
   zz = ZERO;

   ii = null; do
   {
      xx += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][null] );
      yy += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][ONE] );
      zz += ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][TWO] );
      ii++;
   } while ( ii < TWO );
   do
   {
      xx -= ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][null] );
      yy -= ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][ONE] );
      zz -= ( spt->ppt->cpt->c[spt->tpt->cm[cell][vtx[ii]]][TWO] );
      ii++;
   } while ( ii < FOUR );

   xx /= TWO;
   yy /= TWO;
   zz /= TWO;

/* the port vector components */

   crd.px = xx;
   crd.py = yy;
   crd.pz = zz;

   crd.pm = sqrt( xx*xx + yy*yy + zz*zz );
   
   return str;
}
/*============================================================================*/
/*********************** end of function cllcrd(*) ****************************/
