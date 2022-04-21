/* [ file: triads.c ] */
# define DO_TRIADS "triads(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function triads(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   returns cross product                                                      *
*                          p[2] = p[0] ^ p[1]                                  *
*                                                                              *
*        p[2][i] = p[0][i+1]*p[1][i+2] - p[1][i+1]*p[0][i+2] (i+n mod(3))      *
*                                                                              *
*   for given vectors                                                          *
*                   p[0] = ( p[0][0], p[0][1], p[0][2] )                       *
*                   p[1] = ( p[1][0], p[1][1], p[1][2] )                       *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define _POSIX_SOURCE 0 /* set to 1: if POSIX.1 standard headers will be used*/
/*----------------------------------------------------------------------------*/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
# include <time.h>             /* cf. time( ),ctime( ),asctime( ),localtime( )*/
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# define TRD_PRECISION ( 1.e-277 )
# define TRD_RNDOFFBND ( 5.*TRD_PRECISION )
# define TRD_RNDOFF 1
/*----------------------------------------------------------------------------*/
# include "../math/consts.h"
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
/*============================================================================*/

TRIADS *\
lndpdt( TRIADS *tpt )
{
   TRIADS 
     *rtp = &trd;

   static short
      ii = null,
      jj = null;
/*............................................................................*/
   if(( tpt->opt ) == 's' )
   {
      ii = ( tpt->rtn );
      ( rtp->n[ii] ) = ZERO;

      jj = null; do
      {
         ( rtp->vu[ii][jj] ) = ZERO;
         ( rtp->uv[ii][jj] ) = ZERO;
      } while(( ++jj ) < THREE );
   };
   return rtp;
}
/*============================================================================*/

TRIADS *\
triads( TRIADS *tpt )
{
/* declarations: */

   TRIADS *rtp;
   
   static short
      ii = null,
      jj = null,
      kk = null;

   static double
      ss = ZERO,
      nrm = ZERO,
      pr[THREE] = {ZERO};

/* prototypes: - */

/*----------------------------------------------------------------------------*/
/* initialize: */

   if ( tpt == null )
   {
      rtp = &trd;

     init:

      ( rtp->opt ) = null;
      ( rtp->det ) = ZERO;

      ii = null; do
      {
         ( rtp->n[ii] ) = ZERO;
         jj = null; do
         {
            ( rtp->v[ii][jj] ) = ZERO;
            ( rtp->uv[ii][jj] ) = ZERO;
            ( rtp->vu[ii][jj] ) = ZERO;
         } while(( ++jj ) < THREE );
      } while(( ++ii ) < THREE );

      ( rtp->rtn ) = null;
      return rtp;
   }; 
/*............................................................................*/
   rtp = tpt;

   switch( tpt->opt )
   {
/*............................................................................*/
     case 'i': /* option: 'i'nitialize */
     case 'I':

      ( tpt->opt ) = 'i';

      goto init;

      break; /* end ( opt == 'i' ); initialize */
/*............................................................................*/
     case 'c': /* option: 'c'ross [wedge] product */
     case 'C':

      ( tpt->opt ) = 'c';

      ii = null; do
      {
         jj = (( ii + ONE ) % THREE );
         kk = (( ii + TWO ) % THREE );

         ( rtp->v[TWO][ii] ) = \
            (( tpt->v[null][jj] )*( tpt->v[ONE][kk] ) - \
             ( tpt->v[ONE][jj] )*( tpt->v[null][kk] ));
      } while (( ++ii ) < THREE );

      break; /* end ( opt == 'c' ); cross product */
/*............................................................................*/
     case 'd': /* determinant */
     case 'D':

      ( rtp->opt ) = 'd';

      ( rtp->det ) = ZERO;
      ii = null; do /* determinant [ Sarrus' scheme ] */
      {
         jj = (( ii + ONE ) % THREE );
         kk = (( ii + TWO ) % THREE );

         ( rtp->det ) += \
         (( tpt->v[ii][null] )*( tpt->v[jj][ONE] )*( tpt->v[kk][TWO] ));
         ( rtp->det ) -= \
         (( tpt->v[kk][null] )*( tpt->v[jj][ONE] )*( tpt->v[ii][TWO] ));
      } while(( ++ii ) < THREE );

      break; /* end ( opt == 'd' ); determinant */
/*............................................................................*/
     case 's': /* Gram Schmidt orthogonalization: */
     case 'S':
     case 'o':
     case 'O':

      ( rtp->opt ) = 's';

      nrm = ZERO;
      jj = null; do
      {
         ss = tpt->v[null][jj];
         ss *= ss;
         nrm += ss;
      } while(( ++jj ) < THREE );
      nrm = sqrt( nrm );

      if( TRD_PRECISION < nrm )
      {
         ( rtp->n[null] ) = nrm;
         ( rtp->vu[null][null] ) = nrm;

         jj = null; do
         {
            ( rtp->uv[null][jj] ) = ( tpt->v[null][jj] ) / nrm;

            if ( null < jj )
               ( rtp->vu[null][jj] ) = ZERO;
         } while(( ++jj ) < THREE );
      }
      else
      {
         ( rtp->rtn ) = ii;
/*............................................................................*/
         rtp = lndpdt( rtp );        /* linearly dependent vectors v[ii][*]   */
/*.................................*/
      };

      ( rtp->det ) = ( rtp->vu[null][null] );

      ii = ONE; do
      {
         jj = null; do
         {
            ( rtp->uv[ii][jj] ) = ( tpt->v[ii][jj] );
         } while(( ++jj ) < THREE );

         kk = null; do
         {
            pr[kk] = ZERO; /* projection of vector v onto uv[kk] */
            jj = null; do
            {
               pr[kk] += (( tpt->v[ii][jj] )*( rtp->uv[kk][jj] ));
            } while(( ++jj ) < THREE );
            ( rtp->vu[ii][kk] ) = pr[kk];

            jj = null; do
            {
               ( rtp->uv[ii][jj] ) -= ( pr[kk] * ( tpt->uv[kk][jj] ));
            } while(( ++jj ) < THREE );
         } while(( ++kk ) < ii );

         nrm = ZERO;
         jj = null; do
         {
            ss = ( rtp->uv[ii][jj] );
            ss *= ss;
            nrm += ss;
         } while(( ++jj ) < THREE );
         nrm = sqrt( nrm );

         if ( TRD_PRECISION < nrm )
         {
            jj = null; do
            {
               ( rtp->uv[ii][jj] ) /= nrm;
            } while(( ++jj ) < THREE );
            ( rtp->vu[ii][ii] ) = nrm;
         }
         else
         {
            ( rtp->rtn ) = ii;
/*............................................................................*/
            rtp = lndpdt( rtp );     /* linearly dependent vectors v[ii][*]   */
/*.................................*/
         };

         ( rtp->det ) *= ( rtp->vu[ii][ii] );

         jj = ( ++ii );
         while( jj < THREE )
         {
            ( rtp->vu[ii][jj] ) = ZERO;
            jj++;
         };
      } while( ii < THREE );

      break; /* end ( opt == 's' ); Gram-Smith orthonormalization */
/*............................................................................*/

     default: 

      fprintf( stderr, "\n\n Error message from function %s:", DO_TRIADS );
      fprintf( stderr, "\n Unknown or unspecified option !!!" ); 
      fprintf( stderr, "\n [ Usage: opt = 'c' [cross_product] "
         "or opt = 'd' [determinant], e.g.]\n" );

      exit( EXIT_FAILURE );

      break;
   }; /* end switch( opt ) */
/*............................................................................*/
# if TRD_RNDOFF == 1

   if( fabs( rtp->det ) < TRD_RNDOFFBND )
      ( rtp->det ) = ZERO;

   ii = null; do
   {
      if( fabs( rtp->n[ii] ) < TRD_RNDOFFBND )
         ( rtp->n[ii] ) = ZERO;

      jj = null; do
      {
         if( fabs( rtp->uv[ii][jj] ) < TRD_RNDOFFBND )
            ( rtp->uv[ii][jj] ) = ZERO;
         if( fabs( rtp->vu[ii][jj] ) < TRD_RNDOFFBND )
            ( rtp->vu[ii][jj] ) = ZERO;
      } while(( ++jj ) < THREE );
   } while(( ++ii ) < THREE );

# endif

   ( rtp->rtn ) = null;
   return rtp;
}
/*============================================================================*/
# undef TRD_PRECISION
# undef TRD_RNDOFFBND
# undef TRD_RNDOFF
/************************* end of function triads(*) **************************/
