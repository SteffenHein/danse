/* [ file: cshape.c ] */
# define DO_CSHAPE "cshape(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function cshape(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given eigth points in 3-space, crp->c[i][] ( i=0,...,7 ), this function    *
*   returns geometric ['cell shape'] data, such als port and face [vectors],   *
*   their lengths and spatial orientations and positions, characterizing the   *
*   hexahedral cell with these vertex points crp->c[i]                         *
*                                                                              *
*   All parameters are transferred to and returned from this function in       *
*   struct cll of type CSHAPE.                                                 *
*                                                                              *
*   (C) SHEIN; Bad Aibling, February 2007                     Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <string.h>
# include <stdio.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# include "../CONFIG.H"
# include "../former/FORMER.CONF"
/*----------------------------------------------------------------------------*/
# define CLL_INCLUDE 1 /* 1: include parameter transfer struct of type CSHAPE */
# define CLL_DEBUG   0 /*>0: activate various printing functions [ for debug- */
                       /*    ging purposes, mainly ] */
/*----------------------------------------------------------------------------*/
/* operation marks: */

# define CLL_FULPIVT 1       /* 1: fully pivoted Gauss-Jordan algorithm       */
                             /*    [ more precise but CPU time expensive ]    */

# define CLL_FCEMODE 1       /* must be 1 or 2:                               */
                             /* 1: compute face vectors from cell ports       */
                             /* 2: compute face vectors from cell edges       */

# define CLL_DSPINTM 0       /* CLL_DSPINTM 1: display intermediate results   */
/*----------------------------------------------------------------------------*/
/* computational parameters, regularization factors, bounds etc.: */

# define CLL_PRECISION ( double )( 1.000e-15 )
# define CLL_TRIVIAL   ( double )( 1.e-277 ) /* trivial media (ke,km)         */
/*----------------------------------------------------------------------------*/
/* array dimensions: */

# define DIMNS  3
# define CRNRS  8
# define SORDR  6
# define PORTS 12
# define FACES  6
/*----------------------------------------------------------------------------*/
/* natural constants - to be changed only with the help of God:               */
/* [ all units are international units (mks), if not specified otherwise.]    */

# define EPS_VAC      ( 8.8541878170e-12 ) /* vac. permittivity [A*sec/(V*m)] */
# define MY_VAC_      ( 1.2566370614e-06 ) /* "    permeability [V*sec/(A*m)] */
# define ELECTR_CHRGE (-1.6021773349e-19 ) /* electron charge [Coulomb=A*sec] */
# define ELECTR_MASS_ ( 9.1093897540e-31 ) /* electron mass [kg]              */
# define BOHRs_MAGNTN ( 1.1654071500e-29 ) /* Bohr's magneton [Volt*sec*m]    */
# define PLANCKs_CNST ( 6.6260755400e-34 ) /* Planck's constant [Joule*sec]   */
/*----------------------------------------------------------------------------*/
# include "../former/formtp.h"
# include "../math/trdstp.h"
# include "../former/gssjtp.h"
# include "../former/jacbtp.h"
/*---------------------------------------------------------------------------*/
# if CLL_INCLUDE == 1
   # include "cshptp.h"
# else
/*-----------------------------------------------------------------------------
*//* geometrical cell parameter transfer structure, release 5.1.
*//* All units are international units (mks), if not otherwise specified.
------------------------------------------------------------------------------*/
typedef struct
{
   signed char /* any return character [ may be used as an error flag, e.g. ] */
      rtn;

   signed char /* any option */
      opt;

   char
      skew; /* geometric skew indicator [ no skew:0 / skew:1] */

/*............................................................................*/
/* input parameters: */

   double
      c[CRNRS][DIMNS]; /* corner point coordinates */

   long
      cell;

   char
      face, port;

   TOPOLOGY *tpt;

   PARAMETERS *ppt;
/*............................................................................*/
/* computed geometric parameters: */

   double
      vol, /* cell volume [ m^3 ] */
      skv, /* cell skew volume [ m^3 ] */
      xn, yn, zn;

   double
      px[PORTS],
      py[PORTS],
      pz[PORTS],
      pm[PORTS],

      xp[PORTS],
      yp[PORTS],
      zp[PORTS],

      fx[FACES],
      fy[FACES],
      fz[FACES],
      fm[FACES],

      xf[FACES],
      yf[FACES],
      zf[FACES],

      vp[FACES]; /* face pyramide volume: vp[jj] = ( f[jj] | b[jj/2] )/6 */

   double
      e[PORTS][DIMNS], /* edge vectors */
     eu[PORTS][DIMNS], /* edge vectors relative to cell basis ub[] */

      p[PORTS][DIMNS], /* port vectors */
     pu[PORTS][DIMNS], /* port vectors relative to cell basis ub[] */

      f[FACES][DIMNS], /* face vectors; inner or outer - depending on */
                       /* face label !!! */
     nf[FACES][DIMNS], /* outer face normal vectors */
     fu[FACES][DIMNS], /* face vectors relative to cell basis ub[] */
     fa[DIMNS][DIMNS], /* opposite face vector arithmetic means */

      a[DIMNS][DIMNS], /* area vector matrix A[i][j] = adj((arv[i])[j]) */
     au[DIMNS][DIMNS], /* area vector matrix relative to cell basis ub[] */
     ai[DIMNS][DIMNS], /* Inverse are vector matrix A^-1 */

      b[DIMNS][DIMNS], /* node vector matrix B[i][j] = adj((ndv[i])[j]) */
     ub[DIMNS][DIMNS], /* Gram-Schmidt orthonormalzd node vectors [cell basis]*/
     bu[DIMNS][DIMNS], /* normalized Node vector matrix relative to ub[] */
     bi[DIMNS][DIMNS], /* inverse node vector basis B^-1 */

     cs[DIMNS][DIMNS], /* cell shape skew vectors */
     cu[DIMNS][DIMNS], /* cell shape skew vectors relative to cell basis */

     uv[DIMNS][DIMNS], /* any orthonormal basis */
     vu[DIMNS][DIMNS]; /* any coordinate vectors [ with respect to on ] */

} CSHAPE;
/*----------------------------------------------------------------------------*/
# endif /* CLL_INCLUDE != 1 */
/*----------------------------------------------------------------------------*/
/* system function prototypes: */

void abort( void );
int printf( const char *format,...);
int scanf( const char *format,...);
/*----------------------------------------------------------------------------*/
/* macros: */
/*----------------------------------------------------------------------------*/
/* complex number product: UU+i*VV = (AA+i*BB)*(CC+i*DD): */
# define CPRODUCT( AA, BB, CC, DD, UU, VV ) \
{ \
   (UU) = (AA)*(CC) - (BB)*(DD); \
   (VV) = (BB)*(CC) + (AA)*(DD); \
}

/* complex number quotient: UU+i*VV = (AA+i*BB)/(CC+i*DD): */
# define CQUOTIENT( AA, BB, CC, DD, UU, VV ) \
{ \
   (VV) = (CC)*(CC) + (DD)*(DD); \
   (UU) = ((AA)*(CC) + (BB)*(DD))/(VV); \
   (VV) = ((BB)*(CC) - (AA)*(DD))/(VV); \
}
/*----------------------------------------------------------------------------*/
/* compute mesh cell center: */
# define CLL_CENTER( ) \
{ \
   xx = ZERO; \
   yy = ZERO; \
   zz = ZERO; \
 \
   ii = null; do \
   {  \
      xx += ( rtp->c[ii][null] ); \
      yy += ( rtp->c[ii][ONE] ); \
      zz += ( rtp->c[ii][TWO] ); \
   } while (( ++ii ) < CRNRS ); \
 \
/* the node [viz. mesh cell center] coordinates: */ \
 \
   ( rtp->xn ) = xx/CRNRS; \
   ( rtp->yn ) = yy/CRNRS; \
   ( rtp->zn ) = zz/CRNRS; \
} /* end of macro CLL_CENTER( ) */
/*----------------------------------------------------------------------------*/
/* compute mesh cell port vectors from vertex points: */
# define CLL_PORTS( ) \
{ \
   ii = null; do \
   { \
      xx = ZERO; \
      yy = ZERO; \
      zz = ZERO; \
 \
      kk = null; do \
      { \
         ( rtp->p[ii][kk] ) = ZERO; \
      } while (( ++kk ) < THREE ); \
 \
      jj = null; do \
      { \
         kk = null; do \
         { \
            ( rtp->p[ii][kk] ) += .5*( rtp->c[(short)pvx[ii][jj]][kk] ); \
            ( rtp->p[ii][kk] ) -= .5*( rtp->c[(short)pvx[ii][jj+TWO]][kk] ); \
         } while (( ++kk ) < THREE ); \
 \
         xx += ( rtp->c[(short)pvx[ii][jj]][null] ); \
         xx += ( rtp->c[(short)pvx[ii][jj+TWO]][null] ); \
         yy += ( rtp->c[(short)pvx[ii][jj]][ONE] ); \
         yy += ( rtp->c[(short)pvx[ii][jj+TWO]][ONE] ); \
         zz += ( rtp->c[(short)pvx[ii][jj]][TWO] ); \
         zz += ( rtp->c[(short)pvx[ii][jj+TWO]][TWO] ); \
 \
      } while (( ++jj ) < TWO ); \
 \
      ( rtp->xp[ii] ) = xx*QUART; \
      ( rtp->yp[ii] ) = yy*QUART; \
      ( rtp->zp[ii] ) = zz*QUART; \
 \
      ( rtp->px[ii] ) = ( rtp->p[ii][null] ); \
      ( rtp->py[ii] ) = ( rtp->p[ii][ONE] ); \
      ( rtp->pz[ii] ) = ( rtp->p[ii][TWO] ); \
 \
      zz = ZERO; \
      jj = null; do \
      { \
         xx = ( rtp->p[ii][jj] ); \
         zz += ( xx*xx ); \
      } while (( ++jj ) < THREE ); \
      ( rtp->pm[ii] ) = sqrt( zz ); \
   } while (( ++ii ) < PORTS ); \
} /* end of macro CLL_PORTS( ) */
/*----------------------------------------------------------------------------*/
/* compute mesh cell edge vectors from vertex points: */

# define CLL_EDGES( ) \
{ \
   jj = null; do \
   { \
      ( rtp->e[0][jj] ) = ( rtp->c[3][jj] - rtp->c[2][jj] ); \
      ( rtp->e[1][jj] ) = ( rtp->c[7][jj] - rtp->c[6][jj] ); \
      ( rtp->e[2][jj] ) = ( rtp->c[5][jj] - rtp->c[4][jj] ); \
      ( rtp->e[3][jj] ) = ( rtp->c[1][jj] - rtp->c[0][jj] ); \
 \
      ( rtp->e[4][jj] ) = ( rtp->c[6][jj] - rtp->c[4][jj] ); \
      ( rtp->e[5][jj] ) = ( rtp->c[7][jj] - rtp->c[5][jj] ); \
      ( rtp->e[6][jj] ) = ( rtp->c[3][jj] - rtp->c[1][jj] ); \
      ( rtp->e[7][jj] ) = ( rtp->c[2][jj] - rtp->c[0][jj] ); \
 \
      ( rtp->e[8][jj] ) = ( rtp->c[5][jj] - rtp->c[1][jj] ); \
      ( rtp->e[9][jj] ) = ( rtp->c[7][jj] - rtp->c[3][jj] ); \
      ( rtp->e[10][jj] ) = ( rtp->c[6][jj] - rtp->c[2][jj] ); \
      ( rtp->e[11][jj] ) = ( rtp->c[4][jj] - rtp->c[0][jj] ); \
 \
   } while (( ++jj ) < DIMNS ); \
} /* end of macro CLL_EDGES( ) */
/*----------------------------------------------------------------------------*/
/* orthonormalize node vectors b[] to a 'cell basis' ub[], and transform node */
/* vectors into cell basis coordinates; this yields vectors bu[]: */
# define CLL_ORTNORM( ) \
{ \
   ii = null; do \
   { \
      jj = null; do \
      { \
         ( trp->v[ii][jj] ) = ( rtp->b[ii][jj] ); \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < THREE ); \
 \
   ( trp->opt ) = 'o'; /* option: 'o'rthonormalize */\
/*..........................................................................*/ \
   trp = triads( trp );         /* call triads(*), option 'o'rthonormlze */ \
/*............................*/ \
   ii = null; do \
   { \
      jj = null; do \
      { \
         ( rtp->ub[ii][jj] ) = ( trp->uv[ii][jj] ); \
         ( rtp->bu[ii][jj] ) = ( trp->vu[ii][jj] ); \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < THREE ); \
}
/*----------------------------------------------------------------------------*/
/* There are different ways to compute the cell face vectors [FCEMODE 1,2]: */
# if CLL_FCEMODE == 1
# define CLL_CLFCES( ) \
{ \
   ii = null; do \
   { \
      xx = ZERO; \
      yy = ZERO; \
      zz = ZERO; \
 \
      switch(ii) \
      { \
        case 0: \
         kk = 7; \
         ll = 10; \
         break ; \
 \
        case 1: \
         kk = 5; \
         ll = 8; \
         break ; \
 \
        case 2: \
         kk = 11; \
         ll = 2; \
         break ; \
 \
        case 3: \
         kk = 9; \
         ll = 0; \
         break ; \
 \
        case 4: \
         kk = 3; \
         ll = 6; \
         break ; \
 \
        case 5: \
         kk = 1; \
         ll = 4; \
         break ; \
      }; \
 \
      jj = null; do \
      { \
         ( trp->v[null][jj] ) = ( rtp->p[kk][jj] ); \
         ( trp->v[ONE][jj] ) = ( rtp->p[ll][jj] ); \
      } while (( ++jj ) < THREE ); \
 \
      ( trp->opt ) = 'c'; /* option: 'c'ross product */\
/*..........................................................................*/ \
      trp = triads( trp );      /* call triads(*) in option 'c'ross_product */ \
/*............................*/ \
      jj = null; do \
      { \
         ( rtp->f[ii][jj] ) = ( trp->v[TWO][jj] ); \
      } while (( ++jj ) < THREE ); \
 \
      ( rtp->fx[ii] ) = ( rtp->f[ii][null] ); \
      ( rtp->fy[ii] ) = ( rtp->f[ii][ONE] ); \
      ( rtp->fz[ii] ) = ( rtp->f[ii][TWO] ); \
 \
      jj = null; do \
      { \
         xx += ( rtp->c[(short)fvx[ii][jj]][null] ); \
         yy += ( rtp->c[(short)fvx[ii][jj]][ONE] ); \
         zz += ( rtp->c[(short)fvx[ii][jj]][TWO] ); \
      } while (( ++jj ) < FOUR ); \
 \
      ( rtp->xf[ii] ) = xx*QUART; \
      ( rtp->yf[ii] ) = yy*QUART; \
      ( rtp->zf[ii] ) = zz*QUART; \
 \
      zz = ZERO; \
      jj = null; do \
      { \
         xx = ( rtp->f[ii][jj] ); \
         zz += ( xx*xx ); \
      } while (( ++jj ) < THREE ); \
      zz = sqrt( zz ); \
      ( rtp->fm[ii] ) = zz; \
 \
      jj = null; do \
      { \
         ( rtp->nf[ii][jj] ) = \
            (( double )( TWO*( ii%TWO ) - ONE )*( rtp->f[ii][jj] )/zz ); \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < FACES ); \
} /* end of macro CLL_CLFCES */
/*----------------------------------------------------------------------------*/
# elif CLL_FCEMODE == 2
# define CLL_CLFCES( ) \
{ \
   ii = null; do \
   { \
      xx = ZERO; \
      yy = ZERO; \
      zz = ZERO; \
 \
      switch(ii) \
      { \
        case 0: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->e[10][jj] ); \
            edge[1][jj] = - ( rtp->e[4][jj] ); \
            edge[2][jj] = - ( rtp->e[11][jj] ); \
            edge[3][jj] = ( rtp->e[7][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
 \
        case 1: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->e[9][jj] ); \
            edge[1][jj] = - ( rtp->e[5][jj] ); \
            edge[2][jj] = - ( rtp->e[8][jj] ); \
            edge[3][jj] = ( rtp->e[6][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
 \
        case 2: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->e[2][jj] ); \
            edge[1][jj] = - ( rtp->e[8][jj] ); \
            edge[2][jj] = - ( rtp->e[3][jj] ); \
            edge[3][jj] = ( rtp->e[11][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
 \
        case 3: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->e[1][jj] ); \
            edge[1][jj] = - ( rtp->e[9][jj] ); \
            edge[2][jj] = - ( rtp->e[0][jj] ); \
            edge[3][jj] = ( rtp->e[10][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
 \
        case 4: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->e[6][jj] ); \
            edge[1][jj] = - ( rtp->e[0][jj] ); \
            edge[2][jj] = - ( rtp->e[7][jj] ); \
            edge[3][jj] = ( rtp->e[3][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
 \
        case 5: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->e[5][jj] ); \
            edge[1][jj] = - ( rtp->e[1][jj] ); \
            edge[2][jj] = - ( rtp->e[4][jj] ); \
            edge[3][jj] = ( rtp->e[2][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
      }; \
 \
      jj = null; do \
      { \
         ( rtp->f[ii][jj] ) = ZERO; \
      } while (( ++jj ) < THREE ); \
 \
      kk = null; do \
      { \
         jj = null; do \
         { \
            ll = (( kk+ONE )%FOUR ); \
            ( trp->v[null][jj] ) = edge[kk][jj]; \
            ( trp->v[ONE][jj]  ) = edge[ll][jj]; \
         } while (( ++jj ) < THREE ); \
 \
         ( trp->opt ) = 'c'; /* option: 'c'ross product */\
/*..........................................................................*/ \
         trp = triads( trp );   /* call triads(*) in option 'c'ross_product */ \
/*............................*/ \
         jj = null; do \
         { \
            ( rtp->f[ii][jj] ) += ( trp->v[TWO][jj] ); \
         } while (( ++jj ) < THREE ); \
 \
         xx += ( rtp->c[fvx[ii][kk]][null] ); \
         yy += ( rtp->c[fvx[ii][kk]][ONE] ); \
         zz += ( rtp->c[fvx[ii][kk]][TWO] ); \
 \
      } while (( ++kk ) < FOUR ); \
 \
      jj = null; do \
      { \
         ( rtp->f[ii][jj] ) *= QUART; \
      } while (( ++jj ) < THREE ); \
 \
      ( rtp->fx[ii] ) = ( rtp->f[ii][null] ); \
      ( rtp->fy[ii] ) = ( rtp->f[ii][ONE] ); \
      ( rtp->fz[ii] ) = ( rtp->f[ii][TWO] ); \
 \
      ( rtp->xf[ii] ) = xx*QUART; \
      ( rtp->yf[ii] ) = yy*QUART; \
      ( rtp->zf[ii] ) = zz*QUART; \
 \
      zz = ZERO; \
      jj = null; do \
      { \
         xx = ( rtp->f[ii][jj] ); \
         zz += ( xx*xx ); \
      } while (( ++jj ) < THREE ); \
      ( rtp->fm[ii] ) = sqrt( zz ); \
 \
   } while (( ++ii ) < FACES ); \
} /* end of macro CLL_CLFCES */
# endif /* CLL_FCEMODE == 2 */
/*----------------------------------------------------------------------------*/
# define CLL_CLSKEW( ) \
{ \
   ii = null; do \
   { \
      switch(ii) \
      { \
        case 0: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->p[0][jj] ); \
            edge[1][jj] = ( rtp->p[1][jj] ); \
            edge[2][jj] = ( rtp->p[2][jj] ); \
            edge[3][jj] = ( rtp->p[3][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
 \
        case 1: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->p[4][jj] ); \
            edge[1][jj] = ( rtp->p[5][jj] ); \
            edge[2][jj] = ( rtp->p[6][jj] ); \
            edge[3][jj] = ( rtp->p[7][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
 \
        case 2: \
         jj = null; do \
         { \
            edge[0][jj] = ( rtp->p[8][jj] ); \
            edge[1][jj] = ( rtp->p[9][jj] ); \
            edge[2][jj] = ( rtp->p[10][jj] ); \
            edge[3][jj] = ( rtp->p[11][jj] ); \
         } while (( ++jj ) < THREE ); \
         break; \
      }; \
 \
      jj = null; do \
      { \
         ( rtp->cs[ii][jj] ) = ZERO; \
      } while (( ++jj ) < THREE ); \
 \
      kk = null; do \
      { \
         jj = null; do \
         { \
            ll = (( kk+ONE )%FOUR ); \
            ( trp->v[null][jj] ) = edge[kk][jj]; \
            ( trp->v[ONE][jj]  ) = edge[ll][jj]; \
         } while (( ++jj ) < THREE ); \
 \
         ( trp->opt ) = 'c'; /* option: 'c'ross product */\
/*..........................................................................*/ \
         trp = triads( trp );   /* call triads(*) in option 'c'ross_product */ \
/*............................*/ \
         jj = null; do \
         { \
            ( rtp->cs[ii][jj] ) += ( trp->v[TWO][jj] ); \
         } while (( ++jj ) < THREE ); \
      } while (( ++kk ) < FOUR ); \
 \
      jj = null; do \
      { \
         ( rtp->cs[ii][jj] ) *= QUART; \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < THREE ); \
} /* end of macro CLL_CLSKEW */
/*----------------------------------------------------------------------------*/
/* The cell skew volume is defined as (1/3)*trace(B~*S), where B = ( b[i][j] )*/
/* and S = ( cs[i][j] ) are the node and skew vektor matrices, respectively, */
/* of the mesh cell: */

# define CLL_SKWVOL( ) \
{ \
   ( rtp->skv ) = ZERO; \
   ii = null; do \
   { \
      jj = null; do \
      { \
         ( rtp->skv ) += ( rtp->b[ii][jj] )*( rtp->cs[ii][jj] ); \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < THREE ); \
   ( rtp->skv ) /= THREE; \
} /* end of macro CLL_SKVOL */
/*----------------------------------------------------------------------------*/
/* There are different ways to compute the mesh cell volume, viz. by formulae */
/* (case NN=0)      vol = det(B) + (1/3) trace (B~*S)                         */
/* (case NN=1)      vol = (1/3) trace(B~*FM )                                 */
/* (case NN=2)      vol = (1/6) sum{j=1,2,3} ( b[j][]*( f[2j][]+f[2j+1][] ))  */
/* with the matrices B = (b[i][j]), S = ( cs[i][j] ), and F = ( fa[i][j] )    */

# define CLL_VOLUME(NN) \
if ((NN) == null ) \
{ \
   ii = null; do \
   { \
      jj = null; do \
      { \
         ( trp->v[ii][jj] ) = ( rtp->b[ii][jj] ); \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < THREE ); \
 \
   ( trp->opt ) = 'd'; /* option: 'd'eterminant */\
/*..........................................................................*/ \
   trp = triads( trp );         /* call triads(*) in option 'd'eterminant */ \
/*............................*/ \
 \
   ( rtp->vol ) = ( trp->det ) + ( rtp->skv ); \
} \
else if ((NN) == 1 ) \
{ \
   ( rtp->vol ) = ZERO; \
   ii = null; do \
   { \
      jj = null; do \
      { \
         ( rtp->vol ) += ( rtp->b[ii][jj] )*( rtp->fa[ii][jj] ); \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < THREE ); \
   ( rtp->vol ) /= THREE; \
} \
else /* NN == 2 */ \
{ \
   ( rtp->vol ) = ZERO; \
   jj = null; do \
   { \
      ( rtp->vp[jj] ) = ZERO; \
      ii = null; do \
      { \
         ( rtp->vp[jj] ) += ( rtp->b[( jj/2 )][ii] )*( rtp->f[jj][ii] ); \
      } while (( ++ii ) < THREE ); \
      ( rtp->vp[jj] ) /= 6.; \
      ( rtp->vol ) += ( rtp->vp[jj] ); \
   } while (( ++jj ) < FACES ); \
} /* end of macro CLL_VOLUME */
/*----------------------------------------------------------------------------*/
static GAUSS_JRD gss = { null };
/*
static JACOBI_EV jac = { null };
*/
/*============================================================================*/

CSHAPE *
cshape( CSHAPE *crp )
{
/* allusions: */
/*
   extern GAUSS_JRD gss;
   extern JACOBI_EV jac;
*/
/*----------------------------------------------------------------------------*/
/* declarations: */

   static CSHAPE
      cmx = {null},
     *rtp = &cmx;

   static GAUSS_JRD *gjp = &gss;
/*
   static JACOBI_EV *jev = &jac;

   static COMPLEX 
      cc = { null },
     *cpt = &cc;
*/
   static TRIADS 
      trd = { null },
     *trp = &trd;
/*
      ind = null,
*/
   static double
      xx = ZERO,
      yy = ZERO,
      zz = ZERO;

   static char
      fvx[FACES][FOUR] = {{ null }},
      pvx[PORTS][FOUR] = {{ null }};

# if CLL_DEBUG != 0
   static char
      ptr[STS_SIZE] = {null};
# endif

   static short 
      ii = null, 
      jj = null,
      kk = null,
      ll = null,
      mm = null,
      nn = null,
      pp = null,
      qq = null;

   static long
      vtx = null;

   static double
      edge[FOUR][DIMNS] = {{ZERO}}; /* edge vectors */

# if CLL_DEBUG == 3
   static double
      mtx[DIMNS][DIMNS] = {{ZERO}};
# endif
/*----------------------------------------------------------------------------*/
/* math. function prototypes: */

   double  
      fabs( double x ),
      sqrt( double x ),
      cos( double x ),
      sin( double x ),
     atan( double x ),
    atan2( double x, double y );

   int 
      abs( int h );

# if CLL_FULPIVT == 1
   GAUSS_JRD 
      *gssjpv( GAUSS_JRD *gjp );
# else
   GAUSS_JRD 
      *gssjrd( GAUSS_JRD *gjp );
# endif

   TRIADS
      *triads( TRIADS *trp );

   JACOBI_EV 
      *jacobi( JACOBI_EV *jcp );

   COMPLEX 
      *expc( COMPLEX *p );
/*----------------------------------------------------------------------------*/
/* initialize/reset: */

   if ( crp == NULL )
   {
      rtp = &cmx;

# if CLL_DEBUG != 0
      ptr[null] = null;
# endif

      ( rtp->rtn ) = null;
      ( rtp->opt ) = null;

      ( rtp->cell ) = null;
      ( rtp->face ) = null;
      ( rtp->port ) = null;

      ( rtp->skew ) = null;

      ( rtp->xn ) = ZERO;
      ( rtp->yn ) = ZERO;
      ( rtp->zn ) = ZERO;

      ( rtp->vol ) = ZERO;
      ( rtp->skv ) = ZERO;

      jj = null; do
      {
         ii = null; do /* clear vertex points */
         {
            ( rtp->c[ii][jj] ) = ZERO;
         } while (( ++ii ) < CRNRS );

         ii = null; do /* clear port vectors */
         {
            ( rtp->p[ii][jj] ) = ZERO;
            ( rtp->e[ii][jj] ) = ZERO;
         } while (( ++ii ) < PORTS );

         ii = null; do /* clear face vectors */
         {
            ( rtp->f[ii][jj] ) = ZERO;
         } while (( ++ii ) < FACES );

         ii = null; do
         {
            ( rtp->fa[ii][jj] ) = ZERO;

            ( rtp->a[ii][jj] ) = ZERO;
            ( rtp->ai[ii][jj] ) = ZERO;
            ( rtp->au[ii][jj] ) = ZERO;

            ( rtp->b[ii][jj] ) = ZERO;
            ( rtp->bi[ii][jj] ) = ZERO;
            ( rtp->bu[ii][jj] ) = ZERO;
            ( rtp->ub[ii][jj] ) = ZERO;

            ( rtp->cs[ii][jj] ) = ZERO;
            ( rtp->cu[ii][jj] ) = ZERO;

         } while (( ++ii ) < DIMNS );
      } while (( ++jj ) < DIMNS );

      ( rtp->rtn ) = ONE; /* normal return for structure initialization */
      return rtp;
   };
/*............................................................................*/
   rtp = crp;
/*............................................................................*/
/* now, do the requested job [ option ( crp->opt ) ]: */

   if ( null <= ( crp->opt ))
   {
/* overtake cell vertex coordinates [ for any given cell index > 0 ] */
      if ( null < ( crp->cell ))
      {
         ii = null; do
         {
            vtx = ( rtp->tpt->cm[( rtp->cell )][ii] );
            jj = null; do
            {
               ( rtp->c[ii][jj] ) = ( rtp->ppt->cpt->c[vtx][jj] );
            } while (( ++jj ) < THREE );
         } while (( ++ii ) < CRNRS );
      };
/*............................................................................*/
/* face vertices: */

      fvx[0][null] = 0;
      fvx[0][ONE] = 2;
      fvx[0][TWO] = 6;
      fvx[0][THREE] = 4;

      fvx[1][null] = 1;
      fvx[1][ONE] = 3;
      fvx[1][TWO] = 7;
      fvx[1][THREE] = 5;

      fvx[2][null] = 0;
      fvx[2][ONE] = 4;
      fvx[2][TWO] = 5;
      fvx[2][THREE] = 1;

      fvx[3][null] = 2;
      fvx[3][ONE] = 6;
      fvx[3][TWO] = 7;
      fvx[3][THREE] = 3;

      fvx[4][null] = 0;
      fvx[4][ONE] = 1;
      fvx[4][TWO] = 3;
      fvx[4][THREE] = 2;

      fvx[5][null] = 4;
      fvx[5][ONE] = 5;
      fvx[5][TWO] = 7;
      fvx[5][THREE] = 6;
/*............................................................................*/
/* port vertices: */

      pvx[0][null] = 3;
      pvx[0][ONE] = 7;
      pvx[0][TWO] = 2;
      pvx[0][THREE] = 6;

      pvx[1][null] = 5;
      pvx[1][ONE] = 7;
      pvx[1][TWO] = 4;
      pvx[1][THREE] = 6;

      pvx[2][null] = 1;
      pvx[2][ONE] = 5;
      pvx[2][TWO] = 0;
      pvx[2][THREE] = 4;

      pvx[3][null] = 1;
      pvx[3][ONE] = 3;
      pvx[3][TWO] = 0;
      pvx[3][THREE] = 2;

      pvx[4][null] = 6;
      pvx[4][ONE] = 7;
      pvx[4][TWO] = 4;
      pvx[4][THREE] = 5;

      pvx[5][null] = 3;
      pvx[5][ONE] = 7;
      pvx[5][TWO] = 1;
      pvx[5][THREE] = 5;

      pvx[6][null] = 2;
      pvx[6][ONE] = 3;
      pvx[6][TWO] = 0;
      pvx[6][THREE] = 1;

      pvx[7][null] = 2;
      pvx[7][ONE] = 6;
      pvx[7][TWO] = 0;
      pvx[7][THREE] = 4;

      pvx[8][null] = 5;
      pvx[8][ONE] = 7;
      pvx[8][TWO] = 1;
      pvx[8][THREE] = 3;

      pvx[9][null] = 6;
      pvx[9][ONE] = 7;
      pvx[9][TWO] = 2;
      pvx[9][THREE] = 3;

      pvx[10][null] = 4;
      pvx[10][ONE] = 6;
      pvx[10][TWO] = 0;
      pvx[10][THREE] = 2;

      pvx[11][null] = 4;
      pvx[11][ONE] = 5;
      pvx[11][TWO] = 0;
      pvx[11][THREE] = 1;
   };

   if ( null <= ( crp->opt ))
   {
/*............................................................................*/
      ( rtp->skew ) = null;
/*............................................................................*/
/* Mesh cell center ['node'] determination:

   Given 8 vertex points rtp->c[i][] ( i=0,...,7 ) nodal coordinates are
   computed as geometrical mesh cell center coordinates:
*/
      CLL_CENTER( );

/*............................................................................*/
/* Port vector determination:

   Given 8 vertex points rtp->c[i][] ( i=0,...,7 ) port vectors are
   computed in to the 'parcel twines' scheme ( viz. the port vectors
   interconnect the midpoints of the mesh edges , cf. [1] ).
*/
      CLL_PORTS( );

/*............................................................................*/
/* Compute cell edges: */

      CLL_EDGES( );

/*............................................................................*/
/* Compute mesh cell skew vectors */

      CLL_CLSKEW( );

/*............................................................................*/
/* Compute face vectors [ from cell port or edge vectors ]: */

      CLL_CLFCES( );

/*............................................................................*/
/* Compute opposite face vector [ arithmetic ] means: */

      ii = null; do
      {
         jj = null; do
         {
            ( rtp->fa[jj][ii] ) = \
               .5*(( rtp->f[2*jj][ii] ) + ( rtp->f[2*jj+ONE][ii]));
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* Area and node vector determination:

   Area and node vectors are computed as port vector cycles and sums
   of port vectors, respectively. Node and area [ line- ] vectors
   ( with respect to the canonical basis of euklidean 3-space, e.g.)
   constitute matrices A and B, respectively, such that  A *(B^-1) is 
   selfadjoint ( This follows from the 'parcel twines' scheme, cf.[1].)
   For every positive oriented node vector basis the matrix  A *(B^-1)
   then is positive definite. 
   Positive selfadjointness of A *(B^-1) will be checked in the following.
*/
/*............................................................................*/
/* Compute area vectors [ as port vector cycles ]: */

      ii = null; do
      {
         pp = (( ii+ONE ) % THREE );
         qq = (( ii+TWO ) % THREE );

         jj = null; do
         {                                       
            ( rtp->a[jj][ii] ) = ZERO;

            kk = ( 4*( jj + ONE )) % 12;
            ll = ( 4*( jj + TWO )) % 12;

            mm = null; do
            {
               nn = ONE; do
               {
                  ( rtp->a[jj][ii] ) += \
                   (( rtp->p[kk+mm][pp] )*( rtp->p[ll+nn][qq] ) - \
                   ( rtp->p[kk+mm][qq] )*( rtp->p[ll+nn][pp] ));
                  nn += TWO;
               } while ( nn < FOUR );

               mm += TWO;
            } while ( mm < THREE );
            ( rtp->a[jj][ii] ) *= QUART; /* j-th area vector: */
                                         /* ( rtp->a[j][*] ) */
                                         /* A := ( rtp->a[*][*] ) */
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* Compute node vectors: */

      ii = null; do
      {
         jj = null; do
         {                                       
            ( rtp->b[jj][ii] ) = ZERO;
            kk = FOUR*jj;

            ll = null; do
            {             
               ( rtp->b[jj][ii] ) += ( rtp->p[kk+ll][ii] );
                                           /* j-th node vector: */
            } while (( ++ll ) < FOUR );    /* ( rtp->b[j][*] ) */
            ( rtp->b[jj][ii] ) *= QUART;   /* B := ( rtp->b[*][*] ) */
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* Orthormalize node vectors b[j], yields intrinsic cell [ON] basis ub[]: */

      CLL_ORTNORM( );
/*............................................................................*/
/* Transform area, face, port, edge, and skew vectors into cell coordinates   */
/* [ i.e. with respect to basis ub[] ]: */

      ii = null; do
      {
         jj = null; do
         {
            ( rtp->pu[ii][jj] ) = ZERO;
            ( rtp->eu[ii][jj] ) = ZERO;
            ( rtp->cu[ii][jj] ) = ZERO;
            ( rtp->au[ii][jj] ) = ZERO;
            ( rtp->fu[ii][jj] ) = ZERO;

            kk = null; do
            {
               ( rtp->pu[ii][jj] ) += ( rtp->p[ii][kk] )*( rtp->ub[jj][kk] );
               ( rtp->eu[ii][jj] ) += ( rtp->e[ii][kk] )*( rtp->ub[jj][kk] );
               ( rtp->cu[ii][jj] ) += ( rtp->cs[ii][kk] )*( rtp->ub[jj][kk] );
               ( rtp->au[ii][jj] ) += ( rtp->a[ii][kk] )*( rtp->ub[jj][kk] );
               ( rtp->fu[ii][jj] ) += ( rtp->f[ii][kk] )*( rtp->ub[jj][kk] );
            } while (( ++kk ) < THREE );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );
      do
      {
         jj = null; do
         {
            ( rtp->pu[ii][jj] ) = ZERO;
            ( rtp->eu[ii][jj] ) = ZERO;
            ( rtp->fu[ii][jj] ) = ZERO;

            kk = null; do
            {
               ( rtp->pu[ii][jj] ) += ( rtp->p[ii][kk] )*( rtp->ub[jj][kk] );
               ( rtp->eu[ii][jj] ) += ( rtp->e[ii][kk] )*( rtp->ub[jj][kk] );
               ( rtp->fu[ii][jj] ) += ( rtp->f[ii][kk] )*( rtp->ub[jj][kk] );
            } while (( ++kk ) < THREE );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < FACES );
      do
      {
         jj = null; do
         {
            ( rtp->pu[ii][jj] ) = ZERO;
            ( rtp->eu[ii][jj] ) = ZERO;

            kk = null; do
            {
               ( rtp->pu[ii][jj] ) += ( rtp->p[ii][kk] )*( rtp->ub[jj][kk] );
               ( rtp->eu[ii][jj] ) += ( rtp->e[ii][kk] )*( rtp->ub[jj][kk] );
            } while (( ++kk ) < THREE );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < PORTS );
/*............................................................................*/
/* Compute mesh cell skew: skv = (1/3)*trace(B~*S). If ( skv != ZERO ) then   */
/* mark for skew in the cell ( skew := ONE ) */

      CLL_SKWVOL( ); /* skwvol := (1/3)*trace(B~*S) */

/*............................................................................*/
/* Compute mesh cell volume */
/* [ N=0: vol = det(B) + (1/3)*trace(B~*S) */
/*   N=1: vol = (1/3)*trace(B~*F) ]: */

      CLL_VOLUME(2);

      if(( fabs(( rtp->skv )/( rtp->vol ))) < ( 33.*CLL_PRECISION ))
         ( rtp->skew ) = null;
      else
         ( rtp->skew ) = ONE;
/*............................................................................*/
/* Debugging section [ display geometric cell characteristics ] */
/*............................................................................*/
# if CLL_DEBUG == 1
/* Display cell volume */

      if( *ptr != 48 ) /* dec 48 = ASCII char '0' */
      {
         printf( "\n det(B)=% .5e", ( trp->det ));
         printf( "\n vol(B)=% .5e", ( rtp->vol ));

         printf( "\n\n please acknowledge "
            "[ enter any character / Escape: enter 0 ]:" );
         scanf( "%s", ptr );

         if( *ptr == 48 )
            printf( "\n " );
      };
/*............................................................................*/
# elif CLL_DEBUG == 2
/* Display area vectors, opposite face vector means and skew vectors: */

      if( *ptr != 48 ) /* dec 48 = ASCII char '0' */
      {
         if( null != ( rtp->skew ))
         {
            printf( "\n\n Warning from function %s :", __func__ );
            printf( "\n Geometric skew detected !!!" );
            printf( "\n [ Area vectors differ from opposite "
               "face vector means.]" );

            printf( "\n\n Area vectors [ port vector cycles ] "
               "- the k-th vector is a[k][*]:" );

            ii = null; do
            {
               printf( "\n |" );
               jj = null; do
               {
                  printf( " a[%2d][%2d]=% .5e|", jj, ii, ( rtp->a[jj][ii] ));
               } while (( ++jj ) < THREE );
            } while (( ++ii ) < THREE );

            printf( "\n\n Opposite face vector means "
               "- the k-th vector is fa[k][*]:" );

            ii = null; do
            {
               printf( "\n |" );
               jj = null; do
               {
                  printf( "fa[%2d][%2d]=% .5e|", jj, ii, ( rtp->fa[jj][ii] ));
               } while (( ++jj ) < THREE );
            } while (( ++ii ) < THREE );

            printf( "\n\n Cell geometry skew vektors "
               "- the k-th vector is cs[k][*]:" );

            ii = null; do
            {
               printf( "\n |" );
               jj = null; do
               {
                  printf( "cs[%2d][%2d]=% .5e|", jj, ii, ( rtp->cs[jj][ii] ));
               } while (( ++jj ) < THREE );
            } while (( ++ii ) < THREE );

            printf( "\n\n please acknowledge "
               "[ enter any character / Escape: enter 0 ]:" );
            scanf( "%s", ptr );

            if( *ptr == 48 )
               printf( "\n " );
         };
      }; /* end if ptr[null] == null */
/*............................................................................*/
# elif CLL_DEBUG == 3
/* Display node and area vectors in exterior, and cell intrinsic coordinates */

      if( *ptr != 48 ) /* dec 48 = ASCII char '0' */
      {
         printf( "\n\n Node vectors "
            "- the k-th vector is b[k][*]:" );

         ii = null; do
         {
            printf( "\n |" );
            jj = null; do
	    {
	       printf( " b[%2d][%2d]=% .5e|", jj, ii, ( rtp->b[jj][ii] ));
	    } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         printf( "\n\n Orthonormalized node vectors [ node vector basis ]"
	    "\n - the k-th vector is ub[k][*]:" );

         ii = null; do
         {
	    printf( "\n |" );
	    jj = null; do
	    {
	       printf( "ub[%2d][%2d]=% .5e| ", jj, ii, ( rtp->ub[jj][ii] ));
	    } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         printf( "\n\n Node vectors in cell coordinates "
	    "[ relative to cell basis ub[] ]:" );

         ii = null; do
         {
	    printf( "\n |" );
	    jj = null; do
	    {
	       printf( "bu[%2d][%2d]=% .5e|", jj, ii, ( rtp->bu[jj][ii] ));
	    } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         printf( "\n\n The product of the former two "
	    "[ yields the node vectors, again ]:" );

         ii = null; do
         {
	    jj = null; do
	    {
	       ( rtp->mtx[ii][jj] ) = ZERO;
	       kk = null; do
	       { 
	          ( rtp->mtx[ii][jj] ) += \
                      ( rtp->bu[ii][kk] )*( rtp->ub[kk][jj] );
	       } while (( ++kk ) < THREE );
	    } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         ii = null; do
         {
	    printf( "\n |" );
	    jj = null; do
	    {
	       printf( "mtx[%2d][%2d]=% .5e|", jj, ii, ( rtp->mtx[jj][ii] ));
	    } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         printf( "\n\n Area vectors [ port vector cycles ]"
	    "\n - the k-th vector is a[k][*]:" );

         ii = null; do
         {
            printf( "\n |" );
	    jj = null; do
	    {
	       printf( " a[%2d][%2d]=% .5e|", jj, ii, ( rtp->a[jj][ii] ));
	    } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         printf( "\n\n Area vectors in cell coordinates "
            "[ with respect to cell basis ub[] ]:" );

         ii = null; do
         {
            printf( "\n |" );
            jj = null; do
            {
               printf( "au[%2d][%2d]=% .5e|", jj, ii, ( rtp->au[jj][ii] ));
            } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         printf( "\n\n please acknowledge "
            "[ enter any character / Escape: enter 0 ]:" );
         scanf( "%s", ptr );

         if( *ptr == 48 )
            printf( "\n " );

      }; /* end if ( *ptr != 48 ) */

# endif /* CLL_DEBUG == 3 */
/*............................................................................*/
/* compute (A^-1): */

      ii = null; do        
      {
         jj = null; do
         {
            ( gjp->mr[ii][jj] ) = ( rtp->a[ii][jj] );
            ( gjp->mi[ii][jj] ) = ZERO;
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: A |-> (A^-1) */

      ( gjp->rank ) = DIMNS;
      ( gjp->neqs ) = DIMNS;
      ( gjp->opt ) = 'i';

# if CLL_FULPIVT == 1
      gjp = gssjpv( gjp );
# else
      gjp = gssjrd( gjp );
# endif
/*............................................................................*/
      if ( gjp == NULL )
      { 
         printf( "\n\n Message from function %s :", __func__ );
         printf( "\n Error on area vector matrix "
            "inversion A |-> (A^-1) !!!\n " );

         ( rtp->rtn ) = null;
         return rtp;
      };
/*............................................................................*/
/* copy (A^-1), then compute (B^-1): */

      ii = null; do           
      {
         jj = null; do
         {
            ( rtp->ai[ii][jj] ) = ( gjp->zr[ii][jj] ); /* copy (A^-1) */
            ( gjp->mr[ii][jj] ) = ( rtp->b[ii][jj] );
            ( gjp->mi[ii][jj] ) = ZERO;
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: B |-> (B^-1): */

      ( gjp->rank ) = DIMNS;
      ( gjp->neqs ) = DIMNS;
      ( gjp->opt ) = 'i';

# if CLL_FULPIVT == 1
      gjp = gssjpv( gjp );
# else
      gjp = gssjrd( gjp );
# endif
/*............................................................................*/
      if ( gjp == NULL )
      {
         printf( "\n\n Message from function %s :", __func__ );
         printf( "\n Error on node vector matrix "
            "inversion B |-> (B^-1) !!!\n " );

         ( rtp->rtn ) = null;
         return rtp;
      };

      if (( gjp->dtr ) < 1.e-277 )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Degenerate cell !!!\n " );
         printf( "\n [ Non-positive node vector orientation.]\n " );

         ( rtp->rtn ) = null;
         return rtp;
      }
      else if (( gjp->dtr ) <= ZERO )
      {
         printf( "\n\n Error message from function %s :", __func__ );
         printf( "\n Degenerate cell !!!\n " );
         printf( "\n [ Non-positive volume.]\n " );

         ( rtp->rtn ) = null;
         return rtp;
      };
/*............................................................................*/
/* copy (B^-1) into CSHAPE *rtp: */

      ii = null; do
      {
         jj = null; do
         {
            ( rtp->bi[ii][jj] ) = ( gjp->zr[ii][jj] ); /* copy (B^-1) */
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
   }; /* end if ( crp->opt ) == ... */
/*............................................................................*/
   if ((( crp->opt ) == THREE )
     ||(( crp->opt ) == 'f' )
     ||(( crp->opt ) == 'F' ))
   {
/*............................................................................*/
/* cell face parameters: */
      kk = ( crp->face );

      jj = null; do
      {
         switch( kk )
         {
           case 0:
            ( trp->v[null][jj] ) = ( rtp->p[7][jj] );
            ( trp->v[ONE][jj] ) = ( rtp->p[10][jj] );
            break;

           case 1:
            ( trp->v[null][jj] ) = ( rtp->p[5][jj] );
            ( trp->v[ONE][jj] ) = ( rtp->p[8][jj] );
            break;

           case 2:
            ( trp->v[null][jj] ) = ( rtp->p[11][jj] );
            ( trp->v[ONE][jj] ) = ( rtp->p[2][jj] );
            break;

           case 3:
            ( trp->v[null][jj] ) = ( rtp->p[9][jj] );
            ( trp->v[ONE][jj] ) = ( rtp->p[0][jj] );
            break;

           case 4:
            ( trp->v[null][jj] ) = ( rtp->p[3][jj] );
            ( trp->v[ONE][jj] ) = ( rtp->p[6][jj] );
            break;

           case 5:
            ( trp->v[null][jj] ) = ( rtp->p[1][jj] );
            ( trp->v[ONE][jj] ) = ( rtp->p[4][jj] );
            break;
         };

         ( trp->v[TWO][jj] ) = ( rtp->f[kk][jj] );

      } while (( ++jj ) < THREE );

      ( trp->opt ) = 'o'; /* option: 'o'rthonormalize */
/*............................................................................*/
      trp = triads( trp );     /* call triads(*), option 'o'rthonormlze       */
/*...........................*/

      if ( 1.E-277 < fabs( trp->det ))
      {
         ii = null; do
         {
            jj = null; do
            {                                             /* uv[] is a [local */
               ( rtp->uv[ii][jj] ) = ( trp->uv[ii][jj] ); /* face] ON basis   */
               ( rtp->vu[ii][jj] ) = ( trp->vu[ii][jj] ); /* vu[] is the new  */
                                         /* port vector matrix V [composed of */
            } while (( ++jj ) < THREE ); /* the columns of port vector coor-  */
         } while (( ++ii ) < THREE );    /* dinates relative to uv[] ]        */
/*............................................................................*/
/* compute adj(V^-1); first copy adj(V)[jj][ii] = vu[ii][jj] into struct *gjp:*/

         ii = null; do
         {
            jj = null; do
            {
               ( gjp->mr[jj][ii] ) = ( rtp->vu[ii][jj] );
               ( gjp->mi[jj][ii] ) = ZERO;
            } while (( ++jj ) < THREE ); /* sufficient: TWO */
         } while (( ++ii ) < THREE ); /* sufficient: TWO */
/*............................................................................*/
/* ... then perform matrix inversion; adj(V) |-> adj(V^-1) = ( adj(V) )^-1: */

         ( gjp->rank ) = THREE; /* sufficient: TWO */
         ( gjp->neqs ) = THREE; /* sufficient: TWO */
         ( gjp->opt ) = 'i'; /* [ matrix 'i'nversion ] */
/*............................................................................*/
# if CLL_FULPIVT == 1
         gjp = gssjpv( gjp );
# else
         gjp = gssjrd( gjp );
# endif
/*............................................................................*/
         if ( gjp == NULL )
         { 
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n Error on port vector matrix inversion "
            "(V*) |-> ((V*)^-1) !!!\n " );

            ( rtp->rtn ) = ONE;
            return rtp;
         };
/*............................................................................*/
/* ... finally copy adj(V^-1) as the [thus redefined] matrix vu[ii][jj] := */ 
/* := adj(V^-1)[ii][jj] */

         ii = null; do
         {
            jj = null; do
            {
               ( rtp->vu[ii][jj] ) = ( gjp->zr[ii][jj] ); /* copy ((V*)^-1) */
            } while (( ++jj ) < THREE ); /* sufficient: TWO */
         } while (( ++ii ) < THREE ); /* sufficient: TWO */
/*............................................................................*/
      }
      else /* if |( trp->det )| =~ ZERO [case of degenerate face] */
      {
         ii = null; do           
         {
            jj = null; do
            {
               ( rtp->vu[ii][jj] ) = ZERO;
            } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );
      };
   };

   ( rtp->rtn ) = null;
   return rtp;
} 
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/*----------------------------------------------------------------------------*/
# ifdef CLL_CLFCES
   # undef CLL_CLFCES
# endif
# ifdef CLL_CLSKEW
   # undef CLL_CLSKEW
# endif
# ifdef CLL_SKWVOL
   # undef CLL_SKWVOL
# endif
# undef CLL_CENTER
# undef CLL_PORTS
# undef CLL_EDGES
# undef CLL_ORTNORM
# undef CLL_VOLUME
/*----------------------------------------------------------------------------*/
# undef CLL_DSPINTM
# undef CLL_FULPIVT 
# undef CLL_FCEMODE
# undef CLL_PRECISION
# undef CLL_TRIVIAL
# undef CLL_INCLUDE
# undef CLL_DEBUG
/*----------------------------------------------------------------------------*/
# undef CPRODUCT
# undef EPS_VAC
# undef MY_VAC_
# undef ELECTR_CHRGE
# undef ELECTR_MASS_
# undef BOHRs_MAGNTN
# undef PLANCKs_CNST
/************************ end of function 'cshape(*)' *************************/
