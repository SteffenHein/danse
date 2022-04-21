/* [ file: hcrsmx.c ] */
/*******************************************************************************
*                                                                              *
*   ISO/ANSI C function hcrsmx(*)                                              *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   DSC heat current s-matrix generation function supporting non-orthogonal    *
*   cell.                                                                      *
*                                                                              *
*   Given eigth points in 3-space, hsp->c[i][] ( i=0,...,7 ), a DSC time       *
*   step, hsp->dt [seconds], and material parameters, hsp->kh, hsp->cv,        *
*   which respectively denote the heat current conductivity [W/(K*m)]          *
*   and [volume specific] heat capacity [Joule/(K*m^3)] in the cell,           *
*   this function returns the nodal heat and fluid propagation parameters,     *
*   the form vectors                                                           *
*             s  =  ( hsp->s[j][k] ) j,k = 0,...,5   ,                         *
*                                                                              *
*   as well as the area and node vector matrices adj(A) = ( hsp->a[][] )       *
*   and adj(B) = ( hsp->b[][] ) [with their inverses] of the mesh cell         *
*   cell with vertex points hsp->c[i][] [ opt = 'r'], needed by the DSC        *
*   heat and fluid propagation algorithm.                                      *
*                                                                              *
*   Also [ in option opt = 't' ], given the vertex points hsp->c[i][],         *
*   and material parameters, an approximate stability upper bound for          *
*   the DSC time step is is returned as hsp->dt .                              *
*                                                                              *
*   All parameters are transferred with structure cmx [ of type  HCRSMX ]      *
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
# if defined ( OPTIMIZE )
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize this header for FORMER.C configuration: */
# include "../former/FORMER.CONF" /* Edit and customize this configuration  */
                                  /* header for the preprocessing routine */
/*----------------------------------------------------------------------------*/
# define HMX_INCLUDE 1 /* 1: include heat current s-parameters transfer struct*/
# define HMX_DEBUG   0 /*>0: activate various printing functions [ for debug- */
                       /*    ging purposes, mainly ] */
/*----------------------------------------------------------------------------*/
/* operation marks: */

# define HMX_FULPIVT 1     /* 1: fully [line and column] pivoted Gauss-Jordan */
                           /* elimination */

# define HMX_FCEMODE 1     /* 1 or 2: */
                           /* 1: compute face vectors from cell ports */
                           /* 2: compute face vectors from cell edges */
# if DSC_FLDMDE != 0
   # define HMX_CLLCRDS 0  /* 0: use global coordinates */
# else
   # define HMX_CLLCRDS 1  /* 1: use local [ cell internal ] coordinates */
# endif                    /* [ repeater may identify more equivalent cells */
                           /* -- this option is illeagal with unisotropic */
                           /* media, such as fluids, e.g. ] */
# define HMX_DSPINTM 0     /* HMX_DSPINTM 1: display intermediate results   */
/*----------------------------------------------------------------------------*/
/* computational parameters, regularization factors, bounds etc.: */

# ifdef PRECISION
   # define HMX_PRECISION PRECISION
# else
   # define HMX_PRECISION ( 1.e-15 ) 
# endif
                                 /* lower bounds for non-trivial check: */
# define HMX_TRIVIAL ( 1.e-277 ) /* heat current conductivity (kh) */
# define HMX_INITSTP ( 1.e+277 ) /* the initial time step */
# define HMX_1_MINUS ( .999999 ) /* 'secure' 1 [ for process stability ] */
/*----------------------------------------------------------------------------*/
/* array dimensions: */
# ifndef CRNRS
   # define CRNRS  8
# endif
# ifndef SORDR
   # define SORDR  6
# endif
# ifndef PORTS
   # define PORTS 12
# endif
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
# define STEFAN_BOLTZ ( 5.6705100000e-08 ) /* Stefan-Boltzmann [W/(K^4*m^2)]  */
/*----------------------------------------------------------------------------*/
# include "../math/trdstp.h"
# include "../former/gssjtp.h"
# include "../former/jacbtp.h"
/*---------------------------------------------------------------------------*/
# if HMX_INCLUDE == 1
   # include "hcrstp.h"
# else
/*-----------------------------------------------------------------------------
*//* heat current s-matrix data transfer structure, release 6.1.
*//* All units are international units (mks), if not otherwise specified.
------------------------------------------------------------------------------*/
typedef struct
{
   signed char /* any return character [ may be used as an error flag, e.g. ] */
      rtn;

   char 
      opt; /* hcrsmx(*) options: 't' time_step, 'r': real s-parameters [ in */
           /* time-domain ], 'c': complex s-parameters [ in frequency-dmn ] */
   char
      orth, /* orthogonal cell indicator [ non-orth:0 / orth:1 ] */
      skew, /* geometric skew indicator [ no skew:0 / skew:1 ] */
      loss, /* loss indicator [ lossless cell:0 / lossy cell:1 ] */
      isotrop; /* media isotropy indicator [ unisotropic: 0 / isotropic: 1 ] */ 

   short 
      med; /* media label */
/*............................................................................*/
/* input parameters: */

   double
      dt,  /* time step [sec] */
      dp,  /* phase shift [rad] */
      ke,  /* electric conductivity [A/(V*m) = S/m] */
      km,  /* magnetic conductivity [V/(A*m) = Ohm/m] */
      kh,  /* heat current conductivity [W/(K*m)] */
      cv;  /* heat capacity [J/(K*m^3)] */
/*............................................................................*/
# if DSC_FLDMDE != 0
   double
      rm, /* mean mass density [Kg/(m^3)] */
      tm, /* mean temperature [C] */
      bm, /* mean expansion coefficient [1/K] */
      cm, /* mean [ adiabatic ] compression coefficient [1/Pa] */
      ny, /* dynamic viscosity [Kg/(sec*m)] */
      q1, /* Cp/Cv - 1 [dimensionless] */
      td, /* [dissipation] time constant [sec] */
      LL; /* characteristic length [in Prandtl turbulence model, e.g.] */

   double
      gr[THREE], /* gravitational acceleration [m/sec^2] */ 
      gp[THREE]; /* pressure gradient [Kg/((sec*m)^2)] */ 

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
 double
      c[CRNRS][THREE]; /* corner point coordinates */
/*............................................................................*/
/* computed geometric parameters: */

   double
      vol, /* cell volume [ m^3 ] */
      skv; /* cell skew volume [ m^3 ] */

   double
     fm[FACES],        /* size of face */
     vp[FACES];        /* vp[jj] pyramide volume ( F[jj], b[jj/2] )/6 */

   double
      e[PORTS][THREE], /* edge vectors */
     eu[PORTS][THREE], /* edge vectors relative to cell basis ub[] */

      p[PORTS][THREE], /* port vectors */
     pu[PORTS][THREE], /* port vectors relative to cell basis ub[] */

      f[FACES][THREE], /* face vectors */
     fu[FACES][THREE], /* face vectors relative to cell basis ub[] */
     fa[THREE][THREE], /* opposite face vector arithmetic means */

      a[THREE][THREE], /* area vector matrix A = arv[i][j] */
     au[THREE][THREE], /* area vector matrix relative to cell basis ub[] */
     ai[THREE][THREE], /* Inverse are vector matrix A^-1 */

      b[THREE][THREE], /* node vector matrix B = ndv[i][j] */
     bu[THREE][THREE], /* normalized Node vector matrix relative to ub[] */
     ub[THREE][THREE], /* Gram-Schmidt orthonormalzd node vectors [cell basis]*/
     bi[THREE][THREE], /* inverse node vector basis B^-1 */

     cs[THREE][THREE], /* cell shape skew vectors */
     cu[THREE][THREE]; /* cell shape skew vectors relative to cell basis */
/*............................................................................*/
/* computed [ temperature updating ] s-parameters [-> elf.smx<N>] */

   double
     hdt, /* heat diffusion time step [sec] */
      ct; /* heating coefficient [ dt/cv ] */
/*............................................................................*/
# if DSC_FLDMDE != 0
   double
     fdt, /* fluid [ viscous ] diffusion time step [sec] */
      ft; /* fluid updating coefficient [ fdt/rm ] */
# endif
/*............................................................................*/
/* face form vectors [ (F)*(B^-1) ] */
   double
     s[FACES][THREE];  /* face form vectors [ (F)*(B^-1) ] */

   char
     name[STS_SIZE],
     text[STS_SIZE], 
     ttyp[STS_SIZE]; 

} HCRSMX;
/*----------------------------------------------------------------------------*/
# endif /* HMX_INCLUDE != 1 */
/*----------------------------------------------------------------------------*/
/* system function prototypes: */

void abort( void );
int printf( const char *format,...);
int scanf( const char *format,...);
/*----------------------------------------------------------------------------*/
/* macros: */
/*----------------------------------------------------------------------------*/
/* product of complex numbers: UU+i*VV = (AA+i*BB)*(CC+i*DD): */
# define CPRODUCT( AA, BB, CC, DD, UU, VV ) \
{ \
   (UU) = (AA)*(CC) - (BB)*(DD); \
   (VV) = (BB)*(CC) + (AA)*(DD); \
}
/*............................................................................*/
/* quotient of complex numbers: UU+i*VV = (AA+i*BB)/(CC+i*DD): */
# define CQUOTIENT( AA, BB, CC, DD, UU, VV ) \
{ \
   (VV) = (CC)*(CC) + (DD)*(DD); \
   (UU) = ((AA)*(CC) + (BB)*(DD))/(VV); \
   (VV) = ((BB)*(CC) - (AA)*(DD))/(VV); \
}
/*----------------------------------------------------------------------------*/
/* compute mesh cell port vectors from vertex points: */
# define HMX_PORTS( ) \
{ \
   jj = null; do \
   { \
      ( rtp->p[0][jj] ) = .5*( rtp->c[3][jj] + rtp->c[7][jj] - \
                               rtp->c[2][jj] - rtp->c[6][jj] ); \
      ( rtp->p[1][jj] ) = .5*( rtp->c[5][jj] + rtp->c[7][jj] - \
                               rtp->c[4][jj] - rtp->c[6][jj] ); \
      ( rtp->p[2][jj] ) = .5*( rtp->c[1][jj] + rtp->c[5][jj] - \
                               rtp->c[0][jj] - rtp->c[4][jj] ); \
      ( rtp->p[3][jj] ) = .5*( rtp->c[1][jj] + rtp->c[3][jj] - \
                               rtp->c[0][jj] - rtp->c[2][jj] ); \
      ( rtp->p[4][jj] ) = .5*( rtp->c[7][jj] + rtp->c[6][jj] - \
                               rtp->c[4][jj] - rtp->c[5][jj] ); \
      ( rtp->p[5][jj] ) = .5*( rtp->c[3][jj] + rtp->c[7][jj] - \
                               rtp->c[1][jj] - rtp->c[5][jj] ); \
      ( rtp->p[6][jj] ) = .5*( rtp->c[3][jj] + rtp->c[2][jj] - \
                               rtp->c[0][jj] - rtp->c[1][jj] ); \
      ( rtp->p[7][jj] ) = .5*( rtp->c[2][jj] + rtp->c[6][jj] - \
                               rtp->c[0][jj] - rtp->c[4][jj] ); \
      ( rtp->p[8][jj] ) = .5*( rtp->c[5][jj] + rtp->c[7][jj] - \
                               rtp->c[1][jj] - rtp->c[3][jj] ); \
      ( rtp->p[9][jj] ) = .5*( rtp->c[7][jj] + rtp->c[6][jj] - \
                               rtp->c[3][jj] - rtp->c[2][jj] ); \
      ( rtp->p[10][jj] ) = .5*( rtp->c[4][jj] + rtp->c[6][jj] - \
                                rtp->c[0][jj] - rtp->c[2][jj] ); \
      ( rtp->p[11][jj] ) = .5*( rtp->c[4][jj] + rtp->c[5][jj] - \
                                rtp->c[0][jj] - rtp->c[1][jj] ); \
   } while (( ++jj ) < THREE ); \
}
/*----------------------------------------------------------------------------*/
/* compute mesh cell edge vectors from vertex points: */

# define HMX_EDGES( ) \
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
   } while (( ++jj ) < THREE ); \
}
/*----------------------------------------------------------------------------*/
/* orthonormalize node vectors b[] to a 'cell basis' ub[], and transform node */
/* vectors into cell basis coordinates; this yields vectors bu[]: */
# define HMX_ORTNORM( ) \
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
/* There are different ways to compute the cell face vectors [ FCEMODE 1,2 ]: */
# if HMX_FCEMODE == 1
# define HMX_CLFCES( ) \
{ \
   ii = null; do \
   { \
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
      xx = ZERO; \
      jj = null; do \
      { \
         yy = ( trp->v[TWO][jj] ); \
         ( rtp->f[ii][jj] ) = yy; \
         xx += ( yy*yy ); \
      } while (( ++jj ) < THREE ); \
      ( rtp->fm[ii] ) = sqrt( xx ); \
   } while (( ++ii ) < FACES ); \
} /* end of macro HMX_CLFCES */
/*----------------------------------------------------------------------------*/
# elif HMX_FCEMODE == 2
# define HMX_CLFCES( ) \
{ \
   ii = null; do \
   { \
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
      } while (( ++kk ) < FOUR ); \
 \
      xx = ZERO; \
      jj = null; do \
      { \
         ( rtp->f[ii][jj] ) *= QUART; \
         xx += (( rtp->f[ii][jj] )*( rtp->f[ii][jj] )); \
      } while (( ++jj ) < THREE ); \
      ( rtp->fm[ii] ) = sqrt( xx ); \
   } while (( ++ii ) < FACES ); \
} /* end of macro HMX_CLFCES */
# endif /* HMX_FCEMODE == 2 */
/*----------------------------------------------------------------------------*/
# define HMX_CLSKEW( ) \
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
} /* end of macro HMX_CLSKEW */
/*----------------------------------------------------------------------------*/
/* The cell skew volume is defined as (1/3)*trace(B~*S), where B = ( b[i][j] )*/
/* and S = ( cs[i][j] ) are the node and skew vektor matrices, respectively, */
/* of the mesh cell: */

# define HMX_SKWVOL( ) \
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
} /* end of macro HMX_SKVOL */
/*----------------------------------------------------------------------------*/
/* There are different ways to compute the mesh cell volume, viz. by formulae */
/* case NN=0:      vol = det(B) + (1/3) trace(B~*S)                           */
/* case NN=1:      vol = (1/3) trace(B~*FM)                                   */
/* case NN=2:      vol = (1/6) sum{j=1,2,3} ( b[j][]*( f[2j][]+f[2j+1][] ))   */
/* with the matrices B = ( b[i][j] ), S = ( cs[i][j] ), and F = ( fa[i][j] )  */

# define HMX_VOLUME(NN) \
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
 \
   ( rtp->vol ) /= THREE; \
 \
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
      ( rtp->vp[jj] ) /= SIX; \
      ( rtp->vol ) += ( rtp->vp[jj] ); \
   } while (( ++jj ) < FACES ); \
} /* end of macro HMX_VOLUME */
/*----------------------------------------------------------------------------*/
static GAUSS_JRD gss = { null };
/*
static JACOBI_EV jac = { null };
*/
/*============================================================================*/

HCRSMX *
hcrsmx( HCRSMX *hsp )
{
/* allusions: */
/*
   extern GAUSS_JRD gss;
   extern JACOBI_EV jac;
*/
/* declarations: */

   static HCRSMX
      cmx = {null},
     *rtp = &cmx;

   static GAUSS_JRD 
      *gjp = &gss;
/*
   static JACOBI_EV
      *jev = &jac;

   static COMPLEX 
      cc = { null },
     *cpt = &cc;
*/
   static TRIADS 
      trd = { null },
     *trp = &trd;

   static signed char
      opt = null;
/*
      ind = null,
*/
   static char
     *trv = "trivial_TH",
     *hcr = "hcurr";

# if HMX_DEBUG != 0
   static char
      ptr[STS_SIZE] = {null};
# endif
/*............................................................................*/
# if DSC_FLDMDE != 0
   static char
     *fld = "fluid";
# endif
/*............................................................................*/
   static short 
      ii = null, 
      jj = null,
      kk = null,
      ll = null,
      mm = null,
      nn = null,
      pp = null,
      qq = null;

   static double
      xx = ZERO,
      yy = ZERO,
      edge[FOUR][THREE] = {{ZERO}}; /* edge vectors */

# if HMX_DEBUG == 3
   static double
      mtx[THREE][THREE] = {{ZERO}};
# endif
/*
   struct emhcr
   {
      struct emhcr *ecp;

      double
         ec, mc,

         tue[THREE][THREE],
          ve[THREE][THREE],
          pe[THREE][THREE],
          qe[THREE][THREE],
         eue[THREE][THREE],
         tum[THREE][THREE],
          vm[THREE][THREE],
          pm[THREE][THREE],
          qm[THREE][THREE],
         eum[THREE][THREE];
   };
   static struct emhcr emc = {null};
*/
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

# if HMX_FULPIVT == 1
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
/* initialize, reset: */

   if ( hsp == NULL )
   {
      rtp = &cmx;

     init:

# if HMX_DEBUG != 0
      ptr[null] = null;
# endif

      ( rtp->rtn ) = null;
      ( rtp->opt ) = null;

      ( rtp->loss ) = null;
      ( rtp->skew ) = null;
      ( rtp->isotrop ) = null;

      ii = null; do
      {
         ( rtp->name[ii] ) = null;
         ( rtp->text[ii] ) = null;
      } while (( ++ii ) < STS_SIZE );

      ( rtp->hdt ) = ZERO;
      ( rtp->dt ) = HMX_INITSTP;
      ( rtp->ct ) = ZERO;
      ( rtp->kh ) = ZERO;
      ( rtp->ke ) = ZERO;
      ( rtp->km ) = ZERO;
      ( rtp->cv ) = ZERO;
      ( rtp->vol ) = ZERO;
      ( rtp->skv ) = ZERO;
/*............................................................................*/
# if DSC_FLDMDE != 0
      ( rtp->fdt ) = ZERO;
      ( rtp->ft ) = ZERO;
      ( rtp->rm ) = ZERO;
      ( rtp->tm ) = ZERO;
      ( rtp->bm ) = ZERO;
      ( rtp->cm ) = ZERO;
      ( rtp->ny ) = ZERO;
      ( rtp->q1 ) = ZERO;
# endif
/*............................................................................*/
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
         } while (( ++ii ) < THREE );
/*............................................................................*/
# if DSC_FLDMDE != 0
         ( rtp->gr[jj] ) = ZERO;
         ( rtp->gp[jj] ) = ZERO;
# endif
/*............................................................................*/
      } while (( ++jj ) < THREE );

      jj = null; do
      {
         ( rtp->fm[jj] ) = ZERO;
         ii = null; do
         {
            ( rtp->s[jj][ii] ) = ZERO;
         } while (( ++ii ) < THREE );
      } while (( ++jj ) < FACES );

      strncpy(( rtp->ttyp ), trv, STS_SIZE );

      ( rtp->rtn ) = ONE; /* normal return for structure initialization */
      return rtp;
   };
/*............................................................................*/
   rtp = hsp;

   switch(( hsp->opt ))
   {
     case 0:
      goto init;
      break;

     case 't':
     case 'T':
     case 'p':
     case 'P':

      opt = 't';
      break;

     case 's':
     case 'S':
     case 'r':
     case 'R':

      opt = 's';
      break;

     default:
      fprintf( stderr, "\n\n Error message from function %s:", __func__ );
      fprintf( stderr, "\n Illegal option '%c' in function call !!!", opt );
      fprintf( stderr, "\n Allowed options are [ case insensitive ]:" );
      fprintf( stderr, "\n opt = 't' - compute stability time step" );
      fprintf( stderr, "\n opt = 's' - compute s-parameters" );
      fprintf( stderr, "\n [ Check function call in calling program.]\n " );

      ( rtp->rtn ) = null;
      return rtp;
      break;
   };
/*............................................................................*/
   if (( rtp->med ) < TWO )
   {
      strcpy(( rtp->ttyp ), trv );
      goto init; /* reset parameters and */
                 /* return: trivial cell */
   }
   else if (( fabs( hsp->kh ) < HMX_TRIVIAL )
          &&( fabs( hsp->cv ) < HMX_TRIVIAL ))
   {
      strcpy(( rtp->ttyp ), trv );
      goto init; /* reset parameters and */
                 /* return: trivial cell */
   }
   else /* non trivial media parameters */
   {
      ( rtp->kh ) = fabs( hsp->kh );
      ( rtp->cv ) = fabs( hsp->cv );
/*............................................................................*/
# if DSC_FLDMDE != 0

      if (( fabs( hsp->bm ) < HMX_TRIVIAL )
        ||( fabs( hsp->ny ) < HMX_TRIVIAL ))
      {
         strcpy(( rtp->ttyp ), hcr );
         ( rtp->isotrop ) = ONE;

         ( hsp->rm ) = ZERO;
         ( hsp->tm ) = ZERO;
         ( hsp->bm ) = ZERO;
         ( hsp->cm ) = ZERO;
         ( hsp->ny ) = ZERO;
         ( hsp->q1 ) = ZERO;

         jj = null; do
         {
            ( hsp->gr[jj] ) = ZERO;
            ( hsp->gp[jj] ) = ZERO;
         } while (( ++jj ) < THREE );
      }
      else /* fluid type */
      {
         strcpy(( rtp->ttyp ), fld );
         ( rtp->isotrop ) = null;

         ( rtp->rm ) = fabs( hsp->rm );
         ( rtp->tm ) = fabs( hsp->tm );
         ( rtp->bm ) = fabs( hsp->bm );
         ( rtp->cm ) = fabs( hsp->cm );
         ( rtp->ny ) = fabs( hsp->ny );
         ( rtp->q1 ) = fabs( hsp->q1 );
      };
/*............................................................................*/
# else /* if DSC_FLDMDE == 0 */
      strcpy(( rtp->ttyp ), hcr );
      ( rtp->isotrop ) = ONE;
# endif
/*............................................................................*/
      if ( HMX_TRIVIAL <= fabs( hsp->ke ))
      {
         strcat(( rtp->ttyp ), "_E" );
         ( rtp->loss ) = ONE;

         ( rtp->ke ) = fabs( hsp->ke );

         if ( HMX_TRIVIAL <= fabs( hsp->km ))
         {
            strcat(( rtp->ttyp ), "&M" );
            ( rtp->km ) = fabs( hsp->km );
         }
	 else
            ( hsp->km ) = ZERO;
      }
      else if ( HMX_TRIVIAL <= fabs( hsp->km ))
      {
         strcat(( rtp->ttyp ), "_M" );
         ( rtp->loss ) = ONE;

         ( rtp->km ) = fabs( hsp->km );
         ( hsp->ke ) = ZERO;
      }
      else /* non lossy electric media */
      {
         strcat(( rtp->ttyp ), "__" );
         ( rtp->loss ) = null;

         ( hsp->ke ) = ZERO;
         ( hsp->km ) = ZERO;
      };
   }; /* non trivial media parameters */
/*............................................................................*/
/* Compute port vectors:

   Given 8 vertex points rtp->c[i][] ( i=0,...,7 ) port vectors are
   computed in to the 'parcel twines' scheme ( viz. the port vectors
   interconnect the midpoints of the mesh edges , cf. [1] ).
*/
   HMX_PORTS( );
/*............................................................................*/
/* Cell edges: */

   HMX_EDGES( );
/*............................................................................*/
/* Mesh cell skew: */

   HMX_CLSKEW( );
/*............................................................................*/
/* face vectors [ computed from cell port or edge vectors ]: */

   HMX_CLFCES( );
/*............................................................................*/
/* Opposite face vector [ arithmetic ] means: */

   ii = null; do
   {
      jj = null; do
      {
         ( rtp->fa[jj][ii] ) = \
            .5*(( rtp->f[2*jj][ii] ) + ( rtp->f[2*jj+ONE][ii]));
      } while (( ++jj ) < THREE );
   } while (( ++ii ) < THREE );
/*............................................................................*/
/* Area and node vectors:

   Area and node vectors are computed as port vector cycles and sums
   of port vectors, respectively. Node and area [ line- ] vectors
   ( with respect to the canonical basis of euklidean 3-space, e.g.)
   constitute matrices A and B, respectively, such that  A*(B^-1) is 
   selfadjoint ( This follows from the 'parcel twines' scheme, cf.[1].)
   For every positive oriented node vector basis the matrix A*(B^-1)
   is then positive definite. 
   Positive selfadjointness of A*(B^-1) will be checked in the following.
*/
/*............................................................................*/
/* Area vectors [ computed as port vector cycles ]: */

   ii = null; do
   {
      pp = (( ii + ONE )%THREE );
      qq = (( ii + TWO )%THREE );

      jj = null; do
      {                                       
         ( rtp->a[jj][ii] ) = ZERO;

         kk = (( 4*( jj + ONE )) % 12 );
         ll = (( 4*( jj + TWO )) % 12 );

         mm = null; do
         {
            nn = ONE; do
            {
               ( rtp->a[jj][ii] ) += \
                (( rtp->p[kk+mm][pp] )*( rtp->p[ll+nn][qq] ) - \
                ( rtp->p[kk+mm][qq] )*( rtp->p[ll+nn][pp] ));
               nn += TWO;
            } while( nn < FOUR );

            mm += TWO;
         } while ( mm < THREE );
         ( rtp->a[jj][ii] ) *= QUART; /* j-th area vector: */
      } while (( ++jj ) < THREE );    /* ( rtp->a[j][*] )  */
   } while (( ++ii ) < THREE );       /* ( A* ) = adj(A) := ( rtp->a[*][*] ) */
/*............................................................................*/
/* Node vectors: */

   ii = null; do
   {
      jj = null; do
      {                                       
         ( rtp->b[jj][ii] ) = ZERO;
         kk = FOUR*jj;

         ll = null; do
         {             
            ( rtp->b[jj][ii] ) += ( rtp->p[kk+ll][ii] );
         } while (( ++ll ) < FOUR );            
         ( rtp->b[jj][ii] ) *= QUART; /* j-th node vector: */
      } while (( ++jj ) < THREE );    /* ( rtp->b[j][*] ) */
   } while (( ++ii ) < THREE );       /* ( B* ) = adj(B) := ( rtp->b[*][*] ) */
/*............................................................................*/
/* Orthormalize node vectors b[j], yields intrinsic cell [ON] basis ub[]: */

   HMX_ORTNORM( );
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
/* Compute mesh cell skew: skv = (1/3)*trace(B*S). If ( skv != ZERO ) then   */
/* mark for skew in the cell ( skew := ONE ) */

   HMX_SKWVOL( ); /* skwvol := (1/3)*trace(B*S) */
/*............................................................................*/
/* Compute mesh cell volume */
/* [ N=0: vol = det(B) + (1/3)*trace(B*S) */
/*   N=1: vol = (1/3)*trace(B~*F) ]: */
/*   N=2: vol = (1/6) sum{j=1,2,3} ( b[j][]*( f[2j][]+f[2j+1][] )) */

   HMX_VOLUME(2);

   if ( fabs( rtp->vol ) < HMX_TRIVIAL )
   {
      fprintf( stderr, "\n\n Error message from function %s:", __func__ );
      fprintf( stderr, "\n Vanishing cell volume !!!" );

      ( rtp->rtn ) = null; /* unnormal return on error */
      return rtp;
   }
   else if(( fabs(( rtp->skv )/( rtp->vol ))) < ( 33.*HMX_PRECISION ))
      ( rtp->skew ) = null;
   else
      ( rtp->skew ) = ONE;
/*............................................................................*/
/* Debugging section [ display geometric cell characteristics ] */
/*............................................................................*/
# if HMX_DEBUG == 1
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
# elif HMX_DEBUG == 2
/* Display area vectors, opposite face vector means and skew vectors: */

   if( *ptr != 48 ) /* dec 48 = ASCII char '0' */
   {
      if( null != ( rtp->skew ))
      {
         printf( "\n\n Warning from function %s : ", __func__ );
         printf( "\n Geometric skew detected !!!" );
         printf( "\n [ Area vectors differ from opposite face vector means.]" );

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
# elif HMX_DEBUG == 3
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
	       ( rtp->mtx[ii][jj] ) += ( rtp->bu[ii][kk] )*( rtp->ub[kk][jj] );
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

# endif /* HMX_DEBUG == 3 */
/*............................................................................*/
# if HMX_CLLCRDS == 1
/* Use node and area vectors in cell coordinates [ with respect to basis ub ] */
/* This option yield incorrect results with anisotropic media !!! */

   if (( rtp->isotrop ) == ONE )
   {
      ii = null; do
      {
         jj = null; do
         {
            ( rtp->a[ii][jj] ) = ( rtp->au[ii][jj] );
            ( rtp->b[ii][jj] ) = ( rtp->bu[ii][jj] );
            ( rtp->p[ii][jj] ) = ( rtp->pu[ii][jj] );
            ( rtp->e[ii][jj] ) = ( rtp->eu[ii][jj] );
            ( rtp->cs[ii][jj] ) = ( rtp->cu[ii][jj] );
            ( rtp->f[ii][jj] ) = ( rtp->fu[ii][jj] );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );
      do
      {
         jj = null; do
         {
            ( rtp->p[ii][jj] ) = ( rtp->pu[ii][jj] );
            ( rtp->e[ii][jj] ) = ( rtp->eu[ii][jj] );
            ( rtp->f[ii][jj] ) = ( rtp->fu[ii][jj] );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < SIX );
      do
      {
         jj = null; do
         {
            ( rtp->p[ii][jj] ) = ( rtp->pu[ii][jj] );
            ( rtp->e[ii][jj] ) = ( rtp->eu[ii][jj] );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < PORTS );
   }; /* end if (( rtp->isotrop ) == 1 ) */
# endif /* HMX_CLLCRDS == 1 */
/*............................................................................*/
/* compute adj(A^-1) = ( A* )^-1: */

   ii = null; do        
   {
      jj = null; do
      {
         ( gjp->mr[ii][jj] ) = ( rtp->a[ii][jj] );
         ( gjp->mi[ii][jj] ) = ZERO;
      } while (( ++jj ) < THREE );
   } while (( ++ii ) < THREE );
/*............................................................................*/
/* matrix inversion: adj(A) |-> adj(A^-1): */

   ( gjp->rank ) = THREE;
   ( gjp->neqs ) = THREE;
   ( gjp->opt ) = 'i';

# if HMX_FULPIVT == 1
   gjp = gssjpv( gjp );
# else
   gjp = gssjrd( gjp );
# endif
/*............................................................................*/
   if ( gjp == NULL )
   { 
      printf( "\n\n Message from function %s : ", __func__ );
      printf( "\n Error on area vector matrix inversion "
         "adj(A) |-> adj(A^-1) !!!\n " );

      ( rtp->rtn ) = null;
      return rtp;
   };
/*............................................................................*/
/* copy adj(A^-1), then compute adj(B^-1): */

   ii = null; do           
   {
      jj = null; do
      { /* copy adj(A^-1), then enter adj(B) */
         ( rtp->ai[ii][jj] ) = ( gjp->zr[ii][jj] );
         ( gjp->mr[ii][jj] ) = ( rtp->b[ii][jj] );
         ( gjp->mi[ii][jj] ) = ZERO;
      } while (( ++jj ) < THREE );
   } while (( ++ii ) < THREE );
/*............................................................................*/
/* matrix inversion: adj(B) |-> adj(B^-1): */

   ( gjp->rank ) = THREE;
   ( gjp->neqs ) = THREE;
   ( gjp->opt ) = 'i';

# if HMX_FULPIVT == 1
   gjp = gssjpv( gjp );
# else
   gjp = gssjrd( gjp );
# endif
/*............................................................................*/
   if ( gjp == NULL )
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Error on node vector matrix inversion "
         "adj(B) |-> adj(B^-1) !!!\n " );

      ( rtp->rtn ) = null;
      return rtp;
   };
   if (( gjp->dtr ) < 1.e-77 )
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Non-positive node vector orientation !!!\n " );

      ( rtp->rtn ) = null;
      return rtp;
   }
   else if (( gjp->dtr ) <= ZERO )
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Degenerate node [ volume equals zero ] !!!\n " );

      ( rtp->rtn ) = null;
      return rtp;
   };
/*............................................................................*/
   ii = null; do
   {
      jj = null; do
      { /* copy adj(B^-1) */
         ( rtp->bi[ii][jj] ) = HMX_1_MINUS*( gjp->zr[ii][jj] );
      } while (( ++jj ) < THREE );
   } while (( ++ii ) < THREE );

/* here ends the option inedependent part ... */
/*............................................................................*/
/* ... applied in all the following options: */ 
/*............................................................................*/
/* option ( hsp->opt ) == 't': stability upper bound for DSC time step */
/* [ accurate for heat conduction and necessarily approximate for fluid flow  */
/*   - note that stability time step depends also on external sources ] */

   if ( opt == 't' ) /* compute stability upper bound for DSC time step */
   { 
      if ( *( rtp->ttyp ) == 't' ) /* trivial_e cell */
      {
         ( rtp->dt ) = HUGE_VALF;
         goto init;
      };
/*............................................................................*/
/* the stability time step: min ( ...dt ) */
/*...........................................................................*/
# if HMX_FCEMODE == 0     
      printf( "\n\n Error message from function %s :", __func__ );
      printf( "\n Can't compute thermal time step in HMX_FCEMODE = null !!!" );
      ( rtp->rtn ) = null; /* normal return: trivial cell */
      return rtp;
# else
/*...........................................................................*/
      ( rtp->hdt ) = HMX_INITSTP;
/*...........................................................................*/
# if DSC_FLDMDE != 0
      ( rtp->fdt ) = HMX_INITSTP;
# endif /* DSC_FLDMDE != 0 */
/*...........................................................................*/
      ii = null; do
      {
/* F*adj(B^-1): */

         xx = ZERO;
         jj = null; do
         {
            yy = ZERO;
            kk = null; do
            {
               yy += ( rtp->f[ii][kk] )*( rtp->bi[kk][jj] );
            } while (( ++kk ) < THREE );
            xx += ( yy*yy );
         } while (( ++jj ) < THREE );

         xx = sqrt( xx );

         if ( HMX_TRIVIAL < xx )
         {
            if ( HMX_TRIVIAL < ( rtp->kh ))
            {
               yy = (( rtp->cv )*( rtp->vol ))/( xx*( rtp->kh ));

               if( yy < ( rtp->hdt ))
                  ( rtp->hdt ) = yy;
            };
/*...........................................................................*/
# if DSC_FLDMDE != 0

            if ( *( rtp->ttyp ) == 'f' ) /* fluid type */
            {
               if ( HMX_TRIVIAL < ( rtp->ny ))
               {
                  yy = (( rtp->rm )*( rtp->vol ))/( xx*( rtp->ny ));

                  if( yy < ( rtp->fdt ))
                     ( rtp->fdt ) = yy;
               };
            };
# endif /* DSC_FLDMDE != 0 */
/*...........................................................................*/
         }; /* end if ( HMX_TRIVIAL < xx ) */
      } while (( ++ii ) < FACES );

      ( rtp->hdt ) *= QUART;
      ( rtp->dt ) = ( rtp->hdt );
/*...........................................................................*/
# if DSC_FLDMDE != 0
      ( rtp->fdt ) *= QUART;

      if (( rtp->fdt ) < ( rtp->dt ))
         ( rtp->dt ) = ( rtp->fdt );

# endif /* DSC_FLDMDE != 0 */
/*...........................................................................*/
      ( rtp->rtn ) = ONE; /* normal return: */
      return rtp;         /* stability time step determination terminated */
# endif /* HMX_FCEMODE != 0 */
/*...........................................................................*/
   } /* end if opt == 't'ime step */                
/*............................................................................*/
/* option ( hsp->opt ) == 's' */
/* heat and fluid flow s-parameters: */
/*............................................................................*/
   else if ( opt == 's' ) 
   {
      if ( *( rtp->ttyp ) == 't' ) /* trivial_e cell */
      {
         ( rtp->rtn ) = ONE; /* normal return: trivial cell */
         return rtp;
      };
/*............................................................................*/
# if HMX_FCEMODE == 0     
      printf( "\n\n Error message from function %s :", __func__ );
      printf( "\n Can't compute heat current s-parameters with "
         "HMX_FCEMODE = 0 !!!" );

      ( rtp->rtn ) = null; /* normal return: trivial cell */
      return rtp;
# else
      ( rtp->ct ) = ( rtp->dt )/( rtp->cv ); /* normalized temperature */
/*............................................................................*/
# if DSC_FLDMDE != 0
      if ( *( rtp->ttyp ) == 'f' )
         ( rtp->ft ) = ( rtp->dt )/( rtp->rm ); /* normalized fluid velocity */
      else                                      /* updating coefficient */
         ( rtp->ft ) = ZERO;
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      ii = null; do
      {
/* compute s[ii] = F[ii]*adj(B^-1): */

         jj = null; do
         {
            xx = ZERO;
            kk = null; do
            {
               xx += ( rtp->f[ii][kk] )*( rtp->bi[kk][jj] );
            } while (( ++kk ) < THREE );
            ( rtp->s[ii][jj] ) = xx;
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < FACES );
/*............................................................................*/
      ( rtp->rtn ) = ONE; /* normal return: */
      return rtp;         /* current s-parameter determination terminated */

# endif /* HMX_FCEMODE != 0 */
/*............................................................................*/
   }; /* end of option 's': nodal heat current s_parameters */
/*............................................................................*/
/* else: Error - This statement can never be reached in any legal option: */ 

   ( rtp->rtn ) = null;
   return rtp;
} 
/*============================================================================*/
# if defined ( OPTIMIZE )
   # pragma OPTIMIZE OFF
# endif
/*----------------------------------------------------------------------------*/
# if defined ( HMX_CLFCES )
   # undef HMX_CLFCES
# endif
# if defined ( HMX_CLSKEW )
   # undef HMX_CLSKEW
# endif
# if defined ( HMX_SKWVOL )
   # undef HMX_SKWVOL
# endif
# undef HMX_PORTS
# undef HMX_EDGES
# undef HMX_ORTNORM
# undef HMX_VOLUME
/*----------------------------------------------------------------------------*/
# undef HMX_CLLCRDS
# undef HMX_DSPINTM
# undef HMX_FULPIVT 
# undef HMX_FCEMODE
# undef HMX_PRECISION
# undef HMX_TRIVIAL
# undef HMX_INITSTP
# undef HMX_INCLUDE
# undef HMX_DEBUG
/*----------------------------------------------------------------------------*/
# undef CPRODUCT
# undef EPS_VAC
# undef MY_VAC_
# undef ELECTR_CHRGE
# undef ELECTR_MASS_
# undef BOHRs_MAGNTN
# undef PLANCKs_CNST
# undef STEFAN_BOLTZ
/***** end of DSC heat current s-parameters generation function hcrsmx(*) *****/
