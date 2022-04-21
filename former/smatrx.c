/* [ file: smatrx.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function smatrx(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   DSC S-matrix generation function supporting non-orthogonal hexahedral      *
*   mesh cell.                                                                 *
*                                                                              *
*   Given eigth points in 3-space, spt->c[i][] ( i=0,...,7 ), a DSC            *
*   time step, spt->dt, and material tensors spt->ep[][], spt->my[][],         *
*   spt->ke[][], spt->km[][] which respectively represent                      *
*   the relative permittivity, relative permeability, an the electric          *
*   and magnetic conductivities in the DSC cell,                               *
*                                                                              *
*   this function returns the nodal S-matrix                                   *
*                                                                              *
*      S  =  ( spt->se[j][k] ) x ( spt->sm[j][k] ) ; j,k = 0,...,5             *
*                                                                              *
*   and the set of form operators spt->a[][], spt->b[][],... for the           *
*   DSC mesh cell with vertex points spt->c[i][] [ option opt = 'r'].          *
*                                                                              *
*   Also [ in option opt = 't' ], given the vertex points spt->c[i][],         *
*   material tensors  spt->ep ,..., spt->km [ and optionally further           *
*   media parameters ] the stability upper bound for the DSC time step         *
*   of the cell with these vertex points is returned as spt->dt.               *
*                                                                              *
*   In option 's' , the S-matrix is returned in  block-diagonalized            *
*   form, viz. in the representation of the port voltage vector space          *
*   described in author's paper :                                              *
*   [1] 'Finite-Difference Time-Domain Approximation of Maxwell's              *
*   Equations with Non-orthogonal Condensed TLM Mesh' ,                        *
*   International Journal of Numerical Modelling, 1993                         *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
/*----------------------------------------------------------------------------*/
# if defined ( OPTIMIZE )
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2
# endif
/*----------------------------------------------------------------------------*/
# include "../math/consts.h"   /* various constants [null, ZERO, ONE,...]: */
/*----------------------------------------------------------------------------*/
# define SMX_INCLUDE 1      /* 1: include s-parameters transfer struct type */
# define SMX_DEBUG   0      /* >0: activate various printing functions */
                            /* [ mainly for debugging purposes ] */
/*----------------------------------------------------------------------------*/
/* operation marks:                                                           */

# define SMX_FULPIVT 1      /* 1: fully pivoted Gauss-Jordan algorithm        */
                            /*    [ more precise but CPU time expensive ]     */

# define SMX_FCEMODE 1      /* 0: don't compute any face vectors              */
                            /* 1: compute face vectors from cell ports        */
                            /* 2: compute face vectors from cell edges        */

/* At least one of the following three macros must be set to 1 */
/* to enable TRIVIAL CELL identification */ 

# define SMX_MEDCHCK 0      /* positivity and symmetry of material parameter  */
                            /* tensors is checked with MEDCHCK set to 1.      */
                            /* Also then, TRIVIAL CELLS are identified.       */

# define SMX_CLLCRDS 0      /* 1: forces cell specific node and area vectors. */
                            /* Also, identifies TRIVIAL CELLS. This option is */
                            /* inconsistent with unisotropic media */

# define SMX_ISOCHCK 1      /* 1: check media isotropy                        */
                            /* Also identifies TRIVIAL CELLS.                 */

# define SMX_UNICHCK 1      /* 1: check unitarity of frequency domain         */
                            /* S-matrix for lossless media.                   */

# define SMX_DSPINTM 0      /* SMX_DSPINTM  1: display intermediate results   */
/*----------------------------------------------------------------------------*/
# include "../math/CPRODUCT.H" /* complex product [ macro ] */
# include "../math/trdstp.h"   /* structure used in function triads(*) */
# include "../former/gssjtp.h" /* structure used in function gssjrd(*) */
# include "../former/jacbtp.h" /* structure used in function jacobi(*) */
/*----------------------------------------------------------------------------*/
/* computational parameters , regularization factors, bounds etc.:            */

# ifndef SMX_PRECISION
   # define SMX_PRECISION ( 1.e-14 )
# endif
# ifndef SMX_SQRPRECIS
   # define SMX_SQRPRECIS ( 1.e-08 )
# endif
# ifndef SMX_GIANTVAL
   # define SMX_GIANTVAL ( 1.e+277 )
# endif

# define SMX_DIAGBND ( 33.3*SMX_PRECISION )   /* bound for diagonality check  */
# define SMX_UNITBND ( 1.11*SMX_SQRPRECIS )   /* bound for unitary check      */

# define SMX_QRT_MINUS ( .25-1.*SMX_PRECISION ) /*'Secure'.25 used as limiter */
# define SMX_ONE_PLUS  ( 1.+77.*SMX_PRECISION ) /*'Secure' 1. used as bounds  */
# define SMX_ONE_MINUS ( 1.-77.*SMX_PRECISION ) /* for Hilbert norm of S-mtrx */
                                              /* block N, for instance,       */
                                              /* cf. DSC stability time step  */
                                              /* computation [ option 't' ].  */
# define SMX_TRIVBND ( 1.e-277 )              /* trivial cell recognitn bound */
# define SMX_ECRRBND ( 1.e-277 )              /* bound for recognition of     */
                                              /* electric current cell        */
# define SMX_MCRRBND ( 1.e-277 )              /* bound for recognition of     */
                                              /* magnetic current cell        */
/*----------------------------------------------------------------------------*/
/* natural constants - to be changed only with the help of God.               */
/* [ all units are international units (mks), if not specified otherwise;     */
/*   - not all these constants are used in this particular function version.] */

# define EPS_VAC      ( 8.8541878170e-12 ) /* vac. permittivity [A*sec/(V*m)] */
# define MY_VAC_      ( 1.2566370614e-06 ) /* "    permeability [V*sec/(A*m)] */
# define ELECTR_CHRGE (-1.6021773349e-19 ) /* electron charge [Coulomb=A*sec] */
# define ELECTR_MASS_ ( 9.1093897540e-31 ) /* electron mass [kg]              */
# define BOHRs_MAGNTN ( 1.1654071500e-29 ) /* Bohr's magneton [Volt*sec*m]    */
# define PLANCKs_CNST ( 6.6260755400e-34 ) /* Planck's constant [Joule*sec]   */
/*----------------------------------------------------------------------------*/
/* array dimensions:                                                          */

# define DIMNS  3
# define CRNRS  8
# define SORDR  6
# define PORTS 12
# define FACES  6
/*---------------------------------------------------------------------------*/
# if SMX_INCLUDE == 1
   # include "smxtyp.h"
# else
/*-----------------------------------------------------------------------------
*//* s-matrix data transfer structure [ in header 'smxctp.h' ], Release 5.2.
*//* All units are international units (mks), if not otherwise specified
*//* [ not all variables are used in this particular function version ]
------------------------------------------------------------------------------*/
typedef struct   
{
   signed char /* any return character [ may be used as an error flag, e.g. ] */
      rtn;

   char 
      opt; /* smatrx(*) options: 't' time_step, 'r': real s-parameters [ in */
           /* time-domain ], 'c': complex s-parameters [ in frequency-dmn ] */
   char
      loss, /* loss indicator [ lossless cell:0 / lossy cell:1 ] */
      skew, /* geometric skew indicator [ no skew:0 / skew:1] */
      isotrop; /* media isotropy indicator [ unisotropic: 0 / isotropic: 1 ] */ 

   short 
      med; /* media label */

   double
      adm, dt, dp, fr, omega, no, tr, tg, ld, vol, skv;

                       /* dt  = DSC time step                 [sec^-1       ] */
                       /* adm = y = field admittance          [Ohm^-1       ] */
                       /* no  = plasma electron density       [1/m^3        ] */
                       /* tr  = "      " relaxation time      [sec^-1       ] */
                       /* tg  = gyrmagn. relaxation time      [sec^-1       ] */
                       /* ld  = LANDE factor                  [dimensionless] */
   double
     mi[DIMNS],        /* mi = int.mag.fl.density,plasma[Tesla]               */
     ms[DIMNS],        /* ms = saturat. magnetization   [Tesla]               */
     hg[DIMNS],        /* hg = internal (static) magn.fld.[A/m]               */

     vp[FACES];        /* vp[jj] pyramide volume ( F[jj], b[jj/2] )/6         */

   double
    epr[DIMNS][DIMNS], /* ep = rel. permittivity tensor                       */
    epi[DIMNS][DIMNS],
    myr[DIMNS][DIMNS], /* my = rel. permeability "                            */
    myi[DIMNS][DIMNS],
     ke[DIMNS][DIMNS], /* ke = el.  conductivity "                            */
     km[DIMNS][DIMNS], /* km = magn conductivity "                            */

      c[CRNRS][DIMNS], /* corner point coordinates */

      e[PORTS][DIMNS], /* edge vectors */
     eu[PORTS][DIMNS], /* edge vectors relative to cell basis ub[]            */

      p[PORTS][DIMNS], /* port vectors */
     pu[PORTS][DIMNS], /* port vectors relative to cell basis ub[]            */

      f[FACES][DIMNS], /* face vectors */
     fu[FACES][DIMNS], /* face vectors relative to cell basis ub[]            */
     fa[DIMNS][DIMNS], /* opposite face vector means                          */

      a[DIMNS][DIMNS], /* area vector matrix A = arv[i][j]                    */
     au[DIMNS][DIMNS], /* area vector matrix relative to cell basis ub[]      */
     ai[DIMNS][DIMNS], /* Inverse are vector matrix A^-1                      */

      b[DIMNS][DIMNS], /* node vector matrix B = ndv[i][j]                    */
     bu[DIMNS][DIMNS], /* normalized Node vector matrix relative to ub[]      */
     ub[DIMNS][DIMNS], /* Gram-Schmidt orthonormalzd node vectors [cell basis]*/
     bi[DIMNS][DIMNS], /* inverse node vector basis B^-1                      */

     cs[DIMNS][DIMNS], /* cell geometric skew vectors                         */
     cu[DIMNS][DIMNS]; /* cell geometric skew vectors relative to cell basis  */

   double
     se[SORDR][SORDR], /* S-parameters [ electric type, time domain ]         */
    ser[DIMNS][DIMNS], /* S-parameters [ electric type, frq. dmn, real part ] */
    sei[DIMNS][DIMNS], /* S-parameters [ electric type, frq. dmn, imag.part ] */

     sm[SORDR][SORDR], /* S-parameters [ magnetic type ]                      */
    smr[DIMNS][DIMNS], /* S-parameters [ magnetic type, frq. dmn, real part ] */
    smi[DIMNS][DIMNS], /* S-parameters [ magnetic type, frq. dmn, imag.part ] */

     ge[DIMNS][DIMNS], /* gyroelectric bias [ not used in this function ]     */
     gm[DIMNS][DIMNS], /* gyromagnetic bias [ not used in this function ]     */

    te0[DIMNS][DIMNS], /* A*(ke/2+(EP+GE)/dt)*(B^-1)/4y                       */
    te1[DIMNS][DIMNS], /* A*(ke/2-(EP+GE)/dt)*(B^-1)/4y                       */
    tei[DIMNS][DIMNS], /* TE0^-1                                              */
    tet[DIMNS][DIMNS], /* (TE0^-1)*TE1/y                                      */
    tm0[DIMNS][DIMNS], /* y*A*(km/2+(MY+GM)/dt)*(B^-1)/4                      */
    tm1[DIMNS][DIMNS], /* y*A*(km/2-(MY+GM)/dt)*(B^-1)/4                      */
    tmi[DIMNS][DIMNS], /* TM0^-1                                              */
    tmt[DIMNS][DIMNS]; /* y*(TM0^-1)*TM1                                      */

   char 
     name[STS_SIZE], /* DSC system model name [ not used in this function ] */
     text[STS_SIZE], /* any comment [ not used in this function ] */
     etyp[SHS_SIZE], /* electric type indicator ['diagonal','symmetric', etc.]*/
     mtyp[SHS_SIZE], /* magnetic type indicator ['diagonal','symmetric', etc.]*/
     getp[SHS_SIZE], /* gyroelectric type identifier [ not used in this func] */
     gmtp[SHS_SIZE];

} S_MATRIX;
# endif /* SMX_INCLUDE != 1 */
/*----------------------------------------------------------------------------*/
/* system function prototypes:                                                */

void abort( void );
int fprintf( FILE *stream, const char *format,...);
int scanf( const char *format,...);
/*----------------------------------------------------------------------------*/
/* macros: */

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
/* compute mesh cell port vectors from vertex points: */
# define SMX_PORTS( ) \
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
   } while (( ++jj ) < DIMNS ); \
}
/*----------------------------------------------------------------------------*/
/* compute mesh cell edge vectors from vertex points: */

# define SMX_EDGES( ) \
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
   } while (( ++jj ) < DIMNS ); \
}
/*----------------------------------------------------------------------------*/
/* orthonormalize node vectors b[] to an intrinsic 'cell basis' ub[], and     */
/* transform node vectors into cell basis coordinates; this yields vectors    */
/* bu[]: */
# define SMX_ORTNORM( ) \
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
# if SMX_FCEMODE == 1
# define SMX_CLFCES( ) \
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
      jj = null; do \
      { \
         ( rtp->f[ii][jj] ) = ( trp->v[TWO][jj] ); \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < FACES ); \
} /* end of macro SMX_CLFCES */
/*----------------------------------------------------------------------------*/
# elif SMX_FCEMODE == 2
# define SMX_CLFCES( ) \
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
      jj = null; do \
      { \
         ( rtp->f[ii][jj] ) *= QUART; \
      } while (( ++jj ) < THREE ); \
   } while (( ++ii ) < FACES ); \
} /* end of macro SMX_CLFCES */
# endif /* SMX_FCEMODE == 2 */
/*----------------------------------------------------------------------------*/
# define SMX_CLSKEW( ) \
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
} /* end of macro SMX_CLSKEW */
/*----------------------------------------------------------------------------*/
/* The cell skew volume is defined as (1/3)*trace(B~*S), where B = ( b[i][j] )*/
/* and S = ( cs[i][j] ) are the node and skew vektor matrices, respectively, */
/* of the mesh cell: */

# define SMX_SKWVOL( ) \
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
} /* end of macro SMX_SKVOL */
/*----------------------------------------------------------------------------*/
/* There are different ways to compute the mesh cell volume, viz. by formulae */
/* (case NN=0)      vol = det(B) + (1/3) trace (B~*S)                         */
/* (case NN=1)      vol = (1/3) trace(B~*FM )                                 */
/* (case NN=2)      vol = (1/6) sum{j=1,2,3} ( b[j][]*( f[2j][]+f[2j+1][] ))  */
/* with the matrices B = (b[i][j]), S = ( cs[i][j] ), and F = ( fa[i][j] )    */

# define SMX_VOLUME(NN) \
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
} /* end of macro SMX_VOLUME */
/*----------------------------------------------------------------------------*/
static GAUSS_JRD gss = { null };
static JACOBI_EV jac = { null };
/*============================================================================*/

S_MATRIX *
smatrx( S_MATRIX *spt )
{
/* allusions: */
/*
   extern GAUSS_JRD gss;
   extern JACOBI_EV jac;
*/
/* declarations: */

   static S_MATRIX
      smx = {null},
     *rtp = &smx;

   static GAUSS_JRD *gjp = &gss;
   static JACOBI_EV *jev = &jac;

   static COMPLEX 
      cc = { null },
     *cpt = &cc;

   static TRIADS 
      trd = { null },
     *trp = &trd;

   static signed char
      opt = null;

   static char
     ptr[STS_SIZE] = {null},
     *trv_e = "trivial_E",
     *dgn_e = "diagonal_E",
     *sym_e = "symmetric_E",
     *asy_e = "asymmetric_E",
     *trv_m = "trivial_H",
     *dgn_m = "diagonal_H",
     *sym_m = "symmetric_H",
     *asy_m = "asymmetric_H";

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
      nrm1   = ZERO, 
      nrm2   = ZERO,
      hnorm1 = ZERO,
      hnorm2 = ZERO,
      bound1 = ZERO,
      bound2 = ZERO,
      dte    = ZERO,
      dtm    = ZERO,
      aa     = ZERO,
      rr     = ZERO,
      ss     = ZERO,
      uu     = ZERO,
      vv     = ZERO,
      xx     = ZERO,
      yy     = ZERO;

   static double
      edge[FOUR][DIMNS] = {{ZERO}},
      mtx1[DIMNS][DIMNS] = {{ZERO}},
      mtx2[DIMNS][DIMNS] = {{ZERO}};

/* [ so far not required; outcommented 2022.03.18 ]

   struct emcrr
   { 
      struct emcrr *ecp;

      double
         ec, mc,

         tue[DIMNS][DIMNS], 
          ve[DIMNS][DIMNS],
          pe[DIMNS][DIMNS],
          qe[DIMNS][DIMNS],
         eue[DIMNS][DIMNS],
         tum[DIMNS][DIMNS],
          vm[DIMNS][DIMNS],
          pm[DIMNS][DIMNS],
          qm[DIMNS][DIMNS],
         eum[DIMNS][DIMNS];
   };
   static struct emcrr emc = {null}; 
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

# if SMX_FULPIVT == 1  /* GAUSS JORDAN eleimination functions */
   GAUSS_JRD           /* [ used for matrix inversion, e.g. */
     *gssjpv( GAUSS_JRD *gjp ); /* full pivoting version */
# else
   GAUSS_JRD 
     *gssjrd( GAUSS_JRD *gjp ); /* line pivoting version */
# endif

   TRIADS
     *triads( TRIADS *trp ); /* various 3-vector operations [wedge prd,.e.g.] */

   JACOBI_EV 
     *jacobi( JACOBI_EV *jcp ); /* JACOBI eigenvalue determination function */
/*----------------------------------------------------------------------------*/
/* initialize/reset: */

   if ( spt == NULL )
   {
      rtp = &smx;

     init:

      ptr[null] = null;

      ( rtp->opt ) = null;
      ( rtp->loss ) = null;
      ( rtp->skew ) = null;
      ( rtp->isotrop ) = null;

      ii = null; do
      {
         ( rtp->name[ii] ) = null;
         ( rtp->text[ii] ) = null;
      } while (( ++ii ) < STS_SIZE );

      ( rtp->dt ) = SMX_GIANTVAL;
      ( rtp->fr ) = ZERO;
      ( rtp->no ) = ZERO;
      ( rtp->tr ) = ZERO;
      ( rtp->tg ) = ZERO;
      ( rtp->ld ) = ZERO;
      ( rtp->adm ) = ZERO;
      ( rtp->vol ) = ZERO;
      ( rtp->skv ) = ZERO;
      ( rtp->adm ) = sqrt( EPS_VAC/MY_VAC_ ); /* vacuum admittance [A/V] */
      ( rtp->omega ) = ZERO;

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
            ( rtp->fu[ii][jj] ) = ZERO;
         } while (( ++ii ) < FACES );

         ii = null; do
         {
            ( rtp->epr[ii][jj] ) = ZERO;
            ( rtp->epi[ii][jj] ) = ZERO;
            ( rtp->ke[ii][jj] ) = ZERO;
            ( rtp->myr[ii][jj] ) = ZERO;
            ( rtp->myi[ii][jj] ) = ZERO;
            ( rtp->km[ii][jj] ) = ZERO;

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

            ( rtp->te0[ii][jj] ) = ZERO;
            ( rtp->te1[ii][jj] ) = ZERO;
            ( rtp->tei[ii][jj] ) = ZERO;
            ( rtp->tet[ii][jj] ) = ZERO;

            ( rtp->tm0[ii][jj] ) = ZERO;
            ( rtp->tm1[ii][jj] ) = ZERO;
            ( rtp->tmi[ii][jj] ) = ZERO;
            ( rtp->tmt[ii][jj] ) = ZERO;

         } while (( ++ii ) < DIMNS );

         ii = null; do
         {
            ( rtp->ser[ii][jj] ) = ZERO;
            ( rtp->sei[ii][jj] ) = ZERO;
            ( rtp->smr[ii][jj] ) = ZERO;
            ( rtp->smi[ii][jj] ) = ZERO;
         } while (( ++ii ) < DIMNS );

         ( rtp->mi[jj] ) = ZERO;
         ( rtp->ms[jj] ) = ZERO;
         ( rtp->hg[jj] ) = ZERO;

         ii = null; do
         {
            ( rtp->ge[ii][jj] ) = ZERO;
            ( rtp->gm[ii][jj] ) = ZERO;
         } while (( ++ii ) < DIMNS );
      } while (( ++jj ) < DIMNS );

      ii = null; do
      {
         jj = null; do
         {
            ( rtp->se[ii][jj] ) = ZERO;
            ( rtp->sm[ii][jj] ) = ZERO;
         } while (( ++jj ) < SORDR );
      } while (( ++ii ) < SORDR );

      strncpy(( rtp->etyp ), trv_e, SHS_SIZE );
      strncpy(( rtp->mtyp ), trv_m, SHS_SIZE );
      strncpy(( rtp->getp ), trv_e, SHS_SIZE );
      strncpy(( rtp->gmtp ), trv_m, SHS_SIZE );

      ( rtp->rtn ) = ONE; /* normal return for structure initialization */
      return rtp;
   };
/*............................................................................*/
   rtp = spt;

   switch(( spt->opt ))
   {
     case 0:
      goto init;
      break;

     case 't':
     case 'T':
     case 'p':
     case 'P':

      opt = 't';
      ( rtp->dt ) = HUGE_VALF;
      break;

     case 's':
     case 'S':
     case 'r':
     case 'R':

      opt = 's';
      break;

     case 'c':
     case 'C':
     case 'f':
     case 'F':

      opt = 'f';
      break;

     default:
      fprintf( stdout, "\n\n Error message from function %s :", __func__ );
      fprintf( stdout,
         "\n Unknown or unspecified option on function call !!!" );
      fprintf( stdout, "\n Legal options are [ case insensitive ]:" );
      fprintf( stdout, "\n opt = 't' - compute stability time step" );
      fprintf( stdout, "\n opt = 'p' - compute stability phase shift" );
      fprintf( stdout,
         "\n opt = 's' - compute [ real ] time domain s-parameters" );
      fprintf( stdout, "\n opt = 'c' - compute [ complex ] frequency domain "
         "s-parameters" );
      fprintf( stdout, "\n [ Check function call in calling program.]\n " );

      ( rtp->rtn ) = null;
      return rtp;
   };
      
   ( rtp->adm ) = sqrt( EPS_VAC / MY_VAC_ );      /* vacuum admittance [A/V] */

/* [ so far not required; outcommented 2022.03.18 ]
   emc.ec = ELECTR_CHRGE/ELECTR_MASS_;          [ Coulomb/kg] 
   emc.mc = - 2.*PI*BOHRs_MAGNTN/PLANCKs_CNST;  [ m/Coulomb ] 
*/
   strcpy( rtp->etyp, trv_e );
   strcpy( rtp->getp, trv_e );
   strcpy( rtp->mtyp, trv_m );
   strcpy( rtp->gmtp, trv_m );

   ii = null; do
   { 
      jj = null; do
      {
         ( rtp->ge[ii][jj] ) = ZERO;
         ( rtp->gm[ii][jj] ) = ZERO;
      } while (( ++jj ) < DIMNS );
   } while (( ++ii ) < DIMNS );

   ( rtp->loss ) = null;
   ( rtp->skew ) = null;
   ( rtp->isotrop ) = null;
/*............................................................................*/
# if (( SMX_CLLCRDS == 1 )\
    ||( SMX_ISOCHCK == 1 ))

/* Media isotropy check: */

   ( rtp->isotrop ) = ONE;

   xx = fabs(( rtp->epr[null][null] ));
   yy = fabs(( rtp->ke[null][null] ));

   ( rtp->epr[null][null] ) = xx; /* [ must be non-negative ] */
   ( rtp->ke[null][null] ) = yy;  /* [ must be non-negative ] */

   ii = ONE; do
   {
      if ( SMX_TRIVBND < fabs( xx - ( rtp->epr[ii][ii] )))
      {
         ( rtp->isotrop ) = null;
# if SMX_CLLCRDS == 1
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Unisotropic dielectric media !!!" );
         fprintf( stderr,
            "\n This is incompatible with option SMX_CLLCRDS == 1 "
            "in function %s.", __func__ );
         fprintf( stderr, "\n [ Set function configuration macro "
            "SMX_CLLCRDS to null.]\n" );
         exit( EXIT_FAILURE );
# endif
      };

      if ( SMX_TRIVBND < fabs( yy - ( rtp->ke[ii][ii] )))
      {
         ( rtp->isotrop ) = null;
# if SMX_CLLCRDS == 1
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Unisotropic electric conductivity !!!" );
         fprintf( stderr,
            "\n This is incompatible with option SMX_CLLCRDS == 1 "
            "in function %s.", __func__ );
         fprintf( stderr, "\n [ Set function configuration macro "
            "SMX_CLLCRDS to null.]\n" );
         exit( EXIT_FAILURE );
# endif
      };

      jj = null; do
      {
         if( jj != ii )
         {
            if (( SMX_TRIVBND < fabs( rtp->epr[ii][jj] ))
              ||( SMX_TRIVBND < fabs( rtp->epr[jj][ii] )))
            {
               ( rtp->isotrop ) = null;
# if SMX_CLLCRDS == 1
               fprintf( stderr,
                  "\n\n Error message from function %s :", __func__ );
               fprintf( stderr, "\n Unisotropic dielectric media !!!" );
               fprintf( stderr,
                  "\n This is incompatible with option SMX_CLLCRDS == 1 "
                  "in function %s.", __func__ );
               fprintf( stderr, "\n [ Set function configuration macro "
                  "SMX_CLLCRDS to null.]\n" );
               exit( EXIT_FAILURE );
# endif
            };
            if (( SMX_TRIVBND < fabs( rtp->ke[ii][jj] ))
              ||( SMX_TRIVBND < fabs( rtp->ke[jj][ii] )))
            {
               ( rtp->isotrop ) = null;

# if SMX_CLLCRDS == 1
               fprintf( stderr,
                  "\n\n Error message from function %s :", __func__ );
               fprintf( stderr, "\n Unisotropic electric conductivity !!!" );
               fprintf( stderr,
                  "\n This is incompatible with option SMX_CLLCRDS == 1 "
                  "in function %s.", __func__ );
               fprintf( stderr, "\n [ Set function configuration macro "
                  "SMX_CLLCRDS to null.]\n" );
               exit( EXIT_FAILURE );
# endif
            };
         };
      } while (( ++jj ) < DIMNS );
   } while (( ++ii ) < DIMNS );

   if ( SMX_TRIVBND < yy )    /* loss = ONE: skip unitarity check */
      ( rtp->loss ) = ONE;    /* [ in the presence of losses ] */

   if ( xx < SMX_TRIVBND )    /* singular permittivity tensor rtp->epr */
   { 
      strcpy( rtp->etyp, trv_e );
      goto read_myr1;
   }
   else
      strcpy( rtp->etyp, dgn_e );

  read_myr1:

   xx = fabs(( rtp->myr[null][null] ));
   yy = fabs(( rtp->km[null][null] ));

   ( rtp->myr[null][null] ) = xx; /* [ must be non-negative ] */
   ( rtp->km[null][null] ) = yy;  /* [ must be non-negative ] */

   ii = ONE; do
   {
      if ( SMX_TRIVBND < fabs( xx - ( rtp->myr[ii][ii] )))
      {
         ( rtp->isotrop ) = null;
# if SMX_CLLCRDS == 1
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Unisotropic magnetic media !!!" );
         fprintf( stderr,
            "\n This is incompatible with option SMX_CLLCRDS == 1 "
            "in function %s.", __func__ );
         fprintf( stderr, "\n [ Set function configuration macro "
            "SMX_CLLCRDS to null.]\n" );
         exit( EXIT_FAILURE );
# endif
      };

      if ( SMX_TRIVBND < fabs( yy - ( rtp->km[ii][ii] )))
      {
         ( rtp->isotrop ) = null;
# if SMX_CLLCRDS == 1
         fprintf( stderr, "\n\n Error message from function %s :", __func__ );
         fprintf( stderr, "\n Unisotropic magnetic conductivity !!!" );
         fprintf( stderr,
            "\n This is incompatible with option SMX_CLLCRDS == 1 "
            "in function %s.", __func__ );
         fprintf( stderr, "\n [ Set function configuration macro "
            "SMX_CLLCRDS to null.]\n" );
         exit( EXIT_FAILURE );
# endif
      };

      jj = null; do
      {
         if( jj != ii )
         {
            if (( SMX_TRIVBND < fabs( rtp->myr[ii][jj] ))
              ||( SMX_TRIVBND < fabs( rtp->myr[jj][ii] )))
            {
               ( rtp->isotrop ) = null;
# if SMX_CLLCRDS == 1
               fprintf( stderr,
                  "\n\n Error message from function %s :", __func__ );
               fprintf( stderr, "\n Unisotropic magnetic media !!!" );
               fprintf( stderr,
                  "\n This is incompatible with option SMX_CLLCRDS == 1 "
                  "in function %s.", __func__ );
               fprintf( stderr, "\n [ Set function configuration macro "
                  "SMX_CLLCRDS to null.]\n" );
               exit( EXIT_FAILURE );
# endif
            };
            if (( SMX_TRIVBND < fabs( rtp->km[ii][jj] ))
              ||( SMX_TRIVBND < fabs( rtp->km[jj][ii] )))
            {
               ( rtp->isotrop ) = null;
# if SMX_CLLCRDS == 1
               fprintf( stderr,
                  "\n\n Error message from function %s :", __func__ );
               fprintf( stderr, "\n Unisotropic magnetic conductivity !!!" );
               fprintf( stderr,
                  "\n This is incompatible with option SMX_CLLCRDS == 1 "
                  "in function %s.", __func__ );
               fprintf( stderr, "\n [ Set function configuration macro "
                  "SMX_CLLCRDS to null.]\n" );
               exit( EXIT_FAILURE );
# endif
            };
         };
      } while (( ++jj ) < DIMNS );
   } while (( ++ii ) < DIMNS );

   if ( SMX_TRIVBND < yy )    /* loss = ONE: skip unitarity check */
      ( rtp->loss ) = ONE;    /* [ in the presence of losses ] */

   if ( xx < SMX_TRIVBND )    /* singular permeability tensor rtp->my */
   {
      if ( *( rtp->etyp ) == 't' )
         goto init; /* reset all parameters to ZERO and */
                    /* return for trivial mesh cell */
   }
   else
      strcpy( rtp->mtyp, dgn_m );

# endif /* if (( SMX_CLLCRDS == 1 )\
             ||( SMX_ISOCHCK == 1 )) */
/*............................................................................*/
# if SMX_MEDCHCK == 1
/* Check positivity and symmetry of material tensors, rtp->ep [ relative */
/* permittivity ] and rtp->my [ relative permeability ]. Also, this option */
/* identifies trivial cells */

   if (( rtp->isotrop ) == null )
   {
/* read_eps: */

      xx = ZERO;
      yy = ZERO;

      ii = null; do
      {
         jj = null; do
         {
            xx += fabs(( rtp->epr[ii][jj] ));
            yy += fabs(( rtp->ke[ii][jj] ));

            jac.h.r[ii][jj] = ( rtp->epr[ii][jj] );
            jac.h.i[ii][jj] = ( rtp->epi[ii][jj] );
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );

      if ( SMX_TRIVBND < yy )    /* loss = ONE: skip unitarity check */
         ( rtp->loss ) = ONE;    /* [ in the presence of losses ] */

      if ( xx < SMX_TRIVBND )    /* singular permittivity tensor rtp->epr */
         goto read_myr2;
      else 
         strcpy( rtp->etyp, dgn_e );
/*............................................................................*/
/* check permittivity tensor ( rtp->ep[][] ): */

      ( jev->rank ) = DIMNS;
/*............................................................................*/
      jev = jacobi( jev );     /* selfadjointness check, rtp->ep[][]          */
/*...........................*/
      if ( jev == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout, "\n Error returned from jacobi eigenvalue algorithm"
            "\n on computing spectrum of matrix [(ep)*(ke)] !!!" );
         ( rtp->rtn ) = null;

         return rtp;
      };

      ii = null; do
      {   
         if ( jac.h.r[ii][ii] < ZERO )
         {
            fprintf( stdout, "\n\n Message from function %s :", __func__ );
            fprintf( stdout,
               "\n Negative eigenvalue in matrix [(ep)*(ke)] !!!" );
            fprintf( stdout,
               "\n [ Check permittivity and conductivity tensors.]" );
            ( rtp->rtn ) = null;

            return rtp;
         };
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
 read_myr2:

      xx = ZERO;
      yy = ZERO;

      ii = null; do
      {
         jj = null; do
         {
            xx += fabs(( spt->myr[ii][jj] ));
            yy += fabs(( spt->km[ii][jj] ));

            jac.h.r[ii][jj] = ( spt->myr[ii][jj] );
            jac.h.i[ii][jj] = ( spt->myi[ii][jj] );
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );

      if ( SMX_TRIVBND < yy )    /* loss = ONE: skip unitarity check */
         ( rtp->loss ) = ONE;    /* [ in the presence of losses ] */

      if ( xx < SMX_TRIVBND )    /* singular permeability tensor rtp->my */
      {
         if ( *( rtp->etyp ) == 't' ) 
            goto init; /* reset all parameters to ZERO and */
                       /* return for trivial mesh cell */
      }
      else
         strcpy( rtp->mtyp, dgn_m );

/*............................................................................*/
/* check permeability tensor ( rtp->my[][] ): */

      ( jev->rank ) = DIMNS;
/*............................................................................*/
      jev = jacobi( jev );     /* selfadjointness check, rtp->my[][]          */
/*...........................*/
      if ( jev == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout, "\n Error returned from jacobi eigenvalue algorithm"
            "\n on computing spectrum of matrix [(my)*(km)] !!!" );
         ( rtp->rtn ) = null;
         return rtp;
      };

      ii = null; do
      {
         if ( jac.h.r[ii][ii] < ZERO )
         {
            fprintf( stdout, "\n\n Message from function %s :", __func__ );
            fprintf( stdout,
               "\n Negative eigenvalue in matrix [(my)*(km)] !!!" );
            ( rtp->rtn ) = null;
            return rtp;
         };
      } while (( ++ii ) < DIMNS );
   }; /* end if ( isotrop == null ) */

# else /* if ( SMX_MEDCHCK == null ) */

   strcpy( rtp->etyp, dgn_e );
   strcpy( rtp->mtyp, dgn_m );

# endif /* if ( SMX_MEDCHCK == null ) */
/*............................................................................*/
/* Port vector determination:

   Given 8 vertex points rtp->c[i][] ( i=0,...,7 )  port vectors  are
   computed  along  the 'parcel twines' scheme ( viz. the port vectors
   interconnect the midpoints of the mesh edges , cf. [1] ).
*/
   SMX_PORTS( );

/*............................................................................*/
/* Compute cell edges: */

   SMX_EDGES( );

/*............................................................................*/
/* Compute geometric skew: */

   SMX_CLSKEW( );

/*............................................................................*/
/* Compute face vectors [ from cell ports ]: */
# if SMX_FCEMODE != 0

   SMX_CLFCES( );

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
# endif /* SMX_FCEMODE != 0 */
/*............................................................................*/
/* Area and node vector determination:

   Area and node vectors are computed as port vector cycles and sums
   of port vectors, respectively. Node and area [line-] vectors
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
         } while( mm < THREE );
         ( rtp->a[jj][ii] ) *= QUART;            /* j-th area vector: */
      } while (( ++jj ) < DIMNS );               /* rtp->a[j][*]      */
   } while (( ++ii ) < DIMNS );                  /* A := rtp->a[*][*] */
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
         } while (( ++ll ) < FOUR );             /* j-th node vector: */
         ( rtp->b[jj][ii] ) *= QUART;            /* rtp->b[j][*]      */
      } while (( ++jj ) < DIMNS );               /* B := rtp->b[*][*] */
   } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* Orthormalize node vectors b[j], yields intrinsic cell [ON] basis ub[]: */

   SMX_ORTNORM( );

/*............................................................................*/
/* Transform area, face, port, and edge vectors into cell coordinates         */
/* [ i.e. with respect to basis ub[] ]: */

   ii = null; do
   {
      jj = null; do
      {
         ( rtp->au[ii][jj] ) = ZERO;
         ( rtp->pu[ii][jj] ) = ZERO;
         ( rtp->eu[ii][jj] ) = ZERO;
         ( rtp->cu[ii][jj] ) = ZERO;

# if SMX_FCEMODE != 0
         ( rtp->fu[ii][jj] ) = ZERO;
# endif
         kk = null; do
         {
            ( rtp->au[ii][jj] ) += ( rtp->a[ii][kk] )*( rtp->ub[jj][kk] );
            ( rtp->pu[ii][jj] ) += ( rtp->p[ii][kk] )*( rtp->ub[jj][kk] );
            ( rtp->eu[ii][jj] ) += ( rtp->e[ii][kk] )*( rtp->ub[jj][kk] );
            ( rtp->cu[ii][jj] ) += ( rtp->cs[ii][kk] )*( rtp->ub[jj][kk] );

# if SMX_FCEMODE != 0
            ( rtp->fu[ii][jj] ) += ( rtp->f[ii][kk] )*( rtp->ub[jj][kk] );
# endif
         } while (( ++kk ) < THREE );
      } while (( ++jj ) < THREE );
   } while (( ++ii ) < THREE );
   do
   {
      jj = null; do
      {
         ( rtp->pu[ii][jj] ) = ZERO;
         ( rtp->eu[ii][jj] ) = ZERO;

# if SMX_FCEMODE != 0
         ( rtp->fu[ii][jj] ) = ZERO;
# endif
         kk = null; do
         {
            ( rtp->pu[ii][jj] ) += ( rtp->p[ii][kk] )*( rtp->ub[jj][kk] );
            ( rtp->eu[ii][jj] ) += ( rtp->e[ii][kk] )*( rtp->ub[jj][kk] );
# if SMX_FCEMODE != 0
            ( rtp->fu[ii][jj] ) += ( rtp->f[ii][kk] )*( rtp->ub[jj][kk] );
# endif
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
/* Compute mesh cell skew volume  */
/* skv := (1/3)*trace(B~*S): */

   SMX_SKWVOL( ); /* skv := (1/3) trace(B~*S) */

/*............................................................................*/
/* Compute mesh cell volume using                                             */
/*           vol = det(B) + (1/3)* trace(B~*S)      [ case SMX_FCEMODE == 0 ] */
/*           vol = (1/3)*trace(B~*F),               [ case SMX_FCEMODE != 0 ] */
/* where B, S, and F denote the node, skew, and mean opposite face vector     */
/* matrices of the cell, respectively:                                        */

   SMX_VOLUME(SMX_FCEMODE);

   if(( fabs(( rtp->skv )/( rtp->vol ))) < ( 33.*SMX_PRECISION ))
      ( rtp->skew ) = null;
   else
      ( rtp->skew ) = ONE;
/*............................................................................*/
/* Debugging section [ display geometric cell characteristics ] */
/*............................................................................*/
# if SMX_DEBUG == 1
/* Display cell volume */

   if( *ptr != 48 ) /* dec 48 = ASCII char '0' */
   {
      fprintf( stdout, "\n det(B)=% .5e", ( trp->det ));
      fprintf( stdout, "\n vol(B)=% .5e", ( rtp->vol ));

      fprintf( stdout, "\n\n please acknowledge "
         "[ enter any character / Escape: enter 0 ]:" );
      scanf( "%s", ptr );

      if( *ptr == 48 )
         fprintf( stdout, "\n " );
   };
/*............................................................................*/
# elif SMX_DEBUG == 2
# if SMX_FCEMODE == 0
/* void */

   if( *ptr != 48 ) /* dec 48 = ASCII char '0' */
   {
      fprintf( stdout, "\n\n No cell geometry skew vektors computed"
         "\n [ Inactive face mode: SMX_FCEMODE = 0 ]" );

      fprintf( stdout, "\n\n please acknowledge "
         "[ enter any character ]:" );
      scanf( "%s", ptr );

      ptr[null] = 48;
      fprintf( stdout, "\n ");
   };
/*............................................................................*/
# elif SMX_FCEMODE != 0
/* Display area vectors, opposite face vector means and skew vectors: */

   if( *ptr != 48 ) /* dec 48 = ASCII char '0' */
   {
      if(( rtp->skew ) == ONE )
      {
         fprintf( stdout, "\n\n Warning from function %s :", __func__ );
         fprintf( stdout, "\n Geometric skew detected !!!" );
         fprintf( stdout,
            "\n [ Area vectors differ from opposite face vector means.]" );

         fprintf( stdout, "\n\n Area vectors [ port vector cycles ] "
            "- the k-th vector is a[k][*]:" );

         ii = null; do
         {
            fprintf( stdout, "\n |" );
            jj = null; do
            {
               fprintf( stdout,
                  " a[%2d][%2d]=% .5e|", jj, ii, rtp->a[jj][ii] );
            } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         fprintf( stdout, "\n\n Opposite face vector means "
            "- the k-th vector is fa[k][*]:" );

         ii = null; do
         {
            fprintf( stdout, "\n |" );
            jj = null; do
            {
               fprintf( stdout,
                  "fa[%2d][%2d]=% .5e|", jj, ii, ( rtp->fa[jj][ii] ));
            } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         fprintf( stdout, "\n\n Cell geometry skew vektors "
            "- the k-th vector is cs[k][*]:" );

         ii = null; do
         {
            fprintf( "\n |" );
            jj = null; do
            {
               fprintf( stdout,
                  "cs[%2d][%2d]=% .5e|", jj, ii, ( rtp->cs[jj][ii] ));
            } while (( ++jj ) < THREE );
         } while (( ++ii ) < THREE );

         fprintf( stdout, "\n\n please acknowledge "
            "[ enter any character / Escape: enter 0 ]:" );
         scanf( "%s", ptr );

         if( *ptr == 48 )
            fprintf( stdout, "\n " );
      };
   }; /* end if ptr[null] == null */
# endif /* SMX_FCEMODE != 0 */
/*............................................................................*/
# elif SMX_DEBUG == 3
/* Display node and area vectors in exterior end cell intrinsic coordinates */

   if( *ptr != 48 ) /* dec 48 = ASCII char '0' */
   {
      fprintf( stdout, "\n\n Node vectors "
         "- the k-th vector is b[k][*]:" );

      ii = null; do
      {
         fprintf( stdout, "\n |" );
         jj = null; do
	 {
	    fprintf( stdout,
               " b[%2d][%2d]=% .5e|", jj, ii, rtp->b[jj][ii] );
	 } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );

      fprintf( stdout, "\n\n Orthonormalized node vectors [ node vector basis ]"
	 "\n - the k-th vector is ub[k][*]:" );

      ii = null; do
      {
	 fprintf( stdout, "\n |" );
	 jj = null; do
	 {
	    fprintf( stdout,
               "ub[%2d][%2d]=% .5e| ", jj, ii, rtp->ub[jj][ii] );
	 } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );

      fprintf( stdout, "\n\n Node vectors in cell coordinates "
	 "[ relative to cell basis ub[] ]:" );

      ii = null; do
      {
	 fprintf( stdout, "\n |" );
	 jj = null; do
	 {
	    fprintf( stdout,
               "bu[%2d][%2d]=% .5e|", jj, ii, rtp->bu[jj][ii] );
	 } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );

      fprintf( stdout, "\n\n The product of the former two "
	 "[ yields the node vectors, again ]:" );

      ii = null; do
      {
	 jj = null; do
	 {
	    ( rtp->te0[ii][jj] ) = ZERO;
	    kk = null; do
	    { 
	       ( rtp->te0[ii][jj] ) += ( rtp->bu[ii][kk] )*( rtp->ub[kk][jj] );
	    } while (( ++kk ) < THREE );
	 } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );

      ii = null; do
      {
	 fprintf( stdout, "\n |" );
	 jj = null; do
	 {
	    fprintf( stdout,
               "te[%2d][%2d]=% .5e|", jj, ii, rtp->te0[jj][ii] );
	 } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );

      fprintf( stdout, "\n\n Area vectors [ port vector cycles ]"
	 "\n - the k-th vector is a[k][*]:" );

      ii = null; do
      {
         fprintf( stdout, "\n |" );
	 jj = null; do
	 {
	    fprintf( stdout,
               " a[%2d][%2d]=% .5e|", jj, ii, rtp->a[jj][ii] );
	 } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );

      fprintf( stdout, "\n\n Area vectors in cell coordinates "
         "[ with respect to cell basis ub[] ]:" );

      ii = null; do
      {
         fprintf( stdout, "\n |" );
         jj = null; do
         {
            fprintf( stdout,
               "au[%2d][%2d]=% .5e|", jj, ii, rtp->au[jj][ii] );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );

      fprintf( stdout, "\n\n please acknowledge "
         "[ enter any character / Escape: enter 0 ]:" );
      scanf( "%s", ptr );

      if( *ptr == 48 )
         fprintf( stdout, "\n " );

   }; /* end if ( *ptr != 48 ) */

# endif /* SMX_DEBUG == 3 */
/*............................................................................*/
# if (( SMX_CLLCRDS == 1 )\
    ||( SMX_ISOCHCK == 1 ))
/* Use node and area vectors in cell coordinates [ with respect to basis ub ] */
/* This option yields incorrect results with anisotropic media !!! */

   if (( rtp->isotrop ) == ONE )
   {
      ii = null; do
      {
         jj = null; do
         {
            ( rtp->a[ii][jj] ) = ( rtp->au[ii][jj] );
            ( rtp->b[ii][jj] ) = ( rtp->bu[ii][jj] );
# if SMX_FCEMODE != 0
            ( rtp->e[ii][jj] ) = ( rtp->eu[ii][jj] );
            ( rtp->f[ii][jj] ) = ( rtp->fu[ii][jj] );
# endif
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < THREE );
# if SMX_FCEMODE != 0
      do
      {
         jj = null; do
         {
            ( rtp->e[ii][jj] ) = ( rtp->eu[ii][jj] );
            ( rtp->f[ii][jj] ) = ( rtp->fu[ii][jj] );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < SIX );
      do
      {
         jj = null; do
         {
            ( rtp->f[ii][jj] ) = ( rtp->fu[ii][jj] );
         } while (( ++jj ) < THREE );
      } while (( ++ii ) < PORTS );
# endif /* SMX_FCEMODE != 0 */
   }; /* end if ( rtp->isotrop == 1 ) */
# endif /* ((SMX_CLLCRDS == 1 )||( SMX_ISOCHCK == 1 )) */
/*............................................................................*/
/* compute (A^-1): */

   ii = null; do        
   {
      jj = null; do
      {
         gss.mr[ii][jj] = ( rtp->a[ii][jj] );
         gss.mi[ii][jj] = ZERO;
      }  while (( ++jj ) < DIMNS );
   } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: A |-> (A^-1) */

   ( gjp->rank ) = DIMNS; /* the rank of the coeff matrix */
   ( gjp->neqs ) = DIMNS; /* the number of [simultaneously solved] eq sytems */
   ( gjp->opt ) = 'i';    /* option: matrix inversion */

# if SMX_FULPIVT == 1
   gjp = gssjpv( gjp );
# else
   gjp = gssjrd( gjp );
# endif
/*............................................................................*/
   if ( gjp == NULL )
   { 
      fprintf( stdout, "\n\n Message from function %s :", __func__ );
      fprintf( stdout,
         "\n Error on area vector matrix inversion A |-> (A^-1) !!!\n " );
      ( rtp->rtn ) = null;
      return rtp;
   };
/*............................................................................*/
/* copy matrix B into structure gss: */

   ii = null; do           
   {
      jj = null; do
      {
         ( rtp->ai[ii][jj] ) = gss.zr[ii][jj]; /* copy matrix (A^-1) */
         gss.mr[ii][jj] = ( rtp->b[ii][jj] );
         gss.mi[ii][jj] = ZERO;
      } while (( ++jj ) < DIMNS );
   } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: B |-> (B^-1) */

   ( gjp->rank ) = DIMNS; /* the rank of the coeff matrix */
   ( gjp->neqs ) = DIMNS; /* the number of [simultaneously solved] eq sytems */
   ( gjp->opt ) = 'i';    /* option: matrix inversion */

# if SMX_FULPIVT == 1
   gjp = gssjpv( gjp );
# else
   gjp = gssjrd( gjp );
# endif
/*............................................................................*/
   if ( gjp == NULL )
   {
      fprintf( stdout, "\n\n Message from function %s :", __func__ );
      fprintf( stdout,
         "\n Error on node vector matrix inversion B |-> (B^-1) !!!\n " );
      ( rtp->rtn ) = null;
      return rtp;
   };
   if ( gss.dtr < 1.e-77 )
   {
      fprintf( stdout, "\n\n Message from function %s :", __func__ );
      fprintf( stdout, "\n Negative node vector orientation !!!\n " );
      ( rtp->rtn ) = null;
      return rtp;
   }
   else if ( gss.dtr <= ZERO )
   {
      fprintf( stdout, "\n\n Message from function %s :", __func__ );
      fprintf( stdout, "\n Degenerate node [ volume equals zero ] !!!\n " );
      ( rtp->rtn ) = null;
      return rtp;
   };

   ii = null; do
   {
      jj = null; do
      {
         ( rtp->bi[ii][jj] ) = gss.zr[ii][jj]; /* copy matrix (B^-1) */
      } while (( ++jj ) < DIMNS );
   } while (( ++ii ) < DIMNS );
/*-------------------- end of general ( common ) part ------------------------*/
/* options: */

   if ( opt == 't' )
/*.............................................................................
   option 't' :   compute stability upper bound for DSC time step
.............................................................................*/
   {
/* tstep_e: */

      if ( *( rtp->etyp ) == 't' ) /* trivial_e cell: In the main program     */
      {                            /* the 'trivial cell' is encapsulated by   */
         dte = HUGE_VALF;         /* electric walls. The scattering cycle    */
         goto tstep_m;             /* overrides that cell, while the connec-  */
      };                           /* tion cycle encounters the electric walls*/
/*............................................................................*/
/* copy matrix ep into structure gss: */

      if (( rtp->isotrop ) == null )
      {
         ii = null; do
         {
            jj = null; do
            {
               gss.mr[ii][jj] = ( rtp->epr[ii][jj] );
               gss.mi[ii][jj] = ZERO;
            } while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: ep |-> (ep^-1) */

         ( gjp->rank ) = DIMNS; /* the rank */
         ( gjp->neqs ) = DIMNS; 
         ( gjp->opt ) = 'i';    /* otion inversion */

# if SMX_FULPIVT == 1
         gjp = gssjpv( gjp );
# else
         gjp = gssjrd( gjp );
# endif
/*............................................................................*/
         if ( gjp == NULL )
         {
            fprintf( stdout, "\n\n Message from function %s :", __func__ );
            fprintf( stdout,
               "\n Error on matrix inversion ep |-> (ep^-1) !!!\n " );
            ( rtp->rtn ) = null;
            return rtp;
         };

         ii = null; do /* matrix product: (ep^-1)*(A^-1) */
         {
            jj = null; do
            {
               mtx1[ii][jj] = ZERO;
               kk = null; do
               {
                  mtx1[ii][jj] += ( gss.zr[ii][kk] * ( rtp->ai[kk][jj] ));
               } while (( ++kk ) < DIMNS );
            } while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
      } /* end if ( isotrop == null ) */
      else /* if ( isotrop != null ) */
      {
         ii = null; do /* matrix product: (ep^-1)*(A^-1) */
         {
            jj = null; do
            {
               mtx1[ii][jj] = ( rtp->ai[ii][jj] )/( rtp->epr[null][null] );
            } while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
      }; /* end if ( isotrop != null ) */
/*............................................................................*/
/*  Compute matrix product [ B * (ep^-1) * (A^-1) ]: */

      ii = null; do
      {
         jj = null; do
         {
            jac.h.r[ii][jj] = ZERO;
            jac.h.i[ii][jj] = ZERO;

            kk = null; do
            {
               jac.h.r[ii][jj] += (( rtp->b[ii][kk] ) * mtx1[kk][jj] );
            } while (( ++kk ) < DIMNS );
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );

      ( jev->rank ) = DIMNS;
/*............................................................................*/
      jev = jacobi( jev );     /*  returns the eigenvalues of matrix          */
/*...........................*//*  T = B*(ep^-1)*(A^-1)                       */
      if ( jev == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout,
            "\n Error returned from jacobi eigenvalue algorithm"
            "\n on computing spectrum of matrix B*(ep^-1)*(A^-1) !!!\n " );
         ( rtp->rtn ) = null;
         return rtp;
      };

      ii = null; do
      {
         if ( jac.h.r[ii][ii] <= ZERO )
         {
            fprintf( stdout, "\n\n Message from function %s :", __func__ );
            fprintf( stdout, "\n Non-positive eigenvalue "
               "in matrix B*(ep^-1)*(A^-1) !!!" );
            fprintf( stdout, "\n [ Check permittivity tensor.]\n " );
            ( rtp->rtn ) = null;
            return rtp;
         };
      } while (( ++ii ) < DIMNS );

      hnorm1 = jac.hn;
/*
      fprintf( stdout, "\n\n Hilbert norm :  H = %.16le\n ", hnorm1 );
*/
      dte = EPS_VAC/( rtp->adm )/hnorm1/2.;
/*............................................................................*/
/* compute matrix (my^-1): */

     tstep_m:

      if ( *( rtp->mtyp ) == 't' ) /* trivial_m cell */
      {
         dtm = HUGE_VALF;
         goto tstep;
      };
/*............................................................................*/
/* copy matrix my into structure gss: */

      if (( rtp->isotrop ) == null )
      {
         ii = null; do
         {
            jj = null; do
            {
               gss.mr[ii][jj] = ( rtp->myr[ii][jj] );
               gss.mi[ii][jj] = ZERO;
            } while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: my |-> (my^-1) */

         ( gjp->rank ) = DIMNS; /* the rank */
         ( gjp->neqs ) = DIMNS; 
         ( gjp->opt ) = 'i';    /* inversion */

# if SMX_FULPIVT == 1
         gjp = gssjpv( gjp );
# else
         gjp = gssjrd( gjp );
# endif
/*............................................................................*/
         if ( gjp == NULL )
         {
            fprintf( stdout, "\n\n Message from function %s :", __func__ );
            fprintf( stdout,
               "\n Error on matrix inversion my |-> (my^-1) !!!\n " );
            ( rtp->rtn ) = null;
            return rtp;
         };
/*............................................................................*/
/* matrix product (my^-1)*(A^-1): */

         ii = null; do
         {
            jj = null; do
            {
               mtx1[ii][jj] = ZERO;
               kk = null; do
               {
                  mtx1[ii][jj] += ( gss.zr[ii][kk] * ( rtp->ai[kk][jj] ));
               } while (( ++kk ) < DIMNS );
            }  while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
      } /* end if ( isotrop == null ) */
      else /* if ( isotrop == ONE ) */
      {
/*............................................................................*/
/* matrix product (my^-1)*(A^-1): */

         ii = null; do
         {
            jj = null; do
            {
               mtx1[ii][jj] = ( rtp->ai[ii][jj] )/( rtp->myr[null][null] );
            } while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
      }; /* end if ( isotrop == ONE ) */
/*............................................................................*/
/* matrix product B*(my^-1)*(A^-1): */

      ii = null; do
      {
         jj = null; do
         {
            jac.h.r[ii][jj] = ZERO;
            jac.h.i[ii][jj] = ZERO;

            kk = null; do
            {
               jac.h.r[ii][jj] += (( rtp->b[ii][kk] ) * mtx1[kk][jj] );
            } while (( ++kk ) < DIMNS );
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );

      ( jev->rank ) = DIMNS;
/*............................................................................*/
      jev = jacobi( jev );     /*  returns the eigenvalues of matrix          */
/*...........................*//*  T = B*(my^-1)*(A^-1)                       */
      if ( jev == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout, "\n Error returned from jacobi eigenvalue algorithm"
            "\n on computing spectrum of matrix B*(my^-1)*(A^-1) !!!\n " );
         ( rtp->rtn ) = null;
         return rtp;
      };

      ii = null; do
      {
         if ( jac.h.r[ii][ii] <= ZERO )
         {
            fprintf( stdout, "\n\n Message from function %s :", __func__ );
            fprintf( stdout, "\n Non-positive eigenvalue "
               "in matrix B*(my^-1)*(A^-1) !!!" );
            fprintf( stdout, "\n [ Check permeability tensor.]\n " );
            ( rtp->rtn ) = null;
            return rtp;
         };
      } while (( ++ii ) < DIMNS );

      hnorm2 = jac.hn;
/*
      fprintf( stdout, "\n\n Hilbert norm :  H = %.16le\n ", hnorm2 );
*/
      dtm = MY_VAC_*( rtp->adm )/hnorm2/2.;
/*............................................................................*/
/* the stability time step equals min ( dte , dtm )                           */

     tstep: 

      if ( dtm < dte ) /* rtp->dt = stability upper bound for DSC time step */
         ( rtp->dt ) = dtm;
      else
         ( rtp->dt ) = dte;
       
      ( rtp->dt ) *= SMX_ONE_MINUS;

      ( rtp->rtn ) = ONE;
      return rtp;         /* stability time step computation terminated */

   } /* end if opt == 't'ime step */                

   else if ( opt == 's' ) /* option: compute 's'-parameters */
   {
/*.............................................................................
      1. S-parameters derived from the generalized Amprere's law
..............................................................................*/
/* Amperes_law: */

      if ( *( rtp->etyp ) == 't' ) /* trivial_e cell */
         goto Faradays_law;
      else
      {
         strcpy( rtp->getp, trv_e );

         ii = null; do
         {
            jj = null; do
            {
               ( rtp->ge[ii][jj] ) = ZERO;
            } while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
      }; /* end if ( |rtp->bi| <= SMX_ECRRBND ) */
/*............................................................................*/
/* compute matrix TE0(1) = A*( ke/2 +(-)(EP/dt + GE))*(B^-1)/y/4: */

      aa = EPS_VAC/( rtp->dt );

      ii = null; do
      {
         jj = null; do
         {
            uu = ZERO;
            vv = ZERO;
            kk = null; do
            {
               xx = ( rtp->ke[kk][jj] )/2. + aa * ( rtp->epr[kk][jj] );
               xx += ( rtp->ge[kk][jj] );
               uu += (( rtp->a[ii][kk] ) * xx );
               xx = ( rtp->ke[kk][jj] )/2. - aa * ( rtp->epr[kk][jj] );
               xx -= ( rtp->ge[kk][jj] );
               vv += (( rtp->a[ii][kk] ) * xx );
            } while (( ++kk ) < DIMNS );
            mtx1[ii][jj] = uu / ( rtp->adm ) / 4.;
            mtx2[ii][jj] = vv / ( rtp->adm ) / 4.;
         } while (( ++jj ) < DIMNS );

         jj = null; do
         {
            uu = ZERO;
            vv = ZERO;

            kk = null; do
            {
               uu += ( mtx1[ii][kk] * ( rtp->bi[kk][jj] ));
               vv += ( mtx2[ii][kk] * ( rtp->bi[kk][jj] ));
            } while (( ++kk ) < DIMNS );
            ( rtp->te0[ii][jj] ) = uu;
            ( rtp->te1[ii][jj] ) = vv;
            gss.mr[ii][jj] = uu; 
            gss.mi[ii][jj] = ZERO;
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/*  matrix inversion: TE0 |-> (TE0^-1) */

         ( gjp->rank ) = DIMNS; /* the rank */
         ( gjp->neqs ) = DIMNS; 
         ( gjp->opt ) = 'i';    /* inversion */

# if SMX_FULPIVT == 1
         gjp = gssjpv( gjp );
# else
         gjp = gssjrd( gjp );
# endif
/*............................................................................*/
      if ( gjp == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s : ", __func__ );
         fprintf( stdout, "\n Error on matrix inversion T |-> (T^-1) !!!\n " );
         ( rtp->rtn ) = null;
         return rtp;
      };

      nrm1 = ZERO;
      nrm2 = ZERO;

      ii = null; do  
      {
         jj = null; do
         {
            ( rtp->tei[ii][jj] ) = gss.zr[ii][jj];

            xx = ( rtp->tei[ii][jj] ) * ( rtp->tei[ii][jj] );

            if ( nrm1 < xx )
               nrm1 = xx;

         } while (( ++jj ) < DIMNS );  

         jj = null; do
         {
            ( rtp->tet[ii][jj] ) = ZERO; 

            kk = null; do 
            {
               ( rtp->tet[ii][jj] ) += \
                           (( rtp->tei[ii][kk] ) * ( rtp->te1[kk][jj] ));
            } while (( ++kk ) < DIMNS );

            xx = ( rtp->tet[ii][jj] ) * ( rtp->tet[ii][jj] );

            if ( nrm2 < xx )
               nrm2 = xx;

         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );

      nrm1   = sqrt( nrm1 );
      bound1 = nrm1*SMX_DIAGBND;
      nrm2   = sqrt( nrm2 );
      bound2 = nrm2*SMX_DIAGBND;
/*............................................................................*/
/* check if (TE0^-1) or (TE0^-1)*TE1 diagonal : */

      strcpy( rtp->etyp, dgn_e );            

      ii = ONE; do            
      {
         jj = null; do 
         {
            if ( *rtp->etyp  == 'd' )
            {
               xx = fabs( rtp->tei[ii][jj] ) + fabs( rtp->tei[jj][ii] );
               yy = fabs( rtp->tet[ii][jj] ) + fabs( rtp->tet[jj][ii] );

               if (( bound1 < xx )
                 ||( bound2 < yy ))
                  strcpy(( rtp->etyp ), sym_e );
            };
            if ( *( rtp->etyp ) == 's' )
            {
               xx = fabs(( rtp->tei[ii][jj] ) - ( rtp->tei[jj][ii] ));
               yy = fabs(( rtp->tet[ii][jj] ) - ( rtp->tet[jj][ii] ));

               if (( bound1 < xx )
                 ||( bound2 < yy ))
                  strcpy(( rtp->etyp ), asy_e );
            };
         } while (( ++jj ) < ii );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* compute NE = mat1 = -(TE0^-1) - TET; ME = mat2 = (Id + mat1 )*(TE0^-1)    */

      ii = null; do
      {                          
         jj = null; do
         {
            mtx1[ii][jj] = - (( rtp->tei[ii][jj] ) + ( rtp->tet[ii][jj] )); 
         } while (( ++jj ) < DIMNS );

         jj = null; do
         {
            mtx2[ii][jj] = ZERO;

            kk = null; do
            {
               mtx2[ii][jj] += \
                   ((( double )( ii == kk ) + mtx1[ii][kk] ) *
                    ( rtp->tei[kk][jj] ));
            } while (( ++kk ) < DIMNS );
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*----------------------------------------------------------------------------*/
/* Write E-type nodal S-parameters: */

      bound1 = 33.*SMX_PRECISION;

      ii = null; do 
      {
         jj = null; do
         {
/*............................................................................*/
/* block KE: */
            uu = ( rtp->tei[ii][jj] );                              
            vv = fabs( uu );

            if( vv < bound1 )
               ( rtp->se[ii][jj] ) = - ( double )( ii == jj );
            else if( fabs( 1. - vv ) < bound1 )
               ( rtp->se[ii][jj] ) = uu/vv - ( double )( ii == jj );
            else
               ( rtp->se[ii][jj] ) = uu - ( double )( ii == jj );
/*............................................................................*/
/* block LE:                                                                  */
            ( rtp->se[ii][jj+DIMNS] ) = ( double )( ii == jj );            
/*............................................................................*/
/* block ME:                                                                  */
            uu = mtx2[ii][jj];  
            vv = fabs( uu );   

            if( vv < bound1 )
               ( rtp->se[ii+DIMNS][jj] ) = ZERO;
            else if( fabs( 1. - vv ) < bound1 )
               ( rtp->se[ii+DIMNS][jj] ) = uu/vv;
            else
               ( rtp->se[ii+DIMNS][jj] ) = uu;
/*............................................................................*/
/* block NE:                                                                  */
            uu = mtx1[ii][jj];                                
            vv = fabs( uu );   

            if( vv < bound1 )
               ( rtp->se[ii+DIMNS][jj+DIMNS] ) = ZERO;
            else if( fabs( 1. - vv ) < bound1 )
               ( rtp->se[ii+DIMNS][jj+DIMNS] ) = uu/vv;
            else
               ( rtp->se[ii+DIMNS][jj+DIMNS] ) = uu;
/*............................................................................*/
/* Hilbert norm:                                                              */
            jac.h.r[ii][jj] = ZERO;
            jac.h.i[ii][jj] = ZERO;

            kk = null; do
            {
               jac.h.r[ii][jj] += ( mtx1[ii][kk] * mtx1[jj][kk] );
            } while (( ++kk ) < DIMNS ); 
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );

      ( jev->rank ) = DIMNS;
/*............................................................................*/
      jev = jacobi( jev );     /*  returns the Hilbert [spectral] norm of NE  */
/*...........................*//*                                             */
      if ( jev == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout, "\n Error on Hilbert norm computation "
            "in matrix NE !!!\n ");
         ( rtp->rtn ) = null;
         return rtp;
      };
      if ( SMX_ONE_PLUS < sqrt( jac.hn ) )
      {  
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout, "\n Sup norm of electric mode S-matrix"
            " ( block NE ) = %.16e > 1 .", sqrt( jac.hn ) );
         fprintf( stdout, "\n - This may yield unstable DSC process !" );
         fprintf( stdout, "\n [ Should be remedied "
            "by chosing a smaller DSC time step.]\n " );
         ( rtp->rtn ) = null;
         return rtp;
      };           
/*.............................................................................
   2. S-parameters derived from Faraday's law
..............................................................................*/
     Faradays_law:

      if ( *( rtp->mtyp ) == 't' ) /* trivial_m cell */
      {
            goto init; /* reset all parameters to ZERO and */
                       /* return for trivial mesh cell */
      }
      else
      {
         strcpy(( rtp->gmtp ), trv_e );

         ii = null; do
         {
            jj = null; do
            {
               ( rtp->gm[ii][jj] ) = ZERO;
            } while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
      }; /* end if ( |rtp->ms| <= SMX_MCRRBND ) */
/*............................................................................*/
/* compute matrix TM0(1) = y*A*( km/2 +(-)(MY/dt + GM))*(B^-1)/4 : */

      aa = MY_VAC_/( rtp->dt );

      ii = null; do
      {
         jj = null; do
         {
            uu = ZERO;
            vv = ZERO;

            kk = null; do
            {
               xx = ( rtp->km[kk][jj] )/2. + aa * ( rtp->myr[kk][jj] );
               xx += ( rtp->gm[kk][jj] );
               uu += (( rtp->a[ii][kk] ) * xx );
               xx = ( rtp->km[kk][jj] )/2. - aa * ( rtp->myr[kk][jj] );
               xx -= ( rtp->gm[kk][jj] );
               vv += (( rtp->a[ii][kk] ) * xx );
            } while (( ++kk ) < DIMNS );
            mtx1[ii][jj] = uu * ( rtp->adm ) / 4.;
            mtx2[ii][jj] = vv * ( rtp->adm ) / 4.; 
         } while (( ++jj ) < DIMNS );

         jj = null; do
         {
            uu = ZERO;
            vv = ZERO;

            kk = null; do
            {
               uu += ( mtx1[ii][kk] * ( rtp->bi[kk][jj] ));
               vv += ( mtx2[ii][kk] * ( rtp->bi[kk][jj] ));
            }  while (( ++kk ) < DIMNS );
            ( rtp->tm0[ii][jj] ) = uu;
            ( rtp->tm1[ii][jj] ) = vv;
            gss.mr[ii][jj] = uu; 
            gss.mi[ii][jj] = ZERO;
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: TM0 |-> (TM0^-1) */

      ( gjp->rank ) = DIMNS; /* the rank */
      ( gjp->neqs ) = DIMNS; 
      ( gjp->opt ) = 'i';    /* matrix inversion */

# if SMX_FULPIVT == 1
      gjp = gssjpv( gjp );
# else
      gjp = gssjrd( gjp );
# endif
/*............................................................................*/
      if ( gjp == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout,
            "\n Error on matrix inversion TM0 |-> (TM0^-1) !!!\n " );
         ( rtp->rtn ) = null;
         return rtp;
      };

      nrm1 = ZERO;
      nrm2 = ZERO; 

      ii = null; do  
      {
         jj = null; do
         {
            ( rtp->tmi[ii][jj] ) = gss.zr[ii][jj];

            xx = ( rtp->tmi[ii][jj] ) * ( rtp->tmi[ii][jj] );

            if ( nrm1 < xx )
               nrm1 = xx;

         } while (( ++jj ) < DIMNS );   

         jj = null; do
         {
            ( rtp->tmt[ii][jj] ) = ZERO; 

            kk = null; do
            {
               ( rtp->tmt[ii][jj] ) += 
                      (( rtp->tmi[ii][kk] ) * ( rtp->tm1[kk][jj] ));
            } while (( ++kk ) < DIMNS );

            xx = ( rtp->tmt[ii][jj] ) * ( rtp->tmt[ii][jj] );

            if ( nrm2 < xx )
               nrm2 = xx;

         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );

      nrm1   = sqrt( nrm1 );
      bound1 = nrm1*SMX_DIAGBND;
      nrm2   = sqrt( nrm2 );
      bound2 = nrm2*SMX_DIAGBND;
/*............................................................................*/
/* check if (TM0^-1) or (TM0^-1)*TM1 diagonal: */

      strcpy(( rtp->mtyp ), dgn_m );

      ii = ONE; do            
      {
         jj = null; do
         {
            if ( *rtp->mtyp  == 'd' )
            {
               xx = fabs( rtp->tmi[ii][jj] ) + fabs( rtp->tmi[jj][ii] );
               yy = fabs( rtp->tmt[ii][jj] ) + fabs( rtp->tmt[jj][ii] );

               if (( bound1 < xx )
                 ||( bound2 < yy )) 
                  strcpy(( rtp->mtyp ), sym_m );
            };
            if ( *rtp->mtyp == 's' )
            {
               xx = fabs(( rtp->tmi[ii][jj] ) - ( rtp->tmi[jj][ii] ));
               yy = fabs(( rtp->tmt[ii][jj] ) - ( rtp->tmt[jj][ii] ));

               if (( bound1 < xx )
                 ||( bound2 < yy ))
                  strcpy(( rtp->mtyp ), asy_m );
            };
         } while (( ++jj ) < ii );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* compute NH = mat1 = -(TM0^-1) - TMT; MH = mat2 = - (Id + mat1 )*(TM0^-1)   */

      ii = null; do
      {                          
         jj = null; do
         {
            mtx1[ii][jj] = - (( rtp->tmi[ii][jj] ) + ( rtp->tmt[ii][jj] )); 
         } while (( ++jj ) < DIMNS );

         jj = null; do
         {
            mtx2[ii][jj] = ZERO;

            kk = null; do
            {
               mtx2[ii][jj] -= \
                (((double)( ii == kk ) + mtx1[ii][kk] ) * ( rtp->tmi[kk][jj] ));
            } while (( ++kk ) < DIMNS );
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*----------------------------------------------------------------------------*/
/* Write H-type nodal S-parameters: */

      bound1 = 33.*SMX_PRECISION;

      ii = null; do 
      {
         jj = null; do
         {
/*............................................................................*/
/* block KH: */

            uu = - ( rtp->tmi[ii][jj] );                              
            vv = fabs( uu );

            if( vv < bound1 )
               ( rtp->sm[ii][jj] ) = ( double )( ii == jj );
            else if( fabs( 1. - vv ) < bound1 )
               ( rtp->sm[ii][jj] ) = uu/vv + ( double )( ii == jj );
            else
               ( rtp->sm[ii][jj] ) = uu + ( double )( ii == jj );
/*............................................................................*/
/* block LH: */

            ( rtp->sm[ii][jj+DIMNS] ) = ( double )( ii == jj );
/*............................................................................*/
/* block MH: */

            uu = mtx2[ii][jj];                         
            vv = fabs( uu );

            if ( vv < bound1 )
               ( rtp->sm[ii+DIMNS][jj] ) = ZERO;
            else if( fabs( 1. - vv ) < bound1 )
               ( rtp->sm[ii+DIMNS][jj] ) = uu/vv;
            else
               ( rtp->sm[ii+DIMNS][jj] ) = uu;
/*............................................................................*/
/* block NH: */

            uu = mtx1[ii][jj];                              
            vv = fabs( uu );

            if( vv < bound1 )
               ( rtp->sm[ii+DIMNS][jj+DIMNS] ) = ZERO;
            else if( fabs( 1. - vv ) < bound1 )
               ( rtp->sm[ii+DIMNS][jj+DIMNS] ) = uu / vv;
            else
               ( rtp->sm[ii+DIMNS][jj+DIMNS] ) = uu;
/*............................................................................*/
/* Hilbert [ spectral ] norm: */

            jac.h.r[ii][jj] = ZERO;
            jac.h.i[ii][jj] = ZERO;

            kk = null; do
            {
               jac.h.r[ii][jj] += ( mtx1[ii][kk]*mtx1[jj][kk] );
            } while (( ++kk ) < DIMNS ); 
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );

      ( jev->rank ) = DIMNS;
/*............................................................................*/
      jev = jacobi( jev );     /*  returns the Hilbert [spectral] norm of NH  */
/*...........................*//*                                             */
      if ( jev == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout, "\n Error on Hilbert norm computation "
            "in matrix NH !!!\n ");
         ( rtp->rtn ) = null;
         return rtp;
      };
      if ( SMX_ONE_PLUS < sqrt( jac.hn ))
      {  
         fprintf( stdout, "\n\n Warning from function %s :", __func__ );
         fprintf( stdout, "\n Sup-norm of magnetic mode S-matrix "
            "( block NH ) = %.16e > 1 .", sqrt( jac.hn ) );
         fprintf( stdout, "\n - This may yield unstable DSC process !" );
         fprintf( stdout, "\n [ Should be remedied "
            "by chosing a smaller DSC time step.]\n " );
         ( rtp->rtn ) = null;
         return rtp;
      };           
/*............................................................................*/
      hnorm2 = sqrt( jac.hn );

# if SMX_DSPINTM == 1
      fprintf( stdout, "\n Hilbert norm: jac.hn = %.16le   \n", hnorm2 );
# endif
      if ( SMX_ONE_PLUS < hnorm2 )
      {
         fprintf( stdout, "\n\n Warning from function %s :", __func__ );
         fprintf( stdout, "\n Sup-norm of magnetic mode S-matrix "
            "equals %.16e > 1.", hnorm2 );
         fprintf( stdout, "\n - This may yield unstable DSC process !" );
         fprintf( stdout, "\n [ Should be remedied "
            "by chosing a smaller DSC time step.]" );
         fprintf( stdout,
            "\n\n Please acknowledge ! [ enter any character ]:" );
         scanf( "%s", ptr );
         ( rtp->rtn ) = null;
         return rtp;
      };
/*............................................................................*/
      if( fabs( rtp->fr ) < SMX_TRIVBND )
      {
         ( rtp->fr ) = ZERO;
         ( rtp->dp ) = ZERO;

         ii = null; do
         {
            jj = null; do
            {
               ( rtp->ser[ii][jj] ) = ZERO;
               ( rtp->sei[ii][jj] ) = ZERO;
               ( rtp->smr[ii][jj] ) = ZERO;
               ( rtp->smi[ii][jj] ) = ZERO;
            } while (( ++jj ) < DIMNS );
         } while (( ++ii ) < DIMNS );
         
         return rtp;
      };
/*------------------------------------------------------------------------------

      Transformation of S-parameters into the frequency domain

------------------------------------------------------------------------------*/
      bound1 = SMX_PRECISION;

      ( rtp->omega ) = 2.*PI*( rtp->fr );
      ( rtp->dp ) = ( rtp->omega ) * ( rtp->dt );

      xx = cos( rtp->dp );
      yy = sin( rtp->dp );

      ii = null; do
      {
         jj = null; do
         {
            gss.mr[ii][jj] = ( double ) ( ii == jj ) * xx;
            gss.mr[ii][jj] -= ( rtp->se[ii+DIMNS][jj+DIMNS] );
            gss.mi[ii][jj] = ( double ) ( ii == jj )*yy;
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: ( exp(*) - NE ) |-> ( exp(*) - NE )^-1 */

      ( gjp->rank ) = DIMNS; /* the rank */
      ( gjp->neqs ) = DIMNS; 
      ( gjp->opt ) = 'i';    /* inversion */

# if SMX_FULPIVT == 1
      gjp = gssjpv( gjp );
# else
      gjp = gssjrd( gjp );
# endif
/*............................................................................*/
      if ( gjp == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout,
            "\n Error on matrix inversion ( exp(*) - NE )^-1 !!!\n " );
         ( rtp->rtn ) = null;
         return rtp;
      };
/*............................................................................*/
      ii = null; do
      {
         jj = null; do
         {
            ( cpt->r ) = ZERO;
            ( cpt->i ) = ZERO;

/* matrix (( exp(*) - NE )^-1 * ME ): */

            kk = null; do
            {
               CPRODUCT( gss.zr[ii][kk], gss.zi[ii][kk],
                         ( rtp->se[kk+DIMNS][jj] ), 0., uu, vv );

               ( cpt->r ) += uu;
               ( cpt->i ) += vv;

            } while (( ++kk ) < DIMNS );
            ( cpt->r ) += ( rtp->se[ii][jj] ); /* add matrix KE */ 

            CPRODUCT( xx, -yy, ( cpt->r ), ( cpt->i ), uu, vv );

            if( fabs( uu ) < bound1 )
               ( rtp->ser[ii][jj] ) = ZERO;
            else
               ( rtp->ser[ii][jj] ) = uu;

            if( fabs( vv ) < bound1 )
               ( rtp->sei[ii][jj] ) = ZERO;
            else
               ( rtp->sei[ii][jj] ) = vv;

         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/* 
      S-matrix SE[][] ready: 
      SE[][] = (( rtp->ser[][] ) + j*( rtp->sei[][] ))
*/
/*............................................................................*/
      ii = null; do
      {
         jj = null; do
         {
            gss.mr[ii][jj] = ( double ) ( ii == jj ) * xx;
            gss.mr[ii][jj] -= ( rtp->sm[ii+DIMNS][jj+DIMNS] );
            gss.mi[ii][jj] = ( double ) ( ii == jj )*yy;
         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
/* matrix inversion: ( exp(*) - NH ) |-> ( exp(*) - NH )^-1 : */

      ( gjp->rank ) = DIMNS; /* the rank */
      ( gjp->neqs ) = DIMNS; 
      ( gjp->opt ) = 'i';    /* inversion */

# if SMX_FULPIVT == 1
      gjp = gssjpv( gjp );
# else
      gjp = gssjrd( gjp );
# endif
/*............................................................................*/
      if ( gjp == NULL )
      {
         fprintf( stdout, "\n\n Message from function %s :", __func__ );
         fprintf( stdout,
            "\n Error on matrix inversion ( exp(*) - NH )^-1 !!!\n " );
         ( rtp->rtn ) = null;
         return rtp;
      };
/*............................................................................*/
      ii = null; do
      {
         jj = null; do
         {
            ( cpt->r ) = ZERO;
            ( cpt->i ) = ZERO;

/* matrix (( exp(*) - NH )^-1 * MH ): */

            kk = null; do       
            {
               CPRODUCT( gss.zr[ii][kk], gss.zi[ii][kk],
                         ( rtp->sm[kk+DIMNS][jj] ), 0., uu, vv );

               ( cpt->r ) += uu;
               ( cpt->i ) += vv;

            } while (( ++kk ) < DIMNS );
            ( cpt->r ) += ( rtp->sm[ii][jj] ); /* add matrix KH */

            CPRODUCT( xx, -yy, ( cpt->r ), ( cpt->i ), uu, vv );

            if( fabs( uu ) < bound1 )
               ( rtp->smr[ii][jj] ) = ZERO;
            else
               ( rtp->smr[ii][jj] ) = uu;

            if( fabs( vv ) < bound1 )
               ( rtp->smi[ii][jj] ) = ZERO;
            else
               ( rtp->smi[ii][jj] ) = vv;

         } while (( ++jj ) < DIMNS );
      } while (( ++ii ) < DIMNS );
/*  
      S-matrix SH[][] ready:
      SH[][] = (( rtp->smr[][] ) + j*( rtp->smi[][] ))
*/
/*............................................................................*/
/* symmetry check: */

      bound1 = SMX_DIAGBND;

      strcpy( rtp->etyp, dgn_e );
      strcpy( rtp->mtyp, dgn_m );

      ii = ONE; do
      {
         jj = null; do 
         {
            if ( *( rtp->etyp ) == 'd' )
            {
               uu = ( fabs( rtp->ser[ii][jj] ) + fabs( rtp->ser[jj][ii] ));
               vv = ( fabs( rtp->sei[ii][jj] ) + fabs( rtp->sei[jj][ii] ));

               if (( bound1 < uu )
                 ||( bound1 < vv ))
                  strcpy(( rtp->etyp ), sym_e );
            };
            if ( *rtp->etyp == 's' )
            {
               uu = fabs(( rtp->ser[ii][jj] ) - ( rtp->ser[jj][ii] ));
               vv = fabs(( rtp->sei[ii][jj] ) - ( rtp->sei[jj][ii] ));

               if (( bound1 < uu )
                 ||( bound1 < vv ))
                  strcpy(( rtp->etyp ), asy_e );
            };
            if ( *( rtp->mtyp ) == 'd' )
            {
               uu = ( fabs( rtp->smr[ii][jj] ) + fabs( rtp->smr[jj][ii] ));
               vv = ( fabs( rtp->smi[ii][jj] ) + fabs( rtp->smi[jj][ii] ));

               if (( bound1 < uu )
                 ||( bound1 < vv ))
                  strcpy(( rtp->mtyp ), sym_m );
            };
            if ( *rtp->mtyp == 's' )
            {
               uu = fabs(( rtp->smr[ii][jj] ) - ( rtp->smr[jj][ii] ));
               vv = fabs(( rtp->smi[ii][jj] ) - ( rtp->smi[jj][ii] ));

               if (( bound1 < uu )
                 ||( bound1 < vv ))
                  strcpy(( rtp->mtyp ), asy_m );
            };
         } while (( ++jj ) < ii );
      } while (( ++ii ) < DIMNS );
/*............................................................................*/
# if SMX_UNICHCK == 1
/* unitarity check [ only in the lossless case, loss = null ]: */

      if (( rtp->loss ) == null ) /* loss = null: lossless mesh cell */
      {
         nrm1 = ZERO;
         nrm2 = ZERO;

         bound1 = SMX_UNITBND;
         bound1 *= bound1;

         bound2 = 1. + 3.*SMX_PRECISION;

         ii = null; do
         {
            jj = null; do
            {
               uu = ZERO;
               vv = ZERO;

               kk = null; do
               {
                  CPRODUCT(( rtp->ser[ii][kk]), ( rtp->sei[ii][kk] ),
                     ( rtp->ser[jj][kk] ), - ( rtp->sei[jj][kk] ), rr, ss );
                  uu += rr;
                  vv += ss;
               } while (( ++kk ) < DIMNS );

               if ( ii != jj )
               {
                  rr = uu*uu + vv*vv;
                  if ( bound1 < fabs( rr ) )
                  {
                     fprintf( stdout,
                        "\n\n Error message from function %s :", __func__ );
                     fprintf( stdout, "\n Non - unitary S-matrix !!!\n " );
                     ( rtp->rtn ) = null;
                     return rtp;
                  };
               }
               else /* if ( ii == jj ) */
               {
                  rr = uu*uu + vv*vv;

                  if ( nrm1 < rr )
                     nrm1 = rr;

                  uu -= 1.;
                  rr = uu*uu + vv*vv;
                  if ( bound1 < fabs( rr ) )
                  {
                     fprintf( stdout,
                        "\n\n Error message from function %s :", __func__ );
                     fprintf( stdout, "\n Non - unitary S-matrix !!!\n " );
                     ( rtp->rtn ) = null;
                     return rtp;
                  };
               };

               uu = ZERO;
               vv = ZERO;

               kk = null; do
               {
                  CPRODUCT(( rtp->smr[ii][kk] ), ( rtp->smi[ii][kk] ),
                     ( rtp->smr[jj][kk] ), - ( rtp->smi[jj][kk] ), rr, ss );
                  uu += rr;
                  vv += ss;
               } while (( ++kk ) < DIMNS );
     
               if ( ii != jj )
               {
                  rr = uu*uu + vv*vv;
                  if ( bound1 < fabs( rr ) )
                  {
                     fprintf( stdout,
                        "\n\n Error message from function %s :", __func__ );
                     fprintf( stdout, "\n Non - unitary S-matrix !!!\n " );
                     ( rtp->rtn ) = null;
                     return rtp;
                  };
               }
               else /* if ( ii == jj ) */
               {
                  rr = uu*uu + vv*vv;
                  if ( nrm2 < rr )
                     nrm2 = rr;

                  uu -= 1.;
                  rr = uu*uu + vv*vv;
                  if ( bound1 < fabs( rr ) )
                  {
                     fprintf( stdout,
                        "\n\n Error message from function %s :", __func__ );
                     fprintf( stdout, "\n Non - unitary S-matrix !!!\n " );
                     ( rtp->rtn ) = null;
                     return rtp;
                  }; 
               };
            } while (( ++jj ) < DIMNS );
         } while  (( ++ii ) < DIMNS );

         if ( bound2 < nrm1 )
         {
            nrm1 = sqrt( sqrt( nrm1 ));

            ii = null; do
            {
               jj = null; do
               {
                  ( rtp->ser[ii][jj] ) /= nrm1;
                  ( rtp->sei[ii][jj] ) /= nrm1;
               } while (( ++jj ) < DIMNS );
            } while (( ++ii ) < DIMNS );
         };

         if ( bound2 < nrm2 )
         {
            nrm2 = sqrt( sqrt( nrm2 ));

            ii = null; do
            {
               jj = null; do
               {
                  ( rtp->smr[ii][jj] ) /= nrm2;
                  ( rtp->smi[ii][jj] ) /= nrm2;
               } while (( ++jj ) < DIMNS );
            } while (( ++ii ) < DIMNS );
         };
      }; /* end if ( rtp->loss ) == null: unitarity check for lossless cell */
# endif
/*............................................................................*/
      ( rtp->rtn ) = ONE;
      return rtp;

   }; /* end of option 's': nodal s_parameters */

   ( rtp->rtn ) = null;
   return rtp; /* Error: statement never reached in any legal option */ 
} 
# if defined ( OPTIMIZE )
   # pragma OPTIMIZE OFF
# endif
/*----------------------------------------------------------------------------*/
# if defined ( SMX_CLFCES )
   # undef SMX_CLFCES
# endif
# if defined ( SMX_CLSKEW )
   # undef SMX_CLSKEW
# endif
# if defined ( SMX_SKWVOL )
   # undef SMX_SKWVOL
# endif
# undef SMX_PORTS
# undef SMX_EDGES
# undef SMX_ORTNORM
# undef SMX_VOLUME
/*----------------------------------------------------------------------------*/
# undef SMX_MEDCHCK
# undef SMX_CLLCRDS
# undef SMX_ISOCHCK
# undef SMX_UNICHCK
# undef SMX_DSPINTM
# undef SMX_FULPIVT 
# undef SMX_FCEMODE
# undef SMX_DIAGBND
# undef SMX_UNITBND
# undef SMX_QRT_MINUS
# undef SMX_ONE_PLUS
# undef SMX_ONE_MINUS
# undef SMX_TRIVBND
# undef SMX_ECRRBND
# undef SMX_MCRRBND
# undef SMX_PRECISION
# undef SMX_SQRPRECIS
# undef SMX_INCLUDE
/*----------------------------------------------------------------------------*/
# undef CPRODUCT
# undef EPS_VAC
# undef MY_VAC_
# undef ELECTR_CHRGE
# undef ELECTR_MASS_
# undef BOHRs_MAGNTN
# undef PLANCKs_CNST
/*----------------------------------------------------------------------------*/
# undef SMX_DEBUG /* should be removed after program check */
/************* end of DSC S-matrix generation function smatrx(*) **************/
