/* [ file: gradnt.c ] */
# define DO_GRADNT "gradnt(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function gradnt(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   This macro returns the nodal gradient, as well as the node vector          *
*   finite differences [representing the gradient in canonical coordinates]    *
*   of any quantity within a CLUSTER structure, pointed to by cls *CLUSTER,    *
*   hh being the index the pertinent (heat&fluid) s-parameter set.             *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include <string.h>
# include <stdio.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# include "../CONFIG.H"
# include "./SOLVER.CONF"
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*============================================================================*/

CLUSTER *gradnt( CLUSTER *cls, long hh )
{
/* declarations: */

   register signed char
      kk = null,
      ll = null;

   static CLUSTER
     *rtp = NULL;
/*----------------------------------------------------------------------------*/
   rtp = cls;
/*............................................................................*/
/* restore nodal gradient, */
/* in cell coordinates [ with respect to basis b[] ] */
/* df[k] = < b[k] | grad(*) >, */
/* and in canonical coordinates */
/* df[k] = < e[k] | grad(*) > = adj(B^-1 )*db ]: */
/*............................................................................*/
   kk = null; do
   {
      ( rtp->db[kk] ) = ( cls->face[TWO*kk+ONE] ) - ( cls->face[TWO*kk] );
   } while(( ++kk ) < THREE );

   kk = null; do
   {
      ( rtp->grad[kk] ) = ZERO;
      ll = null; do
      {
         ( rtp->grad[kk] ) += (( cls->hsp->bi[hh][kk][ll] )*( rtp->db[ll] ));
      } while(( ++ll ) < THREE );
   } while(( ++kk ) < THREE );

   ( rtp->hh ) = hh; /* remember pertinent s-parameter set index */
/*............................................................................*/
   return rtp;
} 
/*============================================================================*/
/************************* end of function gradnt(*) **************************/

/*******************************************************************************
*                                                                              *
*   ANSI C function gradfc(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   This function returns the face gradients of any quantity in a CLUSTER      *
*   structure pointed to by cls *CLUSTER, hh being the index of the perti-     *
*   nent (heat & fluid) s-parameter set.                                       *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: 23 March 2011 ]                       <contact@steffen-hein.org> *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

CLUSTER *gradfc( CLUSTER *cls, long hh )
{
/* declarations: */

   register signed char
      fc = null,
      kk = null,
      ll = null,
      sgn = ONE;

   static CLUSTER
     *rtp = NULL;
/*----------------------------------------------------------------------------*/
   rtp = cls;
/*............................................................................*/
/* restore face temperature gradients */
/* [ dtdf[k] = adj(B^-1)*dtdb ] */
/*............................................................................*/
   sgn = ONE;

   ( rtp->dvgr ) = ZERO;
   fc = null; do
   {
      kk = null; do
      {
         ll = ( short )( fc/TWO );
         ( rtp->grdf[fc][kk] ) = 2.*sgn*\
            ( cls->hsp->bi[hh][kk][ll] )*(( cls->node ) - ( cls->face[fc] ));

         ++ll; ll%=THREE; /* ll := ( ll+ONE ) mod THREE */
         ( rtp->grdf[fc][kk] ) +=\
            (( cls->hsp->bi[hh][kk][ll] )*( cls->db[ll] ));
         ++ll; ll%=THREE;
         ( rtp->grdf[fc][kk] ) +=\
            (( cls->hsp->bi[hh][kk][ll] )*( cls->db[ll] ));
/* div grad: */
         ( rtp->dvgr ) += (( rtp->grdf[fc][kk] )*( cls->hsp->f[hh][fc][kk] ));
      } while(( ++kk ) < THREE );
      sgn *= ( -ONE );
   } while(( ++fc ) < FACES );

   ( rtp->hh ) = hh; /* remember pertinent s-parameter set index */
/*............................................................................*/
   return rtp;
} 
/*============================================================================*/
# endif /* DSC_HCRMDE != 0 */
/************************* end of function gradfc(*) **************************/
/************************** end of file gradnt.c ******************************/
