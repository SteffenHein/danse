/* [ file: grdnts.c ] */
# define DO_GRDNTS "grdnts(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function gradnt(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0r3                                                                *
*                                                                              *
*   This function returns the nodal gradient, as well as the node vector       *
*   finite differences [which represent the gradient in nodal coordinates],    *
*   of any quantity in a CLUSTER structure pointed to by cls *CLUSTER;         *
*   binvrs denotes the inverse adjoint node vector matrix.                     *
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
# include "SOLVER.CONF"
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*============================================================================*/

CLUSTER *grdnde( CLUSTER *cls, double binvrs[THREE][THREE] )
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
/* restore node temperature gradient, */
/* in cell coordinates [ with respect to basis b[] ] */
/* dtdb[k] = < b[k] | grad(T) >, */
/* and in canonical coordinates */
/* dt[k] = < e[k] | grad(T) > = adj(B^-1 )*dtdb ]: */
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
         ( rtp->grad[kk] ) += (( binvrs[kk][ll] )*( rtp->db[ll] ));
      } while(( ++ll ) < THREE );
   } while(( ++kk ) < THREE );
/*............................................................................*/
   return rtp;
} 
/*============================================================================*/
/************************* end of function grdnde(*) **************************/

/*******************************************************************************
*                                                                              *
*   ANSI C function grdfce(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0r3                                                              *
*                                                                              *
*   This function returns the face gradients of any quantity in a CLUSTER      *
*   structure pointed to by cls *CLUSTER, binvrs being the inverse adjoint     *
*   node vector matrix and face the face vector coordinates.                   *
*   nent (heat & fluid) s-parameter set.                                       *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: 23 March 2011 ]                       <contact@steffen-hein.org> *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

CLUSTER *grdfce( CLUSTER *cls, double binvrs[THREE][THREE], \
   double face[SIX][THREE] )
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
/*   dtdf[k] = adj(B^-1)*dtdb */
/*............................................................................*/
   sgn = ONE;
   ( rtp->dvgr ) = ZERO;
   fc = null; do
   {
      kk = null; do
      {
         ll = ( short )( fc/TWO );
         ( rtp->grdf[fc][kk] ) = 2.*sgn*\
            ( binvrs[kk][ll] )*(( cls->node ) - ( cls->face[fc] ));

         ++ll; ll%=THREE; /* ll := ( ll+ONE ) mod THREE */
         ( rtp->grdf[fc][kk] ) +=\
            (( binvrs[kk][ll] )*( cls->db[ll] ));
         ++ll; ll%=THREE;
         ( rtp->grdf[fc][kk] ) +=\
            (( binvrs[kk][ll] )*( cls->db[ll] ));
/* div grad: */
         ( rtp->dvgr ) += (( rtp->grdf[fc][kk] )*( face[fc][kk] ));
      } while(( ++kk ) < THREE );
      sgn *= ( -ONE );
   } while(( ++fc ) < FACES );
/*............................................................................*/
   return rtp;
} 
/*============================================================================*/
# endif /* DSC_HCRMDE != 0 */
/************************* end of function grdfce(*) **************************/
/************************** end of file grdnts.c ******************************/
