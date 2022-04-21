/* [ file: scattr.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function scattr(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC scattering algorithm based on non-orthogonal condensed node            *
*   with tensor material parameters and gyrotropic structure.                  *
*   [ cf. 'Finite Difference Time-Domain Approximation of Maxwell's            *
*   Equations with Non-Orthogonal Condensed TLM Mesh', International           *
*   Journal of Numerical Modelling, vol.7, 179-188 (1994) and 'TLM             *
*   Numerical Solution of Bloch's Equations for Magnetized Gyrotropic          *
*   Media', Applied Mathematical Modelling, vol.22, 221-229 (1997).]           *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
/*----------------------------------------------------------------------------*/
# if defined ( OPTIMIZE )
   # pragma OPT_LEVEL 3
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize solver configuration header SOLVER.CONF */
# include "SOLVER.CONF"
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
/* the configuration header */
# include "SCATFLD.CONF"
/*============================================================================*/
/* the function body */
# include "scatfld-1.0r3.a"
/*============================================================================*/
/************* end of Maxwell field updating function scatfld(*) **************/





/*******************************************************************************
*                                                                              *
*   ANSI C function scathcr(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC nodal temperatures updating function                                   *
*   [ Stable version, fluid flows separately updated in function sctflow(*) ]  *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: 23 March 2011 ]                       <contact@steffen-hein.org> *
*                                                                              *
*******************************************************************************/
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
/* the configuration header */
# include "SCATHCR.CONF"
/*============================================================================*/
/* the function body */
# include "scathcr-1.0r3.a"
/*============================================================================*/
# endif /* DSC_HCRMDE != 0 */
/***************** end of thermal updating function scathcr(*) ****************/





/*******************************************************************************
*                                                                              *
*   ANSI C function sctflow(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC nodal fluid flow updating function                                     *
*   [ Stable version, diffusive temperature updates in function scathcr(*) ]   *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: 23 March 2011 ]                       <contact@steffen-hein.org> *
*                                                                              *
*******************************************************************************/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE != 0
/*----------------------------------------------------------------------------*/
/* the configuration header */
# include "SCTFLOW.CONF"
/*============================================================================*/
/* the function body */
# include "sctflow-1.0r3.a"
/*============================================================================*/
/************** end of fluid flow updating function sctflow(*) ****************/





/*******************************************************************************
*                                                                              *
*   Function fluidstate(*)                                                     *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   Returns nodal pressure as a function of temperature and density            *
*   [ Boyle-Mariotte law of ideal gas ]                                       *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: 23 March 2011 ]                       <contact@steffen-hein.org> *
*                                                                              *
*******************************************************************************/
# if (( CMPRSSBL != 0 )\
    &&( BSQAPRX == 0 ))
/*============================================================================*/
# include "fluidstate.b"
/*============================================================================*/
# endif /* (( CMPRSSBL != 0 )\
          &&( BSQAPRX == 0 )) */
/********************** end of function fluidstate(*) *************************/
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*************************** end of file scattr.c *****************************/
