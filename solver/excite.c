/* [ file: excite.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function excfld(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   TLM/DSC Maxwell field excitation function                                  *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPT_LEVEL 2
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize solver configuration header SOLVER.CONF */
# include "SOLVER.CONF"
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
/* the following [two] macros should be defined in "../CONFIG.H"              */
# ifndef DSC_ADJCLL      /* assign neighbouring cell index top.mn[][k] thus:  */
   # define DSC_ADJCLL 0 /* 0: to faces [ k is a face index; 0 <= k < 6 ]     */
# endif                  /* 1: to ports [ k is a port index; 0 <= k < 12 ]    */
/*----------------------------------------------------------------------------*/
# ifndef DSC_FCELBL      /* label neighbouring face index top.fn[][] thus:    */
   # define DSC_FCELBL 0 /* 0: unsigned, started with index null              */
# endif                   /* 1: unsigned, started with index ONE              */
                          /* 2: signed, started with index ONE                */
                          /* [ the negative/positive sign indicating opposed/ */
                          /*   aligned orientation to neighbouring cell face ]*/
/*----------------------------------------------------------------------------*/
/* There are different ways to interlace excited with incident quantities:    */
     
# ifndef EXC_IMPOSE
   # define EXC_IMPOSE 0 /* 0: additive superposition [ recommended ]         */
# endif                  /* 1: direct imposition                              */
                         /* 2: imposition with scattered quantities subtracted*/
                         /*    [ yields fixed total quantities ].             */
# ifndef EXC_MAXFLD
   # define EXC_MAXFLD 1 /* 1: write maximum excited field into ( ept->mx )   */
# endif
/*----------------------------------------------------------------------------*/
/* constants, stabilizing ONE [ to be modified with care ]:                   */

# ifndef PRECISION
   # define PRECISION ( double )( 1.000e-15 )
# endif
/*----------------------------------------------------------------------------*/
# if EXC_IMPOSE != 0
   # define ONE_MINUS ( 1. - 33.*PRECISION ) 
# endif
/*============================================================================*/

DSC_FIELDS *\
excfld( struct solverstat *state )
{
/* allusions: [void] */
/* declarations: */
/*............................................................................*/
# if DSC_DOMAIN != 0
   FILE 
     *display = stdout;
# endif
/*............................................................................*/
   static long
      mm = null;

   static short 
      ii = null;

   static signed char
      pp = null;
/*
   static char
      ptr[STS_SIZE] = {null};
*/
   static struct
      excitation *ept;

   static DSC_FIELDS
      *ret = NULL;

/* prototypes: */

   char 
     *lotos ( long lngint, char length );
/*----------------------------------------------------------------------------*/
   if(( state->dmn ) == 't' )
   {
/*............................................................................*/
# if DSC_DOMAIN == 2
      fprintf( display, "\n\n Error message from function %s :",
         __func__ );
      fprintf( display, "\n Solver is compiled for frequency "
         "domain operation 'DSC_DOMAIN = 2'" );
      fprintf( display, "\n - yet run in time domain mode 1 !!!" );
      fprintf( display, "\n [ Change time/frequency domain macro "
         "DSC_DOMAIN in SOLVER.CONF from 2 to" );
      fprintf( display, "\n   1: time domain, or          "
         "\n   0: time & frequency domain. " );
      fprintf( display, "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
         "in_mode_DSC_DOMAIN=2] ", LGS_SIZE );
      return null; /* abnormal return */
/*............................................................................*/
# else /* if DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
      ept = ( state->ept );
      ret = ( state->inc );

      ii = null;
      while( ii < ( ept->ne ))
      {
         mm = ( ept->me[ii] );
         pp = ( ept->pe[ii] ) - ONE;
/*............................................................................*/
# if EXC_IMPOSE == 1
         ( ret->r[mm][pp] ) = ZERO;
# elif EXC_IMPOSE == 2
         ( ret->r[mm][pp] ) = - ONE_MINUS*( ret->i[mm][pp] );
# endif
/*............................................................................*/
         ( ret->r[mm][pp] ) += ( state->swing )*( ept->er[ii] );

         ii++ ;
      };

      ii = null; 
      while( ii < ( ept->nh ))
      {
         mm = ( ept->mh[ii] );
         pp = ( ept->ph[ii] ) - ONE;
/*............................................................................*/
# if EXC_IMPOSE == 1
         ( ret->r[mm][pp] ) = ZERO;
# elif EXC_IMPOSE == 2
         ( ret->r[mm][pp] ) = ONE_MINUS*( ret->i[mm][pp] );
# endif
/*............................................................................*/
         ( ret->r[mm][pp] ) += ( state->swing )*( ept->hr[ii] );

         ii++ ; 
      };

      return ret; /* normal return */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
   }
   else /* if (( state->dmn ) != 't' ): frequency domain */
   {
/*............................................................................*/
# if DSC_DOMAIN == 1
      fprintf( display, "\n\n Error message from function %s :", __func__);
      fprintf( display, "\n Solver is compiled for time domain "
         "operation, 'DSC_DOMAIN = 1'" );
      fprintf( display, "\n - yet run in frequency domain mode 2 !!!" );
      fprintf( display, "\n [ Change time/frequency domain macro "
         "DSC_DOMAIN in SOLVER.CONF from 1 to" );
      fprintf( display, "\n   2: frequency domain, or     "
         "\n   0: time & frequency domain. " );
      fprintf( display, "\n   Then re-compile and restart program. ]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
         "in_mode_DSC_DOMAIN=1] ", LGS_SIZE );
      return NULL;
/*............................................................................*/
# else /* DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
      ept = ( state->ept );
      ret = ( state->inc );

      ii = null;
      while( ii < ( ept->ne ))
      {
         mm = ( ept->me[ii] );
         pp = ( ept->pe[ii] ) - ONE;
/*............................................................................*/
# if EXC_IMPOSE == 1
         ( ret->r[mm][pp] ) = ZERO;
         ( ret->i[mm][pp] ) = ZERO;
# elif EXC_IMPOSE == 2
         ( ret->r[mm][pp] ) = - ONE_MINUS*( ret->r[mm][pp] );
         ( ret->i[mm][pp] ) = - ONE_MINUS*( ret->i[mm][pp] );
#endif
/*............................................................................*/
         ( ret->r[mm][pp] ) += ( state->swing )*( ept->er[ii] );
         ( ret->i[mm][pp] ) += ( state->swing )*( ept->ei[ii] );

         ii++ ;
      };

      ii = null;
      while( ii < ( ept->nh ))
      {
         mm = ( ept->mh[ii] );
         pp = ( ept->ph[ii] ) - ONE;

# if EXC_IMPOSE == 1
         ( ret->r[mm][pp] ) = ZERO;
         ( ret->i[mm][pp] ) = ZERO;
# elif EXC_IMPOSE == 2
         ( ret->r[mm][pp] ) = ONE_MINUS*( ret->r[mm][pp] );
         ( ret->i[mm][pp] ) = ONE_MINUS*( ret->i[mm][pp] );
# endif
         ( ret->r[mm][pp] ) += ( state->swing )*( ept->hr[ii] );
         ( ret->i[mm][pp] ) += ( state->swing )*( ept->hi[ii] );

         ii++ ;
      };

      return ret; /* normal return */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
   }; /* end if (( state->dmn ) != 't' ), frequency domain */
} 
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
# if EXC_IMPOSE != 0
   # undef ONE_MINUS
# endif
# undef EXC_IMPOSE
/************** end of DSC system excitation function excfld(*) ***************/



/*----------------------------------------------------------------------------*/
/* [ file: exchcr.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function exchcr(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   DSC temperature distribution and heat current excitation function          *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPT_LEVEL 3
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
/*
# include "../math/consts.h"
*/
/*----------------------------------------------------------------------------*/
/* structure excitation exc: cf. function scattr.h                            */
/*----------------------------------------------------------------------------*/
/* There are different ways to interlace excited with incident quantities:    */
/*----------------------------------------------------------------------------*/
/* constants, stabilizing ONE [ to be modified with care ]:                   */

# ifndef PRECISION
   # define PRECISION ( 1.e-15 )
# endif
/*============================================================================*/

DSC_HCRRTS * \
exchcr( struct solverstat *state )
{
/* allusions: */
/*
   extern struct solverstat solver;
   extern struct excitation exc;
*/
/* declarations: */
/*............................................................................*/
# if DSC_DOMAIN != 0
   FILE 
     *display = stdout;
# endif
/*............................................................................*/
   static long
      mm = null;

   static short 
      fc = null,
      ii = null,
      cc = null;
/*
   static char
      ptr[STS_SIZE] = {null};
*/
   static struct
     excitation *ept = NULL;

   static DSC_HCRRTS
      *hre = NULL;

/* prototypes: - */
/*----------------------------------------------------------------------------*/
   if(( state->dmn ) == 't' )
   {
/*............................................................................*/
# if DSC_DOMAIN == 2
      fprintf( display, "\n\n Error message from function %s :", __func__ );
      fprintf( display, "\n Solver is compiled for frequency "
         "domain operation 'DSC_DOMAIN = 2'" );
      fprintf( display, "\n - yet run in time domain mode 1 !!!" );
      fprintf( display, "\n [ Change time/frequency domain macro "
         "DSC_DOMAIN in SOLVER.CONF from 2 to" );
      fprintf( display, "\n   1: time domain, or          "
         "\n   0: time & frequency domain. " );
      fprintf( display, "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
                            "in_mode_DSC_DOMAIN=2] ", LGS_SIZE );
      return NULL; /* abnormal return */
/*............................................................................*/
# else /* if DSC_DOMAIN == 0, 1 */

      ept = ( state->ept );
      cc = ( state->hclbl );
      hre = ( state->hci[cc] );

      ii = null;
      while( ii < ( ept->nhc ))
      {
         mm = ( ept->mhc[ii] );
         fc = ( ept->fhc[ii] );
         ( hre->ic[mm][fc] ) = ( state->hcswg )*( ept->hc[ii] );
         ii++ ;
      };

      ii = null;
      while( ii < ( ept->ntf ))
      {
         mm = ( ept->mtf[ii] );
         fc = ( ept->ftf[ii] );
         ( hre->tf[mm][fc] ) = ( state->hcswg )*( ept->tf[ii] );

         ii++ ;
      };

      ii = null;
      while( ii < ( ept->ntn ))
      {
         mm = ( ept->mtn[ii] );
         ( hre->tn[mm] ) = ( state->hcswg )*( ept->tn[ii] );

         ii++ ;
      };

      return hre; /* normal return */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 1 */
/*............................................................................*/
   }
   else /* if (( state->dmn ) != 't' ): frequency domain */
   {
/*............................................................................*/
# if DSC_DOMAIN == 1
      fprintf( display, "\n\n Error message from function %s :", __func__);
      fprintf( display, "\n Solver is compiled for time domain "
         "operation, 'DSC_DOMAIN = 1'" );
      fprintf( display, "\n - yet run in frequency domain mode 2 !!!" );
      fprintf( display, "\n [ Change time/frequency domain macro "
         "DSC_DOMAIN in SOLVER.CONF from 1 to" );
      fprintf( display, "\n   2: frequency domain, or     "
         "\n   0: time & frequency domain." );
      fprintf( display, "\n   Then re-compile and restart program.]\n " );

      strncpy(( state->fcterr ), __func__, SHS_SIZE );
      strncpy(( state->errmsg ), "domain_!!!_[program_compiled_"
         "in_mode_DSC_DOMAIN=1] ", LGS_SIZE );

      return NULL; /* abnormal return */
/*............................................................................*/
# else /* if DSC_DOMAIN == 0 or 2 */

      ept = ( state->ept );
      cc = ( state->hclbl );
      hre = ( state->hci[cc] );

      ii = null;
      while( ii < ( ept->nhc ))
      {
         mm = ( ept->mhc[ii] );
         fc = ( ept->fhc[ii] );
         ( hre->ic[mm][fc] ) = ( state->hcswg )*( ept->hc[ii] );
         ii++ ;
      };

      ii = null;
      while( ii < ( ept->ntf ))
      {
         mm = ( ept->mtf[ii] );
         fc = ( ept->ftf[ii] );
         ( hre->tf[mm][fc] ) = ( state->hcswg )*( ept->tf[ii] );

         ii++ ;
      };

      ii = null;
      while( ii < ( ept->ntn ))
      {
         mm = ( ept->mtn[ii] );
         ( hre->tn[mm] ) = ( state->hcswg )*( ept->tn[ii] );

         ii++ ;
      };

      return hre; /* normal return */
/*............................................................................*/
# endif /* DSC_DOMAIN == 0 or 2 */
/*............................................................................*/
   }; /* end if (( state->dmn ) != 't' ), frequency domain */
}
/*============================================================================*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif
/******** end of DSC heat current and temp distrib function exchcr(*) *********/
# endif /* DSC_HCRMDE != 0 */
/********* end of DSC system excitation functions excfld.c, exchcr.c **********/
