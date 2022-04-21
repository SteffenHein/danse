/* [ file: spline.c ] */
# define DO_SPLINE "spline(*)"
/*******************************************************************************
*                                                                              *
*   ANSI-C function spline(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   Cubical spline interpolation routine                                       *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <time.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H"
/*----------------------------------------------------------------------------*/
/* Edit and customize this header for POSTER configuration: */
# include "../poster/POSTER.CONF"
/*----------------------------------------------------------------------------*/
/*
# include "../math/SPL.CONF"
*/
/*----------------------------------------------------------------------------*/
# ifndef SPL_SPNTS
   # define SPL_SPNTS  500  /* maximum number of support points               */
# endif
# ifndef SPL_INTPL
   # define SPL_INTPL 5000  /* maximum number of interpolation points         */
# endif
# ifndef SPL_MAXIT
   # define SPL_MAXIT 2000  /* maximum number of iterations                   */
# endif
# ifndef SPL_BOUND
   # define SPL_BOUND ( 1.e-07 )
# endif
/*----------------------------------------------------------------------------*/
typedef struct
{
   signed char
      rtn;

   long
      mm, /* number of support points */
      nn; /* number of interpolated points */

   double
      intgr, fmin, fmax,
      vct[SPL_SPNTS][TWO], /* given points */
      dmn[SPL_INTPL],      /* interpolated points */
      fct[SPL_INTPL],
      drv[SPL_INTPL];

} SPLINES;
/*============================================================================*/

SPLINES *
spline( SPLINES *spt )
{
/* allusions: */

/* declarations: */

   static double
      ss[SPL_SPNTS] = {ZERO},
      gg[SPL_SPNTS] = {ZERO},
      ww[SPL_SPNTS] = {ZERO};

   static double
      xi    = ZERO,
      xim1  = ZERO,
      xip1  = ZERO,
      xd    = ZERO,
      x1    = ZERO,
      yi    = ZERO,
      yim1  = ZERO,
      yip1  = ZERO,
      hh    = ZERO,
      s1    = ZERO,    
      tt    = ZERO,
      tw    = ZERO,
      uu    = ZERO,
      vv    = ZERO,
      bound = SPL_BOUND,
      newbd = ZERO;
     
   char *ptr;

   char **endp = null;

   long 
      m_ = null,
      ii = null,
      jj = null;

/*============================================================================*/

   ptr = ( char *) calloc ( SHS_SIZE, ONE );

   if ( SPL_SPNTS <= (spt->mm) )
   {
      fprintf( stderr, "\n\n Error in function %s:", __func__ );
      fprintf( stderr, "\n\n Too many support points transferred "
         "in ( spt->vct[][] ) !!!" );
      fprintf( stderr, "\n [ The maximum number is %d = macro SPL_SPNTS "
         "in function %s.", SPL_SPNTS, __func__ );
      fprintf( stderr, "\n - Change macro only in compliance "
         "with memory resources.]\n " );

      exit( EXIT_FAILURE );
   };

   if ( SPL_INTPL <= (spt->nn) )
   {
      fprintf( stderr, "\n\n Error in function %s:", __func__ );
      fprintf( stderr, "\n\n Too many interpolation points requested !!!" );
      fprintf( stderr, "\n [ The maximum number is %d = macro SPL_INTPL "
         "in function %s.", SPL_INTPL, __func__ );
      fprintf( stderr, "\n - Change macro only in compliance "
         "with memory resources.]\n " );

      exit( EXIT_FAILURE );
   };

   m_ = (spt->mm) - ONE;
   ( spt->fmin ) = 1.e+277;
   ( spt->fmax ) = -( spt->fmin );

   if (( spt->vct[null][ONE] ) < ( spt->fmin ))
       ( spt->fmin ) = ( spt->vct[null][ONE] );
   
   if (( spt->vct[m_][ONE] ) < ( spt->fmin ))
       ( spt->fmin ) = ( spt->vct[m_][ONE] );

   if (( spt->fmax ) < ( spt->vct[null][ONE] ))
       ( spt->fmax ) = ( spt->vct[null][ONE] );

   if (( spt->fmax ) < ( spt->vct[m_][ONE] ))
       ( spt->fmax ) = ( spt->vct[m_][ONE] );  

   ii = ONE;
   while ( ii < m_ )
   {
      xi   = ( spt->vct[ii][null] );
      xim1 = ( spt->vct[ii-ONE][null] );
      xip1 = ( spt->vct[ii+ONE][null] );
      yi   = ( spt->vct[ii][ONE] );
      yim1 = ( spt->vct[ii-ONE][ONE] );
      yip1 = ( spt->vct[ii+ONE][ONE] );
      xd   = xi - xim1;
      hh   = xip1 - xim1;

      ww[ii] = xd/hh/2.;
      tt = (( yip1 - yi )/( xip1 - xi ) - ( yi - yim1 ) / xd ) / hh;
      ss[ii] = 2.*tt;
      gg[ii] = 3.*tt;
      ii++ ;
   };

   ss[null] = ZERO;
   ss[(( spt->mm )-ONE )] = ZERO;

   tw = 8. - 4.*sqrt( 3. );

   jj = null;

  uu_0:
   
   uu = ZERO;
   jj++;

   if ( SPL_MAXIT < jj )
   {
      fprintf( stdout, "\n\n Message from function %s: ", __func__ );
      fprintf( stdout, "\n\n Too many iterations !!!" );
      fprintf( stdout, "\n [ You may try to proceed by entering a new bound"
              " > %.7e ] ...: ", bound );
      scanf( "%s", ptr );
      newbd = strtod( ptr, endp );

      if ( bound < newbd )
      {
         bound = newbd;
         jj = null; 
      }
      else
      {
         fprintf( stderr, "\n + program stop + program stop + program stop "
            "+ program stop + program stop +\n " );
         exit( EXIT_FAILURE );
      };
   };

   ii = ONE;
   while ( ii < m_ )
   {
      tt = tw*(-ss[ii] - ww[ii]*ss[ii-ONE] - (.5 - ww[ii])*ss[ii+ONE]+gg[ii] );
      hh = fabs( tt );

      if ( uu < hh )
         uu = hh;

      ss[ii] += tt;
      ii++ ;
   };

   if ( bound < uu )
      goto uu_0;

/*............................................................................*/

   ii = null; 
   while ( ii < m_ )
   {
      gg[ii] = ss[ii+ONE] - ss[ii];
      vv = (( spt->vct[ii+ONE][null] ) - ( spt->vct[ii][null] ));

      if ( 1.e-277 < fabs( vv ))
         gg[ii] /= vv;
      else
      {
         fprintf( stdout, "\n\n Error message from function %s:", __func__ );
         fprintf( stdout, "\n\n division by ZERO !!!" );
         fprintf( stdout,
            "\n [ illegal equal arguments vct[%ld][0], vct[%ld][0] ]",
            ( long ) ( ii+ONE ), ( long ) ii );
         ( spt->rtn ) = null;

         return spt;
      };
         
      ii++ ;
   };

   if (( spt->nn ) == null )
      goto integral;

   jj = null; 
   while ( jj < ( spt->nn ))
   {
      tt = ( spt->dmn[jj] );

      if (( tt < ( spt->vct[null][null] ))||(( spt->vct[m_][null] ) < tt ))
      {
        error:       

         fprintf( stderr, "\n\n Error message from function %s:", __func__ );
         fprintf( stderr, "\n\n Arguments out of bounds !!!" );
         fprintf( stderr, "\n [ x0(vct) = %.15e, x%ld(vct) = %.15e,", 
            ( spt->vct[null][null] ), ( spt->mm ), ( spt->vct[m_][null] ));
         fprintf( stderr,
            "\n   x%ld(dmn) = %.15e ] \n\n", jj, ( spt->dmn[jj] ));

         exit( EXIT_FAILURE );
      };
     
      ii = null; do
      {
         if ((( spt->mm )-ONE ) < ( ++ii ))
            goto error;

      } while(( spt->vct[ii][null] ) < tt );
      ii-- ;

      hh  = ( spt->dmn[jj] ) - ( spt->vct[ii][null] );
      tt  = ( spt->dmn[jj] ) - ( spt->vct[ii+ONE][null] );

      x1  = hh*tt;
      s1  = ss[ii] + hh*gg[ii];
      vv  = 1./6.;
      uu  = vv*( ss[ii]+ss[ii+ONE] + s1 );

      tw = (( spt->vct[ii+ONE][null] ) - ( spt->vct[ii][null] ));
      if ( 1.e-277 < fabs( tw ))
      { 
         tw = ( 1./ tw );
         tw *= ( spt->vct[ii+ONE][ONE] ) - ( spt->vct[ii][ONE] );
      }
      else
      {
         fprintf( stdout, "\n\n Error message from function %s:", __func__ );
         fprintf( stdout, "\n\n division by ZERO !!!" );
         fprintf( stdout,
            "\n [ illegal equal arguments vct[%ld][0], vct[%ld][0] ]",
            ( long ) ( ii+ONE ), ( long ) ii );
         ( spt->rtn ) = null;

         return spt;
      };

      ( spt->fct[jj] ) = tw*hh + ( spt->vct[ii][ONE] ) + x1*uu;
      
      if (( spt->fct[jj] ) < ( spt->fmin ))
         ( spt->fmin ) = ( spt->fct[jj] );

      if (( spt->fmax ) < ( spt->fct[jj] ))
          ( spt->fmax ) = ( spt->fct[jj] );

      ( spt->drv[jj] ) = tw + ( hh + tt )*uu + vv*x1*gg[ii];

      jj++ ;
   }; /* while( jj < ( spt->nn )) */

  integral:

   ( spt->intgr ) = ZERO;

   ii = null;
   while ( ii < m_ )
   {
      hh = ( spt->vct[ii+ONE][null] ) - ( spt->vct[ii][null] );
      ( spt->intgr ) += \
          .5*hh*(( spt->vct[ii][ONE] ) + ( spt->vct[ii+ONE][ONE] ));
      ( spt->intgr ) -= \
          (( 1./24. ) * ( hh*hh*hh ) * ( ss[ii] + ss[ii+ONE] ));

      ii++ ;
   }; 

   return spt;
}
/*======================== end of function spline(*) =========================*/
/*************************** end of file spline.c *****************************/
