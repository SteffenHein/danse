/* [ file: setmed.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function setmed(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Switches media index of mesh cells ( mdp->ci ),...,( mdp->cf ) to          *
*   medium indexed ( mdp->idx ), if the latter is yet defined, or defines      *
*   a new medium, indexed ( mpt->idx ), if that index is not yet assigned      *
*   to a formerly defined medium.                                              *
*   In the second case, the function reads the material parameters stored      *
*   in struct med into the media parameter structure &( stp->ppr->mpt )        *
*   and assigns index ( mdp->idx ) to the thereby defined new medium.          *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: April 14, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
# include "formtp.h"
*/
/*----------------------------------------------------------------------------*/
# define INCLDE_MED 0
/*............................................................................*/
# if INCLDE_MED == 1
# include "./setmpt.h"
# elif INCLDE_MED == 2
typedef struct
{            
   signed char
      rtn, /* any returned character */
      opt; /* any returned option */

   short
      idx; /* media index */

   long
      ci, /* initial cell, to that media label <idx> shall be assigned */
      cf; /* final cell, to that ... */

/* E/H field relevant media parameters: */
   double
      eps, /* rel. permittivity tensor */
      myr, /* rel. permeability tensor */
      ke,  /* electric conductivity tns. [A/(V*m)] */
      km;  /* magnetic conductivity tns. [V/(A*m)] */
/*............................................................................*/
# if DSC_HCRMDE != 0
/* thermal [diffusion] media parameters: */

   double
      cv, /* specific heat [J/(K*m^3)] */
      kh; /* heat current conductivity [W/(K*m)] */
/*............................................................................*/
# if DSC_FLDMDE != 0
/* fluid media parameters: */

   short
      cnn; /* fluid connected component index; POSITIVE integer */

   double
      rm, /* mean mass density [Kg/m^3] */
      tm, /* mean temperature [C] */
      bm, /* mean thermal expansion coefficient [1/K] */
      cm, /* mean [ adiabatic ] compression coefficient [1/Pa] */
      ny, /* dynamic viscosity [Kg/(sec*m)] */
      q1, /* Cp/Cv - 1 [dimensionless] */
      td, /* dissipation time constant [sec] */
      LL; /* characteristic length [in Prandtl turbulence model, e.g.] */

   double
      gp[THREE], /* pressure gradient [Pa/m]*/
      gr[THREE]; /* gravitation acceleration [m/sec^2]*/

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* Gyrotropic media parameters: */

   double
      no,        /* plasma electron density  [1/m^3] */
      tr,        /* plasma curr.relax.time [seconds] */
      tg,        /* spin-spin  relax.time  [seconds] */
      ld,        /* LANDE factor */
      mi[DIMNS], /* int. magn. flux, plasma  [Tesla] */
      ms[DIMNS], /* saturation magnetization [Tesla] */
      hg[DIMNS]; /* internal magnetic field    [A/m] */

   char
      type[SHS_SIZE];
} MEDIUM;
# endif /* INCLDE_MED == 2 */
/*----------------------------------------------------------------------------*/
static MEDIUM medium = {null};
/*----------------------------------------------------------------------------*/

/*============================================================================*/

MEDIUM *setmed( MEDIUM *mdp )
{ 
/* allusions: */
/*
   extern FORMSTATE *spt;
*/
/* declarations: */

   static long
      ll = null;

   static short
      ii = null;

   static signed char
      pp = null;

   static struct media
     *mpt = NULL;
   
   static MEDIUM
     *rtp = &medium;

# if DSC_HCRMDE != 0
# if DSC_FLDMDE != 0
   static struct fconnect
     *fcp = NULL;
# endif
# endif /* DSC_HCRMDE != 0 */
/*----------------------------------------------------------------------------*/
   mpt = ( spt->ppt->mpt );
   rtp = &medium;
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE != 0
   fcp = ( spt->ppt->fcp );
# endif
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   if ( mdp == NULL ) /* initialize structure med [ type MEDIUM ] */
   {
      ( rtp->opt ) = null;

      ( rtp->rtn ) = null;
      ( rtp->idx ) = null;

      ( rtp->ci ) = null;
      ( rtp->cf ) = null;

      ( rtp->eps ) = 1.;
      ( rtp->myr ) = 1.;
      ( rtp->ke ) = ZERO;
      ( rtp->km ) = ZERO;
/*............................................................................*/
# if DSC_HCRMDE != 0
      
      ( rtp->kh ) = ZERO;
      ( rtp->cv ) = ZERO;
/*............................................................................*/
# if DSC_FLDMDE != 0
      ( rtp->cnn ) = null; 

      ( rtp->rm ) = ZERO;
      ( rtp->tm ) = ZERO;
      ( rtp->bm ) = ZERO;
      ( rtp->cm ) = ZERO;
      ( rtp->ny ) = ZERO;
      ( rtp->q1 ) = ZERO;
      ( rtp->td ) = ZERO;

      ii = null; do           
      {
	 ( rtp->gr[ii] ) = ZERO;
         ( rtp->gp[ii] ) = ZERO;
      } while(( ++ii ) < THREE );
/*............................................................................*/
# if ( TURBMOD == 1 )||( TURBMOD == 2 ) /* Prandtl turbulence models */
     ( rtp->LL ) = ZERO;
# endif /* TURBMOD ! null */
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      ii = null; do           
      {
         ( rtp->type[ii] ) = null;
      } while(( ++ii ) < SHS_SIZE );

      return rtp;
   }
   else if ( blc.cov != null )
      return rtp;
   else if ((( mdp->opt ) == 'c' )
          &&(( mpt->idx[null] ) < (( mdp->idx ) - ONE )))
   {
      fprintf( stderr, "\n\n ERROR message from media switching function %s(*)",
         __func__ );
      fprintf( stderr, "\n Illegal medium index %d !!!", ( mdp->idx ));
      fprintf( stderr, "\n Medium (indexed) %d has not yet been defined.",
         (( mdp->idx ) - ONE ));
      fprintf( stderr, "\n [ Media must be labelled in natural order" );
      fprintf( stderr, "\n   without gaps - i.e. no numbers omitted.]\n\n" );

      exit( EXIT_FAILURE );
   };
/*............................................................................*/
   if ( null < ( mdp->ci ))
   {
      for( ll=( mdp->ci ); ll<=( mdp->cf ); ll++ )
      {
         if(( mdp->opt ) == 'c' ) /* option: 'c'oordinates */
            ( mpt->idx[ll] ) = ( mdp->idx );

         if(( mdp->idx ) == null ) /* trivial cell */
         { /* encapsulate cell with electric walls */
/*............................................................................*/
# if DSC_ADJCLL == 0
            for ( pp=null; pp<FACES; pp++ )
               ( spt->tpt->mn[ll][(int)pp] ) = -ONE;
# elif DSC_ADJCLL == 1
            for ( pp=null; pp<PORTS; pp++ )
               ( spt->tpt->mn[ll][(int)pp] ) = -ONE;
# endif /* DSC_ADJCLL == 1 */
/*............................................................................*/
         };
      };
   };
/*............................................................................*/
/* [ mpt->idx[null]: the number ( = maximum index ) of yet defined media ] */

   if (( mdp->opt ) == 'c' )
   {
      if (( mpt->idx[null] ) < ( mdp->idx ))
      {
         ( mpt->idx[null] ) = ( mdp->idx );
         ( mpt->tg[( mdp->idx )] ) = ZERO;

         if ( fabs( mdp->eps ) < 1.e-277 )
            strcpy(( mpt->type[( mdp->idx )] ), "trivial_E" );
         else if ( fabs( mdp->myr ) < 1.e-277 )
            strcpy(( mpt->type[( mdp->idx )] ), "trivial_M" );
         else
            strcpy(( mpt->type[( mdp->idx )]), "non-trv_EM" );
/*............................................................................*/
# if DSC_HCRMDE != 0
         if ( fabs( mdp->kh ) < 1.e-277 )
            strcat(( mpt->type[( mdp->idx )] ), "_trv_hcr" );
         else if ( fabs( mdp->cv ) < 1.e-277 )
            strcat(( mpt->type[( mdp->idx )] ), "_trv_hcr" );
         else
         {
            ( mpt->kh[( mdp->idx )] ) = ( mdp->kh );
            ( mpt->cv[( mdp->idx )] ) = ( mdp->cv );
/*............................................................................*/
# if DSC_FLDMDE != 0
            if ( fabs( mdp->rm ) < 1.e-277 )
            {
               strcat(( mpt->type[( mdp->idx )] ), "_hcurr" );
	       ( mpt->rm[( mdp->idx )] ) = ZERO;
            }
            else
            {
               strcat(( mpt->type[( mdp->idx )] ), "_fluid" );
	       ( mpt->rm[( mdp->idx )] ) = ( mdp->rm );
	       ( mpt->tm[( mdp->idx )] ) = ( mdp->tm );
	       ( mpt->bm[( mdp->idx )] ) = ( mdp->bm );
	       ( mpt->cm[( mdp->idx )] ) = ( mdp->cm );
	       ( mpt->ny[( mdp->idx )] ) = ( mdp->ny );
	       ( mpt->q1[( mdp->idx )] ) = ( mdp->q1 );
	       ( mpt->td[( mdp->idx )] ) = ( mdp->td );

               ii = null; do
               {
	          ( mpt->gr[( mdp->idx )][ii] ) = ( mdp->gr[ii] );
	          ( mpt->gp[( mdp->idx )][ii] ) = ( mdp->gp[ii] );
               } while(( ++ii ) < THREE );
/*............................................................................*/
# if ( TURBMOD == 1 )||( TURBMOD == 2 ) /* Prandtl turbulence models */
	       ( mpt->LL[( mdp->idx )] ) = ( mdp->LL );
# endif /* TURBMOD ! null */
/*............................................................................*/
            };
# else /* DSC_FLDMDE == 0 */
/*............................................................................*/
            strcat(( mpt->type[( mdp->idx )] ), "_hcurr" );
# endif /* DSC_FLDMDE == 0 */
/*............................................................................*/
         };
# endif /* if DSC_HCRMDE != 0 */
/*............................................................................*/
         ii = null;
         do
         {
            ( mpt->ep[( mdp->idx )][ii] ) = ( mdp->eps );
            ( mpt->my[( mdp->idx )][ii] ) = ( mdp->myr );
            ( mpt->ke[( mdp->idx )][ii] ) = ( mdp->ke );
            ( mpt->km[( mdp->idx )][ii] ) = ( mdp->km );
            ( mpt->ms[( mdp->idx )][ii] ) = ZERO;
            ( mpt->hg[( mdp->idx )][ii] ) = ZERO;
         } while(( ++ii ) < THREE );
         do
         {
            ( mpt->ep[( mdp->idx )][ii] ) = ZERO;
            ( mpt->my[( mdp->idx )][ii] ) = ZERO;
            ( mpt->ke[( mdp->idx )][ii] ) = ZERO;
            ( mpt->km[( mdp->idx )][ii] ) = ZERO;
         } while(( ++ii ) < SIX );
      }; /* end if (( mpt->idx[null] ) < ( mdp->idx )) */
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE != 0

      if ( null != strstr(( mpt->type[( mdp->idx )]), "fluid" ))
      {
         if (( mdp->cnn ) < ONE )
         {
            fprintf( stderr, "\n\n ERROR message from media switching "
               "function %s(*):", __func__ );
            fprintf( stderr, "\n Illegal connected component label cnn =" );
            fprintf( stderr, "\n = %d assigned to fluid medium no %d !!!",
               ( mdp->cnn ), ( mdp->idx ));
            fprintf( stderr, "\n - Assign a POSITIVE index cnn <= %d", NFCNN );
            fprintf( stderr, "\n [ = macro NFCNN ] to fluid medium no %d.\n\n",
               ( mdp->idx ));

            exit( EXIT_FAILURE );
         }
         else if ( NFCNN < ( mdp->cnn ))
         {
            fprintf( stderr, "\n\n ERROR message from media switching "
               "function %s(*):", __func__ );
            fprintf( stderr, "\n Too many fluid components defined in "
               "DSC model !!!" );
            fprintf( stderr, "\n [ Number %d exceeds maximum number %d "
               "= macro NFCNN, ", ( mdp->cnn ), NFCNN );
            fprintf( stderr, "fixed in file FORMER.CONF." );
            fprintf( stderr, "\n   - Change macro only in compliance with" );
            fprintf( stderr, "\n   memory resources.]\n\n" );

            exit( EXIT_FAILURE );
         };

         for( ll=( mdp->ci ); ll<=( mdp->cf ); ll++ )
            ( fcp->cnn[ll] ) = ( mdp->cnn );

         if (( fcp->cnn[null] ) < ( mdp->cnn ))
            ( fcp->cnn[null] ) = ( mdp->cnn );
      }
      else /* non fluid medium */
      {
         for( ll=( mdp->ci ); ll<=( mdp->cf ); ll++ )
            ( fcp->cnn[ll] ) = null;
      };
# endif
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
   }; /* end if (( mpt->opt ) == 'c' ) */

   return rtp;
}
/*============================================================================*/
/*********************** end of function setmed(*) ****************************/
