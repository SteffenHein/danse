/* [ file: storepar.h ] */
/*******************************************************************************
*                                                                              *
*   DSC model parameter storage function store_params(*)                       *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: April 14, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# ifndef LINLEN
   # define LINLEN 53
# endif
/*----------------------------------------------------------------------------*/
# define STORE_PAR( FN, MD ) \
{ \
   fprintf((FN), spformat, ( spt->ppt->name )); \
   fprintf((FN), spformat, ( spt->ppt->text )); \
 \
   if( *( spt->ppt->domain ) == 'f' ) \
      fprintf((FN), "%s\n\n", "FREQUENCY_DOMAIN" ); \
   else \
      fprintf((FN), "%s\n\n", "TIME_DOMAIN" ); \
 \
   parameters = ( short ) trf.s[null]; \
 \
  /* if ( null == strcmp( fleptr, tmpfle )) */ \
   if (( MD ) == 't' ) \
   { \
      fprintf((FN), "%-3d %s", \
         parameters, trf.stx[null] ); \
 \
      fprintf((FN), "\nINTERNATIONAL_UNITS_[if_not_otherwise_specified]\n" ); \
 \
      ii = null; do \
      { \
         fprintf((FN), "%c", '-' ); \
      } while(( ii++ ) < ( LINLEN+25 )); \
 \
      for ( ii=ONE; ii<=parameters; ii++ ) \
         fprintf((FN), "\n%+.8e <--- s%03d=%s", \
         trf.s[ii], ii, trf.stx[ii] ); \
 \
      fprintf((FN), "\n" );\
   } \
   else \
   { \
      fprintf((FN), "%-3d %s", parameters, trf.stx[null] ); \
 \
      fprintf((FN), "\nINTERNATIONAL_UNITS_[if_not_otherwise_specified]\n" ); \
 \
      ii = null; do \
      { \
         fprintf((FN), "%c", '-' ); \
      } while(( ii++ ) < ( LINLEN+25 )); \
 \
      for ( ii=ONE; ii<=parameters; ii++ ) \
         fprintf((FN), "\ns%03d=%s: % .12e", \
         ii, trf.stx[ii], trf.s[ii] ); \
 \
      fprintf((FN), "\n\nwaveguide_type" ); \
      ii = null; do \
      { \
         fprintf((FN), "%c", '_' ); \
      } while(( ii++ ) < ( LINLEN-10 )); \
      fprintf((FN), ":  %s\n", trf.wgtype ); \
      fprintf((FN), "\n[ operation parameters at f = " \
         "%.16e GHz\n", ( trf.fr*1.e-9 )); \
 \
      if ( null == strncmp( trf.wgtype, "coax", ( size_t ) THREE )) \
         wve = wvepar( trf.wgtype, 2.*trf.ra, 2.*trf.ri, 1., 1., trf.fr ); \
      else if ( null == strncmp( trf.wgtype, "elliptic", ( size_t ) THREE )) \
         wve = wvepar( trf.wgtype, 2.*trf.ra, 2.*trf.ri, 1., 1., trf.fr ); \
      else \
         wve = wvepar( trf.wgtype, trf.a, trf.b, 1., 1., trf.fr ); \
 \
      if (( wve->q ) < 1. ) \
      { \
         fprintf((FN), "  free_wavelength.............[m]: " \
            "%.15e\n", ( wve->w0 )); \
         fprintf((FN), "  waveguide_wavelength........[m]: " \
            "%.15e\n", ( wve->wg )); \
         fprintf((FN), "  TEmn_mode_line_impedance__[Ohm]: " \
            "%.15e\n", ( wve->zh )); \
         fprintf((FN), "  phase_velocity............[m/s]: " \
            "%.15e\n", ( wve->vp )); \
         fprintf((FN), "  group_velocity............[m/s]: " \
            "%.15e ]\n", ( wve->vg )); \
         fprintf((FN), "\ncritical_frequency" ); \
         ii = null; do \
         { \
            fprintf((FN), "%c", '_' ); \
         } while(( ii++ ) < ( LINLEN-18 )); \
         fprintf((FN), "[Hz]:  %.12e\n", ( wve->fc )); \
      } \
      else \
      { \
         fprintf((FN), " ------------ frequency below cutoff !!! " \
            "------------ ]\n" ); \
      }; \
   }; \
}
/*============================================================================*/

short store_params( char *filename, char mode )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
*/
/* declarations: */

   static FILE 
     *paramtrs = NULL;

   static long 
      ll = null;

   static short
      ii = null,
      parameters = null;

   static char
     *timeptr = null,
      ptr[STS_SIZE] = {null},
      fleptr[STS_SIZE] = {null};

   static const char
     *scformat = "%80s",
     *spformat = "%s\n";

   time_t 
      nseconds = null,
     *timer = null;

   static WVGDPAR
     *wve;

   WVGDPAR
     *wvepar( char *type, double aa, double bb,
              double eps, double my, double ff );
/*----------------------------------------------------------------------------*/

   strcpy( fleptr, filename );
   paramtrs = fopen( fleptr, "w+" );  /* save coordinates on log file */
/*............................................................................*/
/* set stream unbuffered: */

/* setbuf( paramtrs, NULL ); */

   ii = null;
   while ( ii < IPT_MAXLBL )
   {
      fscanf( paramtrs, scformat, ptr );

      if ( null == strncmp( ptr, "PARAMETERS", ( size_t ) SEVEN ))
      {
         ll = ftell( paramtrs ) + ONE;
         break;
      }
      else
      {
         ii++;
         if ( ii == IPT_MAXLBL )
         {
            if ( mode == 't' ) /* temporary file */
            { 
               fprintf( paramtrs, spformat, "The actually charged" );
               fprintf( paramtrs, spformat, "---> PARAMETERS" );
            }
            else
               fprintf( paramtrs, "%s%s\n", "PARAMETERS-", ( spt->flbl ));

            ll = ftell( paramtrs );
         };
      };
   };

   fseek( paramtrs, ll, SEEK_SET );

   STORE_PAR( paramtrs, mode );

/* append material parameters: */

   if ( mode != 't' ) /* eventually append complete configuration files */
   {
/* append material parameters file: */

      fprintf( paramtrs, "\n==============================="
         "=================================================\n" );
      fprintf( paramtrs, "\n%s%s\n", "DIVISIONS-", ( spt->flbl ));

      STORE_DIV( paramtrs, 'p' );

/* append operation parameters file: */

      fprintf( paramtrs, "\n==============================="
         "=================================================\n" );
      fprintf( paramtrs, "\n%s%s\n", "OPERATIONS-", ( spt->flbl ));

      STORE_OPR( paramtrs, 'p' );
   };

   nseconds = time( timer );
   timeptr = ctime( &nseconds );

   fprintf( paramtrs, "\nDSC model parameter input file %s", fleptr );
   fprintf( paramtrs, " created:\n%24s", timeptr );
   fprintf( paramtrs, "%c", EOF );

   fclose( paramtrs );

   return null;
}
/*============================================================================*/
/********************** end of function store_params(*) ***********************/
