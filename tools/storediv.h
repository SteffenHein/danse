/* [ file: storediv.h ] */
/*******************************************************************************
*                                                                              *
*   DSC model topology storage function store_divsns(*)                        *
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
# define STORE_DIV( FN, MD ) \
{ \
   fprintf((FN), spformat, ( spt->tpt->name )); \
   fprintf((FN), spformat, ( spt->tpt->text )); \
 \
  /* if ( null == strcmp( fleptr, tmpfle )) */ \
   if (( MD ) == 't' ) \
   { \
      fprintf((FN), "\n%-3d          %s", trf.n[null], trf.ntx[null] ); \
 \
      fprintf((FN), "\n" ); \
      ii = null; do \
      { \
         fprintf((FN), "%c", '-' ); \
      } while(( ii++ ) < ( LINLEN+17 )); \
 \
      for ( ii=ONE; ii<=trf.n[null]; ii++ ) \
         fprintf((FN), "\n%-4d <---    n%03d=%s", \
         trf.n[ii], ii, trf.ntx[ii] ); \
 \
      fprintf((FN), "\n\n%-3d          %s", blc.m[null], trf.mtx[null] ); \
 \
      fprintf((FN), "\n" ); \
      ii = null; do \
      { \
         fprintf((FN), "%c", '-' ); \
      } while(( ii++ ) < ( LINLEN+17 )); \
 \
      for ( ii=ONE; ii<=blc.m[null]; ii++ ) \
         fprintf((FN), "\n%-4d <---    m%03d=%s", \
         blc.m[ii], ii, trf.mtx[ii] ); \
   } \
   else \
   { \
      fprintf((FN), "\n%-3d |%s", trf.n[null], trf.ntx[null] ); \
 \
      fprintf((FN), "\n" ); \
      ii = null; do \
      { \
         fprintf((FN), "%c", '-' ); \
      } while(( ii++ ) < ( LINLEN+9 )); \
 \
      for ( ii=ONE; ii<=trf.n[null]; ii++ ) \
         fprintf((FN), "\nn%03d=%s: %3d", ii, trf.ntx[ii], trf.n[ii] ); \
 \
      fprintf((FN), "\n\n%-3d |%s", blc.m[null], trf.mtx[null] ); \
 \
      fprintf((FN), "\n" ); \
      ii = null; do \
      { \
         fprintf((FN), "%c", '-' ); \
      } while(( ii++ ) < ( LINLEN+9 )); \
 \
      for ( ii=ONE; ii<=blc.m[null]; ii++ ) \
         fprintf((FN), "\nm%03d=%s: %3d", ii, trf.mtx[ii], blc.m[ii] ); \
   }; \
   fprintf((FN), "\n" ); \
}
/*============================================================================*/

short store_divsns( char *filename, char mode )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
*/
/* declarations: */

   static FILE 
     *divisns = NULL;

   static long 
      ll = null;

   static short
      ii = null;

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
/*----------------------------------------------------------------------------*/
   strcpy( fleptr, filename );
   divisns = fopen( fleptr, "w+" );
/*............................................................................*/
/* set stream unbuffered: */

   setbuf( divisns, NULL );

   ll = null;
   ii = null;
   while ( ii < IPT_MAXLBL )
   {
      fscanf( divisns, scformat, ptr );

      if ( null == strncmp( ptr, "DIVISIONS", ( size_t) SEVEN ))
      {
         ll = ftell( divisns ) + ONE;
            break;
      }
      else
      {
         ii++;
         if( ii == IPT_MAXLBL )
         {
            if ( mode == 't' ) /* temporary file */
            {
               fprintf( divisns, spformat, "The actually charged" );
               fprintf( divisns, spformat, "---> DIVISIONS" );
            }
            else
               fprintf( divisns, "%s%s\n", "DIVISIONS-", ( spt->flbl ));

            ll = ftell( divisns );
         };
      };
   };

   fseek( divisns, ll, SEEK_SET );

   STORE_DIV( divisns, mode );

   nseconds = time( timer );
   timeptr = ctime( &nseconds );

   fprintf( divisns, "\nDSC model divisions logfile %s ", fleptr );
   fprintf( divisns, "created:\n%24s", timeptr );
   fprintf( divisns, "%c", EOF );

   fclose( divisns );

   return null;
}
/*============================================================================*/
/******************** end of function store_divsns(*) *************************/
