/* [ file: readdiv.h ] */
/*******************************************************************************
*                                                                              *
*   Divisions file [type 'div.log<N>'] reader function rread_divsns(*)         *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

short rread_divsns( char *filename, char mode )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern struct blcstruc blc;
   extern struct transfer trf;
*/
/* declarations: */

   static FILE
      *divisns = NULL;

   static short 
      ii = null,
      jj = null;

   static char
      ptr[STS_SIZE] = {null},
      fleptr[STS_SIZE] = {null};

   static const char
     *scformat = "%80s";
/*
   static char 
     *buf = NULL;
*/
   static char 
    **endp = NULL;

   static TXCNSL 
     *csp; /* cf. txcnsl(*) */

   TXCNSL 
     *txcnsl( TXCNSL *csp );
/*----------------------------------------------------------------------------*/
/* initialize text console */

   csp = txcnsl( null ); /* initialize the text console */
/*............................................................................*/

   strcpy( fleptr, filename );

   while( *fleptr != '0' )
   {
      divisns = fopen( fleptr, "r" );

      while ( divisns == NULL )
      {

         printf( "\n Divisions file %s not found "
            "in present directory\n ", fleptr );
         strcpy(( csp->rqfrm ), "points" );
         strcpy(( csp->rqstr ), "Re-enter filename [Escape: enter null]" );
      /* strcpy(( csp->dfstr ), "0" ); */
/*............................................................................*/
         csp = txcnsl( csp );       /*                                        */
/*................................*/
         strcpy( fleptr, ( csp->instr ));

         if ( *fleptr == '0' )
            return ONE;

         divisns = fopen( fleptr, "r" );
      };
/*............................................................................*/
/* set stream to unbuffered: */

    /* setbuf( divisns, buf ); */
/*............................................................................*/

      ii = null;
      while ( ii < MAX_LABEL ) 
      {
         fscanf( divisns, scformat, ptr );

         if ( null == strncmp( ptr, "DIVISIONS", FIVE ))
         { 
            fscanf( divisns, scformat, ptr );
            strncpy(( spt->tpt->name ), DSC_MODEL, SHS_SIZE );
            fscanf( divisns, scformat, ptr );
            strncpy(( spt->tpt->text ), ptr, STS_SIZE );

            fscanf( divisns, scformat, ptr );
            trf.n[null] = strtol( ptr, endp, DEC );
            fscanf( divisns, scformat, ptr );
            fscanf( divisns, scformat, ptr ); /* underscore line */

            for ( jj=ONE; jj<=trf.n[null]; jj++)
            {
               fscanf( divisns, scformat, ptr );

               if ( mode == 't' )
               {
                  trf.n[jj] = strtol( ptr, endp, DEC );
                  fscanf( divisns, scformat, ptr ); /* string "<---" */
               };

               fscanf( divisns, scformat, ptr );

               if ( mode != 't' )
                  trf.n[jj] = strtol( ptr, endp, DEC );
            };

            fscanf( divisns, scformat, ptr );
            blc.m[null] = strtol( ptr, endp, DEC );
            fscanf( divisns, scformat, ptr );
            fscanf( divisns, scformat, ptr ); /* underscore line */

            for ( jj=ONE ; jj<=blc.m[null] ; jj++)
            {
               fscanf( divisns, scformat, ptr );

               if ( mode == 't' )
               {
                  blc.m[jj] = strtol( ptr, endp, DEC );
                  fscanf( divisns, scformat, ptr ); /* string "<---" */
               };

               fscanf( divisns, scformat, ptr );

               if ( mode != 't' )
                  blc.m[jj] = strtol( ptr, endp, DEC );
            };

            fclose( divisns );
            return null; 
            break;
         }
         else
         {
            ii++;

            if( ii == MAX_LABEL )
            {
               printf( "\n Divisions not found "
                  "in file `%s':\n ", fleptr );
               printf( "\n Please re-enter filename [ Escape: "
                  "enter null ] >----> " );
               scanf( "%s", fleptr );
            };
         };/* end if( ptr != "DIVISIONS" ) */
      }; /* end while( ii < MAX_LABEL ) */
   }; /* end while( *fleptr != '0' */

   fclose( divisns );
   return null;
}
/*============================================================================*/
/********************* end of function rread_divsns(*) ************************/
