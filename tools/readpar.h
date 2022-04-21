/* [ file: readpar.h ] */
/*******************************************************************************
*                                                                              *
*   Parameter file [type 'par.log<N>'] reader function rread_params(*)         *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

short rread_params( char *filename, char mode )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern struct transfer trf;
*/
/* declarations: */

   static FILE
      *paramtrs = NULL;

   static short 
      ii = null,
      jj = null,
      parameters;

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
      paramtrs = fopen( fleptr, "r" );

      while ( paramtrs == NULL )
      {

         printf( "\n Parameter file %s not found "
            "in present directory\n ", fleptr );
         strcpy(( csp->rqfrm ), "points" );
         strcpy(( csp->rqstr ), "Re-enter filename [Escape: enter null]" );
      /* strcpy(( csp->dfstr ), "0" ); */
/*............................................................................*/
         csp = txcnsl( csp );       /*                                        */
/*................................*/
         strcpy( fleptr, ( csp->instr ));

         if ( *fleptr == '0' ) 
         {
            fclose( paramtrs );
            return null;
         };
         paramtrs = fopen( fleptr, "r" );
      };
/*............................................................................*/
/* set stream to unbuffered: */

   /* setbuf( paramtrs, buf ); */
/*............................................................................*/

      ii = null;
      while ( ii < MAX_LABEL )
      {
         fscanf( paramtrs, scformat, ptr );

         if ( null == strncmp( ptr, "PARAMETERS", ( size_t ) FIVE ))
         {
            fscanf( paramtrs, scformat, ptr );
            strncpy(( spt->ppt->name ), DSC_MODEL, SHS_SIZE );
            fscanf( paramtrs, scformat, ptr );
            strncpy(( spt->ppt->text ), ptr, STS_SIZE );

            fscanf( paramtrs, scformat, ptr ); /* time/frequency domain idntf. */

            if (( null == strncmp( ptr, "FREQUENCY_DOMAIN", ( size_t ) ONE ))||
               ( null == strncmp( ptr, "frequency_domain", ( size_t ) ONE )))
               strcpy(( spt->ppt->domain ), "frequency_domain" );
            else
               strcpy(( spt->ppt->domain ), "time_domain" );

            fscanf( paramtrs, scformat, ptr );
            parameters = strtol( ptr , endp, DEC );
            trf.s[null] = parameters; 
            fscanf( paramtrs, scformat, ptr ); /* string "parameters", e.g. */
            fscanf( paramtrs, scformat, ptr ); /* string "INTERNATIONAL ..." */
            fscanf( paramtrs, scformat, ptr ); /* underscore line */

            for ( jj=ONE; jj<=parameters; jj++ )
            {
               fscanf( paramtrs, scformat, ptr );

               if ( mode == 't' )
               {
                  trf.s[jj] = strtod( ptr , endp );
                  fscanf( paramtrs, scformat, ptr ); /* string "<---" */
               };

               fscanf( paramtrs, scformat, ptr );

               if ( mode != 't' )
                  trf.s[jj] = strtod( ptr , endp );
            };

            fclose( paramtrs );
            return null;
            break;
         }
         else
         {
            ii++;
            if ( ii == MAX_LABEL )
            {
               printf( "\n parameters not found "
                  "in file `%s':\n ", fleptr );
               printf( "\n Please re-enter filename [ Escape: "
                  "enter null ] >----> " );
               scanf( "%s", fleptr );
            }
         };/* end if( ptr != "PARAMETERS" ) */
      }; /* end while( ii < MAX_LABEL ) */
   }; /* end while( *flelbl == '0' ) */

   fclose( paramtrs );
   return null;
}
/*============================================================================*/
/********************* end of function rread_params(*) ************************/
