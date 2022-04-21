/* [ file: readopr.h ] */
/*******************************************************************************
*                                                                              *
*   Operations file [type 'opr.log<N>'] reader function rread_operts(*)        *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/

short rread_operts( char *filename, char mode )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern struct blcstruc blc;
   extern struct transfer trf;
*/
/* declarations: */

   static FILE
      *operats = NULL;

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
      operats = fopen( fleptr, "r+" );

      while ( operats == NULL )
      {

         printf( "\n Operations file %s not found "
            "in present directory\n ", fleptr );
         strcpy(( csp->rqfrm ), "points" );
         strcpy(( csp->rqstr ), "Re-enter filename [ Escape: enter null ]" );
      /* strcpy(( csp->dfstr ), "0" ); */
/*............................................................................*/
         csp = txcnsl( csp );       /*                                        */
/*................................*/
         strcpy( fleptr, ( csp->instr ));

         if ( *fleptr == '0' )
            return ONE;

         operats = fopen( fleptr, "r+" );
      };
/*............................................................................*/
/* set stream to unbuffered: */

   /* setbuf( operats, buf ); */
/*............................................................................*/

      ii = null;
      while ( ii < MAX_LABEL ) 
      {
         fscanf( operats, scformat, ptr );

         if ( null == strncmp( ptr, "OPERATIONS", FIVE ))
         { 
            fscanf( operats, scformat, ptr );
            strncpy(( spt->tpt->name ), DSC_MODEL, SHS_SIZE );
            fscanf( operats, scformat, ptr );
            strncpy(( spt->tpt->text ), ptr, STS_SIZE );

            fscanf( operats, scformat, ptr );
            trf.c[null] = strtol( ptr, endp, DEC );
            fscanf( operats, scformat, ptr );
            fscanf( operats, scformat, ptr ); /* underscore line */

            for ( jj=ONE; jj<=trf.c[null]; jj++)
            {
               fscanf( operats, scformat, ptr );

               if ( mode == 't' )
               { 
                  trf.c[jj] = strtol( ptr, endp, DEC );
                  fscanf( operats, scformat, ptr ); /* string "<---" */
               };

               fscanf( operats, scformat, ptr );

               if ( mode != 't' )
                  trf.c[jj] = strtol( ptr, endp, DEC );
            };

            fclose( operats );
            return null; 
            break;
         }
         else
         {
            ii++;

            if( ii == MAX_LABEL )
            {
               printf( "\n Operations not found "
                  "in file `%s':\n ", fleptr );
               printf( "\n Please re-enter filename [ Escape: "
                  "enter null ] >----> " );
               scanf( "%s", fleptr );
            };
         };/* end if( ptr != "OPERATIONS" ) */
      }; /* end while( ii < MAX_LABEL ) */
   }; /* end while( *fleptr != '0' */

   fclose( operats );
   return null;
}
/*============================================================================*/
/********************* end of function rread_operts(*) ************************/
