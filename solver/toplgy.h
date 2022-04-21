/* [ file: toplgy.h ] */
# define DO_TOPLGY "toplgy(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function toplgy(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   This subroutine  enters the  topological structure of the DSC system       *
*   ( top.name ) from file  DSC_PRFX.TOPOLOGY_FILE<n>.   The S-matrix          *
*   of the referenced system is stored in file DSC_PRFX.S_MATRIX_FILE<n>       *
*   and transferred to the main program by function smtrix(*).                 *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* the following macro should be defined in "../CONFIG.H"              */
# ifndef DSC_ADJCLL      /* assign neighbouring cell index top.mn[][k] so:    */
   # define DSC_ADJCLL 0 /* 0: to faces [ k is a face index; 0 <= k < 6 ]     */
# endif                  /* 1: to ports [ k is a port index; 0 <= k < 12 ]    */
/*----------------------------------------------------------------------------*/
static short toplbl[DSC_JOBS] = {null};
/*============================================================================*/

short toplgy( const short jj )
{
/* allusions: */
/*
   extern struct solverstat solver;
   extern struct topology top;
   extern short toplbl[];
*/
/* declarations: */

   static FILE 
      *tplgfle = NULL;

   static char
      ptr[STS_SIZE] = { null },
     *sysptr = TOPOLOGY_FILE,               /* cf. main program solver.c  */
    **endp = NULL;

   static const char 
     *scformat = "%80s";

   static signed char 
      kk = null;

   static short 
      ind = null;

   static long
      ii = null,
      meshi  = null;
/*............................................................................*/
# if DSC_ADJCLL == 0
   static const char /* 1st and 2nd ports pertinent to faces 0,...,5: */
      prt1[FACES] = { 7, 5, 11, 9, 3, 1 }, /* 1st port */
      prt2[FACES] = { 10, 8, 2, 0, 6, 4 }; /* 2nd port */
# elif DSC_ADJCLL == 1
   static const char /* faces pertinent to ports 0,1,...,11: */
      fce[PORTS] = { 3, 5, 2, 4, 5, 1, 4, 0, 1, 3, 0, 2 };
# endif
/*............................................................................*/
/* prototypes: */

   char
     *lotos ( long lngint, char length );
/*----------------------------------------------------------------------------*/
/* memory allocations: - [void] */
/*............................................................................*/
# if DSC_LNGNAMES == 1
   strcpy( solver.top, solver.prfx );
   strcat( solver.top, sysptr );
# else
   strncpy( solver.top, solver.prfx, VSS_SIZE );
   strncat( solver.top, sysptr, SHS_SIZE - VSS_SIZE - THREE );
# endif
/*............................................................................*/

   strcat( solver.top, lotos( toplbl[jj], null ));

   tplgfle = fopen( solver.top, "r" );
     
   if ( tplgfle == null ) 
   {
      fprintf( display, "\n\n Error on opening topology "
         "file %s\n", solver.top );
      strncpy( solver.fcterr, DO_TOPLGY, SHS_SIZE );
      strncpy( solver.errmsg, solver.top, SHS_SIZE );
      strncat( solver.errmsg, ":unable_to_open_file ", 
         LGS_SIZE-SHS_SIZE );
      return null; /* abnormal return */
   };
/*............................................................................*/
        
   fscanf( tplgfle, scformat, ptr );
   strncpy(( top.name ), ptr, SHS_SIZE - ONE ); 
   fscanf( tplgfle, scformat, ptr );
   strncpy(( top.text ), ptr, STS_SIZE - ONE );

   fscanf( tplgfle, scformat, ptr ); /* string "______________..." */
   fscanf( tplgfle, scformat, ptr ); /* string "first_mesh_cell_index_:" */
   fscanf( tplgfle, scformat, ptr ); /* long integer string */
   meshi = strtol( ptr, endp, DEC ); /* meshi = initial mesh index  */
                                     /* [ should be ONE ]           */
   if ( meshi != ONE )
   {
      fclose( tplgfle );
      fprintf( display, "\n Truncated system in topology file %s ",
                                                                  solver.top );
      fprintf( display, "\n lowest mesh index  %ld != ONE ", meshi ); 
      strncpy( solver.fcterr, DO_TOPLGY, SHS_SIZE ); 
      strncpy( solver.errmsg, solver.top, SHS_SIZE );
      strncat( solver.errmsg, ":truncated_DSC_mesh ", 
         LGS_SIZE-SHS_SIZE );
      return null; /* abnormal return */
   };
/*............................................................................*/
   fscanf( tplgfle, scformat, ptr ); /* string "last_mesh_cell_index__:" */
   fscanf( tplgfle, scformat, ptr ); /* long integer string */
   top.n = strtol( ptr, endp, DEC ); /* ( top.n ) = final cell index    */
                                     /* [ must be equal to number of cells ]  */
   if ( NODES <  top.n )
   {
      fprintf( display, "\n\n Message from function %s:", DO_TOPLGY );
      fprintf( display, "\n\n Too many nodes defined "
         "in DSC mesh %s !!!", top.name );
      fprintf( display, "\n [ Maximum number is %ld = macro NODES "
         "in SOLVER.CONF or %s.", ( long ) NODES, DO_TOPLGY );
      fprintf( display, "\n - Change macro only in compliance "
         "with memory resources !\n " );
      strncpy( solver.fcterr, DO_TOPLGY, SHS_SIZE );
      strncpy( solver.errmsg, solver.top, SHS_SIZE );
      strncat( solver.errmsg, ":too_many_nodes_defined_", 
         ( LGS_SIZE - SHS_SIZE ));
      return null; /* abnormal return */
   };
/*............................................................................*/
   for ( ii=ONE; ii<=top.n; ii++)
   {
      fscanf( tplgfle, scformat, ptr ); /* mesh index [ even if in general not used ] */
/*............................................................................*/
# if DSC_ADJCLL == 0
      kk=null; do
      {
         fscanf( tplgfle, scformat, ptr );
         top.m[ii][kk] = strtol( ptr, endp, DEC );

         if ( null < top.m[ii][kk] )
         {
            fscanf( tplgfle, scformat, ptr );
            top.p[ii][(int)prt1[kk]] = strtol( ptr, endp, DEC );
            fscanf( tplgfle, scformat, ptr );
            top.p[ii][(int)prt2[kk]] = strtol( ptr, endp, DEC );
            fscanf( tplgfle, scformat, ptr );
            top.f[ii][kk] = strtol( ptr, endp, DEC );
         };
      } while (( ++kk ) < FACES );
# elif DSC_ADJCLL == 1
      kk=null; do
      {
         fscanf( tplgfle, scformat, ptr );
         top.m[ii][kk] = strtol( ptr, endp, DEC );

         if ( null < top.m[ii][kk] )
         {
            fscanf( tplgfle, scformat, ptr );
            top.p[ii][kk] = strtol( ptr, endp, DEC );
            fscanf( tplgfle, scformat, ptr );
            top.f[ii][(int)fce[kk]] = strtol( ptr, endp, DEC );
         };
      } while(( ++kk ) < PORTS );
# endif
/*............................................................................*/
   };

   fclose( tplgfle );
   fprintf( display, "\n\n ====================================="
                     "=========================================" );
   if ( solver.rmfle != ONE )
   {
      fprintf( display, "\n entered: topology file %s", solver.top );
   }
   else if ( solver.rmfle == ONE )
   {
      for ( kk=jj+ONE; kk<solver.nbrjbs; kk++ )
      {
         if ( toplbl[kk] == toplbl[jj] )
         {
            fprintf( display, "\n entered: topology file %s", solver.top );
            goto disp;
         };
      };
      ind = remove( solver.top );
      if ( ind == null )
         fprintf( display, "\n entered and removed: "
            "topology file %s ", solver.top );
   };
     
  disp:  
/*............................................................................*/
# if DSC_DISP == 1

   fprintf( display, "\n\n identifier ( system name ).............: " );

   ind = strlen( top.name );
   if ( 35 < ind )
      fprintf( display, "\n %s\n", top.name );
   else
      fprintf( display, "%s", top.name );

   fprintf( display, "\n text...................................: " );

   ind = strlen( top.text );
   if ( 35 < ind )  
      fprintf( display, "\n %s\n", top.text ); 
   else
      fprintf( display, "%s", top.text );

   fprintf( display, "\n number of mesh cells ..................: %ld ",
      top.n );
   fprintf( display, "\n ------------------------------"
      "------------------------------------------------" );
# endif
/*............................................................................*/
   return ONE; /* normal return */
}
/*============================================================================*/
/*********** end of DSC system topology input function toplgy(*) **************/
