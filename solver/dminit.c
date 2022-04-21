/* [ file: dminit.c ] */
# define DO_DMINIT "dminit(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function dminit(*)                                                  *
*                                                                              *
*   This function permits terminal free ['daemon'] operation of a program      *
*   in a Unix related system environment [ BSD, Linux etc.]. The program       *
*   can thus be automatically re-started at boot time, after power chute       *
*   and subsequent reboot, for instance. Obviously, this option is particu-    *
*   larily well suited for extensive batch file operation, in the absence      *
*   of a system operator [ during weekend, e.g.].                              *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <fcntl.h>
/*----------------------------------------------------------------------------*/
/* another system header: */
# include "../math/unxsys.h"
/*----------------------------------------------------------------------------*/
static DMINIT dmn = { null };
/*============================================================================*/

DMINIT *\
dminit( DMINIT *dpt )
{
   static FILE 
     *dmstrm;

   static short
      ii = null;

   static DMINIT 
     *rpt = &dmn;
/*----------------------------------------------------------------------------*/
/* initialize structure DMINIT dmn: */

   if( dpt == null )
   {
      ( rpt->uid ) = null;
      ( rpt->pid ) = null;

      ii = null;
      while( ii < STS_SIZE )
      {
         ( rpt->dmndir )[ii] = null;
         ( rpt->logfle )[ii] = null;
         ( rpt->pidfle )[ii] = null;
         ii++ ;
      };

      ( rpt->rtn ) = null;
      return rpt;
   };
/*............................................................................*/
   rpt = dpt;

/*............................................................................*/
/* check permission: */

   if (( rpt->uid ) == null )
   {
      ( rpt->uid ) = ( long ) getuid( );

      if ( 3 < ( rpt->uid ))
      {
         fprintf( stderr, "only superuser can start solver daemon\n" );
         exit( EXIT_FAILURE );
      };
   };
/*............................................................................*/
/* fork to child process: */

   if ((( rpt->pid ) = ( long ) fork( )) < null )
   {
      fprintf( stderr, "can't fork\n" );
      ( rpt->rtn ) = -1;
      return( rpt );
   }
   else if (( rpt->pid ) != null )
      exit( EXIT_SUCCESS ); /* parent process says good-bye */
/*...........................................................................*/
   setsid( ); /* becomes controling process */

   chdir( dpt->dmndir ); /* change into directory where daemon resides */

   umask( null ); /* clear file protection bit mask */
/*...........................................................................*/
/* store pid: */

   ( rpt->pid ) = ( long ) getpid( ); /* new process id */

   dmstrm = fopen( dpt->pidfle, "w+" );
   if ( null == dmstrm )
   {
      fprintf( stderr, "can't write pid into %s\n", dpt->pidfle );
      exit( EXIT_FAILURE );
   }
   else
   {
      fprintf( dmstrm, "%ld", rpt->pid );
      fclose( dmstrm );
   };
/*...........................................................................*/
/* check log file: */

   dmstrm = fopen( dpt->logfle, "a+" );
   if ( null == dmstrm )
   {
      fprintf( stderr, "can't open logfile %s\n", dpt->logfle );
      exit( EXIT_FAILURE );
   }
   else
      fclose( dmstrm );
/*...........................................................................*/
   ( rpt->rtn ) = null;
   return( rpt );
}
/*============================================================================*/
/************************* end of function dminit(*) **************************/
