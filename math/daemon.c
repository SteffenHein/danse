/* [ file: daemon.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function dminit(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   This function permits terminal free ['daemon'] operation of a program      *
*   in a Unix related system environment [ BSD, Linux etc.]. The program       *
*   can thus be automatically re-started at boot time, after power chute       *
*   and subsequent reboot for instance. Obviously, this option is particu-     *
*   larily well suited for extensive batch file operation, in the absence      *
*   a system operator [ during wheekends, e.g.].                               *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <fcntl.h>
/*----------------------------------------------------------------------------*/
/* This header may be used to implement some functions [ daemons, e.g.] on    */
/* UNIX related systems, such as any BSD version or Linux, for instance.      */
/* Included after other headers it is suitable to harmonize some differences  */
/* between those unix like systems.                                           */
/*                                                                            */
/* (C) SHEIN, Bad Aibling                                    Steffen Hein     */
/* [ Update: October 24, 2005 ]                            <contact@sfenx.de> */
/*----------------------------------------------------------------------------*/
# ifndef __unxsys_h /* ... then mark here the presence of header unxsys.h */
# define __unxsys_h
/*----------------------------------------------------------------------------*/
# include <sys/types.h>	/* required to define some prototypes */
# include <stdio.h>     /* the following just to make things easier ... */
# include <stdlib.h>
# include <string.h>
# include <unistd.h>
/*----------------------------------------------------------------------------*/
# ifndef null
   # define null 0
# endif

# ifndef STS_SIZE
   # define STS_SIZE 80
# endif
 
# define MAXLINE 4096	/* maximum line length */

# define FILE_MODE	(S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)
			/* standard permissions for new files */
# define DIR_MODE	(S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)
			/* standard permissions for new directories */

typedef void Sigfunc(int);
			/* for signal processing functions */
			/* SIG_ERR is not defined in 4.3BSD Reno <signal.h> */
# if	defined(SIG_IGN) && !defined(SIG_ERR)
   # define SIG_ERR ((Sigfunc *)-1)
# endif

# define	min(a,b)	((a) < (b) ? (a) : (b))
# define	max(a,b)	((b) < (a) ? (a) : (b))

/* structure type of daemon initializer function dminit(*) */

typedef struct
{
   signed char
      rtn;
   
   char
      dmndir[STS_SIZE], /* daemon directory [ where all input files reside ] */
      logfle[STS_SIZE], /* the log file [ e.g. '/var/log/elf.log'] */
      pidfle[STS_SIZE]; /* the process number file [ e.g. '/var/run/elf.pid'] */

/* pid_t, uid_t are not known by some older HPUX machines !!! */

   pid_t        
      pid;

   uid_t
      uid;

} DMINIT;

/* the daemom initializer function prototype: */
DMINIT *\
   dminit( DMINIT *dpt );

# endif /* end ifndef __unxsys_h */
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

   if(( rpt->uid ) == null )
   {
      ( rpt->uid ) = getuid( );

      if ( 3 < ( rpt->uid ))
      {
         fprintf( stderr, "only superuser can start DANSE daemon\n" );
         exit( EXIT_FAILURE );
      };
   };
/*............................................................................*/
/* fork to child process: */

   if ((( rpt->pid ) = fork( )) < null )
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

   ( rpt->pid ) = getpid( ); /* new process id */

   dmstrm = fopen( dpt->pidfle, "w+" );
   if ( null == dmstrm )
   {
      fprintf( stderr, "can't write pid into %s\n", dpt->pidfle );
      exit( EXIT_FAILURE );
   }
   else
   {
      fprintf( dmstrm, "%d", rpt->pid );
      fclose( dmstrm );
   };
/*............................................................................*/
/* check log file: */

   dmstrm = fopen( dpt->logfle, "a+" );
   if ( null == dmstrm )
   {
      fprintf( stderr, "can't open logfile %s\n", dpt->logfle );
      exit( EXIT_FAILURE );
   }
   else
      fclose( dmstrm );
/*............................................................................*/
   ( rpt->rtn ) = null;
   return( rpt );
}
/************************* end of function dminit(*) **************************/
