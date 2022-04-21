/* [ file: unxsys.h ] */
/*----------------------------------------------------------------------------*/
/* This header may be used to implement some functions [e.g. daemons] on      */
/* UNIX related systems [ such as BSD versions, or Linux, for instance ].     */
/* Included along with other headers [ viz. after them ], it is suited to     */
/* harmonize some differences between those unix-type systems.                */
/*                                                                            */
/* (C) SHEIN; Bad Aibling, February 2007                      Steffen Hein    */
/* [ Update: May 01, 2007 ]                                <contact@sfenx.de> */
/*                                                                            */
/*----------------------------------------------------------------------------*/
# ifndef __unxsys_h /* this condition embraces the entire content of this file*/

# define __unxsys_h
/*----------------------------------------------------------------------------*/
# define _POSIX_SOURCE 1 /* some headers of the POSIX.1 standard will be used */
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

/* the following doesn't work on some older HPUX machines !!! */
/*
   pid_t        
      pid;

   uid_t
      uid;
*/
/* thus it has been replaced by: */

   long
      pid;

   long
      uid;

} DMINIT;

/* the daemom initializer function prototype: */

DMINIT *\
   dminit( DMINIT *dpt );

# endif /* end ifndef __unxsys_h */
/*************************** end of file unxsys.h *****************************/
