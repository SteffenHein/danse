/* [ file: TIMEFORM.M ] */
/*******************************************************************************
*                                                                              *
*  Macro TIMEFORM(*)                                                           *
*                                                                              *
*  Gives the string pointed to by function ctime a conveniently readable form  *
*                                                                              *
*  (C) SHEIN; Bad Aibling, August 2007                         Steffen Hein    *
*  [ Update: April 13, 2022 ]                               <contact@sfenx.de> *
*                                                                              *
*******************************************************************************/
# define TIMEFORM(NSTR,OSTR) \
{ \
/* clear strings: */ \
\
   kk = null; do \
   { \
      (NSTR)[kk] = null; \
   } while(( ++kk ) < STS_SIZE ); \
\
/* day of the week, e.g. "Wed" : */ \
   strncpy( (NSTR), (OSTR), THREE ); \
\
   strncat( (NSTR), " ", ONE ); \
/* day of the month: */ \
   if ( null == strncmp( (OSTR) + 8 , " ", ONE )) \
   { \
/* replace trailing space with '0': */ \
      strncat( (NSTR), "0", 1 ); \
      strncat( (NSTR), (OSTR) + 9, ONE ); \
   } \
   else \
      strncat( (NSTR), (OSTR) + 8, TWO ); \
\
   strncat( (NSTR), " ", ONE ); \
/* month, e.g. "Aug" : */ \
   strncat( (NSTR), (OSTR) + 4, THREE ); \
\
   strncat( (NSTR), " ", ONE ); \
/* year, e.g. "2007" : */ \
   strncat( (NSTR), (OSTR) + 20, FOUR ); \
\
   strncat( (NSTR), "; ", 4 ); \
/* time, e.g. "15:12:54" : */ \
   strncat( (NSTR), (OSTR) + 11, EIGHT ); \
} \
/************************ end of macro TIMEFORM(*) ****************************/
