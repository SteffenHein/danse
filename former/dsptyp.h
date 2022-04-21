/* [ file: dsptyp.h ] */
/* Update: 13 April 2007 */
/*----------------------------------------------------------------------------*/
typedef struct
{
   char 
      rtn;

   char 
      option,
      messge[LGS_SIZE];

   short
      rcsps; /* relative curser position */
   
   long
      fleps, /* file position pointer [ fleps = ftell(*) etc. ] */
      state,
      range;

   FILE 
     *display;

} DSPLAY;
/*************************** end of file dsptyp.h *****************************/
