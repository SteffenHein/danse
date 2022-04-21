/* [ file: grphtp.h ] */
/* Update: April 14, 2007 */
/*--------------------------------------------------------------------------*/
/* graphics data transfer structure [ used in function graphp(*), e.g.]: */
# define TYPE_GRAPHICS 1
typedef struct
{
   signed char
      rtn, /* return operation mark: 0: returm with error */
      dsp; /* display operation mark: 1 display some file saving messages */

   char
      name[STS_SIZE],
      text[STS_SIZE];

   char
      file[STS_SIZE],
      curve[GPH_CURVES][STS_SIZE],
      format[SHS_SIZE];

   char
      xunit[SHS_SIZE],
      yunit[SHS_SIZE];

   short
      nc; /* nc = number of graphics [ 'curves' ] */

   long
      nn,
      np[GPH_CURVES]; /* np = number of sample points */

   double
      xmin,
      xmax,
      ymin,
      ymax,
      vct[GPH_POINTS][GPH_CURVES+ONE];

} GRAPHICS;
/******************************* end of file grphtp.h *************************/
