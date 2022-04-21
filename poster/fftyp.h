/* [ file: fftyp.h ] */
/* Update: April 14, 2007 */
/*----------------------------------------------------------------------------*/
typedef struct
{                     
   signed char
      rtn;

   char
      opt[SHS_SIZE];

   short 
      p, q,
      mult[FTR_NMBR+ONE];

   long
      ttlg[FTR_NMBR+ONE], stlg[FTR_NMBR+ONE];

/* distributions, real and imaginary parts: */
   double
      r[FTR_NMBR+ONE][FTR_SIZE+ONE],
      i[FTR_NMBR+ONE][FTR_SIZE+ONE],

      t[FTR_NMBR+ONE], tt[FTR_NMBR+ONE], dt[FTR_NMBR+ONE],
      s[FTR_NMBR+ONE], ss[FTR_NMBR+ONE], ds[FTR_NMBR+ONE];

/* normalization constant: */
   double
      nor;

} FFT;
/*************************** end of file fftyp.h ******************************/
