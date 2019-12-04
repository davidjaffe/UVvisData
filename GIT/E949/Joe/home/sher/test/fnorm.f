      REAL FUNCTION fnorm()
      REAL
     +BR      ,LKLO    ,CLLO    ,LK      ,CL      ,LKHI    ,
     +CLHI     
*
      LOGICAL         CHAIN
      CHARACTER*128   CFILE
*
      COMMON /PAWCHN/ CHAIN, NCHEVT, ICHEVT
      COMMON /PAWCHC/ CFILE
*
      COMMON/PAWIDN/IDNEVT,OBS(13),
     +BR      ,LKLO    ,CLLO    ,LK      ,CL      ,LKHI    ,
     +CLHI     
*

      sumlo = sumlo + lklo
      sum = sum + lk
      sumhi = sumhi + lkhi

      write(6,*) br,sumlo,sum,sumhi
      
      fnorm=1.
      END
