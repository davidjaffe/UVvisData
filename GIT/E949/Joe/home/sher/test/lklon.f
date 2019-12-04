      REAL FUNCTION lklon()
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
      real sumlo,sum,sumhi,trlo,tr,trhi
      data sumlo,sum,sumhi/7607.255  ,   15419.00  ,   59828.86/
      data trlo,tr,trhi/0.000243,0.000516,0.00086/

*      sumlo = sumlo + lklo
*      sum = sum + lk
*      sumhi = sumhi + lkhi
*      write(6,*) br,sumlo,sum,sumhi
      
      lklon = lklo/(sumlo/(1.-trlo))
      
      END
