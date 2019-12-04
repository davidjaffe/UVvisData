      REAL FUNCTION width
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

      vector wid(501),fcl(501),fclh(501),fcll(501),fbr(501)
      
*      sumlo = sumlo + lklo
*      sum = sum + lk
*      sumhi = sumhi + lkhi
*      write(6,*) br,cllo,cl,clhi

      i = i + 1
      wid(i) = max(cllo,cl,clhi)-min(cllo,cl,clhi)
      fclh(i) = max(cllo,cl,clhi) - cl
      fcll(i) = cl - min(cllo,cl,clhi)
      fcl(i) = cl
      fbr(i) = br

      write(6,*) i,br,fcl(i),fcll(i),fclh(i)
      
      return
      end
      








