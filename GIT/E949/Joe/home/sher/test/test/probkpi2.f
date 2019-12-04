      function probkpi2(iflag,bkg,acc,rdev,edev,pdev)
c implicit statement goes here
      real probkpi2,rdev,edev,pdev
*
* iflag = 0   Returns number of expected kpi2 kinematic bkg in bkg, 
*              acceptance in acc
* iflag = 1   Returns the acceptance given the number of expected bkg
*              as specified by bkg.
* iflag = 2   Returns the expected bkg given the acceptance as
*              specified by acc.
*
* Acceptance and background are scaled by the acceptance and
* background levels of the 1995 TRIUMF analysis on the entire 1995
* data set. 
*
* The function value is 1 if everything is ok, 0 otherwise.
*
********************************************************************
* RTOT,ETOT AND PTOT NEED TO BE SCALED BEFORE CALLING THIS ROUTINE *
********************************************************************
*
* Modified:
*  10-Jul-98 (GR) Removed setup cuts. We now calculate a function
*                  value for any event.
*                 Modified inside-the-box function.
*                 Dvxpi cut is already applied before measuring the
*                  function.
*                 Take cos3d dependence of R,E,P peaks and sigmas
*                  into account.
*                 R-P and E cut positions are no longer optimized;
*                  we set the cut positions in advance, and just
*                  measure the performance.
*                 We now tighten R,E and P simultaneously, rather
*                  than tightening the two cuts R-P and E separately.
*                 We no longer apply the cut on nhz at the end, for
*                  fear of low statistics bias.
*  21-Sep-98 (GR) Changed ascale from 1.056 to 1.036 to take out
*                  acceptance gain that is correlated with kmu2.
*                  Removed hard cut at R=40, E=135, P=229.
*  02-Oct-98 (GR) Recalibrated with summary ntuples (w/o tracking fix
*                  and w/o empirical correction).  Added EIC cut
*                  and bad run cut when calibrating. Fixed an error
*                  in the normalization.
*  18-Oct-98 (GR) Recalibrated with summary ntuples after tracking fix.
*  22-Oct-98 (GR) "Final" recalibration after finalizing EIC and DVXPI 
*                   cuts
*  23-Oct-98 (GR) One more iteration, this time applying EIC cut properly.
*  21-Nov-98 (GR) Updated rdev,edev,pdev calculation for 96-97.
*
***************************************************************************
*
*--   Enter user code here
*
      real boxt(5,35),dnum(5),dden(5),delta(5)
      real rejt(35),acct(35),accint,rpeak,epeak,ppeak
      real rscale,ascale,nz,rscale2,bkg,acc,dacc,drej
      real esig,psig,rsig,dmin,rdevl,edevl,pdevl,rej

      integer nrowt,iflag,i,j,imin
      logical lfirst,ldebug

      data lfirst/.true./
      data ldebug/.false./
      data nrowt/35/

* Note that the numbers in the following table get scaled by
* rscale, rscale2 and ascale to get numbers relative to the 1995 TRIUMF
* analysis.

*                 range   energy   mom      Acc      Nexp
*                 ------  ------   ------   ------   ------
      data boxt/
     &1.600,1.250,1.400,1.350,1.00000,
     &1.650,1.300,1.450,1.336,0.79457,
     &1.700,1.350,1.500,1.320,0.59360,
     &1.750,1.400,1.550,1.307,0.45791,
     &1.800,1.450,1.600,1.293,0.34057,
     &1.850,1.500,1.650,1.276,0.25493,
     &1.900,1.550,1.700,1.264,0.19455,
     &1.950,1.600,1.750,1.251,0.14745,
     &2.000,1.650,1.800,1.238,0.11058,
     &2.050,1.700,1.850,1.224,0.07708,
     &2.100,1.750,1.900,1.205,0.05960,
     &2.150,1.800,1.950,1.190,0.04363,
     &2.200,1.850,2.000,1.169,0.03282,
     &2.250,1.900,2.050,1.151,0.02394,
     &2.300,1.950,2.100,1.133,0.01928,
     &2.350,2.000,2.150,1.117,0.01343,
     &2.400,2.050,2.200,1.100,0.00909,
     &2.450,2.100,2.250,1.086,0.00682,
     &2.500,2.150,2.300,1.070,0.00486,
     &2.550,2.200,2.350,1.051,0.00398,
     &2.600,2.250,2.400,1.034,0.00300,
     &2.650,2.300,2.450,1.017,0.00205,
     &2.700,2.350,2.500,1.000,0.00147,
     &2.750,2.400,2.550,0.984,0.00110,
     &2.800,2.450,2.600,0.969,0.00065,
     &2.850,2.500,2.650,0.951,0.00055,
     &2.900,2.550,2.700,0.935,0.00045,
     &2.950,2.600,2.750,0.918,0.00031,
     &3.000,2.650,2.800,0.902,0.00021,
     &3.050,2.700,2.850,0.888,0.00014,
     &3.100,2.750,2.900,0.870,0.00010,
     &3.150,2.800,2.950,0.852,0.00006,
     &3.200,2.850,3.000,0.835,0.00004,
     &3.250,2.900,3.050,0.818,0.00003,
     &3.300,2.950,3.100,0.802,0.00002/
*-- Initialization

      if(lfirst) then
        lfirst = .false.
        do i=1,nrowt
          if(i.le.14) then
            rejt(i) = boxt(5,i) 
          else
            rejt(i) = boxt(5,i) 
          endif
          acct(i) = boxt(4,i) 
        enddo
      endif

      probkpi2 = 1.

c      call repdev(rdev,edev,pdev)

      if(iflag.eq.0) then

        if(rdev.gt.1.60.and.edev.gt.1.25.and.pdev.gt.1.40) then

          j = 1
          do i=1,nrowt
            if(rdev.gt.boxt(1,i).and.
     &         edev.gt.boxt(2,i).and.
     &         pdev.gt.boxt(3,i)) then
              j = i
            endif
          enddo
          if(ldebug) write(*,*) 'Row number is ',j

          if(j.lt.nrowt) then
            dnum(1) = rdev-boxt(1,j)
            dnum(2) = edev-boxt(2,j)
            dnum(3) = pdev-boxt(3,j)

            dmin = 100.
            imin = 0
            do i=1,3
              dden(i) = boxt(i,j+1)-boxt(i,j)
              if(ldebug) then
                write(*,*) 'dnum,dden = ',dnum(i),dden(i)
              endif
              if(dden(i).ne.0.) then
                delta(i) = dnum(i)/dden(i)
              else
                delta(i) = 100.
              endif
              if(delta(i).lt.dmin) then
                dmin = delta(i)
                imin = i
              endif
            enddo
            if(ldebug) write(*,*) 'Closest quantity is ',imin

            drej = rejt(j+1)-rejt(j)
            dacc = acct(j+1)-acct(j)
            bkg = rejt(j)+dmin*drej
            acc = acct(j)+dmin*dacc
          else
            bkg = 0.
            acc = 0.
          endif
          return
        else
          bkg = 999.
          acc = 999.
          if(ldebug) write(*,*) 'R,E,P too low'
          return
        endif

      elseif(iflag.eq.1) then

        rej = bkg
        if(rej.gt.rejt(1)) then
          acc = 999.
        elseif(rej.lt.rejt(nrowt)) then
          acc = 0.
        else
          acc = accint(acct,rejt,nrowt,rej)
        endif

      elseif(iflag.eq.2) then

        if(acc.gt.acct(1)) then
          bkg = 999.
        elseif(acc.lt.acct(nrowt)) then
          bkg = 0.
        else
          bkg = accint(rejt,acct,nrowt,acc)
        endif

      else
        bkg = 999.
        acc = 999.       
        probkpi2 = 0.
        if(ldebug) write(*,*) 'Unknown value of iflag'
      endif

      return
      end
c needs : accint     
c needs : repdev
      include 'accint.f'

