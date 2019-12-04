      function probkmu2(ityp,iflag,bkg,acc,ptot,etot,rtot,rngmom)
c implicit statement goes here
      real probkmu2,ptot,etot,rtot,rngmom
* ityp: == 1 for tail events
*       == 2 for band events
* iflag = 0   Returns number of expected kmu2 kinematic bkg in bkg, 
*              acceptance in acc
* iflag = 1   Returns the acceptance given the number of expected bkg
*              as specified by bkg.
* iflag = 2   Returns the expected bkg given the acceptance as
*              specified by acc.
*
* The function value is 1 if everything is ok, 0 otherwise.
*
      integer ityp,iflag,ntail,nband
      parameter (ntail = 24, nband = 32)
      real bkg,acc,linint
      real tail(ntail),acc1(ntail),pcut(ntail)
      real band(nband),acc2(nband),xcut(nband)
      data tail/
     &0.4474419E-04, 
     &0.1566047E-03,0.3355814E-03,0.5145581E-03,0.7382791E-03, 
     &0.1140977E-02,0.1812140E-02,0.2975488E-02,0.4340186E-02, 
     &0.7405163E-02,0.1346800E-01,0.2356900E-01,0.3703700E-01, 
     &0.5387200E-01,0.9090900E-01,0.158249,0.276094,0.424242, 
     &0.572390,0.683501,0.794612,0.895622,0.966329, 1.000000/
      data acc1/0.346,0.426,0.511,0.598,0.684,0.761,0.838,0.897,0.944,
     &0.976,1.000,1.016,1.022,1.024,1.025,1.026,1.026,1.026,1.026,
     &1.026,1.026,1.026,1.026,1.026/
      data pcut/219.,220.,221.,222.,223.,224.,225.,226.,227.,
     &228.,229.,230.,231.,232.,233.,234.,235.,236.,237.,238.,
     &239.,240.,241.,242./ 
      data band/
     &0.7800310E-03,0.1560060E-02,0.2340090E-02,0.3900150E-02,      
     &0.3900150E-02,0.7020280E-02,0.7020280E-02,0.8580340E-02, 
     &0.1092040E-01,0.1326050E-01,0.1638060E-01,0.2106080E-01, 
     &0.2496100E-01,0.2964120E-01,0.3666150E-01,0.4368180E-01, 
     &0.5538230E-01,0.6786280E-01,0.8268340E-01,0.110765, 
     &0.163027,0.227770,0.333854,0.439938,0.564743,0.682528,
     &0.792512,0.876755,0.939157,0.979719,0.996880,1.000000/
      data acc2/
     &0.025,0.053,0.105,0.179,0.279, 
     &0.406,0.541,0.670,0.772,0.860,0.923,0.966,0.988,
     &1.000,1.007,1.010,1.012,1.013,1.014,1.014,1.014,1.014,1.014,
     &1.014,1.014,1.014,1.014,1.014,1.014,1.014,1.015,1.015/
      data xcut/
     &-1.25,-1.00,-0.75,-0.50,-0.25, 
     &0.00,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.25,
     &2.50,2.75,3.00,3.25,3.50,3.75,4.00,4.25,4.50,4.75,
     &5.00,5.25,5.50,5.75,6.00,6.25,6.50/ 

      probkmu2=1.

c The range-tail function involves cutting on ptot, and was trained
c using km2 peak events, which are defined as SKIM2 with TD inverted,
c and all cuts other than the upper-side cuts on rtot, etot, and ptot
c applied. The applied cuts include the rngmom < 2 cut which removes
c muon-band events.
c
c The muon-band function involves cutting on rngmom, and was trained
c using muon-band events, which are defined as SKIM2 with TD inverted
c and all cuts other than the PV, rngmom, and rtot cuts applied,  The
c applied cuts include a tight ptot cut at ptot < 225 which removes
c range-tail events.  The number of muon-band events "inside the box"
c which have rngmom < 2 and ptot < 229 is very small, which is why we
c can cut at ptot < 225 and exclude the ptot = [225,229] region in an
c attempt to suppress range-tail events.
c
c Range tail N values range from 0 (signal-like) to 1 (background-like),
c and are assigned to 999 if the event has ptot very background-like
c and outside the range of calibrated function values (ptot > 242).
c Muon band N values range from 0 (signal-like) to 1 (background-like),
c and are assigned to 999 if the event has rngmom very background-like
c and outside the range of calibrated function values (rngmom > 6.5).
c
c Note: there is no need to initialize bkg and acc for iflag = 0, 1,
c or 2, because bkg and acc are always defined by the following code,
c and/or are passed to this routine (bkg is passed for iflag = 1;
c acc is passed for iflag = 2).
c This routine returns a value of probkmu2 = 1 if bkg and acc are in
c some way meaningful.  probkmu2 = 0 only when there is an illegal
c use of this routine (i.e., iflag not equal to [0,1,2] or ityp not
c equal to 1 or 2).

c for tail events      
      if(ityp.eq.1)then
        if(iflag.eq.0)then
          if(rtot.gt.40.or.etot.gt.135.or.ptot.gt.pcut(ntail)) then
             bkg=999.
             acc=999.
          else if (ptot.lt.pcut(1)) then
             bkg=0.
             acc=0.
          else
             bkg = linint(tail,pcut,ntail,ptot)
             acc = linint(acc1,pcut,ntail,ptot)
          end if
        else if(iflag.eq.1)then     
          if (bkg.gt.tail(ntail)) then
             acc = 999.
          else if (bkg.lt.tail(1)) then
             acc = 0.
          else
             acc = linint(acc1,tail,ntail,bkg)
          endif
        else if(iflag.eq.2)then
          if (acc.gt.acc1(ntail)) then
             bkg = 999.
          else if (acc.lt.acc1(1)) then
             bkg = 0.
          else
             bkg = linint(tail,acc1,ntail,acc) 
          endif
        else
          probkmu2=0.
          write(6,*)'Wrong input iflag!'
          stop
        end if
c for band events
      else if(ityp.eq.2)then
        if(iflag.eq.0)then
          if(rtot.gt.40.or.etot.gt.135.or.rngmom.gt.xcut(nband)) then
             bkg=999.
             acc=999.
          else if (rngmom.lt.xcut(1)) then
             bkg=0.
             acc=0.
          else
             bkg = linint(band,xcut,nband,rngmom)
             acc = linint(acc2,xcut,nband,rngmom)
          end if
        else if(iflag.eq.1)then     
          if (bkg.gt.band(nband)) then
             acc = 999.
          else if (bkg.lt.band(1)) then
             acc = 0.
          else
             acc = linint(acc2,band,nband,bkg)
          endif
        else if(iflag.eq.2)then
          if (acc.gt.acc2(nband)) then
             bkg = 999.
          else if (acc.lt.acc2(1)) then
             bkg = 0.
          else
             bkg = linint(band,acc2,nband,acc) 
          endif
        else
          probkmu2=0.
          write(6,*)'Wrong input iflag!'
          stop
        end if
      else
        probkmu2=0.
        write(6,*)'The input ityp value does not match 1 or 2'
        stop         
      end if
      return
      end
c needs : linint      
      include 'linint.f'

