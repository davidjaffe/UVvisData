*
      real function scale_epr(xdummy)
*
*    File:  scale_epr.f
* Purpose:  scale the etot,ptot,rtot values in the ntuple
*   Usage:  status = scale_epr(0.)
*
* History:  01 Feb  PB  created, using constants from Shaomin
*

      include 'pass3nt.inc'

      real xdummy
      real ppeak,epeak,rpeak

      ppeak=205.14
      epeak=108.55
      rpeak=30.37

      if (run.ge.29000.and.run.le.33000) then
         rtot = rtot*rpeak/( 30.32-0.0333*cos3d-0.0835*cos3d**2)
         rtot = rtot - 0.07
         etot = etot*epeak/(108.82+0.0373*cos3d-0.5807*cos3d**2)
         etot = etot + 0.02
         ptot = ptot*ppeak/(204.98+0.0120*cos3d-1.3443*cos3d**2)
         ptot = ptot - 0.06
         rdev = (rtot-rpeak)/(0.8689+0.0071*cos3d+0.7309*cos3d**2)
         edev = (etot-epeak)/2.9677
         pdev = (ptot-ppeak)/(2.1927-0.0512*cos3d+1.5051*cos3d**2)
      else if(run.ge.33001.and.run.le.35500) then
         rtot = rtot*rpeak/( 30.34-0.0114*cos3d-0.1337*cos3d**2)
         rtot = rtot - 0.03
         etot = etot*epeak/(108.75+0.0374*cos3d-0.0775*cos3d**2)
         etot = etot - 0.04
         ptot = ptot*ppeak/(205.02+0.1998*cos3d-0.8461*cos3d**2)
         ptot = ptot - 0.06
         rdev = (rtot-rpeak)/(0.8508+0.0341*cos3d+0.8567*cos3d**2)
         edev = (etot-epeak)/2.9551
         pdev = (ptot-ppeak)/(2.1353-0.1105*cos3d+1.4936*cos3d**2)          
      else if(run.gt.35500) then
         rtot = rtot*rpeak/( 30.31+0.0363*cos3d-0.3741*cos3d**2)
         etot = etot*epeak/(108.41-0.1775*cos3d-1.2285*cos3d**2)
         ptot = ptot*ppeak/(205.17-0.1369*cos3d-1.1713*cos3d**2)
         ptot = ptot - 0.04 !2.3104 
         rtot = rtot - 0.06 !0.9031
         etot = etot - 0.07 !3.3845
         rdev = (rtot-rpeak)/(0.8608+0.0249*cos3d+0.7955*cos3d**2)
         edev = (etot-epeak)/3.3845
         pdev = (ptot-ppeak)/(2.1877+0.0376*cos3d+1.7232*cos3d**2) 
      endif 

      end
