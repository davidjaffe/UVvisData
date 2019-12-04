*
      real function close_to_box(xdummy)
*
*    File:  close_to_box.f
* Purpose:  require events to fail at least one cut
*   Usage:  PAW> if (close_to_box(0.).gt.0.5) --> event is not in the box
*
* History: 2001 Jul 20  PB  created
*

      implicit none

      include 'pass3nt.inc'

      real xdummy,status,statkp2,stattd,statkm2t,statkm2b
      real scale_epr,probkpi2,tdfunc,tdfunc_patch,probkmu2
      real l2_5_5_8
      real nfpv,nfkp2,nftd,nfkm2t,nfkm2b
*      parameter (nfpv = 1.0, nfkp2 = 0.00318,
*     &           nftd = 1.0,
*     &           nfkm2t = 0.008405, nfkm2b = 0.003284)
      parameter (nfpv = 1.0, nfkp2 = 0.00147,
     &           nftd = 1.0,
     &           nfkm2t = 0.013468, nfkm2b = 0.0296412)
      real bgbm1,sigbm1,bgbm2,sigbm2
      real statbm1,statbm2
      real probbeam
      logical first
      data first /.true./

*
*  Initialize.
*
      close_to_box = 0.

*
*  Re-calculate the functions.
*
      status = scale_epr(0.)
      statkp2 = probkpi2(0,bgkp2,sigkp2,rdev,edev,pdev)
*      status = tdfunc_patch(0.)
*      stattd = tdfunc(0,bgtd,sigtd)
      statkm2t = probkmu2(1,0,bgkm2t,sigkm2t,ptot,etot,rtot,rngmom_val)
      statkm2b = probkmu2(2,0,bgkm2b,sigkm2b,ptot,etot,rtot,rngmom_val)

*
*  Test each cut class for pass/fail.
*
*  (1) PV function
*
*      if (bgpv.gt.nfpv)
*     &   return

*  (2) Kp2 kinematic function
*
      if (bgkp2.gt.1.)
     &   return

*  (3) TD function
*
*      if (bgtd.gt.nftd)
*     &   return

*  (4) Km2 kinematic function
*
      if ((bgkm2t.gt.1.).or.(bgkm2b.gt.1.))
     &    return

*  (5) fiducial cut class (track stays within detector):
*      LAYV4, COS3D, ZFRF, ZUTOUT, LAY14
*
      if (.not.(LAYV4_BIT.and.COS3D_BIT.and.ZFRF_BIT.and.
     &          ZUTOUT_BIT.and.ICODEL14_BIT))
     &   return

*  (6) scatter cut class (cuts sensitive to muon scattering in RS):
*      RSDEDXMAX, RSDEDXCL, RSLIKE, PR_RF, and PRRFZ1_N, these 5 cuts
*      combined into a neural net, PR_RF_AB, PR_RFZ2, UTCQUAL, and
*      the new UTC safety cut UTCFIT
*
      RSDEDXMAX_BIT = (rsdedxmax_val.le.4)
      RSDEDXCL_BIT = (rsdedxcl_val.ge.0.1)
      RSLIKE_BIT = ((rslike_val.le.10.).and.(rslike_val.gt.0.))
      PR_RF_BIT = (pr_rf_val.ge.0.1)
      PRRFZ1_N_BIT = (prrfz1_n_val.ge.-2.5)
      PR_RF_AB_BIT = (.not.((nhz.le.4).and.
     &                      ((pr_rfz1.lt.-7.0).or.(pr_rfz2.lt.-3.30))))
*      if (.not.((l2_5_5_8(0.).ge.0.76).and.PR_RFZ2_BIT.and.
*     &          UTCQUAL_BIT))
*     &   return
*      if (.not.(RSDEDXMAX_BIT.and.RSDEDXCL_BIT.and.RSLIKE_BIT.and.
*     &          PR_RF_BIT.and.PRRFZ1_N_BIT.and.
*     &          (l2_5_5_8(0.).ge.0.76).and.
*     &          PR_RF_AB_BIT.and.PR_RFZ2_BIT.and.
*     &          UTCQUAL_BIT.and.UTCFIT_BIT))
*     &   return
      if (.not.(RSDEDXMAX_BIT.and.RSDEDXCL_BIT.and.RSLIKE_BIT.and.
     &          PR_RF_BIT.and.PRRFZ1_N_BIT.and.
     &          (l2_5_5_8(0.).ge.0.76).and.
     &          PR_RF_AB_BIT.and.PR_RFZ2_BIT.and.
     &          UTCQUAL_BIT))
     &   return

*  (7) tgdedx cut class (target K-pi reconstruction)
*      TGDEDX, PIGAP, TGLIKE, TGB4
*
      if (.not.(TGDEDX_BIT.and.PIGAP_BIT.and.
     &          TGLIKE_BIT.and.TGB4_BIT))
     &   return

*  (8) tgepi cut class (pion-time energy in target):
*      TGCCDPF, EPITG, EPIMAXK, PHIVTX, OPSVETO, TGEDGE
*
      if (.not.(TGCCDPF_BIT.and.EPITG_BIT.and.EPIMAXK_BIT.and.
     &          PHIVTX_BIT.and.OPSVETO_BIT.and.TGEDGE_BIT))
     &   return

*  (9)  tgtrack cut class (pion tracking in target):
*       TGQUALT, TGER, TARGF, DTGTTP, RTDIF, DRP
*
      if (.not.(TGQUALT_BIT.and.TGER_BIT.and.TARGF_BIT.and.
     &          DTGTTP_BIT.and.RTDIF_BIT.and.DRP_BIT))
     &   return

*  (10) timcon cut class (timing consistency):
*       TIMCON, TIC, TGCCD
*
      if (.not.(TIMCON_BIT.and.TIC_BIT.and.TGCCD_BIT))
     &   return

*  (11) ic cut class (IC cuts):
*       EIC, KIC, TGGEO
*
      if (.not.(EIC_BIT.and.KIC_BIT.and.TGGEO_BIT))
     &   return

*  (12) b4ekz cut class (beam correlation cuts):
*       B4EKZ, TGZFOOL
*
      if (.not.(B4EKZ_BIT.and.TGZFOOL_BIT))
     &   return

*  (13) bhtrs cut class (hole counter):  BHTRS
*
      if (.not.(BHTRS_BIT))
     &   return

*  (14) Beam function cuts (used at bottom of normalization and
*       rejection branches for estimation of beam background):
*       DELC, CKTRS, CKTAIL, PBNRS, B4DEDX,
*       BWTRS, CPITRS, CPITAIL B4TRS, B4TD
*
*      if (.not.(DELC_BIT.and.CKTRS_BIT.and.CKTAIL_BIT.and.
*     &          PBNRS_BIT.and.B4DEDX_BIT.and.BWTRS_BIT.and. 
*     &          CPITRS_BIT.and.CPITAIL_BIT.and.
*     &          B4TRS_BIT.and.B4TD_BIT))
*     &   return
      if (.not.(CKTRS_BIT.and.CKTAIL_BIT.and.
     &          PBNRS_BIT.and.B4DEDX_BIT.and.
     &          CPITRS_BIT.and.CPITAIL_BIT.and.
     &          B4TD_BIT))
     &   return

      close_to_box = 1.

      if (first) then
         first = .false.
         open(71,file='close_to_box.dat',status='new')
      endif

      statbm1 = probbeam(1,0,bgbm1,sigbm1)
      statbm2 = probbeam(2,0,bgbm2,sigbm2)

      write (71,11) run,event,bgpv,bgkp2,bgtd,bgkm2t,bgkm2b,bgbm1,
     &              bgbm2,delc_val,bwtrs_val,bw1hrs,bw2hrs,utcfit_bit
 11   format(I7,I9,F9.4,G12.4,F11.4,2(G12.4),6(F10.4),I3)

      end

      include 'scale_epr.f'
      include 'probkpi2.f'
*      include 'tdfunc_patch.f'
*      include 'tdfunc.forpaw'
      include 'probkmu2.f'
      include 'l2_5_5_8.f'
      include 'probbeam.f'

