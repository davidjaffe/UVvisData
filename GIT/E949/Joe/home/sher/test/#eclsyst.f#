c-------------------------------------------------------------------------------

c confidence level calculator:  -- using the Likelihood Ratio
c estimator.   See A. Read, Optimal Statistical Analysis of Search
c Results based on the Likelihood Ratio and its Application to the
c Search for the MSM Higgs Boson at roots=161 and 172 GeV
c DELPHI 97-158 PHYS 737

c Author:  Tom Junk,  Tom.Junk@cern.ch
c modification history:
c 2 Nov 1998 -- change -- sort the inputs by s/b when combining them.
c               this requires allocating local storage for each input
c               channel -- so give up on passing in scratch arrays and
c               have a fixed maximum channel count defined in ecl.inc
c               name:  maxbin -- increase if you need more channels.
c               note: now e_calcsb and e_calcb now need a sort list
c               in their arguments.  Note -- technical detail.  All
c               use of the SAVE statement here is simply to keep large arrays
c               from being created on the stack, causing stack overflows.
c               SAVE puts them in the heap.
c 10 Nov 1998   fix a bug in the Gaussian approximation
c 10 Nov 1998   Add e_simadd to add together the contents of similar
c               channels.
c 2 Dec  1998   Introduce conservativeness to the binning procedure -- take
c               the smallest estimator and the biggest probability within a bin.
c 18 Jan 1999   Compute expected limits more efficiently from the pdf functions
c               of the test statistic rather than with a toy MC
c 19 Jan 1999   Introduce the effect of systematic uncertainty in the expected
c               s and b -- new file -- eclsyst.f    
c 22 Jan 1999   Fix bug in computation of expected limits
c 16 Dec 1999   Add an n95 calculator
c 20 Jun 2001   Copy in e_avcls_f2 from the non-syst. version -- does
c               median and +-1,2 sigma expecteds on the median, and e_tstatmed
c               and e_tstatmed2 for additional diagnostic info which can be
c               gotten from the probability lists.
c               Also copy in the s95 expectation routines from ecl.f
c 29 Jun 2001   Fix problem in expected cls normalization: in e_avcls_f2
c               it assumed a flat 1-CLb distribution which is not true for
c               very small stats.
c------------------------------------------------------------------------------
c INPUTS
c integer n: number of channels
c real*8 s(2,n): expected signal events for each channel -- first entry: value
c                                                           second entry: error
c real* b(2,n): expected background events for each channel -- first entry: value
c                                                           second entry: error
c integer d(n): observed data candidate count 
c 
c OUTPUTS:
c clsb: probability that an experiment will have fewer candidates
c     than the actual experiment assuming that s+b events were expected. real*8
c clb   probability that an experiment will have fewer candidates
c     than the actual experiment assuming that b events were expected. real*8
c cls = clsb/clb    (real*8)
c  By convention, the confidence level for excluding a signal is 1-cls.
c  So if cls is 0.049, the s+b hypothesis is excluded at 95% over the background-only
c  hypothesis.  If cls is 0.051, it is not.
c-------------------------------------------------------------------------------

      subroutine e_cls(n,s,b,d,cls,clsb,clb)
      implicit none
      include 'ecl.inc'
      integer n,d(*)
      real*8 s(2,*),b(2,*),cls,clsb,clb
      integer isort(maxbin)
      save isort
      call e_insort(n,s,b,isort)
      call e_calcsb(n,s,b,d,clsb,isort)
      call e_calcb(n,s,b,d,clb,isort)
      if (clb.gt.1.0D-10) then
        cls = clsb/clb
      else
        cls = 1.0
      endif
      return
      end

      subroutine e_s95aux(n,s,b,d,rt,sf)
      implicit none
      include 'ecl.inc'
      integer n,d(*)
      real*8 s(2,*),b(2,*),cls
      real*8 sf
c local signal storage for scaled signals
      real*8 sl(2,maxbin)
      integer i,j
      logical foundit
      real*8 clh,cll,cla,clsf,sfh,sfl
      real*8 sflast
      external rt

      sf = 1.0
      clsf = 0.0
      cll = 0.0
      cla = 0.0
      clh = 0.0
      foundit = .false.

      do j=1,32
        if (.not.foundit) then

          do i=1,n
            sl(1,i) = s(1,i)*sf
            sl(2,i) = s(2,i)*sf
          enddo
          call rt(n,sl,b,d,cls)
c the previous version when we were only doing this for cls
c          call e_cls(n,sl,b,d,cls,clsb,clb)
          if (j.eq.1) then
            clsf = cls
            cla = cls
          endif

          if (cls.lt.0.05) then
            if (clsf.gt.0.05) then
              sfh = sf
              clh = cls
              sfl = sf/2.0
              cll = cla
              foundit = .true.
            endif
            sf = sf/2.0
          elseif (cls.gt.0.05) then
            if (clsf.lt.0.05) then
              sfl = sf
              cll = cls
              sfh = sf*2.0
              clh = cla
              foundit = .true.
            endif
            sf = sf*2.0
          else
            return
          endif
          cla = cls
        endif
      enddo

      sf = sfh
      if (.not.foundit) then
        write(*,*) 'e_s95: could not find s95'//
     >             ' within 2**32 of original guess'
        sf = 0.0
      else
        do j=1,30
          sflast = sf
          if (cll.le.0.0 .or. clh.le.0.0 .or. cll.eq.clh) then
            sf = (sfh+sfl)/2.0
          else
            sf = sfl + (log(0.05)-log(cll))*
     >                   (sfl-sfh)/(log(cll)-log(clh))
          endif
          do i=1,n
            sl(1,i) = s(1,i)*sf
            sl(2,i) = s(2,i)*sf
          enddo
          call rt(n,sl,b,d,cls)
c the previous version when we were only doing this for cls
c          call e_cls(n,sl,b,d,cls,clsb,clb)
          if (cls.gt.0.05) then
            sfl = sf
            cll = cls
          elseif (cls.lt.0.05) then
            sfh = sf
            clh = cls
          endif
          if (abs(cls-0.05).lt.0.001 .or.
     >      abs((sflast/sf) - 1.0).lt.0.001 ) then
c             write(*,*) 'N95 iterations: ',j
             return
          endif
        enddo
      endif

      return
      end


c compute "n95", or rather, the scaling factor on the siganl required
c to make cls=0.05 -- requires several calls to e_cls.
c This is done first with a factor-of-two search, and then a Newton's
c method search to get the answer.  Inputs: n,s,b,d as above.  Output: sf,
c the scale factor on s needed to get cls=0.05

c Generalize this to compute observed, median expected, and +-1,2 sigma
c 95% limit scale factors

      subroutine e_s95(n,s,b,d,sf)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*)
      real*8 sf
      external e_cls_rt
      call e_s95aux(n,s,b,d,e_cls_rt,sf)
      return
      end

c median expected s95 -- "d" is ignored here, just a dummy array
c to be passed in (so as not to waste memory inside)

      subroutine e_s95med(n,s,b,d,sf)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*)
      real*8 sf
      external e_clsmed_rt
      call e_s95aux(n,s,b,d,e_clsmed_rt,sf)
      return
      end

c compute all the s95's:  the observed, median expected, and +-1, 2 sigma
c versions.
c sf is the observed s95, and sfexp(5) are the
c +2sigma +1sigma, median, -1sigma, -2sigma expected s95 values.

      subroutine e_s95_withexpect(n,s,b,d,sf,sfexp)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*)
      real*8 sf,sfexp(5)
      external e_cls_rt,e_clsmed_rt
      external e_clsm68_rt,e_clsm95_rt
      external e_clsp68_rt,e_clsp95_rt
      call e_s95aux(n,s,b,d,e_cls_rt,sf)
      call e_s95aux(n,s,b,d,e_clsm95_rt,sfexp(1))
      call e_s95aux(n,s,b,d,e_clsm68_rt,sfexp(2))
      call e_s95aux(n,s,b,d,e_clsmed_rt,sfexp(3))
      call e_s95aux(n,s,b,d,e_clsp68_rt,sfexp(4))
      call e_s95aux(n,s,b,d,e_clsp95_rt,sfexp(5))
      return
      end


c and a little standard-format routine to be passed in.  Compute
c CLs and nothing else.
      subroutine e_cls_rt(n,s,b,d,cls)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*),cls,clsb,clb
      call e_cls(n,s,b,d,cls,clsb,clb)
      return
      end

c and for the median expected CLs
      subroutine e_clsmed_rt(n,s,b,d,clsmed)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*),cls,clsb,clb
      real*8 aclsp95,aclsp68,clsmed,aclsm68,aclsm95
      real*8 clbm
      call e_cls(n,s,b,d,cls,clsb,clb)
      call e_avcls_f2(n,s,b,
     >                aclsp95,aclsp68,clsmed,aclsm68,aclsm95,clbm)
      return
      end

c and for the median -1 sigma
      subroutine e_clsm68_rt(n,s,b,d,aclsm68)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*),cls,clsb,clb
      real*8 aclsp95,aclsp68,clsmed,aclsm68,aclsm95
      real*8 clbm
      call e_cls(n,s,b,d,cls,clsb,clb)
      call e_avcls_f2(n,s,b,
     >                aclsp95,aclsp68,clsmed,aclsm68,aclsm95,clbm)
      return
      end

c and for the median -2 sigma
      subroutine e_clsm95_rt(n,s,b,d,aclsm95)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*),cls,clsb,clb
      real*8 aclsp95,aclsp68,clsmed,aclsm68,aclsm95
      real*8 clbm
      call e_cls(n,s,b,d,cls,clsb,clb)
      call e_avcls_f2(n,s,b,
     >                aclsp95,aclsp68,clsmed,aclsm68,aclsm95,clbm)
      return
      end

c and for the median +1 sigma
      subroutine e_clsp68_rt(n,s,b,d,aclsp68)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*),cls,clsb,clb
      real*8 aclsp95,aclsp68,clsmed,aclsm68,aclsm95
      real*8 clbm
      call e_cls(n,s,b,d,cls,clsb,clb)
      call e_avcls_f2(n,s,b,
     >                aclsp95,aclsp68,clsmed,aclsm68,aclsm95,clbm)
      return
      end

c and for the median +2 sigma
      subroutine e_clsp95_rt(n,s,b,d,aclsp95)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*),cls,clsb,clb
      real*8 aclsp95,aclsp68,clsmed,aclsm68,aclsm95
      real*8 clbm
      call e_cls(n,s,b,d,cls,clsb,clb)
      call e_avcls_f2(n,s,b,
     >                aclsp95,aclsp68,clsmed,aclsm68,aclsm95,clbm)
      return
      end


c calc confidence for s+b.  Note -- requires now a sort order isort.
c this is provided by the call from e_cls.

      subroutine e_calcsb(n,s,b,d,cl,isort)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*),cl
      integer isort(*)
      call e_calc(n,s,b,d,cl,.false.,isort)
      return
      end

c calc confidence for b.  Note -- requires now a sort order isort.
c this is provided by the call from e_cls.

      subroutine e_calcb(n,s,b,d,cl,isort)
      implicit none
      integer n,d(*)
      real*8 s(2,*),b(2,*),cl
      integer isort(*)
      call e_calc(n,s,b,d,cl,.true.,isort)
      return
      end


c---------------------------------------------------------------
c    Create a sort order for the input array -- put the
c    ones with the best s/b at the end of the list

      subroutine e_insort(n,s,b,isort)
      implicit none
      include 'ecl.inc'
      integer n,isort(*)
      real*8 s(2,*),b(2,*)
      real*8 sdb(maxbin)
      integer i
      do i=1,n
        if (b(1,i).gt.0.0D0) then
          sdb(i) = s(1,i)/b(1,i)
        else
          sdb(i) = 1.0E20
        endif
      enddo
      call sortzv(sdb,isort,n,1,0,0)
      return
      end

c---------------------------------------------------------------
c    compute a s+b or b confidence level
      
      subroutine e_calc(n,s,b,d,cl,bgonly,isort)
      implicit none
      include 'ecl.inc'
      include 'eclplist.inc'
      include 'eclout.inc'
      integer n,d(*),isort(*)
      real*8 s(2,*),b(2,*),cl,x
      logical bgonly
      integer i,nlist,nacc
      real*8 accprob(maxacc),accest(maxacc)
      real*8 prob(maxlist),est(maxlist)
      real*8 dest
      real*8 zeros(2)
      data zeros /0.0D0,0.0D0/

      cl = 0.0D0
      if (check_inputs) call e_check(n,s,b,d)

c start with a null channel and then combine the others one by one onto it.

      call e_estlist(zeros,zeros,accprob,accest,nacc,bgonly)
      do i=1,n
        call e_estlist(s(1,isort(i)),
     >                 b(1,isort(i)),
     >                 prob,est,nlist,bgonly)
        call e_combine(accprob,accest,nacc,prob,est,nlist)
      enddo

c save the pdf of the estimator for future use (see e_avcls_f below)

      if (bgonly) then
        naccb = nacc
        do i=1,nacc
          baccp(i) = accprob(i)
          baccest(i) = accest(i)
        enddo
      else
        naccsb = nacc
        do i=1,nacc
          sbaccp(i) = accprob(i)
          sbaccest(i) = accest(i)
        enddo
      endif

c find the estimator corresponding to the observed data point.

      call e_dest(n,s,b,d,dest,x)
      xobs = exp(x)

c find all points in the list with an estimator bigger than
c the data estimator and report the probability to see an experiment
c with a smaller estimator as the confidence level.
c (sometimes you get some points in the list with equal values of the
c estimator)

      cl = 0.0D0
      do i=1,nacc
        if (accest(i).le.dest) then
          cl = cl + accprob(i)
        endif
      enddo
      if (cl.gt.0.5D0) then
        cl = 0.0D0
        do i=1,nacc
          if (accest(i).ge.dest) then
            cl = cl + accprob(i)
          endif
        enddo
        cl = 1 - cl
      endif

 10   continue

      return
      end

c----------------------------------------------------------------------
c compute average expected confidence levels
c need only do the combination once -- then read off the CL's
c from the tables in little toy runs 

c INPUTS
c n: number of channels
c s: expected signal events for each channel (real*8, of dimension n)
c b: expected background events for each channel (real*8, of dimension n)
c OUTPUTS:
c aclsb: average expected cls if many experiments were run, and only background
c        was genuinely present. real*8
c aclb   average expected clb if many experiments were run, and only background
c        was genuinely present.  Should be equal to 0.5 except in degenerate cases,
c        such as zero expected background.
c acls = <clsb/clb>    (real*8)

c simplify this because we often need to calculate both confidence levels and
c expected ones -- do not need to recombine all the experimental inputs in
c this case.

      subroutine e_avcls(n,s,b,acls,aclsb,aclb)
      implicit none
      include 'ecl.inc'
      integer n
      real*8 s(2,*),b(2,*),acls,aclsb,aclb
      real*8 cls,clsb,clb
      integer iwork(maxbin)
      save iwork

c data is unimportant here, we just need to compute the pdf's of the
c estimators assuming s and b -- this is done in e_cls

      call vzero(iwork,n)
      call e_cls(n,s,b,iwork,cls,clsb,clb)
      call e_avcls_f(n,s,b,acls,aclsb,aclb)
      return
      end

c-----------------------------------------------------
c e_avcls_f is the same thing as e_avcls except it assumes that
c e_cls has already been called, filling in the pdf of the
c estimator for signal and signal+background

      subroutine e_avcls_f(n,s,b,acls,aclsb,aclb)
      implicit none
      include 'ecl.inc'
      include 'eclplist.inc'
      real*8 s(2,*),b(2,*),acls,aclsb,aclb
      integer i,j,n
      real*8 psumb,psumsb

      aclb = 0.0D0
      aclsb = 0.0D0
      acls = 0.0D0

c compute aclb, aclsb, and acls analytically

      psumb = 0.0D0
      psumsb = 0.0D0
      j = 1
      do i=1,naccb
        psumb = psumb + baccp(i)
 10     continue
        if (baccest(i).ge.sbaccest(j)) then
          psumsb = psumsb + sbaccp(j)
          if (j.lt.naccsb) then
            j = j + 1
            goto 10
          endif
        endif
        if (baccp(i).ne.0.0D0) then
          aclb = aclb + baccp(i)*psumb
          aclsb = aclsb + baccp(i)*psumsb
          acls = acls + baccp(i)*psumsb/psumb
        endif
      enddo

      return
      end


c-----------------------------------------------------
c e_avcls_f2 computes median and +- 1 and 2 sigma CLs's.

      subroutine e_avcls_f2(n,s,b,
     >   aclsp95,aclsp68,aclsm,aclsm68,aclsm95,clbm)
      implicit none
      include 'ecl.inc'
      include 'eclplist.inc'
      real*8 s(2,*),b(2,*),aclsp95,aclsp68,aclsm,aclsm68,aclsm95
      integer i,j,n
      real*8 psumb,psumsb
      real*8 clbm

      aclsp95 = 0.0
      aclsp68 = 0.0
      aclsm = 0.0
      aclsm68 = 0.0
      aclsm95 = 0.0
      clbm = 0.0

      psumb = 0.0
      do i=1,naccb
        psumb = psumb + baccp(i)
        if (psumb.gt.0.025.and.aclsp95.eq.0.0) then
           do j=1,naccsb
             if (baccest(i).ge.sbaccest(j)) then
               aclsp95 = aclsp95 + sbaccp(j)/psumb
             endif
           enddo
        endif
        if (psumb.gt.0.16.and.aclsp68.eq.0.0) then
           do j=1,naccsb
             if (baccest(i).ge.sbaccest(j)) then
               aclsp68 = aclsp68 + sbaccp(j)/psumb
             endif
           enddo
        endif
        if (psumb.gt.0.5.and.aclsm.eq.0.0) then
           do j=1,naccsb
             if (baccest(i).ge.sbaccest(j)) then
               aclsm = aclsm + sbaccp(j)/psumb
             endif
           enddo
        endif
        if (psumb.gt.0.84.and.aclsm68.eq.0.0) then
           do j=1,naccsb
             if (baccest(i).ge.sbaccest(j)) then
               aclsm68 = aclsm68 + sbaccp(j)/psumb
             endif
           enddo
        endif
        if (psumb.gt.0.975.and.aclsm95.eq.0.0) then
           do j=1,naccsb
             if (baccest(i).ge.sbaccest(j)) then
               aclsm95 = aclsm95 + sbaccp(j)/psumb
             endif
           enddo
        endif
      enddo

      psumsb = 0.0
      do i=1,naccsb
        psumsb = psumsb + sbaccp(i)
        if (psumsb.gt.0.50.and.clbm.eq.0.0) then
          do j=1,naccb
            if (baccest(j).ge.sbaccest(i)) then
              clbm = clbm + baccp(j)
            endif
          enddo
        endif
      enddo

      return
      end


c-----------------------------------------------------
c e_avclb_lr -- compute median and +-1, 2 sigma expectations on 1-CLb
c for the signal hypothesis.  You must have called e_cls first in order
c to have all the PDFs computed, and to have done it with the likelihood
c ratio test statistic (just the sum of d*log(1+s/b)'s).  Other test statistics
c will not work with this routine as it uses a property of the likelihood ratio.
c 4 July 2001 -- moved to eclsyst 18 July 2001

      subroutine e_avclb_lr(n,s,b,asigclb)
      implicit none
      include 'ecl.inc'
      include 'eclplist.inc'
      real*8 s(2,*),b(2,*),asigclb(5)
      integer i,j,n
      real*8 psumsb,ssum

c to compute the full LR, need the signal sum..

      ssum = 0.0D0
      do i=1,n
        ssum = ssum + s(1,i)
      enddo

      do i=1,5
        asigclb(i) = 0.0D0
      enddo
      psumsb = 0.0D0

      do i=1,naccsb
        psumsb = psumsb + sbaccp(i)
        if (psumsb.gt.0.025D0.and.asigclb(1).eq.0.0D0) then
           do j=i,naccsb
             asigclb(1) = asigclb(1)+sbaccp(j)*exp(ssum-sbaccest(j))
           enddo
        endif
        if (psumsb.gt.0.16D0.and.asigclb(2).eq.0.0D0) then
           do j=i,naccsb
             asigclb(2) = asigclb(2)+sbaccp(j)*exp(ssum-sbaccest(j))
           enddo
        endif
        if (psumsb.gt.0.5D0.and.asigclb(3).eq.0.0D0) then
           do j=i,naccsb
             asigclb(3) = asigclb(3)+sbaccp(j)*exp(ssum-sbaccest(j))
           enddo
        endif
        if (psumsb.gt.0.84D0.and.asigclb(4).eq.0.0D0) then
           do j=i,naccsb
             asigclb(4) = asigclb(4)+sbaccp(j)*exp(ssum-sbaccest(j))
           enddo
        endif
        if (psumsb.gt.0.975D0.and.asigclb(5).eq.0.0D0) then
           do j=i,naccsb
             asigclb(5) = asigclb(5)+sbaccp(j)*exp(ssum-sbaccest(j))
           enddo
        endif
      enddo

      return
      end

c---------------------------------------------------------------------
c a better attempt at computing small 1-CLb's:
c still some limitations:  One must use the likelihood ratio and call e_cls
c to compute the pdf's before calling this routine.  And if 1-CLb is near 1,
c it may not be accurate (deficit situation).  This routine uses the signal
c PDF multiplied by the LR to get a better computation of the bg PDF.  So if
c the data are much much more signal-like than even the bulk of the signal, then
c the computation may also lose a bit of precision -- about 1% of the signal
c PDF should be more signal-like than the data to get a good computation of
c 1-CLb (this is an estimate -- not rigorously defined).
c inputs: n, s, b, d -- #bins, signal, background, data.  
c output: clbm1 -- 1-clb  (double precision)
c 18 July 2001

      subroutine e_clb_lr(n,s,b,d,clbm1)
      implicit none
      include 'ecl.inc'
      include 'eclplist.inc'
      real*8 s(2,*),b(2,*),clbm1,x
      integer i,n,d(*)
      real*8 ssum,dest

      ssum = 0.0D0
      do i=1,n
        ssum = ssum + s(1,i)
      enddo
      call e_dest(n,s,b,d,dest,x)
      clbm1 = 0.0D0

      do i=1,naccsb
        if (sbaccest(i).gt.dest) then
          clbm1 = clbm1 + sbaccp(i)*exp(ssum - sbaccest(i))
        endif
      enddo

      return
      end


c----------------------------------------------------------------------
c compute the median test statistic in the bg hypothesis, as well as
c the median in the signal hypothesis and +- 1 and 2 sigma on it.

      subroutine e_tstatmed(tsbmed,tss)
      implicit none
      include 'ecl.inc'
      include 'eclplist.inc'
      real*8 tsbmed,tss(5)
      integer i
      real*8 psumb,psumsb
      tsbmed = 0.0
      do i=1,5
        tss(i) = 0.0
      enddo

      psumsb = 0.0
      do i=1,naccsb
        psumsb = psumsb + sbaccp(i)
        if (psumsb.gt.0.025.and.tss(1).eq.0.0) then
          tss(1) = sbaccest(i)
        endif
        if (psumsb.gt.0.16.and.tss(2).eq.0.0) then
          tss(2) = sbaccest(i)
        endif
        if (psumsb.gt.0.5.and.tss(3).eq.0.0) then
          tss(3) = sbaccest(i)
        endif
        if (psumsb.gt.0.84.and.tss(4).eq.0.0) then
          tss(4) = sbaccest(i)
        endif
        if (psumsb.gt.0.975.and.tss(5).eq.0.0) then
          tss(5) = sbaccest(i)
        endif
      enddo

      psumb = 0.0
      do i=1,naccb
        psumb = psumb + baccp(i)
        if (psumb.gt.0.50.and.tsbmed.eq.0.0) then
          tsbmed = baccest(i)
        endif
      enddo

      return
      end

c----------------------------------------------------------------------
c compute the median test statistic in the background hypothsis
c as well as +-1, +-2 sigma on it, as well as
c the median in the signal hypothesis and +- 1 and 2 sigma on it.
c tsb(1:5) = -2, -1, median, +1, +2 on the background tstat
c      (sum log(1+s/b))
c tss(1:5) = -2, -1, median, +1, +2 on the background tstat

      subroutine e_tstatmed2(tsb,tss)
      implicit none
      include 'ecl.inc'
      include 'eclplist.inc'
      real*8 tsb(5),tss(5)
      integer i
      real*8 psumb,psumsb
      do i=1,5
        tss(i) = 0.0
        tsb(i) = 0.0
      enddo

      psumsb = 0.0
      do i=1,naccsb
        psumsb = psumsb + sbaccp(i)
        if (psumsb.gt.0.025.and.tss(1).eq.0.0) then
          tss(1) = sbaccest(i)
        endif
        if (psumsb.gt.0.16.and.tss(2).eq.0.0) then
          tss(2) = sbaccest(i)
        endif
        if (psumsb.gt.0.5.and.tss(3).eq.0.0) then
          tss(3) = sbaccest(i)
        endif
        if (psumsb.gt.0.84.and.tss(4).eq.0.0) then
          tss(4) = sbaccest(i)
        endif
        if (psumsb.gt.0.975.and.tss(5).eq.0.0) then
          tss(5) = sbaccest(i)
        endif
      enddo

      psumb = 0.0
      do i=1,naccb
        psumb = psumb + baccp(i)
        if (psumb.gt.0.025.and.tsb(1).eq.0.0) then
          tsb(1) = baccest(i)
        endif
        if (psumb.gt.0.16.and.tsb(2).eq.0.0) then
          tsb(2) = baccest(i)
        endif
        if (psumb.gt.0.5.and.tsb(3).eq.0.0) then
          tsb(3) = baccest(i)
        endif
        if (psumb.gt.0.84.and.tsb(4).eq.0.0) then
          tsb(4) = baccest(i)
        endif
        if (psumb.gt.0.975.and.tsb(5).eq.0.0) then
          tsb(5) = baccest(i)
        endif
      enddo

      return
      end

c----------------------------------------------------------------------
c    ESTIMATOR COMPUTATION

c compute the estimator function for the observed distribution
c of candidates in the data d (a list for all channels).

      subroutine e_dest(n,s,b,d,dest,x)
      implicit none
      include 'ecl.inc'
      integer n,d(*)
      real*8 s(2,*),b(2,*),dest,x
      real*8 estim,xi
      integer i
      dest = 0.0D0
      x = 0.0D0
      do i=1,n
        call e_estimator(s(1,i),b(1,i),d(i),estim,xi)
        if (estim.ge.estbig) then
          dest = estbig
          x = estbig
        else
          dest = dest + estim
          x = x + xi
        endif
      enddo
      return
      end

c single channel test statistic computer -- compute for all possible
c outcomes 

      subroutine e_estlist(s,b,prob,est,nlist,bgonly)
      implicit none
      include 'ecl.inc'
      real*8 s(2),b(2),prob(maxlist),est(maxlist),xi
      real*8 psum,x,dx,xlo,xhi,cvm
      integer nlist
      integer i
      logical bgonly

c start the list at 1 but the event count starts at zero.
c doesn't real*8ly matter to the combiner because all it wants
c is the probability of each outcome, however indexed.
c put a cap on the end for the integral off to infinity.

c special case -- expect no events (channel switched off)

      if (s(1).le.0.0D0.and.b(1).le.0.0D0) then
        nlist = 1
        prob(1) = 1.0D0
        est(1) = 0.0D0
        return
      endif

      if (bgonly) then
        cvm = b(1)
      else
        cvm = s(1) + b(1)
      endif

c if we expect lots and lots of events then produce something
c with a gaussian estimated prob. distribution. -- Use only 30
c bins for this out to +-5 sigma for now.  Keep computation time down.
c sqrt(30) is a little more than 5, so going 5 sigma below average
c won't take us negative.

      if (cvm.gt.30) then
        xlo = cvm - 5.*sqrt(cvm)
        xhi = cvm + 5.*sqrt(cvm)
        dx = (xhi - xlo)/float(30)
        psum = 0.0D0
        do i=1,30
          x = xlo + dx*(i-1)
          call e_poiss_e(s,b,bgonly,nint(x),prob(i))
          call e_estimator(s,b,nint(x),est(i),xi)
          psum = psum + prob(i)
        enddo
        do i=1,30
          prob(i) = prob(i)/psum
        enddo
        nlist = 30
        return
      endif

c otherwise do an explicit listing of all possible final outcomes

      psum = 0.0D0
      do i=1,maxlist-1
        nlist = i
        call e_estimator(s,b,i-1,est(i),xi)
        if (est(i).ge.estbig) then
          prob(i) = 1.0D0-psum
        else
          call e_poiss_e(s,b,bgonly,i-1,prob(i))
        endif
        psum = psum + prob(i)
        if (psum.gt.pmaxcut) goto 10
      enddo
 10   continue
 
      if (psum.lt.1.0D0) then
        nlist = nlist + 1
        est(nlist) = est(nlist-1)+1.0D0
        prob(nlist) = 1.0D0 - psum
      elseif (psum.gt.1.00000001D0) then
        write(*,*) 'Invalid psum in e_estlist: ',psum
        stop
      endif

      return
      end

c single-channel estimator with systematic uncertainties
      subroutine e_estimator(s,b,d,estim,xi)
      implicit none
      include 'ecl.inc'
      integer d
      real*8 s(2),b(2),pi,estim,xi
      real*8 dst,dbt,dst1,dbt1,dsum,ssum
      real*8 bsig,ssig,stmp,btmp
      parameter (pi=3.14159265358979)

c handle gracefully cases with zero error also.  Do others numerically
c within +-3 sigma

      if (s(2).le.0.0D0 .and. b(2).le.0.0D0) then
        call e_estimatorn(s(1),b(1),d,estim,xi)
      elseif (s(2).le.0.0D0) then
        dbt1 = 1.0D0/(sqrt(2.0D0*pi)*b(2))
        dsum = 0.0D0
        ssum = 0.0D0
        stmp = s(1)
        do bsig=-systwid,systwid,dsystwid
          btmp = b(1) + bsig*b(2)
          if (btmp.gt.0.0D0) then
            dbt = dbt1*exp(-bsig**2/2.0D0)
            dsum = dsum + dbt
            call e_estimatorn(stmp,btmp,d,estim,xi)
            ssum = ssum + estim*dbt
          endif
        enddo
        estim = ssum/dsum
      elseif (b(2).le.0.0D0) then
        dst1 = 1.0D0/(sqrt(2.0D0*pi)*s(2))
        dsum = 0.0D0
        ssum = 0.0D0
        btmp = b(1)
        do ssig=-systwid,systwid,dsystwid
          stmp = s(1) + ssig*s(2)
          if (stmp.gt.0.0D0) then
            dst = dst1*exp(-ssig**2/2.0D0)
            dsum = dsum + dst
            call e_estimatorn(stmp,btmp,d,estim,xi)
            ssum = ssum + estim*dst
          endif
        enddo
        estim = ssum/dsum
      else
        dst1 = 1.0D0/(sqrt(2.0D0*pi)*s(2))
        dbt1 = 1.0D0/(sqrt(2.0D0*pi)*b(2))
        dsum = 0.0D0
        ssum = 0.0D0
        do ssig=-systwid,systwid,dsystwid
          stmp = s(1) + ssig*s(2)
          if (stmp.ge.0.0D0) then
            dst = dst1*exp(-ssig**2/2.0D0)
            do bsig=-systwid,systwid,dsystwid
              btmp = b(1) + bsig*b(2)
              if (btmp.gt.0.0D0) then
                dbt = dbt1*exp(-bsig**2/2.0D0)
                dsum = dsum + dst*dbt
                call e_estimatorn(stmp,btmp,d,estim,xi)
                ssum = ssum + estim*dst*dbt
              endif
            enddo
          endif
        enddo
        estim = ssum/dsum
      endif

      return
      end


c--------------------------------------------------------------------------
c single channel estimator -- change this subroutine
c to change the overall estimator
c function
c IMPORTANT -- estimators are now considered additive in the calling
c routines -- the estimator for a combination of channels is the
c sum of the estimators of each channel.  This means that if you have
c a multiplicative estimator, please use the logarithm!
c (this makes life easier for testing other additive estimators)

c 13 November 1998 TRJ -- added in ALEPH and L3 estimator functions
c note -- this subroutine does not handle systematic uncertainty --
c e_estimator calls it many times to estimate the average value of
c the estimator given systematic uncertainties on s and b

      subroutine e_estimatorn(s,b,d,estim,xi)
      implicit none
      include 'ecl.inc'
      real*8 s,b,estim,xi
      integer d
      real*8 e_factorial,e_psi,e_pbinom
      external e_factorial,e_psi,e_pbinom
      real*8 prob
      logical pbock,lratio,l3estimator,alephestim
      parameter (pbock       = .false.)
      parameter (lratio      = .true. )
      parameter (l3estimator = .false.)
      parameter (alephestim  = .false.)
      real*8 cpbock
      parameter (cpbock = 1.0D0)
      real*8 sumn,sumd,sbp,bp
      integer i

      if (pbock) then
        if (s.eq.0.0) then
          estim = 0.0D0
        else
          estim = dble(d)/(cpbock + b/s)
        endif

      elseif (lratio) then

c the ratio of Poisson probabilities of s and b -- can be simplified,
c see below.
c        estim = 1.0
c        call e_poiss(s+b,d,prob)
c        estim = estim*prob
c        call e_poiss(b,d,prob)
c        if (prob.eq.0.0D0) then
c          estim = estbig
c        else
c          estim = log(estim/prob)
c        endif

        if (b.le.0.0D0 .and. d.gt.0) then
          estim = estbig
          xi = estbig
        elseif (b.eq.0.0D0.and.s.eq.0.0D0) then
          estim = 0.0D0
          xi = 0.0D0
        elseif (b.eq.0.0D0.and.d.eq.0) then
          estim = 0.0D0
          xi = 0.0D0
        else
          estim = dble(d)*log(1.0 + s/b)
          xi = -s + dble(d)*log(1.0 + s/b)
        endif

      elseif (l3estimator) then
        if (b.le.0.0D0 .and. d.gt.0) then
          estim = estbig
        elseif (b.eq.0.0D0.and.s.eq.0.0D0) then
          estim = 0.0D0

c this case is defined fine, it is just separated here in order
c not to divide zero by zero.

        elseif (b.eq.0.0D0.and.s.gt.0.0D0) then
          sumd = e_factorial(d)
          sumn = 0.0
          sbp = (s+b)**d
          do i=d,0,-1
            sumn = sumn + sbp
            sbp = sbp*dble(i)/(s+b)
          enddo
          estim = log(sumn/sumd) - s
        else
          sumn = 0.0
          sumd = 0.0
          sbp = (s+b)**d
          bp = b**d
          do i=d,0,-1
            sumn = sumn + sbp
            sbp = sbp*dble(i)/(s+b)
            sumd = sumd + bp
            bp = bp*dble(i)/b
          enddo
          estim = log(sumn/sumd) - s
        endif

      elseif (alephestim) then
        if (b.le.0.0D0 .and. d.gt.0) then
          estim = estbig
        elseif (s.eq.0.0D0) then
          estim = 0.0D0
        else
          sumn = 0.0D0
          do i=0,d
            call e_poiss(s+b,i,prob)
            sumn = sumn + prob
c            sumn = sumn + prob*(s/(s+b))**i
c            sumn = sumn + prob*e_psi(i,(s/(s+b))**i)
c this one real*8ly doesn't work
c            sumn = sumn + prob*e_pbinom(i,d,(s/(s+b)))
          enddo
          estim = log(sumn)
        endif
      endif

      return
      end

c an auxiliary function for the ALEPH estimator (maybe do not need it now
c because we bin the inputs)

      real*8 function e_psi(k,z)
      implicit none
      integer k
      real*8 z
      real*8 x,sum,term
      integer j

      x = -log(z)
      term = 1.0D0
      sum = 1.0D0
      do j=1,k-1
        term = term*x/dble(j)
        sum = sum + term
      enddo
      e_psi = z*sum
      return
      end

c another aux. function for ALEPH estimator tests

      real*8 function e_pbinom(j,n,f)
      integer j,n
      real*8 f
      real*8 sum,e_factorial
      external e_factorial
      sum = 0.0
      do i=j,j
        sum = sum + e_factorial(n)/(e_factorial(i)*e_factorial(n-i))*
     >              (f**i)*(1.0D0-f)**(n-i)
      enddo
      e_pbinom = sum
      return
      end

c-------------------------------------------------------------------------

c   COMBINATION

c combine two channels -- put the results in the first one.
c keep a sparseish list of the
c estimators and probabilities so the storage and CPU are manageable.

      subroutine e_combine(accprob,accest,nacc,prob,est,nlist)
      implicit none
      include 'ecl.inc'
      real*8 accprob(maxacc),accest(maxacc)
      real*8 prob(maxlist),est(maxlist)
      real*8 est12(maxlist*maxacc)
      real*8 accp(maxacc),accestl(maxacc),psum
      integer nacc,nlist
      integer i,j,k,ix1,ix2,ibin,naccs
      integer iso(maxacc*maxlist),ix(maxacc)
      save est12,iso,accp,accestl,ix

c add up the probabilities to see less than or equal to
c the estimator at each particular point.  Do not do a 2D
c sum here -- put it all in a linear array, sort it, and
c do linear sums -- take advantage of all the sums done before.

c put a ceiling on big estimators -- effectively an "infinity"
c which can be expressed as a floating point number for later
c sorting

c the save statements above are just to keep the big arrays out
c of the stack.

      k = 0
      do i=1,nlist
        do j=1,nacc
          k = k + 1
          if (est(i).ge.estbig.or.accest(j).ge.estbig) then
            est12(k) = estbig
          else
            est12(k) = min(estbig,est(i)+accest(j))
          endif
        enddo
      enddo
      call sortzv(est12,iso,k,1,0,0)

c bin the possible final outcomes by their cumulative probability
c conservative approach -- take the largest accumulated probability
c in a bin and the smallest test-statistic.

      psum = 0.0D0
      call vzero(ix,maxacc)
      naccs = nacc
      nacc = 0
      do i=1,k
        ix1 = 1 + (iso(i)-1)/naccs
        ix2 = iso(i) - naccs*(ix1-1)
        psum = psum + prob(ix1)*accprob(ix2)
        call e_probbin(psum,ibin)
        if (ix(ibin).eq.0) then
          ix(ibin) = 1
          nacc = nacc + 1
          accestl(nacc) = est(ix1) + accest(ix2)
        endif
        accp(nacc) = psum
      enddo

c differential probability -- get from the sums
      accprob(1) = accp(1)
      do i=2,nacc
        accprob(i) = accp(i)-accp(i-1)
      enddo
      do i=1,nacc
        accest(i) = accestl(i)
      enddo
      return
      end

c probability bin calculator -- make sure there is plenty of coverage
c at the small probability end -- some people make cl curves with
c probabilities all the way down to 1E-7  -- let's go down to 1E-10
c cover the small probabilities in logarithmic steps and larger ones
c (bigger than 1%) linearly.
c one linear section -- from 0.01 to 1.0 so as to cover 5% a little
c better - always get within 0.1% of 5%. 
c Also bin 1062 is populated by
c prob .eq. 1.0 identically.  Be sure that there is always a bin reserved
c for exactly zero and one reserved for exactly one.

      subroutine e_probbin(prob,ibin)
      implicit none
      include 'ecl.inc'
      real*8 prob,probl
      integer ibin

      if (prob.lt.0.0D0) then
        write(*,*) 'Negative probability sent to e_probbin ',prob
        stop
      elseif (prob.gt.1.00000000001D0) then
        write(*,*) 'Probability bigger than 1 sent to e_probbin ',prob
        stop
      elseif (prob.eq.0.0D0) then
        ibin = 1
      else
        probl = log10(prob)
        if (probl.lt.-10.0D0) then
          ibin = 1
        elseif (probl.lt.-2.0D0) then
          ibin = 2 + int(20.0D0*(probl+10.0D0))
c          ibin = 2 + int(50.0D0*(probl+10.0D0))
        else
          ibin = 162 + int(scaleplin*(prob-0.01D0))        
c          ibin = 402 + int(scaleplin*(prob-0.01D0))        
        endif
      endif

      return
      end

c---------------------------------------------------------
c   CHECKING INPUTS

      subroutine e_check(n,s,b,d)
      implicit none
      include 'ecl.inc'
      integer n,d(*)
      real*8 s(2,*),b(2,*)
      integer i

      if (n.gt.maxbin) then
        write(*,*) 'Channel count ',n,
     >    ' Exceeds internal limit: ',maxbin
        write(*,*) 'Please increase maxbin in ecl.inc and recompile.'
        write(*,*) 'Fatal error.'
        stop
      endif

      do i=1,n
        if (b(1,i).lt.0.0D0) then
          write(*,*) 'Negative expected background: ',i,b(1,i)
        endif
        if (s(1,i).lt.0.0D0) then
          write(*,*) 'Negative expected signal: ',i,s(1,i)
        endif
        if (d(i).lt.0) then
          write(*,*) 'Negative candidate count: ',i,d(i)
        endif
        if (b(1,i).eq.0.0D0.and.d(i).gt.0.and.b(2,i).le.0.0D0) then
          write(*,*) 'No BG or error but have a candidate ',
     >                i,b(1,i),'+-',b(2,i),d(i)
        endif
c        if (b(1,i)+s(1,i).gt.130.) then
c          write(*,*) 'Very large expected s+b -- caution so far'
c          write(*,*) ' i,s,b: ',i,s(1,i),b(1,i)
c        endif
      enddo

      return
      end

c --------------------------------------------------------
c add together contents of channels with similar s/b
c warning -- this preprocessing stage writes over the input arrays.
c to be called by the program assembling the inputs as an option to
c speed up the combination computation

      subroutine e_simadd(n,s,b,d)
      include 'ecl.inc'
      integer n,d(*)
      real*8 s(2,*),b(2,*)
      integer isort(maxbin)
      real*8 rtmp(maxbin),x,xmin
      integer itmp(maxbin)
      integer i,nb

c sort the inputs and put them in the sort order

      call e_insort(n,s,b,isort)
      do j=1,2
        do i=1,n
          rtmp(i) = s(j,isort(i))
        enddo
        do i=1,n
          s(j,i) = rtmp(i)
        enddo
        do i=1,n
          rtmp(i) = b(j,isort(i))
        enddo
        do i=1,n
          b(j,i) = rtmp(i)
        enddo
      enddo

      do i=1,n
        itmp(i) = d(isort(i))
      enddo
      do i=1,n
        d(i) = itmp(i)
      enddo

c collect inputs -- step through the array -- nb
c will never step faster than i, so can write over
c s, b, and d at nb.

      nb = 1
      call e_simadda(s(1,1),b(1,1),xmin)
      do i=2,n
        call e_simadda(s(1,i),b(1,i),x)
        if (abs(x-xmin).gt.addwid) then
          nb = nb + 1
          s(1,nb) = s(1,i)
          b(1,nb) = b(1,i)
          s(2,nb) = s(2,i)
          b(2,nb) = b(2,i)
          d(nb) = d(i)
          xmin = x
        else
          s(1,nb) = s(1,nb) + s(1,i)
          b(1,nb) = b(1,nb) + b(1,i)
          s(2,nb) = sqrt(s(2,nb)**2 + s(2,i)**2)
          b(2,nb) = sqrt(b(2,nb)**2 + b(2,i)**2)
          d(nb) = d(nb) + d(i)
        endif
      enddo
      n = nb

      return
      end

c something to group by (no systematic errors here -- they get combined
c above in e_simadd)

      subroutine e_simadda(s,b,x)
      real*8 s,b,x
      if (b.le.0.0D0) then
        x = 1.0D20
      elseif (s.le.0.0D0) then
        x = 0.0D0
      else
        x = log(1.0D0 + s/b)
      endif
      return
      end

c utilities ----------------------------------------------
c Poisson probabilities and factorials.

      subroutine e_poiss(u,n,prob)
      implicit none
      integer n
      real*8 u,prob,e_factorial,pi
      external e_factorial
      parameter (pi=3.14159265358979)
c make sure 0^0=1 in this case.  Zero expected events gives 100%
c probability to zero observed events.
      if (n.lt.100) then
        if (u.lt.1.0D-6.and.n.eq.0) then
          prob = 1.0D0
        else
          prob = (u**n)*exp(-u)/e_factorial(n)
        endif
      else
        prob = exp(-(dble(n)-u)**2/(2.0D0*u))/sqrt(2.0D0*pi*u)
      endif
      return
      end


c Poisson probability with a systematic uncertainty on s and b
c Do the calculation numerically --
c calculation courtesy of C. Giunti, hep-ex/9901015 -- There is an
c exact calculation for one syst. error, but if there are two that can be
c big, do it numerically.  Maybe do this with adaptive quadrature later,
c but for now, this ought to do.

      subroutine e_poiss_e(s,b,bgonly,n,prob)
      implicit none
      include 'ecl.inc'
      integer n,i
      logical bgonly
      real*8 s(2),b(2),prob,e_factorial,pi
      real*8 dst,dbt,dst1,dbt1,u,ssave(2),dsum,ssum
      real*8 bsig,ssig,stmp,btmp
      external e_factorial
      parameter (pi=3.14159265358979)

c only background?

      do i=1,2
        ssave(i) = s(i)
      enddo
      if (bgonly) then
        do i=1,2
          s(i) = 0.0D0
        enddo
      endif

c handle gracefully cases with zero error also.  Do others numerically
c within +-3 sigma

      if (s(2).le.0.0D0 .and. b(2).le.0.0D0) then
        u = s(1) + b(1)
        call e_poiss(s(1)+b(1),n,prob)
      elseif (s(2).le.0.0D0) then
        dbt1 = 1.0D0/(sqrt(2.0D0*pi)*b(2))
        dsum = 0.0D0
        ssum = 0.0D0
        stmp = s(1)
        do bsig=-systwid,systwid,dsystwid
          btmp = b(1) + bsig*b(2)
          if (btmp.gt.0.0D0) then
            dbt = dbt1*exp(-bsig**2/2.0D0)
            dsum = dsum + dbt
            u = stmp + btmp
            call e_poiss(u,n,prob)
            ssum = ssum + prob*dbt
          endif
        enddo
        prob = ssum/dsum
      elseif (b(2).le.0.0D0) then
        dst1 = 1.0D0/(sqrt(2.0D0*pi)*s(2))
        dsum = 0.0D0
        ssum = 0.0D0
        btmp = b(1)
        do ssig=-systwid,systwid,dsystwid
          stmp = s(1) + ssig*s(2)
          if (stmp.gt.0.0D0) then
            dst = dst1*exp(-ssig**2/2.0D0)
            dsum = dsum + dst
            u = stmp + btmp
            call e_poiss(u,n,prob)
            ssum = ssum + prob*dst
          endif
        enddo
        prob = ssum/dsum
      else
        dst1 = 1.0D0/(sqrt(2.0D0*pi)*s(2))
        dbt1 = 1.0D0/(sqrt(2.0D0*pi)*b(2))
        dsum = 0.0D0
        ssum = 0.0D0
        do ssig=-systwid,systwid,dsystwid
          stmp = s(1) + ssig*s(2)
          if (stmp.ge.0.0D0) then
            dst = dst1*exp(-ssig**2/2.0D0)
            do bsig=-systwid,systwid,dsystwid
              btmp = b(1) + bsig*b(2)
              if (btmp.gt.0.0D0) then
                dbt = dbt1*exp(-bsig**2/2.0D0)
                dsum = dsum + dst*dbt
                u = stmp + btmp
                call e_poiss(u,n,prob)
                ssum = ssum + prob*dst*dbt
              endif
            enddo
          endif
        enddo
        prob = ssum/dsum
      endif

      do i=1,2
        s(i) = ssave(i)
      enddo

      return
      end

      real*8 function e_factorial(n)
      implicit none
      integer n,i
      real*8 fl(0:100)
      save fl
      real*8 e_f2
      external e_f2
      logical first /.true./
      save first
      if (first) then
        first = .false.
        do i=0,100
          fl(i) = e_f2(i)
        enddo
      endif
      if (n.gt.100) then
        e_factorial = e_f2(n)
      else
        e_factorial = fl(n)
      endif
      return
      end

      real*8 function e_f2(n)
      implicit none
      integer i,n
      e_f2 = 1.0D0
      do i=1,n
        e_f2 = e_f2 * i
      enddo
      return
      end
