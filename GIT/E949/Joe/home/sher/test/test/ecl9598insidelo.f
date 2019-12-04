      program ecl9598insidelo
*
*    File:  ecl9598insidelo.f
* Purpose:  single-channel tester for eclsyst.f -- find CL
*           for different values of the K -> pnn BR
*   Usage:  > make -f Makefile9598insidelo
*           > ecl9598insidelo.exe >& ecl9598insidelo.dat
*           > awk -f ecl.awk ecl9598insidelo.dat >& ecl9598insidelo.vect
*           > eclresults
*
      implicit none
      include 'eclout.inc'

      integer nbins
      parameter (nbins = 8000)
      real*8 s(nbins),b(nbins),cls,clb,clsb,acls,aclb,aclsb
      real*8 clbm1,asigclb(5)
      real*8 q(nbins),scale,ss(nbins)
      integer d(nbins),n,ns95,i
      integer iwork(nbins)
      real*8 clstoy,clbtoy,clsbtoy
      real*8 xobs_cls
c  Total background levels and their uncertainties
c  inside the 1995-7 and 1998 boxes
      real*8 boxbg9597,dboxbg9597,boxbg98,dboxbg98
      parameter (boxbg9597 = 0.0804, dboxbg9597 = 0.0201)
      parameter (boxbg98 = 0.0657, dboxbg98 = 0.0438)
c  1/single-event-sensitivities
c  error on signal ~= 5%  (from UMC nuclear interaction code -- see E787 TN365)
      real*8 sens9597,sens98,relerrsens
      parameter (sens9597 = 1.496, sens98 = 1.892, relerrsens = 0.05)

c
c  Initialize.
c
      do i = 1,nbins
         s(i) = 0
         b(i) = 0
         ss(i) = 0
         s(i) = 0
         b(i) = 0
         ss(i) = 0
         q(i) = 0
         d(i) = 0
      enddo

      n = 0

c
c  Get the 1998 expected signal and background values
c  (inside the box only).
c  s(i) summed over inside-the-box bins = 1 (see bgtotal.c).
c  b(i) summed over inside-the-box bins = 1 (see bgtotal.c).
c
      open(71,file='bgtotal.inside.dat',status='old')
      read(71,*,end=999) n
      do i=1,n
        read(71,*) q(i),q(i),b(i),s(i),q(i),
     &             q(i),q(i),q(i),q(i),q(i),q(i),q(i)
      enddo

c
c  Get the 1995-7 expected signal and background values
c  (2 bins inside the 1995-7 box only).
c  s(i) summed over inside-the-box bins = 1*(1998 s.e.s./1995-7 s.e.s.).
c  b(i) summed over inside-the-box bins = 1.
c
      n = n + 1
      ns95 = n
      b(n) = 0.1
      s(n) = 0.55*sens98/sens9597
      n = n + 1
      b(n) = 0.9
      s(n) = 0.45*sens98/sens9597

c
c  Convert the b(i) into the absolute expected number of background events.
c
      do i=1,ns95-1
        b(i) = b(i)*(boxbg98+dboxbg98)
      enddo
      do i=ns95,n
        b(i) = b(i)*(boxbg9597+dboxbg9597)
      enddo

c
c  Observed number of events in each bin:
c
      d(ns95) = 1   ! 1995 signal event
      d(64) = 1     ! 1998 signal event

c
c  Loop over scaling factors (i.e., K -> pnn branching ratios)
c  which convert the s(i) into the absolute expected number of signal events,
c  and calculate the confidence level for each scaling factor.
c
      do scale = 0.0,5.0,0.01
c      do scale = 10,30,1  ! large expected signal for CLb calculation
        do i=1,n
          ss(i) = s(i)*(1+relerrsens)*scale
        enddo

        call e_cls(n,ss,b,d,cls,clsb,clb)
        xobs_cls = xobs
        call e_avcls(n,ss,b,acls,aclsb,aclb)
        call e_clb_lr(n,ss,b,d,clbm1)
        call e_avclb_lr(n,ss,b,asigclb)
c        call e_clstoymc(n,ss,b,d,iwork,clstoy,clsbtoy,clbtoy)

        write(*,*) 'signal scaling: ',scale
        write(*,*) 'Xobs: ',xobs_cls
        write(*,*) 'CL(s): ',cls
        write(*,*) 'CL(s+b): ',clsb
        write(*,*) 'CL1(b): ',clb
        write(*,*) 'CL2(b): ',1-clbm1

c        write(*,*) 'CL(s) toyMC: ',clstoy
c        write(*,*) 'CL(s+b) toyMC: ',clsbtoy
c        write(*,*) 'CL(b) toyMC: ',clbtoy

        write(*,*) '<CL(s)>: ',acls
        write(*,*) '<CL(s+b)>: ',aclsb
        write(*,*) '<CL1(b)>: ',aclb
        write(*,*) '<CL2(b)>: ',asigclb(3)
      enddo

 999  continue
      end
