      program ecl9598inside
*
*    File:  ecl9598inside.f
* Purpose:  single-channel tester for eclsyst.f -- find CL
*           for different values of the K -> pnn BR
*   Usage:  > make -f Makefile9598inside
*           > ecl9598inside.exe >& ecl9598inside.dat
*           > awk -f ecl.awk ecl9598inside.dat >& ecl9598inside.vect
*           > eclresults
*
* History:  2001 Sep  PB  created K->pnn version of Tom Junk's testecl.f
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
      real*8 xobs_cls,bgscale
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
*      open(71,file='bgtotal.inside.dat',status='old')
*      read(71,*,end=999) n
*      do i=1,n
*        read(71,*) q(i),q(i),b(i),s(i),q(i),
*     &             q(i),q(i),q(i),q(i),q(i),q(i),q(i)
*      enddo


c
c  Get the 1995-7 expected signal and background values
c  (2 bins inside the 1995-7 box only).
c  s(i) summed over inside-the-box bins = 1*(1998 s.e.s./1995-7 s.e.s.).
c  b(i) summed over inside-the-box bins = 1.
c
*      n = n + 1
*      ns95 = n
*      b(n) = 0.1
*      s(n) = 0.55*sens98/sens9597
*      n = n + 1
*      b(n) = 0.9
*      s(n) = 0.45*sens98/sens9597
* 1. Put all of 98 stuff in one bin
      n = n + 1
      b(n) = 1.
      s(n) = 1.
* 2. Put all of 95-97 results in 1 bin
      n = n + 1
      ns95 = n
      b(n) = 1.
      s(n) = sens98/sens9597

      open(71,file='bgscale.dat',status='old')
      read(71,*) bgscale
      close (71)
      write(6,*) bgscale
      open(72,file='ecl9598inside.dat',status='new')
c
c  Convert the b(i) into the absolute expected number of background events.
c
      do i=1,ns95-1
        b(i) = b(i)*boxbg98*bgscale
        write(6,*) ' BG for bin ',i,b(i)
      enddo
      do i=ns95,n
        b(i) = b(i)*boxbg9597*bgscale
        write(6,*) ' BG for bin ',i,b(i)
      enddo

c
c  Observed number of events in each bin:
c
*      d(ns95) = 1   ! 1995 signal event
*      d(64) = 1     ! 1998 signal event
      d(ns95) = 1   ! 1995 signal event
      d(1) = 1     ! 1998 signal event

c
c  Loop over scaling factors (i.e., K -> pnn branching ratios)
c  which convert the s(i) into the absolute expected number of signal events,
c  and calculate the confidence level for each scaling factor.
c
      do scale = 0.0,5.0,0.01
        do i=1,n
          ss(i) = s(i)*scale
        enddo

        call e_cls(n,ss,b,d,cls,clsb,clb)
        xobs_cls = xobs
        call e_avcls(n,ss,b,acls,aclsb,aclb)
        call e_clb_lr(n,ss,b,d,clbm1)
        call e_avclb_lr(n,ss,b,asigclb)
c        call e_clstoymc(n,ss,b,d,iwork,clstoy,clsbtoy,clbtoy)

        write(72,*) 'signal scaling: ',scale
        write(72,*) 'Xobs: ',xobs_cls
        write(72,*) 'CL(s): ',cls
        write(72,*) 'CL(s+b): ',clsb
        write(72,*) 'CL1(b): ',clb
        write(72,*) 'CL2(b): ',1-clbm1

c        write(72,*) 'CL(s) toyMC: ',clstoy
c        write(72,*) 'CL(s+b) toyMC: ',clsbtoy
c        write(72,*) 'CL(b) toyMC: ',clbtoy

        write(72,*) '<CL(s)>: ',acls
        write(72,*) '<CL(s+b)>: ',aclsb
        write(72,*) '<CL1(b)>: ',aclb
        write(72,*) '<CL2(b)>: ',asigclb(3)
      enddo

 999  continue

      close (71)
      
      end
