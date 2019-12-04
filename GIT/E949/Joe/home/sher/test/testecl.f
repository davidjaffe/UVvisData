c single-channel tester for ecl.f

      program testecl
      implicit none
      integer nbins
      parameter (nbins = 7600)
      real*8 s(nbins),b(nbins),cls,clb,clsb,acls,aclb,aclsb
      real*8 clbm1,asigclb(5)
      real*8 x(nbins),scale,ss(nbins),bg
      integer d(nbins),n,ns95,i,j
      integer iwork(nbins)
      real*8 clstoy,clbtoy,clsbtoy
c  vary the bin in which an inside-the-box event is put
      integer test(9)
      data test /1,7,20,37,132,235,372,406,486/

c
c  Get the inputs
c
      n = 0

c  1998 expected signal and background
      open(71,file='bgtotal.dat',status='old')
c      open(71,file='bgtotal.inside.dat',status='old')
      read(71,*,end=999) n
      write(*,*) 'Number of channels: ',n
      do i=1,n
        read(71,*) x(i),x(i),b(i),s(i),x(i),
     &             x(i),x(i),x(i),x(i),x(i),x(i),x(i)
      enddo

c  1995 expected signal and background
c      n = n + 1
c      ns95 = n
c      b(n) = 0.55
c      s(n) = 0.1
c      n = n + 1
c      b(n) = 0.45
c      s(n) = 0.9

c
c  Observed number of events in each bin:
c
c      do j = 1,9
c      do bg = 0.02,0.2,0.02

      do i=1,n
        d(i) = 0
      enddo

c  for inside-the-box function:
c      d(ns95) = 1   ! 1995 signal event
c      d(37) = 1     ! 1998 signal event
c  for total function:
c      d(ns95) = 1   ! 1995 signal event
      d(47) = 1     ! 1998 signal event
      d(461) = 1    ! signal-like event (bgtd = 3.1)
c      d(5378) = 1   ! Km2 L/R ambiguity event

c      d(test(j)) = 1
c      write (*,*) 'event in bin ',test(j)
c      write (*,*) 'background in box = ',bg

c
c  N and A are scaled such that N = 1, A = 1 corresponds to the final box.
c  Scale s(i) such that it corresponds to the expected number of events
c  in the box (i.e., the K -> pnn branching ratio), and scale b(i) such
c  that it corresponds to the background level in the box (see bgtotal.c).
c
c      do i=1,ns95-1
      do i=1,n
        b(i) = b(i) * 0.0828             ! 1998 background level
      enddo
c      do i=ns95,n
c        b(i) = b(i) * 0.0675             ! 1995-7 background level
c      enddo
c      do i = 1,n
c        b(i) = b(i) * 0.0828             ! 1998 background level
c        b(i) = b(i) * 0.0675             ! 1995-7 background level
c        b(i) = b(i) * bg
c        b(i) = b(i) * (0.0828 + 0.0675)  ! 1995-7 + 1998 background level
c      enddo
      do scale = 0.05,5.50,0.05
        do i=1,n
          ss(i) = s(i) * scale
        enddo


c
c  Calculate confidence interval
c
        call e_cls(n,ss,b,d,cls,clsb,clb)
        call e_avcls(n,ss,b,acls,aclsb,aclb)
        call e_clb_lr(n,ss,b,d,clbm1)
        call e_avclb_lr(n,ss,b,asigclb)
c        call e_clstoymc(n,ss,b,d,iwork,clstoy,clsbtoy,clbtoy)

        write(*,*) 'signal scaling: ',scale
        write(*,*) 'CL(s): ',cls
        write(*,*) 'CL(s+b): ',clsb
c        write(*,*) 'CL(b): ',clb
        write(*,*) 'CL(b): ',1-clbm1

c        write(*,*) 'CL(s) toyMC: ',clstoy
c        write(*,*) 'CL(s+b) toyMC: ',clsbtoy
c        write(*,*) 'CL(b) toyMC: ',clbtoy

        write(*,*) '<CL(s)>: ',acls
        write(*,*) '<CL(s+b)>: ',aclsb
c        write(*,*) '<CL(b)>: ',aclb
        write(*,*) '<CL(b)>: ',asigclb(3)
      enddo
c      enddo

 999  continue
      end
