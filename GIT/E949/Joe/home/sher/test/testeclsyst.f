c single-channel tester for eclsyst.f

      program testeclsyst
      implicit none
      real*8 s(2,500),b(2,500),cls,clb,clsb,acls,aclb,aclsb
      integer d(500),n,i
c      integer iwork(500)
c      real*8 clstoy,clbtoy,clsbtoy

 10   continue
      write(*,*) 'Number of channels: '
      read(*,*,end=999) n
      if (n.gt.10) then
        do i=1,n
          s(1,i) = 3.0/n
          s(2,i) = 0.0
          b(1,i) = 0.0
          b(2,i) = 0.0
          d(i) = 0
        enddo
      else
        do i=1,n
          write(*,*) 'Expected signal and error in channel ',i,':'
          read(*,*) s(1,i),s(2,i)
          write(*,*) 'Expected background and error in channel ',i,':'
          read(*,*) b(1,i),b(2,i)
          write(*,*) 'Observed candidates in channel ',i,':'
          read(*,*) d(i)
        enddo
      endif

      call e_cls(n,s,b,d,cls,clsb,clb)
      write(*,*) 'CL(s): ',cls
      write(*,*) 'CL(s+b): ',clsb
      write(*,*) 'CL(b): ',clb

c      call e_clstoymc(n,s,b,d,iwork,clstoy,clsbtoy,clbtoy)
c      write(*,*) 'CL(s) toyMC: ',clstoy
c      write(*,*) 'CL(s+b) toyMC: ',clsbtoy
c      write(*,*) 'CL(b) toyMC: ',clbtoy

      call e_avcls_f(n,s,b,acls,aclsb,aclb)

      write(*,*) '<CL(s)>: ',acls
      write(*,*) '<CL(s+b)>: ',aclsb
      write(*,*) '<CL(b)>: ',aclb

      goto 10
 999  continue
      end
