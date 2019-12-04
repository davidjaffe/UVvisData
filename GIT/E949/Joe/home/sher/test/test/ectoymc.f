c toy MC to compute the same things ecl.f computes.  Test of
c validity:

      subroutine e_clstoymc(n,s,b,d,iwork,cls,clsb,clb)
      implicit none
      integer n,d(*),iwork(*)
      real*8 s(*),b(*),cls,clsb,clb
      real*8 dest,dtest
      integer i,j,ierr
      integer nexpt,nsb,nb
      parameter (nexpt=100000)
      real rsingle

      call e_dest(n,s,b,d,dest)

      nsb = 0
      nb = 0
      do i=1,nexpt

c s+b
        do j=1,n
          rsingle = s(j) + b(j)
          if (rsingle.le.0.0) then
            iwork(j) = 0
          else
            call rnpssn(rsingle,iwork(j),ierr)
          endif
          if (ierr.ne.0) then
            write(*,*) 'RNPSSN problem ',ierr
          endif
        enddo
        call e_dest(n,s,b,iwork,dtest)
        if (dtest.le.dest) then
          nsb = nsb + 1
        endif

c b
        do j=1,n
          rsingle = b(j)
          if (rsingle.le.0.0) then
            iwork(j) = 0
          else
            call rnpssn(rsingle,iwork(j),ierr)
          endif
          if (ierr.ne.0) then
            write(*,*) 'RNPSSN problem ',ierr
          endif
        enddo
        call e_dest(n,s,b,iwork,dtest)
        if (dtest.le.dest) then
          nb = nb + 1
        endif

      enddo
        
      clb = dble(nb)/dble(nexpt)
      clsb = dble(nsb)/dble(nexpt)
      if (clb.gt.0.0D0) then
        cls = clsb/clb
      else
        cls = 1.0D0
      endif

      return
      end

