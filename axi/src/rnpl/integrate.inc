cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     integrate: manually update (integrate) the total charge and
c                energy. Note that RNPL supplies the routine header,
c                declaration of some routine arguments, and RETURN/END
c                statements
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer i, j
      real*8  Qtotsum, Etotsum, Eemsum
      integer unum, stat

      real*8      pi
      parameter ( pi = 3.14159265358979d0 )

      real*8, allocatable :: jac(:,:)
      allocate(jac(g1_NR,g1_NX))

      do i=1, g1_NR
      do j=1, g1_NX

        ! linear coordinates
c        jac(i,j) = R(i)

        ! compactified coordinates
        jac(i,j) =
     #  c ** 2 * d ** 3 * (exp(c * X(j)) + exp(-c * X(j))) * (exp(c 
     #* R(i)) + exp(-c * R(i))) * (exp(c * R(i)) - exp(-c * R(i)))

      end do
      end do

      ! total charge
      unum = 20
      open(unit=unum, file='Qtot.dat', status='unknown',
     &     action='write', position='append', iostat=stat)

      Qtotsum = 0.0d0

      do i=1, g1_NR-1
      do j=1, g1_NX-1
        Qtotsum = Qtotsum + dR*dX/4.0d0*
     &  (   jac(i,j)*Qden_n(i,j)     + jac(i+1,j)*Qden_n(i+1,j)
     &    + jac(i,j+1)*Qden_n(i,j+1) + jac(i+1,j+1)*Qden_n(i+1,j+1) )
      end do
      end do

      Qtotsum = 2d0*pi*Qtotsum

      write(unum,*,iostat=stat) t, Qtotsum
      close(unit=unum,   status='keep', iostat=stat)

      ! total energy
      unum = 20
      open(unit=unum, file='Etot.dat', status='unknown',
     &     action='write', position='append', iostat=stat)

      Etotsum = 0.0d0

      do i=1, g1_NR-1
      do j=1, g1_NX-1
        Etotsum = Etotsum + dR*dX/4.0d0*
     &  (   jac(i,j)*Eden_n(i,j)     + jac(i+1,j)*Eden_n(i+1,j)
     &    + jac(i,j+1)*Eden_n(i,j+1) + jac(i+1,j+1)*Eden_n(i+1,j+1) )
      end do
      end do

      Etotsum = 2d0*pi*Etotsum

      write(unum,*,iostat=stat) t, Etotsum
      close(unit=unum,   status='keep', iostat=stat)

      ! EM field energy
      unum = 20
      open(unit=unum, file='Eem.dat', status='unknown',
     &     action='write', position='append', iostat=stat)

      Eemsum = 0.0d0

      do i=1, g1_NR-1
      do j=1, g1_NX-1
        Eemsum = Eemsum + dR*dX/4.0d0*
     &  (   jac(i,j)*Eem_n(i,j)     + jac(i+1,j)*Eem_n(i+1,j)
     &    + jac(i,j+1)*Eem_n(i,j+1) + jac(i+1,j+1)*Eem_n(i+1,j+1) )
      end do
      end do

      Eemsum = 2d0*pi*Eemsum

      write(unum,*,iostat=stat) t, Eemsum
      close(unit=unum,   status='keep', iostat=stat)

      dummy_np1(:,:) = 0d0
