      program readdata

      implicit none

      integer nx, ns
      parameter (nx = nxnum, ns = nsnum)

      real*4 dataarr(nx, ns)
      real*8 x(nx)
      real*4 r4_x
      real*8 xmin, xmax
      real*8 xb(0:nx), xc(nx), dx
      integer i, j, ll
      real*4 maxdata,mindata

      character*4 arg
      character*8 file1

      call getarg(1, arg)
      file1 = arg//"0001"

      xmin = xmin1
      xmax = xmax1

      dx = ( xmax - xmin) / nx
      xb ( 0) = xmin

      DO i = 1, nx
        xb(i) = xb(i-1)+dx
        xc(i) = 0.5*(xb(i)+xb(i-1))
      END DO

      open (12, file='filelist', status='new')

      do ll = 1, nslice

      open (10, file=file1, status='old')
      open (11, file=file1//'_rewrite', status='new')

      read (10, *) ((dataarr(i, j), i = 1, nx), j = 1, ns)

      maxdata = dataarr (1, 1)
      mindata = dataarr (1, 1)
      do j = 1, ns
      do i = 1, nx
        r4_x = xc(i)
        write(11,*) r4_x,dataarr(i, j)
        if (dataarr(i, j).gt.maxdata) maxdata = dataarr(i, j)
        if (dataarr(i, j).lt.mindata) mindata = dataarr(i, j)
      enddo
        write(11,*) ' '
      enddo

      close(10)
      close(11)
      write(12,*) file1

      CALL INCNAME_8 ( file1, 1)
      enddo

      if (abs(mindata-maxdata).le.1.e-16) then
        if ( abs(mindata).le.1.e-16) then
          mindata = -1.e-16
          maxdata = 1.e-16
        else
          mindata = mindata - abs ( mindata) * 1.e-16
          maxdata = maxdata + abs ( maxdata) * 1.e-16
        endif
      endif

      open (10, file='min', status='new')
      open (11, file='max', status='new')
      write(10,*) mindata
      write(11,*) maxdata
      close(10)
      close(11)
      close(12)


      end program readdata


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INCNAME_8 ( filname, inc)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL MYSTOP

      CHARACTER*8 filname

      CHARACTER*1 unitchar ( 8)

      INTEGER unitint ( 8), zeroint, value, inc, i, j

! ----------------------------------------------------------------------

      unitchar ( 1) = filname ( 1:1)
      unitchar ( 2) = filname ( 2:2)
      unitchar ( 3) = filname ( 3:3)
      unitchar ( 4) = filname ( 4:4)
      unitchar ( 5) = filname ( 5:5)
      unitchar ( 6) = filname ( 6:6)
      unitchar ( 7) = filname ( 7:7)
      unitchar ( 8) = filname ( 8:8)

      zeroint = ichar ( '0')

      unitint ( 8) = ichar ( unitchar ( 8) ) - zeroint
      unitint ( 7) = ichar ( unitchar ( 7) ) - zeroint
      unitint ( 6) = ichar ( unitchar ( 6) ) - zeroint
      unitint ( 5) = ichar ( unitchar ( 5) ) - zeroint

      value = 1000 * unitint ( 5) + 100 * unitint ( 6)
     1              + 10 * unitint ( 7) + unitint ( 8)

      value = value + inc


      unitint ( 5) = value / 1000

      value = value - 1000 * unitint ( 5)

      unitint ( 6) = value / 100

      value = value - 100 * unitint ( 6)

      unitint ( 7) = value / 10

      value = value - 10 * unitint ( 7)

      unitint ( 8) = value

      unitint ( 5) = unitint ( 5) + zeroint
      unitint ( 6) = unitint ( 6) + zeroint
      unitint ( 7) = unitint ( 7) + zeroint
      unitint ( 8) = unitint ( 8) + zeroint

      unitchar ( 5) = char ( unitint ( 5) )
      unitchar ( 6) = char ( unitint ( 6) )
      unitchar ( 7) = char ( unitint ( 7) )
      unitchar ( 8) = char ( unitint ( 8) )

      filname ( 5:5) = unitchar ( 5)
      filname ( 6:6) = unitchar ( 6)
      filname ( 7:7) = unitchar ( 7)
      filname ( 8:8) = unitchar ( 8)

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INCNAME_8

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INCNAME_SP ( filname, inc)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL MYSTOP

      CHARACTER*8 filname

      CHARACTER*1 unitchar ( 8)

      INTEGER unitint ( 8), zeroint, value, inc, i, j

! ----------------------------------------------------------------------

      unitchar ( 1) = filname ( 1:1)
      unitchar ( 2) = filname ( 2:2)
      unitchar ( 3) = filname ( 3:3)
      unitchar ( 4) = filname ( 4:4)
      unitchar ( 5) = filname ( 5:5)
      unitchar ( 6) = filname ( 6:6)
      unitchar ( 7) = filname ( 7:7)
      unitchar ( 8) = filname ( 8:8)

      zeroint = ichar ( '0')

      unitint ( 4) = ichar ( unitchar ( 4) ) - zeroint
      unitint ( 3) = ichar ( unitchar ( 3) ) - zeroint

      value = 10 * unitint ( 3) + unitint ( 4)

      value = value + inc

      unitint ( 3) = value / 10

      value = value - 10 * unitint ( 3)

      unitint ( 4) = value

      unitint ( 3) = unitint ( 3) + zeroint
      unitint ( 4) = unitint ( 4) + zeroint

      unitchar ( 3) = char ( unitint ( 3) )
      unitchar ( 4) = char ( unitint ( 4) )

      filname ( 3:3) = unitchar ( 3)
      filname ( 4:4) = unitchar ( 4)

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INCNAME_SP

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
