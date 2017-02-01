      real function compare_frames(rot)

      use param, only: synimg, M, nx, ny

      real, intent(in)  :: rot(4) 
      real              :: f

      call projection(rot)
      
      f = 0.0

!$OMP  PARALLEL DO SCHEDULE(GUIDED)
!$OMP$ DEFAULT(SHARED)
!$OMP$ PRIVATE(i,j)
!$OMP$ REDUCTION(+:f)
      do j=1,ny
        do i=1,nx
          f = f + real(abs(M(i,j) - synimg(i,j)))
        enddo
      enddo
!$OMP  END PARALLEL DO

c-----Non-OpenMP version would be:
c-----f    = real(sum(abs(M(i,j) - synimg(i,j)))))

      compare_frames = f
      return

      end function
