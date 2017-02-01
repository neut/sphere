      subroutine projection(rot)

      use param, only: x, y, r, nx, ny, synimg

      integer :: i, j, s
      real    :: theta, phi
      real    :: rot(4) 

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP$ DEFAULT(SHARED)
!$OMP$ PRIVATE(i,j,phi,theta,s)
      do j=1,ny      
        do i=1,nx       
          synimg(j,i) = 1

          if (r(i,j) < 1.0d0) then 
            phi         = acos(y(j))
            theta       = asin(x(i)/sin(phi))

            call rotated_pattern(theta, phi, rot, s)

            synimg(j,i) = s
          endif
        enddo
      enddo
!$OMP  END PARALLEL DO
      
      end subroutine
