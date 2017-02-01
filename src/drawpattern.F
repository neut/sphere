      subroutine drawpattern(draw)

      use param
      implicit none
      
      integer :: i, j, draw
      
      if (draw == 1) then
        do i=1,ny
          do j=1,nx
            write(*,'(I1)',advance='no') synimg(ny+1-j,i)
          enddo
          write(*,*) ''
         enddo
      endif
      
      end subroutine drawpattern
