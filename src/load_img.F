      subroutine load_img(filename)

      use param
      integer        :: i, j
      character(128) :: filename
    
      filename = trim(filename)

      open(101,file=filename,status='unknown')
      rewind(101)

      read(101,*) ny
      read(101,*) nx

      read(101,*) ((M(j,i), i=1,nx), j=1,ny)
      close(101)
 
      end subroutine
