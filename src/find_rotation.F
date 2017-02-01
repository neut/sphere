      program find_rotation
      ! main program
      ! 
      ! governs all subroutines. should open an 1D concatentated image 
      ! from a file, save this image into a 2-dimensional array and use
      ! nelmin() to solve compare_frames() for rot such that f -> 0
      
      use param
      use omp_lib

      implicit none

      integer         :: mode, draw, fcalls=0, n, i, j, num, numres
      real            :: rot(4), rot_ming(4)
      real            :: fmin=1e6
      real(kind=8)    :: t_start, t_end 
      character(128)  :: path, filename, listfile, writefile, f_list
      character(128)  :: f_out, nstr, ifilestr

c-----mode = 1 ! 0: slow. 1: fast

      path     = 'data/'
      f_list   = 'n_050/img_list.txt'
      f_out    = 'n_050/output.txt'

      listfile  = trim(path)//trim(f_list)
      writefile = trim(path)//trim(f_out)

c-----find number of images to load
      open(102,file=listfile,status='unknown')
      rewind(102)

      read(102,*) mode
      read(102,*) draw
      read(102,*) num
      read(102,*) lnx
      read(102,*) lny

c-----allocate largest image (mostly images will be same size)
c-----for M and synimg, the array will not be filled completely
c-----for x and y, x(1) = 0 o and x(nx) = 1 at all times. see
c-----projection.f90

      allocate(M(lny,lnx))
      allocate(synimg(lny,lnx))
      allocate(x(lnx))
      allocate(y(lny))
      allocate(r(lny,lnx))
      
c-----allocate and fill x and y values of the gridpoints on the
c-----projection grid

      do j=1,lny
        y(j) = (2.0*j - lny - 1.0)/(lny - 1.0)
        do i=1,lnx
          x(i) = (2.0*i - lnx - 1.0)/(lnx - 1.0)
          r(i,j) = sqrt(x(i)*x(i) + y(j)*y(j))
        enddo
      enddo

      rot = (/0, 0, 0, 0/)

3     format(4(F20.10))
4     format(5(F10.3), I5)

      open(111,file=writefile,status='unknown')
    
c-----start solving

c-----SLOW MODE-------------------------------------------------
      if (mode == 0) then
        write(*,*) '-----------SLOW MODE-----------'
       do i=1,num
        read(102,*) filename
c-------reset fmin and fcalls
        fmin = 1e6
        fcalls = 0

        call load_img(filename)
        call solver(0, rot, rot_ming, fmin, fcalls, numres)

        call drawpattern(draw)

        write(111,3) rot_ming(1), rot_ming(2), rot_ming(3), 
     &               mod(rot_ming(4)/pi*180,360.0)
      
        write(*,*) '------------------------------------------------'
        write(*,4) rot_ming(1), rot_ming(2), rot_ming(3),
     &             mod(rot_ming(4)/pi*180,360.0), i
        write(*,*) '------------------------------------------------'
       
       enddo
c-----FAST MODE-------------------------------------------------
      elseif (mode == 1) then
        write(*,*) '-----------FAST MODE----------'
        do i = 1,num ! iterate through each frame
          read(102,*) filename  ! load each frames filename
          call load_img(filename) ! and load it into memory
 
          if (i == 1) then ! if at the first frame, do a 'slow' solve
            call solver(0, rot, rot_ming, fmin, fcalls, numres)

            rot = rot_ming 

c-----------write to file   
            write(111,3) rot_ming(1), rot_ming(2), rot_ming(3),
     &                   mod(rot_ming(4)/pi*180,360.0)
c-----------write to screen
            write(*,*) '-----------------------------------------------'
            write(*,*) rot_ming(1), rot_ming(2), rot_ming(3),
     &                 mod(rot_ming(4)/pi*180,360.0), fcalls
            write(*,*) '-----------------------------------------------'

            call drawpattern(draw)
c-----------etime is called after the calculation of the first frame      
            t_start = omp_get_wtime()

          else ! for the rest of the frames, solve causally using the previous solution
            call solver(1, rot, rot_ming, fmin, fcalls, numres)
 
            call drawpattern(draw)
      
            write(111,3) rot_ming(1), rot_ming(2), rot_ming(3),
     &                   mod(rot_ming(4)/pi*180,360.0)

            write(*,4) rot_ming(1), rot_ming(2), rot_ming(3),
     &                 mod(rot_ming(4)/pi*180,360.0), fmin, i

             rot = rot_ming
          endif
        enddo

        t_end = omp_get_wtime()

      endif

      close(102)
      close(111)

      write(*,*) 'Elapsed time fast frames only', t_end - t_start 

      end program find_rotation
