      subroutine solver(mode, rot_guess, rot_ming, fmin, fcalls, numres)

      use param
      implicit none

      integer            :: numres, ifault, icount
      integer            :: i, j, mode, fcalls
      real               :: rnorm, rot(4), rot_guess(4)
      real               :: rot_min(4), rot_ming(4)
      real               :: step(4), f, fmin
      real               :: u_axes(26,3), u_angles(8)

      real, parameter    :: reqmin       = 1
      real, parameter    :: minthreshold = 1
      integer, parameter :: konvge       = 1
      integer, parameter :: kcount       = 500

      real, external     :: compare_frames

c-----6 main axes, plus and minus
      u_axes(1,:)  = (/1, 0, 0/)
      u_axes(2,:)  = (/0, 1, 0/)
      u_axes(3,:)  = (/0, 0, 1/)

      u_axes(4,:)  = -(/1, 0, 0/)
      u_axes(5,:)  = -(/0, 1, 0/)
      u_axes(6,:)  = -(/0, 0, 1/)

c-----12 in between each of the main axes (4 for each xy, xz, yz surface)
      u_axes(7,:)  = (/1, 1, 0/)
      u_axes(8,:)  = (/-1, 1, 0/)
      u_axes(9,:)  = (/-1, -1, 0/)
      u_axes(10,:)  = (/1, -1, 0/)

      u_axes(11,:)  = (/1, 0, 1/)
      u_axes(12,:)  = (/-1, 0, 1/)
      u_axes(13,:)  = (/-1, 0, -1/)
      u_axes(14,:)  = (/1, 0, -1/)

      u_axes(15,:)  = (/0, 1, 1/)
      u_axes(16,:)  = (/0, -1, 1/)
      u_axes(17,:)  = (/0, -1, -1/)
      u_axes(18,:)  = (/0, 1, -1/)

c-----8 in the center of each octant
      u_axes(19,:) = (/0.5, 0.5, 0.7071/)
      u_axes(20,:) = (/0.5, 0.5, -0.7071/)
      u_axes(21,:) = (/0.5, -0.5, 0.7071/)
      u_axes(22,:) = (/0.5, -0.5, -0.7071/)

      u_axes(23,:) = -(/0.5, 0.5, 0.7071/)
      u_axes(24,:) = -(/0.5, 0.5, -0.7071/)
      u_axes(25,:) = -(/0.5, -0.5, 0.7071/)
      u_axes(26,:) = -(/0.5, -0.5, -0.7071/)

      u_angles(:) = (/360, 405, 450, 495, 540, 585, 630, 675/) 
c-----u_angles(:) = (/360, 450, 540, 630/) 

c-----Initial simplex determination 
      if (mode == 0) then ! slow
c-----DONT CHANGE THIS AREA, THIS WORKS--------------------------------
       step(1) = 1 
       step(2) = 1
       step(3) = 1
       step(4) = pi     
c-----END DONT CHANGE THIS AREA, THIS WORKS----------------------------

c-----  step(1) = 0.5 
c-----  step(2) = 0.5
c-----  step(3) = 0.5     
c-----  step(4) = 60*pi/180     
      elseif (mode == 1) then ! fast
c-----DONT CHANGE THIS AREA, THIS WORKS--------------------------------
       step(1) = 0.1 
       step(2) = 0.1
       step(3) = 0.1    
       step(4) = 15*pi/180      
c-----END DONT CHANGE THIS AREA, THIS WORKS----------------------------

c----- step(1) = 1
c----- step(2) = 1
c----- step(3) = 1
c----- step(4) = pi
      endif

4     format(5(F10.3))

      if (mode == 1) then
c-----EDUCATED GUESS MODE-----------------------------------------------
       rot = 1.00*rot_guess
       call nelmin(compare_frames, 4, rot, rot_min, f, reqmin, step,
     &             konvge, kcount, icount, numres, ifault)
c----- write(*,*) rot_min/pi*180-360, f
       fmin = f
       rot_ming = rot_min
       fcalls = icount
      elseif (mode == 0) then
c-----AUTOMATIC MODE----------------------------------------------------
       do i = 1,26! loop over axis
         rot(1) = u_axes(i,1)
         rot(2) = u_axes(i,2) 
         rot(3) = u_axes(i,3)

         do j = 1,8 ! loop over angle
           rot(4) = u_angles(j)/180*pi 
 
           call nelmin(compare_frames,4,rot,rot_min,f,reqmin,step,
     &                 konvge, kcount, icount, numres, ifault)
  
           fcalls = fcalls + icount

           if (f < fmin) then
             write(*,4) rot_min(1), rot_min(2), rot_min(3),
     &                  mod(rot_min(4)/pi*180,360.0), f
             fmin = f
             rot_ming = rot_min
            
c-----For testing puposes, the automatic mode breaks for low threshold
             if (fmin < minthreshold) then
               goto 11 
             endif
c-----End of testing purposes

           endif
         enddo
       enddo

11    endif

      rnorm = sqrt(rot_ming(1)**2 + rot_ming(2)**2 + rot_ming(3)**2)

      rot_ming(1) = rot_ming(1)/rnorm
      rot_ming(2) = rot_ming(2)/rnorm
      rot_ming(3) = rot_ming(3)/rnorm

      call projection(rot_ming)

      end subroutine
