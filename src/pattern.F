c-----FINAL VERSION----------------------------------------------------
      subroutine pattern(theta, phi, s)
c-----pattern() outputs the color of the sphere surface as function of
c-----the spherical coordinates given by theta and phi

      use param, only: pi, pi2, x1, x2, y1, num1, num2, num3, num4,
     &                 num5, num6, num7, num8, num9, num10, num11, 
     &                 num12
      
      implicit none

      integer          :: s
      real             :: theta, phi, y, z
      real, parameter  :: bw = 0.353553390593274 ! sqrt(2)/4

      theta = mod(mod(theta,2.0*pi)+2.0*pi,2.0*pi)

c-----Initialize output to zero
      s = 0

      if (theta >= pi .and. phi >= pi2) then
        s = 1
      elseif (theta <= pi .and. phi <= pi2) then
        s = 1
      endif
      

      if ( (phi >= num6) .and. (phi <= num7) .and. 
     &     (theta >= num12) ) then
        s = 0
      endif

      if ( (phi >= num5 .and. phi <= num6) .and.    
     &     (theta >= num10 .and. theta <= num11) ) then
        s = 1
      endif

      y = sin(phi)*sin(theta)
      z = cos(phi)

      if ( (abs(sin(phi)*cos(theta)) <= (bw + y)) .and. (theta > pi) )
     &then
        s = 1
      endif

      if ( (abs(sin(phi)*cos(theta)) <= (bw - z)) .and. (theta < pi) ) 
     &then
        s = 0
      endif

      end 


