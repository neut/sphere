      subroutine rotated_pattern(theta, phi, rot, s)
c-----! rotated_pattern() outputs the synthetic color of some rotated coordinates
c-----! Uses Rodriges' formula to rotate along axis rot(1:3) with angle rot(4)

      real,intent(in)     :: theta, phi, rot(4)
      real                :: k(3), kxv(3), kdv,  alpha, knorm 
      real                :: v(3), vnew(3)
      real                :: theta_new, phi_new
      integer,intent(out) :: s

      v(1) = sin(phi)*cos(theta)
      v(2) = sin(phi)*sin(theta)
      v(3) = cos(phi)
           
      knorm = sqrt(rot(1)**2 + rot(2)**2 + rot(3)**2)
      k(1)  = rot(1)/knorm
      k(2)  = rot(2)/knorm
      k(3)  = rot(3)/knorm
      alpha = rot(4)

!!    Cross product
      kxv(1) = k(2)*v(3) - k(3)*v(2)
      kxv(2) = k(3)*v(1) - k(1)*v(3)
      kxv(3) = k(1)*v(2) - k(2)*v(1)
        
!!    Dot product
      kdv = k(1)*v(1) + k(2)*v(2) + k(3)*v(3)

!     Rodrigues formula 
      vnew = v*cos(alpha) + kxv*sin(alpha) + k*kdv*(1 - cos(alpha))

      theta_new = atan2(vnew(2), vnew(1))
      phi_new   = acos(vnew(3))

      call pattern(theta_new,phi_new,s) 
      
      end subroutine     




