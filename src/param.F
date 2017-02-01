      module param
      
      implicit none
      integer,allocatable :: M(:,:)
      integer,allocatable :: synimg(:,:)
      integer             :: nx,ny, firstframeloaded = 0
      integer             :: lnx, lny
      real, allocatable   :: x(:), y(:), r(:,:)

      real, parameter     :: pi    = 3.14159265358979323
      real, parameter     :: pi2   = 0.5*3.14159265358979323
!     real, parameter     :: num1  = 20.0/180.0*pi
      real, parameter     :: num2  = 22.5/180.0*pi
      real, parameter     :: num3  = 45.0/180.0*pi
      real, parameter     :: num4  = 66.5/180.0*pi
!     real, parameter     :: num5  = 70.0/180.0*pi
!     real, parameter     :: num6  = 90.0/180.0*pi
!     real, parameter     :: num7  = 110.0/180.0*pi
      real, parameter     :: num8  = 113.5/180.0*pi
      real, parameter     :: num9  = 135.0/180.0*pi
!     real, parameter     :: num10 = 180.0/180.0*pi
!     real, parameter     :: num11 = 200.0/180.0*pi
!     real, parameter     :: num12 = 340.0/180.0*pi

      real, parameter     :: num1  = 30.0/180*pi
      real, parameter     :: num5  = 60.0/180.0*pi
      real, parameter     :: num6  = 90.0/180.0*pi
      real, parameter     :: num7  = 120.0/180.0*pi
      real, parameter     :: num10 = 180.0/180.0*pi
      real, parameter     :: num11 = 210.0/180.0*pi
      real, parameter     :: num12 = 330.0/180.0*pi

      real, parameter     :: x1    = sin(num4)*cos(num4)
      real, parameter     :: x2    = sin(num4)*cos(num8)
      real, parameter     :: y1    = sin(num2)*sin(num6)

      end module param
