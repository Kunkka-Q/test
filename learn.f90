module Grid_Area

   IMPLICIT NONE
   PRIVATE

   !public
   PUBLIC :: set_geometry
   public :: print_geometry
   REAL(8), ALLOCATABLE, PUBLIC :: show(:, :)
   PUBLIC :: Area_Calculation

   !  private
   REAL(8), private :: inR
   REAL(8), private :: outR
   REAL(8), private :: BladeLength
   REAL(8), private :: TankSize
   REAL(8), private :: w
   REAL(8), private :: area_dt              
   REAL(8), private :: rotate_dt           
   REAL(8), private :: tmax                 
   REAL(8), private :: centrePoint(2)
   REAL(8), private :: grid_dx
   REAL(8), private :: grid_dy
   REAL(8), private :: previous_startangle1   
   REAL(8), private :: previous_startangle2   
   REAL(8), private :: angle1
   REAL(8), private :: angle2
   !grid
   PRIVATE :: generate_mesh
   REAL(8), ALLOCATABLE, PRIVATE :: Xgrid(:, :), Ygrid(:, :)
contains
   SUBROUTINE generate_mesh()
      IMPLICIT NONE
      INTEGER :: nx, ny, i, j
      REAL(8), ALLOCATABLE :: XX(:), YY(:)

      nx = INT(TankSize/grid_dx + 1.0D0)
      ny = INT(TankSize/grid_dy + 1.0D0)

      ALLOCATE (XX(nx), YY(ny))
      ALLOCATE (Xgrid(ny, nx), Ygrid(ny, nx))

      DO i = 1, nx
         XX(i) = (i - 1)*grid_dx
      END DO

      DO j = 1, ny
         YY(j) = (j - 1)*grid_dy
      END DO

      DO j = 1, ny
         DO i = 1, nx
            Xgrid(j, i) = XX(i)
            Ygrid(j, i) = YY(j)
         END DO
      END DO
      ALLOCATE (show(ny - 1, nx - 1))
      show = 0.0D0

   END SUBROUTINE generate_mesh

   SUBROUTINE set_geometry(rinR, routR, rBlade, rTank, rw, rAreaDT, rRotateDT, rtmax, rCentre, rdx, rdy, rAngle1, rAngle2)
      REAL(8), INTENT(IN) :: rinR, routR, rBlade, rTank, rw, rAreaDT, rRotateDT, rtmax
      REAL(8), INTENT(IN) :: rCentre(2), rdx, rdy, rAngle1, rAngle2

      inR = rinR
      outR = routR
      BladeLength = rBlade
      TankSize = rTank
      w = rw
      area_dt = rAreaDT
      rotate_dt = rRotateDT
      tmax = rtmax
      centrePoint = rCentre
      grid_dx = rdx
      grid_dy = rdy
      previous_startangle1 = rAngle1
      previous_startangle2 = rAngle2
      angle1 = rAngle1 + rAreaDT*rw
      angle2 = rAngle2 + rAreaDT*rw
      CALL generate_mesh()
   END SUBROUTINE set_geometry

   SUBROUTINE print_geometry()
      PRINT *, "==== Current Geometry Parameters ===="
      PRINT *, "Inner Radius (inR): ", inR
      PRINT *, "Outer Radius (outR): ", outR
      PRINT *, "Blade Length: ", BladeLength
      PRINT *, "Tank Size: ", TankSize
      PRINT *, "Angular Velocity (w): ", w
      PRINT *, "Area Time Step (area_dt): ", area_dt
      PRINT *, "Rotation Time Step (rotate_dt): ", rotate_dt
      PRINT *, "Total Rotation Time (tmax): ", tmax
      PRINT *, "Center Position: ", centrePoint
      PRINT *, "Grid dx: ", grid_dx
      PRINT *, "Grid dy: ", grid_dy
      PRINT *, "Blade1 Initial Angle: ", previous_startangle1
      PRINT *, "Blade2 Initial Angle: ", previous_startangle2
      PRINT *, "==== Current Mesh Grid ====", Xgrid
   END SUBROUTINE print_geometry

   function trapezoid_area(points) result(area)
      implicit none
      real(8), intent(in) :: points(4, 2)
      real(8) :: area
      real(8) :: sorted_points(4, 2)
      real(8) :: top_points(2, 2), bottom_points(2, 2)
      real(8) :: top_base, bottom_base, height
      integer :: i, j, max_idx, k
      real(8) :: y_values(4)

      do i = 1, 4
         do j = 1, 2
            sorted_points(i, j) = points(i, j)
         end do
         y_values(i) = points(i, 2)
      end do

      do i = 1, 3
         max_idx = i
         do j = i + 1, 4
            if (sorted_points(j, 2) > sorted_points(max_idx, 2)) then
               max_idx = j
            end if
         end do
         if (max_idx /= i) then
            do k = 1, 2
               call swap(sorted_points(i, k), sorted_points(max_idx, k))
            end do
         end if
      end do

      top_points = sorted_points(1:2, :)
      bottom_points = sorted_points(3:4, :)

      top_base = abs(top_points(2, 1) - top_points(1, 1))
      bottom_base = abs(bottom_points(2, 1) - bottom_points(1, 1))

      height = (top_points(1, 2) + top_points(2, 2))/2.0d0 - (bottom_points(1, 2) + bottom_points(2, 2))/2.0d0

      area = 0.5d0*(top_base + bottom_base)*height

   contains
      subroutine swap(a, b)
         real(8), intent(inout) :: a, b
         real(8) :: temp
         temp = a
         a = b
         b = temp
      end subroutine swap
   end function trapezoid_area

   function find_adjacent(p1, p2) result(adjacent_points)
      implicit none
      integer, intent(in) :: p1, p2
      integer :: grid(2, 2) = reshape([1, 2, 3, 4], [2, 2])
      integer :: adjacent_points(2)
      integer :: i, j, idx

      idx = 1

      do i = 1, 2
         do j = 1, 2
            if (grid(i, j) /= p1 .and. grid(i, j) /= p2) then
               adjacent_points(idx) = grid(i, j)
               idx = idx + 1
            end if
         end do
      end do

   end function find_adjacent

   function triangle_area(point1, point2, point3) result(area)
      implicit none
      real(8), intent(in) :: point1(2), point2(2), point3(2)
      real(8) :: area, base, height

      base = sqrt((point2(1) - point1(1))**2 + (point2(2) - point1(2))**2)

      height = sqrt((point2(1) - point3(1))**2 + (point2(2) - point3(2))**2)

      area = 0.5d0*base*height
   end function triangle_area

   function line_circle_intersection(x1, y1, x2, y2, a, b, r) result(points)
      use, intrinsic :: ieee_arithmetic
      implicit none

      real(8), intent(in) :: x1, y1, x2, y2, a, b, r
      real(8), allocatable :: points(:, :)
      real(8) :: dx, dy, coef_a, coef_b, coef_c, disc, sqrt_disc, t1, t2
      real(8), allocatable :: temp_points(:, :)
      integer :: count
      real(8), parameter :: tol = 1.0d-10

      dx = x2 - x1
      dy = y2 - y1

      coef_a = dx*dx + dy*dy

      ! 线段退化为点
      if (coef_a < tol) then
         if (abs((x1 - a)**2 + (y1 - b)**2 - r**2) < tol) then
            allocate (points(1, 2))
            points(1, 1) = x1
            points(1, 2) = y1
         else
            allocate (points(1, 2))
            points = reshape([ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)], [1, 2])

         end if
         return
      end if

      coef_b = 2.0d0*(dx*(x1 - a) + dy*(y1 - b))
      coef_c = (x1 - a)**2 + (y1 - b)**2 - r**2
      disc = coef_b**2 - 4.0d0*coef_a*coef_c

      if (disc < -tol) then
         allocate (points(1, 2))
         points = reshape([ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)], [1, 2])

         return
      end if

      allocate (temp_points(2, 2))
      count = 0

      if (disc >= -tol) then
         sqrt_disc = sqrt(max(disc, 0.0d0))
         t1 = (-coef_b - sqrt_disc)/(2.0d0*coef_a)
         t2 = (-coef_b + sqrt_disc)/(2.0d0*coef_a)

         if (t1 >= -tol .and. t1 <= 1.0d0 + tol) then
            count = count + 1
            temp_points(count, 1) = x1 + t1*dx
            temp_points(count, 2) = y1 + t1*dy
         end if
         if (abs(t2 - t1) > tol .and. t2 >= -tol .and. t2 <= 1.0d0 + tol) then
            count = count + 1
            temp_points(count, 1) = x1 + t2*dx
            temp_points(count, 2) = y1 + t2*dy
         end if
      end if

      if (count == 0) then
         allocate (points(1, 2))
         points = reshape([ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)], [1, 2])

      else
         allocate (points(count, 2))
         points = temp_points(1:count, :)
         call remove_duplicates(points, tol)
      end if

   end function line_circle_intersection

   subroutine remove_duplicates(points, tol)
      implicit none
      real(8), dimension(:, :), intent(inout) :: points
      real(8), intent(in) :: tol
      integer :: i, j, n
      real(8) :: diff

      n = size(points, 1)
      do i = 1, n - 1
         do j = i + 1, n
            diff = maxval(abs(points(i, :) - points(j, :)))
            if (diff < tol) then
               points(j, :) = points(n, :)
               points = points(1:n - 1, :)
               n = n - 1
               exit
            end if
         end do
      end do
   end subroutine remove_duplicates

   function lineSegmentIntersection(x1, y1, x2, y2, x3, y3, x4, y4) result(intersection)
      use, intrinsic :: ieee_arithmetic
      implicit none
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(8) :: dx1, dy1, dx2, dy2, denom, t, u, px, py
      real(8), allocatable :: intersection(:)

      dx1 = x2 - x1
      dy1 = y2 - y1
      dx2 = x4 - x3
      dy2 = y4 - y3

      denom = dx1*dy2 - dy1*dx2

      if (abs(denom) > 1e-10) then  ! Avoid division by zero
         t = ((x3 - x1)*dy2 - (y3 - y1)*dx2)/denom
         u = ((x3 - x1)*dy1 - (y3 - y1)*dx1)/denom

         if (t >= 0.0d0 .and. t <= 1.0d0 .and. u >= 0.0d0 .and. u <= 1.0d0) then
            px = x1 + t*dx1
            py = y1 + t*dy1
            allocate (intersection(2))
            intersection(1) = px
            intersection(2) = py
         else
            allocate (intersection(2))
            intersection = [ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)]
         end if
      else
         allocate (intersection(2))
         intersection = [ieee_value(0.0_8, ieee_quiet_nan), ieee_value(0.0_8, ieee_quiet_nan)]
      end if
   end function lineSegmentIntersection

   logical function PointInSector(px, py, centre, innerR, outerR, startAngle, endAngle)
      implicit none
      real(8), intent(in) :: px, py
      real(8), intent(in) :: centre(2)
      real(8), intent(in) :: innerR, outerR, startAngle, endAngle
      real(8) :: distance, angle
      real(8) :: modStartAngle, modEndAngle

      modStartAngle = mod(startAngle, 2.0d0*3.141592653589793d0)
      modEndAngle = mod(endAngle, 2.0d0*3.141592653589793d0)

      distance = sqrt((px - centre(1))**2 + (py - centre(2))**2)

      angle = atan2(py - centre(2), px - centre(1))

      if (angle < 0.0d0) then
         angle = angle + 2.0d0*3.141592653589793d0
      end if

      if (distance >= innerR .and. distance <= outerR .and. &
          ((modStartAngle <= modEndAngle .and. angle >= modStartAngle .and. angle <= modEndAngle) .or. &
           (modStartAngle > modEndAngle .and. (angle >= modStartAngle .or. angle <= modEndAngle)))) then
         PointInSector = .true.
      else
         PointInSector = .false.
      end if

   end function PointInSector

   subroutine Area_Calculation(step)
      real(8)::sector1_startAngle
      real(8)::sector1_endAngle
      real(8)::sector2_startAngle
      real(8)::sector2_endAngle
      real(8)::add_angle
      real(8)::t
      REAL(8) :: innerStartX1, innerStartY1
      REAL(8) :: innerEndX1, innerEndY1
      REAL(8) :: OuterStartX1, OuterStartY1
      REAL(8) :: OuterEndX1, OuterEndY1
      REAL(8) :: innerStartX2, innerStartY2
      REAL(8) :: innerEndX2, innerEndY2
      REAL(8) :: OuterStartX2, OuterStartY2
      REAL(8) :: OuterEndX2, OuterEndY2
      integer::step
      INTEGER :: i, j
      INTEGER :: ni, nj
      t = step*rotate_dt
      add_angle = w*t
      show = 0.0d0
      !sectror1
      sector1_startAngle = previous_startangle1 + add_angle
      sector1_endAngle = angle1 + add_angle
      !sector2
      sector2_startAngle = previous_startangle2 + add_angle
      sector2_endAngle = angle2 + add_angle
      !endPoints-sector1
      innerStartX1 = centrePoint(1) + inR*cos(sector1_startAngle); 
      innerStartY1 = centrePoint(2) + inR*sin(sector1_startAngle); 
      innerEndX1 = centrePoint(1) + inR*cos(sector1_endAngle); 
      innerEndY1 = centrePoint(2) + inR*sin(sector1_endAngle); 
      OuterStartX1 = centrePoint(1) + outR*cos(sector1_startAngle); 
      OuterStartY1 = centrePoint(2) + outR*sin(sector1_startAngle); 
      OuterEndX1 = centrePoint(1) + outR*cos(sector1_endAngle); 
      OuterEndY1 = centrePoint(2) + outR*sin(sector1_endAngle); 
      !endPoints-sector2
      innerStartX2 = centrePoint(1) + inR*cos(sector2_startAngle); 
      innerStartY2 = centrePoint(2) + inR*sin(sector2_startAngle); 
      innerEndX2 = centrePoint(1) + inR*cos(sector2_endAngle); 
      innerEndY2 = centrePoint(2) + inR*sin(sector2_endAngle); 
      OuterStartX2 = centrePoint(1) + outR*cos(sector2_startAngle); 
      OuterStartY2 = centrePoint(2) + outR*sin(sector2_startAngle); 
      OuterEndX2 = centrePoint(1) + outR*cos(sector2_endAngle); 
      OuterEndY2 = centrePoint(2) + outR*sin(sector2_endAngle); 
      !loop every grid
      ni = size(Xgrid, 2)
      nj = size(Xgrid, 1)
      DO i = 1, nj - 1
         DO j = 1, ni - 1

         end do
      end do

   end subroutine Area_Calculation

end module Grid_Area

PROGRAM learn
   USE Grid_Area
   IMPLICIT NONE

   REAL(8) :: rinR, routR, rBlade, rTank, rw
   REAL(8) :: rAreaDT, rRotateDT, rtmax
   REAL(8) :: rCentre(2), rdx, rdy, rAngle1, rAngle2
   integer ::q, nsteps
   ! user input
   PRINT *, "Enter inner radius:"
   READ (*, *) rinR
   PRINT *, "Enter outer radius:"
   READ (*, *) routR
   PRINT *, "Enter blade length:"
   READ (*, *) rBlade
   PRINT *, "Enter tank size:"
   READ (*, *) rTank
   PRINT *, "Enter angular velocity:"
   READ (*, *) rw
   PRINT *, "Enter area time step:"
   READ (*, *) rAreaDT
   PRINT *, "Enter rotate time step:"
   READ (*, *) rRotateDT
   PRINT *, "Enter total rotate time:"
   READ (*, *) rtmax
   PRINT *, "Enter center point coordinates (x y):"
   READ (*, *) rCentre(1), rCentre(2)
   PRINT *, "Enter grid dx:"
   READ (*, *) rdx
   PRINT *, "Enter grid dy:"
   READ (*, *) rdy
   PRINT *, "Enter blade 1 initial angle:"
   READ (*, *) rAngle1
   PRINT *, "Enter blade 2 initial angle:"
   READ (*, *) rAngle2

   CALL set_geometry(rinR, routR, rBlade, rTank, rw, rAreaDT, rRotateDT, rtmax, &
                     rCentre, rdx, rdy, rAngle1, rAngle2)

   CALL print_geometry
   !calculate grid
   nsteps = INT(aint(rtmax/rRotateDT))
   do q = 0, nsteps
      call Area_Calculation(q)
   end do

END PROGRAM learn

