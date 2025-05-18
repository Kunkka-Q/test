program swept_area_simulation
   implicit none
   ! background
   real(8), parameter :: tanksize = 15.7  ! Tank size
   real(8), parameter :: bl = 6.4        ! Blade length
   real(8), parameter :: innerR = 1.45   ! Inner radius
   real(8), parameter :: outR = 7.85     ! Outer radius
   real(8), parameter :: w = 2.83        ! Angular velocity
   real(8), parameter :: dt = 0.4        ! Time for forming swept area
   real(8), parameter :: DT = 0.1        ! Rotate dt
   real(8), parameter :: TMAX = 1.5      ! Total rotate time
   ! User input(change to read--------------------------------------------------------------)
   real(8) :: centre(2)
   real(8) :: X(100, 100), Y(100, 100)
   real(8) :: show(100, 100)
   real(8) :: previous_angle1, previous_angle2, angle1, angle2
   real(8) :: rotation_angle
   integer :: i, j
   integer :: nX, nY
 
   ! Initialize
   centre = [tanksize / 2, tanksize / 2]  ! Centre location
   nX = nint(tanksize / 1.0) + 1  ! Size for X grid
   nY = nint(tanksize / 1.0) + 1  ! Size for Y grid
 
   ! Initialize show matrix with zeros
   show = 0.0
 
   ! Initialize swept angle
   previous_angle1 = 0.0
   previous_angle2 = 3.14159  ! Pi
   angle1 = previous_angle1 + w * dt
   angle2 = angle1 + 3.14159  ! Pi
 
   ! Rotation angle initialization
   rotation_angle = 0.0
 
   ! Generate the X and Y meshgrid
   call generate_grid(X, Y, nX, nY, tanksize)
 
   ! Main loop to simulate swept area
   do while (rotation_angle < TMAX)
      ! Update the rotation angle
      rotation_angle = rotation_angle + DT
 
      ! Calculate the new angle positions
      angle1 = previous_angle1 + w * dt
      angle2 = angle1 + 3.14159  ! Pi
 
      ! Update previous angles
      previous_angle1 = angle1
      previous_angle2 = angle2
 
      ! Add logic for updating show matrix here, based on the new angle positions
 
      ! Print or visualize the results if needed (e.g., print matrix)
      print *, 'Rotation Angle: ', rotation_angle
 
   end do
 
 contains
 
   ! Subroutine to generate meshgrid
   subroutine generate_grid(X, Y, nX, nY, size)
     implicit none
     integer, intent(in) :: nX, nY
     real(8), intent(in) :: size
     real(8), intent(out) :: X(nX, nY), Y(nX, nY)
     integer :: i, j
     real(8) :: dx, dy
     real(8) :: xx, yy
 
     dx = size / (nX - 1)
     dy = size / (nY - 1)
 
     do i = 1, nX
        xx = (i - 1) * dx
        do j = 1, nY
           yy = (j - 1) * dy
           X(i, j) = xx
           Y(i, j) = yy
        end do
     end do
   end subroutine generate_grid
 
 end program swept_area_simulation
 
 !function
 !calculate trapezoid area
 FUNCTION trapezoid_area(points) RESULT(area)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: points(4, 2)
    REAL(8) :: area
    REAL(8) :: sorted_points(4, 2)
    REAL(8) :: top_points(2, 2), bottom_points(2, 2)
    REAL(8) :: top_base, bottom_base, height
    INTEGER :: i, j, max_idx
    REAL(8) :: temp(2)

    sorted_points = points
    DO i = 1, 3
       max_idx = i
       DO j = i + 1, 4
          IF (sorted_points(j, 2) > sorted_points(max_idx, 2)) THEN
             max_idx = j
          END IF
       END DO
       IF (max_idx /= i) THEN
          temp = sorted_points(i, :)
          sorted_points(i, :) = sorted_points(max_idx, :)
          sorted_points(max_idx, :) = temp
       END IF
    END DO

    top_points = sorted_points(1:2, :)
    bottom_points = sorted_points(3:4, :)

    top_base = ABS(top_points(2, 1) - top_points(1, 1))
    bottom_base = ABS(bottom_points(2, 1) - bottom_points(1, 1))

    height = (top_points(1, 2) + top_points(2, 2))/2.0 - &
             (bottom_points(1, 2) + bottom_points(2, 2))/2.0

    area = 0.5*(top_base + bottom_base)*height

 END FUNCTION trapezoid_area
!adjacent point
 FUNCTION find_adjacent(p1, p2) RESULT(adjacent_points)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: p1, p2
    INTEGER :: adjacent_points(2)
    INTEGER :: grid(4) = [1, 2, 3, 4]
    INTEGER :: temp(2), idx, i

    idx = 0
    DO i = 1, 4
       IF (grid(i) /= p1 .AND. grid(i) /= p2) THEN
          idx = idx + 1
          temp(idx) = grid(i)
       END IF
    END DO

    adjacent_points = temp

 END FUNCTION find_adjacent
!calculate the area of a right triangle
 FUNCTION triangle_area(point1, point2, point3) RESULT(area)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: point1(2), point2(2), point3(2)
    REAL(8) :: base, height, area

    base = SQRT((point2(1) - point1(1))**2 + (point2(2) - point1(2))**2)

    height = SQRT((point2(1) - point3(1))**2 + (point2(2) - point3(2))**2)

    area = 0.5*base*height

 END FUNCTION triangle_area
!intersection point
 SUBROUTINE line_circle_intersection(x1, y1, x2, y2, a, b, r, points, npoints)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: x1, y1, x2, y2, a, b, r
    REAL(8), INTENT(OUT) :: points(2, 2)
    INTEGER, INTENT(OUT) :: npoints

    REAL(8) :: dx, dy, seg_length_sq
    REAL(8) :: A, B, C, D, sqrtD
    REAL(8) :: t1, t2, t
    REAL(8), PARAMETER :: tolerance = 1.0D-10
    INTEGER :: i
    LOGICAL :: used(2)

    npoints = 0
    points(:, :) = 0.0D0/0.0D0
    used = .FALSE.

    dx = x2 - x1
    dy = y2 - y1
    seg_length_sq = dx**2 + dy**2

    IF (seg_length_sq < 1.0D-20) THEN
       IF (ABS((x1 - a)**2 + (y1 - b)**2 - r**2) < tolerance) THEN
          points(1, :) = [x1, y1]
          npoints = 1
       END IF
       RETURN
    END IF

    A = dx**2 + dy**2
    B = 2.0D0*(dx*(x1 - a) + dy*(y1 - b))
    C = (x1 - a)**2 + (y1 - b)**2 - r**2
    D = B**2 - 4.0D0*A*C

    IF (D < -tolerance) THEN
       RETURN
    END IF

    IF (ABS(D) <= tolerance) THEN
       sqrtD = 0.0D0
    ELSE
       sqrtD = SQRT(D)
    END IF

    t1 = (-B - sqrtD)/(2.0D0*A)
    t2 = (-B + sqrtD)/(2.0D0*A)

    IF (t1 >= -tolerance .AND. t1 <= 1.0D0 + tolerance) THEN
       t = t1
       points(1, :) = [x1 + t*dx, y1 + t*dy]
       npoints = npoints + 1
       used(1) = .TRUE.
    END IF

    IF (t2 >= -tolerance .AND. t2 <= 1.0D0 + tolerance) THEN
       IF (.NOT. used(1) .OR. ABS(t2 - t1) > tolerance) THEN
          t = t2
          points(npoints + 1, :) = [x1 + t*dx, y1 + t*dy]
          npoints = npoints + 1
          used(2) = .TRUE.
       END IF
    END IF

 END SUBROUTINE line_circle_intersection

