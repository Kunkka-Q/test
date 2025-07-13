MODULE Grid_Area
   USE core_mod
   IMPLICIT NONE(type, external)
   PRIVATE

   !public
   PUBLIC :: init_Grid_Area, finish_Grid_Area, Area_Calculation
   REAL(realk), ALLOCATABLE, PRIVATE :: show(:, :)

   !  private
   REAL(realk), private :: rCentre(2)
   REAL(realk), private :: rTank
   REAL(realk), private :: rw
   REAL(realk), private :: rtmax
   LOGICAL, PROTECTED :: has_Grid_Area = .FALSE.
   REAL(realk), private :: inR
   REAL(realk), private :: outR
   REAL(realk), private :: BladeLength
   REAL(realk), private :: TankSize
   REAL(realk), private :: w
   REAL(realk), private :: area_dt
   REAL(realk), private :: rotate_dt
   REAL(realk), private :: tmax
   REAL(realk), private :: centrePoint(2)
   REAL(realk), private :: previous_startangle1
   REAL(realk), private :: previous_startangle2
   REAL(realk), private :: angle1
   REAL(realk), private :: angle2
   !grid - 移除内部网格生成，使用外部网格
   !PRIVATE :: generate_mesh
   !REAL(realk), ALLOCATABLE, PRIVATE :: Xgrid(:, :), Ygrid(:, :)
contains
   SUBROUTINE init_Grid_Area()

      ! leaving inactive if no parameters specified
      has_Grid_Area = .FALSE.
      IF (.NOT. fort7%exists("/flow/Grid_Area_Calculation")) RETURN   

      ! retrieving rotation rate vector from parameters.json
      !centre position
       CALL fort7%get_array("/flow/Grid_Area_Calculation/rCentre", rCentre)
      !tank size
       CALL fort7%get_array("/flow/Grid_Area_Calculation/rTank", rTank)
      !angular velocity
       CALL fort7%get_array("/flow/Grid_Area_Calculation/rw", rw)
      !total simulation time
       CALL fort7%get_array("/flow/Grid_Area_Calculation/rtmax", rtmax)
      !inner radius
       CALL fort7%get_array("/flow/Grid_Area_Calculation/inR", inR)
      !outer radius
       CALL fort7%get_array("/flow/Grid_Area_Calculation/outR", outR)
      !blade length
       CALL fort7%get_array("/flow/Grid_Area_Calculation/BladeLength", BladeLength)
      !area calculation time step
       CALL fort7%get_array("/flow/Grid_Area_Calculation/area_dt", area_dt)
      !rotation time step
       CALL fort7%get_array("/flow/Grid_Area_Calculation/rotate_dt", rotate_dt)
      !initial angles
       CALL fort7%get_array("/flow/Grid_Area_Calculation/previous_startangle1", previous_startangle1)
       CALL fort7%get_array("/flow/Grid_Area_Calculation/previous_startangle2", previous_startangle2)
      !calculate angles
       angle1 = previous_startangle1 + area_dt*rw
       angle2 = previous_startangle2 + area_dt*rw
      ! display obtained parameters
      IF (myid == 0) THEN
            WRITE(*, '("Calculate_Grid_Area TERM:")')
            WRITE(*, '(2X, "CentrePosition: ", 2(G0, 1X))') rCentre
            WRITE(*, '(2X, "TankSize: ", 1(G0, 1X))') rTank
            WRITE(*, '(2X, "AngularVelocity: ", 1(G0, 1X))') rw
            WRITE(*, '(2X, "total_simulation_time: ", 1(G0, 1X))') rtmax
            WRITE(*, '(2X, "InnerRadius: ", 1(G0, 1X))') inR
            WRITE(*, '(2X, "OuterRadius: ", 1(G0, 1X))') outR
            WRITE(*, '(2X, "BladeLength: ", 1(G0, 1X))') BladeLength
            WRITE(*, '(2X, "AreaTimeStep: ", 1(G0, 1X))') area_dt
            WRITE(*, '(2X, "RotateTimeStep: ", 1(G0, 1X))') rotate_dt
            WRITE(*, '(2X, "InitialAngle1: ", 1(G0, 1X))') previous_startangle1
            WRITE(*, '(2X, "InitialAngle2: ", 1(G0, 1X))') previous_startangle2
            WRITE(*, '()')
        END IF

      ! set active
      has_Grid_Area = .TRUE.

   END SUBROUTINE init_Grid_Area

   SUBROUTINE finish_Grid_Area

      ! revoking activity
      has_Grid_Area = .FALSE.

      RETURN

   END SUBROUTINE finish_Grid_Area

   ! 移除generate_mesh子程序，改为使用外部网格信息



   
   FUNCTION trapezoid_area(points) result(area)
      implicit none
      real(realk), intent(in) :: points(4, 2)
      real(realk) :: area
      real(realk) :: sorted_points(4, 2)
      real(realk) :: top_points(2, 2), bottom_points(2, 2)
      real(realk) :: top_base, bottom_base, height
      integer(intk) :: i, j, max_idx, k
      real(realk) :: y_values(4)

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

      height = (top_points(1, 2) + top_points(2, 2))/2.0_realk - (bottom_points(1, 2) + bottom_points(2, 2))/2.0_realk

      area = 0.5_realk*(top_base + bottom_base)*height

   contains
      SUBROUTINE swap(a, b)
         real(realk), intent(inout) :: a, b
         real(realk) :: temp
         temp = a
         a = b
         b = temp
      END SUBROUTINE swap
   END FUNCTION trapezoid_area

   FUNCTION find_adjacent(p1, p2) result(adjacent_points)
      implicit none
      integer(intk), intent(in) :: p1, p2
      integer(intk) :: grid(2, 2) = reshape([1, 2, 3, 4], [2, 2])
      integer(intk) :: adjacent_points(2)
      integer(intk) :: i, j, idx

      idx = 1

      do i = 1, 2
         do j = 1, 2
            if (grid(i, j) /= p1 .and. grid(i, j) /= p2) then
               adjacent_points(idx) = grid(i, j)
               idx = idx + 1
            end if
         end do
      end do

   END FUNCTION find_adjacent

   FUNCTION triangle_area(point1, point2, point3) result(areaa)
      implicit none
      real(realk), intent(in) :: point1(2), point2(2), point3(2)
      real(realk) :: areaa, base, height

      base = sqrt((point2(1) - point1(1))**2 + (point2(2) - point1(2))**2)

      height = sqrt((point2(1) - point3(1))**2 + (point2(2) - point3(2))**2)

      areaa = 0.5_realk*base*height
   END FUNCTION triangle_area

   FUNCTION line_circle_intersection(x1, y1, x2, y2, a, b, r) result(points)
      use, intrinsic :: ieee_arithmetic
      implicit none

      real(realk), intent(in) :: x1, y1, x2, y2, a, b, r
      real(realk), allocatable :: points(:, :)
      real(realk) :: dx, dy, coef_a, coef_b, coef_c, disc, sqrt_disc, t1, t2
      real(realk), allocatable :: temp_points(:, :)
      integer(intk) :: count
      real(realk), parameter :: tol = 1.0_realk-10

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
            points = reshape([ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)], [1, 2])

         end if
         return
      end if

      coef_b = 2.0_realk*(dx*(x1 - a) + dy*(y1 - b))
      coef_c = (x1 - a)**2 + (y1 - b)**2 - r**2
      disc = coef_b**2 - 4.0_realk*coef_a*coef_c

      if (disc < -tol) then
         allocate (points(1, 2))
         points = reshape([ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)], [1, 2])

         return
      end if

      allocate (temp_points(2, 2))
      count = 0

      if (disc >= -tol) then
         sqrt_disc = sqrt(max(disc, 0.0_realk))
         t1 = (-coef_b - sqrt_disc)/(2.0_realk*coef_a)
         t2 = (-coef_b + sqrt_disc)/(2.0_realk*coef_a)

         if (t1 >= -tol .and. t1 <= 1.0_realk + tol) then
            count = count + 1
            temp_points(count, 1) = x1 + t1*dx
            temp_points(count, 2) = y1 + t1*dy
         end if
         if (abs(t2 - t1) > tol .and. t2 >= -tol .and. t2 <= 1.0_realk + tol) then
            count = count + 1
            temp_points(count, 1) = x1 + t2*dx
            temp_points(count, 2) = y1 + t2*dy
         end if
      end if

      if (count == 0) then
         allocate (points(1, 2))
         points = reshape([ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)], [1, 2])

      else
         allocate (points(count, 2))
         points = temp_points(1:count, :)
      end if

   END FUNCTION line_circle_intersection

   FUNCTION lineSegmentIntersection(x1, y1, x2, y2, x3, y3, x4, y4) result(intersection)
      use, intrinsic :: ieee_arithmetic
      implicit none
      real(realk), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(realk) :: dx1, dy1, dx2, dy2, denom, t, u, px, py
      real(realk), allocatable :: intersection(:)

      dx1 = x2 - x1
      dy1 = y2 - y1
      dx2 = x4 - x3
      dy2 = y4 - y3

      denom = dx1*dy2 - dy1*dx2

      if (abs(denom) > 1e-10) then  ! Avoid division by zero
         t = ((x3 - x1)*dy2 - (y3 - y1)*dx2)/denom
         u = ((x3 - x1)*dy1 - (y3 - y1)*dx1)/denom

         if (t >= 0.0_realk .and. t <= 1.0_realk .and. u >= 0.0_realk .and. u <= 1.0_realk) then
            px = x1 + t*dx1
            py = y1 + t*dy1
            allocate (intersection(2))
            intersection(1) = px
            intersection(2) = py
         else
            allocate (intersection(2))
            intersection = [ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)]
         end if
      else
         allocate (intersection(2))
         intersection = [ieee_value(0.0_realk, ieee_quiet_nan), ieee_value(0.0_realk, ieee_quiet_nan)]
      end if
   END FUNCTION lineSegmentIntersection

   LOGICAL FUNCTION PointInSector(px, py, centre, innerR, outerR, startAngle, endAngle)
      implicit none
      real(realk), intent(in) :: px, py
      real(realk), intent(in) :: centre(2)
      real(realk), intent(in) :: innerR, outerR, startAngle, endAngle
      real(realk) :: distance, angle
      real(realk) :: modStartAngle, modEndAngle

      modStartAngle = mod(startAngle, 2.0_realk*3.141592653589793_realk)
      modEndAngle = mod(endAngle, 2.0_realk*3.141592653589793_realk)

      distance = sqrt((px - centre(1))**2 + (py - centre(2))**2)

      angle = atan2(py - centre(2), px - centre(1))

      if (angle < 0.0_realk) then
         angle = angle + 2.0_realk*3.141592653589793_realk
      end if

      if (distance >= innerR .and. distance <= outerR .and. &
          ((modStartAngle <= modEndAngle .and. angle >= modStartAngle .and. angle <= modEndAngle) .or. &
           (modStartAngle > modEndAngle .and. (angle >= modStartAngle .or. angle <= modEndAngle)))) then
         PointInSector = .true.
      else
         PointInSector = .false.
      end if

   END FUNCTION PointInSector

   ! 修改Area_Calculation子程序，使其与cori的网格处理方式一致
   subroutine Area_Calculation(step, show_result)
      use, intrinsic :: ieee_arithmetic
      INTEGER(intk), INTENT(IN) :: step
      REAL(realk), INTENT(OUT), ALLOCATABLE :: show_result(:, :)
      real(realk)::sector1_startAngle
      real(realk)::sector1_endAngle
      real(realk)::sector2_startAngle
      real(realk)::sector2_endAngle
      real(realk)::add_angle
      real(realk)::t
      REAL(realk) :: innerStartX1, innerStartY1
      REAL(realk) :: innerEndX1, innerEndY1
      REAL(realk) :: OuterStartX1, OuterStartY1
      REAL(realk) :: OuterEndX1, OuterEndY1
      REAL(realk) :: innerStartX2, innerStartY2
      REAL(realk) :: innerEndX2, innerEndY2
      REAL(realk) :: OuterStartX2, OuterStartY2
      REAL(realk) :: OuterEndX2, OuterEndY2
      INTEGER(intk) :: i, j, k, igrid
      INTEGER(intk) :: ni, nj, kk, jj, ii
      real(realk) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(realk) :: xpoints(4), ypoints(4)
      integer(intk) :: q(4)
      integer(intk) :: mapping(4, 2)
      integer(intk) :: m
      integer(intk) :: idx
      integer(intk) :: New_m, m_count
      integer(intk) :: New_idx(2)
      integer(intk) :: adjacent_points(2)
      real(realk) :: orthoPoint(2)
      integer(intk) :: calcupoint(2)
      integer(intk) :: New_New_m, New_New_idx
      real(realk) :: New_orthoPoint(2)
      integer(intk) :: New_calcupoint(2)
      real(realk) :: coordinateA(2), coordinateB(2), coordinateC(2), coordinateD(2), coordinateE(2), coordinateF(2)
      real(realk) :: coordinate1(2), coordinate2(2), coordinate3(2), coordinate4(2)
      real(realk) :: coordinate7(2), coordinate8(2), coordinate9(2), coordinate10(2)
      real(realk) :: coordinate13(2), coordinate14(2), coordinate15(2), coordinate16(2)
      real(realk) :: coordinate19(2), coordinate20(2), coordinate21(2), coordinate22(2)
      real(realk) :: coordinate25(2), coordinate26(2), coordinate27(2), coordinate28(2)
      real(realk) :: coordinate31(2), coordinate32(2), coordinate33(2), coordinate34(2)
      real(realk), allocatable :: coordinate5(:, :), coordinate6(:, :)
      real(realk), allocatable :: coordinate11(:, :), coordinate12(:, :)
      real(realk), allocatable :: coordinate17(:, :), coordinate18(:, :)
      real(realk), allocatable :: coordinate23(:, :), coordinate24(:, :)
      real(realk), allocatable :: coordinate29(:, :), coordinate30(:, :)
      real(realk), allocatable :: coordinate35(:, :), coordinate36(:, :)
      real(realk) ::area, New_area
      real(realk), dimension(4, 2) :: trapezoid_Point
      real(realk) :: trapezoid_Area_calculated
      
      ! 网格坐标相关变量
      REAL(realk), ALLOCATABLE :: Xgrid(:, :), Ygrid(:, :)
      REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
      REAL(realk) :: dx, dy
      
      ! 网格字段指针
      TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
      REAL(realk), POINTER, CONTIGUOUS :: dx_ptr(:), dy_ptr(:), dz_ptr(:)
      REAL(realk), POINTER, CONTIGUOUS :: ddx_ptr(:), ddy_ptr(:), ddz_ptr(:)

      !check status
      IF (.NOT. has_Grid_Area) RETURN

      t = step*rotate_dt
      add_angle = w*t
      
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
      
      ! 修改为与cori一致的网格处理方式
      ! 遍历所有网格块（类似cori的方式）
      DO i = 1, nmygrids
         igrid = mygrids(i)
         
         ! 获取网格尺寸（使用MGLET的接口）
         CALL get_mgdims(kk, jj, ii, igrid)
         
         ! 对于2D网格，对应关系：ni=ii, nj=jj
         ni = ii  ! x方向网格点数
         nj = jj  ! y方向网格点数
         
         ! 获取网格边界框
         CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
         
         ! 获取基本网格间距
         CALL get_field(dx_f, "DX")    ! 获取 dx 字段
         CALL get_field(dy_f, "DY")    ! 获取 dy 字段  
         CALL get_field(dz_f, "DZ")    ! 获取 dz 字段

         ! 获取网格间距的倒数
         CALL get_field(ddx_f, "DDX")  ! 获取 ddx 字段
         CALL get_field(ddy_f, "DDY")  ! 获取 ddy 字段
         CALL get_field(ddz_f, "DDZ")  ! 获取 ddz 字段

         ! 获取指针
         CALL dx_f%get_ptr(dx_ptr, igrid)   ! dx_ptr(:) - 基本网格间距
         CALL dy_f%get_ptr(dy_ptr, igrid)   ! dy_ptr(:) - 基本网格间距
         CALL dz_f%get_ptr(dz_ptr, igrid)   ! dz_ptr(:) - 基本网格间距

         CALL ddx_f%get_ptr(ddx_ptr, igrid) ! ddx_ptr(:) - 网格间距倒数
         CALL ddy_f%get_ptr(ddy_ptr, igrid) ! ddy_ptr(:) - 网格间距倒数
         CALL ddz_f%get_ptr(ddz_ptr, igrid) ! ddz_ptr(:) - 网格间距倒数

         ! 分配网格坐标数组
         ALLOCATE(Xgrid(nj, ni), Ygrid(nj, ni))
         
         ! 直接使用从外部获取的 dx, dy
         DO j = 1, nj
            DO k = 1, ni
               Xgrid(j, k) = minx + REAL(k - 1, realk) * dx_ptr(k)  ! 使用 dx_ptr(k)
               Ygrid(j, k) = miny + REAL(j - 1, realk) * dy_ptr(j)  ! 使用 dy_ptr(j)
            END DO
         END DO
         
         ! 分配show数组（如果还没有分配）
         IF (.NOT. ALLOCATED(show)) THEN
            ALLOCATE(show(nj-1, ni-1))
         END IF
         show = 0.0_realk
         
         ! 遍历网格单元（使用与cori一致的逻辑：从第三个网格点开始）
         DO j = 3, nj-2
            DO k = 3, ni-2
               ! 获取网格单元顶点
               x1 = Xgrid(j, k); y1 = Ygrid(j, k)
               x2 = Xgrid(j, k + 1); y2 = Ygrid(j, k + 1)
               x3 = Xgrid(j + 1, k); y3 = Ygrid(j + 1, k)
               x4 = Xgrid(j + 1, k + 1); y4 = Ygrid(j + 1, k + 1)
               
               xpoints = [x1, x2, x3, x4]
               ypoints = [y1, y2, y3, y4]

               q = [0, 0, 0, 0]

               mapping(1, :) = [2, 3]
               mapping(2, :) = [1, 4]
               mapping(3, :) = [1, 4]
               mapping(4, :) = [2, 3]
               Do m = 1, 4
                  if (PointInSector(xpoints(m), ypoints(m), centrePoint, inR, outR, sector1_startAngle, sector1_endAngle) .or. &
                      PointInSector(xpoints(m), ypoints(m), centrePoint, inR, outR, sector2_startAngle, sector2_endAngle)) then
                     q(m) = 1
                  else
                     q(m) = 0
                  end if
               end do
               
               ! 面积计算逻辑保持不变
               !4 points
               if (sum(q) == 4) then
                  show(j, k) = 1
                  !3 points
               elseif (sum(q) == 3) then
                  do m = 1, 4
                     if (q(m) == 0) then
                        idx = m
                        exit
                     end if
                  end do
                  orthoPoint(1) = xpoints(idx)
                  orthoPoint(2) = ypoints(idx)
                  calcupoint = mapping(idx, :)
                  !Calculate A Point
                  coordinate1=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate1(1))) then
                     coordinateA = coordinate1
                  else
                     coordinate2=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate2(1))) then
                        coordinateA = coordinate2
                     else
                        coordinate3=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate3(1))) then
                           coordinateA = coordinate3
                        else
                         coordinate4=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(1)),ypoints(calcupoint(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate4(1))) then
                              coordinateA = coordinate4
                           else
                              coordinate5 = line_circle_intersection(xpoints(calcupoint(1)), ypoints(calcupoint(1)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), inR)
                              if (.not. ieee_is_nan(coordinate5(1, 1))) then
                                 coordinateA = coordinate5(1, :)
                              else
                                 coordinate6 = line_circle_intersection(xpoints(calcupoint(1)), ypoints(calcupoint(1)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate6(1, 1))) then
                                    coordinateA = coordinate6(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !Calculate B Point
                  coordinate7=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate7(1))) then
                     coordinateB = coordinate7
                  else
                     coordinate8=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate8(1))) then
                        coordinateB = coordinate8
                     else
                     coordinate9=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate9(1))) then
                           coordinateB = coordinate9
                        else
                     coordinate10=lineSegmentIntersection(xpoints(idx),ypoints(idx),xpoints(calcupoint(2)),ypoints(calcupoint(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate10(1))) then
                              coordinateB = coordinate10
                           else
                              coordinate11 = line_circle_intersection(xpoints(calcupoint(2)), ypoints(calcupoint(2)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), inR)
                              if (.not. ieee_is_nan(coordinate11(1, 1))) then
                                 coordinateB = coordinate11(1, :)
                              else
                                 coordinate12 = line_circle_intersection(xpoints(calcupoint(2)), ypoints(calcupoint(2)), xpoints(idx), ypoints(idx), centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate12(1, 1))) then
                                    coordinateB = coordinate12(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !calculate cut area
                  area = triangle_area(coordinateA, orthoPoint, coordinateB)
                  show(j, k) = ((dx_ptr(k)*dy_ptr(j)) - area)/(dx_ptr(k)*dy_ptr(j))
                  ! 2 points
               elseif (sum(q) == 2) then
                  m_count = 0
                  do New_m = 1, 4
                     if (q(New_m) == 0) then
                        m_count = m_count + 1
                        New_idx(m_count) = New_m
                        if (m_count == 2) exit
                     end if
                  end do
                  adjacent_points = find_adjacent(New_idx(1), New_idx(2))
                  !Calculate C Point
                  coordinate13=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate13(1))) then
                     coordinateC = coordinate13
                  else
                     coordinate14=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate14(1))) then
                        coordinateC = coordinate14
                     else
                     coordinate15=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate15(1))) then
                           coordinateC = coordinate15
                        else
                     coordinate16=lineSegmentIntersection(xpoints(New_idx(1)),ypoints(New_idx(1)),xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate16(1))) then
                              coordinateC = coordinate16
                           else
                         coordinate17=line_circle_intersection(xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),xpoints(New_idx(1)),ypoints(New_idx(1)),centrePoint(1), centrePoint(2), inR)
                              if (.not. ieee_is_nan(coordinate17(1, 1))) then
                                 coordinateC = coordinate17(1, :)
                              else
                             coordinate18=line_circle_intersection(xpoints(adjacent_points(1)),ypoints(adjacent_points(1)),xpoints(New_idx(1)),ypoints(New_idx(1)),centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate18(1, 1))) then
                                    coordinateC = coordinate18(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !Calculate D Point
                  coordinate19=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate19(1))) then
                     coordinateD = coordinate19
                  else
                     coordinate20=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate20(1))) then
                        coordinateD = coordinate20
                     else
                     coordinate21=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate21(1))) then
                           coordinateD = coordinate21
                        else
                     coordinate22=lineSegmentIntersection(xpoints(New_idx(2)),ypoints(New_idx(2)),xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate22(1))) then
                              coordinateD = coordinate22
                           else
                         coordinate23=line_circle_intersection(xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),xpoints(New_idx(2)),ypoints(New_idx(2)),centrePoint(1), centrePoint(2), inR)
                              if (.not. ieee_is_nan(coordinate23(1, 1))) then
                                 coordinateD = coordinate23(1, :)
                              else
                             coordinate24=line_circle_intersection(xpoints(adjacent_points(2)),ypoints(adjacent_points(2)),xpoints(New_idx(2)),ypoints(New_idx(2)),centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate24(1, 1))) then
                                    coordinateD = coordinate24(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !Calculate trapzoid area
                  trapezoid_Point(1, 1) = xpoints(New_idx(1))
                  trapezoid_Point(1, 2) = ypoints(New_idx(1))
                  trapezoid_Point(2, 1) = xpoints(New_idx(2))
                  trapezoid_Point(2, 2) = ypoints(New_idx(2))
                  trapezoid_Point(3, 1) = coordinateC(1)
                  trapezoid_Point(3, 2) = coordinateC(2)
                  trapezoid_Point(4, 1) = coordinateD(1)
                  trapezoid_Point(4, 2) = coordinateD(2)
                  trapezoid_Area_calculated = trapezoid_area(trapezoid_Point)
                  show(j, k) = ((dx_ptr(k)*dy_ptr(j)) - trapezoid_Area_calculated)/(dx_ptr(k)*dy_ptr(j))
               elseif (sum(q) == 1) then
                  do New_New_m = 1, 4
                     if (q(New_New_m) == 1) then
                        New_New_idx = New_New_m
                        exit
                     end if
                  end do
                  New_orthoPoint(1) = xpoints(New_New_idx)
                  New_orthoPoint(2) = ypoints(New_New_idx)
                  New_calcupoint = mapping(New_New_idx, :)
                  !Calculate E point
                  coordinate25=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
                  if (.not. ieee_is_nan(coordinate25(1))) then
                     coordinateE = coordinate25
                  else
                      coordinate26=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                     if (.not. ieee_is_nan(coordinate26(1))) then
                        coordinateE = coordinate26
                     else
                     coordinate27=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                        if (.not. ieee_is_nan(coordinate27(1))) then
                           coordinateE = coordinate27
                        else
                     coordinate28=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                           if (.not. ieee_is_nan(coordinate28(1))) then
                              coordinateE = coordinate28
                           else
                    coordinate29=line_circle_intersection(xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)),xpoints(New_New_idx),ypoints(New_New_idx), centrePoint(1), centrePoint(2), inR) 
                              if (.not. ieee_is_nan(coordinate29(1, 1))) then
                                 coordinateE = coordinate29(1, :)
                              else
                                  coordinate30=line_circle_intersection(xpoints(New_calcupoint(1)),ypoints(New_calcupoint(1)), xpoints(New_New_idx),ypoints(New_New_idx),centrePoint(1), centrePoint(2), outR)
                                 if (.not. ieee_is_nan(coordinate30(1, 1))) then
                                    coordinateE = coordinate30(1, :)
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
                  !Calculate F point
                  coordinate31=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerStartX1,innerStartY1,OuterStartX1,OuterStartY1)
               if (.not. ieee_is_nan(coordinate31(1))) then
                  coordinateF = coordinate31
               else
                   coordinate32=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerEndX1,innerEndY1,OuterEndX1,OuterEndY1)
                  if (.not. ieee_is_nan(coordinate32(1))) then
                     coordinateF = coordinate32
                  else
                  coordinate33=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerStartX2,innerStartY2,OuterStartX2,OuterStartY2)
                     if (.not. ieee_is_nan(coordinate33(1))) then
                        coordinateF = coordinate33
                     else
                  coordinate34=lineSegmentIntersection(xpoints(New_New_idx),ypoints(New_New_idx),xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),innerEndX2,innerEndY2,OuterEndX2,OuterEndY2)
                        if (.not. ieee_is_nan(coordinate34(1))) then
                           coordinateF = coordinate34
                        else
                 coordinate35=line_circle_intersection(xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)),xpoints(New_New_idx),ypoints(New_New_idx), centrePoint(1), centrePoint(2), inR) 
                           if (.not. ieee_is_nan(coordinate35(1, 1))) then
                              coordinateF = coordinate35(1, :)
                           else
                               coordinate36=line_circle_intersection(xpoints(New_calcupoint(2)),ypoints(New_calcupoint(2)), xpoints(New_New_idx),ypoints(New_New_idx),centrePoint(1), centrePoint(2), outR)
                              if (.not. ieee_is_nan(coordinate36(1, 1))) then
                                 coordinateF = coordinate36(1, :)
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
                  New_area = triangle_area(coordinateE, New_orthoPoint, coordinateF)
                  show(j, k) = (New_area)/(dx_ptr(k)*dy_ptr(j))
               end if
            END DO
         END DO
          
         ! 释放网格坐标内存
         DEALLOCATE(Xgrid, Ygrid)
      END DO
      
      ! 返回结果
      IF (ALLOCATED(show_result)) DEALLOCATE(show_result)
      ALLOCATE(show_result, SOURCE=show)

   END SUBROUTINE Area_Calculation

END MODULE Grid_Area 