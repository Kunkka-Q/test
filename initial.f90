    SUBROUTINE init_Grid_Area()

        ! leaving inactive if no parameters specified
        has_Grid_Area = .FALSE.
        IF (.NOT. fort7%exists("/flow/Grid_Area")) RETURN

        ! retrieving rotation rate vector from parameters.json
        

        ! display obtained parameters
        IF (myid == 0) THEN
            WRITE(*, '(start calculate Grid)')
        END IF

        ! set active
        has_Grid_Area = .TRUE.

    END SUBROUTINE init_Grid_Area