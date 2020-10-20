MODULE degfun
    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = KIND(1.d0)
    
    REAL, PARAMETER :: PI  = ATAN(1.d0) * 4.d0
    REAL, PARAMETER :: RAD = (180.d0/PI)

    CONTAINS
        !Time functions
        SUBROUTINE print_hours(deg)
            REAL :: m, s
            REAL, INTENT(IN) :: deg

            m = (deg-FLOOR(deg))*60.d0
            s = (m-FLOOR(m))*60.d0

            10 FORMAT(I3, "h ", I3, "m ", F4.1, "s")
            WRITE(*,10) FLOOR(deg), FLOOR(m), s
        END SUBROUTINE print_hours
        
        FUNCTION from_hours(h, m, s)
            REAL :: from_hours 
            INTEGER  :: h, m, s

            from_hours = h + (m/60.d0) + (s/3600.d0)
        END FUNCTION from_hours
        
        !General math functions
        FUNCTION normalize_ang(x) !Keeps degree value within proper range
           REAL :: x, normalize_ang
           normalize_ang = x - FLOOR(x/360.d0)*360.d0
       END FUNCTION normalize_ang

        FUNCTION cbrt(x) !Cube root
            REAL :: x, cbrt
            IF (x .GT. 0) THEN
                cbrt = EXP(LOG(x)/3.d0)
            ELSE
                cbrt = -EXP(LOG(ABS(x))/3.d0)
            END IF
        END FUNCTION cbrt

        FUNCTION to_rad(x) !Convert to radians
            REAL :: x, to_rad
            to_rad = x * RAD
        END FUNCTION to_rad

        FUNCTION to_deg(x) !Convert to degrees
            REAL :: x, to_deg
            to_deg = x * (PI/180.d0)
        END FUNCTION to_deg

        !Trig functions altered for degrees
        FUNCTION SIND(x)
            REAL :: x, SIND
            SIND = SIN(x * (1.d0/RAD))
        END FUNCTION SIND

        FUNCTION COSD(x)
            REAL :: x, COSD
            COSD = COS(x * (1.d0/RAD))
       END FUNCTION COSD

       FUNCTION ATAND(x)
           REAL :: x, ATAND
           ATAND = ATAN(x) * RAD
       END FUNCTION ATAND

       FUNCTION ATAN2D(y, x)
           REAL :: y, x, ATAN2D
           ATAN2D = RAD * ATAN2(y, x)
       END FUNCTION ATAN2D
END MODULE degfun
