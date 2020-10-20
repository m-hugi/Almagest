PROGRAM almagest
    USE degfun

    IMPLICIT NONE

    REAL :: UT, t
    REAL :: lat, lon
    INTEGER  :: hh, mm, ss, M, d, Y, UTC_offset

    PRINT *, "Enter lat/lon: "
    READ  *, lat, lon

    PRINT *, "Enter Month/Day/Year: "
    READ  *, M, d, Y

    PRINT *, "Enter hour/minute/second: "
    READ  *, hh, mm, ss

    PRINT *, "Enter UTC offset: "
    READ  *, UTC_offset

    UT = from_hours(hh-UTC_offset, mm, ss)
    t  = getJDN(Y, M, d, UT)

    call computeSolarPosition(lat, lon, t, UT)

    CONTAINS
        SUBROUTINE computeSolarPosition(lat, lon, t, UT)
            REAL, INTENT(IN) :: lat, lon, t, UT
            
            REAL :: ST, GMST0, HA
            REAL :: EA, EL, RA, Decl
            REAL :: x, xe, y, ye, z, ze, r, v
            REAL :: w, a, e, M, L, ML, o, alt, azm
            
            w = 282.9404d0 + 4.70935E-5     * t    !(longitude of perihelion)
            a = 1.d0                               !(mean distance, a.u.)
	    e = 0.016709d0 - 1.151E-9       * t    !(eccentricity)
            M = 356.0470d0 + 0.9856002585d0 * t    !(mean anomaly)
            o = 23.4393d0 - 3.563E-7        * t

            M  = normalize_ang(M)     !Normalize mean anomaly
            EA = calculateEA0(e, M)   !Calculate eccentric anomaly
            ML = normalize_ang(w + M) !Sun's mean longitude

            !Sun's coordinates in ecliptic plane
            x = COSD(EA) - e
            y = SIND(EA) * SQRT(1.d0 - e**2)

            r = SQRT(x**2 + y**2)
            v = ATAN2D(y, x)

            !Longitude of the Sun
            L = normalize_ang(w + v)

            !Sun's ecliptic rectuangular coordinates
            x = r * COSD(L)
            y = r * SIND(L)
            z = 0.d0

            !Rotate coordinates
            xe = x
            ye = y * COSD(o) - z * SIND(o)
            ze = y * SIND(o) + z * COSD(o)

            !Convert to Geocentric RA/Decl/Dist
            RA   = ATAN2D(ye, xe)
            Decl = ATAN2D(ze, SQRT(xe**2 + ye**2))
            r    = SQRT(xe**2 + ye**2 + ze**2)

            GMST0 = ML/15.d0 + 12.d0       !Sidereal Time at Greenwich
            ST    = GMST0 + UT + lon/15.d0 !Local Sidereal Time

            IF (ST .LT. 0.d0) THEN
                ST =  ST + 24.d0
            ELSE IF (ST .GT. 24.d0) THEN
                ST = ST - 24.d0
            END IF

            HA = (ST - (RA/15.d0)) * 15    !Hour angle

            !Convert HA and Decl to rectangular coordinates
            x = COSD(HA) * COSD(Decl) ! -> South celestial equator
            y = SIND(HA) * COSD(Decl) ! -> Western horizon
            z = SIND(Decl)            ! -> North celestial pole

            !Rotate coordinates along an East-West axis
            xe = x * SIND(lat) - z * COSD(lat) ! -> North-South
            ye = y                             ! -> East-West
            ze = x * COSD(lat) + z * SIND(lat) ! -> Zenith

            !Convert to Topocentric azimuth/altitude
            azm = ATAN2D(ye, xe) + 180.d0
            alt = ASIND(ze)

            !Account for Atmospheric reftaction
            IF ((alt .LT. 85.d0) .AND. (alt .GT. 5.d0)) THEN
                alt = alt + 1.d0/3600.d0 * ((58.1d0/TAND(alt)) - (0.07d0/(TAND(alt)**3)) + (0.000086d0/(TAND(alt)**5)))
            ELSE IF ((alt .LT. 5.d0) .AND. (alt .GT. -0.575d0)) THEN
                alt = alt + 1.d0/3600.d0 * (1735.d0 - 518.2d0*TAND(alt) + 12.79d0*(alt**3) + 0.711d0*(alt**4))
            ELSE IF (alt .LT. -0.575d0) THEN
                alt = alt + 1.d0/3600.d0 * (-20.774d0/TAND(alt))
            END IF

            PRINT *, "Azimuth: ", azm
            PRINT *, "Altitude: ", alt
        END SUBROUTINE computeSolarPosition

        FUNCTION calculateEA0(e, M)
            REAL :: e, M, calculateEA0

            calculateEA0 = M + to_rad(e) * SIND(M) * (1.0 + e * COSD(M))
        END FUNCTION calculateEA0

        FUNCTION calculateEA(e, M, accu)
            REAL :: e, M, accu, EA0, EA1, calculateEA

            EA0 = calculateEA0(e, M)
            DO WHILE (.TRUE.)
                EA1 = EA0 - (EA0 - to_deg(e) * SIND(EA0)-M) / (1 - e * COSD(EA0))
                IF (ABS(EA1 - EA0) .LE. accu) THEN
                    STOP
                ELSE
                    EA0 = EA1
                END IF
            END DO
            calculateEA = EA1
        END FUNCTION calculateEA

        !Valid for all Gregorian calendar dates after November 23, 4713 B.C.
        FUNCTION getJDN(Y, M, D, UT)
            REAL :: getJDN, UT
            INTEGER  :: Y, M, D, JD

            JD = 367*Y - 7*(Y+(M+9.)/12)/4 - 3*((Y+(M-9)/7)/100 + 1)/4 + 275*M/9 + D-730515
            getJDN = JD + UT/24.d0
        END FUNCTION getJDN
END PROGRAM almagest
