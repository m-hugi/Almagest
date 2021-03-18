!TODO: Might need to normalize all angles
!Sufficiently test for accuracy

MODULE almagest
  IMPLICIT NONE

  REAL*8, PARAMETER :: PI  = ATAN(1.d0) * 4.d0
  REAL*8, PARAMETER :: RAD = (180.d0/PI)

CONTAINS
  !Valid for all Gregorian calendar dates after November 23, 4713 B.C.
  SUBROUTINE computeSolarTimeValues(m, d, y, hh, mm, ss, utc_offset, lon, sma, sml, ooe, LST, JDN)
    IMPLICIT NONE

    REAL*8,  INTENT(IN)  :: lon
    INTEGER, INTENT(IN)  :: m, d, y, hh, mm, ss, utc_offset
    REAL*8,  INTENT(OUT) :: sma, sml, ooe, LST, JDN

    INTEGER :: UT
    REAL*8  :: lo, ma

    !Initial time values
    UT = (hh-utc_offset) + (mm/60) + (ss/3600)
    JDN = REAL((367*Y - 7*(Y+(M+9)/12)/4 - 3*((Y+(M-9)/7)/100 + 1)/4 + 275*M/9 + D-730515) + UT/24)

    !Solar values
    lo  = 282.9404d0 + 4.70935E-5     * JDN !Longitude of the perihelion
    ma  = 356.0470d0 + 0.9856002585d0 * JDN !Mean anomaly
    ooe = 23.4393d0  - 3.563E-7       * JDN !Obliquity of the Ecliptic

    sma = normAng(ma)     !Normalized mean anomaly; might need to be normalized everywhere
    sml = normAng(lo + ma) !Solar mean longitude

    !Sidereal Time at Greenwich + UT + lon = Local SideREAL*8 Time
    LST = (sml/15.d0 + 12.d0) + REAL(UT) + (lon/15.d0)

    !Normalize time value
    IF (LST .LT. 0.d0) THEN
       LST =  LST + 24.d0
    ELSE IF (LST .GT. 24.d0) THEN
       LST = LST - 24.d0
    END IF
  END SUBROUTINE computeSolarTimeValues

  FUNCTION computeSolarPosition(lat, sma, ooe, LST, JDN) RESULT(altazm) !=> [Azm, Alt]
    REAL*8 :: LST, JDN
    REAL*8 :: lat, sma, ooe 

    REAL*8, DIMENSION(2) :: altazm

    REAL*8 :: w, e, ea, L
    REAL*8 :: x, xe, y, ye, z, ze, r, v
    REAL*8 :: HA, RA, Decl, alt, azm

    w = 282.9404d0 + 4.70935E-5     * JDN !Longitude of perihelion)
    e = 0.016709d0 - 1.151E-9       * JDN !Eccentricity)

    ea = calculateEA0(e, sma)   !Calculate eccentric anomaly

    !Sun's coordinates in ecliptic plane
    x = COSD(ea) - e
    y = SIND(ea) * SQRT(1.d0 - e**2)

    r = SQRT(x**2 + y**2)
    v = ATAN2D(y, x)

    !Longitude of the Sun
    L = normAng(w + v)

    !Sun's ecliptic rectuangular coordinates
    x = r * COSD(L)
    y = r * SIND(L)
    z = 0.d0

    !Rotate coordinates
    xe = x
    ye = y * COSD(ooe) - z * SIND(ooe)
    ze = y * SIND(ooe) + z * COSD(ooe)

    !Convert to Geocentric RA/Decl/Dist
    RA   = ATAN2D(ye, xe)
    Decl = ATAN2D(ze, SQRT(xe**2 + ye**2))
    r    = SQRT(xe**2 + ye**2 + ze**2)

    HA = (LST - (RA/15.d0)) * 15    !Hour angle

    !Convert HA and Decl to rectangular coordinates
    x = COSD(HA) * COSD(Decl) ! -> South celestial equator
    y = SIND(HA) * COSD(Decl) ! -> Western horizon
    z = SIND(Decl)            ! -> North celestial pole

    !Rotate coordinates along an East-West axis
    xe = x * SIND(lat) - z * COSD(lat) ! -> North-South
    ye = y                             ! -> East-West
    ze = x * COSD(lat) + z * SIND(lat) ! -> Zenith

    !Convert to Topocentric azimuth/altitude
    alt = ASIND(ze)
    azm = ATAN2D(ye, xe) + 180.d0

    !Account for Atmospheric reftaction
    IF ((alt .LT. 85.d0) .AND. (alt .GT. 5.d0)) THEN
       alt = alt + 1.d0/3600.d0 * ((58.1d0/TAND(alt)) - (0.07d0/(TAND(alt)**3)) + (0.000086d0/(TAND(alt)**5)))
    ELSE IF ((alt .LT. 5.d0) .AND. (alt .GT. -0.575d0)) THEN
       alt = alt + 1.d0/3600.d0 * (1735.d0 - 518.2d0*TAND(alt) + 12.79d0*(alt**3) + 0.711d0*(alt**4))
    ELSE IF (alt .LT. -0.575d0) THEN
       alt = alt + 1.d0/3600.d0 * (-20.774d0/TAND(alt))
    END IF

    altazm = [alt, azm]
  END FUNCTION computeSolarPosition

  FUNCTION computeLunarPosition(olat, sma, sml, ooe, LST, JDN) RESULT(radecl) !=> [RA, Decl]
    REAL*8 :: LST, JDN
    REAL*8 :: olat, sma, sml, ooe 

    REAL*8, DIMENSION(2) :: radecl

    REAL*8 :: n, i, w, a, e, m, d, l, f, ea
    REAL*8 :: RA, Decl, glat, HA, rho, g, lp, lat, lon 
    REAL*8 :: x, xe, xeq, y, ye, yeq, ze, zeq, r, v


    !Calculate lunar orbital elements
    i = 5.1454d0   !Inclination
    e = 0.0549d0   !Eccentricity
    a = 60.2666d0  !Mean distance

    n = normAng(125.1228d0 - 0.0529538083d0  * JDN)  !Longitude of the ascending node
    w = normAng(318.0634d0 + 0.1643573223d0  * JDN)  !Argument of the perigee
    m = normAng(115.3654d0 + 13.0649929509d0 * JDN)  !Mean anomaly

    l = n + w + m !Mean longitude
    d = l - sml   !Mean elongation
    f = l - n     !Argument of Latitude

    !Compute eccentric anomaly
    ea = calculateEA(e, m, 0.005d0)

    !Compute rectanglar coordinates
    x = a * (COSD(ea) - e)
    y = a * SQRT(1 - e**2) * SIND(ea)

    !Convert to polar coordinatesx
    r = SQRT(x**2 + y**2)  !Distance
    v = ATAN2D(y, x)       !True anomaly

    !Compute ecliptic coordinates
    xe = r * (COSD(n)   * COSD(v+w) - SIND(n) * SIND(v+w) * COSD(i))
    ye = r * (SIND(n)   * COSD(v+w) + COSD(n) * SIND(v+w) * COSD(i))
    ze = r * SIND(v+w) * SIND(i)

    !Convert to lat, lon, dist
    lon = normAng(ATAN2D(ye, xe))
    lat = ATAN2D(ze, SQRT(xe**2 + ye**2))
    r   = SQRT(xe**2 + ye**2 + ze**2)

    !Account for lunar perturbations (Might need to be SIND)
    lon = lon -1.274d0 * SIND(m - 2*d)       & !(Evection)
              +0.658d0 * SIND(2*d)           & !(Variation)
              -0.186d0 * SIND(sma)           & !(Yearly equation)
              -0.059d0 * SIND(2*m - 2*d)     &
              -0.057d0 * SIND(m - 2*d + sma) &
              +0.053d0 * SIND(m + 2*d)       &
              +0.046d0 * SIND(2*d - sma)     &
              +0.041d0 * SIND(m - sma)       &
              -0.035d0 * SIND(d)             & !(Parallactic equation)
              -0.031d0 * SIND(m + sma)       &
              -0.015d0 * SIND(2*F - 2*D)     &
              +0.011d0 * SIND(m - 4*D)      

    lat = lat -0.173d0 * SIND(f - 2*d)       &
              -0.055d0 * SIND(m - F - 2*d)   &
              -0.046d0 * SIND(m + F - 2*d)   &
              +0.033d0 * SIND(d + 2*d)       &
              +0.017d0 * SIND(2*m + f)       

    r   =   r -0.58d0 * COSD(m - 2*d)        &
              -0.46d0 * COSD(2*d)            

    !Compute rectangular ecliptic coordinates
    xe = r * COSD(lon) * COSD(lat)
    ye = r * SIND(lon) * COSD(lat)
    ze = r * SIND(lat)

    !Rotate coordinates to equatorial
    xeq = xe
    yeq = ye * COSD(ooe) - ze * SIND(ooe)
    zeq = ye * SIND(ooe) + ze * COSD(ooe)

    !Compute geocentric Right Acension and Declination
    RA = normAng(ATAN2D(yeq, xeq))
    Decl = ATAN2D(zeq, SQRT(xeq**2 + yeq**2))
    

    !Transalte RA, Decl into topocentic values
    lp = ASIND(1.d0/r) !Lunar parallax

    !Account for flattening of the Earth
    glat = olat - 0.1924d0 * SIND(2*olat)        !Geocentric latitude
    rho  = 0.99833d0 + 0.00167d0 * COSD(2*olat) !Distance from the center of the Earth

    HA = normAng((LST * 15) - RA) !Hour Angle (LST * 15 converts to degrees)
    g  = ATAND(TAND(glat) / COSD(HA)) ! Auxillary angle

    !Convert geocentric RA,Decl to topocentric RA,Decl
    RA   = RA   - lp * rho * COSD(glat) * SIND(HA) / COSD(Decl)
    Decl = Decl - lp * rho * SIND(glat) * SIND(g - Decl) / SIND(g)

    radecl = [RA, Decl]
  END FUNCTION computeLunarPosition

  FUNCTION calculateEA0(e, M) !Calculate inital eccentric anomaly
    REAL*8 :: e, M, calculateEA0

    calculateEA0 = M + to_rad(e) * SIND(M) * (1.0 + e * COSD(M))
  END FUNCTION calculateEA0

  FUNCTION calculateEA(e, M, accu) !Iteratively calculate a more accurate eccentric anomaly
    REAL*8  :: e, M, EA0, EA1, accu, calculateEA

    EA0 = calculateEA0(e, M)
    DO WHILE (.TRUE.)
       EA1 = EA0 - (EA0 - to_rad(e) * SIND(EA0) - M) / (1 - e * COSD(EA0))
       IF (ABS(EA1 - EA0) .LE. accu) THEN
          EXIT
       ELSE
          EA0 = EA1
       END IF
    END DO

    calculateEA = EA0
  END FUNCTION calculateEA

  FUNCTION fromHours(h, m, s)
    INTEGER :: fromHours
    INTEGER  :: h, m, s

    fromHours = h + (m/60) + (s/3600)
  END FUNCTION fromHours

  !General math functions
  FUNCTION normAng(x) !Keeps degree value within proper range
    REAL*8 :: x, normAng
    
    normAng = x - FLOOR(x/360.d0)*360.d0
  END FUNCTION normAng

  FUNCTION to_rad(x) !Convert to radians
    REAL*8 :: x, to_rad

    to_rad = x * RAD
  END FUNCTION to_rad

  FUNCTION to_deg(x) !Convert to degrees
    REAL*8 :: x, to_deg

    to_deg = x * (PI/180.d0)
  END FUNCTION to_deg
END MODULE almagest