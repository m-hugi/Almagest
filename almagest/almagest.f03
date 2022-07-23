!TODO:
! Might need to normalize all angles
! Use computeHeliocentricCoords everywhere
! Sufficiently test for accuracy

module almagest
  implicit none

  real*8, parameter :: PI  = ATAN(1.d0) * 4.d0
  real*8, parameter :: RAD = (180.d0/PI)

  real*8 :: lat, lon
  real*8 :: LST, JDN
  real*8 :: olat, sma, sml, ooe 
  real*8 :: sol_x, sol_y, sol_z

  private planet
  type planet
    real*8 :: n ! Longitude of the ascending node
    real*8 :: i ! Inclination
    real*8 :: w ! Argument of perihelion
    real*8 :: a ! Semi-major axis
    real*8 :: e ! Eccentricity
    real*8 :: m ! Mean anomaly
    real*8 :: lon_adj, lat_adj ! Perturbations
  end type
contains
  !Valid for all Gregorian calendar dates after November 23, 4713 B.C.
  subroutine computeSolarTimeValues(latitiude, longitude, m, d, y, hh, mm, ss, utc_offset)
    implicit none

    real*8,  intent(in)  :: latitiude, longitude 
    integer, intent(in)  :: m, d, y, hh, mm, ss, utc_offset

    integer :: UT
    real*8  :: lo, ma

    lat = latitiude
    lon = longitude

    !Initial time values
    UT = (hh-utc_offset) + (mm/60) + (ss/3600)
    JDN = real((367*Y - 7*(Y+(M+9)/12)/4 - 3*((Y+(M-9)/7)/100 + 1)/4 + 275*M/9 + D-730515) + UT/24)

    !Solar values
    lo  = 282.9404d0 + 4.70935E-5     * JDN !Longitude of the perihelion
    ma  = 356.0470d0 + 0.9856002585d0 * JDN !Mean anomaly
    ooe = 23.4393d0  - 3.563E-7       * JDN !Obliquity of the Ecliptic

    sma = normAng(ma)     !Normalized mean anomaly; might need to be normalized everywhere
    sml = normAng(lo + ma) !Solar mean longitude

    !Sidereal Time at Greenwich + UT + lon = Local Sidereal Time
    LST = (sml/15.d0 + 12.d0) + real(UT) + (lon/15.d0)

    !Normalize time value
    if (LST .lt. 0.d0) then
       LST =  LST + 24.d0
    else if (LST .gt. 24.d0) then
       LST = LST - 24.d0
    end if
  end subroutine computeSolarTimeValues

  subroutine computeSolarPosition(azm, alt)
    implicit none

    real*8, intent(out) :: alt, azm

    real*8 :: w, e, ea, L
    real*8 :: x, xe, y, ye, z, ze, r, v
    real*8 :: HA, RA, Decl

    w = 282.9404d0 + 4.70935E-5     * JDN !Longitude of perihelion)
    e = 0.016709d0 - 1.151E-9       * JDN !Eccentricity)

    ea = calculateEA0(e, sma)   !Calculate eccentric anomaly

    !Sun's coordinates in ecliptic plane
    x = cosd(ea) - e
    y = sind(ea) * sqrt(1.d0 - e**2)

    r = sqrt(x**2 + y**2)
    v = atan2d(y, x)

    !Longitude of the Sun
    L = normAng(w + v)

    !Sun's ecliptic rectuangular coordinates`
    sol_x = r * cosd(L)
    sol_y = r * sind(L)
    sol_z = 0.d0

    !Rotate coordinates
    xe = sol_x
    ye = sol_y * cosd(ooe) - sol_z * sind(ooe)
    ze = sol_y * sind(ooe) + sol_z * cosd(ooe)

    !Convert to Geocentric RA/Decl/Dist
    RA   = atan2d(ye, xe)
    Decl = atan2d(ze, sqrt(xe**2 + ye**2))
    r    = sqrt(xe**2 + ye**2 + ze**2)

    HA = (LST - (RA/15.d0)) * 15    !Hour angle

    !Convert HA and Decl to rectangular coordinates
    x = cosd(HA) * cosd(Decl) ! -> South celestial equator
    y = sind(HA) * cosd(Decl) ! -> Western horizon
    z = sind(Decl)            ! -> North celestial pole

    !Rotate coordinates along an East-West axis
    xe = x * sind(lat) - z * cosd(lat) ! -> North-South
    ye = y                             ! -> East-West
    ze = x * cosd(lat) + z * sind(lat) ! -> Zenith

    !Convert to Topocentric azimuth/altitude
    alt = asind(ze)
    azm = atan2d(ye, xe) + 180.d0

    !Account for Atmospheric reftaction
    if ((alt .lt. 85.d0) .and. (alt .gt. 5.d0)) then
       alt = alt + 1.d0/3600.d0 * ((58.1d0/tand(alt)) - (0.07d0/(tand(alt)**3)) + (0.000086d0/(tand(alt)**5)))
    else if ((alt .lt. 5.d0) .and. (alt .gt. -0.575d0)) then
       alt = alt + 1.d0/3600.d0 * (1735.d0 - 518.2d0*tand(alt) + 12.79d0*(alt**3) + 0.711d0*(alt**4))
    else if (alt .lt. -0.575d0) then
       alt = alt + 1.d0/3600.d0 * (-20.774d0/tand(alt))
    end if
  end subroutine computeSolarPosition

  subroutine computeLunarPosition(RA, Decl)
    implicit none

    real*8, intent(out) :: RA, Decl

    real*8 :: n, i, w, a, e, m, d, l, f, ea
    real*8 :: mlat, mlon, glat, HA, rho, g, lp
    real*8 :: x, xe, xeq, y, ye, yeq, ze, zeq, r, v


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
    x = a * (cosd(ea) - e)
    y = a * sqrt(1 - e**2) * sind(ea)

    !Convert to polar coordinatesx
    r = sqrt(x**2 + y**2)  !Distance
    v = atan2d(y, x)       !True anomaly

    !Compute ecliptic coordinates
    xe = r * (cosd(n)   * cosd(v+w) - sind(n) * sind(v+w) * cosd(i))
    ye = r * (sind(n)   * cosd(v+w) + cosd(n) * sind(v+w) * cosd(i))
    ze = r * sind(v+w) * sind(i)

    !Convert to lat, lon, dist
    mlon = normAng(atan2d(ye, xe))
    mlat = atan2d(ze, sqrt(xe**2 + ye**2))
    r   = sqrt(xe**2 + ye**2 + ze**2)

    !Account for lunar perturbations (Might need to be sind)
    mlon = mlon -1.274d0 * sind(m - 2*d)       & !(Evection)
                +0.658d0 * sind(2*d)           & !(Variation)
                -0.186d0 * sind(sma)           & !(Yearly equation)
                -0.059d0 * sind(2*m - 2*d)     &
                -0.057d0 * sind(m - 2*d + sma) &
                +0.053d0 * sind(m + 2*d)       &
                +0.046d0 * sind(2*d - sma)     &
                +0.041d0 * sind(m - sma)       &
                -0.035d0 * sind(d)             & !(Parallactic equation)
                -0.031d0 * sind(m + sma)       &
                -0.015d0 * sind(2*F - 2*D)     &
                +0.011d0 * sind(m - 4*D)      

    mlat = mlat -0.173d0 * sind(f - 2*d)       &
                -0.055d0 * sind(m - F - 2*d)   &
                -0.046d0 * sind(m + F - 2*d)   &
                +0.033d0 * sind(d + 2*d)       &
                +0.017d0 * sind(2*m + f)       

    r    =  r -0.58d0 * cosd(m - 2*d)        &
              -0.46d0 * cosd(2*d)            

    !Compute rectangular ecliptic coordinates
    xe = r * cosd(mlon) * cosd(mlat)
    ye = r * sind(mlon) * cosd(mlat)
    ze = r * sind(mlat)

    !Rotate coordinates to equatorial
    xeq = xe
    yeq = ye * cosd(ooe) - ze * sind(ooe)
    zeq = ye * sind(ooe) + ze * cosd(ooe)

    !Compute geocentric Right Acension and Declination
    RA = normAng(atan2d(yeq, xeq))
    Decl = atan2d(zeq, sqrt(xeq**2 + yeq**2))

    !Transalte RA, Decl into topocentic values
    lp = Asind(1.d0/r) !Lunar parallax

    !Account for flattening of the Earth
    glat = lat - 0.1924d0 * sind(2*lat)        !Geocentric latitude
    rho  = 0.99833d0 + 0.00167d0 * cosd(2*lat) !Distance from the center of the Earth

    HA = normAng((LST * 15) - RA) !Hour Angle (LST * 15 converts to degrees)
    g  = atand(tand(glat) / cosd(HA)) ! Auxillary angle

    !Convert geocentric RA,Decl to topocentric RA,Decl
    RA   = RA   - lp * rho * cosd(glat) * sind(HA) / cosd(Decl)
    Decl = Decl - lp * rho * sind(glat) * sind(g - Decl) / sind(g)
  end subroutine computeLunarPosition

  subroutine computePlanetaryPositions(radecls)
    implicit none

    integer :: i

    real*8 :: plat, plon, pr, HA, glat, g, rho, pp
    real*8 :: x, y, z, xeq, yeq, zeq, RA, Decl, r

    real*8, dimension(7,2), INTENT(OUT) :: radecls

    type(planet) :: p
    type(planet) :: Mercury, Venus, Mars, Jupiter, Saturn, Uranus, Neptune
    type(planet), dimension(7) :: planets

    Mercury%N =  48.3313d0 + 3.24587E-5   * JDN 
    Mercury%i =   7.0047d0 + 5.00E-8      * JDN 
    Mercury%w =  29.1241d0 + 1.01444E-5   * JDN 
    Mercury%a = 0.387098d0                            
    Mercury%e = 0.205635d0     + 5.59E-10         * JDN 
    Mercury%M = 168.6562d0 + 4.0923344368d0 * JDN

    Venus%N =  76.6799d0 + 2.46590E-5   * JDN
    Venus%i =   3.3946d0 + 2.75E-8      * JDN
    Venus%w =  54.8910d0 + 1.38374E-5   * JDN
    Venus%a = 0.723330d0
    Venus%e = 0.006773d0     - 1.302E-9         * JDN
    Venus%M =  48.0052d0 + 1.6021302244d0 * JDN

    Mars%N =  49.5574d0 + 2.11081E-5   * JDN
    Mars%i =   1.8497d0 - 1.78E-8      * JDN
    Mars%w = 286.5016d0 + 2.92961E-5   * JDN
    Mars%a = 1.523688d0
    Mars%e = 0.093405d0     + 2.516E-9         * JDN
    Mars%M =  18.6021d0 + 0.5240207766d0 * JDN

    Jupiter%N = 100.4542d0 + 2.76854E-5   * JDN
    Jupiter%i =   1.3030d0 - 1.557E-7     * JDN
    Jupiter%w = 273.8777d0 + 1.64505E-5   * JDN
    Jupiter%a = 5.20256d0
    Jupiter%e = 0.048498d0     + 4.469E-9  * JDN
    Jupiter%M = 19.8950d0 + 0.0830853001d0 * JDN

    Saturn%N = 113.6634d0 + 2.38980E-5   * JDN
    Saturn%i =   2.4886d0 - 1.081E-7     * JDN
    Saturn%w = 339.3939d0 + 2.97661E-5   * JDN
    Saturn%a = 9.55475d0
    Saturn%e = 0.055546d0     - 9.499E-9         * JDN
    Saturn%M = 316.9670d0 + 0.0334442282d0 * JDN

    Uranus%N =  74.0005d0 + 1.3978E-5    * JDN
    Uranus%i =   0.7733d0 + 1.9E-8       * JDN
    Uranus%w =  96.6612d0 + 3.0565E-5    * JDN
    Uranus%a = 19.18171d0 - 1.55E-8        * JDN
    Uranus%e = 0.047318d0 + 7.45E-9        * JDN
    Uranus%M = 142.5905d0 + 0.011725806d0  * JDN

    Neptune%N = 131.7806d0 + 3.0173E-5    * JDN
    Neptune%i =   1.7700d0 - 2.55E-7      * JDN
    Neptune%w = 272.8461d0 - 6.027E-6     * JDN
    Neptune%a = 30.05826d0 + 3.313E-8     * JDN
    Neptune%e = 0.008606d0 + 2.15E-9      * JDN
    Neptune%M = 260.2471d0 + 0.005995147d0  * JDN

    ! Calculate perturbations
    Jupiter%lon_adj = -0.332d0 * sind(2*Jupiter%M - 5*Saturn%M - 67.6d0) &
                      -0.056d0 * sind(2*Jupiter%M - 2*Saturn%M + 21)     &
                      +0.042d0 * sind(3*Jupiter%M - 5*Saturn%M + 21)     &
                      -0.036d0 * sind(Jupiter%M - 2*Saturn%M)            &
                      +0.022d0 * cosd(Jupiter%M - Saturn%M)              &
                      +0.023d0 * sind(2*Jupiter%M - 3*Saturn%M + 52)     &
                      -0.016d0 * sind(Jupiter%M - 5*Saturn%M - 69)       

    Saturn%lon_adj = +0.812d0 * sind(2*Jupiter%M - 5*Saturn%M - 67.6d0) &
                     -0.229d0 * cosd(2*Jupiter%M - 4*Saturn%M - 2d0)    &
                     +0.119d0 * sind(Jupiter%M - 2*Saturn%M - 3d0)      &
                     +0.046d0 * sind(2*Jupiter%M - 6*Saturn%M - 69d0)   &
                     +0.014d0 * sind(Jupiter%M - 3*Saturn%M + 32d0)     

    Saturn%lat_adj = -0.020d0 * cosd(2*Jupiter%M - 4*Saturn%M - 2d0) &
                     +0.018d0 * sind(2*Jupiter%M - 6*Saturn%M - 49d0)

    Uranus%lon_adj = +0.040d0 * sind(Saturn%M - 2*Uranus%M + 6d0)  &
                     +0.035d0 * sind(Saturn%M - 3*Uranus%M + 33d0) &
                     -0.015d0 * sind(Jupiter%M - Uranus%M + 20d0)

    planets = [Mercury, Mars, Venus, Jupiter, Saturn, Uranus, Neptune]

    do i = 1, size(planets)
      p = planets(i)
      call computeHeliocentricCoords(p%n, p%i, p%w, p%a, p%e, p%M, plat, plon, pr)

      ! Apply perturbation calculations
      plon = plon + p%lon_adj
      plat = plat + p%lat_adj

      ! Precess to epoch 2000
      plon = plon + 3.82394E-5 * (0-JDN)

      ! Conversion to ecliptic system
      x = pr * cosd(plon) * cosd(plat)
      y = pr * sind(plon) * cosd(plat)
      z = pr * sind(plat)

      ! Add solar rectangular coordinates to convert to geocentric
      x = x + sol_x
      y = y + sol_y
      z = z + sol_z

      ! Rotate to equatorial
      xeq = x
      yeq = y * cosd(ooe) - z * sind(ooe)
      zeq = y * sind(ooe) + z * cosd(ooe)

      !Compute geocentric Right Acension and Declination
      RA = normAng(atan2d(yeq, xeq))
      Decl = atan2d(zeq, sqrt(xeq**2 + yeq**2))
      r = sqrt(xeq**2 + yeq**2 + zeq**2)

      glat = lat - 0.1924d0 * sind(2*lat)        !Geocentric latitude
      rho  = 0.99833d0 + 0.00167d0 * cosd(2*lat) !Distance from the center of the Earth

      HA = normAng((LST * 15) - RA) !Hour Angle (LST * 15 converts to degrees)
      g  = atand(tand(glat) / cosd(HA)) ! Auxillary angle

      ! Calculate parllax angle
      pp = (8.794d0/3600) / r

      !Convert geocentric RA,Decl to topocentric RA,Decl
      RA   = RA   - pp * rho * cosd(glat) * sind(HA) / cosd(Decl)
      Decl = Decl - pp * rho * sind(glat) * sind(g - Decl) / sind(g)

      radecls(i, 1) = RA
      radecls(i, 2) = Decl
    end do
  end subroutine computePlanetaryPositions

  subroutine computeHeliocentricCoords(n, i, w, a, e, m, o_lat, o_lon, o_r)
    real*8, intent(in) :: n, i, w, a, e, m
    real*8, intent(out) :: o_lat, o_lon, o_r

    real*8 :: ea, x, y, xe, ye, ze, r, v

    ea = calculateEA(e, m, 0.005d0)

    !Compute rectanglar coordinates
    x = a * (cosd(ea) - e)
    y = a * sqrt(1 - e**2) * sind(ea)

    !Convert to polar coordinatesx
    r = sqrt(x**2 + y**2)  !Distance
    v = atan2d(y, x)       !True anomaly

    !Compute ecliptic coordinates
    xe = r * (cosd(n)   * cosd(v+w) - sind(n) * sind(v+w) * cosd(i))
    ye = r * (sind(n)   * cosd(v+w) + cosd(n) * sind(v+w) * cosd(i))
    ze = r * sind(v+w) * sind(i)

    !Convert to lat, lon, dist
    o_lon = normAng(atan2d(ye, xe))
    o_lat = atan2d(ze, sqrt(xe**2 + ye**2))
    o_r   = sqrt(xe**2 + ye**2 + ze**2)
  end subroutine computeHeliocentricCoords

  function calculateEA0(e, M) !Calculate inital eccentric anomaly
    real*8 :: e, M, calculateEA0

    calculateEA0 = M + to_rad(e) * sind(M) * (1.0 + e * cosd(M))
  end function calculateEA0

  function calculateEA(e, M, accu) !Iteratively calculate a more accurate EA
    real*8  :: e, M, EA0, EA1, accu, calculateEA

    EA0 = calculateEA0(e, M)
    DO WHILE (.TRUE.)
      EA1 = EA0 - (EA0 - to_rad(e) * sind(EA0) - M) / (1 - e * cosd(EA0))
      if (abs(EA1 - EA0) .le. accu) then
        exit
      else
        EA0 = EA1
      end if
    end do

    calculateEA = EA0
  end function calculateEA

  ! Currently unused
  !function fromHours(h, m, s)
  !  integer :: fromHours
  !  integer  :: h, m, s
  !
  !  fromHours = h + (m/60) + (s/3600)
  !end function fromHours

  !General math functions
  function normAng(x) !Keeps degree value within proper range
    real*8 :: x, normAng
    
    normAng = x - FLOOR(x/360.d0)*360.d0
  end function normAng

  function to_rad(x) !Convert to radians
    real*8 :: x, to_rad

    to_rad = x * RAD
  end function to_rad

  ! Currently unused
  !function to_deg(x) !Convert to degrees
  !  real*8 :: x, to_deg
  !
  !  to_deg = x * (PI/180.d0)
  !end function to_deg
end MODULE almagest