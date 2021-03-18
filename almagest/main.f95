PROGRAM main
  USE almagest

  IMPLICIT NONE

  INTEGER  :: hh, mm, ss, M, d, Y, utc_offset  
  REAL(8)  :: lat, lon, sma, sml, ooe, LST, JDN

  REAL(8), DIMENSION(2) :: solp, lunp

  PRINT *, "Enter lat/lon: "
  READ  *, lat, lon

  PRINT *, "Enter Month/Day/Year: "
  READ  *, M, d, Y

  PRINT *, "Enter hour/minute/second: "
  READ  *, hh, mm, ss

  PRINT *, "Enter UTC offset: "
  READ  *, utc_offset

  call computeSolarTimeValues(M, d, Y, hh, mm, ss, utc_offset, lon, sma, sml, ooe, LST, JDN)

  solp = computeSolarPosition(lat, sma, ooe, LST, JDN)
  lunp = computeLunarPosition(lat, sma, sml, ooe, LST, JDN)

  PRINT *, "Solar Alt: ", solp(1)
  PRINT *, "Solar Azm: ", solp(2)

  PRINT *, "Lunar RA:   ", lunp(1)
  PRINT *, "Lunar Decl: ", lunp(2)
END PROGRAM main
