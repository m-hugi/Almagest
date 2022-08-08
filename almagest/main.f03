program main
  use almagest

  implicit none

  integer  :: hh, mm, ss, M, d, Y, utc_offset  

  real*8, dimension(7, 2) :: radecls
  real*8, dimension(2) :: solar_pos, lunar_pos

  !print *, "Enter lat/lon: "
  !read  *, lat, lon

  !print *, "Enter Month/Day/Year: "
  !read  *, M, d, Y

  !print *, "Enter hour/minute/second: "
  !read  *, hh, mm, ss

  !print *, "Enter UTC offset: "
  !read  *, utc_offset

  ! Test Values
  lat = 60
  lon = 15
  M = 4
  d = 19
  Y = 1990
  hh = 0
  mm = 0
  ss = 00
  utc_offset = 0

  call computeSolarTimeValues(lat, lon, M, d, Y, hh, mm, ss, utc_offset)
  
  solar_pos = computeSolarPosition()
  lunar_pos = computeLunarPosition()
  radecls = computePlanetaryPositions()

  print *, "Solar Topo Alt: ", solar_pos(1)
  print *, "Solar Topo Azm: ", solar_pos(2)
  print "(A)"
  print *, "Lunar Topo RA:   ", toRAHours(lunar_pos(1))
  print *, "Lunar Topo Decl: ", toArcTime(lunar_pos(2))
  print "(A)"
  print *, "Mercury Topo RA: ", toRAHours(radecls(1, 1))
  print *, "Mercury Topo Decl: ", toArcTime(radecls(1, 2))
  print "(A)"
  print *, "Venus Topo RA: ", toRAHours(radecls(2, 1))
  print *, "Venus Topo Decl: ", toArcTime(radecls(2, 2))
  print "(A)"
  print *, "Mars Topo RA: ", toRAHours(radecls(3, 1))
  print *, "Mars Topo Decl: ", toArcTime(radecls(3, 2))
  print "(A)"
  print *, "Jupiter Topo RA: ", toRAHours(radecls(4, 1))
  print *, "Jupiter Topo Decl: ", toArcTime(radecls(4, 2))
  print "(A)"
  print *, "Saturn Topo RA: ", toRAHours(radecls(5, 1))
  print *, "Saturn Topo Decl: ", toArcTime(radecls(5, 2))
  print "(A)"
  print *, "Uranus Topo RA: ", toRAHours(radecls(6, 1))
  print *, "Uranus Topo Decl: ", toArcTime(radecls(6, 2))
  print "(A)"
  print *, "Neptune Topo RA: ", toRAHours(radecls(6, 1))
  print *, "Neptune Topo Decl: ", toArcTime(radecls(6, 2))
end program main