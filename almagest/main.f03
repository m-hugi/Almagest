program main
  use almagest

  implicit none

  integer  :: hh, mm, ss, M, d, Y, utc_offset  
  real(8) :: alt, amz, RA, Decl
  real*8, dimension(7, 2) :: radecls

  print *, "Enter lat/lon: "
  read  *, lat, lon

  print *, "Enter Month/Day/Year: "
  read  *, M, d, Y

  print *, "Enter hour/minute/second: "
  read  *, hh, mm, ssd

  print *, "Enter UTC offset: "
  read  *, utc_offset

  ! Test Values
  ! lat = 60
  ! lon = 15
  ! M = 4
  ! d = 19
  ! Y = 1990
  ! hh = 0
  ! mm = 0
  ! ss = 00
  ! utc_offset = 2

  call computeSolarTimeValues(lat, lon, M, d, Y, hh, mm, ss, utc_offset)
  call computeSolarPosition(alt, amz)
  call computePlanetaryPositions(radecls)
  call computeLunarPosition(RA, Decl)

  print *, "Solar Topo Alt: ", alt
  print *, "Solar Topo Azm: ", amz
  print "(A)"
  print *, "Lunar Topo RA:   ", RA
  print *, "Lunar Topo Decl: ", Decl
  print "(A)"
  print *, "Mercury Topo RA: ", radecls(1, 1)
  print *, "Mercury Topo Decl: ", radecls(1, 2)
  print "(A)"
  print *, "Venus Topo RA: ", radecls(2, 1)
  print *, "Venus Topo Decl: ", radecls(2, 2)
  print "(A)"
  print *, "Mars Topo RA: ", radecls(3, 1)
  print *, "Mars Topo Decl: ", radecls(3, 2)
  print "(A)"
  print *, "Jupiter Topo RA: ", radecls(4, 1)
  print *, "Jupiter Topo Decl: ", radecls(4, 2)
  print "(A)"
  print *, "Saturn Topo RA: ", radecls(5, 1)
  print *, "Saturn Topo Decl: ", radecls(5, 2)
  print "(A)"
  print *, "Uranus Topo RA: ", radecls(6, 1)
  print *, "Uranus Topo Decl: ", radecls(6, 2)
  print "(A)"
  print *, "Neptune Topo RA: ", radecls(6, 1)
  print *, "Neptune Topo Decl: ", radecls(6, 2)
end program main
