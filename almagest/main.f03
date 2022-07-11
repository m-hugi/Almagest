program main
  use almagest

  implicit none

  integer  :: hh, mm, ss, M, d, Y, utc_offset  
  real(8) :: alt, amz, RA, Decl

  print *, "Enter lat/lon: "
  read  *, lat, lon

  print *, "Enter Month/Day/Year: "
  read  *, M, d, Y

  print *, "Enter hour/minute/second: "
  read  *, hh, mm, ss

  print *, "Enter UTC offset: "
  read  *, utc_offset

  call computeSolarTimeValues(lat, lon, M, d, Y, hh, mm, ss, utc_offset)

  call computeSolarPosition(alt, amz)
  call computeLunarPosition(RA, Decl)

  print *, "Solar Alt: ", alt
  print *, "Solar Azm: ", amz

  print *, "Lunar RA:   ", RA
  print *, "Lunar Decl: ", Decl
end program main