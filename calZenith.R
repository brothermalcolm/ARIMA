library(insol)


# function to calculate sun position and extraterrestial irradiance
calZen <- function(Tm, lat, lon, tz, tilt=0, orientation=0) {
  # tilt angle: between 0 and 90
  # orientation: between 0 and 360, north 0, positive in clockwise
  
  jd = JD(Tm) # Computes Julian Day from dates as POSIXct object. 
  sunv = sunvector(jd, lat, lon, tz) # Calculates a unit vector in the direction of the sun from the observer position. 
  # sunpos() Returns a matrix of azimuth and zenith angles of the sun given the unit vectors from the observer to the direction of the sun.
  azi = round(sunpos(sunv)[,1],3) # azimuth of the sun
  zen = round(sunpos(sunv)[,2],3) # zenith angle
  surface.norm = normalvector(tilt, orientation)
  inc = round(as.numeric(degrees(acos(sunv%*% as.vector(surface.norm)))),3)
  dec = declination(jd)*pi/180
  re = 1.000110+0.034221*cos(dec)+0.001280*sin(dec)+0.00719*cos(2*dec)+0.000077*sin(2*dec)
  Io = round(1362*re,3) # extraterrestrial direct normal irradiance
  Ioh = round(1362*re*cos(zen*pi/180)) # extraterrestrial horizontal irradiance
  Ic = round(0.8298*1362*re*(cos(zen*pi/180))^1.3585*exp(-0.00135*(90-zen)), 1)
  Ioh = ifelse(zen>=90, 0, Ioh)
  Ic = ifelse(zen>=90, 0, Ic)
  out = list(zen, azi, inc, Io, Ioh, Ic)
  names(out) = c("zenith", "azimuth", "incidence angle", "Io", "Ioh", "Ic")
  out
}