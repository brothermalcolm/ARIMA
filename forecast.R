###########################################################
# This file is run periodically every min (live forecast) #
###########################################################

#Clear all workspace
rm(list = ls(all = TRUE))

library(insol, quietly=TRUE)
library(forecast, quietly=TRUE)

source("/home/seris/ARIMA/db_util.R")
source("/home/seris/ARIMA/calZenith.R")

stationdata <- read.csv("/home/seris/ARIMA/stations.csv")
stationdata <- stationdata[,c("ID_VIRTUAL","LATITUDE","LONGITUDE")]
stationdata[nrow(stationdata)+1,] <- list("_AVG", 1.3005, 103.771)  # append the Lat and Lng for the center of Singapore


###########################################################################
# a general forecast function                                             #
# reading aggregated dataset from DB => detrend => forecast => predicted  #
###########################################################################
fcst <- function(station="_AVG", model="ARIMA", now=NULL, interval=60, nr.periods) {
  
  print(paste0("Forecasting for station", station, " with ", model, " model, ", interval, "min ahead"))
  # get the coord of the station for Solar Position calculations
  lat <- stationdata[stationdata$ID_VIRTUAL==station,]$LATITUDE
  lng <- stationdata[stationdata$ID_VIRTUAL==station,]$LONGITUDE
  tm.zone <- 8 # time zone 8
  
  tm.fmt <- "%Y-%m-%d %H:%M:%S"
  
  if(is.null(now)){
    now <- round(Sys.time(),"mins")
  }
  
  tm.begin <- strftime(now-60*interval*nr.periods, tm.fmt)
  tm.end <- strftime(now, tm.fmt)
  
  # a series of time stamps for clearsky index, corresponding to the middle of time binning intervals
  tm.ci.begin <- strftime(now-60*interval*(nr.periods-0.5), tm.fmt)
  tm.ci.seq <- seq(from=as.POSIXct(tm.ci.begin, tz="GMT"),
                by = paste0(interval," mins"), length.out=nr.periods)
  
  data.aggregated <- db.read.avg(paste0("SIN", station), par="AvgGSi00", tm.begin, tm.end, interval)
  data.aggregated$timekey <- nr.periods - data.aggregated$timekey # reverse the order of time key
  
  
  # === calc clear sky index ===
  sol.pos <- calZen(tm.ci.seq, lat, lng, tm.zone)
  ci <- data.aggregated$Avg / sol.pos$Ic[data.aggregated$timekey]
  
  angle.limit.Z <- 87 # limiting zentih angle
  select <- which(sol.pos$zenith[data.aggregated$timekey] < angle.limit.Z)
  
  ci <- ci[select]
  
  # =============================================================
  # === just in case the zenith angle filtering is not enough ===
  ci[abs(ci)==Inf] <- 1
  ci[is.nan(ci)] <- 1
  ci[is.na(ci)] <- 1
  ci[ci>=(1.36/0.9)] <- 1
  ci[ci<0] <- 0
  #==============================================================
  
  
  model <- toupper(model)
  if(identical(model, "ARIMA")) {
    fit <- auto.arima(ci, seasonal = F)
  } else if(identical(model, "ETS")){
    fit <- ets(ci, model = "ZZN") # model=error type, trend type, season type
  }
  
  out <- forecast(fit, h=1) # forecast the next period (one-step-ahead)
  pred.ci <- out$mean[1]
  
  # calc clear sky irradiance and forecast irradiance at forecast time stamp
  tm.fcst <- strftime(now+60*interval*0.5, tm.fmt)
  fcst.irr.clearsky <- calZen(as.POSIXct(tm.fcst, tz="GMT"), lat, lng, tm.zone)$Ic
  pred.irr <- round(pred.ci*fcst.irr.clearsky, 3)
  
  # write results to DB
  tm.fcst <- strftime(round(now+60*interval*0.5,"min"), tm.fmt) # round the forecast time stamp to the nearest minute
  sql.statement <- paste0("INSERT INTO SIN", station,"_FCST ", 
                          "(Tm, ", model, "_", interval, "m) VALUES ('", tm.fcst, "', ", pred.irr, ")",
                          " ON DUPLICATE KEY UPDATE ", model, "_", interval, "m=", pred.irr, ";")
  db.query(sql.statement)
  print(sql.statement)
  
  # also write the average of measured irradiances into DB, at middle of last interval
  tm.last.avg <- strftime(round(now-60*interval*0.5,"min"), tm.fmt) # round the forecast time stamp to the nearest minute
  sql.statement <- paste0("INSERT INTO SIN", station,"_FCST ", 
                          "(Tm, MEAS_", interval, "m) VALUES ('", tm.last.avg, "', ", tail(data.aggregated$Avg, n=1), ")",
                          " ON DUPLICATE KEY UPDATE MEAS_", interval, "m=", tail(data.aggregated$Avg, n=1), ";")
  db.query(sql.statement)
  
  # write to text file
  output.dir <- "/home/seris/ARIMA/output/SRA/"
  output.dir <- paste0(output.dir, station, "/", substr(tm.fcst, 1, 4), "/", substr(tm.fcst, 6, 7), "/")
  dir.create(file.path(output.dir), showWarnings = FALSE, recursive = TRUE) # create the directory if doesn't exist
  output.file <- paste0(output.dir, "SIN", station, "_", model, "_", interval, "min", " ", substr(tm.fcst, 1, 10), ".txt")
  
  if(!file.exists(output.file)){
    write("TimeStamp, Predicted Irradiance", file=output.file)
  }
  write(paste(tm.fcst,pred.irr, sep = ", "), file=output.file, append=TRUE)
  
  
  list(tm=tm.fcst, irr=pred.irr)
}



#===============================================================================
now <- round(Sys.time(),"mins")

if(now$hour>=7 && now$hour<=18) {
  stations <- c(401:425) # list of stations ("_AVG" corresponds to the avg of all 25 stn)

  # 60min ahead, 30min ahead and 15min ahead forecast with ETS model 
  pred.1h.ets <- sapply(stations, fcst, "ETS", now, 60, 24*7)
  pred.30min.ets <- sapply(stations, fcst, "ETS", now, 30, 2*24*3)
  pred.15min.ets <- sapply(stations, fcst, "ETS", now, 15, 4*24*3)
  
  # 60min ahead, 30min ahead and 15min ahead forecast with ARIMA model
  pred.1h.arima <- sapply(stations, fcst, "ARIMA", now, 60, 24*7)
  pred.30m.arima <- sapply(stations, fcst, "ARIMA", now, 30, 2*24*3)
  pred.15m.arima <- sapply(stations, fcst, "ARIMA", now, 15, 4*24*3) 
  
  # 60min ahead, 30min ahead and 15min ahead forecast with all models, averaged over all stations
  algorithms <- c(ls(pattern = '^pred.+?'))
  avg.pred.irr <- 0
  for(i in 1:length(algorithms)) {
    avg.pred.irr[i] <- mean(as.numeric(get(algorithms[i])['irr',]))
  }
  df <- data.frame(intervals=c(15,15,60,60,30,30), models=c("ARIMA","ETS","ARIMA","ETS","ARIMA","ETS"), avg.pred.irr, row.names = c(algorithms))
  for(algorithm in rownames(df)) {
    interval = df[algorithm, "intervals"]
    model = df[algorithm, "models"]
    pred.irr = df[algorithm, "avg.pred.irr"]
    tm.fmt <- "%Y-%m-%d %H:%M:%S"
    tm.fcst <- strftime(round(now+60*interval*0.5,"min"), tm.fmt) # round the forecast time stamp to the nearest minute
    sql.statement <- paste0("INSERT INTO SIN_AVG_FCST ", 
                            "(Tm, ", model, "_", interval, "m) VALUES ('", tm.fcst, "', ", pred.irr, ")",
                            " ON DUPLICATE KEY UPDATE ", model, "_", interval, "m=", pred.irr, ";")
    db.query(sql.statement)
    print(sql.statement)  
  }
              
}



