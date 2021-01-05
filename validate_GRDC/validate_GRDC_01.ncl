;Filter GRDC data for country(DE) and year of data(2008>)
;Discharge Units:  m3/s 
; Average last 24 hours of model data to daily average
; Compute day of year and extract GRDC data for that day

;Outputs discharge data both for model and measurements
;Additionally taylor statistics are also computed annually
; fplot flag generates time-series plots for comparison

; Generally, since the modeled river widths are around 1132 m, does not 
; make sense to calculate discharge for small rivers in catchment, whereas
; the big rivers miss boundary flows

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$PROM/functions/readGRDC.ncl"

begin
  dir_GRDC = "/p/project/chbn33/hbn331/database/bonnRadar/validation/NC_GRDC"
  dir_TSMP = "./data_prom_"
  opt       = False   ;set to False if comparing with daily output
  opt@freq  = 5 
  nens      = 20
  ;year      = 2015
  inst      = 0       ;to extract day of year only
; 
  stid_GRDC = dir_GRDC + "/GRDC_level_3.txt"
  fils_TSMP = dir_TSMP + year + "/debugTSMP_" + year + "_" + inst + ".nc"

  nc_TSMP   = addfile(fils_TSMP,"r")
  lat2d     = nc_TSMP->lat2d
  lon2d     = nc_TSMP->lon2d
  time      = nc_TSMP->time
  ntim      = dimsizes(time)
  ntims     = ntim-24

  print("Obtaining day of year ...")

  utc_date  = cd_calendar(time(ntims), 0)
  year      = tointeger(utc_date(:,0))    ; Convert to integer for
  month     = tointeger(utc_date(:,1))    ; use sprinti 
  day       = tointeger(utc_date(:,2))
  doy       = day_of_year(year,month,day)
  print(year + " " +  month  + " " +  day + " " + doy)

  ;

  stid      = asciiread(stid_GRDC,-1,"string")
  ndim      = dimsizes(stid)
  
  ;Initialize arrays
  ;For annual discharge comparison output
  data      = new((/2,nens,ndim/),"float")
  data@_FillValue = -999
  data@long_name  = "streamflow"

  ; do not change loop below
  do nf = 0, dimsizes(stid)-1
    temp   = str_split(stid(nf),";")
    gridno = temp(0)
    river  = temp(1)
    jy     = tointeger(temp(6))
    ix     = tointeger(temp(7))

    print("Debug " + nf + " " + gridno + " " + river + " " + jy + " " + ix)
    fil_GRDC = dir_GRDC + "/" + gridno + ".nc" 

    qGRDC  = readGRDC(fil_GRDC,year,opt)
    if (qGRDC@status.eq.0) then
      data(0,0,nf) = (/qGRDC(doy)/)
      do ni = 0, nens-1   ;ensemble number
        fil_TSMP = dir_TSMP + "/debugTSMP_" + year + "_" + ni + ".nc"
        nc_TSMP  = addfile(fil_TSMP,"r")
        qTSMP    = dim_avg_n((nc_TSMP->QO(ntims:ntim-1,:,:,:)),0)  ; daily average discharge
        qom  = -1.*(qTSMP(0,jy,ix) + qTSMP(1,jy,ix))   ;m3/s
        data(1,ni,nf) = (/tofloat(qom)/)
      end do   
    end if
    delete([/qGRDC/])
  end do

  ; 
  ;opt      = True
  ;opt@fout = "discharge_" + year
  ;opt@title= "Model,Observation"
  ;fmtf     = "f8.2,1x,f8.2"
  ;write_matrix (data, fmtf, opt)

  fname  = "discharge_" + year + ".nc"
  system("rm " + fname)
  fout = addfile(fname,"c")
  fout->Band1 = data
end
