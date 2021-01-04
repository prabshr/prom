;Plot the plan veiw of model variables at 
;a) Model Level
;b) Interpolated data at average model mass point height (a.m.s.l) for that model level
; Also overlaid in plot a) is the spatial pattern of model mass point height at that level

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$prom_root/functions/cosmo_time_nclf.ncl"

begin
 diri        = "/daten04/pshrestha/shared-data/bonnRadar_20150705/"
 itype_gscp  = "HET2483"
 mvar        = "W"
 inst        = 8
 iz          = 45 ;21(low level) ; 47(mid level) , 41? ;model level to plot
 fn         = 15;20    ;file to read
 ;
 glon_id     = (/7.071663/) ;
 glat_id     = (/50.73052/) ;
 ext         = 90

;;;;;;;;;;

 subevent    = "runE_" + inst + "_" + itype_gscp
 diri_cosI   = diri + "/" + itype_gscp + "/" + subevent + "/cosout/ivr/"
 diri_cos    = diri + "/" + itype_gscp + "/" + subevent + "/cosout/"

 nc_cosI     = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat         = nc_cosI->lat
 lon         = nc_cosI->lon
 topo        = nc_cosI->HSURF(0,:,:)

 ;Need to reverse the vertical levels
 hhl         = nc_cosI->HHL(0,:,:,:)
 ndimhhl     = dimsizes(hhl)
 nz          = ndimhhl(0)-1
 ny          = ndimhhl(1)
 nx          = ndimhhl(2)
 zmodel      = 0.5*(hhl(0:nz-1,:,:) + hhl(1:nz,:,:))  ;m
;
 ;Average height at iz
 zmodel      = zmodel(::-1,:,:) 
 mhgt        = zmodel(iz,:,:)
 avghgt      = avg(mhgt)
 print("Average model mass point height at iz =" +iz + " is " + avghgt + " m a.m.s.l")

 fils_cos    = systemfunc("ls " + diri_cos + "lf*")
 nc_cos      = addfile(fils_cos(fn),"r")
 time        = nc_cos->time
 str_time    = cosmo_time_nclf(time,1)  ; 0 for date
 print(str_time + " " + fils_cos(fn))

 vartemp     = nc_cos->$(mvar)$ 
 if (mvar.eq."W") then
   vartemp2    = 0.5*(vartemp(0,0:nz-1,:,:) + vartemp(0,1:nz,:,:))
   ;Need to reverse the vertical levels
   varh        = vartemp2(::-1,:,:)    
 else
   varh = vartemp(0,::-1,:,:)
 end if

 ;Make comparable with radar extent
 nm          = getind_latlon2d (lat,lon, glat_id, glon_id)
 sclp        = tointeger((/nm(0,0) - ext, nm(0,1) - ext/))
 eclp        = tointeger((/nm(0,0) + ext, nm(0,1) + ext/))

 ;Variable to plot
 spol          = varh(iz,:,:)         ;
 spol!0        = "Y"
 spol!1        = "X"
 spol&Y        = ispan(0,ny-1,1)  ;~ km
 spol&X        = ispan(0,nx-1,1)  ;~ km

 print(min(spol) + " " + max(spol))

 ;Now, interpolate also the variable to average model mass point height
 print("Interpolating " + mvar + " data at " + avghgt + " m")
 spol_ip       = wrf_user_intrp3d(varh,zmodel,"h",avghgt,0.,False)
 spol_ip@_FillValue = varh@_FillValue
 copy_VarCoords(spol,spol_ip)
 copy_VarCoords(spol,mhgt)


 print("Plotting ...")

;Differences
 spol_1d  = ndtooned(spol)
 spol_ip_1d  = ndtooned(spol_ip)
 rms      = dim_rmsd(spol_1d,spol_ip_1d)
 print("RMSD " + rms)

 plot     = new(2,"graphic")
 res      = True
 wks      = gsn_open_wks("x11", "Figure")

 res@gsnDraw    = False
 res@gsnFrame   = False
 res@cnFillOn   = True
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 ;res@cnLevelSelectionMode       = "ExplicitLevels"
 ;res@cnLevels                   =  ipol_levs

 res@trXMinF = sclp + 60 
 res@trYMinF = sclp + 60
 res@trXMaxF = eclp - 40
 res@trYMaxF = eclp - 40

 res@tiMainString = "Model level " + iz
 plot(0)  = gsn_csm_contour(wks,spol,res)
 res@tiMainString = "Avg. model height " + avghgt + " m"
 plot(1)  = gsn_csm_contour(wks,spol_ip,res)

 delete(res@cnLevelSelectionMode)
 delete(res@cnLevels)
 res@cnFillOn = False
 res@cnLinesOn = True
 res@cnLineLabelsOn = False
 res@cnLineColor    = "Gray"
 dum      = gsn_csm_contour(wks,mhgt,res)
 overlay(plot(0),dum)
 
 gsn_panel(wks,plot,(/1,2/),False)


end
