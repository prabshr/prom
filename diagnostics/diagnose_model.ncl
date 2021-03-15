;Plot the plan and elevation veiw of model variables at 

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$prom_root/functions/cosmo_time_nclf.ncl"

begin
 diri        = "/daten04/pshrestha/shared-data/bonnRadar_20150705/"
 itype_gscp  = "HET2483"
 finterp     = False
 mvar        = "QG"
 inst        =  8 
; iz          = 44 ;21(low level) ; 47(mid level) , 41? ;model level to plot
; ix          = 90
; iy          = 90
 ;fn         = 10;20    ;file to read

 ;Bonn Radar
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
 vartemp     = (nc_cos->QR) + (nc_cos->QI) + (nc_cos->QS) + (nc_cos->QG)+ (nc_cos->QH)
 if (mvar.eq."W") then
   vartemp2    = 0.5*(vartemp(0,0:nz-1,:,:) + vartemp(0,1:nz,:,:))
   ;Need to reverse the vertical levels
   varh        = vartemp2(::-1,:,:)    
 else if (mvar.eq."T") then
   varh = vartemp(0,::-1,:,:) -273.14
   cnLevels = fspan(min(varh),max(varh),50)
 else
   varh = 1000.*vartemp(0,::-1,:,:)
   varh = mask(varh,varh.gt.0, True)
   print(mvar + " " + min(varh) + " " + max(varh))
  cnLevels = (/0.0001,0.005,0.01,0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0,10./)
 end if
 end if

 ;Crop model domain - Make comparable with radar extent
 nm          = getind_latlon2d (lat,lon, glat_id, glon_id)
 sclp        = tointeger((/nm(0,0) - ext, nm(0,1) - ext/))
 eclp        = tointeger((/nm(0,0) + ext, nm(0,1) + ext/))

 print(sclp + " " + eclp)

 varhc       = varh(:,sclp(0):eclp(0),sclp(1):eclp(1)) 
 lonhc       = lon(sclp(0):eclp(0),sclp(1):eclp(1))
 lathc       = lat(sclp(0):eclp(0),sclp(1):eclp(1))
 zmodelc     = zmodel(:,sclp(0):eclp(0),sclp(1):eclp(1))

 ndim        = dimsizes(varhc)
 nz          = ndim(0)
 ny          = ndim(1)
 nx          = ndim(2)

; xval2y    = conform_dims((/nz,nx/),lonc(iy,:),1)
;xval2x    = conform_dims((/nz,ny/),latc(:,ix),1)
 xval2x    = conform_dims((/nz,ny/),ispan(0,ny-1,1),1)
; zval2y    = zmodelc(:,iy,:)
 zval2x    = zmodelc(:,:,ix)

 ;Variable to plot
 spol          = varhc(iz,:,:)         ;
 spol!0        = "Y"
 spol!1        = "X"
 spol&Y        = ispan(0,ny-1,1)  ;~ km
 spol&X        = ispan(0,nx-1,1)  ;~ km

 print(min(spol) + " " + max(spol))

 ;
 spol1          = varhc(:,:,ix)         ;
 spol1!0        = "Z"
 spol1!1        = "Y"

 print("Plotting ...")

 plot     = new(2,"graphic")
 res      = True
 wks      = gsn_open_wks("x11", "Figure")

 res@gsnDraw    = False
 res@gsnFrame   = False
 res@cnFillOn   = True
 res@cnFillMode = "CellFill"
 res@cnLinesOn  = False
 res@cnLevelSelectionMode       = "ExplicitLevels"
 res@cnLevels   = cnLevels 

 res@trXMinF = 60
 res@trYMinF = 60
 res@trXMaxF = 140
 res@trYMaxF = 140

 res1 = res

 res@tiMainString = "Model level " + iz
 plot(0)  = gsn_csm_contour(wks,spol,res)

 ; Trick to linearize the plot
 res1@trGridType = "TriangularMesh"
 res1@sfXArray   = xval2x
 res1@sfYArray   = zval2x

 delete(res1@trYMinF)
 res1@trYMaxF = 15000 

 res1@tiMainString = "CS at X = " + ix
 plot(1)  = gsn_csm_contour(wks,spol1,res1)

 res@cnFillOn = False
 res@cnLinesOn = True
 res@cnLineLabelsOn = False
 res@cnLineColor    = "Gray"
 dum      = gsn_csm_contour(wks,mhgt,res)
 overlay(plot(0),dum)

 pres                  = True
 pres@gsLineThicknessF = 2.0 
 dum1 = gsn_add_polyline(wks,plot(0),(/ix,ix/), (/0,ny/), pres)
 dum2 = gsn_add_polyline(wks,plot(1),(/0,ny/), (/zval2x(iz,:)/), pres)

 gsn_panel(wks,plot,(/1,2/),False)


end
