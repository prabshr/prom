;Plot the vertically streched co-ordinates in COSMO
;Top in native, Bottom in linearized plot
;This is important for exploring RHI of polarimetric variables

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
  diri     = "/daten04/pshrestha/shared-data/bonnRadar_20150705"
  diri_ivr = diri + "/HET2483/runE_0_HET2483/cosout/ivr/"

  nc_ivr   = addfile(diri_ivr + "/lfff00000000c.nc","r")
  hhl      = nc_ivr->HHL(0,::-1,:,:)

  ;lets look at jy  and ix
  jy       = 173    ;y-cs
  ix       = 137    ;x-cs
  nz       = 40 

  ;Estimate dz
  hhl2dy    = hhl(1:nz+1,jy,:) - hhl(0:nz,jy,:)
  hhl2y     = hhl(0:nz,jy,:)
  ydims     = dimsizes(hhl2y)
  hhl2dx    = hhl(1:nz+1,:,ix) - hhl(0:nz,:,ix)
  hhl2x     = hhl(0:nz,:,ix)
  xdims     = dimsizes(hhl2x)
;
  xval2y    = conform_dims(ydims,hhl2y&rlon,1) 
  xval2x    = conform_dims(xdims,hhl2x&rlat,1) 
  zval2y    = hhl(0:nz,jy,:)
  zval2x    = hhl(0:nz,:,ix) 

  plot     = new(4,"graphic")
  res      = True
  wks      = gsn_open_wks("x11", "Figure")

  res@gsnDraw    = False
  res@gsnFrame   = False

  ;res@cnCellFillEdgeColor = "gray" 
  res@cnFillOn   = True
  res@cnFillMode = "CellFill"
  res@cnLinesOn  = False
  ;res@cnLevelSelectionMode = "ExplicitLevels"
  ;res@cnLevels   = (/0.,50.,100.,300.,500.,700.,1000.,3000.,5000.,7000/)
  plot(0)  = gsn_csm_contour(wks,hhl2dy,res)
  plot(1)  = gsn_csm_contour(wks,hhl2y,res)

;
; Trick to linearize the plot
  res@trGridType = "TriangularMesh"
  res@sfXArray   = xval2y
  res@sfYArray   = zval2y

  plot(2)  = gsn_csm_contour(wks,hhl2dy,res)
  plot(3)  = gsn_csm_contour(wks,hhl2y,res)
 
  gsn_panel(wks,plot,(/2,2/),False)

end
