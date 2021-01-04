;Plot the plan veiw of elevation at different model levels in COSMO
;This is important for exploring PPI of polarimetric variables
;Usage: need to specify iz (model level, 0 -surface)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin
  diri     = "/daten04/pshrestha/shared-data/bonnRadar_20150705"
  diri_ivr = diri + "/HET2483/runE_0_HET2483/cosout/ivr/"

  nc_ivr   = addfile(diri_ivr + "/lfff00000000c.nc","r")
  hhl      = nc_ivr->HHL(0,::-1,:,:)

  iz       = 21 ;40 

  ;Estimate dz
  hhldz     = hhl(iz+1,:,:) - hhl(iz,:,:)
  hhl2d     = hhl(iz,:,:)
  copy_VarCoords(hhl2d,hhldz)

  plot     = new(2,"graphic")
  res      = True
  wks      = gsn_open_wks("x11", "Figure")

  res@gsnDraw    = False
  res@gsnFrame   = False
  res@cnFillOn   = True
  res@cnFillMode = "CellFill"
  res@cnLinesOn  = False

  ext            = 1.3
  res@trXMinF    = -ext
  res@trXMaxF    = ext

  res@trYMinF    = -ext
  res@trYMaxF    = ext
  plot(0)  = gsn_csm_contour(wks,hhldz,res)
  plot(1)  = gsn_csm_contour(wks,hhl2d,res)
 
  gsn_panel(wks,plot,(/1,2/),False)

end
