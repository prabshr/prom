undef("pm25ART")
;Estimate PM25 chemical composition of aerosol
;Vertiacal profiles of 10 differentiated species are returned
;AQMEII-2 project setup

function pm25ART(fil_cos,opt)

begin
  iy    = opt@iy
  ix    = opt@ix
  favg  = opt@favg
  fdry  = opt@fdry

  nc_cos   = addfile(fil_cos,"r")

; Do not touch below
;---------------------------------------------------------------------------------
  if (favg) then
    print("Processing " + fil_cos + " domain average")
  else
    print("Processing " + fil_cos + " at " + iy + " " + ix)
  end if

  ;Read vertical co-ordinate
  hhl = nc_cos->vcoord
  nz  = dimsizes(hhl)-1
  mhgt = new(nz,"float")
  mhgt = (/0.5*(hhl(0:nz-1) + hhl(1:nz))/)   ;agl

  ;1.Secondary Inorganic
  pm25_so4 = (nc_cos->VSO4I) + (nc_cos->VSO4J) + (nc_cos->VSO4Im) + (nc_cos->VSO4Jm) + \
             (nc_cos->VSO4SEASA) + 0.5*(nc_cos->VSO4SEASB)
  pm25_no3 = (nc_cos->VNO3I) + (nc_cos->VNO3J) + (nc_cos->VNO3Im) + (nc_cos->VNO3Jm)
  pm25_nh4 = (nc_cos->VNH4I) + (nc_cos->VNH4J) + (nc_cos->VNH4Im) + (nc_cos->VNH4Jm)
  ;2.Soot
  pm25_ec  = (nc_cos->VSOOTI) + (nc_cos->VSOOTJ) + (nc_cos->VSOOT)
  ;3.Primary organic
  pm25_poc = (nc_cos->VORGPAI) + (nc_cos->VORGPAJ) + (nc_cos->VORGPAIm) + (nc_cos->VORGPAJm)
  ;4.Primary inorganic
  pm25_inorg = (nc_cos->VP25AI) + (nc_cos->VP25AJ)
  ;5.Seconadry organics
  pm25_soa = (nc_cos->VORG1I) + (nc_cos->VORG1J) + (nc_cos->VORG10I) + (nc_cos->VORG10J) + \
             (nc_cos->VORG100I) + (nc_cos->VORG100J) + (nc_cos->VORG1000I) + (nc_cos->VORG1000J)
  ;6.Water vapor
  pm25_h2o = (nc_cos->VH2OI) + (nc_cos->VH2OJ) + (nc_cos->VH2OIm) + (nc_cos->VH2OJm) + \
             (nc_cos->VH2OSEASA) + 0.5*(nc_cos->VH2OSEASB)
  ;7.Crustal
  pm25_crust = (nc_cos->VSOILA)
  ;8.Sea salt
  pm25_seas  = (nc_cos->VSEASA) + 0.5*(nc_cos->VSEASB)

  ;Number of classified pm25 species
   if (fdry) then
     nc = 3*1 + 6
   else
     nc = 3*1 + 7
   end if

   pm25     = new((/nc,nz/),"float")
   pm25!0   = "ic"
   pm25!1   = "Z"
   pm25&ic  = ispan(0,nc-1,1)
   pm25&Z   = mhgt
   pm25@units    = pm25_crust@units
   pm25@long_name     = "PM2.5"
   if (fdry) then
      pm25@colors   = (/"goldenrod1", "darkseagreen", "palevioletred1", "grey29", "mistyrose", \
                     "mediumorchid2","powderblue","sandybrown","skyblue4"/)
     pm25@names    = (/"SO4", "NO3", "NH4", "EC", "POC", "INORG", "SOA", "CRUST", "SEAS"/)
   else
      pm25@colors   = (/"goldenrod1", "darkseagreen", "palevioletred1", "grey29", "mistyrose", \
                     "mediumorchid2","powderblue","sandybrown","skyblue4","royalblue4" /)
     pm25@names    = (/"SO4", "NO3", "NH4", "EC", "POC", "INORG", "SOA", "CRUST", "SEAS","H2O"/)
   end if

   if (favg) then
     pm25(0,:) = dim_avg_n(dim_avg_n(pm25_so4(0,:,:,:),2),1)
     pm25(1,:) = dim_avg_n(dim_avg_n(pm25_no3(0,:,:,:),2),1)
     pm25(2,:) = dim_avg_n(dim_avg_n(pm25_nh4(0,:,:,:),2),1)
     pm25(3,:) = dim_avg_n(dim_avg_n(pm25_ec(0,:,:,:),2),1)
     pm25(4,:) = dim_avg_n(dim_avg_n(pm25_poc(0,:,:,:),2),1)
     pm25(5,:) = dim_avg_n(dim_avg_n(pm25_inorg(0,:,:,:),2),1)
     pm25(6,:) = dim_avg_n(dim_avg_n(pm25_soa(0,:,:,:),2),1)
     if (.not.(fdry)) then
     pm25(9,:) = dim_avg_n(dim_avg_n(pm25_h2o(0,:,:,:),2),1)
     end if
     pm25(7,:) = dim_avg_n(dim_avg_n(pm25_crust(0,:,:,:),2),1)
     pm25(8,:) = dim_avg_n(dim_avg_n(pm25_seas(0,:,:,:),2),1)
   else
     pm25(0,:) = (/pm25_so4(0,:,iy,ix)/)
     pm25(1,:) = (/pm25_no3(0,:,iy,ix)/)
     pm25(2,:) = (/pm25_nh4(0,:,iy,ix)/)
     pm25(3,:) = (/pm25_ec(0,:,iy,ix)/)
     pm25(4,:) = (/pm25_poc(0,:,iy,ix)/)
     pm25(5,:) = (/pm25_inorg(0,:,iy,ix)/)
     pm25(6,:) = (/pm25_soa(0,:,iy,ix)/)
     if (.not.(fdry)) then
     pm25(9,:) = (/pm25_h2o(0,:,iy,ix)/)
     end if
     pm25(7,:) = (/pm25_crust(0,:,iy,ix)/)
     pm25(8,:) = (/pm25_seas(0,:,iy,ix)/)
   end if

   return(pm25(:,::-1))   
end
