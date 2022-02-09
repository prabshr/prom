;Plot aerosol size distribution for select modes
;opt@favg = True, returns spatial mean and standard deviation of each mode
; Data are ouput at specified heights in 4 different subplots
; fprcs = REad from model or porcessed data
load "$prom_root/functions/lndistART.ncl"

begin
;-----------------------------------------------------------------------------------------
;#1 USER SETTINGS
;-----------------------------------------------------------------------------------------
  pname = "EMPA2_PROM1";
  inst  = 6
  lnd   = "HET"
  itype = "ART"
  diri  = "/p/largedata/hbn33/shrestha1/"

  heights     = (/10., 1000., 2000., 4000., 7000./)
  ;heights     = (/2000.,4000./) ;
  colors      = (/"lightpink3","forestgreen","steelblue","gray20","gray"/)
  dashs       = (/0,0,0,0,0/)
  nf          = 32
  glon_idB    = 7.071663 + 13*0.01  ; Bonn Radar
  glat_idB    = 50.73052
  aerofloor   = 1.E-30   ; m3 (minimum concentration)
  fprcs       = False ;True
;--------------------------------------------------------------------------------
  itype_gscp  = lnd + itype
  sfold       = "runE_" + inst + "_" + itype_gscp
  expname     = lnd + "_" + pname
  subevent    = expname + "/" + itype_gscp + "/" + sfold
  diri_cos    = diri + "/" + subevent + "/a_cosout/"
  diri_cosI   = diri + "/" + subevent + "/a_cosout/ivr/"

  fils_cos = systemfunc("ls " + diri_cos + "*.nc")
  fil_cos  = fils_cos(nf)

;;--------------------------------------------------------------------------------
  ;Identify geo-location and extract model height
  nc_cosI = addfile(diri_cosI + "/lfff00000000c.nc","r")

  lat     = nc_cosI->lat
  lon     = nc_cosI->lon

  ;
  nm      = getind_latlon2d (lat,lon, glat_idB, glon_idB)
  yloc    = nm(0,0)
  xloc    = nm(0,1)

  print("Location " + yloc + " " + xloc)

;--------------------------------------------------------------------------------
  ;Geomeric standard deviations of lognormal modes
  sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

  ;Obtain aggregate of 12 log-normal distribution modes
  opt      = True
  opt@iy   = yloc 
  opt@ix   = xloc
  opt@favg = True  ;Returns average and std. deviation of each modes
  opt@fout = 0   ;Output ASD
  opt@nbdl = 20  ; Crop
  opt@spread = False ; Plot the spread

  smodes = ispan(1,12,1) 


 if (fprcs) then
;-----------------------------------------------------------------------------------------
;#2 Extract lognormal ASD for each mode 
;-----------------------------------------------------------------------------------------

  do im = 0, dimsizes(smodes)-1
    opt@smode = smodes(im)
    if (im.eq.0) then
      temp = lndistART(fil_cos,opt)
      ndim = dimsizes(temp)
      nlDp = new((/12,ndim(0),ndim(1),ndim(2)/),"float")
      units = temp@units
    end if
    nlDp(im,:,:,:) = lndistART(fil_cos,opt)
  end do
  Dp   = nlDp&Dp
  Z    = nlDp&Z
  nchan = nlDp@dlnDp

  ;Extract average or location value of dN/dlog(Dp)p
  ;Average is additive , E(X+Y) = E(X) + E(Y)

  nid = 0   ;avg
  nlDp1     = dim_sum_n(nlDp(0:3,nid,:,:),0)    ;nuc-acc mode ;small
;  nlDp2     = dim_sum_n(nlDp(4:5,nid,:,:),0)    ;soot +coarse
  nlDp2     = nlDp(4,nid,:,:)   ; soot only
    nlDp2  =   where(nlDp2.eq.0, aerofloor ,nlDp2)  ;Remove zero values with minimum
  nlDp3     = dim_sum_n(nlDp(6:8,nid,:,:),0)    ;dust
    nlDp3  =   where(nlDp3.eq.0, aerofloor ,nlDp3)    ;Remove zero values with minimum
  nlDp4     = dim_sum_n(nlDp(9:11,nid,:,:),0)  ;sea

  ;Extract standard deviation of dN/dlog(Dp)p
  ;Std[X+y] = SQRT[VAR(X) + VAR)Y)] assuming X and Y are independent, 
  ; but this might be a very strong assumption here !

  if (opt@favg) then
    nid = 1   ;stddev
 
    snlDp1     = sqrt((/nlDp(0,nid,:,:)^2 + nlDp(1,nid,:,:)^2 + nlDp(2,nid,:,:)^2 + nlDp(3,nid,:,:)^2/))    ;nuc-acc mode ;small
    ;;snlDp2     = dim_sum_n(nlDp(4:5,nid,:,:),0)    ;soot+coarse
    snlDp2     = nlDp(4,nid,:,:) ; soot
      snlDp2  =   where(nlDp2.eq.0, 0 ,snlDp2)    ;Set s.d = 0
    snlDp3     =  sqrt((/nlDp(6,nid,:,:)^2 + nlDp(7,nid,:,:)^2 + nlDp(8,nid,:,:)^2/))    ;dust
      snlDp3  =   where(nlDp3.eq.0, 0 ,snlDp3)      ;Set s.d = 0
    snlDp4     =  sqrt((/nlDp(9,nid,:,:)^2 + nlDp(10,nid,:,:)^2 + nlDp(11,nid,:,:)^2/))   ;sea
  end if
 
 ;Extract data at specific levels 
  nlevs       = dimsizes(heights)
  izlevs      = new(nlevs,"integer")
  do iz       = 0, nlevs-1
    izlevs(iz) = closest_val(heights(iz),Z)
  end do
  znlDp1       = transpose(nlDp1(:,izlevs))
  znlDp2       = transpose(nlDp2(:,izlevs))
  znlDp3       = transpose(nlDp3(:,izlevs))
  znlDp4       = transpose(nlDp4(:,izlevs))
  ;
  if (opt@favg) then
    zsnlDp1       = transpose(snlDp1(:,izlevs))
    zsnlDp2       = transpose(snlDp2(:,izlevs))
    zsnlDp3       = transpose(snlDp3(:,izlevs))
    zsnlDp4       = transpose(snlDp4(:,izlevs))
  end if

  ;Debug Output
  foutname = "asd_" + pname + "_ens" + inst + "_t" + nf + ".nc"
  system("rm " + foutname)
  fout = addfile(foutname,"c")
  fout->Dp = (/Dp/)
  fout->znlDp1 = (/znlDp1/)
  fout->znlDp2 = (/znlDp2/)
  fout->znlDp3 = (/znlDp3/)
  fout->znlDp4 = (/znlDp4/)
  if (opt@favg) then
    fout->zsnlDp1 = (/zsnlDp1/)
    fout->zsnlDp2 = (/zsnlDp2/)
    fout->zsnlDp3 = (/zsnlDp3/)
    fout->zsnlDp4 = (/zsnlDp4/)
  end if
  dummy = -999
  dummy@units = units
  dummy@nlevs = nlevs
  dummy@nchan = nchan
  fout->Band1 = dummy

  delete(dummy)
 else
   foutname = "asd_" + pname + "_ens" + inst + "_t" + nf + ".nc"
   fout = addfile("./data/" + foutname,"r")
   Dp=fout->Dp
   znlDp1=fout->znlDp1 
   znlDp2=fout->znlDp2 
   znlDp3=fout->znlDp3 
   znlDp4=fout->znlDp4 
   if (opt@favg) then
     zsnlDp1=fout->zsnlDp1 
     zsnlDp2=fout->zsnlDp2 
     zsnlDp3=fout->zsnlDp3 
     zsnlDp4=fout->zsnlDp4 
   end if
   dval  = fout->Band1
   units = dval@units
   nlevs = dval@nlevs
   nchan = dval@nchan
 end if    ;fprcs

 ; Diagnostics Print
  N_cn1 = ceil(dim_sum_n(znlDp1,1)) *nchan
  N_cn2 = ceil(dim_sum_n(znlDp2,1)) *nchan
  N_cn3 = ceil(dim_sum_n(znlDp3,1)) *nchan
  N_cn4 = ceil(dim_sum_n(znlDp4,1)) *nchan

  N_cn1 = decimalPlaces(N_cn1,2,True)
  N_cn2 = decimalPlaces(N_cn2,2,True)
  N_cn3 = decimalPlaces(N_cn3,2,True)
  N_cn4 = decimalPlaces(N_cn4,2,True)

  legend1 = heights + " m : N_x = " + N_cn1 + " " + units
  legend2 = heights + " m : N_x = " + N_cn2 + " " + units
  legend3 = heights + " m : N_x = " + N_cn3 + " " + units
  legend4 = heights + " m : N_x = " + N_cn4 + " " + units

;  print(N_cn1 + " " + N_cn2 + " " + N_cn3 + " " + N_cn4)

;-----------------------------------------------------------------------------------------
;#3 Plotting Resource
;-----------------------------------------------------------------------------------------

  plot = new(4,"graphic")

  wks_type = "png"
  fname_out = "Figure_ASD_" + nf + "_" + opt@nbdl ;+ sprinti("%0.2i",iz)
  wks = gsn_open_wks(wks_type,fname_out)

  fncps_no = ":F21:"
  fncps_no = 21
  fncps_hgt = 0.035

 res                            = True                 ; plot mods desired
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 res@gsnRightString             = ""
 res@gsnMaximize                = False
 res@gsnLeftStringFontHeightF   = fncps_hgt
 res@gsnRightStringFontHeightF  = fncps_hgt
 ;res@gsnLeftString              = Z(iz) + " m"

 res@tmXBLabelFont              = fncps_no
 res@tmYLLabelFont              = fncps_no
 res@tmYLLabelFontHeightF       = fncps_hgt
 res@tmXBLabelFontHeightF       = fncps_hgt
 res@tmYROn                     = False
 res@tmXTOn                     = False
 res@tmXBLabelFontColor         = "black"
 res@tmYLLabelFontColor         = "black"
 res@tmXBLabelFontThicknessF    = 2.0
 res@tmYLLabelFontThicknessF    = 2.0

 res@tiXAxisFont                = fncps_no
 res@tiYAxisFont                = fncps_no
 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt

 res@xyLineThicknessF           = 3.
 res@xyDashPattern              = 0
 res@xyLineColors                = colors
 res@tiXAxisString              = "Dp [nm]"
 res@tiYAxisString             = "dN/dlnDp [" + units + "]" 

 ;res@trXMinF = 1.   ;nm
 res@trYMinF = 1.E-3;dN/dlogDp, cm-3 
 res@trYMaxF = 1.E6 ;

 res@trYLog  = True
 res@trXLog = True
 res@xyComputeXMin = True

 res@gsnLeftString = "nuc./acc. (pure + mixed)"
 plot(0)  =  gsn_csm_xy (wks,Dp,znlDp1,res) ; create plot
 res@gsnLeftString = "soot" ; + coarse"
 plot(1)  =  gsn_csm_xy (wks,Dp,znlDp2,res)
 res@gsnLeftString = "dust"
 plot(2)  =  gsn_csm_xy (wks,Dp,znlDp3,res)

 res@gsnLeftString = "sea-salt"
 ;CPS plot(3)  =  gsn_csm_xy (wks,Dp,znlDp4,res)

 ; Add legend
 lgres                    = True
 lgres@lgLineColors       = colors
 lgres@lgLineThicknessF   = 3.
 lgres@lgDashIndexes      = dashs 
 lgres@lgItemType         = "Lines"
 lgres@lgLabelFontHeightF = fncps_hgt*17
 lgres@vpWidthF           = 0.36
 lgres@vpHeightF          = 0.20
 lgres@lgPerimOn          = False
 lbid1 = gsn_create_legend(wks,nlevs,legend1, lgres)
 lbid2 = gsn_create_legend(wks,nlevs,legend2, lgres)
 lbid3 = gsn_create_legend(wks,nlevs,legend3, lgres)
 lbid4 = gsn_create_legend(wks,nlevs,legend4, lgres)

 amres = True
 amres@amParallelPosF   = 0.15
 amres@amOrthogonalPosF = -0.35
 annoid1 = gsn_add_annotation(plot(0),lbid1,amres)
 amres@amParallelPosF   = 0.15
 annoid2 = gsn_add_annotation(plot(1),lbid2,amres)
 amres@amParallelPosF   = -0.15
 annoid3 = gsn_add_annotation(plot(2),lbid3,amres)
 amres@amParallelPosF   = -0.25
 ;CPS annoid4 = gsn_add_annotation(plot(3),lbid4,amres)


 if (opt@spread) then
   ng = 4*nlevs
   dummy = new(ng,"graphic")
   ictr=-1
   do ip = 0,3
   nDp = dimsizes(Dp)
   xp = new((/2*nDp/),"float")
   yp = new((/2*nDp/),"float")
   if (ip.eq.0) then
     znlDp = (/znlDp1/)
     zsnlDp = (/zsnlDp1/)
   end if
   if (ip.eq.1) then
     znlDp = (/znlDp2/)
     zsnlDp = (/zsnlDp2/)
   end if
   if (ip.eq.2) then
     znlDp = (/znlDp3/)
     zsnlDp = (/zsnlDp3/)
   end if
   if (ip.eq.3) then
     znlDp = (/znlDp4/)
     zsnlDp = (/zsnlDp4/)
   end if 
   do iz = 0, nlevs - 1
   ictr = ictr + 1
   do k=0,nDp-1
     yp(k)          = znlDp(iz,k) + zsnlDp(iz,k)
     xp(k)          = Dp(k)
     xp(2*nDp-1-k)  = Dp(k)
     yp(2*nDp-1-k)  = znlDp(iz,k) - zsnlDp(iz,k) 
   end do 
   print(min(xp) + " " + max(xp) + " " + min(yp) + " " + max(yp))

   gsres                   = True                        ; poly res

   gsres@tfPolyDrawOrder   = "Predraw"                   ; draw this first
   gsres@gsFillColor       = colors(iz) ;"SlateBlue"                 ; color chosen
   gsres@gsFillOpacityF     = 0.2
   dummy(ictr) = gsn_add_polygon (wks,plot(ip),xp,yp,gsres)
   end do
   end do ;ip
 end if

 resP  = True
 ;resP@gsnPanelMainString = nlDp@time 
 gsn_panel(wks,plot,(/2,2/), resP)


end
