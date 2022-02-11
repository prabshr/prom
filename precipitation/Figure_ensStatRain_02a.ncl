load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 

; First use dataExtract_rainfall_01.ncl to extract data
begin
;fenscolor =True gives different colors to each itype_gscp
;nes,and nee represents the start and end of ensemble members for each itype_gscp

 prom = "empa2_prom" + expid

 ;lcov       = "CASE I :  20150704 0300 - 20150706 0000"
; itype_gscp = (/"HET2483","HET2563","HET2683","TGR2483","TGR2563","RADOLAN"/)
; itype_gscp = (/"HET2483","HET2563","HET2683","RADOLAN"/)
; itype_gscp = (/"HET2483","TGR2483","CUR2483","RADOLAN"/)
; itype_gscp = (/"HET2483","HET2563","HET2683","RADOLAN"/)

 itype_gscp = (/"HETART","RADOLAN"/)
; itype_gscp = (/"HET2483","HET2463","HET2493","HET2583","HET2563","HET2593","HET2683","HET2663","HET2693", \
;                "EHD2483","EHD2463","EHD2493","EHD2583","EHD2563","EHD2593","EHD2683","EHD2663","EHD2693", \
;                "RHD2483","RHD2463","RHD2493","RHD2583","RHD2563","RHD2593","RHD2683","RHD2663","RHD2693", \
;                "RADOLAN"/)
 nes        = 5 ;0 
 nee        = 9 ;19 
 fenscolor  = False

 ;PDF parameters
 nbin         = 20
 opt          = True
 opt@bin_min  = 0.0     ;mm
 opt@bin_max  = 100.0    ;mm
 
 ;Open Plot resources
 fncps_no     = ":F21:"
 fncps_no     = 21
 fncps_hgt    = 0.02

 plot         = new(2,"graphic")

 fnameout     = "Figure_histo_" + expid
 wks          = gsn_open_wks("png",fnameout)
 res          = True

;gsn Resources
 ;-------------------------------------------
 res@gsnDraw                    = False
 res@gsnFrame                   = False
 ;res@gsnLeftString              = str_time
 res@gsnMaximize                = False
 ;res@gsnLeftStringFontHeightF   = fncps_hgt
 ;res@gsnRightStringFontHeightF  = fncps_hgt

 ;tm Resources
 ;-------------------------------------------
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

 res@tiXAxisFontHeightF         = fncps_hgt
 res@tiYAxisFontHeightF         = fncps_hgt
 ires = res

 ;vp Resources
 ;-------------------------------------------
 res@vpWidthF                   = 0.75
 res@vpHeightF                  = 0.50

 ;xy Resources
 ;-------------------------------------------
 res@xyMarkLineMode             = "MarkLines"
 res@xyMarkers                  = 16
 res@xyDashPatterns             = 3 
 res@xyLineThicknessF           = 2.0
 
 ;tr Resources
 ;-------------------------------------------
 res@trYLog                     = True
 res@trYMinF                    = 0.1    ; %0.1
 res@trYMaxF                    = 100.0  ; %
 res@trXMaxF                    = 60.    ; mm

 ;ti Resources
 ;-------------------------------------------
 res@tiXAxisString              = "Accumulated Precipitation (mm)"
 res@tiYAxisString              = "Frequency (%)"

 ctr = 0     ;counter for statistics
 do ityp  = 0, dimsizes(itype_gscp)-1           ;BIG LOOP
   res@xyDashPattern = ityp 
   if (ityp.eq.0) then
     nmax   = (nee-nes+1)*(dimsizes(itype_gscp)-1) + 1  ;Assume that we have 1 observation data
     pstat  = new((/4, nmax/),"float")
     pstat  = pstat@_FillValue
     barcolors = new(nmax,"string")
     barcolors = barcolors@_FillValue
   end if
   if (itype_gscp(ityp).eq."RADOLAN") then
     nes = 0
     nee = 0
   else
   end if

   do nf  = nes,nee ;0, nens-1                   ;EPS Loop
     filname    = "data_" + prom + "/accprecip_" + itype_gscp(ityp) + "_" + nf + ".nc"
     fout       = addfile(filname,"r")
     pvar       = fout->ACCRAIN
     pvar1d     = pdfx(ndtooned(pvar),nbin,opt) 
     pstat(0,ctr)  = ctr
     pstat(1,ctr)  = avg(pvar)
     pstat(2,ctr)  = stddev(pvar)
     pstat(3,ctr)  = pvar1d@nUse
     print(nf + " " + filname + " " + pstat(1,ctr) + " " + pstat(2,ctr) + " " + sum(pvar1d) + " " + pstat(3,ctr))
     ;  
     pvar1d      = pvar1d + 0.001
     ;pvar1d     = cumsum(pvar1d,2)
     if (fenscolor) then
      if (nf.le.4) then
        res@xyLineColor = "gray"
        res@xyMarkerColor = "gray"
        barcolors(ctr)  = "gray"
      end if
      if (nf.gt.4 .and. nf.le.9) then
        res@xyLineColor = "forestgreen"
        res@xyMarkerColor = "forestgreen"
        barcolors(ctr) = "forestgreen"
      end if 
      if (nf.gt.9 .and. nf.le.14) then
        res@xyLineColor = "gray"
        res@xyMarkerColor = "gray"
        barcolors(ctr) = "gray"
      end if
      if (nf.gt.14 .and. nf.le.19) then
        res@xyLineColor = "gray"
        res@xyMarkerColor = "gray"
        barcolors(ctr) = "gray"
      end if
     else
      res@xyLineColor = "gray"
      res@xyMarkerColor = "gray"
      barcolors(ctr) = "gray"
     end if

     if (itype_gscp(ityp).eq."RADOLAN") then
       res@xyLineColor = "black"
       res@xyMarkerColor = "black"
       res@xyLineThicknessF = 4.
       res@xyDashPattern  = 0
       barcolors(ctr) = "black"
     end if

     if (nf.eq.nes .and. ityp.eq.0) then
       res@tmXBMode                   = "Explicit"
       res@tmXBValues                 = pvar1d@bin_center + 2.5
       res@tmXBLabels                 = pvar1d@bin_center + 2.5
       ;res@gsnRightString             = lcov
       plot(0) = gsn_csm_xy(wks,pvar1d@bin_center, pvar1d,res)
     else
      ;delete(res@gsnRightString)
      dum = gsn_csm_xy(wks,pvar1d@bin_center, pvar1d,res)
      overlay(plot(0),dum)
     end if
     delete([/pvar,pvar1d/])

     ctr = ctr + 1
   end do  ;EPS LOOP
 end do  ;BIG LOOP

 ires = True
 ires@vpWidthF = 0.4
 ires@vpHeightF = 0.2
 ires@vpXF = 0.54
 ires@vpYF = 0.79

 ires@gsnXYBarChart = True            ; Create bar plot
 ires@gsnYRefLine    = 0
 ires@gsnXYBarChartColors = barcolors; 
 print(barcolors)

 ires@trYMinF = 0.
 ires@trYMaxF = 20.

 ires@tiYAxisString = "Avg. Precip (mm)"
 ires@tiXAxisString = "Ens. Members / RADOLAN"
 ;ires@tiXAxisString = "HET - EHD - RHD - RADOLAN"
 plot(1) = gsn_csm_xy (wks,pstat(0,:),pstat(1,:),ires)

 ;add stddev
  polyres                   = True                       ; marker resources
  polyres@gsMarkerIndex     = 1                          ; index for circle
  polyres@gsMarkerSizeF     = .02                        ; size 
  
  error_bar = new(ctr,graphic)
  centers   = new(ctr,graphic)

; loop through the points  
  do t=0,ctr-1
    centers(t)   = gsn_add_polymarker(wks,plot(1),pstat(0,t),pstat(1,t),polyres)
    error_bar(t) = gsn_add_polyline(wks,plot(1),(/pstat(0,t),pstat(0,t)/),(/pstat(1,t)+pstat(2,t),pstat(1,t)-0./),\
    polyres)
  end do


 draw(plot)
 frame(wks)

end

