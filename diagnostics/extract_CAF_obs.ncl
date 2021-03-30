;This is a iterative process,
;First generate the CAF and debug ouput
;Based on debug ouput, mask the CAF for non-valid times (e.g.,when storm is out of extent)
;Re-run script
;See STEP EVAL section

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
 dbz_th    = 40.  ;dBZ
 dbz_min   = 0.  ;dBZ
 hgt_thr   = 2.  ;km
 ydate     = "20150705"
 radar     = "JuXPol"
 diri      = "/daten01/z4/database/PROM/" + ydate + "/" + radar + "/GRIDS/"

 gfils     = systemfunc("ls " + diri + "/*.nc")
 ndim      = dimsizes(gfils)

 time      = new(ndim,"float")
 do nf = 0, dimsizes(gfils)-1
   nc_grid  = addfile(gfils(nf),"r")
   zh_grid  = nc_grid->ZH_grid
   if (nf.eq.0)
     height = nc_grid->z
     x      = nc_grid->x
     y      = nc_grid->y
     iz_ind = ind(height.eq.hgt_thr)
     print("Height : " + height(iz_ind) + " km")
     Zh     = new((/ndim,dimsizes(y),dimsizes(x)/),"float")
     Zh@_FillValue = 999
   end if


   Zh(nf,:,:)    = nc_grid->ZH_grid(iz_ind,:,:)
   Rho           = nc_grid->RHO_grid(iz_ind,:,:)
   Vh            = nc_grid->VH_grid(iz_ind,:,:)

   ;Clutter Removal
   Zh(nf,:,:)    = where(((Rho .lt. 0.9) .and. (Zh(nf,:,:) .lt. 20) .and.(abs(Vh).lt.6)),Zh@_FillValue,Zh(nf,:,:))   
   tstr      = str_split(nc_grid@ScanStart,":")
   time(nf)  = tofloat(tstr(0)) + tofloat(tstr(1))/60.   ;hh + min/60.
   print(nf + " " + time(nf))
 end do

 if (radar.eq."BoXPolXXX") then
 ; STEP EVAL 2 Clutter Removal
 ix= (/(/101,305/),(/108,312/),(/115,321/),(/128,326/),(/134,329/),\
       (/156,334/),(/166,352/),(/178,332/),(/185,333/),(/190,344/)/)
 do nf = 10, 19 
  it = nf-10
  Zh(nf,:,ix(it,0):ix(it,1))   = Zh@_FillValue
 end do

 delete(ix)
 ix = (/(/0,198/),(/0,196/),(/0,198/),(/0,210/),(/113,239/),(/80,271/),\
        (/84,280/),(/94,290/),(/95,300/)/)
 iy = (/(/0,200/),(/0,210/),(/0,217/),(/0,223/),(/0,203/),(/0,202/),\
        (/0,207/),(/0,217/),(/0,227/)/)
 do nf = 27, 35
  it = nf-27
  Zh(nf,:,ix(it,0):ix(it,1))   = Zh@_FillValue
  Zh(nf,iy(it,0):iy(it,1),:)   = Zh@_FillValue
 end do
 end if

 varh       = new(ndim,"float")
 varh1      = new(ndim,"float")
 varh!0     = "time"
 varh&time  = time
 ;
 varh1!0    = "time"
 varh1&time = time
 ;

 do nf = 0, ndim-1
   val1d = ndtooned(Zh(nf,:,:))
   carea = tofloat(num(val1d.ge.dbz_th))
   tarea = tofloat(num(val1d.gt.dbz_min))
   if (tarea.gt.0) then    ;Avoid division by zero
     varh(nf) = carea/tarea
     varh1(nf) = tarea
   end if
   print("Convective area fraction " + carea + " " + tarea)
 end do

 ;;; STEP EVAL2 
 print("Generate debug output to diagnose extent of valid data")
 system("rm debug_" + radar + ".nc")
 fout = addfile("debug_" + radar + ".nc","c")
 fout->Band1 = Zh

 if (ydate .eq. "20150705" .and. radar .eq. "JuXPol") then
  varh(16:19) = varh@_FillValue 
  varh(29:) = varh@_FillValue
  varh1(29:) = varh1@_FillValue
  varh1(16:19) = varh1@_FillValue
 end if

 if (ydate .eq. "20150705" .and. radar .eq. "BoXPol") then
  a = array_append_record (ispan(0,9,1), ispan(20,26,1), 0)
  imsk = array_append_record (a, ispan(35,ndim-1,1),0)
  varh(imsk) = varh@_FillValue
  varh1(imsk) = varh1@_FillValue
 end if

 ;Plot
 ;---------------------------------------------------------------
 wks_type = "x11"
 wks = gsn_open_wks(wks_type,"Figure")

 res                            = True                 ; plot mods desired

 res@xyMarkLineMode = "MarkLines"
 plot    =  gsn_csm_xy (wks,time,varh,res) ; create plot


 foutname = "CAF_" + radar + "_" + ydate + ".nc"
 system("rm " + foutname)
 fout = addfile(foutname,"c")
 fout@date = ydate
 fout->Band1 = varh 
 fout->Band2 = varh1
 fout->time = time

end
