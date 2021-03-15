load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
 dbz_th    = 40.  ;dBZ
 dbz_min   = 0.  ;dBZ
 hgt_thr   = 2.  ;km
 ydate     = "20150705"
 diri      = "/daten01/z4/database/PROM/" + ydate + "/GRIDS/"

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
   Zh(nf,:,:)  = nc_grid->ZH_grid(iz_ind,:,:)
   tstr      = str_split(nc_grid@ScanStart,":")
   time(nf)  = tofloat(tstr(0)) + tofloat(tstr(1))/60.   ;hh + min/60.
   print(nf + " " + time(nf))

 end do

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

 ;Plot
 ;---------------------------------------------------------------
 wks_type = "x11"
 wks = gsn_open_wks(wks_type,"Figure")

 res                            = True                 ; plot mods desired

 res@xyMarkLineMode = "MarkLines"
 plot    =  gsn_csm_xy (wks,time,varh,res) ; create plot


 foutname = "CAF_BoXPol_" + ydate + ".nc"
 system("rm " + foutname)
 fout = addfile(foutname,"c")
 fout@date = ydate
 fout->Band1 = varh 
 fout->Band2 = varh1
 fout->time = time

end
