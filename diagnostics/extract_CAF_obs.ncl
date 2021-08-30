;This is a iterative process,
;First generate the CAF and debug ouput with fcrop= False, fcrop2=False
;Based on debug ouput, mask the CAF for non-valid times (e.g.,when storm is out of extent)
;Copy the debug with fcrop=True, fcrop2=True to the scratch folder, as this mask is necessary for CFADS
;Re-run script
;See STEP EVAL section

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
 dbz_th    = 40.  ;dBZ
 dbz_min   = 0.  ;dBZ
 hgt_thr   = 2.  ;km
 pid       = 2 
 ydates    = (/"20150705","20160513","20170706"/)
 ydate     = ydates(pid)
 radar     = "BoXPol"
 diri      = "/daten01/z4/database/PROM/" + ydate + "/" + radar + "/GRIDS/"
 fcrop     = True   ;spatial filter
 fcrop2    = True  ;temporal filter

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

 if (fcrop) then
   print("Manual storm selection ....")
   nzh   = dimsizes(Zh)
   nt    = nzh(0)
   ny    = nzh(1)
   nx    = nzh(2)
   datac = new((/nt,ny,nx/),typeof(Zh))
   datac = 999
   maskc = datac

;PROM1
   if (pid.eq.0) then

   if (radar.eq."JuXPol") then
     yx = (/(/0,110,0,61/),(/0,124,0,62/), (/0,138,0,74/), (/0,157,0,104/),      \
        (/0,176,0,121/), (/0,190,0,138/), (/0,185,0,151/), (/0,194,0,160/) , \
        (/0,215,0,165/), (/0,230,0,178/), (/0,231,0,183/), (/0,238,0,191/),  \
        (/0,240,0,200/), (/0,258,78,214/),  (/0,243,40,218/), (/0,246,70,228/),  \
        (/0,266,110,269/), (/0,270,136,271/), (/128,276,156,278/), (/111,276,166,282/), \
        (/124,277,177,277/), (/127,285,194,289/), (/114,283,203,302/),(/119,289,208,315/), \
        (/139,289,230,327/), (/146,294,239,341/),(/156,305,245,353/), (/159,296,251,362/), \
        (/166,330,258,375/)  /)
   end if

   if (radar.eq."BoXPol") then
     yx  = (/(/0,125,0,50/), (/0,125,0,50/), (/0,125,0,50/), (/0,125,0,50/), \
         (/0,125,0,50/), (/0,139,0,60/), (/0,228,0,68/), (/0,235,0,80/), \
         (/0,250,0,85/), (/0,264,0,96/), (/0,273,0,99/), (/0,275,0,109/), \
         (/0,279,0,117/), (/0,290,0,125/), (/0,299,0,135/), (/0,305,0,157/), \
         (/0,299,0,170/), (/154,319,0,179/), (/158,317,0,186/), (/159,333,0,189/), \
         (/159,322,107,195/), (/165,300,130,202/), (/170,305,140,212/), (/175,319,147,236/), \
         (/185,334,166,249/), (/193,330,170,261/), (/195,361,180,273/), (/204,370,190,287/), \
         (/214,385,200,302/), (/218,390,210,308/), (/218,390,232,312/),(/218,390,243,nx-1/), \
         (/214,385,259,nx-1/), (/214,385,269,nx-1/) /)
   end if
  
   ntt = dimsizes(yx)
   do nf = 0, ntt(0)-1
     datac(nf,yx(nf,0):yx(nf,1),yx(nf,2):yx(nf,3)) = Zh(nf,yx(nf,0):yx(nf,1),yx(nf,2):yx(nf,3))
     maskc(nf,yx(nf,0):yx(nf,1),yx(nf,2):yx(nf,3)) = 1
   end do

   end if

;PROM2
   if (pid.eq.1 .and. radar.eq."BoXPol") then
     yx = (/0,ny-1,0,250/)
     datac   = Zh
     maskc   = 1
     do nf = 0, 26 
       datac(nf,yx(0):yx(1),yx(2):yx(3)) = 999 
       maskc(nf,yx(0):yx(1),yx(2):yx(3)) = 999
     end do
     yx = (/175,ny-1,0,230/) 
     do nf = 27,nt-1
       datac(nf,yx(0):yx(1),yx(2):yx(3)) = 999
       datac(nf,0:175,0:145) = 999
     end do
   end if

   if (pid.eq.2 ) then
     datac = Zh
     maskc = 1
   end if

   Zh = datac
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
 if (fcrop) then
   fout->Band2 = maskc
 end if

 if (fcrop2) then
 if (ydate .eq. "20150705" .and. radar .eq. "JuXPol") then
  varh(18:21) = varh@_FillValue 
  varh(29:) = varh@_FillValue
  varh1(29:) = varh1@_FillValue
  varh1(18:21) = varh1@_FillValue
 end if

 if (ydate .eq. "20150705" .and. radar .eq. "BoXPol") then
  a = array_append_record (ispan(0,8,1), ispan(19,25,1), 0)
  imsk = array_append_record (a, ispan(34,ndim-1,1),0)
  varh(imsk) = varh@_FillValue
  varh1(imsk) = varh1@_FillValue
 end if

 if (ydate .eq. "20170706" .and. radar .eq. "JuXPol") then
  varh(0:7) = varh@_FillValue
  varh1(0:7) = varh1@_FillValue
  varh(12:22) = varh@_FillValue
  varh1(12:22) = varh1@_FillValue
 end if

 if (ydate .eq. "20170706" .and. radar .eq. "BoXPol") then
  varh(0:9)   = varh@_FillValue
  varh1(0:9)  = varh@_FillValue
  varh(14:15) = varh@_FillValue
  varh1(14:15) = varh1@_FillValue
 end if


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
