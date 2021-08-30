;Plot the time-series of convective area fraction for all ensemble members
;at specified  model level
; CAF is the fractional area of storm where the ZH > 40 dBZ

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
external EXT1 "$prom_root/functions/readdim.so"
external EXT2 "$prom_root/functions/readascii.so"

begin
; diri_dBz    = "/daten04/pshrestha/shared-data/bonnRadar_20150705/Synth/"
 diri_dBz    = "/daten04/pshrestha/shared-data/bonnRadar_20160513/Synth/"
; diri_dBz    = "/daten04/pshrestha/shared-data/bonnRadar_20170706/Synth/"
; diri_dBz    = "/daten/daten04/pshrestha/shared-data/bonnRadar_20150705/Synth_track/"

 expid        = "HET2483"
 pvars       = (/"Zh","Zdr","KDP","RHV"/)
 punits      = (/"dBZ","dB","deg/km ", "-"/)
 iz          = 34 ;21(low level) ; 47(mid level) , 41? ;model level to plot
 dbz_th      = 40.  ;dBZ
 dbz_min     =  0.  ; -5.  ;dBZ
;;;;;;;;;;
 fn0         = 0 ; use only 1 elevation file

 diri_input  = diri_dBz + "/" + expid + "_" + inst + "/input"
 diri_output = diri_dBz + "/" + expid + "_" + inst +"/output"

 fils_var    = systemfunc("ls " + diri_output + "/dbz-tm-1m_*")
 fils_elv    = systemfunc("ls " + diri_input + "/elevation_*")
 fil_lat     = systemfunc("ls " + diri_input + "/lat_*")
 fil_lon     = systemfunc("ls " + diri_input + "/lon_*")

 ;Read header information to extract nx,ny,nz
 array = new(3,integer)
 EXT1::readdim(array(0),array(1),array(2),fils_elv(fn0))

 latPV = new((/1,1,array(1),array(0)/),"float")
 lonPV = new((/1,1,array(1),array(0)/),"float")
 hgtPV = new((/1,array(2),array(1),array(0)/),"float")
 EXT2::readascii(array(0),array(1),1,latPV,1,1,fil_lat)
 EXT2::readascii(array(0),array(1),1,lonPV,1,1,fil_lon)
 EXT2::readascii(array(0),array(1),array(2),hgtPV,1,1,fils_elv(fn0))   ;Model mass point heights
 
 ;Average height at iz
 ;Need to reverse the vertical levels
 hgtPV  = hgtPV(:,::-1,:,:) 
 mhgt   = hgtPV(0,iz,:,:) -  hgtPV(0,0,:,:)
 avghgt = avg(mhgt)
 print("Average model mass point height at iz =" +iz + " is " + avghgt + " m a.g.l")  ;~aprox 

 nfils        = dimsizes(fils_var)
 varh         = new(nfils,"float")    ;Initialize variable for plotting
 tvar         = new(nfils,"float")
 timeh        = new(nfils,"float")
 varh@_FillValue = -99
 varPV        = new((/8,array(2),array(1),array(0)/),"float")
 varPV@_FillValue = -99
 varD         = new((/nfils,array(1),array(0)/),"float")
 
 do fn = 0, nfils - 1    ;BIG LOOP 
   fname        = fils_var(fn)
   print(" " + fname)

   EXT2::readascii(array(0),array(1),array(2),varPV,8,1,fname)
   ;Need to reverse the vertical levels
   Zh  = varPV(1,::-1,:,:)
   val1d = ndtooned(Zh(iz,:,:))
   carea = tofloat(num(val1d.ge.dbz_th))
   tarea = tofloat(num(val1d.ge.dbz_min))
   if (tarea.gt.0) then    ;Avoid division by zero
     varh(fn) = carea/tarea
   end if
   tvar(fn) = tarea
   varD(fn,:,:) = Zh(iz,:,:)
   delete([/Zh, val1d/]) 
   print("Convective area fraction " + carea + " " + tarea)
 end do


 ;Debug output
 ;fdname = "debug_model_" + inst + ".nc"
 ;system("rm " + fdname)
 ;fout = addfile(fdname,"c")
 ;fout->Band1 = varD

 ;CAF output
 varh@model_level = iz
 varh@model_height = avghgt

 foname = "CAF_" + expid + "_" + inst + ".nc"
 system("rm " + foname)
 fout = addfile(foname,"c")
 fout->Band1 = varh
 fout->Band2 = tvar
end
