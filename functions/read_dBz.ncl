load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"
external EXT1 "/home/pshrestha/prom/functions/readdim.so"
external EXT2 "/home/pshrestha/prom/functions/readascii.so"

function read_dBz(diri,varid,nf,fopt)

begin

;User Settings
;---------------------------------------------------------------------------
 ;diri         = "/daten04/pshrestha/shared-data/bonnRadar_20150705/Synth/"
 dataFillVal  = fopt@FillValue
 asciiType    = fopt@asciiType

 varname      = (/"DBZ","zh" ,"zdr" ,"kdp" ,"rhv" ,"dlt" ,"ah", "ade"/)
 long_name    = (/"DBZ", "Z~B~H", "Z~B~DR", "K~B~DP", "~F5~r~B~~F21~HV~N~", "DLT", "AH","ADE"/)
 units        = (/"dBZ", "dBZ", "dB", "deg/km", " "/)


 diri_input   = diri + "/input"
 diri_output  = diri + "/output"

 fils_var     = systemfunc("ls " + diri_output + "/dbz-tm-1m_*")
 fils_elv     = systemfunc("ls " + diri_input + "/elevation_*")

 ndata        = dimsizes(fils_var)
 begTime = get_cpu_time()

 ;Read header information to extract nx,ny,nz
 array = new(3,integer)
 EXT1::readdim(array(0),array(1),array(2),fils_elv(nf))

 ;Identify index of file to read
 hgt3d = new((/array(2),array(1),array(0)/),"float")
 EXT2::readascii(array(0),array(1),array(2),hgt3d,1,1,fils_elv(nf))

 if (fopt@hgt) then
    return(hgt3d)
 else

 ;Parse time
 filestr = str_split(fils_var(nf),"_")
 ndstr   = dimsizes(filestr)
 datestr = stringtochar(filestr(ndstr-3))
 timestr = stringtochar(filestr(ndstr-2))
 time = tofloat(tostring(timestr(0:1)))*24. + \
          tofloat(tostring(timestr(2:3))) + tofloat(tostring(timestr(4:5)))/60.
 time@units = "hours since " + tostring(datestr(0:3)) + "-" + \
                tostring(datestr(4:5)) + "-" + tostring(datestr(6:7)) + \
                " " + tostring(datestr(8:9)) + ":00:00"
 stime  = ut_string(time,"%Y%N%D%H%M")

 ;Read data
 fname   = fils_var(nf)
 print("Reading PV data ..")

 EXT1::readdim(array(0),array(1),array(2),fname)
 var_3d = new((/array(2),array(1),array(0)/),"float")
 varidf = varid+1    ;Fortran wrapper index starts with 1
 if (varid.eq.2) then
   EXT2::readascii(array(0),array(1),array(2),var_3d,8,varid,fname)
   var_3d = where(var_3d.ne.dataFillVal,var_3d,var_3d@_FillValue)
   varh = var_3d
   EXT2::readascii(array(0),array(1),array(2),var_3d,8,varidf,fname)
   var_3d = where(var_3d.ne.dataFillVal,var_3d,var_3d@_FillValue)
   varv = var_3d
   var_3d = varh - varv
 else
    EXT2::readascii(array(0),array(1),array(2),var_3d,8,varidf,fname)
 end if

;;
  var_3d@long_name = varname(varid)
  var_3d@units     = units(varid)
  var_3d@time      = stime

;   ;Convective Fraction Estimate based on Intensity only (Rickenbach_Rutledge1998)
;   if (nf.eq.0) then
;     varout = new((/ndata,array(2),array(1),array(0)/),"float")
;     convf = new(ndata,"float")
;   end if
;
;   varout(nf,:,:,:) = var_3d
;   r1d  = ndtooned(varout(nf,:,:,:))
;   ncf  = dimsizes(ind(r1d.ge.40))
;   nca  = dimsizes(ind(r1d.ne.varout@_FillValue))
;   convf(nf) = tofloat(ncf)/tofloat(nca)

   print("Extracted " + nf + " " + time + (get_cpu_time() - begTime) + " seconds")
   return(var_3d)
 end if

end
