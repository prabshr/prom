;Plot the c-s veiw of hydrometeors alond Y-axis
; specified model grid in X using linearized 

;P. Shrestha
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$PROM/functions/shade_terrain.ncl"
load "$prom_root/functions/cosmo_time_nclf.ncl"
load "$prom_root/functions/calc_rho_nclf.ncl"
load "$prom_root/functions/convol.ncl"


begin
;STEP 1 User Settings
;-----------------------------------------------------------------------------------------
 pid       = 0 
 opt       = 1   ;0-NX, 1-QX , 2-DX

 ydates    = (/"20150705","20160513","20170706"/)
 insts     = (/6, 8, 11/)
 fns       = (/ 9, 11,  16/)
 fne       = (/18, 20,  19/)

 itype     = "HET2483"
 diri      = "/daten04/pshrestha/shared-data/bonnRadar_" + ydates(pid)
 subevent  = itype + "/runE_" + insts(pid) + "_" + itype
 diri_cos  = diri + "/" + "/" + subevent + "/cosout/"
 diri_cosI = diri + "/" + "/" + subevent + "/cosout/ivr/"

;STEP 2 Set variables to read and specify scales or offsets
;-------------------------------------------------------------------------------------
 hydrometeors_qx = (/"QC","QR","QI","QS","QG","QH"/)
 hydrometeors_qn = (/"NCCLOUD","NCRAIN","NCICE","NCSNOW","NCGRAUPEL","NCHAIL"/)

 if (opt.eq.0) then
   hydrometeors = hydrometeors_qn
   scale_varh   = 1/1.E3
   qcnlevs      = (/1.1, 1.E1, 1.E2, 1.E3, 1.E4, 1.E5, 1.E6, 1.E7/)
   hydunits     = "L~S~-1~N~"
   hydname      = "Number Concentration"
 else if (opt.eq.1) then
   hydrometeors = hydrometeors_qx
   scale_varh   = 1000.  ;kg/kg to g/kg
   thr_min      = 0.01
   qcnlevs      = (/0.01,0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0/)
   hydunits     = "g/m~S~3~N~"
   hydname      = "Sp. Hydrometeor Content"
 else if (opt.eq.2) then
   x3           = 1./3.
   hydrometeors = hydrometeors_qx
   a_geo        = (/0.124    , 0.124  , 0.835 , 2.4    , 0.142     , 0.1366/)  ;mkg-b
   b_geo        = (/x3       , x3     , 0.39  , 0.455  , 0.314     , x3/)      ;-
   xmin         = (/4.2E-15  , 2.6E-10, 1.E-12, 1.E-10, 1.E-9, 2.6E-9/)        ;kg
   xmax         = (/2.6E-10  , 3.0E-6 , 1.E-6 , 2.E-5 , 5.E-4, 5.E-4/)         ;g
   qmax         = 30.
   qmin         = 0.001
   qcnlevs       = (/0.001, 0.1, 1.0, 5.0, 10., 15., 20., 25., 30./)
   scale_varh   = 1000.   ; m to  mm
   hydunits     = "mm"
   hydname      = "Max. Volume Diameter"
 else
   print("code not written for opt = " + opt)
   exit
 end if
 end if
 end if
;

;STEP 3  READ NC FILES
;-------------------------------------------------------------------------------------
 nc_cosI   = addfile(diri_cosI +"lfff00000000c.nc","r")
 lat       = nc_cosI->lat
 lon       = nc_cosI->lon
 topo      = nc_cosI->HSURF(0,:,:)

 hhl       = nc_cosI->HHL(0,::-1,:,:)
 ndimhhl   = dimsizes(hhl)
 nz        = ndimhhl(0)-1
 ny        = ndimhhl(1)
 nx        = ndimhhl(2)
 zmodel    = 0.5*(hhl(0:nz-1,:,:) + hhl(1:nz,:,:))/1000.  ;km
;
 fils_cos  = systemfunc("ls " + diri_cos + "lf*")

 do fn = fns(pid),fne(pid)

 nc_cos    = addfile(fils_cos(fn),"r")
 time      = nc_cos->time
 ;
 str_time  = cosmo_time_nclf(time,1)  ; 0 for date
 print("Processing ..." + fn + "  " + str_time)

;STEP 4 Prognostic state variables
;------------------------------------------------------------------------------
 tair            = nc_cos->T(0,::-1,:,:)     ;K
 pressure        = nc_cos->P(0,::-1,:,:)     ;Pa
 qv              = nc_cos->QV(0,::-1,:,:)    ;kg/kg
 ux              = nc_cos->U(0,::-1,:,:)
 uy              = nc_cos->V(0,::-1,:,:)
 uz              = nc_cos->W(0,::-1,:,:)
 uz_hhl          = uz(1:80,:,:)*0.
 do kk = 0,79
   uz_hhl(kk,:,:)  = 0.5*(uz(kk,:,:) + uz(kk+1,:,:))
 end do

;STEP 5 Extract Hydrometeors
;-------------------------------------------------------------------------------------

; Compute ambient air density
 qc            = nc_cos->QC(0,::-1,:,:)    ;kg/kg
 qr            = nc_cos->QR(0,::-1,:,:)    ;kg/kg
 qi            = nc_cos->QI(0,::-1,:,:)    ;kg/kg
 qs            = nc_cos->QS(0,::-1,:,:)    ;kg/kg
 qg            = nc_cos->QG(0,::-1,:,:)    ;kg/kg
 qh            = nc_cos->QH(0,::-1,:,:)    ;kg/kg
 tlc           = qc + qr + qi + qs + qg + qh    ;Precipitation loading
 rho           = calc_rho_nclf(tair,pressure,qv,tlc) ;air density (kgm-3)

 nh           = dimsizes(hydrometeors)
 varh         = new((/nh,nz,ny,nx/),"float")    ;Initialize variable for plotting
 varh@units   = hydunits
 varh@long_name = hydname
 varh@_FillValue = -999.
 do ih = 0, nh-1
   if (opt.eq.0) then
     varh(ih,:,:,:) = (/nc_cos->$(hydrometeors_qn(ih))$(0,::-1,:,:)/)
     varh(ih,:,:,:) = (/scale_varh * varh(ih,:,:,:) * rho/)    ;(kg-1) to (L-1)
           ;Nt    = 1/1000. * Nd * rho            ; number density to number concentration
   else if (opt.eq.1) then
     varh(ih,:,:,:) = (/nc_cos->$(hydrometeors_qx(ih))$(0,::-1,:,:)/)
     varh(ih,:,:,:) = (/scale_varh * varh(ih,:,:,:) * rho/)    ;(kg/kg) to (g/m3)
     ;clipping
     varh(ih,:,:,:) = where(varh(ih,:,:,:).gt.thr_min,varh(ih,:,:,:),0) 
          ;Qt    = 1000. * Qd * rho
   else if (opt.eq.2) then
     qtemp      = (/nc_cos->$(hydrometeors_qx(ih))$(0,::-1,:,:)/)
     Ntemp      = (/nc_cos->$(hydrometeors_qn(ih))$(0,::-1,:,:)/)
     ;qtemp      = where(qtemp.lt.0.0001, 0., qtemp)
     Ntemp      = where(Ntemp.eq.0, varh@_FillValue, Ntemp)
     xmass      = (/qtemp/Ntemp/)   ; kg
     ;Clipping
     xmass      = where(xmass.gt.xmin(ih),xmass,0)
     xmass      = where(xmass.gt.xmax(ih),xmax(ih),xmass)
     ;diameter mass relationship
     varh(ih,:,:,:) = (/scale_varh* (a_geo(ih) * (xmass)^b_geo(ih))/)   ;mm
   end if
   end if
   end if
 end do

 hmask = dim_sum_n(varh(:,:,58:237,65:244),0)
 if (fn.eq.fns(pid)) then
   nvdim = dimsizes(hmask)
   varo  = new((/fne(pid)-fns(pid)+1,nvdim(0),nvdim(1),nvdim(2)/),"float")
 end if
 varo(fn-fns(pid),:,:,:)  = where(hmask.gt.thr_min, 1, 0)

 end do

 fname = "hydro.nc"
 system("rm " + fname)
 fout = addfile(fname,"c")
 fout->Band1 = varo
end
