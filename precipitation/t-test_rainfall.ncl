;======================================================
; Diagnose difference in spatial mean rainfall between two simulations
; are statistically signficant, paired t-test with same variance assumption
;
; The precipitation is first aggregated to 10deltaX before performing
; the test, 10deltax is used as a scale, where the model actually resolves the
; features

; Diffferent results might be expected if the native resolution is used for comparison
; due to the chaotic nature of the atmospheric model

; Also, 40 grid points are removed from the periphery to remove the affect of boundary
;
; 26.03.2012 Created 
; 20.04.2012 Added Flags readcos, readclm, readpfl
;======================================================
; Prabhakar Shrestha
; Transregional Collaborative Research Centre 32 (TR32)
; Z4: Model Development and Maintenance Services
; Meteorological Institute
; University of Bonn
; Auf dem Hügel 20
; D-53121 Bonn
; E-mail: pshrestha@uni-bonn.de
;=======================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"


function extract_rain(diri_cos,nbdl,varname)
 begin
   fils_cos  = systemfunc("ls " + diri_cos + "lf*")
   nf_cos   = dimsizes(fils_cos)

  ;Read specified variables
  nc_cos1 = addfile(fils_cos(0),"r")
  nc_cos2 = addfile(fils_cos(nf_cos-1),"r")
  time1   = nc_cos1->time
  time2   = nc_cos2->time

  ;print(time1 + " " + time2)
 
  temp    = nc_cos1->$varname$
  ndim    = dimsizes(temp)

  var  = (nc_cos2->$varname$(0,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)) \
          - (nc_cos1->$varname$(0,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1))

  title = "Avg. rainfall " + avg(var) + "  " + stddev(var)
  print(" " + diri_cos)
  print(" " + title)
  return(var)
end

function aggregate(var2d,scale)
begin
  ndim = dimsizes(var2d)
  ny   = ndim(0)/scale
  nx   = ndim(1)/scale
  varnew = new((/ny,nx/),typeof(var2d))

  do iy = 0,ny-1
  sty = iy*scale
  ety = sty + scale-1
  do ix = 0,nx-1
     stx = ix*scale
     etx = stx + scale -1
     varnew(iy,ix) = avg(var2d(sty:ety,stx:etx))
  end do
  end do
  print("Avg. rainfall " + avg(varnew) + " " + stddev(varnew))
  return(varnew)
end 


begin

;-----------------------------------------------------------------------------------------
; USER SETTINGS
;-----------------------------------------------------------------------------------------

 ;;MODEL OUTPUT
 diri        = "/daten04/pshrestha/shared-data/bonnRadar_20141115/"
 case        = "HET2683"
 varname     = "TOT_PREC" ;
 bdr         = 3 + 37  ;exclude boundaries and relaxation zone


 event1      = "graupelhail2test4_satad" 
 event2      = "graupelhail2test4_satad_cpsUB1" 
 event2      = "graupelhail2test4_satad_Blahak1" 
 ;event2 =   "graupelhail2test4_satad_cps2"

 diri_cos1   = diri + "/" + case + "/" + event1 + "/cosout/"
 diri_cos2   = diri + "/" + case + "/" + event2 + "/cosout/"


 vscal = 10
 var2d =  extract_rain(diri_cos1,bdr,varname)
 var1  =  aggregate(var2d,vscal)

 var2d =  extract_rain(diri_cos2,bdr,varname)
 var2  = aggregate(var2d,vscal)

 y1    = ndtooned(var1)
 y2    = ndtooned(var2)

 system("rm temp.nc")
 fout = addfile("temp.nc","c")
 fout->var2d = var2d
 fout->var2 = var2
; t-test
 aveX    = avg (y1)        
 aveY    = avg (y2)         
 varX    = variance (y1)     
 varY    = variance (y2)     
 sX      = dimsizes (y1)    
 sY      = dimsizes (y2)    
 iflag   = True                                   ; population variance similar
 tval    = True                                   ; p-value
 prob = ttest(aveX,varX,sX, aveY,varY,sY, iflag, tval) 
 print(" " + prob)

; Plot
 wks = gsn_open_wks("x11","histo")

 res = True
 res@tiMainString = " " 
 plot = gsn_histogram(wks,y1,res)

end
