;Extract aerosol for model comparison with observations

begin

;STEP 1 User Settings
;-----------------------------------------------------------------------------------------

 pname = "EMPA2_PROM1";
 einst = (/5,6,7,8,9/)
 lnd   = "HET"
 itype = "ART"
 diri  = "/p/largedata/hbn33/shrestha1/"
 vid   = 3   ;Variable ID
 sid   = 1   ;Station ID

;Tracer
 varnames = (/"NO2","HCHO", "VSOILA0", "TAU"/)

 varname  = varnames(vid)

 ;Geo-location
 stname      = (/"BoXPol", "FJZ-JOYCE", "MAINZ"/)
 lonc        = (/7.071663,6.413056,8.300000/); N
 latc        = (/50.73052,50.908333,49.999000/) ;E


; Read Variables
;-----------------------------------------------------------------------
 ninst = dimsizes(einst)    ;ENS LOOP
 print(stname(sid) + " " + varnames(vid))
 do ic = 0, ninst-1 

 itype_gscp  = lnd + itype
 sfold       = "runE_" + einst(ic) + "_" + itype_gscp
 expname     = lnd + "_" + pname
 subevent    = expname + "/" + itype_gscp + "/" + sfold
 diri_cos    = diri + "/" + subevent + "/a_cosout/"
 diri_cosI   = diri + "/" + subevent + "/a_cosout/ivr/"


 fil_ivr    = diri_cosI + "lfff00000000c.nc"
 nc_cosI    = addfile(fil_ivr,"r")
 lon2d      = nc_cosI->lon
 lat2d      = nc_cosI->lat
 rotpole    = nc_cosI->rotated_pole
   pole_lon = rotpole@grid_north_pole_longitude
   pole_lat = rotpole@grid_north_pole_latitude

 topo       = nc_cosI->HSURF(0,:,:)
 hhl       = nc_cosI->HHL(0,::-1,:,:)

 ndimhhl   = dimsizes(hhl)
 nz        = ndimhhl(0)-1
 ny        = ndimhhl(1)
 nx        = ndimhhl(2)

 fils_cos  = systemfunc("ls " + diri_cos + "lf*")
 nf_cos    = dimsizes(fils_cos)-1


 if (ic.eq.0) then
   ;Station location
   nm     = getind_latlon2d (lat2d,lon2d, latc(sid), lonc(sid))
   iys    = nm(0,0)
   ixs    = nm(0,1)
   print("O: " + latc(sid) + " " +  lonc(sid) + " M: " + lat2d(iys,ixs) + " " + lon2d(iys,ixs))

   tracer = new((/ninst,nf_cos/),"float")
   time   = new(nf_cos,"double")
 end if

 print("Processing ensemble ..." + einst(ic))

 do fn = 0,nf_cos-1    ;BIG LOOP for animation

   nc_cos    = addfile(fils_cos(fn),"r")
   time(fn)  = nc_cos->time

   if (varname.eq."TAU") then
     tracer3d = nc_cos->TAU_AERO_A(0,::-1,:,:) + \
              nc_cos->TAU_DUST_A(0,::-1,:,:) + \
                nc_cos->TAU_SEAS_A(0,::-1,:,:)
     tracersum  = dim_sum_n(tracer3d,0)
     tracer(ic,fn) = tracersum(iys,ixs)

     copy_VarAtts(tracer3d,tracer)
     tracer@long_name = "AOD at 550 nm"
     tracer@units     = "-"
   else
     tracer3d  = nc_cos->$varname$(0,::-1,:,:)
   end if
 end do ;BIG LOOP 
 end do ;ENSEMBLE LOOP

 ;Data output
 foutname = stname(sid) + "_" + varname + "_model.nc"
 system("rm " + foutname)
 fout = addfile(foutname,"c")
 fout->Band1 = tracer
 fout->time  = time

end
