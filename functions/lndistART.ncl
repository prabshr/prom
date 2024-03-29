undef("lndistART")
;Estimate median diameter of lognormal distribution for each mode
;Estimate aerosol size distribution of select mode
;Estimates of ASD is provide for entire column at select location
; iy and ix based on entire domain, and favg only applies for fout = 0

;12 modes in COSMO-ART
;4 - nucleation/accumulation (pure and mixed)
;6 - small, medium and biggest (crustal and sea salt)
;2 - soot, coarse particles
;Each mode is approximated by a 3-parameter lognormal distribution 
;with constant chemical composition.
;N_i    - number concentration of mode i [m-3]
;sig_i  - geometric standard devaition of the mode i [-]
;Dpd_i  - median diameter of particles in model i [m]
;While N_i and sig_i are available directly from the model, Dpd_it needs
;to be estimated.
;The kth moment of the distribution can be solved as:
; M_k = N*(Dpd_i)^k [exp((k^2)/2 * (log(sig_i))^2.]
;Ackerman et al. (1998), Whitby 1978;

;--------------------
;art_specis_data includes the difference species of each mode
;eg., mode=NUCm,
;art_aerosol_const has the prescribed parameters of each mode
;--------------------
;Important note:
;This script for the first time step should produce median diameter equivalent to initial diameter 
;prescribed for each modes, except for dust modes, where number concentrations are also initialized !

;nbdl for favg, boundary cropped

function lndistART(fil_cos,opt)

begin
  smode = opt@smode
  favg  = opt@favg   
  iy    = opt@iy
  ix    = opt@ix
  oup   = opt@fout
  nbdl  = opt@nbdl

  if (smode.lt.1 .or. smode.gt.12) then
    print("Mode should be between 1..12")
    exit
  end if

  nc_cos   = addfile(fil_cos,"r")
  time     = nc_cos->time

  utc_date = cd_calendar(time, 0)
 ;
 ; Store return information into more meaningful variables.
 ;
   year   = tointeger(utc_date(:,0))    ; Convert to integer for
   month  = tointeger(utc_date(:,1))    ; use sprinti 
   day    = tointeger(utc_date(:,2))
   hour   = tointeger(utc_date(:,3))
   minute = tointeger(utc_date(:,4))
   second = utc_date(:,5)
 ;
 ; Write out strings in the format "hhZ dd mmm yyyy".
 ;

   date_str =  sprinti("%0.4i", year) + "/" + sprinti("%0.2i", month) +"/"+sprinti("%0.2i", day) + \
               "  " + sprinti("%0.2i", hour) + sprinti("%0.2i", minute) + " UTC"
   print("Processing ART data " + date_str)
;Do not touch below
;---------------------------------------------------------------------------------
  if (favg) then
    print("Processing " + fil_cos + " for domain average")
  else
    print("Processing " + fil_cos + " at " + iy + " " + ix)
  end if

  ;art_aerosol_const.f90 parameters
  ;nucl,acc,nucl-mixed,acc-mixed,soot,coarse,(small,medi,big) dust, 
  mname = (/"nucleation mode","accumulation mode","nucleation mode (mixed)", \
            "accumulation mode (mixed)", "soot mode", "coarse mode", \
            "small dust mode", "medium dust mode", "big dust mode", \
            "small sea-salt mode", "medium sea-salt mode", "big sea-salt mode" /)

  ;Component densities (kg/m3)
  rhoso4  = 1.8E3;
  rhonh4  = 1.8E3
  rhono3  = 1.8E3
  rhoorg  = 1.4E3
  rhosoil = 2.6E3
  rhoseas = 2.2E3
  rhoanth = 2.2E3
  rhosoot = 1.7E3
  rhowater = 1.0E3

  ;Conversion from ug/m3 to m3/m3
  PI  = get_pi("float" ) 
  cfac = 1.E-9 * 6.0/PI

  ;Min. diameter
  dpd_min = 1.E-9 ; m
  ;
  so4fac  = cfac/rhoso4
  nh4fac  = cfac/rhonh4
  no3fac  = cfac/rhono3
  orgfac  = cfac/rhoorg
  soilfac = cfac/rhosoil
  seasfac = cfac/rhoseas
  anthfac = cfac/rhoanth
  sootfac = cfac/rhosoot
  waterfac = cfac/rhowater
 
  ;Geomeric standard deviations of lognormal modes
  sig = (/1.7, 2.0, 1.7, 2.0, 1.4, 2.5, 1.7, 1.6, 1.5, 1.9, 2.0, 1.7/)

  imode = smode-1
  sig_i = sig(imode)

  ;Read vertical co-ordinate
  hhl = nc_cos->vcoord
  nz  = dimsizes(hhl)-1
  mhgt = new(nz,"float")
  mhgt = (/0.5*(hhl(0:nz-1) + hhl(1:nz))/)   ;agl

  if (smode.eq.1) then
    print("Processing " + mname(imode))
    N_i = nc_cos->VNU0

    M_i = so4fac*(nc_cos->VSO4I) + nh4fac*(nc_cos->VNH4I) + no3fac*(nc_cos->VNO3I)           + \  ;secondary inorganic
          orgfac*((nc_cos->VORG1I)+(nc_cos->VORG10I)+(nc_cos->VORG100I)+(nc_cos->VORG1000I)) + \  ;SOA
          anthfac*(nc_cos->VP25AI) + \                                                            ;primary inorganic
          orgfac*(nc_cos->VORGPAI) + \                                                            ;primary organic
          sootfac*(nc_cos->VECI)  +  \                                                            ;soot
          waterfac*(nc_cos->VH2OI)
  end if

  if (smode.eq.2) then
    print("Processing " + mname(imode))
    N_i = nc_cos->VAC0

    M_i = so4fac*(nc_cos->VSO4J) + nh4fac*(nc_cos->VNH4J) + no3fac*(nc_cos->VNO3J)           + \  ;secondary inorganic
          orgfac*((nc_cos->VORG1J)+(nc_cos->VORG10J)+(nc_cos->VORG100J)+(nc_cos->VORG1000J)) + \  ;SOA
          anthfac*(nc_cos->VP25AJ) + \                                                            ;primary inorganic
          orgfac*(nc_cos->VORGPAJ) + \                                                            ;primary organic
          sootfac*(nc_cos->VECJ)   + \                                                            ;soot
          waterfac*(nc_cos->VH2OJ)
  end if

  if (smode.eq.3) then
    print("Processing " + mname(imode))
    N_i = nc_cos->VNU0m

    M_i = so4fac*(nc_cos->VSO4Im) + nh4fac*(nc_cos->VNH4Im) + no3fac*(nc_cos->VNO3Im)           + \  ;secondary inorganic
          orgfac*((nc_cos->VORG1Im)+(nc_cos->VORG10Im)+(nc_cos->VORG100Im)+(nc_cos->VORG1000Im)) + \  ;SOA
          sootfac*(nc_cos->VSOOTI)   +  \                                                           ;soot
          orgfac*(nc_cos->VORGPAIm)  +  \                                                           ;primary organic
          waterfac*(nc_cos->VH2OIm)
  end if

  if (smode.eq.4) then
    print("Processing " + mname(imode))
    N_i = nc_cos->VAC0m

    M_i = so4fac*(nc_cos->VSO4Jm) + nh4fac*(nc_cos->VNH4Jm) + no3fac*(nc_cos->VNO3Jm)           + \  ;secondary inorganic
          orgfac*((nc_cos->VORG1Jm)+(nc_cos->VORG10Jm)+(nc_cos->VORG100Jm)+(nc_cos->VORG1000Jm)) + \  ;SOA
          sootfac*(nc_cos->VSOOTJ)   +  \                                                           ;soot
          orgfac*(nc_cos->VORGPAJm)  +  \                                                            ;primary organic
          waterfac*(nc_cos->VH2OJm)
  end if

  if (smode.eq.5) then
    print("Processing " + mname(imode))
    N_i = nc_cos->VSOOT0

    M_i = sootfac*(nc_cos->VSOOT) + \
          waterfac*(nc_cos->h2o_k) + so4fac*(nc_cos->vso4_k) + nh4fac*(nc_cos->vnh4_k) + no3fac*(nc_cos->vno3_k)  
  end if

  if (smode.eq.6) then
    print("Processing " + mname(imode))
    N_i = nc_cos->VCORN

    M_i = anthfac*(nc_cos->VANTHA) 
  end if

  if (smode.ge.7 .and. smode.le.9) then
    print("Processing " + mname(imode))
    ish = 7
    varnam_m = (/"VSOILA","VSOILB","VSOILC"/)
    varnam_N = (/"VSOILA0","VSOILB0","VSOILC0"/) 
    N_i = nc_cos->$(varnam_N(smode-ish))$

    M_i = soilfac*(nc_cos->$(varnam_m(smode-ish))$)
  end if

  if (smode.ge.10 .and. smode.le.12) then
    print("Processing " + mname(imode))
    ish = 10
    varnam_m1 = (/"VSEASA","VSEASB","VSEASC"/)
    varnam_m2 = (/"VSO4SEASA","VSO4SEASB","VSO4SEASC"/)
    varnam_m3 = (/"VH2OSEASA","VH2OSEASB","VH2OSEASC"/)
    varnam_N = (/"VSEASA0","VSEASB0","VSEASC0"/)
    N_i = nc_cos->$(varnam_N(smode-ish))$

    M_i = seasfac*(nc_cos->$(varnam_m1(smode-ish))$) + so4fac*(nc_cos->$(varnam_m2(smode-ish))$) + \
          waterfac*(nc_cos->$(varnam_m3(smode-ish))$)
  end if

  ;Compute median diameter of the mode
  Dpd_i    = (M_i/(exp(4.5*(log(sig_i))^2.)*N_i))^(1./3.)     ;[m]
  Dpd_i@long_name = "Median diameter"
  Dpd_i@units     = "m"
  Dpd_i@mode      = smode
  Dpd_i@mname     = mname(imode)

 ;Cutoff
  ;Dpd_i=Dpd_i>dpd_min
  print("Mode: " + smode + " Median Diameter: " + 1E9*avg(Dpd_i) + " nm" + " Ni " + avg(N_i)/1.E6 + "cm-3 Mi " + avg(M_i))


  if (oup.eq.0) then     ;Output ASD

  ;Compute lognormal aerosol size distribution
  ;First create size bins with constant dlnDp
  nchan = 1./32.

  p1 = 1.    ;nm
  nDp  = 300
  lp   = new(nDp,"float")

  do ic = 0,nDp-1
    if (ic.eq.0) then
       lp(ic) = log(p1)
    else
       lp(ic) = lp(ic-1) + nchan
    end if
  end do
  Dp = exp(lp)
  ;---------------------------------------------------

  rnlDp   = new((/2,nDp,nz/),"float")
  rnlDp!1 = "Dp"
  rnlDp!2 = "Z"
  rnlDp&Dp = Dp
  rnlDp&Z  = mhgt
  rnlDp@mode      = smode
  rnlDp@mname     = mname(imode)
  rnlDp@time      = date_str
  rnlDp@dlnDp     = nchan

  if (favg) then
    Ncm3   = N_i(0,:,:,:)/1.E6    ;cm-3
    Dpd    = Dpd_i(0,:,:,:)*1.E9  ;nm
    ndim   = dimsizes(Dpd)
    nlDp3d   = new((/nDp,nz,ndim(1),ndim(2)/),"float")
  else
    Ncm3   = N_i(0,:,iy,ix)/1.E6    ;cm-3
    Dpd    = Dpd_i(0,:,iy,ix)*1.E9  ;nm
  end if

  temp1  = Ncm3/(sqrt(2.*PI)*log(sig_i))
  temp3  = 2.*(log(sig_i))^2.0 

  do ib = 0, nDp-1
    temp2  = -1.*(log(Dp(ib)/Dpd))^2.
    if (favg) then
      nlDp3d(ib,:,:,:)  = temp1*exp(temp2/temp3) ; cm-3
    else
      rnlDp(0,ib,:)  = temp1*exp(temp2/temp3) ; cm-3
    end if
  end do

  if (favg) then
    do ib = 0, nDp-1
    print("Averaging for mode :" + ib)
    do iz = 0, nz-1
       if (iy.ne.-99 .and. ix.ne.-99) then
          print("Averaging at " + iy + " " + ix + " with 9x9 box")
          rnlDp(0,ib,iz) = avg(ndtooned(nlDp3d(ib,iz,iy-4:iy+4,ix-4:ix+4)))
          rnlDp(1,ib,iz) = stddev(ndtooned(nlDp3d(ib,iz,iy-4:iy+4,ix-4:ix+4)))
       else  
          rnlDp(0,ib,iz) = avg(ndtooned(nlDp3d(ib,iz,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1))) 
          rnlDp(1,ib,iz) = stddev(ndtooned(nlDp3d(ib,iz,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1))) 
       end if
    end do
    end do
  end if
 
  ;Reverse COSMO vertical levels
  rnlDp = rnlDp(:,:,::-1)
  rnlDp@units = "cm~S~-3~N~"
  return(rnlDp)

  else if (oup.eq.1) then
    ndim   = dimsizes(N_i(0,:,:,:))
    Ncm3   = N_i(0,::-1,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)/1.E6    ;cm-3
    Ncm3@time = date_str
    Ncm3@units = "cm~S~-3~N~"
    return(Ncm3)   ;cm-3
  else if (oup.eq.2) then
    ndim   = dimsizes(M_i(0,:,:,:))
    Mkg    = M_i(0,::-1,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)*1000.   ;m3/m3
    Mkg@time = date_str
    Mkg@units = "m~S~3~N~/m~S~3~N~"
    return(Mkg)
  else if (oup.eq.3) then
    ndim   = dimsizes(Dpd_i(0,:,:,:))
    Dpd    = Dpd_i(0,::-1,nbdl:ndim(1)-nbdl-1,nbdl:ndim(2)-nbdl-1)*1.E9  ;nm
    Dpd@time = date_str
    Dpd@units = "nm"
    return(Dpd)    ;nm
  else
    print("Error ... check opt@fout=0-3")
  end if
  end if
  end if
  end if
end
