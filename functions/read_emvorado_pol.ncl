load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;Extracts data from EVORADARO-POL outputs
;opt@profile = 1 returns vertical profile at opt@iy and opt@ix
;Else returns the whole data
;opt@nf specifies the time snapshhot to process in the folder

function read_emvorado_pol(diri,opt)
begin

  fil  = systemfunc("ls " + diri + "/" + opt@var + "*.dat")
  print("Processing " + fil(opt@nf) + "...." )

  d1   = asciiread(fil(opt@nf), -1, "float")

  if (opt@var.ne."kdp") then
    iind = 2
  else
    iind = 0
  end if

  ny   = tointeger(d1(iind))
  nx   = tointeger(d1(iind+1))
  nz   = tointeger(d1(iind+2))
  data = onedtond(d1(iind+3:),(/nz,ny,nx/))
  data@_FillValue = -999.99

  print("EMVORADO-POL auxillary output " + opt@var + ": Data Size: " + ny + " " + nx + " " + nz)
  print("--------------")

  if (opt@var.eq."zh") then
    ;mask low value 
    data = mask(data,data.lt.10E-9, False)
    data@_longname = "Z~B~H~N~"
    data@units     = "[mm^6 m^−3]"
  end if

  if (opt@var.eq."zv") then
    ;mask low value 
    data = mask(data,data.lt.10E-9, False)
    data@_longname = "Z~B~V~N~"
    data@units     = "[mm^6 m^−3]"
  end if

  if (opt@var.eq."rhv") then
    ;mask low value 
    data = mask(data,data.lt.10E-9, False)
    data = 10*log10(data/1.)
    data@_longname = "~F8~r~B~~F10~rhv~E~"
    data@units     = "[mm^6 m^−3]"
  end if

  if (opt@var.eq."ihv") then
    ;mask low value 
    data = mask(data,data.lt.10E-9, False)
    data@_longname = "~F8~r~B~~F10~ihv~E~"
    data@units     = "[mm^6 m^−3]"
  end if

  if (opt@var.eq."kdp") then
    data = mask(data,data.eq.0, False)
    data@_longname = "K~B~DP~N~"
    data@units     = "[deg/m]"
  end if

  ;Reverse data
  data = data(::-1,:,:)
  if (opt@profile .eq. 1) then
    print("Returning vertical profile at " + opt@iy + " " + opt@ix)
    data1d = data(:,opt@iy,opt@ix)
    return(data1d)
  else
    return(data)
  end if
end
