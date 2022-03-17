;Identify cell location based on initial input

begin
;--User Settings --------------------------------------------------
  pid       = 3            ;Experiment ID
  lnd       = "HET"
  itype     = "ART_SP4"
  ;inst    = 8

  ;First cell
  t  = 0   ;14 
  iy = 120; 128 
  ix = 150 ;157

  ;Search box size
  bs = 40
  pnames   = (/"PROM1","PROM2","PROM3","EMPA2_PROM1"/)

  itype_gscp = lnd + itype   
  fname    = pnames(pid) + "_runE_" + inst + "_" + itype_gscp
;------------------------------------------------------------------

  ;Read nc file
  fout = addfile(fname + ".nc","r")
  dev  = fout->Band2

  ;
  ndim = dimsizes(dev)
  nt   = ndim(0)
  ny   = ndim(1)
  nx   = ndim(2)

  yp   = new(nt,"integer")
  xp   = new(nt,"integer")
  yp   = -999
  xp   = -999

  ;Initialize
  yp(t) = iy
  xp(t) = ix

  do it = t+1, nt-1

    nys  = max((/0,iy-1/))
    nye  = min((/ny-1,iy+bs/))
    nxs  = max((/0,ix-1/))
    nxe  = min((/nx-1,ix+bs/))

    ;print("Processing ..." + it + "X " + nxs + " " + nxe + " Y " + nys + " " + nye)
    sbox = dev(it,nys:nye,nxs:nxe)

    ;Loop over Y to find location with overlapping of storm
    do iyy = nys-nys,nye-nys
      sind =  ind(sbox(iyy,:).eq.2)
      if (.not.all(ismissing(sind))) then
        ixx = sind(0)
        delete(sind)
        break
      else
        delete(sind)
      end if
    end do
    delete(sbox)
    ;Update location by adding offset
    iy = iyy + nys
    ix = ixx + nxs
    print("it: " + it + " iy: " + iy + " ix:  " + ix)
    yp(it) = iy
    xp(it) = ix
  end do


  ;Start writing the file
  csv_filename = "cell_" + fname + ".txt"
  system("rm " + csv_filename)
  hlist = [/"#JJ II"/]
  write_table(csv_filename, "w", hlist, "%s")
  alist = [/yp,xp/]
  write_table(csv_filename, "a", alist, "%i,%i")
end
