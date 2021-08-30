;Edge detection and storm tracking algorithm
;Useful for CAF estimation for selected storm...
;This algorithm works in 2 step iterative fashion
;1) Run the script with opt  = False  , it will create an output with edges
;2) Select the co-ordinates of storm at different times
;3) Enter the co-ordinates with time , opt@jd,opt@id
;4) Run the script again with opt = True 

;Functions:
; convol
; trackedge
; scanfill
;The storm has to be selected manually

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

function convol(data)
;Roberts cross convolution
begin 
  print("Roberts cross convolution to find edges ....")
  ;Kernel
  kx = (/(/1,0/),(/0,-1/)/)
  ky = (/(/0,1/),(/-1,0/)/)
  ;
  ndim = dimsizes(data)
  ny   = ndim(0)
  nx   = ndim(1)
  yy   = data^(0.5)
  GradI = new((/2,ny,nx/),"float")
  do j = 0, ny-2
    jj = j+1
    do i = 0, nx-2
      ii = i + 1
      Gx   = (kx(0,0)*yy(j,i) + kx(1,1)*yy(jj,ii)) \
              + (kx(1,0)*yy(jj,i) + kx(0,1)*yy(j,ii))
      Gy   = (ky(0,0)*yy(j,i) + ky(1,1)*yy(jj,ii)) \ 
              + (ky(1,0)*yy(jj,i) + ky(0,1)*yy(j,ii))
      GradI(0,j,i) = (Gx^2 + Gy^2)^(0.5)
      if (Gx.ne.0) then
        GradI(1,j,i) = atan(Gy/Gx) * 180./3.14159
      end if
    end do
  end do
  GradI = where(GradI.eq.0,0,1)
  return(GradI(0,:,:))
end

function trackedge(jd,id,data)
begin
  print("Edge coordintest tracking for JJ=" + jd + " II=" + id)
  ;Dimensions of data
  ndim = dimsizes(data)
  ny   = ndim(0)
  nx   = ndim(1)
  fdebug = True 

  ;Pad zeros to remove boundary affect
  udata = new((/ny+1,nx+1/),typeof(data))
  udata = 0
  udata(1:ny,1:nx) = data
  ;Copy of data for visited
  ucdata = udata

  ;Find the rightmost edge, starting location
  temp = ind(udata(jd,id:nx-1).eq.1)
  ii0  = id + temp(0)
  jj0  = jd
  print("Start jj ii :" + jj0 +"  " + ii0)

  ;Initialize centre of serach kernel
  jj = 0
  ii = 0
  ctr = -1  ; search counter

  ;Search kernel
  ; |  *  |
  ; *  |  *
  ; |  *  |
  drc     = (/(/1,0/),(/-1,0/),(/0,-1/),(/0,1/)/)

  do while (ctr.lt.5000)
    ctr = ctr + 1
    if (fdebug) then
      print("Count :" + ctr + " jj ii :" + jj +"  " + ii)
    end if
    ;4x4
    if (ctr.eq.0) then
      jjs     = drc(:,0) + jj0
      iis     = drc(:,1) + ii0
    else
      jjs     = drc(:,0) + jj
      iis     = drc(:,1) + ii
    end if
    jjs(0)    = min((/jjs(0),ny-1/))
    jjs(1)    = max((/0, jjs(1)/))
    iis(2)    = max((/0, iis(2)/))
    iis(3)    = min((/iis(3),nx-1/)) 
   
    datac = (/udata(jjs(0),iis(0)),udata(jjs(1),iis(1)),udata(jjs(2),iis(2)),udata(jjs(3),iis(3))/)

    ;Debug print
    ;print("x " + datac(0) +  " x ")
    ;print(datac(2) + "  x " +  datac(3))
    ;print("x " + datac(1) +  " x ")

    ;Edge path
    cind = ind(datac.eq.1)
    if (ctr.eq.0) then
      iloc1 = cind(0) 
    else
      ncdim = dimsizes(cind)
      if (ncdim.eq.2) then
        dp = cind(0)
        if ((jp(ctr-1).eq.jjs(dp)) .and. (ip(ctr-1).eq.iis(dp))) then
          iloc1 = cind(1) 
        else
          iloc1 = dp 
        end if
        if (fdebug) then
          print(ncdim + " " + cind(0) + " " + cind(1) + " " + iloc1 + " jp/ip-1: " + jp(ctr-1) + " " + ip(ctr-1))
          print(jjs(dp) + " " + iis(dp) + " jjs/iis 1: " + jjs(cind(1)) + " " + iis(cind(1)))
        end if
      else if (ncdim.gt.2) then
        if (fdebug) then
          print("Entering thick edge ...")
        end if
        do ic = 0, ncdim-1
           dp = cind(ic)
           if (jp(ctr-1).ne.jjs(dp) .and. ip(ctr-1).ne.iis(dp) .and. ucdata(jjs(dp)+drc(dp,0),iis(dp)+drc(dp,1)).eq.1) then
             iloc1 = dp
           end if

;           if (dp.eq.0 .and. jp(ctr-1).ne.jjs(dp) .and. ip(ctr-1).ne.iis(dp) .and. ucdata(jjs(dp)+1,iis(dp)).eq.1) then
;             iloc1 = dp
;           else if (dp.eq.1 .and. jp(ctr-1).ne.jjs(dp) .and. ip(ctr-1).ne.iis(dp) .and. ucdata(jjs(dp)-1,iis(dp)).eq.1) then
;             iloc1 = dp 
;           else if (dp.eq.2 .and. jp(ctr-1).ne.jjs(dp) .and. ip(ctr-1).ne.iis(dp) .and. ucdata(jjs(dp),iis(dp)-1).eq.1) then
;             iloc1 = dp
;           else if (dp.eq.3 .and. jp(ctr-1).ne.jjs(dp) .and. ip(ctr-1).ne.iis(dp) .and. ucdata(jjs(dp),iis(dp)+1).eq.1) then
;             iloc1 = dp
;           end if
 ;          end if
 ;          end if
 ;          end if
           if (fdebug) then
             ;print("ic = " + ic + " " + dp + "--------------------")
             ;print(ncdim + " " + cind(0) + " " + cind(1) + " " + iloc1 + " jp/ip-1: " + jp(ctr-1) + " " + ip(ctr-1))
             ;print(jjs(dp) + " " + iis(dp) + " jjs/iis 1: " )
           end if
         end do
         sse  = 1    ;reset search extent for iloc1 =-1
      else if (all(ismissing(cind))) then
         iloc1 = -1         
      end if
      end if
      end if
    end if
    delete(cind)

    ;Append co-ordinates
    if (ctr.eq.0) then
      ;Initialize polyline pair
      jp  = jj0
      ip  = ii0
    else
      temp0 =  array_append_record(jp, jj, 0)
      temp1 =  array_append_record(ip, ii, 0)
      delete([/jp, ip/])
      jp = temp0
      ip = temp1
      if (jp(ctr).eq.jj0 .and. ip(ctr).eq.ii0) then
        print("Edge co-ordinate extracted successfully ...")
        print("Total perimeter of strom = " + ctr + " dx")
        break
      end if
      delete([/temp0, temp1/])
    end if

    if (iloc1.ne.-1) then
      ;Visited points
      ucdata(jjs(iloc1),iis(iloc1)) = 2

      jj = (/jjs(iloc1)/)
      ii = (/iis(iloc1)/)
    else 
      ;search increases gradually from min 1 to 5 
      sse    = min((/sse + 1, 5/))   ;max 5x5 box
      sbox   = udata(jj-sse:jj+sse,ii-sse:ii+sse)
      sbox1d = ndtooned(sbox)
      sboxsize = dimsizes(sbox)
      sind   = ind_resolve(ind(sbox1d.eq.1),sboxsize)
       low  = -dimsizes(sind(:,0))*1.0 
       high = dimsizes(sind(:,0))*1.0 
       con  = (high - low) / 32766.0  
       n    = abs(tointeger(low + con * rand()))
       n    = min((/n, dimsizes(sind(:,0))-1/))
       print(sse + " " +  dimsizes(sind(:,0)) + " rand " + n)
     
      jj = tointeger(jj+sind(n,0))
      ii = tointeger(ii+sind(n,1))
      delete([/sbox,sbox1d,sind/])
      ;Visited points
      ;ucdata(jj,ii) = 2
    end if
  end do
 
  ;Concat array for edge co-ordinates
  nji = dimsizes(jp)
  ji  = new((/nji,2/),"integer")
  ji(:,0) = jp -1   ;remove padding
  ji(:,1) = ip -1   ;remove padding

  do ix = 0, nji-1
    udata(ji(ix,0)+1,ji(ix,1)+1) =   2
  end do
  udata = where(udata.eq.2,udata,0)
  udata@jc = ji(:,0)
  udata@ic = ji(:,1)
  return(udata(1:ny,1:nx))
end

function scanfill(data)
begin
  print("Masking selected storm ...")
  print("--------------------------------------------")
  jj    = data@jc
  ii    = data@ic
  udata = (/data/)*0

  ;Crop data surrounding the storm region
  cind  = (/min(jj),max(jj),min(ii),max(ii)/)
  ny    = cind(1) - cind(0) + 1
  nx    = cind(3) - cind(2) + 1

  ;Generate co-ordinates for data
  jj_1d = ispan(cind(0),cind(1),1)
  ii_1d = ispan(cind(2),cind(3),1)

  ;
  jj_2d = conform_dims((/ny,nx/),jj_1d,0) 
  ii_2d = conform_dims((/ny,nx/),ii_1d,1)

  inout = gc_inout(jj_2d, ii_2d, jj, ii)

  gmsk  = where(inout.eq.True,1,0)

  udata(cind(0):cind(1),cind(2):cind(3)) = gmsk 
  return(udata)
end

function bucketfill(jd,id,data,iter,opt)
begin
  directions = tointeger((/(/1,0/),(/-1,0/),(/0,1/),(/0,-1/)/))
  nd         = dimsizes(directions)
  ndim       = dimsizes(data)
  odata      = data
  copt       = opt
  if (jd.ge.ndim(0) .or. id.ge.ndim(1) .or. jd.lt.0 .or. id.lt.0 \
     .or. data(jd,id).eq.1 .or. data(jd,id).eq.2) then
    iter2 = iter + 1
    if (iter2.eq.4) then
      iter2 = 0
    end if
    if (opt@bdr) then
      j = tointeger(copt@jd + directions(iter2,0))
      i = tointeger(copt@id + directions(iter2,1))
    else
      j = tointeger(jd + directions(iter2,0))
      i = tointeger(id + directions(iter2,1))
    end if
    copt@bdr = True 
    copt@jd  = jd
    copt@id  = id
    print("False " + j + " " + i + " " + iter2)
    return(bucketfill(j,i,odata,iter2,copt))
  else
    copt@bdr = False
    jd2 = jd
    id2 = id
    iter2 = iter
    odata(jd,id) = 2 
    j = tointeger(jd2 + directions(iter2,0))
    i = tointeger(id2 + directions(iter2,1))
    print("True " + j + " " + i + " " + iter2)

    return(bucketfill(j,i,odata,iter2,copt))
  end if
end

begin
  nc_radar = addfile("debug_BoXPol.nc","r")
  fopt =  True ;False 
  jd  = (/ 72, 95, 94,102,112,126,165,140,166,176,177,180,160,186,173,182, \
              184,180,208,217,231,232,232,186,203,203,215,216,231/)
  id  = (/ 33, 38, 49, 56, 62, 73, 89,106,100,109,120,131,139,149,167,190,  \
              196,208,201,205,223,229,248,283,288,303,314,325,331/)

  ;;; opt@bdr = False

;
  if (fopt) then
    nc_fil = addfile("data1.nc","r")
    zh     = nc_fil->Band1
    ndim     = dimsizes(zh)
    st = 0
    et = 0 ;dimsizes(jd)-1 
    do it = st,et
      val  = trackedge(jd(it),id(it),zh(it,:,:))
      val  = scanfill(val)
      zh(it,:,:) = (/val/)
    end do
  else
    zh       = nc_radar->Band1
    ndim     = dimsizes(zh)
    delete(zh@_FillValue)
    rnd = random_uniform(0.1,0.5,1)

    zh = where(zh.ne.999,1,rnd)
    st = 0 
    et = ndim(0)-1 ;dimsizes(jd) -1 ; ndim(0)-1 

    do it  = st, et  
      print("Processing " + it + " of " + ndim(0))
      val = convol(zh(it,:,:))
      zh(it,:,:) = (/val/)
      delete(val)
    end do
  end if

  system("rm temp2.nc")
  fout = addfile("temp2.nc","c")
  fout->Band1 = zh(st:et,:,:) 
end
