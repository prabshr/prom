load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; N colors and N-1 ticks, the ticks are the boundaries of the discrete colors, except the lower and upper edges!
; except first 2 background and foreground colors
function load_cfadColors(opt) 
begin
 
 ;Contour Color Map
 colors = (/"white","black","white","gray30","gray40","gray50","gray61","gray75", \
          "gray80","gray90","paleturquoise4","paleturquoise3","orchid4","orchid3","orchid2","orchid1","orchid", \
         "royalblue4","darkturquoise","yellow"/)

 ;Bin for CFADS
 ;Zh, Zdr, Kdp, Rhv
 binw   = (/ 5.0, 0.25,   0.25,   0.05/)
 minv   = (/-10.,  -1.,  -1.0,   0.7 /)
 maxv   = (/ 70.,   7.,   7.,    1.1 /)
 
 rgb = namedcolor2rgb(colors)
 ;write_matrix(rgb,"3f5.2",0)

 colors@binw = binw
 colors@minv = minv
 colors@maxv = maxv

;cn Levels for contour, 18 colors ~ 17 ticks 
 if (opt@YH95) then
   if (opt@vid.eq.0) then
     cnLevels                   = ispan(25,300,25)/10.0 ;2.5 % per dbZ per km
   else if (opt@vid.eq.1) then
     cnLevels                   = ispan(50,350,50)      ;50 % per dbZ per km
   end if
   end if
 else
   ;Lang2007
   cnLevels                   = (/0.01,0.1,1,2,3,4,5,6,7,10,15,20,25,30,40,50,75/)
 end if

 colors@cnLevels = cnLevels
 
 return(colors) 
end
