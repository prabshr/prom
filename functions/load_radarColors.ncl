load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; N colors and N-1 ticks, the ticks are the boundaries of the discrete colors, except the lower and upper edges!
; except first 2 background and foreground colors
function load_radarColors(opt) 
begin
 
 colors = (/"white","black","cyan","deepskyblue2","blue", \
            "chartreuse", "chartreuse3", "chartreuse4", \
            "yellow", "yellow3", "orange", \
            "orangered","orangered3","orangered4", \
            "magenta", "mediumpurple"/) 
 zhh_levs = ispan(0,60,5)*1. -10.
 zdr_levs = (/-1.  ,-0.1 , 0.0,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.8,  1.0,  2.0,   3.0/)
 rho_levs = (/ 0.7 , 0.8 , 0.85, 0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.998/)
 kdp_levs = (/-0.50,-0.10, 0.0,  0.05, 0.10, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 2.00,  20.00/)

 rgb = namedcolor2rgb(colors)
 ;write_matrix(rgb,"3f5.2",0)

 colors@zhh_levs = zhh_levs
 colors@zdr_levs = zdr_levs
 colors@rho_levs = rho_levs
 colors@kdp_levs = kdp_levs

 return(colors) 
end
