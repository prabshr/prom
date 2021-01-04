undef("cosmo_time_nclf")
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
function cosmo_time_nclf(time,topt)
begin
 dtime = calendar_decode2(time,-5)
 year  = dtime(0,0)
 month = dtime(0,1)
 day   = dtime(0,2)
 hh    = dtime(0,3)
 mm    = dtime(0,4)
 ss    = dtime(0,5)

; fix for bug in calendar_decode2
 if (ss.eq.60) then
   mm = mm +1
 end if
;
 str_time  = sprinti("%0.4i",year) + "-" + sprinti("%0.2i",month) + \
              "-" + sprinti("%0.2i",day) + " " + sprinti("%0.2i",hh) + \
               sprinti("%0.2i",mm) + " UTC" ;sprinti("%0.2i",ss) +  " UTC"

 if (topt.eq.0) then
 str_time  = sprinti("%0.4i",year) + "-" + sprinti("%0.2i",month) + \
              "-" + sprinti("%0.2i",day) + " " + sprinti("%0.2i",hh) + \
               sprinti("%0.2i",mm) ;+ " UTC" ;sprinti("%0.2i",ss) +  " UTC"
 else
 str_time  = sprinti("%0.2i",hh) + \
               sprinti("%0.2i",mm) + " UTC"
 end if 

 return(str_time)
end

