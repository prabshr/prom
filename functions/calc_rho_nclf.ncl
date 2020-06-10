undef("calc_rho_nclf")
function calc_rho_nclf(t,p,qv,qrs)
begin
 ;p - pressure (Pa)
 ;qv - specific humidity (kg/kg)
 ;qc - specific cloud water content (kg/kg)
 ;qrs - specific precipitation loading (kg/kg)
 ;rho - air density (kg/m3) 

 Rd     = 287.05   ;Jkg-1K-1
 Rv     = 461.51   ;Jkg-1K-1
 rat    = Rv/Rd -1.0

 rho    = p/(Rd*t*(1 + rat*qv - qrs))

 return(rho)
end

