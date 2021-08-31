undef("calcrho_cos_nclf")
function calcrho_cos_nclf(t,p,qv,qc,qrs)
begin
 ;p - pressure (Pa)
 ;qv - specific humidity (kg/kg)
 ;qc - specific cloud water content (kg/kg)
 ;qrs - specific precipitation loading (kg/kg)
 ;rho - air density (kg/m3) 

 Rd     = 287.05   ;Jkg-1K-1
 Rv     = 461.51   ;Jkg-1K-1
 rat    = Rv/Rd -1.0

 rho    = p/(Rd*t*(1 + rat*qv - qc - qrs))

 return(rho)
end

