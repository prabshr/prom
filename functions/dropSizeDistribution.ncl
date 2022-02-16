function dropSizeDistribution(qx,Nx,opt)
;mass density qx [gm-3], number density Nx [cm-3]
;opt contains shape parameters and power law parameters

begin

  ;Modified Gamma Distribution (MGD)
  ;f(x) = No*(x^mu)*exp(-lam*(x^nu))
  ; x is particle mass and mu,nu are shape parameters

  ag = opt@ag
  bg = opt@bg
  
;Reversed in COSMO-SB2M
  mu = opt@mu
  nu = opt@nu
  ;

;Convert units
  qxu = qx * 1.E-3  ;kg/m3
  Nxu = Nx * 1.E6   ;m-3 
  x   = qxu/Nxu      ;mass [kg] 
  print("Mass x = " + x + " kg")

  lam = ( gamma((mu+1.d0)/nu)/ gamma((mu+2.d0)/nu) * x)^(-nu)
  No  = nu*Nxu*lam^((mu+1.d0)/nu) / gamma((mu+1.d0)/nu) 

  ;particle diameter
  D     = (ispan(1,500,1)/10.)*1.E-6   ;m
  ;modified gamma particle size distribution (PSD)
  nud   = nu/bg
  mud   = ((mu+1.d0)/bg) - 1.d0
  Nod   = No/(bg*ag^((mu+1.)/bg))
  lamd  = lam/(ag^(nud))

  fD    = Nod * (D^(mud)) * exp(-lamd*(D^(nud))) 

  ;Change units
  D     = D *1.E6  ;um
  fD    = fD/1.E6/1.E6  ;cm-3

  D@units = "~F5~m~F21~m"
  fD@units = "cm~S~-3~N~ ~F5~m~F21~m~S~-1~N~"
  fD!0  = "D"
  fD&D  = (/D/)

  return(fD)
end
