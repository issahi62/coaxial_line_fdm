 function [DFX, DGX] = mult_var(Ns, RS)
 
 Nx = Ns(1); 
 Ny = Ns(2); 
 
 xa = RS(1); 
 xb = RS(2);
 
 dx = (xb-xa)/(Nx-1); 
 
 DFX = sparse(Nx, Ny); 
 DFX = spdiags(-1*ones(Nx, 1), 0, DFX); 
 DFX = spdiags(+ones(Nx, 1), 1, DFX); 
 DFX = DFX/dx; 
 
 DGX = -DFX'; 
 

 
