function [DX, DX2] = fdder1(Nx, dx, BC)

%% DX SECTION
DX = sparse(Nx, Nx); 
DX = spdiags(+ones(Nx, 1), 1, DX);
DX = spdiags(-ones(Nx,1), -1, DX); 

%% DX2 SECTION 
DX2 = sparse(Nx, Nx);
DX2 = spdiags(-2*ones(Nx,1), 0, DX2); 
DX2 = spdiags(+ones(Nx, 1), 1, DX2);
DX2 = spdiags(+ones(Nx,1), -1, DX2); 

%% Boundary conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERIODIC BOUNDARY CONDITION
% | -2   1 |
% | 1   -2 | 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   p1 = BC(1); 
   p2 = BC(2); 
if p1 == 1 
   DX2(1, Nx) = p1;
   DX2(Nx, 1) = p2; 
   
elseif p1 == 0   %% Dirichlet BOUNDARY CONDITION
    DX2(1, Nx) =0; 
    DX2(Nx, 1) =0; 
   %% NEUMANN BOUNDARY CONDITION 
elseif p1 == -2
    DX(1, 1) = 2; 
    DX(1, 2) = -2; 
    DX(Nx, Nx-1) = -2; 
    DX(Nx, Nx) = 2;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DX = 1/2/dx*DX; 
DX2 = (1/dx^2)*DX2;
