%% codes for coaxail

%%%% CODES WERE DEVELOPED BY KEN & IBRA


%% INITIALIZATION 
close all 
clc 
clear

%% DASHBOARD

% GRID SIZE
Nx = 512;
Ny = 512;

%% constant values
meters      = 1;
seconds     = 1;
degrees     = pi/180;
F           = 1;
H           = 1;

% CONSTANTS
e0 = 8.85418782e-12 * F/meters;
u0 = 1.25663706e-6 * H/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;



% TRANSMISSION LINE PARAMETERS
erin  = 1 * eye(2,2);     % Permittivity tensor of inside of Coax cable
erout = 1.5 * eye(2,2);     % Permittivity tensot of outside of Coax cable
r1    = 0.35;               % Radius of inner conductor
r2    = 2.5;                % Radius of dielectric core
r3    = 2.7;                % Thickness of cladding


% SPACER REGIONS
BUFF  = 3*r3;
Sx    = BUFF;
Sy    = BUFF;
%% 


% INITIAL GUESS AT RESOLUTION
dx = Sx/Nx;
dy = Sy/Ny;

% COMPUTE 2X GRID
Nx2 = 2*Nx;
dx2 = dx/2;
Ny2 = 2*Ny;
dy2 = dy/2;

%% this>>>>>>>
% GRID AXES
xa = [0:Nx-1]*dx; xa = xa - mean(xa);
ya = [0:Ny-1]*dy; ya = ya - mean(ya);

% 2x GRID AXES
xa2 = [0:Nx2-1]*dx2; xa2 = xa2 - mean(xa2);
ya2 = [0:Ny2-1]*dy2; ya2 = ya2 - mean(ya2);

xc = 0; 
yc = 0; 
rx = 1;
ry = 1;
% CREATE MESH
[Y,X] = meshgrid(ya,xa);
RSQ = (((X-xc)/rx).^2 + ((Y-yc)/ry).^2);

% CREATE 2X MESH
[Y2,X2] = meshgrid(ya2,xa2);
RSQ2 = (((X2-xc)/rx).^2 + ((Y2-yc)/ry).^2);


%% Formulation from slides 
%BUILD INNER CONDUCTOR
CIN = (RSQ<=r1^2);
% BUILD OUTER CONDUCTOR
COUT = (RSQ<r3^2 & RSQ>=r2^2);

er1 = erin(1,1); 
er2 = erout(2,2);
% BUILD DIELECTRIC
r23 = (r2 + r3)/2; %middle of outer conductor
ER2 = ones(Nx2,Ny2); %air
ER2 = ER2 + (er1 - 1)*(RSQ2<r23^2 & X2<0); %dielectric 1
ER2 = ER2 + (er2 - 1)*(RSQ2<r23^2 & X2>=0); %dielectric 2


% EXTRACT ERxx AND ERyy FROM ER2
Erxx = ER2(2:2:Nx2,1:2:Ny2);
Eryy = ER2(1:2:Nx2,2:2:Ny2);
% FORM DIAGONAL PERMITTIVITY MATRICES
ERxx = diag(sparse(Erxx(:)));
ERyy = diag(sparse(Eryy(:)));
% FORM PERMITTIVITY TENSOR
Z = sparse(Nx*Ny,Nx*Ny);
ER = [ ERxx , Z ; Z , ERyy ];


% CALL FUNCTION TO CONSTRUCT DERIVATIVE OPERATORS
NS = [Nx Ny]; %grid size
RES = [dx dy]; %grid resolution

% Construct derivative matrix operators
BC = [-2 -2];
kx = .3; 
ky =0; 
kinc = [kx ky];
[DVX,DVY,DEX,DEY] = yeeder(NS,RES,BC, kinc);

% inhomegeneous for calculation the capacitance of the coaxial ...
% transmission line 
L = [DEX DEY] * ER * [DVX ; DVY];


%% Build the homogenous one for Induction

Lh = [DEX DEY] * [DVX ; DVY];


% FORCE MATRIX
F = CIN | COUT;
F = diag(sparse(F(:)));
% FORCED POTENTIALS
vf = 1*CIN + 0*COUT;
% FORCE KNOWN POTENTIALS
M = Nx*Ny;
I = speye(M,M);
L = (I - F)*L + F;
Lh = (I - F)*Lh + F;
b = F*vf(:);



% COMPUTE POTENTIALS
v = L\b;
vh = Lh\b; % COMPUTE E FIELDS
e = - [ DVX ; DVY ] * v;
eh = - [ DVX ; DVY ] * vh;
% COMPUTE D FIELDS
d = ER*e;
dh = eh;



% DISTRIBUTED CAPACITANCE
V0 = 1;
C = d.'*e*(e0*dx*dy)./V0^2;
% DISTRIBUTED INDUCTANCE
Ch = dh.'*eh*(e0*dx*dy)./V0^2;
L = 1/(c0^2*Ch);
% CHARACTERISTIC IMPEDANCE
Z0 = sqrt(L/C);
% EFFECTIVE REFRACTIVE INDEX
neff = c0*sqrt(L*C);

% SHOW NUMERICAL PARAMETERS ON CONSOLE
disp(['C    = ' num2str(C/1e-12,'%3.5f') ' pF/m']);
disp(['L    = ' num2str(L/1e-09,'%3.5f') ' nH/m']);
disp(['Z0   = ' num2str(Z0) ' Ohms']);
disp(['nEff = ' num2str(neff)]);


% RESHAPE THE FUNCTIONS BACK TO A 2D GRID
Vh = reshape(vh,Nx,Ny);
V = reshape(v,Nx,Ny);

% Obtain fields
Ex = e(1:M);
Ey = e(M+1:2*M);
Ex = reshape(Ex,Nx,Ny);
Ey = reshape(Ey,Nx,Ny);

E = sqrt(abs(Ex).^2 + abs(Ey).^2);

%% Plotting section
% SET SPACING FOR QUIVER
sp = 18;
% VISUALIZE POTENTIAL AND FIELDS
hold on;
imagesc(xa,ya,V');

contour(xa, ya, COUT, 2);
colormap(jet);
set(gca,'FontSize',12,'FontWeight','bold');
colorbar;
axis equal tight;
title('Electric Potential V','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);

figure('Color','w');
imagesc(xa,ya,E'); 

caxis([0 0.75]);
set(gca,'FontSize',12,'FontWeight','bold');
colorbar
colormap(hot);
axis equal tight;
title('|E|','FontSize',14);
xlabel('x (mm)','FontSize',12);
ylabel('y (mm)','FontSize',12);
hold on;

% ADD QUIVER
[Y,X] = meshgrid(ya,xa);
quiver(X(1:sp:Nx,1:sp:Ny),Y(1:sp:Nx,1:sp:Ny),Ex(1:sp:Nx,1:sp:Ny),...
      Ey(1:sp:Nx,1:sp:Ny),'Color','w');
hold off;

