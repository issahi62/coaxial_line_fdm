% DEFINE BOUNDARY VALUES
xa = 0;
xb = 2;
fa = 2;
fb = 0.2;

%% BOUNDARY CONDITIONS
BC = [ 1; 1]; %% periodic BC 
%BC = [-2, -2]; %% Nuemann BC
%BC = [ 0, 0]; % Dirichlet BC


% GRID PARAMETERS
Nx = 10000;
dx = (xb - xa)/(Nx - 1);
x = linspace(xa, xb, Nx); 
% BUILD MATRIX OPERATORS
I = speye(Nx,Nx);

[DX,DX2] = fdder1(Nx,dx, BC);
% CALCULATE [A] AND [b]
A = DX2 + 5*DX + 6*I;
b = sparse(Nx,1);

% INCORPORATE BOUNDARY VALUES
A([1 Nx],:) = 0;
A(1,1) = 1;
A(Nx,Nx) = 1;
b([1 Nx]) = [ fa fb ];
% SOLVE PROBLEM
f = A\b;


%% plotting section 
% PLOT RESULT
subplot(121)
h = plot(x,f,'--b','LineWidth',2);
h2 = get(h,'Parent');
set(h2,'LineWidth',2,'FontSize',18);
xlabel('$x$','Interpreter','LaTex');
ylabel('$f(x)$','Interpreter','LaTex');
T = [0 0.5 1 1.5 2];
L = {'0' '0.5' '1.0' '1.5' '2.0'};
set(gca,'XTick',T,'XTickLabel',L);
T = [0.5 1 1.5 2 2.5 3];
L = {'0.5' '1.0' '1.5' '2.0' '2.5' '3.0'};
set(gca,'YTick',T,'YTickLabel',L);


%% 2D version of the differential equation 
subplot(122)
F = reshape(f, 100, 100);
h = imagesc(x,x, F);
h1 = get(h, 'Parent'); 
set(h1, 'Linewidth', 2, 'FontSize', 18); 
colormap(jet);
shading interp 
xlabel('$x$','Interpreter','LaTex');
ylabel('$f(x)$','Interpreter','LaTex');
title('FDDER1'); 
axis equal tight
