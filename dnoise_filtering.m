% CREATE HIGH FREQUENCY NOISE
Nx = 128;
Ny = Nx;
f = rand(Nx,Ny) - 0.5;
% CALCULATE SPECTRUM
F = fft2(f);
% FILTER SPECTRUM
nx = round(0.05*Nx);
ny = round(0.05*Ny);
F(nx:Nx-nx,:) = 0;
F(:,ny:Ny-ny) = 0;
% RECONSTRUCT LOW FREQUENCY NOISE
f2 = real(ifft2(F));

imagesc(f2)
colormap(jet); 