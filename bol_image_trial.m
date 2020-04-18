%% INITIALIZE 
clc 
close all 
clear 

%% DASHBOARD 
Sx = 1; % physical distance 
Sy = 1; % physical distance in y 

Nx = 312; % grid size 
Ny = 312; % grid size 

dx = Sx/Nx; 
dy = Sy/Ny; 

xa = [0:Nx-1]*dx; xa= xa-mean(xa); 
ya = [0:Ny-1]*dy; ya= ya-mean(ya); 

[Y, X] = meshgrid(ya, xa); 

AA = abs(Y)<= .075; 
BB = abs(X)<= .075; 

ER = AA|BB; % using boolean operator. 

%% Blurring using guassian
r = 0.1; 
E = exp(-(X.^2+Y.^2)/r^2); 

KK = fft2(ER).*fft2(E)/sum(E(:)); % blurring of image

ER = ifftshift(real(ifft2(KK))); 

ER = ER>0.3 ; 

imagesc((real((ER)))); 
colormap(jet); 