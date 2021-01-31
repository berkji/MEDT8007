function plot2dFFT(data, dx, dy, Nfftx, Nffty)

H = fftshift(fft2(data, Nffty, Nfftx)); % Calculate 2D Fourier transform
H = H/max(abs(H(:)));                   % Normalize
Fx = linspace(-1/dx/2, 1/dx/2, Nfftx);  % Frequency axis in x direction 
Fy = linspace(-1/dy/2, 1/dy/2, Nffty);  % Frequency axis in y direction

figure; 
imagesc(Fx,Fy,(abs(H)));
% caxis([-40 0]); % Limiting dynamic range to 40 dB