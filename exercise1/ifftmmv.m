function x=ifftmmv(X, N_FFT)
% x=ifftmmv(X, N_FFT)
%
% Inverse Fourier transform for positive frequency components. Where:
% X is the positive frequency signal and
% N_FFT is the samples in the time signal.
% 
%
% author:  Marco M. Voormolen
% draft:   25 March 2008

% update:  
%
% uses:   
% sub of:  


x=real(ifft([X flipdim(real(X(2:ceil(N_FFT/2))) - imag(X(2:ceil(N_FFT/2)))*1i, 2)])); % Always use 'real' to avoid epsilon size imaginairy results.