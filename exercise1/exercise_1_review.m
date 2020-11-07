% Review-script for exercise 1 in the Simulation Methods in Ultrasound
% Imaging course (MEDT8007).
% 
%
% author:  Marco M. Voormolen
% draft:   25 March 2008

% update:  
%
% uses:   ifftmmv.m and sirmmv.m
% sub of:  


clear all
close all
fclose('all');


% Parameters (see figure 3.16 on page 156)
a=3e-3; % aperture radius in m
ds=0.4e-3; % spatial resolution in m
x_1=-a:ds:a; % x-coordinate of the (candidate) aperture point(s) in m
y_1=-a:ds:a; % y-coordinate of the (candidate) aperture point(s) in m
x_0=5e-3; %-30e-3:0.5e-3:30e-3; % x-coordinate of the observation point(s) in m
y_0=0e-3; % y-coordinate of the observation point(s) in m
z_0=50e-3; % z-coordinate of the observation point(s) in m
f_c=2.5e6; % center frequency of the transducer/aperture in Hz
f_Sample=10*f_c; % sample frequency in Hz
N_FFT=128; % length of the FFT
t_s=0:1/f_Sample:0.8e-6; % non-zero duration of the excitation function in s
s=[sin(0.25*2*pi*f_c*t_s).^2.*sin(2*pi*f_c*t_s) zeros(1, N_FFT - length(t_s))]; % excitation function
c_0=1500; % speed of sound in m/s
df=f_Sample/N_FFT; % frequency resolution in Hz
f=0:df:floor(N_FFT/2)*df; %f_c; % frequency (range) Hz


% Huygens' principle
R=zeros(length(x_1), length(y_1));
H=zeros(length(x_0), length(f)); % complex amplitude function
h_WB=waitbar(0);
tic
for m=1:length(x_1)
  for n=1:length(y_1)
    for q=1:length(x_0)
      if sqrt(x_1(m)^2 + y_1(n)^2)<=a % only accept points that lie within the radius of the aperture
        r=[x_0(q) - x_1(m); y_0 - y_1(n); z_0]; % vector from the aperture point to observation point
        R(m, n)=sqrt(sum(r.^2)); % distance from the aperture point to observation point
        for i=1:length(f)
          Lambda=c_0/f(i); % wave length
          k(i)=(2*pi)/Lambda; % wave number
          H(q, i)=H(q, i) + exp(-1i*k(i)*R(m, n))/R(m, n); % because we are not taking care of the correct amplitude '4*pi' is omitted
        end
      end  
    end  
  end  
  waitbar(m/length(x_1), h_WB);
end
t_1=toc;
close(h_WB)
n=find(R==0);
R(n)=NaN;

figure;
%plot(x_0*1000, abs(H)/max(abs(H)), 'b-', 'LineWidth', 2);
plot(f, abs(H)/max(abs(H)), 'b-', 'LineWidth', 2);
%plot(x_0*1000, 20*log10(abs(H)/max(abs(H))), 'b-', 'LineWidth', 2)
figure;
Theta=atan(x_0/z_0);
%plot(x_0*1000, abs(2*besselj(1, k*a*sin(Theta))./(k*a*sin(Theta))), 'r--', 'LineWidth', 2);
plot(f, abs(2*besselj(1, k*a*sin(Theta))./(k*a*sin(Theta))), 'r--', 'LineWidth', 2);
%plot(x_0*1000, 20*log10(abs(2*bessel(1, k*a*sin(Theta))./(k*a*sin(Theta)))), 'r--', 'LineWidth', 2);
xlabel('x_0 [mm]=')
title(['z_0=' num2str(z_0*1000, '%.2f') 'mm'])


% Acoustic field

% FFT
S=fft(s);
P=S(1:length(H)).*H;
p=ifftmmv(P, N_FFT);
N_Rep=ceil((max(R(:))/c_0 + length(t_s)/f_Sample)/(N_FFT/f_Sample)); % number of repetition of the periodic time signal
p=repmat(p, 1, N_Rep); % periodic pressure signal
t=(0:length(p) - 1)/f_Sample; % time axesk
n_b=floor((min(R(:))/c_0)*f_Sample) - 1; % first sample of the transient time signal
n_e=ceil((max(R(:))/c_0)*f_Sample) + 1; % last sample of the transient time signal
figure
hold on
n=n_b + (length(t_s) + (n_e - n_b + 1) - 1) - 1;
plot(t(n_b:n), p(n_b:n)/max(abs(p(n_b:n))), 'b-', 'LineWidth', 2)

% convolution
h=ifftmmv(H, N_FFT);
h=repmat(h, 1, N_Rep); % periodic SIR
p=conv(s, h(n_b:n_e));
plot(t(n_b) + (0:length(p) - 1)/f_Sample, p/max(abs(p)), 'r--', 'LineWidth', 2)


% SIR
[h, t_0]=sirmmv(a, x_0, y_0, z_0, f_Sample, c_0);
p=conv(s, h);
t_2=toc;
plot(t_0 + (0:length(p) - 1)/f_Sample, p/max(abs(p)), 'm-.', 'LineWidth', 2)


% Field II
field_init(0)
set_field('c', c_0)
set_sampling(f_Sample);
N_Subdiv=10; % number of lateral subdivisions
TxAperture=xdc_piston(a, 2*a/N_Subdiv);
s_R=(2*pi*f_c*(2*a/N_Subdiv)^2)/(2*c_0) % Rayleigh distance (square mathematical elements are assumed)
xdc_excitation(TxAperture, s);
tic
[p_FII, t_FII]=calc_hp(TxAperture, [x_0 y_0 z_0]);
t_3=toc;
plot(t_FII + (0:length(p_FII) - 1)/f_Sample, p_FII/max(abs(p_FII)), 'g:', 'LineWidth', 2)
xdc_free(TxAperture)
field_end
legend('Huygens FFT', 'Huygens convolution', 'SIR', 'Field II')