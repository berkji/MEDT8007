% Starting-point-script for exercise 1 in the Simulation Methods in
% Ultrasound Imaging course (MEDT8007).
% 
% self-study for focused SIR response
% 
%
% author:  JI Kai
% draft:   Nov/2 2020
%%
% 
% update:  
%
% uses:   
% sub of:  

%%
% 
% 
% 
% 
% 
% 


close all
fclose('all');
clear all;


% Parameters (see figure 3.16 on page 156)

height = 5/1000;
width = 1/1000;
kerf = width/4;
no_elements = 64;
focus = [0 0 60]/1000;


ds  = 0.2e-3; % spatial resolution in m
x_p = ((width+kerf)*0.5*no_elements - kerf/2);
x_m = -1*x_p;
x_1 = linspace(x_m,x_p,no_elements); % x-coordinate of the (candidate) aperture point(s) in m
y_1 = 0e-3;       % y-coordinate of the (candidate) aperture point(s) in m
x_a = x_p;

f_c=2.5e6; % center frequency of the transducer/aperture in Hz
f_Sample=4*f_c; % sample frequency in Hz
N_FFT=1024*1; % length of the FFT
t_s=0:1/f_Sample:0.8e-6; % non-zero duration of the excitation function in s
s=[sin(0.25*2*pi*f_c*t_s).^2.*sin(2*pi*f_c*t_s) zeros(1, N_FFT - length(t_s))]; % excitation function
c_0=1500; % speed of sound in m/s
df=f_Sample/N_FFT; % frequency resolution in Hz
f_c = 0:df:floor(N_FFT/2)*df; % 0.1e6:0.1e5:6e6; 
f=f_c; % frequency (range) Hz

steeringAngle = -35; 
focusRange=60e-3;
focus = focusRange*[sin(steeringAngle*pi/180),0,cos(steeringAngle*pi/180)];
x_f = focus(1);
y_f = focus(2);
z_f = focus(3);

%x_0 =-30e-3:0.5e-3:30e-3; % x-coordinate of the observation point(s) in m
x_0 = 0;
z_0 = z_f;

% Huygens' principle
R=zeros(length(x_1), length(y_1));


H=zeros(length(x_0), length(f)); % complex amplitude function
h_WB=waitbar(0);
tic
for m=1:length(x_1)
  for n=1:length(y_1)
    for q=1:length(x_0)
      if sqrt(x_1(m)^2 + y_1(n)^2)<= x_a % only accept points that lie within the aperture
        %<frequency independent variables>
		xf = [x_f 0 z_f];
        xp = [x_1(m) 0 0];
        xpa =[x_a 0 0];
        xma =[-x_a 0 0];
        x0 = [x_0(q) 0 z_0];
        r = x0 - xp;
        R = xf - xma;
        P = xf - xp;
        rpf(m,n) = sqrt(sum(R.*R)) - sqrt(sum(P.*P));
        rr(m,n)  = sqrt(sum(r.*r));
        for i=1:length(f)
		  k(i) = 2*pi*f(i)/c_0;
          H(q, i)=H(q, i) + exp(-1i*k(i)*(rpf(m,n)+rr(m,n)))/(2*pi*rr(m,n));
%         H(q, i)=H(q, i) + exp(-1i*k(i)*(rr(m,n)))/(2*pi*rr(m,n));
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
%plot Normalized Green function frequency response  
%plot(x_0, abs(H(:))/max(abs(H(:))),'b');
%xlabel("x axis (mm)");
%plot(f/1e6, abs(H(:))/max(abs(H(:))),'b');
plot(f/1e6, abs(H(:)),'b');
xlabel("freq (Mhz)");
ylabel("Impulse_Response_Mag");


%impulse response 
t_s=0:1/f_Sample:0.8e-6; % non-zero duration of the excitation function in s
X= fft(s);
Y_Phi = X(1:length(H(1,:))).*H(1,:);
Phi = ifftmmv(Y_Phi, N_FFT);
Phi_a = real(ifft(Y_Phi,N_FFT));
t=(0:length(Phi) - 1)/f_Sample; % time axes
figure;
plot(t*1e6, Phi);
xlabel("time (us)");
hold on;
plot(t*1e6, Phi_a);
xlabel("time (us)");
legend("ifftmmv", "ifft");

n_b = floor(min(rpf(:) + rr(:))/c_0*f_Sample)-1;
n_e = ceil(max(rpf(:) + rr(:))/c_0*f_Sample)+1;
n=n_b + (length(t_s) + (n_e - n_b + 1) - 1) - 1;
n=min(n,N_FFT);

figure;
plot(t,s,'r');
title("original source signal");
xlabel("time (s)");

figure;
plot(t(n_b:n)*1e6, Phi(n_b:n)/max(abs(Phi(n_b:n))), 'b-', 'LineWidth', 2);
title("use impulse response to calcuate pressure signal (zoom in)");
xlabel("time (us)");


%% use mycalc_h to do the SIR 
steeringAngle = -15; 
focusRange=60e-3;
focus = focusRange*[sin(steeringAngle*pi/180),0,cos(steeringAngle*pi/180)];
x_f = focus(1);
y_f = focus(2);
z_f = focus(3);

focus_point = [focus(1) focus(2) focus(3)];

x_0 = 0;
z_0 = z_f;
y_0 = 0;
observation_point = [x_0 y_0 z_0];

f_c = 2.5e6;

Nx = 161; Nz = 30;
x0=linspace(-40e-3,40e-3,Nx)';
z0=linspace( 5e-3,80e-3,Nz)';
[X,Z]=meshgrid(x0,z0);
measure_point = [X(:), zeros(length(X(:)),1),Z(:)];
mysir_array = zeros(Nx*Nz,N_FFT);
for i= 1 : Nx*Nz 
    observation_point = measure_point(i,:);
    [mysir_t, my_tstart]=mycalc_h(no_elements, width, kerf, f_c, focus_point, observation_point);
    mysir_array(i,:)=mysir_t;
end

mysir_array = mysir_array';

figure;
bp = sqrt(mean(mysir_array.^2));
bp = reshape(bp, Nz, Nx);
bp = bp/max(bp(:));
pcolor(x0*1000,z0*1000, 20*log10(bp));
shading interp
title("spatial response [dB]")
xlabel('Azimuth [mm]');
ylabel('Range [mm]');
caxis([-50 0]); % Set dynamic range
colormap(jet(256));
colorbar


% myH = fft(mysir_t, N_FFT);
% 
% figure;
% %plot Normalized Green function frequency response  
% plot(f/1e6, abs(H(:)),'b');
% xlabel("freq (Mhz)");
% ylabel("Impulse_Response_Mag");
% title("mycalc_h");