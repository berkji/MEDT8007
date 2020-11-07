close all;
clear all;

no_elements = 32;
pitch = 0.29e-3;
kerf = 0.025e-3;
width = pitch - kerf;
height=13e-3;
no_sub_x = 5;
no_sub_y = 30;
focus = [0 0 60]/1000;
Rfocus = 60e-3;
c =1540;

field_init(0);

Th = xdc_focused_array(no_elements, width, height, kerf, Rfocus, no_sub_x, no_sub_y, focus);
figure;
show_xdc_geir(Th, 1);
% 
fs = 100e6; %sampling freq (100Mhz)
f0 = 2.5e6; % transducer center freq (2.5Mhz)
t0 = 1/f0; 
dt = 1/fs;  %sampling period
set_sampling(fs);


%% Impulse setup
t_ir = -2/f0:1/fs:2/f0;
Bw = 0.6;
impulse_response = gauspuls(t_ir, f0, Bw);
xdc_impulse (Th, impulse_response);
figure;
excitation = square(2*pi*f0*(0:dt:1.5*t0));
plot(0:dt:1.5*t0, excitation);
xlabel("time (s)");
title("excitation pulse");
xdc_excitation(Th, excitation);
figure;
plot(t_ir, impulse_response);
xlabel("time (s)");
title("impulse response");
figure;
freqz(impulse_response,1,1024,fs);

%% Impulse response from x=-20mm to x=20mm, depth 60mm
N_Points = 100;
x0=linspace(0e-3,0e-3,N_Points);
y0=linspace(0e-3,0e-3,N_Points);
z0=linspace(5e-3,150e-3,N_Points);;
measure_point=[x0', y0', z0'];
figure;
plot3(x0*1e3,y0*1e3,z0*1e3,'o', 'linewidth', 6);
axis tight;
xlabel("x mm"); ylabel("y mm"); zlabel("z mm");

[hp_x0, t_start]=calc_hp(Th, measure_point);
figure;
tAx_hp = t_start+(0:length(hp_x0)-1)/fs;

%plot(t_start+(0:length(h_x0)-1)*dt, h_x0(1));
imagesc(z0*1000, tAx_hp*1e6, hp_x0);
xlabel("z (mm)");
ylabel("t (us)");
stitle = sprintf("pressure at %0.3g~%0.3gmm",min(z0),max(z0));
title(stitle);
cmap = gray(256);
colormap(cmap);
axis tight;

%% Pressure response vs depth
mhp_x0 = sqrt(mean(hp_x0.^2));
figure;
plot(z0*1000,mhp_x0/max(mhp_x0));
xlabel("z (mm)");
ylabel("p (normalized)");
stitle = sprintf("pressure at %0.3g~%0.3gmm",min(z0),max(z0));
title(stitle);


%% estimation of FWHM (Full Width Half Maximum) 
a = no_elements*pitch - kerf;
lambda = c/f0;
halfMaxWidth=7.2*lambda*(focus(3)/a)^2;
s_halfMaxWidth = sprintf("halfMaxWidth = %0.3gmm",halfMaxWidth*1000)

%% diffraction focus
Rdf = a^2/(4*lambda);
s_Rdf = sprintf("Rdf = %0.3gmm",Rdf*1000)







