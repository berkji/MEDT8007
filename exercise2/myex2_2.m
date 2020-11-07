close all;

no_elements = 64;
pitch = 0.29e-3;
kerf = 0.020e-3;
width = pitch - kerf;
height=13e-3;
no_sub_x = 5;
no_sub_y = 15;
focus = [0 0 60]/1000;
Rfocus = 60e-3;
field_init;

Th = xdc_focused_array(no_elements, width, height, kerf, Rfocus, no_sub_x, no_sub_y, focus);
figure;
show_xdc_geir(Th, 1);
% 
fs = 100e6; %sampling freq (100Mhz)
f0 = 2.5e6; % transducer center freq (2.5Mhz)
t0 = 1/f0; 
dt = 1/fs;  %sampling period
set_sampling(fs);

% excitation = sin(2*pi*f0*(0:dt:1.5*t0));
% 
% figure;
% plot(0:dt:1.5*t0, excitation);
% xlabel("time (s)");
% title("excitation pulse");
% 
% xdc_excitation(Th, excitation);
% 
% t_ir = -2/f0:1/fs:2/f0;
% Bw = 0.6;
% impulse_response = gauspuls(t_ir, f0, Bw);
% figure;
% plot(t_ir, impulse_response);
% xlabel("time (s)");
% title("impulse response");
% 
% freqz(impulse_response,1,1024,fs);
% 
%sub problem3
%define a measurement point
x0=0e-3;
y0=0e-3;
z0=60e-3;
measure_point=[x0 y0 z0];
[h_x0, t_start]=calc_h(Th, measure_point);
figure;
plot(t_start+(0:length(h_x0)-1)*dt, h_x0);
xlabel("time (s)");
stitle = sprintf("spatial impulse response at (%2.2d %2.2d %2.2d)", measure_point);
title(stitle);

%sub problem4
%define a measurement point
x0=20e-3;
y0=0e-3;
z0=60e-3;
measure_point=[x0 y0 z0];
[h_x0, t_start]=calc_h(Th, measure_point);
figure;
plot(t_start+(0:length(h_x0)-1)*dt, h_x0);
xlabel("time (s)");
stitle = sprintf("spatial impulse response at (%2.2d %2.2d %2.2d)", measure_point);
title(stitle);

%sub problem4
%define a measurement point
x0=20e-3;
y0=10e-3;
z0=60e-3;
measure_point=[x0 y0 z0];
[h_x0, t_start]=calc_h(Th, measure_point);
figure;
plot(t_start+(0:length(h_x0)-1)*dt, h_x0);
xlabel("time (s)");
stitle = sprintf("spatial impulse response at (%2.2d %2.2d %2.2d)", measure_point);
title(stitle);

% 
% %sub problem3
% no_elements = 1;
% width=18.5e-3;
% height=13e-3;
% kerf = 0;
% no_sub_x = 30;
% no_sub_y = 30;
% focus = [0 0 60]/1000;
% 
% Th = xdc_linear_array(no_elements, width, height, kerf, no_sub_x, no_sub_y, focus);
% figure;
% show_xdc_geir(Th, 1);
% 
% fs = 100e6; %sampling freq (100Mhz)
% set_sampling(fs);
% f0 = 2.5e6; % transducer center freq (2.5Mhz)
% t0 = 1/f0; 
% dt = 1/fs;  %sampling period
% excitation = sin(2*pi*f0*(0:dt:1.5*t0));
% 
% figure;
% plot(0:dt:1.5*t0, excitation);
% xlabel("time (s)");
% title("excitation pulse");
% 
% xdc_excitation(Th, excitation);
% 
% t_ir = -2/f0:1/fs:2/f0;
% Bw = 0.6;
% impulse_response = gauspuls(t_ir, f0, Bw);
% figure;
% plot(t_ir, impulse_response);
% xlabel("time (s)");
% title("impulse response");
% 
% freqz(impulse_response,1,1024,fs);
% 
% %
% %define a measurement point
% x0=0e-3;
% y0=0e-3;
% z0=60e-3;
% measure_point=[x0 y0 z0];
% [h_x0, t_start]=calc_h(Th, measure_point);
% figure;
% plot(t_start+(0:length(h_x0)-1)*dt, h_x0);
% xlabel("time (s)");
% stitle = sprintf("spatial impulse response at (%2.2d %2.2d %2.2d)", measure_point);
% title(stitle);


