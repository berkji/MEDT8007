close all;
clear all;

no_elements = 64;
pitch = 0.29e-3;
kerf = 0.025e-3;
width = pitch - kerf;
height=13e-3;
no_sub_x = 5;
no_sub_y = 30;
focus = [0 0 60]/1000;
Rfocus = 60e-3;
field_init(0);

Th = xdc_focused_array(no_elements, width, height, kerf, Rfocus, no_sub_x, no_sub_y, focus);
%figure;
%show_xdc_geir(Th, 1);
%% setup tx Apod
%txApodWeights = ones(1, no_elements);
txApodWeights = hanning(no_elements)';
%txApodWeights = tukeywin(no_elements, 0.3)';
figure;
stem(txApodWeights);
xdc_apodization(Th, 0, txApodWeights);


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
N_Points = 101;
x0=linspace(-20e-3,20e-3,N_Points);
y0=linspace(0e-3,0e-3,N_Points);
z0=linspace(60e-3,60e-3,N_Points);;
measure_point=[x0', y0', z0'];
figure;
plot3(x0*1e3,y0*1e3,z0*1e3,'o', 'linewidth', 6);
axis tight;
xlabel("x mm"); ylabel("y mm"); zlabel("z mm");


[h_x0, t_start]=calc_h(Th, measure_point);
figure;
tAxhp = t_start+(0:length(h_x0)-1)*dt;

%plot(t_start+(0:length(h_x0)-1)*dt, h_x0(1));
imagesc(x0*1000, tAxhp*1e6, h_x0);
xlabel("x (mm)");
ylabel("t (us)");
stitle = sprintf("SIR at %2.2d~%2.2d mm, depth:%2.2dmm",min(x0),max(x0),min(z0));
title(stitle);
cmap = jet(256);
cmap(1,:)=[0 0 0];
colormap(cmap);


%% Spatial Impulse response from x=-20mm to x=20mm, depth 20mm

N_Points = 101;
x0=linspace(-20e-3,20e-3,N_Points);
y0=linspace(0e-3,0e-3,N_Points);
z0=linspace(20e-3,20e-3,N_Points);;
measure_point=[x0', y0', z0'];
% figure;
% plot3(x0*1e3,y0*1e3,z0*1e3,'o', 'linewidth', 6);
% axis tight;
% xlabel("x mm"); ylabel("y mm"); zlabel("z mm");


[h_x0, t_start]=calc_h(Th, measure_point);
figure;
tAxhp = t_start+(0:length(h_x0)-1)*dt;

%plot(t_start+(0:length(h_x0)-1)*dt, h_x0(1));
imagesc(x0*1000, tAxhp*1e6, h_x0);
xlabel("x (mm)");
ylabel("t (us)");
stitle = sprintf("SIR at %2.2d~%2.2d mm, depth:%2.2dmm",min(x0),max(x0),min(z0));
title(stitle);
cmap = jet(256);
cmap(1,:)=[0 0 0];
colormap(cmap);


%% Spatial Impulse response from x=-20mm to x=20mm, depth 90mm
N_Points = 101;
x0=linspace(-20e-3,20e-3,N_Points);
y0=linspace(0e-3,0e-3,N_Points);
z0=linspace(90e-3,90e-3,N_Points);;
measure_point=[x0', y0', z0'];
% figure;
% plot3(x0*1e3,y0*1e3,z0*1e3,'o', 'linewidth', 6);
% axis tight;
% xlabel("x mm"); ylabel("y mm"); zlabel("z mm");


[h_x0, t_start]=calc_h(Th, measure_point);

figure;
tAxhp = t_start+(0:length(h_x0)-1)*dt;

%plot(t_start+(0:length(h_x0)-1)*dt, h_x0(1));
imagesc(x0*1000, tAxhp*1e6, h_x0);
xlabel("x (mm)");
ylabel("t (us)");
stitle = sprintf("SIR at %2.2d~%2.2d mm, depth:%2.2dmm",min(x0),max(x0),min(z0));
title(stitle);
cmap = jet(256);
cmap(1,:)=[0 0 0];
colormap(cmap);


%% Pressure response from x=-20mm to x=20mm, depth 60mm

%define measurement point matrix from x=-20mm to x=20mm
N_Points = 101;
x0=linspace(-20e-3,20e-3,N_Points);
y0=linspace(0e-3,0e-3,N_Points);
z0=linspace(60e-3,60e-3,N_Points);;
measure_point=[x0', y0', z0'];

[hp_x0, t_start]=calc_hp(Th, measure_point);
figure;
tAx_hp = t_start+(0:length(hp_x0)-1)/fs;

%plot(t_start+(0:length(h_x0)-1)*dt, h_x0(1));
imagesc(x0*1000, tAx_hp*1e6, hp_x0);
xlabel("x (mm)");
ylabel("t (us)");
stitle = sprintf("pressure at %2.2d~%2.2d mm, depth:%2.2dmm",min(x0),max(x0),min(z0));
title(stitle);
cmap = gray(256);
colormap(cmap);
axis tight;

%% BEAM PROFILE
beamprofile = sqrt(mean(hp_x0.^2));
beamprofile_60 = beamprofile;
beampro = beamprofile/max(beamprofile);
figure;
plot(x0*1000,beampro);
title(sprintf('Beamprofile at depth =%0.3gmm', z0(1)*1000))
xlabel('Azimuth position [mm]');


%% Pressure response from x=-20mm to x=20mm, depth 20mm

%define measurement point matrix from x=-20mm to x=20mm
N_Points = 101;
x0=linspace(-20e-3,20e-3,N_Points);
y0=linspace(0e-3,0e-3,N_Points);
z0=linspace(20e-3,20e-3,N_Points);;
measure_point=[x0', y0', z0'];

[hp_x0, t_start]=calc_hp(Th, measure_point);
figure;
tAx_hp = t_start+(0:length(hp_x0)-1)/fs;

%plot(t_start+(0:length(h_x0)-1)*dt, h_x0(1));
imagesc(x0*1000, tAx_hp*1e6, hp_x0);
xlabel("x (mm)");
ylabel("t (us)");
stitle = sprintf("pressure at %2.2d~%2.2d mm, depth:%2.2dmm",min(x0),max(x0),min(z0));
title(stitle);
cmap = gray(256);
colormap(cmap);
axis tight;

%% BEAM PROFILE
beamprofile = sqrt(mean(hp_x0.^2));
beamprofile_20 = beamprofile;
beampro = beamprofile/max(beamprofile);
figure;
plot(x0*1000,beampro);
title(sprintf('Beamprofile at depth =%0.3gmm', z0(1)*1000))
xlabel('Azimuth position [mm]');


%% Pressure response from x=-20mm to x=20mm, depth 90mm

%define measurement point matrix from x=-20mm to x=20mm
N_Points = 101;
x0=linspace(-20e-3,20e-3,N_Points);
y0=linspace(0e-3,0e-3,N_Points);
z0=linspace(90e-3,90e-3,N_Points);;
measure_point=[x0', y0', z0'];

[hp_x0, t_start]=calc_hp(Th, measure_point);
figure;
tAx_hp = t_start+(0:length(hp_x0)-1)/fs;

%plot(t_start+(0:length(h_x0)-1)*dt, h_x0(1));
imagesc(x0*1000, tAx_hp*1e6, hp_x0);
xlabel("x (mm)");
ylabel("t (us)");
stitle = sprintf("pressure at %2.2d~%2.2d mm, depth:%2.2dmm",min(x0),max(x0),min(z0));
title(stitle);
cmap = gray(256);
colormap(cmap);
axis tight;

%% BEAM PROFILE
beamprofile = sqrt(mean(hp_x0.^2));
beamprofile_90 = beamprofile;
beampro = beamprofile/max(beamprofile);
figure;
plot(x0*1000,beampro);
title(sprintf('Beamprofile at depth =%0.3gmm', z0(1)*1000))
xlabel('Azimuth position [mm]');


%% BEAM PROFILE at 20, 60, 90mm
beamprofile_all=[beamprofile_20 beamprofile_60 beamprofile_90];
beampro = beamprofile_20/max(beamprofile_all);
figure;
plot(x0*1000,beampro);
title(sprintf('Beamprofile at depth =20,200,90mm'));
xlabel('Azimuth position [mm]');

beampro = beamprofile_60/max(beamprofile_all);
hold on;
plot(x0*1000,beampro);

beampro = beamprofile_90/max(beamprofile_all);
hold on;
plot(x0*1000,beampro);



%% fourier transform of aperature
aTx = no_elements*pitch - kerf;
c = 1540;
R = 60e-3;
beampro_f = sinc(x0*aTx*f0/(R*c));
hold on;
plot(x0*1000,beampro_f);
legend("20mm","60mm","90mm","aperature fourier");

%
% %define a measurement point
% x0=20e-3;
% y0=0e-3;
% z0=60e-3;
% measure_point=[x0 y0 z0];
% [h_x0, t_start]=calc_h(Th, measure_point);
% figure;
% plot(t_start+(0:length(h_x0)-1)*dt, h_x0);
% xlabel("time (s)");
% stitle = sprintf("spatial impulse response at (%2.2d %2.2d %2.2d)", measure_point);
% title(stitle);
% 
% %sub problem4
% %define a measurement point
% x0=20e-3;
% y0=10e-3;
% z0=60e-3;
% measure_point=[x0 y0 z0];
% [h_x0, t_start]=calc_h(Th, measure_point);
% figure;
% plot(t_start+(0:length(h_x0)-1)*dt, h_x0);
% xlabel("time (s)");
% stitle = sprintf("spatial impulse response at (%2.2d %2.2d %2.2d)", measure_point);
% title(stitle);

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


