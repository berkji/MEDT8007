%%
clear all;
close all;


%% use field II to calculate SIR and comparsion between using xdc_focus
f_c=2.5e6; % center frequency of the transducer/aperture in Hz
f_Sample=4*f_c; % sample frequency in Hz
N_FFT = 1024;
c_0 = 1500;
field_init(0);
set_field('c',c_0);
set_sampling(f_Sample);
height = 5/1000;
width = 1/1000;
kerf = width/4;
no_elements = 64;
focus = [0 0 60]/1000;
TxAperature=xdc_linear_array(no_elements, width, height, kerf, 1, 1, focus);
figure;
show_xdc_geir(TxAperature, 1);
axis equal;
view(3);
x_0 = 0;
y_0 = 0;
z_0 = 60e-3;

xdc_focus(TxAperature, 0, [x_0 y_0 z_0]);

[SIR_4, t_field_sir]=calc_h(TxAperature,[x_0 y_0 z_0]);
figure;
subplot(2,1,1);
plot(t_field_sir + (0:length(SIR_4)-1)/f_Sample, SIR_4/max(SIR_4));
stitle=sprintf("use field-ii to calculate SIR with xdc_focus on [%0.2g %0.2g %0.2g]",[x_0 y_0 z_0]);
title(stitle);
xlabel("time (s)");

subplot(2,1,2);
f_field_sir = linspace(0,f_Sample, N_FFT);
y_sir = fft(SIR_4, N_FFT);
plot(f_field_sir/1e6, abs(y_sir(:)));
xlabel("freq (Mhz)");

steeringAngle = -35; 
focusRange=60e-3;
focus = focusRange*[sin(steeringAngle*pi/180),0,cos(steeringAngle*pi/180)];
xdc_focus(TxAperature, 0, focus);
[SIR_4, t_field_sir]=calc_h(TxAperature,focus);
figure;
subplot(2,1,1);
plot(t_field_sir + (0:length(SIR_4)-1)/f_Sample, SIR_4/max(SIR_4));
stitle=sprintf("use field-ii to calculate SIR with xdc_focus on [%0.2g %0.2g %0.2g]",focus);
title(stitle);
xlabel("time (s)");

subplot(2,1,2);
f_field_sir = linspace(0,f_Sample, N_FFT);
y_sir = fft(SIR_4, N_FFT);
plot(f_field_sir/1e6, abs(y_sir(:)));
xlabel("freq (Mhz)");



%% use field II to calculate SIR in the x-z field
steeringAngle = -35; 
focusRange=60e-3;
focus = focusRange*[sin(steeringAngle*pi/180),0,cos(steeringAngle*pi/180)];
TxAperature=xdc_linear_array(no_elements, width, height, kerf, 1, 1, focus);
xdc_focus(TxAperature,0, focus);
Nx = 161; Nz = 30;
x0=linspace(-40e-3,40e-3,Nx)';
z0=linspace( 5e-3,80e-3,Nz)';
[X,Z]=meshgrid(x0,z0);
measure_point = [X(:), zeros(length(X(:)),1),Z(:)];
[h_tx0, t_start]=calc_h(TxAperature, measure_point);
figure;
bp = sqrt(mean(h_tx0.^2));
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



