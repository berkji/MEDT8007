%% use field II to calculate SIR and comparsion between using xdc_impulse or not 
f_c=2.5e6; % center frequency of the transducer/aperture in Hz
f_Sample=10*4*f_c; % sample frequency in Hz
% Parameters (see figure 3.16 on page 156)
a=3e-3; % aperture radius in m
ds=0.4e-3; % spatial resolution in m
x_1=-a:ds:a; % x-coordinate of the (candidate) aperture point(s) in m
y_1=-a:ds:a; % y-coordinate of the (candidate) aperture point(s) in m
x_0=-30e-3:0.5e-3:30e-3; % x-coordinate of the observation point(s) in m
x_0=0e-3; % 5e-3
y_0=0e-3; % y-coordinate of the observation point(s) in m
z_0=60e-3; % z-coordinate of the observation point(s) in m
f_c=2.5e6; % center frequency of the transducer/aperture in Hz
f_Sample=10*f_c; % sample frequency in Hz
N_FFT=1024; % length of the FFT
c_0=1500; % speed of sound in m/s

field_init(0);
set_field('c',c_0);
set_sampling(f_Sample);
N_div = 32; % number of lateral division
TxAperature=xdc_piston(a, 2*a/N_div);
show_xdc_geir(TxAperature,1);

[SIR_4, t_field_sir]=calc_h(TxAperature,[x_0 y_0 z_0]);
figure;
plot(t_field_sir + (0:length(SIR_4)-1)/f_Sample, SIR_4/max(SIR_4));
title("use field-ii to calculate SIR without xdc_impulse");
xlabel("time (s)");

y_SIR_4 = fft(SIR_4, N_FFT);
f0 = linspace(0, f_Sample, N_FFT);
plot(f0/1e6,abs(y_SIR_4));
xlabel("freq Mhz");