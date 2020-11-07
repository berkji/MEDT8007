% Eksempel på enkelt fasestyrt array i Field II

%% make sure that field_init has been called 
eval('field_init(0)','1;')

%% DEFINE ARRAY 
c = 1540;                   % Speed of sound
f0 = 2.5e6;                 % Transducer center frequency [Hz]
fs = 100e6;                 % Sampling frequency [Hz]
lambda = c/f0;              % Wavelength
element_height = 13/1000;   % Height of element [m] (elevation direction)
element_width = 18.5/1000;       % Element width [m] (azimuth direction)
kerf = 0;
focus = [0 0 60]/1000;      % Fixed emitter focal point [m] (irrelevant for single element transducer)
N_elements = 1;          % Number of physical elements in array
N_sub_x = 90;                % Element sub division in x-direction
N_sub_y = 90;               % Element subdivision in y-direction


%% GENERATE TRANSMIT APERTURE 
emit_aperture = xdc_linear_array (N_elements, element_width, element_height, kerf, N_sub_x, N_sub_y, focus);
eval('close(1)','1;')
figure(1)
show_xdc_geir(emit_aperture, 1); 
axis equal; view(3)
h_txAp = gcf;

%% SET THE IMPULSE RESPONSE AND EXCITATION OF THE TRANSMIT APERTURE 
t_ir = -2/f0:1/fs:2/f0;
Bw = 0.6;
impulse_response=gauspuls(t_ir,f0,Bw);
set_sampling(fs);
xdc_impulse (emit_aperture, impulse_response);

ex_periods = 1.5;
t_ex=(0:1/fs:ex_periods/f0);
excitation=square(2*pi*f0*t_ex);
xdc_excitation (emit_aperture, excitation);

figure(2); 
subplot(211);plot(t_ex*1e6, excitation);ylim([-1.1 1.1]); 
title('Excitation'); xlabel('Time [\mus]');ylabel('Amplitude');
subplot(212);plot(t_ir*1e6, impulse_response);ylim([-1.1 1.1]);
title('Transducer impulse response'); xlabel('Time [\mus]');ylabel('Amplitude');

figure(3);
freqz(impulse_response, 1, 1024, fs)


%% DEFINE MEASUREMENT POINTS
measurement_points = [20,0,60]/1000;
figure(h_txAp); 
hold on; 
plot3(measurement_points(1)*1000,measurement_points(2)*1000,measurement_points(3)*1000, 'o', 'linewidth', 6)
axis tight

%% CALCULATE SPATIAL IMPULSE RESPONSE AND TRANSMIT PRESSURE
[spatImpResp_tx, startTime_tx] = calc_h(emit_aperture, measurement_points);
[pressure_tx, startTime_tx] = calc_hp(emit_aperture, measurement_points);

%% PLOT RESULTS
figure(4);
tAxh = startTime_tx + (0:length(spatImpResp_tx)-1)/fs;
plot(tAxh*1e6, spatImpResp_tx);
title(sprintf('Spatial impulse response at point (x=%0.3gmm, y=%0.3gmm, z=%0.3gmm)', ...
    measurement_points(1)*1000, measurement_points(2)*1000, measurement_points(3)*1000))
xlabel('Time [us]'); 
axis tight

figure(5);
tAxhp = startTime_tx + (0:length(pressure_tx)-1)/fs;
plot(tAxhp*1e6, pressure_tx);
title(sprintf('Transmit pressure field at point (x=%0.3gmm, y=%0.3gmm, z=%0.3gmm)', ...
    measurement_points(1)*1000, measurement_points(2)*1000, measurement_points(3)*1000))
xlabel('Time [us]'); 
axis tight