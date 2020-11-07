%% make sure that field_init has been called 
eval('field_init(0)','1;')

%% DEFINE ARRAY 
c = 1540;                   % Speed of sound
f0 = 2.5e6;                 % Transducer center frequency [Hz]
fs = 100e6;                 % Sampling frequency [Hz]
lambda = c/f0;              % Wavelength
element_height = 13/1000;   % Height of element [m] (elevation direction)
pitch = 0.290/1000;         % Distance between element centers
kerf = 0.025/1000;          % Width of fill material between the ceramic elements
element_width = pitch-kerf; % Element width [m] (azimuth direction)
Rfocus = 60/1000;           % Elevation lens focus (or radius of curvature, ROC) 
focus = [0 0 60]/1000;      % Fixed emitter focal point [m] (irrelevant for single element transducer)
N_elements = 64;            % Number of physical elements in array
N_sub_x = 5;                % Element sub division in x-direction
N_sub_y = 15;               % Element subdivision in y-direction
apodTx = 0;                 % Transmit apodization. 0:boxcar, 1:Hanning, 2:Cosine-tapered 0.3

%% GENERATE TRANSMIT APERTURE 
emit_aperture = xdc_focused_array (N_elements, element_width, element_height, kerf, Rfocus, N_sub_x, N_sub_y, focus);
eval('close(1)','1;')
% figure(1)
% show_xdc_geir(emit_aperture, 1); 
% axis equal; view(3)
% h_txAp = gcf;

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

%% DEFINE APODIZATION FOR THE EMIT APERTURE 
switch apodTx
    case 0
        apo_emit = ones(1,N_elements);   % Rectangular apodization
    case 1
        apo_emit = hanning(N_elements)'; % Hanning apodization on the emit-apperture
    case 2
        apo_emit = tukeywin(N_elements,  0.3)'; % Cosine-tapered apodization
end
xdc_apodization(emit_aperture, 0, apo_emit);
figure(2)
stem(apo_emit);
title('Transmit apodization');
xlabel('Element index');
ylabel('Element amplitude weighting');
axis tight
ylim([0 1]);

%% DEFINE MEASUREMENT POINTS
measDepthStart = 5/1000;    % Start depth along z-axis to place measurement points
measDepthEnd = 150/1000;    % End depth along z-axis to place measurement points
xStart = -15/1000;      % Start position of measurement points in x direction
xEnd = 15/1000;         % End position of measurement points in x direction

Nmpx = 81;
Nmpz = 59;
mx = linspace(xStart,xEnd,Nmpx)';
my = zeros(Nmpx*Nmpz, 1);
mz = linspace(measDepthStart,measDepthEnd,Nmpz)';
[X,Z] = meshgrid(mx,mz);
measurement_points = [X(:),my,Z(:)];
% figure(h_txAp); 
% hold on; 
% plot3(measurement_points(:,1)*1000,measurement_points(:,2)*1000,measurement_points(:,3)*1000, 'o', 'linewidth', 6)
% axis tight

%% CALCULATE TRANSMIT PRESSURE
disp('Calculating pressure');
[pressure_tx, startTime_tx] = calc_hp(emit_aperture, measurement_points);

%% PLOT RESULTS
figure(4)
bp = sqrt(mean(pressure_tx.^2));
bp = reshape(bp, Nmpz, Nmpx);
bp=bp/max(bp(:));
pcolor(mx*1000, mz*1000, 20*log10(bp));
shading interp
title(sprintf('Transmit beam pressure field'))
xlabel('Azimuth [mm]');
ylabel('Range [mm]');
caxis([-35 0]); % Set dynamic range
colormap(jet(256));
colorbar

%% PLOT BP NORMALIZED AT EACH DEPTH
figure(5)
bp= bp./repmat(max(bp')', 1,Nmpx);
pcolor(mx*1000, mz*1000, 20*log10(bp));
shading interp
title(sprintf('Transmit beam pressure field'))
xlabel('Azimuth [mm]');
ylabel('Range [mm]');
caxis([-35 0]); % Set dynamic range
colormap(jet(256));
colorbar
