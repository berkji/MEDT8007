clc
close all
clear

% -- build lookup table
% define sampling grid
axial_lut = [20 100]/1000;
lateral_lut = [0 7]/1000; 

axial_fov = [20 100]/1000;
lateral_fov = [-7 7]/1000; 

FWHM = 2/1000;
f0 = 2.5e6;
fs = 50e6;
c = 1540;
nz = 128;
nx = 32;
z_axis_lut = linspace(axial_lut(1),axial_lut(2),nz);
x_axis_lut = linspace(lateral_lut(1),lateral_lut(2),nx);
[x_meas,z_meas] = meshgrid(x_axis_lut,z_axis_lut); 
y_meas = zeros(size(x_meas));
measurement_points = [x_meas(:), y_meas(:), z_meas(:)];

if ~exist('LUT_1.mat','file')
    %addpath '..\exercise2\fieldii';
    eval('field_init(0)','1;')
    
    lambda = c/f0;              % Wavelength
    element_height = 13/1000;   % Height of element [m] (elevation direction)
    pitch = 0.290/1000;         % Distance between element centers
    kerf = 0.025/1000;          % Width of fill material between the ceramic elements
    element_width = pitch-kerf; % Element width [m] (azimuth direction)
    focus = [0 0 60]/1000;      % Fixed emitter focal point [m] (irrelevant for single element transducer)
    N_elements = 64;            % Number of physical elements in array
    N_sub_x = 5;                % Element sub division in x-direction
    N_sub_y = 15;               % Element subdivision in y-direction

    % GENERATE TRANSMIT AND RECEIVE APERTURE 
    emit_aperture = xdc_linear_array (N_elements, element_width, element_height, kerf, N_sub_x, N_sub_y, focus);
    receive_aperture = xdc_linear_array (N_elements, element_width, element_height, kerf, N_sub_x, N_sub_y, focus); 
   
    % SET THE IMPULSE RESPONSE AND EXCITATION OF THE TRANSMIT AND RECEIVE APERTURE 
    t_ir = -2/f0:1/fs:2/f0;
    Bw = 0.6;
    impulse_response=gauspuls(t_ir,f0,Bw);
    set_sampling(fs);
    xdc_impulse (emit_aperture, impulse_response);
    xdc_impulse (receive_aperture, impulse_response);

    % SET THE EXCITATION OF THE TRANSMIT APERTURE 
    ex_periods = 1.5;
    t_ex=(0:1/fs:ex_periods/f0);
    excitation=square(2*pi*f0*t_ex);
    xdc_excitation (emit_aperture, excitation);
    
    apo_emit = ones(1,N_elements);
    apo_receive = ones(1,N_elements);   % Rectangular apodization
    
    xdc_apodization(emit_aperture, 0, apo_emit);
    xdc_focus(receive_aperture, 0, focus);
    [simData, startTime] = calc_hhp(emit_aperture, receive_aperture, measurement_points);
    
    LUT = reshape(max(simData),nz,nx);
    save LUT.mat LUT
%     simData = reshape(max(simData),nz,nx);
        
else
    load LUT
end

imagesc(x_axis_lut,z_axis_lut,20*log10(LUT))
    
% -- phantom
n_scatt = 20;
start_depth = 30/1000; 
stop_depth = 90/1000; 

% -- probe
Number_of_lines = 30; 

% -- psf 
n_oscill = 3;
t = 0:1/fs:n_oscill/f0;
t = t - mean(t);
s_y = max(t)/3; % in sec!!
s_x = FWHM / (2*sqrt(2*log(2)));
P_y = exp(-.5.*(t/s_y).^2).*cos(2*pi*f0*t + pi/6);

% -- sampled x and y axis
dr = c / (2*fs);
y_axis = axial_fov(1):dr:axial_fov(2);
x_axis = linspace(lateral_fov(1),lateral_fov(2),Number_of_lines); % scan lines

% -- tissue
y_scatt = linspace(start_depth,stop_depth,n_scatt);
x_scatt = zeros(size(y_scatt));
amp = 1;

% -- image
rf_image = zeros(numel(y_axis),numel(x_axis));
cont = 0;
for ii = 1:Number_of_lines
    
    pos_x = x_axis(ii); % -- scan line position
    refl = zeros(size(y_axis)); % reflectivity

    % -- project scatterers on sampled y axis
    
    y_scatt_idx = (y_scatt - axial_fov(1)) / diff(axial_fov) * numel(y_axis-1) + 1;  
    y_scatt_idx_f = floor(y_scatt_idx);
    y_scatt_idx_c = ceil(y_scatt_idx);

    w_f = 1 - abs(y_scatt_idx - y_scatt_idx_f);
    w_c = 1 - abs(y_scatt_idx - y_scatt_idx_c);
    d = abs(x_scatt - pos_x);
    % -- weights
    if ii == 1
        figure
        imagesc(x_axis_lut,z_axis_lut,20*log10(LUT)), colormap gray, axis image, hold on,
        plot(d,y_scatt,'go','markerfacecolor','g')
        
    end
    
    w_amp = interp2(x_axis_lut,z_axis_lut,LUT,d,y_scatt,'cubic',0); 
    figure;
    stem(w_amp);
    title(sprintf("weight at %2g",ii));
    
    refl(y_scatt_idx_f) = w_amp.*w_f;
    refl(y_scatt_idx_c) = w_amp.*w_c;
        
    % -- convolution
    rf_image(:,ii) = conv(refl,P_y,'same');

    if cont == 1
        figure
        stem(y_axis,refl), hold on;
        plot(y_scatt,ones(size(y_scatt)),'ko','markerfacecolor','g');

        figure
        plot(y_axis,rf_image(:,ii)),

    end
            
end

bmode = abs(hilbert(rf_image));
bmode = bmode / max(bmode(:));
bmode = 20*log10(bmode) + 50;
bmode(bmode<0) = 0;
% bmode = imresize(bmode,[512 256]);
% close all
figure
imagesc(x_axis,y_axis,bmode); axis image, colormap gray,
xlabel('[m]')
ylabel('[m]')
