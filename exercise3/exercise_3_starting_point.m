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
focus = [0 0 95]/1000;      % Fixed emitter focal point [m] (irrelevant for single element transducer)
N_elements = 64;            % Number of physical elements in array
N_sub_x = 1;                % Element sub division in x-direction
N_sub_y = 7;                % Element subdivision in y-direction
apodTx = 0;                 % Transmit apodization. 0:boxcar, 1:Hanning, 2:Cosine-tapered 0.3
apodRx = 0;                 % Receive apodization. 0:boxcar, 1:Hanning, 2:Cosine-tapered 0.3
dynamic_receive_focus = 1;  % Enable dynamic receive focusing.
expanding_aperture = 1;     % Enable expanding aperture.
fn_receive = 2.0;           % Receive F-number
show2dFourier = 0;          % Plot 2D Fourier transform of your image
oversamplingratio = 4;      % Oversampling based on Nyquist of the transmit beams


%% GENERATE TRANSMIT AND RECEIVE APERTURE 
emit_aperture = xdc_focused_array (N_elements, element_width, element_height, kerf, Rfocus, N_sub_x, N_sub_y, focus);
receive_aperture = xdc_focused_array (N_elements, element_width, element_height, kerf, Rfocus, N_sub_x, N_sub_y, focus);

eval('close(1)','1;')
figure(1)
show_xdc_geir(emit_aperture, 1); 
axis equal; view(3)
% h_txAp = gcf;

%% SET THE IMPULSE RESPONSE AND EXCITATION OF THE TRANSMIT AND RECEIVE APERTURE 
t_ir = -2/f0:1/fs:2/f0;
Bw = 0.6;
impulse_response=gauspuls(t_ir,f0,Bw);
set_sampling(fs);
xdc_impulse (emit_aperture, impulse_response);
xdc_impulse (receive_aperture, impulse_response);

%% SET THE EXCITATION OF THE TRANSMIT APERTURE 
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
subplot(211)
stem(apo_emit);
title('Transmit apodization');
xlabel('Element index');
ylabel('Element amplitude weighting');
axis tight
ylim([0 1]);

%% DEFINE RECEIVE FOCAL REGIONS
Nfr = 600;
focal_area_start = 0/1000;
focal_area_end = 200/1000;
region_starts = linspace(focal_area_start, focal_area_end, Nfr + 1);
receive_focal_times = (region_starts(1:end-1)/c)';            % Start time of each focal region
receive_focal_width = (focal_area_end - focal_area_start)/Nfr;  % Span of each focal region
receive_focal_distance = region_starts(1:end-1) + receive_focal_width/2;    % Center of each focal region


%% DEFINE EXPANDING APERTURE AND APODIZATION FOR CONSTANT RX F-NUMBER
if (expanding_aperture)
    N_elem_active_rec = round(receive_focal_distance./(fn_receive*(element_width+kerf)));
else
    N_elem_active_rec = ones(1, Nfr)*N_elements;
end

N_pre_receive = 0;      % Initialization
N_post_receive = 0;     % Initialization
apo_matrix_receive = zeros(Nfr, N_elements);    % Initialization

for ii=1:Nfr
    if N_elem_active_rec(ii)>N_elements
        N_elem_active_rec(ii)=N_elements; 
    end
    N_pre_receive(ii) = ceil(N_elements/2  - N_elem_active_rec(ii)/2);              % Number of element set to zero on the LEFT hand side of the apperture
    N_post_receive(ii) = N_elements - N_pre_receive(ii) - N_elem_active_rec(ii);    % Number of element set to zero on the RIGHT hand side of the apperture
    switch apodRx
        case 0
            apo_receive = ones(1,N_elem_active_rec(ii));   % Rectangular apodization
        case 1
            apo_receive = hanning(N_elem_active_rec(ii))'; % Hanning apodization on the emit-apperture
        case 2
            apo_receive = tukeywin(N_elem_active_rec(ii),  0.3)'; % Cosine-tapered apodization
    end

    % Construct a matrix containing the apodization and aperturesize for all the focalzones
    apo_matrix_receive(ii,:)=[zeros(1,N_pre_receive(ii)) apo_receive zeros(1,N_post_receive(ii))];  
end
apo_matrix_receive(1,:) = apo_matrix_receive(2,:);

xdc_apodization(receive_aperture, receive_focal_times, apo_matrix_receive);
figure(2)
subplot(212)
imagesc(1:N_elements, receive_focal_distance*1000, apo_matrix_receive);
title('Receive apodization');
xlabel('Element index');
ylabel('Distance [mm]');
axis tight


%%  DEFINE IMAGING SECTOR 
sector=35 * pi/180;                             % Size of image sector [rad]


%%  GENERATE COMPUTER PHANTOM 
phantom_count = 200;
ph_dim = [44/1000 0 30/1000]; origo = [0/1000 0 90/1000]; 
r_start = 10/1000; r_end = 100/1000; N_rp = 10; phi_r = 0*pi/180;
[phantom_positions, phantom_amplitudes] = ph_rect_gauss(phantom_count, ph_dim, origo);

figure;colordef(gcf, 'black')
scatter(phantom_positions(:,1)*1000, phantom_positions(:,3)*1000,...
    50, abs(phantom_amplitudes), 'filled');axis equal; axis ij;
title('Point scatterers'); xlabel('Azimuth [mm]'); ylabel('Range [mm]')
colormap(hot);


%% DEFINE SPATIAL SAMPLING DENSITY 
aTx = pitch * N_elements - kerf;   % Width of physical aperture
dThetaRayleigh = lambda/aTx;                  % Rayleigh resolution (in radians)
dThetaDesired = dThetaRayleigh/oversamplingratio; % Simulation resolution
N_lines = round(sector/dThetaDesired)+1;      % Number of emit lines in image
d_theta_emit = sector/(N_lines-1);                 % Increment in emit angle for 90 deg. image
d_theta_receive = d_theta_emit;

%% INITIALIZE SIMULATION 
image_data=zeros(800,N_lines);          % Preallocating storage for image data
times = zeros(1, N_lines);              % Preallocating storage for times data
theta_emit = -sector/2;                 % Start angle for transmit aperture
theta_receive = theta_emit;             % Start angle for receive aperture
theta_start = theta_emit; 
h_wb = waitbar(0, 'Generating image...');   % Progress indicator


%% START SCAN LOOP
info_on = 1;
for i=1:N_lines
	
    if(info_on);
        disp(['theta_emit: ' num2str(theta_emit*180/pi) ' degrees']);
    end
    
    % Set the focus for the transmit aperture in this direction
	xdc_focus (emit_aperture, 0, [focus(3)*sin(theta_emit) 0 focus(3)*cos(theta_emit)]);
    
    % Set the focus for the receive aperture in this direction
    if dynamic_receive_focus==1
        xdc_dynamic_focus (receive_aperture, 0, theta_receive, 0);
    else
        xdc_focus (receive_aperture, 0, [focus(3)*sin(theta_receive) focus(2) focus(3)*cos(theta_receive)]);
    end

    [v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
        
    % Store the result
    image_data(1:max(size(v)), i)=v;
    times(i) = t1;
  
    % Increment emit angle
    theta_emit = theta_emit + d_theta_emit;
    theta_receive = theta_emit;
    
    % Update waitbar
    waitbar(i/N_lines, h_wb);
end
close(h_wb)
theta_end = theta_receive-d_theta_emit;


%% ADD TIMECORRECTION
min_sample=round(min(times)*fs);
max_sample=round(max(times)*fs);
[n,m]=size(image_data);
n=n+(max_sample-min_sample);
rf_data=zeros(n, N_lines);

for i=1:N_lines
    rf_temp = [ zeros(round(times(i)*fs)-min_sample,1) ; image_data(:,i) ];
    rf_data(1:size(rf_temp, 1),i)=rf_temp;
end


%% GENERATE IMAGE FROM RF-DATA 
D = 5;              % Decimation factor
dyn = 40;           % Dynamic range
gain = 0;           % Image gain


%% Decimate
rf_d = rf_data(1:D:size(rf_data, 1), :);
fn= fs/D;           % New sampling frequency

% Plot 2D Fourier transform before detection
if show2dFourier
    da = d_theta_receive;   % angular beam spacing
    dr = 1/fn;              % Temporal spacing
    Nfftx = 128;            % Number of points in x dir for the FFT
    Nffty = 256;            % Number of points in y dir for the FFT
    
    plot2dFFT(rf_d, da, dr*1e6, Nfftx, Nffty);
    ylim([-6 6])         % Limit yaxis to -6MHz to 6 MHz
    xlabel('Spatial frequency [rad/s]');
    ylabel('Pulse frequency [MHz]');
    title('2D Fourier transform of RF data');
end



%% Display images
% Create axes
timetxrx = min(times) + (0:size(rf_d,1)-1)/fn;
[max_pulse_val, max_pulse_idx] = max(abs(hilbert(conv(excitation,conv(impulse_response,impulse_response)))));
timetxrx=timetxrx-max_pulse_idx/fs;     % Correct depth axis for duration of excitation and impulse response
disttxrx=(timetxrx)*c/2;
x_axis_angle = linspace(theta_start,theta_end, size(rf_d, 2));

% Plot image in beam space
figure;
h_ph_img = imagesc(x_axis_angle*180/pi, disttxrx*1000, rf_d);
% caxis([-dyn 0]-gain);
xlabel(['Beam angle [deg]'])
ylabel('Range [mm]')
colormap(gray(256));
axis on

% Plot scan converted image
angles_mat = ones(length(disttxrx), 1)*(x_axis_angle+pi/2);
ranges_mat = (disttxrx)'*ones(1, length(x_axis_angle));
[X_im, Y_im] = pol2cart(angles_mat, ranges_mat);

figure; 
pcolor(X_im*1000,Y_im*1000,rf_d); 
% caxis([-dyn 0]-gain);       % Set dynamic range and gain in image
shading interp; colormap(gray(256));
axis ij
axis equal;
xlabel('Azimuth [mm]');ylabel('Depth [mm]');