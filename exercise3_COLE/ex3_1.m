clc
close all
clear

% -- phantom
n_scatt = 20;
start_depth = 30/1000; 
stop_depth = 90/1000; 

% -- probe
axial_fov = [20 100]/1000; 
lateral_fov = [-7 +7]/1000; 
Number_of_lines = 30; 
FWHM = 2/1000;

f0 = 2.5e6;
fs = 50e6;
c0 = 1540;

% -- psf 
n_oscill = 3;
t = 0:1/fs:n_oscill/f0;
t = t - mean(t);
s_y = max(t)/3; % in sec!!
s_x = FWHM / (2*sqrt(2*log(2)));
P_y = exp(-.5.*(t/s_y).^2).*cos(2*pi*f0*t + pi/6);
P_x = @(x) exp(-.5.*(x/s_x).^2);

% -- sampled x and y axis
dr = c0 / (2*fs);
y_axis = axial_fov(1):dr:axial_fov(2);
x_axis = linspace(lateral_fov(1),lateral_fov(2),Number_of_lines); % scan lines

subplot(1,2,1)
plot(t,P_y), title('axial psf'), axis tight
xlabel('time (s)')
x_plot = linspace(-2*s_x,2*s_x,128); % for plot only
subplot(1,2,2)
plot(x_plot,P_x(x_plot)), title('axial psf'), axis tight
xlabel('space (m)')

% -- tissue
y_scatt = linspace(start_depth,stop_depth,n_scatt);
x_scatt = zeros(size(y_scatt));
amp = 1;

figure
plot(zeros(size(y_axis)),y_axis,'b.'), hold on
plot(x_scatt,y_scatt,'ko','markerfacecolor','g'),
ylim(axial_fov)
xlim(lateral_fov)
title('reflectivity')

% -- image
rf_image = zeros(numel(y_axis),numel(x_axis));
cont = 0;
for ii = 1:Number_of_lines
    
    pos_x = x_axis(ii); % -- scan line position

    % -- prefilter (don't consider scatterers far from the scan line)
    isin = abs(x_scatt - pos_x) < 3*s_x;
    
    if sum(isin) > 0
        cont = cont + 1;
        refl = zeros(size(y_axis)); % reflectivity
        
        % -- project scatterers on sampled y axis
        y_scatt_valid = y_scatt(isin);
    
        y_scatt_idx = (y_scatt_valid - axial_fov(1)) / diff(axial_fov) * numel(y_axis-1) + 1;  
        y_scatt_idx_f = floor(y_scatt_idx);
        y_scatt_idx_c = ceil(y_scatt_idx);
        
        w_f = 1 - abs(y_scatt_idx - y_scatt_idx_f);
        w_c = 1 - abs(y_scatt_idx - y_scatt_idx_c);
        
        % -- weights
        x_scatt_valid = x_scatt(isin);
        w_amp = P_x(abs(x_scatt_valid - pos_x)); 
        
        refl(y_scatt_idx_f) = w_amp.*w_f;
        refl(y_scatt_idx_c) = w_amp.*w_c;
        
        % -- convolution
        rf_image(:,ii) = conv(refl,P_y,'same');
        
        if cont >= 1
            figure
            stem(y_axis,refl), hold on;
            plot(y_scatt,ones(size(y_scatt)),'ko','markerfacecolor','g');
            
            figure
            plot(y_axis,rf_image(:,ii)),
            
        end
            
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

