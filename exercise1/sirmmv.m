function [h, t_0]=sirmmv(a, x_0, y_0, z_0, f_Sample, c_0)
% [h, t_0]=sirmmv(a, x_0, y_0, z_0, f_Sample, c_0)
%
% Time domain spatial impulse response (SIR) from a circular (plane)
% aperture. Where:
% h is the SIR in the time domain;
% t_0 is the start time of the SIR (i.e. time of the first sample);
% a is the radius of the aperture with its center at the origin and
% extending in to the x-y plane in [m]; 
% x_0, y_0 and z_0 are the Cartesian coordinates of the observation point
% in [m]; 
% f_Sample is the sample frequency in [Hz] and
% c_0 is the speed of sound in [m/s].
%
%
% author:  Marco M. Voormolen
% draft:   25 March 2008

% update:  
%
% uses:   
% sub of:  


% Parameters
%r=sqrt(x_0^2 + y_0^2 + z_0^2);
r_a=sqrt(x_0^2 + y_0^2);
R_1=sqrt(z_0^2 + (a - r_a)^2);
R_2=sqrt(z_0^2 + (a + r_a)^2);
if a>r_a
  t_SIR=(floor(z_0/c_0*f_Sample) - 2:1:ceil(R_2/c_0*f_Sample) + 1)/f_Sample;
else
  t_SIR=(floor(R_1/c_0*f_Sample) - 2:1:ceil(R_2/c_0*f_Sample) + 1)/f_Sample;
end
t_0=t_SIR(1);


% Spatial impluse response
h=zeros(size(t_SIR));
tic
for i=1:length(t_SIR)
  if (a>r_a && c_0*t_SIR(i)<z_0) || (a<=r_a && c_0*t_SIR(i)<R_1) || c_0*t_SIR(i)>R_2
    h(i)=0;
  else
    if c_0*t_SIR(i)>=r_a && c_0*t_SIR(i)<=R_1
      h(i)=c_0;
    else
      h(i)=c_0/pi*acos((r_a^2 + c_0^2*t_SIR(i)^2 - z_0^2 - a^2)/(2*r_a*sqrt(c_0^2*t_SIR(i)^2 - z_0^2)));
    end  
  end
end
% Directly from Szabo, page 174: Three errors, the 'r' in 'ct<R_1 for a<r'
% should have been 'r_a', the 'r' in 'r_a<r<R_1' should have been 'ct' and
% 'r_0' is 'r_a'. Cobbold, page 151: No errors, 1-0 for Cobbold!
