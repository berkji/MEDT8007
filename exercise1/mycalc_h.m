function [sir_t, t_start]=mycalc_h(no_elements, width, kerf, f_c, focus_point, observation_point)
x_p = ((width+kerf)*0.5*no_elements - kerf/2);
x_m = -1*x_p;
x_1 = linspace(x_m,x_p,no_elements); % x-coordinate of the (candidate) aperture point(s) in m
y_1 = 0e-3;       % y-coordinate of the (candidate) aperture point(s) in m
x_a = x_p;
f_Sample=4*f_c; % sample frequency in Hz
N_FFT=1024*1; % length of the FFT

c_0=1500; % speed of sound in m/s
df=f_Sample/N_FFT; % frequency resolution in Hz
f_cc = 0:df:floor(N_FFT/2)*df; % 0.1e6:0.1e5:6e6; 
f=f_cc; % frequency (range) Hz

x_f = focus_point(1);
y_f = focus_point(2);
z_f = focus_point(3);

x_0 = observation_point(1);
y_0 = observation_point(2);
z_0 = observation_point(3);

R=zeros(length(x_1), length(y_1));
H=zeros(length(x_0), length(f)); % complex amplitude function

for m=1:length(x_1)
  for n=1:length(y_1)
    for q=1:length(x_0)
      if sqrt(x_1(m)^2 + y_1(n)^2)<= x_a % only accept points that lie within the aperture
        %<frequency independent variables>
		xf = [x_f 0 z_f];
        xp = [x_1(m) 0 0];
        xpa =[x_a 0 0];
        xma =[-x_a 0 0];
        x0 = [x_0(q) 0 z_0];
        r = x0 - xp;
        R = xf - xma;
        P = xf - xp;
        rpf(m,n) = sqrt(sum(R.*R)) - sqrt(sum(P.*P));
        rr(m,n)  = sqrt(sum(r.*r));
        for i=1:length(f)
		  k(i) = 2*pi*f(i)/c_0;
          H(q, i)=H(q, i) + exp(-1i*k(i)*(rpf(m,n)+rr(m,n)))/(2*pi*rr(m,n));
%         H(q, i)=H(q, i) + exp(-1i*k(i)*(rr(m,n)))/(2*pi*rr(m,n));
        end
      end  
    end  
  end  
end

n=find(rpf==0);
rpf(n)=NaN;

n=find(rr==0);
rr(n)=NaN;


n_b = floor(min(rpf(:) + rr(:))/c_0*f_Sample)-1;
t_start = n_b*(1/f_Sample);
sir_t = real(ifft(H(1,:), N_FFT));




