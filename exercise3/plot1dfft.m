function plot1dfft(fn_y, rf_d)
rf_freq = fft(rf_d);
figure;
plot(fn_y/1e6, abs(rf_freq));
