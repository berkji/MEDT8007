% Create phatom with gaussian distributed amplitudes
function [positions, amp] = ph_rect_gauss(N, dim, origo)

positions(:,1) = rand(N, 1)*dim(1) - dim(1)/2 + origo(1);
positions(:,2) = rand(N, 1)*dim(2) - dim(2)/2 + origo(2);
positions(:,3) = rand(N, 1)*dim(3) - dim(3)/2 + origo(3);
amp = ones(N,1);

