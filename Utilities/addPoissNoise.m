function [ p_noisy ] = addPoissNoise( p, N_0 )
%ADDPOISSNOISE add Poisson noise to the projections
%   p projection data
%   N_0 photon number. Nomally, 10^5 will be low-dose and 10^6 will be
%   normal dose
%   p_noisy projection data with noise

yi = poissrnd(N_0*exp(-p));             %   yi the measurement along the ith X-ray path
p_noisy = -log((yi/N_0)+0.0000001);  
idx = find(p_noisy(:) < 0);
p_noisy(idx) = 0;
end

