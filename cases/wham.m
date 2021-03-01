function p = wham()

% run with
%
% cold only (via quadratic)
% dispersioneering(@wham,'num_points',300,'use_root_finder',false,'use_cold_eps',true)
%
% cold only (via quadratic and root finder - both using cold eps)
% dispersioneering(@wham,'num_points',30,'use_root_finder',true,'use_cold_eps',eps)
%
% cold and hot
% dispersioneering(@wham,'num_points',30,'use_root_finder',true,'use_cold_eps',false,'kper_min',-100,'kper_max',100,'num_init_kper',5)

phys = constants();

c = phys.('c');
me_amu = phys.('me_amu');
zi = complex(0,1);

% cold parameters

f = 27.1e6;
amu = {me_amu, 2}; % D plasma?
Z   = {-1,1};

B_func = @(x) x.*0 + 0.8; % T

den = @(x) 10.^(x*2+18); % 10.^(x*span+min) where span and min are powers of 10

den_m3_func = {den, den};

% default case - hot only parameters

T_eV_func = {@(x) x.*0 + 600, @(x) x.*0 + 600}; % 100ev electrons, cold ions (600 ev is also an option)
k_par = 15*2*pi*f/c; % n_par = 15

% construct parameters object

p = PARAMETERS(f,amu,Z,B_func,den_m3_func,T_eV_func,k_par);

end