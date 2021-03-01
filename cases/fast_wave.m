function p = fast_wave()

% run with
%
% cold only (via quadratic)
% dispersioneering(@fast_wave,'num_points',300,'use_root_finder',false,'use_cold_eps',true)
%
% cold only (via quadratic and root finder - both using cold eps)
% dispersioneering(@fast_wave,'num_points',30,'use_root_finder',true,'use_cold_eps',eps)
%
% cold and hot
% dispersioneering(@fast_wave,'num_points',30,'use_root_finder',true,'use_cold_eps',false)

phys = constants();

c = phys.('c');
me_amu = phys.('me_amu');
zi = complex(0,1);

% cold parameters

f = 52e6;
amu = {me_amu, 1};
Z   = {-1,1};

B_func = @(x) x.*0 + 2.5; % T

den = @(x) 10.^(x*4+16);

den_m3_func = {den, den};

% default case - hot only parameters

T_eV_func = {@(x) x.*0 + 1e-3, @(x) x.*0 + 10};
k_par = 20;

% construct parameters object

p = PARAMETERS(f,amu,Z,B_func,den_m3_func,T_eV_func,k_par);

end