function p = case1()

% run with
%
% cold only (via quadratic)
% dispersioneering(@case1,'num_points',30,'use_root_finder',true,'use_cold_eps',false)
%
% cold only (via quadratic and root finder - both using cold eps)
% dispersioneering(@case1,'num_points',30,'use_root_finder',true,'use_cold_eps',eps)
%
% cold and hot
% dispersioneering(@case1,'num_points',30,'use_root_finder',true,'use_cold_eps',false)

phys = constants();

c = phys.('c');
me_amu = phys.('me_amu');
zi = complex(0,1);

% cold parameters

f = 7.5e6;
amu = {me_amu, 2};
Z   = {-1,1};

B_func = @(x) x.*0 + 1.2;

den = @(x) 10.^(x*10+12);

den_m3_func = {den, den};

% default case - hot only parameters

T_eV_func = {@(x) x.*0 + 1e-3, @(x) x.*0 + 10};
k_par = 20;

% construct parameters object

p = PARAMETERS(f,amu,Z,B_func,den_m3_func,T_eV_func,k_par);

end