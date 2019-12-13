function [eps,sigma,S,D,P,R,L] = epsilon_hot (f, amu, Z, B, density, T_eV, k_per, k_par)

harmonic_number = 3;
nu_omg = 0;

assert(k_per >= 0);

phys = constants();

c = phys.('c');
u0 = phys.('u0');
eps0 = phys.('eps0');
e = phys.('e');
amu0 = phys.('amu');

if k_per lt 1e-5
    k_per = 1e-5;
end

w = 2*pi*f;
m = amu * amu0;
q = Z * e;
n_par = c * k_par / w;
n_per = c * k_per / w;

num_spec = numel(amu);

Shat = 1 + 0i;
Dhat = 0 + 0i;
Phat = 1 + 0i;

etahat = 0;
tauhat = 0;
epshat = 0;

for alp = 1:num_spec 
    
    wc = q(alp) * B / m(alp);
    wp = sqrt( density(alp) * q(alp)^2 / ( m(alp) * eps0 ) );
    v_th = sqrt( 2 * T_eV(alp) * e / m(alp) );
    
    lambda = k_per^2 * v_th^2 / (2 * wc^2);
    
    Ssum = 0;
    Dsum = 0;
    Psum = 0;
    
    eta_sum = 0;
    tau_sum = 0;
    eps_sum = 0;
    
    for n = -harmonic_number:+harmonic_number
        
        % Brambilla expressions, pg 254-255
        
        x = (w - n*wc) / (k_par * v_th);
        x0 = w / (k_par * v_th);
        
        Zeta_C = (nu_omg * w) / (k_par * v_th);
        [Z,Zp] = z_function(x,Zeta_C);
        
    end
    
end

end

function [Z, Zp] = z_function(x, zeta_C)
    
persistent zeta_re
persistent Z_re
persistent Z_im
persistent Zp_re
persistent Zp_im

if isempty(zeta_re) 
    
    z_function_table = 'zFunction.nc';
    ncid = netcdf.open(z_function_table);
    
    arg_re_id = netcdf.inqVarID(ncid,'arg_re');
    Z_re_id = netcdf.inqVarID(ncid,'Z_re');
    Z_im_id = netcdf.inqVarID(ncid,'Z_im');
    Zp_re_id = netcdf.inqVarID(ncid,'Zp_re');
    Zp_im_id = netcdf.inqVarID(ncid,'Zp_im');
    
    zeta_re = netcdf.getVar(ncid,arg_re_id);
    Z_re = netcdf.getVar(ncid,Z_re_id);
    Z_im = netcdf.getVar(ncid,Z_im_id);
    Zp_re = netcdf.getVar(ncid,Zp_re_id);
    Zp_im = netcdf.getVar(ncid,Zp_im_id);
    
end

end