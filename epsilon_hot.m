function [epsilon,sigma,S,D,P,R,L] = epsilon_hot (f, amu, Z, B, density, T_eV, k_per, k_par)

harmonic_number = 3;
nu_omg = 0;

phys = constants();

c = phys.('c');
u0 = phys.('u0');
eps0 = phys.('eps0');
e = phys.('e');
amu0 = phys.('amu');

% get Z-function interpolant

persistent Z_interp

if isempty(Z_interp)  
    Z_interp = get_z_function_interpolant(); 
end

if abs(k_per) < 1e-5
    k_per = 1e-5*sign(k_per);
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

eta_hat = 0;
tau_hat = 0;
eps_hat = 0;

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
        
        this_Z_re = Z_interp.Z_re(x);
        this_Z_im = Z_interp.Z_im(x);
        this_Zp_re = Z_interp.Zp_re(x);
        this_Zp_im = Z_interp.Zp_im(x);
        
        Z = complex(this_Z_re,this_Z_im);
        Zp = complex(this_Zp_re,this_Zp_im);
        
        if x < Z_interp.min_zeta_re || x > Z_interp.max_zeta_re
            % Analytic limits off the end of the Z function table
            Z = complex(-1./x,0);
            Zp = complex(1/x.^2,0);
        end
        
%         [Z,Zp] = z_function(x,Zeta_C);
        
        In = besseli(n,lambda);
        Inp = besseli_prime(n,lambda);
        
        Ssum = Ssum + n.^2 ./ lambda .* In .* exp(-lambda) .* (-x0.*Z);
        Dsum = Dsum + n .* ( Inp - In ) .* exp(-lambda) .* (-x0 .* Z );
        Psum = Psum + In .* exp(-lambda) .* (x0 .* x .* Zp );
        
        eta_sum = eta_sum + n./lambda .* In .* exp(-lambda) .* (x0.^2 .* Zp );
        tau_sum = tau_sum + ( Inp - In ) .* exp(-lambda) .* (-x0 .* Z );
        eps_sum = eps_sum + ( Inp - In ) .* exp(-lambda) * (x0.^2 .* Zp );
        
    end
    
    Shat = Shat - wp.^2 ./ w.^2 .* Ssum;
    Dhat = Dhat + wp.^2 ./ w.^2 .* Dsum;
    Phat = Phat - wp.^2 ./ w.^2 .* Psum;
    
    eta_hat = eta_hat + wp.^2./(w.*wc) .* v_th.^2./c^2 .* eta_sum;
    tau_hat = tau_hat + wp.^2./wc.^2 .* v_th.^2./c^2 .* tau_sum;
    eps_hat = eps_hat + wp.^2./(w.*wc) .* v_th.^2./c.^2 .* eps_sum;
    
end

% Brambilla

eta_hat = -eta_hat/2.0;
tau_hat = -tau_hat/2.0;
eps_hat = +eps_hat/2.0;

exx = Shat;
exy = -i .* Dhat;
exz = n_par .* n_per .* eta_hat;

eyx = +i .* Dhat;
eyy = Shat - 2.*n_per.^2 .* tau_hat;
eyz = +i .* n_par .* n_per .* eps_hat;

ezx = n_par .* n_per .* eta_hat;
ezy = -i .* n_par .* n_per .* eps_hat;
ezz = Phat;

epsilon_r = zeros(3,3);
epsilon = complex(epsilon_r,0);
sigma = epsilon;

epsilon(1,1) = exx;
epsilon(1,2) = exy;
epsilon(1,3) = exz;

epsilon(2,1) = eyx;
epsilon(2,2) = eyy;
epsilon(2,3) = eyz;

epsilon(3,1) = ezx;
epsilon(3,2) = ezy;
epsilon(3,3) = ezz;

sigma = (epsilon - eye(3)) * w * eps0 / i;

end

function [Ip] = besseli_prime(n,zeta)

Ip = besseli(n-1,zeta) - n./zeta .* besseli(n,zeta);

end

