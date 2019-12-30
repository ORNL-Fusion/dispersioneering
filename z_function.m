function [Z, Zp] = z_function(zeta_re, zeta_re_C)

% zeta_re   = (w +/- wc) /(|kparalllel|*vThermal)
% zeta_re_C = nuCollision/(|kparalllel|*vThermal)

persistent Z_interp

% Read Z function from file if not done so already

if isempty(Z_interp)
    
    Z_interp = get_z_function_interpolant();
    
end

this_Z_re = Z_interp.Z_re(zeta_re);
this_Z_im = Z_interp.Z_im(zeta_re);
this_Zp_re = Z_interp.Zp_re(zeta_re);
this_Zp_im = Z_interp.Zp_im(zeta_re);

Z = complex(this_Z_re,this_Z_im);
Zp = complex(this_Zp_re,this_Zp_im);

N = numel(zeta_re);
for n = 1:N
    if zeta_re(n) < Z_interp.min_zeta_re || zeta_re(n) >Z_interp.max_zeta_re
        % Analytic limits off the end of the Z function table
        Z(n) = complex(-1./zeta_re,0);
        Zp(n) = complex(1/zeta_re.^2,0);
    end
end

factor = complex(1,0) - complex(0,zeta_re_C).*Z;
Z = Z ./ factor;
Zp = Zp ./ (factor.*factor);

end