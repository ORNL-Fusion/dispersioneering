function [eps,sigma,S,D,P,R,L] = epsilon_cold (f, amu, Z, B, n)

num_species = numel(amu);

assert(numel(Z)==num_species);
assert(numel(n)==num_species);

phys = constants();

c = phys.('c');
u0 = phys.('u0');
eps0 = phys.('eps0');
e = phys.('e');
amu0 = phys.('amu');

validateattributes(B,{'numeric'},{'nonnegative'});

w0 = 2 * pi * f;
w = w0 * complex( 1, 0);
m = amu * amu0;
q = Z * e;



% Swanson pg 23 - 24
% ------------------

eps_swan = q./abs(q);
wc_swan = abs(q) * B ./ m;

S = 1;
D = 0;
P = 1;

for s=1:num_species
    
    wp = sqrt( n(s) * q(s)^2 / (m(s) * eps0) );
    
    S = S - wp^2 / (w^2 - wc_swan(s)^2);
    D = D + eps_swan(s) * wc_swan(s) * wp^2 / (w * (w^2 - wc_swan(s)^2) );
    P = P - wp^2 / w^2;
    
end

K1 = S;
K2 = D/1i;
K3 = P;

R = S+D;
L = S-D;

epsilon_swan = complex(zeros(3,3));

epsilon_swan(1,1) = +K1;
epsilon_swan(1,2) = +K2;
epsilon_swan(1,3) = 0;

epsilon_swan(2,1) = -K2;
epsilon_swan(2,2) = +K1;
epsilon_swan(2,3) = 0;

epsilon_swan(3,1) = 0;
epsilon_swan(3,2) = 0;
epsilon_swan(3,3) = +K3;


sigma_swan =  ( epsilon_swan - eye(3) ) * w0 * eps0 / 1i;

sigma = sigma_swan;
eps = epsilon_swan;

end
