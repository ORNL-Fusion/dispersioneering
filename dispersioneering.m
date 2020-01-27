function dispersioneering()

phys = constants();

c = phys.('c');
me_amu = phys.('me_amu');
zi = complex(0,1);

num_points = 1;
x = linspace(0,1,num_points);

% Don's case 1
f = 7.5e6;
k_par = 20;
amu = [me_amu, 2];
Z   = [-1,1];
T_eV = [500,500];
B = x.*0 + 1.2;
den = x.*0 + 1e20;
n = [den; den];
% 
% f = 7.5e6;
% k_par = 20;
% amu = [me_amu, 2];
% Z   = [-1,1];
% T_eV = [500,500];
% B = x.*0 + 1.2;
% den = 10.^linspace(19,20,num_points);
% n = [den; den];

% % Case A
% f = 7.5e6;
% k_par = 20;
% amu = [me_amu, 2];
% Z   = [-1,1];
% T_eV = [50,50];
% B = x.*0 + 1.2;
% den = 10.^linspace(18,19,num_points);
% n = [den; den];

w = 2*pi*f;
n_par = k_par * c / w;

% Produce 2D space over k_per and B

num_points = numel(B);

% solve the cold plasma dispersion relation via the quadratic formula

A1 = zeros(num_points,1);
B1 = zeros(num_points,1);
C1 = zeros(num_points,1);

for i=1:num_points
    
    [eps,sigma,S,D,P,R,L] = epsilon_cold(f, amu, Z, B(i), n(:,i));
    
    % From pg 177 Brambilla
    
    A1(i) = S;
    B1(i) = -1.*(R.*L + P.*S - n_par.^2 .* (P+S)); % note the -sign here
    C1(i) = P.*(n_par.^2-R).*(n_par.^2-L);
    
end

n_per_1 = +sqrt((-B1 + sqrt(B1.^2-4.*A1.*C1))./(2.*A1));
n_per_2 = +sqrt((-B1 - sqrt(B1.^2-4.*A1.*C1))./(2.*A1));
n_per_3 = -sqrt((-B1 + sqrt(B1.^2-4.*A1.*C1))./(2.*A1));
n_per_4 = -sqrt((-B1 - sqrt(B1.^2-4.*A1.*C1))./(2.*A1));

% root finder approach for general hot (or cold) plasma

% k_per range

num_init_k_per = 5;

k_per_min = -1000;
k_per_max = +1000;

n_per_out = complex(zeros(num_points,num_init_k_per^2),zeros(num_points,num_init_k_per^2));

initial_k_pers = linspace(k_per_min,k_per_max,num_init_k_per);

if is_octave()
    options = optimset('Jacobian','false',...
        'TolX',1e-3);
else
    options = optimoptions('fsolve',...
        'Display','iter',...
        'SpecifyObjectiveGradient',false,...
        'CheckGradients',false,...
        'UseParallel',false,...
        'StepTolerance',1e-6,...
        'FunctionTolerance',1e-6,...
        'OptimalityTolerance',1e-6);
end

params.w = w;
params.c = c;
params.f = f;
params.amu = amu;
params.Z = Z;
params.T_eV = T_eV;
params.k_par = k_par;
params.n_par = n_par;

for i=1:num_points
    
    tic;
    B0 = B(i);
    n0 = n(:,i);
    
    % scan the [re,im] space of initial guesses of k_per for the root
    % solver
    
    cnt = 1;
    for k = 1:num_init_k_per
        for j = 1:num_init_k_per
            
            n_per_init_re = initial_k_pers(k) * c / w;
            n_per_init_im = initial_k_pers(j) * c / w;
            
            n_per_init = n_per_init_re + n_per_init_im*zi;
            
            params.B0 = B0;
            params.n0 = n0;
            
            [n_per_out(i,cnt),fval,exitflag,output] = ...
                fsolve(@(x) det_fun(x,params),n_per_init,options);
            
            cnt = cnt + 1;
            
        end
    end
        
    disp(['point ', num2str(i), ' of ', num2str(num_points)]);
    disp(['calculation time : ', num2str(toc)]);
    disp(['hot n_per_out : ', num2str(n_per_out(i))]);
    disp(['cold n_per_1 : ', num2str(n_per_1(i))]);
    disp(['cold n_per_2 : ', num2str(n_per_2(i))]);
    disp(['cold n_per_3 : ', num2str(n_per_3(i))]);
    disp(['cold n_per_4 : ', num2str(n_per_4(i))]);
    disp(['  ']);
    
end

% plot k_per versus the location (B,n)

num_species = numel(Z);

figure()
subplot(3,1,1)
plot(x,B);
subplot(3,1,2)
for s=1:num_species
    semilogy(x,n(s,:));
    hold on
end
hold off
subplot(3,1,3)
plot(x,real(n_per_1).*w./c,'-','Color','black');
ylim([k_per_min,k_per_max]);
hold on
plot(x,real(n_per_2).*w./c,'-','Color','black');
plot(x,real(n_per_3).*w./c,'-','Color','black');
plot(x,real(n_per_4).*w./c,'-','Color','black');
plot(x,imag(n_per_1).*w./c,'-r');
plot(x,imag(n_per_2).*w./c,'-r');
plot(x,imag(n_per_3).*w./c,'-r');
plot(x,imag(n_per_4).*w./c,'-r');
for B1=1:numel(n_per_out(1,:))
    p = plot(x,real(n_per_out(:,B1)*w/c),'o','Color','black','LineWidth',2);
    p = plot(x,imag(n_per_out(:,B1)*w/c),'o','Color','red',  'LineWidth',2);
end
hold off

% figure()
% semilogx(den,real(n_per_sq1)*w/c,'-b');
% hold on
% plot(den,real(n_per_sq2)*w/c,'-b');
% plot(den,real(n_per_sq3)*w/c,'-b');
% plot(den,real(n_per_sq4)*w/c,'-b');
% plot(den,imag(n_per_sq1)*w/c,'-r');
% plot(den,imag(n_per_sq2)*w/c,'-r');
% plot(den,imag(n_per_sq3)*w/c,'-r');
% plot(den,imag(n_per_sq4)*w/c,'-r');
% for bb=1:numel(n_per_out(1,:))
%     p = plot(den,real(n_per_out(:,bb)*w/c),'o','Color','black','LineWidth',2);
%     p = plot(den,imag(n_per_out(:,bb)*w/c),'o','Color','red',  'LineWidth',2);
% end
% hold off
% ylim([k_per_min,k_per_max]);
% symlog();
% disp([' ']);

end

function [det,det_gradient] = det_fun(n_per,params)

w       = params.w;
c       = params.c;
f       = params.f;
amu     = params.amu;
Z       = params.Z;
B0      = params.B0;
n0      = params.n0;
T_eV    = params.T_eV;
k_par   = params.k_par;
n_par   = params.n_par;

use_cold_eps = false; % for testing general determinant with the cold dielectric
use_cold_det = false; % for testing the root solver with the quadratic determinant

if use_cold_eps % Cold dielectric
    
    [eps,sigma,S,D,P,R,L] = epsilon_cold(f, amu, Z, B0, n0);
    
else % Hot dielectric
    
    k_per = n_per .* w / c;
    [eps] = epsilon_hot(f, amu, Z, B0, n0, T_eV, k_per, k_par);
    
end

if use_cold_eps && use_cold_det
    
    % Quadratic form for determinant (works only for cold)
    % from pg 177 of Brambilla
    
    A0 = S;
    B0 = R*L + P*S - n_par^2 * (P+S);
    C0 = P*(n_par^2-R)*(n_par^2-L);
    
    det = A0 .* n_per^4 - B0 .* n_per^2 + C0;
    
else
    
    % Generalized determinant (works for hot or cold)
    
    exx = squeeze(eps(1,1,:));
    exy = squeeze(eps(1,2,:));
    exz = squeeze(eps(1,3,:));
    
    eyx = squeeze(eps(2,1,:));
    eyy = squeeze(eps(2,2,:));
    eyz = squeeze(eps(2,3,:));
    
    ezx = squeeze(eps(3,1,:));
    ezy = squeeze(eps(3,2,:));
    ezz = squeeze(eps(3,3,:));
    
    kx = n_per .* w./c;
    kz = n_par .* w./c;
    k0 = w/c;
    
    sz = size(exx);
    assert(numel(n_per)==numel(exx));
    kx = reshape(kx,sz);
    
    det = -(ezy.*k0.^4.*((exz.*eyx - exx.*eyz).*k0.^2 + eyx.*kx.*kz + eyz.*kz.^2)) + ...
        (-(ezx.*k0.^2) - kx.*kz).* ...
        (exy.*eyz.*k0.^4 + exz.*k0.^2.*(-(eyy.*k0.^2) + kx.^2 + kz.^2) + ...
        kx.*kz.*(-(eyy.*k0.^2) + kx.^2 + kz.^2)) + ...
        (-(ezz.*k0.^2) + kx.^2).* ...
        (-(exy.*eyx.*k0.^4) + (-(exx.*k0.^2) + kz.^2).* ...
        (-(eyy.*k0.^2) + kx.^2 + kz.^2));
    
    if nargout > 1
        
        % specific the gradient of the determinant also
        %
        % note that this does not appear to work well so it is presently
        % disabled via the 'SpecifyObjectiveGradient',false option to
        % fsolve
        
        det_gradient = -((k0.^2.*w.*(c.^3.*kz.*(exy.*eyz.*k0.^2 - eyy.*ezx.*k0.^2 +...
            eyx.*ezy.*k0.^2 + ezx.*kz.^2 + exz.*(-(eyy.*k0.^2) + kz.^2)) + ...
            2.*c.^2.*(exy.*eyx.*k0.^2 + exz.*ezx.*k0.^2 + ezz.*kz.^2 +...
            exx.*(-(eyy.*k0.^2) - ezz.*k0.^2 + kz.^2)).*n_per.*w + ...
            3.*c.*(exz + ezx).*kz.*n_per.^2.*w.^2 + 4.*exx.*n_per.^3.*w.^3))/c.^4);
    end
    
end

end


