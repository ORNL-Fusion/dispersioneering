function output = dispersioneering_run(opts,params)

set(0,'defaultAxesFontSize',20);

phys = constants();
c = phys.('c');
me_amu = phys.('me_amu');
zi = complex(0,1);

x = linspace(0,1,opts.num_points);

B = params.B_func(x);

num_species = numel(params.Z);

for s=1:num_species
    T_eV(s,:) = params.T_eV_func{s}(x);
    n(s,:) = params.den_m3_func{s}(x);
end

w = 2*pi*params.f;
n_par = params.k_par * c / w;

amu = cell2mat(params.amu);
Z = cell2mat(params.Z);

num_points = numel(B);



%% Solve the cold plasma dispersion relation via the quadratic formula

A1 = zeros(num_points,1);
B1 = zeros(num_points,1);
C1 = zeros(num_points,1);

for ii=1:num_points
    
    [eps,sigma,S,D,P,R,L] = epsilon_cold(params.f, amu, Z, B(ii), n(:,ii));
    
    % From pg 177 Brambilla
    
    A1(ii) = S;
    B1(ii) = -1.*(R.*L + P.*S - n_par.^2 .* (P+S)); % note the -sign here
    C1(ii) = P.*(n_par.^2-R).*(n_par.^2-L);    
    
end

n_per_1 = +sqrt((-B1 + sqrt(B1.^2-4.*A1.*C1))./(2.*A1));
n_per_2 = +sqrt((-B1 - sqrt(B1.^2-4.*A1.*C1))./(2.*A1));
n_per_3 = -sqrt((-B1 + sqrt(B1.^2-4.*A1.*C1))./(2.*A1));
n_per_4 = -sqrt((-B1 - sqrt(B1.^2-4.*A1.*C1))./(2.*A1));



%% Plot the cold plasma k_per

figs.profiles = figure();
subplot(3,1,1)
plot(x,B);
ylabel('B(T)');
subplot(3,1,2)
for s=1:num_species
    plot(x,n(s,:));
    hold on
end
hold off
ylabel('Density [m^-3]');
subplot(3,1,3)
plot(x,real(n_per_1).*w./c,'-','Color','black','LineWidth',2);
%ylim([kper_min,kper_max]);
hold on
plot(x,(real(n_per_2)).*w./c,'-','Color','black','LineWidth',2);
plot(x,(real(n_per_3)).*w./c,'-','Color','black','LineWidth',2);
plot(x,(real(n_per_4)).*w./c,'-','Color','black','LineWidth',2);
plot(x,(imag(n_per_1)).*w./c,'-r','LineWidth',2);
plot(x,(imag(n_per_2)).*w./c,'-r','LineWidth',2);
plot(x,(imag(n_per_3)).*w./c,'-r','LineWidth',2);
plot(x,(imag(n_per_4)).*w./c,'-r','LineWidth',2);

ylabel('k_{per} [m^-1]');
title('black=real, red=imag, solid=cold, symbols=root finder');



%% Optinally also use a root finder approach for general hot (or cold) plasma

if opts.use_root_finder
    
    % k_per range
    
    num_init_kper = opts.num_init_kper;
    
    kper_min = opts.kper_min;
    kper_max = opts.kper_max;
    
    n_per_out = complex(zeros(num_points,num_init_kper^2),zeros(num_points,num_init_kper^2));
    
    initial_k_pers = linspace(kper_min,kper_max,num_init_kper);
    
    if is_octave()
        options = optimset('Jacobian','false',...
            'TolX',1e-3);
    else
        options = optimoptions('fsolve',...
            'Display','iter-detailed',...
            'SpecifyObjectiveGradient',false,...
            'CheckGradients',false,...
            'UseParallel',false,...
            'StepTolerance',1e-10,...
            'FunctionTolerance',1e-10,...
            'OptimalityTolerance',1e-10,...
            'Diagnostics','on',...
            'FiniteDifferenceType','central');
    end
    
    point_params.w = w;
    point_params.c = c;
    point_params.f = params.f;
    point_params.amu = amu;
    point_params.Z = Z;
    point_params.T_eV = T_eV;
    point_params.k_par = params.k_par;
    point_params.n_par = n_par;
    point_params.opts = opts;
    
    for ii=1:num_points
        
        tic;
        B0 = B(ii);
        n0 = n(:,ii);
        
        % Create plot data for the scan of the [re,im] space of initial
        % guesses of k_per for the root solve for a particular point
        look_at_point = 9;
        if ii==look_at_point
            num_init_kper2 = 200;
            
            kper_min2 = opts.kper_min;
            kper_max2 = opts.kper_max;
            
            initial_k_pers2 = linspace(kper_min2,kper_max2,num_init_kper2);
            
            cnt = 1;
            for k = 1:num_init_kper2
                for j = 1:num_init_kper2
                    
                    n_per_init_re = initial_k_pers2(k) * c / w;
                    n_per_init_im = initial_k_pers2(j) * c / w;
                    
                    n_per_init(1) = n_per_init_re;
                    n_per_init(2) = n_per_init_im;
                    
                    point_params.B0 = B0;
                    point_params.n0 = n0;
                    
                    det(k,j) = det_fun(n_per_init,point_params);
                    
                    cnt = cnt + 1;
                    
                end
            end          
        end
        
        % Call the root finder for each point in the sampled (re,im) kper
        % space
        cnt = 1;
        for k = 1:num_init_kper
            for j = 1:num_init_kper
                
                n_per_init_re = initial_k_pers(k) * c / w;
                n_per_init_im = initial_k_pers(j) * c / w;
                
                kper_re_start(k,j) = initial_k_pers(k);
                kper_im_start(k,j) = initial_k_pers(j);
                                
                point_params.B0 = B0;
                point_params.n0 = n0;
                
                use_2D = true;             
                if use_2D % 2D using re and im as seperate axes
                    n_per_init(1) = n_per_init_re;
                    n_per_init(2) = n_per_init_im;
                    [n_per_out(ii,1:2,cnt),fval,exitflag,output] = fsolve(@(x) det_fun(x,point_params),n_per_init,options);
                else % 1D using [re,im] as a single input
                    n_per_init = n_per_init_re + n_per_init_im*zi;
                    [n_per_out(i,cnt),fval,exitflag,output] = fsolve(@(x) det_fun(x,point_params),n_per_init,options);
                end
                cnt = cnt + 1;
                
            end
        end

        % Just plot up a [re,im] space of the determinant to see if the
        % roots are found where they should be. 
        if ii==look_at_point
            figure
%             contour(initial_k_pers2,initial_k_pers2,log10(abs(det)))
              contour(initial_k_pers2,initial_k_pers2,det)
            hold on
            if numel(size(n_per_out))>2
                rs_kxre = n_per_out(ii,1,:)*w/c;
                rs_kxim = n_per_out(ii,2,:)*w/c;
            else
                rs_kxre = real(n_per_out(ii,:))*w/c;
                rs_kxim = imag(n_per_out(ii,:))*w/c;
            end
            scatter(rs_kxim,rs_kxre,'SizeData',150,'LineWidth',5)
            scatter(kper_im_start(:),kper_re_start(:));
            xlabel('im(kper)')
            ylabel('re(kper)')
            disp('here');
        end
        
        disp(['point ', num2str(ii), ' of ', num2str(num_points)]);
        disp(['calculation time : ', num2str(toc)]);
        disp(['hot n_per_out : ', num2str(n_per_out(ii))]);
        disp(['cold n_per_1 : ', num2str(n_per_1(ii))]);
        disp(['cold n_per_2 : ', num2str(n_per_2(ii))]);
        disp(['cold n_per_3 : ', num2str(n_per_3(ii))]);
        disp(['cold n_per_4 : ', num2str(n_per_4(ii))]);
        disp(['  ']);
        
    end
    
    % overplot the root finder solution
    figure(figs.profiles);
    hold on
    if numel(size(n_per_out))>2
        for ii=1:numel(n_per_out(1,1,:))
            p = plot(x,n_per_out(:,1,ii)*w/c,'o','Color','black');
            p = plot(x,n_per_out(:,2,ii)*w/c,'o','Color','red');
        end
    else
        
        for ii=1:numel(n_per_out(1,:))
            p = plot(x,real(n_per_out(:,ii)*w/c),'o','Color','black');
            p = plot(x,imag(n_per_out(:,ii)*w/c),'o','Color','red');
        end
    end
    yline(opts.kper_min)
    yline(opts.kper_max)
%     ylim([kper_min,kper_max]);
    
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
% ylim([kper_min,kper_max]);
% symlog();
% disp([' ']);

output = 0;

end

%% Define the determinant for which the root solver will try to find the roots of.

function [det_out,det_gradient] = det_fun(n_per_in,params)


if numel(n_per_in) > 1
    n_per = n_per_in(1) + 1i*n_per_in(2); % the n_per_in is two real numbers (re and im parts) 
else
    n_per = n_per_in; % the n_per_in is a complex number 
end

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
opts    = params.opts;

use_cold_eps = opts.use_cold_eps; % for testing general determinant with the cold dielectric
use_cold_det = opts.use_cold_det; % for testing the root solver with the quadratic determinant


% get dielectric

if use_cold_eps % Cold dielectric
    
    [eps,sigma,S,D,P,R,L] = epsilon_cold(f, amu, Z, B0, n0);
    
else % Hot dielectric
    
    k_per = n_per .* w / c;
    [eps] = epsilon_hot(f, amu, Z, B0, n0, T_eV, k_per, k_par);
    
end

exx = eps(1,1);
exy = eps(1,2);
exz = eps(1,3);

eyx = eps(2,1);
eyy = eps(2,2);
eyz = eps(2,3);

ezx = eps(3,1);
ezy = eps(3,2);
ezz = eps(3,3);

kx = n_per .* w./c; % kx==k_per
kz = n_par .* w./c; % kz==k_par
k0 = w/c;
    

% compute determinant of wave equation    

if use_cold_eps && use_cold_det
    
    % det approach 1 - Quadratic form for determinant (works only for cold)
    % from pg 177 of Brambilla
    
    A0 = S;
    B0 = R*L + P*S - n_par^2 * (P+S);
    C0 = P*(n_par^2-R)*(n_par^2-L);
    
    det_out1 = A0 .* n_per^4 - B0 .* n_per^2 + C0;
    
    % det approach 2 - do the det of abs(matrix) in matlab approach
    
    mat = [...
        -exx*k0^2+kz^2,           -exy*k0^2,  -exz*k0^2-kx*kz; ...
        -eyz*k0^2,      -eyy*k0^2+kx^2+kz^2,        -eyz*k0^2; ...
        -ezx*k0^2-kx*kz,          -ezy*k0^2,    -ezz*k0^2+kx^2  ...
        ];
    
    det_out2 = det(mat);
    
    % det approach 3 - do the det in mathematica in dispersioneering.nb
    
    det_out3 = -(ezy.*k0.^4.*((exz.*eyx - exx.*eyz).*k0.^2 + eyx.*kx.*kz + eyz.*kz.^2)) + ...
        (-(ezx.*k0.^2) - kx.*kz).* ...
        (exy.*eyz.*k0.^4 + exz.*k0.^2.*(-(eyy.*k0.^2) + kx.^2 + kz.^2) + ...
        kx.*kz.*(-(eyy.*k0.^2) + kx.^2 + kz.^2)) + ...
        (-(ezz.*k0.^2) + kx.^2).* ...
        (-(exy.*eyx.*k0.^4) + (-(exx.*k0.^2) + kz.^2).* ...
        (-(eyy.*k0.^2) + kx.^2 + kz.^2));
        
    det_out = det_out1;
    
else
    
    % Generalized determinant (works for hot or cold)
    % here we calculate the det in mathematica and just evaluate in matlab
    
    sz = size(exx);
    assert(numel(n_per)==numel(exx));
    kx = reshape(kx,sz);
    
    det_out = -(ezy.*k0.^4.*((exz.*eyx - exx.*eyz).*k0.^2 + eyx.*kx.*kz + eyz.*kz.^2)) + ...
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

det_out = log10(abs(det_out)+1);

end