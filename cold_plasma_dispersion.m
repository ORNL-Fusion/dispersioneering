function cold_plasma_dispersion()

phys = constants();

c = phys.('c');
me_amu = phys.('me_amu');
zi = complex(0,1);

% Map the space over a variable - here B

f = 7.5e6;
w = 2*pi*f;
k_par = 20;
n_par = k_par * c / w;

% Species

amu = [me_amu, 2];
Z   = [-1,1];
T_eV = [10,10];

num_points = 50;
x = linspace(0,0.1,num_points);

% B field

B = x.*0 + 1.2;

% Density

den = 10.^linspace(18,20,num_points);

n = [den; den];

% Produce 2D space over k_per and B

num_points = numel(B);

A1 = zeros(num_points,1);
B1 = zeros(num_points,1);
C1 = zeros(num_points,1);

for i=1:num_points

    [eps,sigma,S,D,P,R,L] = epsilon_cold(f, amu, Z, B(i), n(:,i));
    
    % From pg 177 Brambilla
    
    A1(i) = S;
    B1(i) = R*L + P*S - n_par^2 * (P+S);
    C1(i) = P*(n_par^2-R)*(n_par^2-L);

end

% k_per range

num_init_k_per = 10;

k_per_min = -1000;
k_per_max = +1000;

k_per = linspace(k_per_min,k_per_max,num_init_k_per);

% Try a root finder approach

n_per_out = complex(zeros(num_points,num_init_k_per^2),zeros(num_points,num_init_k_per^2));

initial_k_pers = linspace(k_per_min,k_per_max,num_init_k_per);


options = optimoptions('fsolve',...
    'Display','none',...
    'SpecifyObjectiveGradient',false,...
    'CheckGradients',false,...
    'UseParallel',false,...
    'StepTolerance',1e-6,...
    'FunctionTolerance',1e-6,...
    'OptimalityTolerance',1e-6);

for i=1:num_points
     
    B0 = B(i);
    n0 = n(:,i);
    
    cnt = 1;
    for k = 1:num_init_k_per
        for j = 1:num_init_k_per
            
            n_per_init_re = initial_k_pers(k) * c / w;
            n_per_init_im = initial_k_pers(j) * c / w;
            
            n_per_init = n_per_init_re + n_per_init_im*zi;
            
            [n_per_out(i,cnt),fval,exitflag,output] = ...
                fsolve(@det_fun,n_per_init,options);
            
            cnt = cnt + 1;
            
        end
    end
    
    aa=A1(i);
    bb=B1(i);
    cc=C1(i);
    
    n_per_sq1(i) = sqrt((-bb + sqrt(bb^2-4*aa*cc))/(2*aa));
    n_per_sq2(i) = sqrt((-bb - sqrt(bb^2-4*aa*cc))/(2*aa));
    n_per_sq3(i) = -sqrt((-bb + sqrt(bb^2-4*aa*cc))/(2*aa));
    n_per_sq4(i) = -sqrt((-bb - sqrt(bb^2-4*aa*cc))/(2*aa));
    
    disp(['point ', num2str(i), ' of ', num2str(num_points)]);
    disp(['n_per_out : ', num2str(n_per_out(i))]);
    disp(['n_per_q1 : ', num2str(n_per_sq1(i))]);
    disp(['n_per_q2 : ', num2str(n_per_sq2(i))]);
    disp(['  ']);
    
end

% density plot

figure()
semilogx(den,real(n_per_sq1)*w/c,'-b');
hold on
plot(den,real(n_per_sq2)*w/c,'-b');
plot(den,real(n_per_sq3)*w/c,'-b');
plot(den,real(n_per_sq4)*w/c,'-b');
plot(den,imag(n_per_sq1)*w/c,'-r');
plot(den,imag(n_per_sq2)*w/c,'-r');
plot(den,imag(n_per_sq3)*w/c,'-r');
plot(den,imag(n_per_sq4)*w/c,'-r');
for bb=1:numel(n_per_out(1,:))
    p = plot(den,real(n_per_out(:,bb)*w/c),'o','Color','black','LineWidth',2);
    p = plot(den,imag(n_per_out(:,bb)*w/c),'o','Color','red',  'LineWidth',2);
end
hold off
ylim([k_per_min,k_per_max]);


figure()
semilogx(den,real(n_per_sq1)*w/c,'-b');
hold on
plot(den,real(n_per_sq2)*w/c,'-b');
plot(den,real(n_per_sq3)*w/c,'-b');
plot(den,real(n_per_sq4)*w/c,'-b');
plot(den,imag(n_per_sq1)*w/c,'-r');
plot(den,imag(n_per_sq2)*w/c,'-r');
plot(den,imag(n_per_sq3)*w/c,'-r');
plot(den,imag(n_per_sq4)*w/c,'-r');
for bb=1:numel(n_per_out(1,:))
    p = plot(den,real(n_per_out(:,bb)*w/c),'o','Color','black','LineWidth',2);
    p = plot(den,imag(n_per_out(:,bb)*w/c),'o','Color','red',  'LineWidth',2);
end
hold off
ylim([k_per_min,k_per_max]);
symlog();
disp([' ']);

    function [det,det_gradient] = det_fun(n_per)       
        
    % Cold dielectric
        
    % [eps,sigma,S,D,P,R,L] = epsilon_cold(f, amu, Z, B0, n0);

    % Quadratic form for determinant (works only for cold)
    % from pg 177 of Brambilla
    %     
    % A1a = S;
    % B1a = R*L + P*S - n_par^2 * (P+S);
    % C1a = P*(n_par^2-R)*(n_par^2-L);
    % 
    % det = A1a .* n_per^4 - B1a .* n_per^2 + C1a; 
    
    % Hot dielectric
    
    k_per = n_per .* w / c;
    [eps] = epsilon_hot(f, amu, Z, B0, n0, T_eV, k_per, k_par);
    
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
        
        det_gradient = -((k0.^2.*w.*(c.^3.*kz.*(exy.*eyz.*k0.^2 - eyy.*ezx.*k0.^2 +...
                eyx.*ezy.*k0.^2 + ezx.*kz.^2 + exz.*(-(eyy.*k0.^2) + kz.^2)) + ...
             2.*c.^2.*(exy.*eyx.*k0.^2 + exz.*ezx.*k0.^2 + ezz.*kz.^2 +...
                exx.*(-(eyy.*k0.^2) - ezz.*k0.^2 + kz.^2)).*n_per.*w + ...
             3.*c.*(exz + ezx).*kz.*n_per.^2.*w.^2 + 4.*exx.*n_per.^3.*w.^3))/c.^4);
    end
    
    end

end