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

num_points = 10;
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

num_k_per = 10;

k_per_min = -1000;
k_per_max = +1000;

k_per = linspace(k_per_min,k_per_max,num_k_per);

det_array = zeros(num_k_per,num_points);

for k=1:num_k_per
    
    n_per = k_per(k) * c / w;
    
    det_array(k,:) = A1 .* n_per^4 - B1 .* n_per^2 + C1;
    
end



% Try a root finder approach

for i=1:num_points
    
    
    B0 = B(i);
    n0 = n(:,i);
    
    %     n_per_init = [0,0];
    
    initial_k_pers = linspace(k_per_min,k_per_max,10);
    num_init_k_per = numel(initial_k_pers);
    
    for k = 1:num_init_k_per
        
        n_per_init_re = initial_k_pers(k) * c / w;
        
        n_per_init = n_per_init_re + 0*zi;
        
        [n_per_out(i,k),fval,exitflag,output] = fsolve(@det_fun,n_per_init);
        
    end
    
    for k = 1:num_init_k_per
        
        n_per_init_im = initial_k_pers(k) * c / w;
        
        n_per_init = 0 + zi*n_per_init_im;
        
        [n_per_out(i,k+num_init_k_per),fval,exitflag,output] = fsolve(@det_fun,n_per_init);
        
    end
    
    aa=A1(i);
    bb=B1(i);
    cc=C1(i);
    
    n_per_sq1(i) = sqrt((-bb + sqrt(bb^2-4*aa*cc))/(2*aa));
    n_per_sq2(i) = sqrt((-bb - sqrt(bb^2-4*aa*cc))/(2*aa));
    n_per_sq3(i) = -sqrt((-bb + sqrt(bb^2-4*aa*cc))/(2*aa));
    n_per_sq4(i) = -sqrt((-bb - sqrt(bb^2-4*aa*cc))/(2*aa));
    
    disp(['n_per_out : ', num2str(n_per_out(i))]);
    disp(['n_per_q1 : ', num2str(n_per_sq1(i))]);
    disp(['n_per_q2 : ', num2str(n_per_sq2(i))]);
    disp([' . ']);
    
end

% % Plot the zero contour of the determinant
% 
% figure()
% subplot(3,1,1)
% semilogy(x,n(1,:))
% subplot(3,1,2)
% plot(x,B(:))
% subplot(3,1,3)
% max_val = max(abs(det_array(:)));
% norm_det_array = det_array./max_val;
% cc = contour(x,k_per,norm_det_array,[-0.01,-0.001,0.0,0.001,+0.01]);
% 
% % spatial plot
% 
% figure()
% plot(x,real(n_per_sq1),'-b');
% hold on
% plot(x,real(n_per_sq2),'-b');
% plot(x,real(n_per_sq3),'-b');
% plot(x,real(n_per_sq4),'-b');
% plot(x,imag(n_per_sq1),'-r');
% plot(x,imag(n_per_sq2),'-r');
% plot(x,imag(n_per_sq3),'-r');
% plot(x,imag(n_per_sq4),'-r');
% for bb=1:numel(n_per_out(1,:))
%     p = plot(x,real(n_per_out(:,bb)),'o','Color','black','LineWidth',2);
%     p = plot(x,imag(n_per_out(:,bb)),'r','Color','black','LineWidth',2);
% end
% hold off
% ylim([k_per_min,k_per_max]*c/w);

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
symlog();
disp([' ']);

%     function det_re_im = det_fun(n_per_re_im)
    function det = det_fun(n_per)
        
%     n_per = complex(n_per_re_im(1),n_per_re_im(2))
          
    [eps,sigma,S,D,P,R,L] = epsilon_cold(f, amu, Z, B0, n0);
%     T_keV = [0.01,0.01];
%     k_per = n_per * w / c;
%     [eps_hot] = epsilon_hot(f, amu, Z, B0, n0, T_keV, k_per, k_par);
    
    
    % From pg 177 Brambilla
    
%     disp([num2str(n_per)]);
    A1a = S;
    B1a = R*L + P*S - n_par^2 * (P+S);
    C1a = P*(n_par^2-R)*(n_par^2-L);
    
    det = A1a .* n_per^4 - B1a .* n_per^2 + C1a;
    
eps = transpose(eps);
    exx = eps(1,1);
    exy = eps(1,2);
    exz = eps(1,3);
    
    eyx = eps(2,1);
    eyy = eps(2,2);
    eyz = eps(2,3);
    
    ezx = eps(3,1);
    ezy = eps(3,2);
    ezz = eps(3,3);
    
    kx = k_per;
    kz = k_par;
    k0 = w/c;
    
    det = ezy.*k0.^4.*(-(exz.*eyx.*k0.^2) + exx.*eyz.*k0.^2 + ...
        kz.*(eyx.*kx + eyz.*kz)) + ...
        (-(ezx.*k0.^2) + kx.*kz).* ...
        (exy.*eyz.*k0.^4 - exz.*k0.^2.*(eyy.*k0.^2 + kx.^2 + kz.^2) + ...
        kx.*kz.*(eyy.*k0.^2 + kx.^2 + kz.^2)) + ...
        (-(ezz.*k0.^2) - kx.^2).* ...
        (-(exy.*eyx.*k0.^4) + (exx.*k0.^2 + kz.^2).*(eyy.*k0.^2 + kx.^2 + kz.^2));
    
%     det = k0.^2.*(exx.*((-(eyz.*ezy) + eyy.*ezz).*k0.^4 - (eyy + ezz).*k0.^2.*kx.^2 + ...
%             kx.^4) + kx.*(-(eyy.*ezx.*k0.^2) + eyx.*ezy.*k0.^2 + ezx.*kx.^2).*kz + ...
%          ((eyz.*ezy - (exx + eyy).*ezz).*k0.^2 + (exx + ezz).*kx.^2).*kz.^2 + ...
%          ezx.*kx.*kz.^3 + ezz.*kz.^4 + ...
%          exy.*k0.^2.*(eyz.*ezx.*k0.^2 - eyx.*ezz.*k0.^2 + eyx.*kx.^2 + ...
%             eyz.*kx.*kz) + exz.*(eyx.*ezy.*k0.^4 - ...
%             eyy.*(ezx.*k0.^4 + k0.^2.*kx.*kz) + ...
%            (ezx.*k0.^2 + kx.*kz).*(kx.^2 + kz.^2)));
%  det/det0
    
%     det_re_im(1) = real(det);
%     det_re_im(2) = imag(det);
    
    end


end