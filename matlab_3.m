% 
% % Signal Model: x[n] = A cos(2πfo n + φ) + w[n]
% % for n = 0,1,2,...N-1
% % where A > 0 and 0 < f0 < 1/2
% % A, f0, and φ are unknown
% 
% %  Signal: x[n] = A cos(2πfo n + φ) + w[n]
% %  let a1 =  Acosφ
% %  and a2 = -Asinφ
% %  -> s[n] = a1*cos2πfon + a2*sin2πfon
% %  in matrix form s = Ha
% 
% figure(1);
% clf
% A = 2;
% phi = 0.5
% f0 = 1
% x = linspace(0,2,100);
% 
% % Signal form 1
% % x[n] = A cos(2πfo n + φ) + w[n] 
% signal_form_1 = A*cos(2*pi*f0.*x + phi);
% % Signal form 2
% %  let a1 =  Acosφ
% %  and a2 = -Asinφ
% %  -> s[n] = a1*cos2πfon + a2*sin2πfon
% 
% a1 = A*cos(phi);
% a2 = -A*sin(phi);
% 
% signal_form_2 = a1*cos(2*pi*f0.*x) + a2*sin(2*pi*f0.*x);
% 
% 
% plot(x, signal_form_1)
% hold on
% plot(x, signal_form_2)
% 
% % LSE for signal forms equals pretty much to zero
% sum((signal_form_1 - signal_form_2).^2)
% %  7.0251e-30
% 
% 
% %
% 
% N = 5;
% syms A f0 phi n;
% x = sym('x_%d', [1 N]);
% 
% symsum((x - A*cos(2*pi*f0*n+phi)).^2, n,0, N-1)
% 

% syms A f0 n phi w N x
% xn_sym = (x - A*cos(2*pi*f0*n + phi)).^2
% fd_xn_sym = gradient(xn_sym, [A,f0,phi]);
% 
% 


% sum(estimator(ns,0.25)')
% estimator = @(n,f0) [
%     pi/2 - 2*pi*f0*n
%           -2*pi*f0*n
%           -2*pi*f0*n
% ];


% 
% N  = 5;
% ns = 0:1:N-1;
% A  = 1;
% f0 = 0.25;
% phi = pi/3;
% sigma_squared = 0.001;
% 
% sn = @(A, n,phi,f0) A*cos(2*pi*f0*n + phi);
% w  = @(n, sigma_squared) randn(1, length(n))*sqrt(sigma_squared)
% xn = @(A, n,phi,f0,sigma_squared) sn(A,n,phi,f0) + w(n,sigma_squared);
% J  = @(A,f0,phi) sum(power(xn(A, ns,phi,f0,sigma_squared) - sn(A,ns,phi,f0),2));
% 





% Signal Model: x[n] = A cos(2πfo n + φ) + w[n]
% for n = 0,1,2,...N-1
% where A > 0 and 0 < f0 < 1/2
% A, f0, and φ are unknown

%  Signal: x[n] = A cos(2πfo n + φ) + w[n]
%  let a1 =  Acosφ
%  and a2 = -Asinφ
%  -> s[n] = a1*cos2πfon + a2*sin2πfon
%  in matrix form s = Ha

% x/cos(phi + 2*pi*f0*n)
% 
% MC = 100000;
% 
% simulation = zeros(1,MC);
% 
% for i = 1:1:MC
%     
%     simulation(i) = sum(power(xn(A, ns,phi,f0,sigma_squared) - sn(A,ns,phi,f0), 2));
%     sum(estimator(ns,0.25)')
% end



N             = 5;
A             = 1;
sigma_squared = 0.001;
f0            = 1/4;
phi           = pi / 3;
MC            = 100000;


f02pi = f0 * 2 * pi
ns    = 0:1:N-1;

c = cos(ns .* f02pi);
s = sin(ns .* f02pi);

a1 =  A*cos(phi);
a2 = -A*sin(phi);

a = [a1 a2]';

xn = @(n) A*cos(2*pi*f0*n + phi) + randn(1, length(n))*sqrt(sigma_squared);

H = [c' s'];

est = zeros(MC, 3);

weird_n = (A^2)/(2*sigma_squared);

CRLB_A = 2*sigma_squared/N;
CRLB_f0 = 12 / ((2*pi)^2*weird_n*N*(N^2 - 1));
CRLB_phi = 2*(2*N-1)/weird_n*N*(N+1);

CRLB = [CRLB_A, CRLB_f0, CRLB_phi];

for i = 1:MC
    
    x = xn(ns)';
    % alpha_hat = inv(H'*H)*H'*x;
    f0_hat    = x'*H * inv(H'*H)*H'*x;
    A_hat     = 2/N * abs(dot(x', exp(-j*2*pi*f0_hat*ns)));
    phi_hat   = atan(-dot(x', sin(ns * 2 * pi * f0_hat))/dot(x', cos(ns * 2 * pi * f0_hat)));
    
    est(i,:) = [A_hat, f0_hat, phi_hat];
    
end







