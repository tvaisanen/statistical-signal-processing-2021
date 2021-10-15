
M             = 2;
A             = 1;
sigma_squared = 0.001;
f0            = 1/4;
phi           = pi / 3;
iterations    = 10;

f02pi = f0 * 2 * pi
ns    = -M:1:M;

h1 = ns*2*pi .* sin(2*pi*f0*ns + phi);
h2 = sin(2*pi*f0*ns + phi);

xn = @(n) A*cos(2*pi*f0*n + phi) + randn(1, length(n))*sqrt(sigma_squared);

f0_k1 = @(f0_k, phi_k) f0_k - (3 / 4*pi*M^3) * sum(ns .* xn(ns) .* sin(2*pi*f0_k*ns + 2*phi_k))
phi_k1 = @(f0_k, phi_k) phi_k - (1/M) * sum(xn(ns) .* sin(2*pi*f0_k*ns + 2*phi_k));

simulated = zeros(iterations,2);


f0_k  = f0;
phi_k = phi;

for i = 1:iterations
    f0_k  = f0_k1(f0_k, phi_k);
    phi_k = phi_k1(f0_k, phi_k);
    
    simulated(i,:) = [f0_k, phi_k];
end

figure(1)
clf
plot(1:10,simulated(:,1))
hold on
plot(1:10,simulated(:,2))

