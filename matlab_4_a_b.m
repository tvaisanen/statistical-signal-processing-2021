%% Declare variables

N = 5;
mu_A = 3;
sigma_A_sqr = 2;
sigma_sqr = 1.25;

n = 0:1:N-1;

%% Check that the theory understood correctly

A = normrnd(mu_A,sigma_A_sqr, 1, 1);
w = normrnd(0, sigma_sqr, 1, 1);

x = A + w;

H = [1];

% Check 10.30 ~ 10.31 for scalar
EAx = mu_A + sigma_A_sqr * H*inv( H*sigma_A_sqr*H'+ sigma_sqr*1)*(x - H*mu_A)
EAx_2 = mu_A + (sigma_A_sqr / (sigma_A_sqr + sigma_sqr / 1))*(x - H*mu_A)

x = A + w;

H = ones(5,1);

% Check 10.30 ~ 10.31 for N
A_hat = mu_A + sigma_A_sqr * H'*inv(H*sigma_A_sqr*H'+ sigma_sqr*eye(N))*(x - H*mu_A)
EAx_2 = mu_A + (sigma_A_sqr / (sigma_A_sqr + sigma_sqr / N))*(mean(x) - mu_A)
% These match

%% Exercise 4.a

%% Use signal model K*A , K = 1
K = 1;
x = @(K) K*normrnd(mu_A,sigma_A_sqr, N, 1) + normrnd(0, sigma_sqr, N, 1);

MC = 10000;
est = zeros(MC,1);
for i = 1:MC
    A_hat = mu_A + (sigma_A_sqr / (sigma_A_sqr + sigma_sqr / N))*(mean(x(K)) - mu_A);
    est(i) = A_hat;
end

figure(1)
clf
plot(1:MC, est)
mean(est)

%% Exercise 4.b
%% Use signal model K*A , K = 2
K = 2;
x = @(K) K*normrnd(mu_A,sigma_A_sqr, N, 1) + normrnd(0, sigma_sqr, N, 1);

MC = 100000;
est = zeros(MC,1);
for i = 1:MC
    A_hat = mu_A + (sigma_A_sqr / (sigma_A_sqr + sigma_sqr / N))*(mean(x(K)) - mu_A);
    est(i) = A_hat;
end

figure(2)
clf
plot(1:MC, est)
mean(est)




