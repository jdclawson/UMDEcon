function [ sigma beta alpha delta M rho tol tol_r itermax rss kss wss ne na al au grid_e tran dist e2 e grid_a a2 a a_prime N Kss css uss] = parasetting( xi,sig_e )

% To set parameters for the model
% Jeff Clawson, Rosa (Heehyun) Lim

% Household
sigma = 2; % % CRRA coefficient
beta = 0.96; % utility discount factor
% Production
alpha = 0.36; % capital share
delta = 0.06; % capital depreciation rate
% Shock
M = 5; % number of possible efficiency units
rho = 0.95; % serial correlation in efficiency shock
tol=1e-5;     % Level of tolerance in VFI
tol_r = 1e-2;
itermax=6000;  % Maximium number of iterations in VFI

% -- Steady State
rss = 1/beta - 1;
kss = (1/alpha*(rss+delta))^(1/(alpha-1)); % per-capita capital
wss = (1-alpha)*kss^alpha;

% -- Grid
% No. of grids
ne = M; % stochastic efficiency
na = 600; % asset
al = -2; % minimum point for e grid (amin = ass * al)
au = 8; % maximumm point for e grid (amax = ass * au)

% Efficiency
[log_e,tran,dist] = tauchen(rho,sig_e,ne);
grid_e = exp(log_e);
e2 = repmat(grid_e, [na 1]); % 2D efficiency values
e = repmat(grid_e, [na 1 na]); % na(asset state) by ne(efficiency) by na(asset choice)

% Asset
abar = -grid_e(1) * wss / rss; % natural borrowing limit
% amin = kss * al; % ass = kss
amin = xi * abar;
amax = kss * au;

grid_a = amin:(amax-amin)/(na-1):amax;
a2 = repmat(grid_a', [1 ne]); % 2D asset values
a = repmat(grid_a', [1 ne na]); % na(asset state) by ne(efficiency) by na(asset choice)
a_prime = zeros(1,1,na);
a_prime(1,1,:) = grid_a';
a_prime = repmat(a_prime, [na ne 1]);

%-- Steady state (continued)
N = grid_e * dist; % exogenously given
Kss = kss * N;
css = rss*a + wss*e;
uss = css.^(1-sigma)./(1-sigma); % steady state utility for initial guess

end

