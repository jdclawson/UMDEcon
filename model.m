%% This script solves 630 HW3 - Similar model with Aiyagari(1994)
function model( xi,sig_e,noisyoutput )


%% 1. Setting

global sigma beta xi alpha delta M rho sig_e

% -- Parameters
%noisyoutput= true;
% Household
sigma = 2; % % CRRA coefficient
beta = 0.96; % utility discount factor
%xi = 0; % borrowing limit coefficient
% Production
alpha = 0.36; % capital share
delta = 0.06; % capital depreciation rate
% Shock
M = 5; % number of possible efficiency units
rho = 0.95; % serial correlation in efficiency shock
%sig_e = 0.05; % st.d of AR(1) process for efficiency

% -- Steady State
rss = 1/beta - 1;
kss = (1/alpha*(rss+delta))^(1/(alpha-1)); % per-capita capital
wss = (1-alpha)*kss^alpha;

% -- Grid
% No. of grids
ne = M; % stochastic efficiency
na = 500; % asset
% al = 0; % minimum point for e grid (amin = ass * al)
au = 8; % maximumm point for e grid (amax = ass * au)
% Efficiency
[log_e,tran,dist] = tauchen(rho,sig_e,ne);
grid_e = exp(log_e);
e2 = repmat(grid_e, [na 1]); % 2D efficiency values
e = repmat(grid_e, [na 1 na]); % na(asset state) by ne(efficiency) by na(asset choice)
% Asset
% amin = kss * al; % ass = kss
abar = -grid_e(1) * wss / rss; % natural borrowing limit
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

% -- Function Set-up 
tol = 10e-9;
itermax = 10000;


%% Capital market clearing for getting r 

% r iteration set-up
tol_r = 10e-3; % takes too long with 10e-9.. 
itermax_r = itermax;
iter_r = 0;
found_r = 1;
dev_r = 10;
rmin = (1/beta - 1) * 0.9;
rmax = (1/beta - 1); 
    % (1/beta - 1) * 0.9: excess saving -6.5133
    % (1/beta - 1) * 0.94: excess saving -1.7532
    % (1/beta - 1) * 0.95: excess saving 1.1041

% Initial guess for r
r0 = (rmin+rmax)/2; % initial guess 
r = r0;

tic;
while dev_r > tol_r    

r = (rmin+rmax)/2;
    
% Get K and w with r guessed
K = N*(1/alpha * (r+delta))^(1/(alpha-1)); % MPK = r + delta
w = (1-alpha)*(K/N)^alpha; % MPL = w

% Consumption & Utility
C = (1+r).*a + w.*e - a_prime; % for each given r and w
C(C<0) = 0; % C>=0 constraint
U = C.^(1-sigma)./(1-sigma); % for each given r and w


%% Value Function Iteration
if noisyoutput==true
disp('Value Function Iteration..')
end

% VFI set-up
tol_vfi = tol;
itermax_vfi = itermax;
iter_vfi = 0;
converge = 1;
dev_vfi = 100;

% Initial guess

if iter_r == 0
    v = uss(:,:,1); % deterministic steady state utility for initial guess
    Ev = uss;
else
    v = V; % better initial guess
    Ev = V; % better initial guess
end

% Value Function Iteration
[V,index,iter_vfi]=vfi(dev_vfi, tol_vfi, beta, U, Ev, v, iter_vfi, itermax_vfi, converge, na, ne, tran,noisyoutput);

% Decision Rule
for i = 1:na
    for j = 1:ne
        a_choice(i,j) = grid_a(index(i,j)); % a_(t+1)

    end
end
c_choice = (1+r)*a2 + w*e2 - a_choice;% c_t

%% Invariant Distibution
if noisyoutput==true
disp(' ')
disp('Finding mu..')
end

% mu function setup
tol_mu = tol;
itermax_mu = itermax;

% Initial guess
% mu0 = ones(na,ne)/(na*ne); % uniform distribution for an initial guess
if iter_r == 0
    mu0 = [ones(na/10,ne)/(na/10*ne);zeros(na*9/10,ne)];
else
    mu0 = mu; % updating the initial guess 
end

% Invariant Distribution
[mu] = musolve(mu0,index,ne,na,tran,tol_mu,itermax_mu);

%% Check capital market clearing 
if noisyoutput==true
disp(' ')
disp('Checking capital market clearing condition..')
end
agg_K = K;
agg_A = N*sum(sum(mu.*a2));

ex_saving = agg_A - agg_K; 
dev_r = abs(ex_saving);
if noisyoutput==true
    disp(' ')
    disp(' ')
    fprintf('ex_saving : %5.7f \n', ex_saving)
    fprintf('r : %5.7f \n', r)
end
    
if ex_saving > 0
    if noisyoutput==true
    disp('We have excess saving; need to lower r')
    end
    rmax = r;
else
    if noisyoutput==true
    disp('We have excess capital demand; need to raise r')
    end
    rmin = r;
end

iter_r = iter_r + 1;
    if iter_r > itermax_r
        found_r = 0;
        break;
    end


    
end

toc;

if noisyoutput==true
disp(' ')
if found_r == 0
    disp('We did not find r')
else
    disp('We found r')
    fprintf('Iteration : %5.0f times \n', iter_r)
end

    disp('Model is solved')
end

%% Euler Error 

[ ee ] = euler_error( c_choice, sigma, ne, na, tran, beta, r );

maxee = max(max(mu.*ee));
minee = min(min(mu.*ee));
meanee = mean(mean(mu.*ee));



%% Welfare Change
% Complete market 
Ncm = 1; % no idiosyncratic risk
Kcm = kss*Ncm; % kss comes from deterministic steady state calculation
Ccm = Kcm^alpha*Ncm^(1-alpha) - delta*Kcm; % aggregate consumption
ufb = Ccm^(1-sigma)/(1-sigma); % aggergate utility level
Wfb = 1/(1-beta)*ufb; % aggregate value level 
wfb = Wfb/(na*ne)*ones(na,ne); % individual value for each (a,e)

lam = (wfb./V).^(1/(1-sigma)) - 1;
Lam = sum(sum(lam.*mu));

% %Maximized steady state value function
% uss_in = max(uss,[],3);
% %Complete asset market will lead to steady state value function
% lam = (uss_in./V).^(1/(1-sigma))-1;
% 
% Lam = sum(sum(lam.*mu));

%% Gini Coefficient
% Distribution of H
H_a = sum(mu,2);
S_a = zeros(1,na+1);

for ii=1:na
    S_a(ii+1)=sum(grid_a(1:ii).*H_a(1:ii)');
end

lag_s = S_a(1:na)+S_a(2:na+1);
num=sum(H_a.*lag_s');

Gini = 1 - num/S_a(na+1);

% Variance of comsumption across agents
var_c = var(reshape(c_choice,na*ne,1));

%% Result 
disp(' ')
disp('***********************************************')
disp('*************** Result Summary ****************')
disp(' ')
fprintf('Time to solve the model : %5.3f seconds', toc)
disp(' ')
fprintf('Interest Rate           : %5.3f\n', r)
fprintf('Capital Stock           : %5.3f\n', agg_K)
disp(' ')
fprintf('Mean Euler (log10)      : %5.6f\n', meanee)
fprintf('Max Euler (log10)       : %5.6f\n', maxee)
fprintf('Min Euler (log10)       : %5.6f\n', minee)
disp(' ')
fprintf('Aggregate Value         : %5.3f\n', sum(sum(mu.*V)))
fprintf('Aggregate Value (CM)    : %5.3f\n', Wfb)
fprintf('Agg Welfare Gain (CM)   : %5.3f\n', Lam)
disp(' ')
fprintf('Gini                    : %5.3f\n', Gini)
fprintf('Variance of Consumption : %5.3f\n', var_c)
disp('***********************************************')
    
%% Graph 

% Equilibrium 
%{
figure(1)
subplot(2,2,1)
mesh(grid_a, grid_e, V')
xlabel('a')
ylabel('e')
axis tight
title('Value Function')

subplot(2,2,2)
mesh(grid_a, grid_e, a_choice')
xlabel('a')
ylabel('e')
axis tight
title('Asset Decision')

subplot(2,2,3)
mesh(grid_a, grid_e, c_choice')
xlabel('a')
ylabel('e')
axis tight
title('Consumption')

subplot(2,2,4)
mesh(grid_a, grid_e, mu')
xlabel('a')
ylabel('e')
axis tight
title('HH Distribution')

% Euler Error
figure(2)
mesh(grid_a, grid_e, ee')
xlabel('a')
ylabel('e')
axis tight
title('Euler Residual')

% Deviation from Complete Market
figure(3)
subplot(2,2,1)
mesh(grid_a, grid_e, V')
xlabel('a')
ylabel('e')
axis tight
title('Value Function')

subplot(2,2,2)
mesh(grid_a, grid_e, wfb')
xlabel('a')
ylabel('e')
axis tight
title('Complete Market Value')

subplot(2,2,3)
mesh(grid_a, grid_e, lam')
xlabel('a')
ylabel('e')
axis tight
title('Lambda Distribution')
%}
end

