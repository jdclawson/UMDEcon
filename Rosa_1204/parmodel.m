%% This script solves 630 HW3 - Similar model with Aiyagari(1994)
% Jeff Clawson, Rosa (Heehyun) Lim
function [Gini,var_c,r,constrained] = parmodel(xi,sig_e,noisyoutput)

%% 1. Setting

% -- Parameters
%noisyoutput= true;
%xi = 0;
%sig_e = .05;
[sigma beta alpha delta M rho tol tol_r itermax rss kss wss ne na al au grid_e tran dist e2 e grid_a a2 a a_prime N Kss css uss] = parasetting( xi,sig_e );

%% Capital market clearing for getting r

% r iteration set-up
tol_r = 1e-2;
itermax_r = itermax;
found_r = 1;
% dev_r = 10;
exsaving=100;
devsaving=abs(exsaving);
iter_r = 0;

rmin = (1/beta - 1) * 0.9;
rmax = (1/beta - 1);
r = (rmin+rmax)/2;

mu0 = (1/(na*ne))*ones(na,ne);
dr0 = 0.0001;

while devsaving >= tol_r && iter_r < itermax_r
    
    exsavingold = exsaving;
    rold = r;
    iter_r = iter_r + 1;
    
    % Get K and w with r guessed
    K = N*(1/alpha * (r+delta))^(1/(alpha-1)); % MPK = r + delta
    w = (1-alpha)*(K/N)^alpha; % MPL = w
    
    % Consumption & Utility
    C = (1+r).*a + w.*e - a_prime; % for each given r and w
    fea = ones(na,ne,na);
    fea(C<0) = -inf;
    % C(C<0) = 0; % C>=0 constraint
    U = C.^(1-sigma)/(1-sigma); % for each given r and w
    C = C.*fea;
    U = U.*fea;
    
    %% Value Function Iteration
    if noisyoutput==true
        disp('Value Function Iteration..')
    end
    
    % Initial guess
    if iter_r == 1
        v = uss(:,:,1); % deterministic steady state utility for initial guess
        Ev = uss;
    else
        v = V; % better initial guess
        Ev = V; % better initial guess
    end
    
    % Value Function Iteration
    [V,index,iter_vfi]=vfi( tol, itermax, beta, U, Ev, v, na, ne, tran,noisyoutput);
    
    % Decision Rule
    for i = 1:na
        for j = 1:ne
            a_choice(i,j) = grid_a(index(i,j)); % a_(t+1)
            c_choice(i,j) = C(i,j,index(i,j));
        end
    end
    
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
    if iter_r == 1
        mu0 = [ones(na/10,ne)/(na/10*ne);zeros(na*9/10,ne)];
        %     mu0 = ones(na,ne)./(na*ne);
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
    exsaving = agg_A - agg_K;
    devsaving = abs(exsaving);
    
    if noisyoutput==true
        disp(' ')
        fprintf('exsaving is %5.6f with interest rate %5.10f \n', exsaving, r)
    end
    
    dr = dr0;
    
    if devsaving > 1
        dr = dr;
    elseif devsaving > 1e-1
        dr = dr*0.1;
    elseif devsaving > 1e-2
        dr = dr*0.01;
    elseif devsaving > 1e-3
        dr = dr*0.001;
    elseif devsaving > 1e-4
        dr = dr*0.0001;
    end
    
    if exsaving > 0
        if noisyoutput==true
            disp(' ')
            disp('Excess saving; LOWER r')
        end
        r = r-dr;
    else
        if noisyoutput==true
            disp(' ')
            disp('Excess capital demand; RAISE r')
        end
        r = r+dr;
    end
end
%
%   figure(1)
%     plot(grid_a, c_choice(:,1)',grid_a, c_choice(:,2)',grid_a, c_choice(:,3)',grid_a, c_choice(:,4)',grid_a, c_choice(:,5)')
%     figure(2)
%     plot(grid_a, a_choice(:,1)',grid_a, a_choice(:,2)',grid_a, a_choice(:,3)',grid_a, a_choice(:,4)',grid_a, a_choice(:,5)')

if noisyoutput==true
    disp(' ')
    if devsaving >= tol_r
        disp(' ')
        disp('We did not find r')
    else
        disp(' ')
        disp('We found r')
        r = rold;
        fprintf('Iteration : %5.0f times \n', iter_r)
    end
    disp(' ')
    disp('Model is solved')
    fprintf('Excess saving: %5.6f \n', exsaving)
    fprintf('Interest rate: %5.6f \n', r)
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
wfb = Wfb*ones(na,ne); % individual value for each (a,e)

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
constrained = sum(mu(1,:));

%% Result
if noisyoutput==true
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
fprintf('Fraction of Constrained : %5.3f\n', constrained)
disp('***********************************************')

else
disp('')
fprintf('Completed for xi: %5.3f and sig_e: %5.3f', xi, sig_e)
disp('')

end

%% Graph

%{
% Equilibrium
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
