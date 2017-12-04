function [ ra, Va, ga, aa, ca, mua, Ka, Aa] = solve( xi, noisyoutput)

% This code solves the model for HW3 (ECON 630, Fall 2017)
% Jeff Clawson, Rosa (Heehyun) Lim

% Parameters
[sigma beta alpha delta M rho sig_e tol tol_r itermax rss kss wss ne na al au grid_e tran dist e2 e grid_a a2 a a_prime N Kss css uss] = parasetting( xi );

itermax_r = itermax;
found_r = 1;
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
%     mu = mu_final;
    
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

ra = r;
Va = V;
ga = index;
aa = a_choice;
ca = c_choice;
mua = mu;
Ka = agg_K;
Aa = agg_A;

end

