%% Shooting algorithm to get a transition between different xi values
% Jeff Clawson, Rosa (Heehyun) Lim

clear
clc
% digits(50); % for shooting
% digits(64);

%% Print option
noisyoutput= true;

%% First steady state (xi = 0)
xi = 0; % borrowing limit coefficient

[ra, Va, ga, aa, ca, mua, Ka, Aa] = solve( xi, noisyoutput);

r_0 = ra;
V_0 = Va;
g_0 = ga;
a_0 = aa;
c_0 = ca;
mu_0 = mua;
K_0 = Ka;
A_0 = Aa;

%% Second steady state (xi = 0.5)
xi = 0.5;

[ra, Va, ga, aa, ca, mua, Ka, Aa] = solve( xi, noisyoutput);

r_1 = ra;
V_1 = Va;
g_1 = ga;
a_1 = aa;
c_1 = ca;
mu_1 = mua;
K_1 = Ka;
A_1 = Aa;

%% Shooting Algorithm
if noisyoutput==true
    disp('Shooting starts..')
end

% Parameters
[sigma beta alpha delta M rho sig_e tol tol_r itermax rss kss wss ne na al au grid_e tran dist e2 e grid_a a2 a a_prime N Kss css uss] = parasetting( xi );

% itermax_r = itermax;
% found_r = 1;
% exsaving=100;
% devsaving=abs(exsaving);
% iter_r = 0;
% mu0 = (1/(na*ne))*ones(na,ne);

% Shooting Set-up
dev = 100;
dev_tr = 100;
tol_tr = 10e-5;
% tol_r = tol;
iter_tr = 0;
T = 30; % try with small T and re-do with bigger T (use the best rpath from the first try)
dr0 = 0.0001;

% r sequence during transition
r_tr = vpa(r_0:(r_1-r_0)/(T-1):r_1); % initial guess
% r_tr = 0.95 * r_1 * ones(T,1); % initial guess
r_tr(T) = r_1;

V_tr = repmat(V_1,[1,1,T]); % V on transition path
g_tr = repmat(g_1,[1,1,T]); % g on transition path
v = V_1;
K_tr = zeros(T,1);
K_tr(T) = K_1;

while dev > tol_tr
    
    if noisyoutput==true
        disp(' ')
        fprintf('Shooting iteration : %5.0f\n', iter_tr)
        fprintf('Deviation (dev_tr) : %5.6f', dev)
    end
    
    
    r_tr_old = vpa(r_tr);
    r_tr_im = zeros(T,1);
    r_tr_im(T) = r_1;
    
    if noisyoutput==true
        disp(' ')
        disp('Solve V and g by backward induction')
    end
    % Solve V and g by backward induction
    for tt = 1:T-1
        % Given r, back-out the economy state
        r = r_tr(T-tt);
        K = vpa(N*(1/alpha * (r+delta))^(1/(alpha-1))); % MPK = r + delta
        K_tr(T-tt) = K;
        w = vpa((1-alpha)*(K/N)^alpha); % MPL = w
        C = vpa((1+r).*a + w.*e - a_prime); % for each given r and w
        C(C<0) = 0; % C>=0 constraint
        U = vpa(C.^(1-sigma)./(1-sigma)); % for each given r and w
        
        % Solve maximization and get the value for this period
        prob_v = vpa(sum(repmat(v,ne,1).*repelem(tran,na,1),2));
        prob_v = reshape(prob_v,[na,ne]);
        Ev = zeros(1,ne,na);
        Ev(1,:,:) = prob_v';
        Ev = repmat(Ev,[na,1,1]);
        RHS = (1-beta)*U + beta*Ev;
        [V, g] = max(RHS, [], 3);
        
        V_tr(:,:,T-tt) = V;
        g_tr(:,:,T-tt) = g;
        
    end
    
    if noisyoutput==true
        disp(' ')
        disp('Backward induction done.')
        disp('Get the mu distribution')
    end
    
    
    mu_tr = zeros(na,ne,T);
    mu_tr(:,:,T) = mu_1;
    mu_prime = zeros(na,ne);
    mu = mu_0;
    
    % Get the mu distribution
    for tt=1:T-1
        
        for i = 1:na
            for j = 1:ne
                mu_prime(g(i,j),:) = mu_prime(g(i,j),:) + mu(i,j)*tran(j,:);
            end
        end
        mu = mu_prime;
        mu_tr(:,:,tt) = mu_prime;
        
    end
    dev_tr = zeros(T-1,1);
    
    if noisyoutput==true
        disp(' ')
        disp('Get the implied r that clears K = A market for each transition period')
    end
    
    % Get the implied r that clears K = A market for each period tt
    for tt = 1:T-1
        %         tt = 1;
        r = r_tr(tt);
        rmin = r_0;
        rmax = r_1;
        devsaving = 100;
        iter_r = 0;
        g = g_tr(:,:,tt);
        
        while devsaving > tol_r
            
            % Check the capital = asset market
            mu = mu_tr(:,:,tt);
            agg_K = vpa(K_tr(tt));
            agg_A = vpa(N*sum(sum(mu.*a2)));
            A_tr(tt) = vpa(agg_A);
            
            ex_saving = vpa(agg_A - agg_K);
            devsaving = abs(ex_saving);
            dev_tr(tt) = devsaving;
            
            
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
        r_tr_im(tt) = vpa(r); % implied r at period tt
    end
    A_tr = A_tr';
    
    % Update r path with both old guess and implied (new) r path
    dev_tr = max(abs(r_tr_old - r_tr_im));
    
    dev = max(dev_tr);
    r_tr = vpa(0.8*r_tr_old + 0.2*r_tr_im);
    iter_tr = iter_tr + 1;
    
    if iter_tr > itermax
        break;
    end
    
    %     fprintf('Shooting iteration : %5.0f\n', iter_tr)
    %     fprintf('Deviation (dev_tr) : %5.6f', dev)
    disp(' ')
end

if noisyoutput == 1
    disp(' ')
    disp('Transition solved')
end
r_path = r_tr_old;

