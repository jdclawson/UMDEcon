function [ mu_final ] = musolve(mu_guess,policy,e_gridpts,a_gridpts,tran_e,tol,itermax)
%{
    Description:
        - This code, given an initial guess for the distribution of assets and
          the policy function calculates the fixed point distribution of assets
          across asset holdings and productivity values
    Inputs:
        - mu_guess       = Array[a_gridpts,e_gridpts], the initial guess of the distribution
                           typically a uniform distribution
        - policy         = Array[a_gridpts,e_gridpts], the policy function INDICIES.
        - params         = Array[7], contains the necessary parameters for setup of
                           the transition kernel.
    
    Other Functions Called:
        - tauchen        = Given rho, sig_e, the number of gridpoints and the dispersion
                           it calculates the grid, transition kernel and invariant
                           distribution for an AR(1) process. See that code for more details


    Outputs:
        - mu_final       = Array[a_gridpts,e_gridpts], fixed point distribution
%}

%Extracts the parameters from the input vector
%     param_list = num2cell(params);
%     [e_gridpts,a_gridpts,tran_e,tol,itermax] = deal(param_list{:});

%     mu_params = [ne,na,tran,tol_mu,itermax_mu];


%     %Uses tauchen to set up the needed transition kernel, the other two are unused
%     e_grid,tran_k,invdist = tauchen(rho,sig_u,e_gridpts,disp);

%Loads the guess to become mu_0
mu_0=mu_guess;

%Initializes distance to be 1
dist = 1;

%Initializes the iteration counter
iternum = 0;

progress = 1;
mu_final=zeros(a_gridpts,e_gridpts);

while dist>tol
    
    %Initializes the next period distribution
    mu_prime = zeros(a_gridpts,e_gridpts);
    
    for ii=1:a_gridpts
        for jj=1:e_gridpts
            
            %For each policy function choice, it takes the distribution from the mu_0
            %multiplies it by transition probability and places into mu_1 the new
            %distribution
            mu_prime(policy(ii,jj),:)=mu_prime(policy(ii,jj),:)+mu_0(ii,jj)*tran_e(jj,:);
                                
        end
        
    end
    
    dist = norm(mu_prime-mu_0);
    
    %Sets mu_0 to be mu_prime if the distance is too large and it tries again. If not,
    %then it saves mu_final to be the last mu_prime
    if dist>tol
        mu_0=mu_prime;
    else
        mu_final=mu_prime;
    end
    
    %Updates the iteration number
    iternum = iternum+1;
    if iternum > itermax
        progress = 0;
        break;
    end
end
%If the progress lever is true, then show the iteration number and distance
% if progress==true
%     print "Mu Iteration Number: ", iternum, "Distance: ",dist
% end

% if progress == 0
%         disp('no invariant Mu found!')
% else
%     disp(' ')
%     disp('Invariant mu is found')
%     fprintf('Iteration : %5.0f times \n', iternum)
% end



end

