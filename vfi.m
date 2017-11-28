function [V,index, iter_vfi ] = vfi( dev_vfi, tol_vfi, beta, U, Ev, v, iter_vfi, itermax_vfi, converge, na, ne, tran,noise)
% VFI to get value function and decision rule (of asset choice)

tic;
while dev_vfi > tol_vfi
    
    RHS = (1-beta)*U + beta*Ev; % normalization
    [Tv, index] = max(RHS, [], 3);
    
    iter_vfi = iter_vfi+1;
    if iter_vfi > itermax_vfi
        converge = 0;
        disp('VF did NOT converge!')
        break;
    end
  
    dev_vfi = max(abs(reshape(Tv,[na*ne,1]) - reshape(v,[ne*na,1])));
    v = Tv; 
                    
    prob_v = sum(repmat(v,ne,1).*repelem(tran,na,1),2);
    prob_v = reshape(prob_v,[na,ne]);
    Ev = zeros(1,ne,na);
    Ev(1,:,:) = prob_v';
    Ev = repmat(Ev,[na,1,1]); % expected value

end
toc;

if converge == 0
    if noise==true
    disp('VF did NOT converge!')
    end
else
    if noise==true
    disp(' ')
    disp('VF converged.')
    fprintf('Elapsed Time : %5.3f seconds \n', toc)
    fprintf('Iteration : %5.0f times \n', iter_vfi)
    end
    V = Tv; % Value conditional for state l and s (= W(l,s))
end

% Value function 
V = V./(1-beta); % un-normalize

end

