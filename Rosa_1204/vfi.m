function [V,index, iter_vfi ] = vfi( tol, itermax, beta, U, Ev, v, na, ne, tran,noise)
% VFI to get value function and decision rule (of asset choice)

    tol_vfi = tol;
    itermax_vfi = itermax;
    iter_vfi = 0;
    converge = 1;
    dev_vfi = 100;

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

if converge == 0
    if noise==true
    disp('VF did NOT converge!')
    toc;
    end
else
    if noise==true
    disp(' ')
    disp('VF converged.')
    fprintf('Elapsed Time : %5.3f seconds \n', toc)
    fprintf('Iteration : %5.0f times \n', iter_vfi) 
    toc;
    end
    V = Tv; % Value conditional for state l and s (= W(l,s))
end


% Value function 
V = V./(1-beta); % un-normalize




% Using HW2 Q1 for reference/checking 
% 
% function [w,index, iter ] = vfi( dev, tol, beta, profit, Ew, w, iter, itermax, converge, nl, ns, tran,noise)
%     % VFI and labor decision
%     while dev > tol
% 
%         RHS = profit + beta * Ew;
%         [Tw, index] = max(RHS, [], 3);
%         
%         iter = iter+1;
%         if iter > itermax
%             converge = 0;
%             disp('VF did NOT converge!')
%             break;
%         end
%         dev = max(abs(reshape(Tw,[nl*ns,1]) - reshape(w,[nl*ns,1])));
%         w = Tw; % nl by ns
%         
%         prob_w = sum(repmat(w,ns,1).*repelem(tran,nl,1),2);
%         prob_w = reshape(prob_w,[nl,ns]);
%         Ew = zeros(1,ns,nl);
%         Ew(1,:,:) = prob_w';
%         Ew = repmat(Ew,[nl,1,1]); % expected value 
%         
%     end
%     
%     toc;
% %     
end

