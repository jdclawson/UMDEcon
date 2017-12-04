function [ ee ] = euler( c_choice, sigma, ne, na, tran, beta, r )
% This function calculates different measures of euler error 

% marginal utility 
marginu = c_choice.^(-sigma); 

% expected marginal utility
ex_marginu = sum(repmat(marginu,ne,1).*repelem(tran,na,1),2);
ex_marginu = reshape(ex_marginu,[na,ne]);

% euler error
ee = log10(abs(1 - beta*(1+r)* ex_marginu./marginu )) ; 

end

