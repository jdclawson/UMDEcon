function [y,p,pi] = tauchen(rho,sig_e,n)

% tauchen.m
% Finds the transition probabilities for a Markov chain that
% approximates an AR(1):
%
%  y[t] = rho*y[t-1]+e[t],  e[t] ~ N(0,sig_e^2).     (1)
%
% where y is the points of approximation, p the transition kernel
% and pi the invariant distribution
%
% References: 
% Tauchen, G. "Finite State Markov-chain approximations
% to Univariate and Vector Autoregressions," Economics 
% Letters 20, pp. 177-181, 1986.
%
% Jesus Fernandez-Villaverde
% Minneapolis, 5-29-2001

	%----------------------------------------------------------------
	% 1. Picking the points in y
	%----------------------------------------------------------------

	y(n)  = 3*sig_e/sqrt(1-rho^2); % Change 3 to increase dispersion
	y(1)  = -y(n);
	w     = 2*y(n)/(n-1);
	y     = y(1):w:y(n);
   
  %----------------------------------------------------------------
	% 2. Computing the transitions
	%----------------------------------------------------------------
   
	p(1:n,1) = normcdf(((y(1)-rho*y(1:n)+(w/2))/sig_e)');
	p(1:n,n) = 1-normcdf(((y(n)-rho*y(1:n)-(w/2))/sig_e)');

	for j = 1:n
      
  		p(j,2:n-1) = normcdf((y(2:n-1)-rho*y(j)+(w/2))/sig_e)- normcdf((y(2:n-1)-rho*y(j)-(w/2))/sig_e);
      
 	end
  
  %----------------------------------------------------------------
	% 3. Computing the invariant distribution
	%----------------------------------------------------------------    
   
  pi    = ones(n,1)/n;
  	piprov = pi+1;
        
  while (norm(pi-piprov)>0.00001)
  		piprov = pi;
	  pi = (pi'*p)';
	end