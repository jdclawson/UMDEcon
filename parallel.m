
clear
clc

xi_vec = [0,.1,.2,.3,.4,.5];
sig_u = [.03,.04,.05,.06,.07,.08];

xilen=length(xi_vec);
siglen=length(sig_u);

Gini = zeros(xilen*siglen);
var_c = zeros(xilen*siglen);
r = zeros(xilen*siglen);
cons = zeros(xilen*siglen);

noisyoutput=false;
tic

parfor ii=1:xilen
    for jj=1:siglen
        fprintf("Running the model with xi=%f and sig_e=%f",xi_vec(ii),sig_u(jj))

        [Gini(ii,jj),var_c(ii,jj),r(ii,jj),cons(ii,jj)]=model(xi_vec(ii),sig_u(jj),noisyoutput);


    end



end
toc
