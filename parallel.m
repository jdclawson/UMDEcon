
clear
clc

xi_vec = [0,.1,.2,.3,.4,.5];
sig_u = [.03,.04,.05,.06,.07,.08];

xilen=length(xi_vec);
siglen=length(sig_u);

noisyoutput=false;
tic

parfor ii=1:xilen
    for jj=1:siglen
        fprintf("Running the model with xi=%f and sig_e=%f",xi_vec(ii),sig_u(jj))

        model(xi_vec(ii),sig_u(jj),noisyoutput);


    end



end
toc
