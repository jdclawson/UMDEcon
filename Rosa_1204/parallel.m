
clear
clc

xi_vec = [0,.1,.2,.3,.4,.5];
sig_u = [.03,.04,.05,.06,.07,.08];
%sig_u = [.05];


xilen=length(xi_vec);
siglen=length(sig_u);
Gini_mat = zeros(xilen,siglen);
varc_mat = zeros(xilen,siglen);
r_mat = zeros(xilen,siglen);
constrained_mat = zeros(xilen,siglen);

noisyoutput=false;
tic;

parfor ii=1:xilen
    for jj=1:siglen
        disp(' ')
        fprintf("Running the model with xi=%f and sig_e=%f",xi_vec(ii),sig_u(jj))
        disp(' ')

        [Gini_mat(ii,jj),varc_mat(ii,jj),r_mat(ii,jj),constrained_mat(ii,jj)]=parmodel(xi_vec(ii),sig_u(jj),noisyoutput);


    end
end
csvwrite('gini.csv',Gini_mat);
csvwrite('varc.csv',varc_mat);
csvwrite('r.csv',r_mat);
csvwrite('constrained.csv',constrained_mat);
toc
