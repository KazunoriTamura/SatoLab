%% test mean-chisq

K = 2;
X = chi2rnd(K,[10000000,1]);
size(X)
mean_chisq(X,K,32)
mean_chisq_v2(X,K,32)