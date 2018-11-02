%% 開始
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% データ
data = datas_all;

%% パラメータ
nCohPoints = 32;

%% コヒーレント積分
data_coh = coherentIntegration(data,nCohPoints);

%% 保存
filename = strcat('data_coh_',num2str(nCohPoints),'pt');
save(filename,'data_coh','rdrdate','nCohPoints','-v7.3');

%% 終了
toc