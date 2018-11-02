%% �J�n
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% �f�[�^
data = datas_all;

%% �p�����[�^
nCohPoints = 32;

%% �R�q�[�����g�ϕ�
data_coh = coherentIntegration(data,nCohPoints);

%% �ۑ�
filename = strcat('data_coh_',num2str(nCohPoints),'pt');
save(filename,'data_coh','rdrdate','nCohPoints','-v7.3');

%% �I��
toc