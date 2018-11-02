%% �J�n
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% �p�����[�^
Tx = 1;
Rx = 1;

%% �����W�v���t�@�C���i�P�`���l���j�̕\��
figure(1)
fig = figure(1);
fig.OuterPosition = [300,300,2103,1300];
T = ['Range Profile','�i Tx',num2str(Tx),'-Rx',num2str(Rx),' �j'];
p = plotRangeProfile(data_coh,'UseTx',Tx,'UseRx',Rx,'Title',T);
set(p,'LineWidth',4);
ax = gca;
ax.XAxis.FontSize = 40;
ax.YAxis.FontSize = 40;
ax.Legend.FontSize = 50;
ax.TitleFontSizeMultiplier = 4;

%% �ۑ�
figName = ['RangeProfile_','Tx',num2str(Tx),'_Rx',num2str(Rx)];
print(['Figures/', figName, '.png'],'-dpng','-r0');

%% �I��
toc