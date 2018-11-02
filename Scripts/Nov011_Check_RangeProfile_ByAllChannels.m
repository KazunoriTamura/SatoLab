%% �J�n
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% �����W�v���t�@�C���i�P�`���l���j�̕\��
figure(2)
fig = gcf;
fig.OuterPosition = [100,100,2103,1300];
T = ['Range Profile','�i16�`���l�����ρj'];
p = plotRangeProfile(data_coh,'UseCh',1:16,'Title',T);
set(p,'LineWidth',4);
ax = gca;
ax.XAxis.FontSize = 40;
ax.YAxis.FontSize = 40;
ax.Legend.FontSize = 50;
ax.TitleFontSizeMultiplier = 4;

%% �ۑ�
figName = 'RangeProfile_16Chs';
print(['Figures/', figName, '.png'],'-dpng','-r0');

%% �I��
toc