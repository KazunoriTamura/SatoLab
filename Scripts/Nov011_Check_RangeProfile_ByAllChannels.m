%% 開始
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% レンジプロファイル（１チャネル）の表示
figure(2)
fig = gcf;
fig.OuterPosition = [100,100,2103,1300];
T = ['Range Profile','（16チャネル平均）'];
p = plotRangeProfile(data_coh,'UseCh',1:16,'Title',T);
set(p,'LineWidth',4);
ax = gca;
ax.XAxis.FontSize = 40;
ax.YAxis.FontSize = 40;
ax.Legend.FontSize = 50;
ax.TitleFontSizeMultiplier = 4;

%% 保存
figName = 'RangeProfile_16Chs';
print(['Figures/', figName, '.png'],'-dpng','-r0');

%% 終了
toc