%% 開始
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% パラメータ
Tx = 1;
Rx = 1;

%% レンジプロファイル（１チャネル）の表示
figure(1)
fig = figure(1);
fig.OuterPosition = [300,300,2103,1300];
T = ['Range Profile','（ Tx',num2str(Tx),'-Rx',num2str(Rx),' ）'];
p = plotRangeProfile(data_coh,'UseTx',Tx,'UseRx',Rx,'Title',T);
set(p,'LineWidth',4);
ax = gca;
ax.XAxis.FontSize = 40;
ax.YAxis.FontSize = 40;
ax.Legend.FontSize = 50;
ax.TitleFontSizeMultiplier = 4;

%% 保存
figName = ['RangeProfile_','Tx',num2str(Tx),'_Rx',num2str(Rx)];
print(['Figures/', figName, '.png'],'-dpng','-r0');

%% 終了
toc