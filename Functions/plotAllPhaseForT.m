function [] = plotAllPhaseForT(data, rangeNo, Tr, varargin)
%% plotAllPhase(data, rangeNo, Tr, <option>) -- サンプリング間隔Trのdataの第rangeNoレンジについて，全チャネルの位相を表示
% 入力データ形式
% data(RANGE, SAMPLE, TX, RX) : 4dims
% rangeNo : 整数
% Tr : double
% <option>'PlotRange', [minAngle maxAngle]
% <option>'Title', 'Graph Title'

fig = gcf;
fig.OuterPosition = [100,100,2500,1300];

% タイトル用レンジ番号の決定
rangeNoForTitle = rangeNo;

% １レンジ分のデータに対する処理
if size(data,1) == 1
    rangeNo = 1;
end

% radに変換
angleOfAllChs = zeros(size(data, 2), size(data, 3), size(data, 4));
for iTx = 1 : size(data, 3)
    for iRx = 1 : size(data, 4)
        angleOfAllChs(:, iTx, iRx) = unwrap(angle(data(rangeNo, :, iTx, iRx)));
    end
end

% 位相の最大/最小値
maxAngle = ceil(max(max(max(angleOfAllChs))));
minAngle = floor(min(min(min(angleOfAllChs))));
defaultPlotRange = [minAngle maxAngle];

%% 引数のパース
p = inputParser;

% 必須引数
addRequired(p,'data',@isnumeric);
addRequired(p,'rangeNo',@isnumeric);
addRequired(p,'Tr',@isnumeric);

% 可変引数
defaultTitle = ['Range:', num2str(rangeNoForTitle)];
addOptional(p, 'PlotRange', defaultPlotRange, @isnumeric)
addOptional(p, 'Title', defaultTitle)

parse(p, data, rangeNo, Tr, varargin{:})


% 時刻
t = (1: size(data, 2)) * Tr;

% 表示する時刻の最大値
maxTime = ceil(t(end));

% 4Ch/16Ch分岐
switch size(data, 3)
    case 4
        subplot(2,2,1)
        plot(t, angleOfAllChs(:, 1, 1),t, angleOfAllChs(:, 1, 2),t, angleOfAllChs(:, 2, 1),t, angleOfAllChs(:, 2, 2))
        ylim(p.Results.PlotRange)
        xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx1-Rx1', 'Tx1-Rx2', 'Tx2-Rx1', 'Tx2-Rx2')
        xlabel('Time (s)')
        ylabel('Phase (rad)')
        ax = gca;
        ax.FontSize = 30;
        ax.XLabel.FontSize = 40;
        ax.YLabel.FontSize = 40;
        
        subplot(2,2,2)
        plot(t, angleOfAllChs(:, 1, 3),t, angleOfAllChs(:, 1, 4),t, angleOfAllChs(:, 2, 3),t, angleOfAllChs(:, 2, 4))
        ylim(p.Results.PlotRange)
        xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx1-Rx3', 'Tx1-Rx4', 'Tx2-Rx3', 'Tx2-Rx4')
        xlabel('Time (s)')
        ylabel('Phase (rad)')
        ax = gca;
        ax.FontSize = 30;
        ax.XLabel.FontSize = 40;
        ax.YLabel.FontSize = 40;
        
        subplot(2,2,3)
        plot(t, angleOfAllChs(:, 3, 1),t, angleOfAllChs(:, 3, 2),t, angleOfAllChs(:, 4, 1),t, angleOfAllChs(:, 4,2))
        ylim(p.Results.PlotRange)
        xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx3-Rx1', 'Tx3-Rx2', 'Tx4-Rx1', 'Tx4-Rx2')
        xlabel('Time (s)')
        ylabel('Phase (rad)')
        ax = gca;
        ax.FontSize = 30;
        ax.XLabel.FontSize = 40;
        ax.YLabel.FontSize = 40;
        
        subplot(2,2,4)
        plot(t, angleOfAllChs(:, 3, 3),t, angleOfAllChs(:, 3, 4),t, angleOfAllChs(:, 4, 3),t, angleOfAllChs(:, 4, 4))
        ylim(p.Results.PlotRange)
        xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx3-Rx3', 'Tx3-Rx4', 'Tx4-Rx3', 'Tx4-Rx4')
        xlabel('Time (s)')
        ylabel('Phase (rad)')
        ax = gca;
        ax.FontSize = 30;
        ax.XLabel.FontSize = 40;
        ax.YLabel.FontSize = 40;
        
        % 全体のタイトルを設定
        axes;
        title({p.Results.Title; ''});
        ax = gca;
        ax.TitleFontSizeMultiplier = 3;
        axis off;
        
        
    case 2
        plot(t, angleOfAllChs(:, 1, 1),t, angleOfAllChs(:, 1, 2),t, angleOfAllChs(:, 2, 1),t, angleOfAllChs(:, 2, 2))
        ylim(p.Results.PlotRange)
        xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx1-Rx1', 'Tx1-Rx2', 'Tx2-Rx1', 'Tx2-Rx2')
        xlabel('Time (s)')
        ylabel('Phase (rad)')
        title(p.Results.Title)
        
end