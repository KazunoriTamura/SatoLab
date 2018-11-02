function [] = plotAllMagForT( data, rangeNo, Tr, varargin)
% 複素時系列の振幅を表示
% 実部と虚部に分けて表示する

% 入力データ形式
% data(RANGE, SAMPLE, TX, RX) : 4dims
% rangeNo : 整数
% Tr : double
% <option>'PlotRange', [minAngle maxAngle]
% <option>'Title', 'Graph Title'

fig = gcf;
fig.OuterPosition = [100,100,2500,1300];

D = ndims(data);
R = size(data,1);

if D ~= 4
    error('入力は４次元で！')
end

%%
rageNoForTitle = rangeNo;

%%
if R == 1
    rangeNo = 1;
end

%% データ抽出
data_r = data(rangeNo,:,:,:);
magnitude = abs(data_r);

%% 位相の最大/最小値
maxMagnitude = ceil(max(max(max(magnitude))));
minMagnitude = floor(min(min(min(magnitude))));
% maxImagMagnitude = ceil(max(max(max(imagData))));
% minImagMagnitude = floor(min(min(min(imagData))));
% maxPlotMagnitude = max(maxRealMagnitude,maxImagMagnitude);
% minPlotMagnitude = min(minRealMagnitude,minImagMagnitude);
defaultPlotRange = [minMagnitude maxMagnitude];

%% 引数のパース
p = inputParser;

%% 必須引数
addRequired(p,'data',@isnumeric);
addRequired(p,'rangeNo',@isnumeric);
addRequired(p,'Tr',@isnumeric);

%% 可変引数
defaultTitle = ['Range: ', num2str(rageNoForTitle)];
addOptional(p, 'PlotRange', defaultPlotRange, @isnumeric)
addOptional(p, 'Title', defaultTitle)

parse(p, data, rangeNo, Tr, varargin{:})

%% 時刻
t = (1: size(data, 2)) * Tr;

%% 4Ch/16Ch分岐
switch size(data, 3)
    case 4
        % 上段
        subplot(2,2,1)
        plot(t,magnitude(:,:,1,1), t,magnitude(:,:,1,2), t,magnitude(:,:,2,1), t,magnitude(:,:,2,2));
        ylim(p.Results.PlotRange)
        %         xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx1-Rx1', 'Tx1-Rx2', 'Tx2-Rx1', 'Tx2-Rx2')
        xlabel('Time (s)')
        ylabel('Magnitude')
        %         title('Real Part')
        ax = gca;
        ax.FontSize = 30;
        ax.XLabel.FontSize = 40;
        ax.YLabel.FontSize = 40;
        
        subplot(2,2,2)
        plot(t,magnitude(:,:,1,3), t,magnitude(:,:,1,4), t,magnitude(:,:,2,3), t,magnitude(:,:,2,4));
        ylim(p.Results.PlotRange)
        %         xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx1-Rx3', 'Tx1-Rx4', 'Tx2-Rx3', 'Tx2-Rx4')
        xlabel('Time (s)')
        ylabel('Magnitude')
        %         title('Real Part')
        ax = gca;
        ax.FontSize = 30;
        ax.XLabel.FontSize = 40;
        ax.YLabel.FontSize = 40;
        
        
        % 下段
        subplot(2,2,3)
        plot(t,magnitude(:,:,3,1), t,magnitude(:,:,3,2), t,magnitude(:,:,4,1), t,magnitude(:,:,4,2));
        ylim(p.Results.PlotRange)
        %         xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx3-Rx1', 'Tx3-Rx2', 'Tx4-Rx1', 'Tx4-Rx2')
        xlabel('Time (s)')
        ylabel('Magnitude')
        %         title('Real Part')
        ax = gca;
        ax.FontSize = 30;
        ax.XLabel.FontSize = 40;
        ax.YLabel.FontSize = 40;
        
        subplot(2,2,4)
        plot(t,magnitude(:,:,3,3), t,magnitude(:,:,3,4), t,magnitude(:,:,4,3), t,magnitude(:,:,4,4));
        ylim(p.Results.PlotRange)
        %         xlim([0 maxTime])
        grid on
        grid minor
        legend('Tx3-Rx3', 'Tx3-Rx4', 'Tx4-Rx3', 'Tx4-Rx4')
        xlabel('Time (s)')
        ylabel('Magnitude')
        %         title('Real Part')
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
        
        %     case 2
        %         plot(t, angleOfAllChs(:, 1, 1),t, angleOfAllChs(:, 1, 2),t, angleOfAllChs(:, 2, 1),t, angleOfAllChs(:, 2, 2))
        %         ylim(p.Results.PlotRange)
        %         xlim([0 maxTime])
        %         grid on
        %         grid minor
        %         legend('Tx1-Rx1', 'Tx1-Rx2', 'Tx2-Rx1', 'Tx2-Rx2')
        %         xlabel('Time (s)')
        %         ylabel('Phase (rad)')
        %         title({p.Results.Title; ''})
end