function [] = plotIqForT(data, rangeNo, varargin)
%% plotIq(data, rangeNo, <option>maxValue) -- 第rangeNoレンジについて、IQ平面で全Ch表示
% 入力データ形式
% data(RANGE, SAMPLE, TX, RX) : 4dims
% rangeNo : 整数
% <option>'MaxValue', maxValue
% <option>'Title', 'Graph Title'
% <option>'MeanDc', meanDc(4, 4)
% <option>'HuDc', huDc(4, 4)
% <option>'ModifiedHuDc', modifiedHuDc(4, 4)
% <option>'MarkerSize', makerSize (default: 2)
% <option>'LineWidth', lineWidth (default: 2)
% 行方向Rx, 列方向Txで表示

fig = gcf;
fig.OuterPosition = [100,100,1600,1300];

nCols = size(data, 4); %行subplot数
nRows = size(data, 3); %列subplot数

% タイトル用レンジ番号の決定
rangeNoForTitle = rangeNo;

% １レンジ分のデータに対する処理
if size(data,1) == 1
    rangeNo = 1;
end

% 振幅の最大値を1E4刻みで取得
maxAmplitude = max(max(max(abs(data(rangeNo, :, :, :))))) / 1e4;
maxAmplitude = ceil(maxAmplitude) * 1e4;
defaultMaxValue = maxAmplitude;

% DC配列初期値
defaultMeanDc = NaN(nCols, nRows);
defaultHuDc = NaN(nCols, nRows);
defaultModifiedHuDc = NaN(nCols, nRows);

%% 引数のパース
p = inputParser;

% 必須引数
addRequired(p,'data',@isnumeric);
addRequired(p,'rangeNo',@isnumeric);

% 可変引数
defaultTitle = ['Range：', num2str(rangeNoForTitle)];
defaultMarkerSize = 2;
defaultLineWidth = 2;
addOptional(p, 'MaxValue', defaultMaxValue, @isnumeric)
addOptional(p, 'Title', defaultTitle)
addOptional(p, 'MeanDc', defaultMeanDc, @isnumeric)
addOptional(p, 'HuDc', defaultHuDc, @isnumeric)
addOptional(p, 'ModifiedHuDc', defaultModifiedHuDc, @isnumeric)
addOptional(p, 'MarkerSize', defaultMarkerSize, @isnumeric)
addOptional(p, 'LineWidth', defaultLineWidth, @isnumeric)

parse(p, data, rangeNo, varargin{:})

meanDc = p.Results.MeanDc;
huDc = p.Results.HuDc;
modifiedHuDc = p.Results.ModifiedHuDc;

%% Plot
for iTx = 1 : nRows
    for iRx = 1 : nCols
        plotNo = (iTx - 1) * nCols + iRx;
        
        subplot(nCols, nRows, plotNo)
        plot(data(rangeNo, :, iTx, iRx), '.')
        hold on
        if ~any(isnan(meanDc))
            plot(meanDc(iTx,iRx), 'ko', 'MarkerSize', p.Results.MarkerSize)
        end
        
        if ~any(isnan(huDc))
            plot(huDc(iTx,iRx), 'ro', 'MarkerSize', p.Results.MarkerSize)
        end

        if ~any(isnan(modifiedHuDc))
            plot(modifiedHuDc(iTx,iRx), 'rx', 'MarkerSize', p.Results.MarkerSize)
        end
        axis square
        xlim([-p.Results.MaxValue p.Results.MaxValue])
        ylim([-p.Results.MaxValue p.Results.MaxValue])
        grid on
        grid minor
        title(strcat('Tx',num2str(iTx),'-Rx',num2str(iRx)))
        ax = gca;
        ax.XLabel.String = 'In-Phase';
        ax.YLabel.String = 'Quadrature-Phase';
        hold off
    end
end

%% 全体のタイトルを設定
axes;
title({p.Results.Title; ''});
axis off;

ax = gca;
ax.TitleFontSizeMultiplier = 2;

end