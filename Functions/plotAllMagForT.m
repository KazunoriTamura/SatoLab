function [] = plotAllMagForT( data, rangeNo, Tr, varargin)
% ���f���n��̐U����\��
% �����Ƌ����ɕ����ĕ\������

% ���̓f�[�^�`��
% data(RANGE, SAMPLE, TX, RX) : 4dims
% rangeNo : ����
% Tr : double
% <option>'PlotRange', [minAngle maxAngle]
% <option>'Title', 'Graph Title'

fig = gcf;
fig.OuterPosition = [100,100,2500,1300];

D = ndims(data);
R = size(data,1);

if D ~= 4
    error('���͂͂S�����ŁI')
end

%%
rageNoForTitle = rangeNo;

%%
if R == 1
    rangeNo = 1;
end

%% �f�[�^���o
data_r = data(rangeNo,:,:,:);
magnitude = abs(data_r);

%% �ʑ��̍ő�/�ŏ��l
maxMagnitude = ceil(max(max(max(magnitude))));
minMagnitude = floor(min(min(min(magnitude))));
% maxImagMagnitude = ceil(max(max(max(imagData))));
% minImagMagnitude = floor(min(min(min(imagData))));
% maxPlotMagnitude = max(maxRealMagnitude,maxImagMagnitude);
% minPlotMagnitude = min(minRealMagnitude,minImagMagnitude);
defaultPlotRange = [minMagnitude maxMagnitude];

%% �����̃p�[�X
p = inputParser;

%% �K�{����
addRequired(p,'data',@isnumeric);
addRequired(p,'rangeNo',@isnumeric);
addRequired(p,'Tr',@isnumeric);

%% �ψ���
defaultTitle = ['Range: ', num2str(rageNoForTitle)];
addOptional(p, 'PlotRange', defaultPlotRange, @isnumeric)
addOptional(p, 'Title', defaultTitle)

parse(p, data, rangeNo, Tr, varargin{:})

%% ����
t = (1: size(data, 2)) * Tr;

%% 4Ch/16Ch����
switch size(data, 3)
    case 4
        % ��i
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
        
        
        % ���i
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
        
        
        
        % �S�̂̃^�C�g����ݒ�
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