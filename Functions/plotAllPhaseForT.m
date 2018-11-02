function [] = plotAllPhaseForT(data, rangeNo, Tr, varargin)
%% plotAllPhase(data, rangeNo, Tr, <option>) -- �T���v�����O�ԊuTr��data�̑�rangeNo�����W�ɂ��āC�S�`���l���̈ʑ���\��
% ���̓f�[�^�`��
% data(RANGE, SAMPLE, TX, RX) : 4dims
% rangeNo : ����
% Tr : double
% <option>'PlotRange', [minAngle maxAngle]
% <option>'Title', 'Graph Title'

fig = gcf;
fig.OuterPosition = [100,100,2500,1300];

% �^�C�g���p�����W�ԍ��̌���
rangeNoForTitle = rangeNo;

% �P�����W���̃f�[�^�ɑ΂��鏈��
if size(data,1) == 1
    rangeNo = 1;
end

% rad�ɕϊ�
angleOfAllChs = zeros(size(data, 2), size(data, 3), size(data, 4));
for iTx = 1 : size(data, 3)
    for iRx = 1 : size(data, 4)
        angleOfAllChs(:, iTx, iRx) = unwrap(angle(data(rangeNo, :, iTx, iRx)));
    end
end

% �ʑ��̍ő�/�ŏ��l
maxAngle = ceil(max(max(max(angleOfAllChs))));
minAngle = floor(min(min(min(angleOfAllChs))));
defaultPlotRange = [minAngle maxAngle];

%% �����̃p�[�X
p = inputParser;

% �K�{����
addRequired(p,'data',@isnumeric);
addRequired(p,'rangeNo',@isnumeric);
addRequired(p,'Tr',@isnumeric);

% �ψ���
defaultTitle = ['Range:', num2str(rangeNoForTitle)];
addOptional(p, 'PlotRange', defaultPlotRange, @isnumeric)
addOptional(p, 'Title', defaultTitle)

parse(p, data, rangeNo, Tr, varargin{:})


% ����
t = (1: size(data, 2)) * Tr;

% �\�����鎞���̍ő�l
maxTime = ceil(t(end));

% 4Ch/16Ch����
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
        
        % �S�̂̃^�C�g����ݒ�
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