%% �J�n
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% �p�����[�^
rangeCr = 7;
rangeCr2 = 8;
rangeMain = 16;
dt = 0.000238;
scalePlotIqForCr = 15e4;
scalePlotIq = 15e4;
c = 299792458;
f_c = 79e9;
% sigma = 1e2;

%% ���M�f�q
Tx = zeros(4,3);
Tx1 = [-12.90e-3,-4.05e-3,0]; % [x,y,z]���W
Tx2 = [-12.90e-3,-1.35e-3,0];
Tx3 = [-12.90e-3,+1.35e-3,0];
Tx4 = [-12.90e-3,+4.05e-3,0];
% Tx1 = [-21.95e-3,-4.05e-3,0]; % [x,y,z]���W
% Tx2 = [-21.95e-3,-1.35e-3,0];
% Tx3 = [-21.95e-3,+1.35e-3,0];
% Tx4 = [-21.95e-3,+4.05e-3,0];
Tx(1,:) = Tx1;
Tx(2,:) = Tx2;
Tx(3,:) = Tx3;
Tx(4,:) = Tx4;

%% ��M�f�q
Rx = zeros(4,3);
Rx1 = [+5.00e-3,0,0]; % [x,y,z]���W
Rx2 = [+7.70e-3,0,0];
Rx3 = [10.40e-3,0,0];
Rx4 = [13.10e-3,0,0];
% Rx1 = [-4.05e-3,0,0]; % [x,y,z]���W
% Rx2 = [-1.35e-3,0,0];
% Rx3 = [+1.35e-3,0,0];
% Rx4 = [+4.05e-3,0,0];
Rx(1,:) = Rx1;
Rx(2,:) = Rx2;
Rx(3,:) = Rx3;
Rx(4,:) = Rx4;

%% �f�q�z��v�Z
r = zeros(16,3);
for tx = 1:4
    for rx = 1:4
        No = txrx2no(tx,rx);
        r(No,:) = Rx(rx,:) - Tx(tx,:);
    end
end

%% �f�[�^
data = data_coh;

%% �O����
dtCoh = dt * nCohPoints;
n_time_max = size(data,2);
n_time_coh = 0:(n_time_max-1);
m_time_coh = n_time_coh * dtCoh;
lambda = c/f_c;
rawFile = dir('rd79_*.mat');
radarId = rawFile.name(25);

%% �P���@�p�ʑ����]����
if radarId == '1'
    data = conj(data);
end

%% ���f�[�^IQ���ʕ\��
% CR�����W
figure(1)
T = ['Raw Data','�iRange�F',num2str(rangeCr),'�j'];
plotIqForT(data,rangeCr,'MaxValue',scalePlotIqForCr,'Title',T);

% CR2�����W
figure(2)
T = ['Raw Data','�iRange�F',num2str(rangeCr2),'�j'];
plotIqForT(data,rangeCr2,'MaxValue',scalePlotIqForCr,'Title',T);

% �l�����W
figure(3)
T = ['Raw Data','�iRange�F',num2str(rangeMain),'�j'];
plotIqForT(data,rangeMain,'MaxValue',scalePlotIq,'Title',T);

%% �ȉ~�␳�i�W���p�����[�^�j
if radarId == '1'
    load('iqparams_1_180803_v4.mat');
elseif radarId == '2'
    load('iqparams_2_180905.mat');
else
    error('�G���[�I');
end
data_balanced = compensateImbalance(data,q,phi);
% data_balanced = data;

%% �����W���o
dataMain = data_balanced(rangeMain,:,:,:);
dataCr = data_balanced(rangeCr,:,:,:);
dataCr2 = data_balanced(rangeCr2,:,:,:);

%% �e��DC�̌v�Z
% CR�����W
meanDc = getDcByMean(dataCr);
huDc = getDcByHu(dataCr);
modifiedHuDc = getDcByModifiedHu(dataCr);

% CR2�����W
meanDc2 = getDcByMean(dataCr2);
huDc2 = getDcByHu(dataCr2);
modifiedHuDc2 = getDcByModifiedHu(dataCr2);

% �l�����W
meanDcMain = getDcByMean(dataMain);


%% �ȉ~�␳���IQ���ʕ\��
% CR�����W
figure(100)
T = ['�ȉ~�␳���IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr),'�j'];
plotIqForT(dataCr,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% CR2�����W
figure(200)
T = ['�ȉ~�␳���IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr2),'�j'];
plotIqForT(dataCr2,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% �l�����W
figure(300)
T = ['�ȉ~�␳���IQ���ʂƊe��DC','�iRange�F',num2str(rangeMain),'�j'];
plotIqForT(dataMain,rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% �ȉ~�␳��̈ʑ��\��
% CR�����W
figure(110)
T = ['�ȉ~�␳��̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
plotAllPhaseForT(dataCr,rangeCr,dtCoh,'Title',T);

% CR2�����W
figure(210)
T = ['�ȉ~�␳��̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
plotAllPhaseForT(dataCr2,rangeCr2,dtCoh,'Title',T);

% �l�����W
figure(310)
T = ['�ȉ~�␳��̈ʑ�','�iRange�F',num2str(rangeMain),'�j'];
plotAllPhaseForT(dataMain,rangeMain,dtCoh,'Title',T);


%% ������



%% �C��DC������
dataCr_partwodc = reduceInputDc(dataCr,modifiedHuDc);
dataCr2_partwodc = reduceInputDc(dataCr2,modifiedHuDc2);

%% ����DC������
% data_wodc = data_balanced - mean(data_balanced,2);
data_wodc = data_balanced;
data_wodc(rangeCr,:,:,:) = dataCr_partwodc;
data_wodc(rangeCr2,:,:,:) = dataCr2_partwodc;

%% �e��DC�̌v�Z�i�āj
% CR�����W
meanDc = getDcByMean(dataCr_partwodc);
huDc = getDcByHu(dataCr_partwodc);
modifiedHuDc = getDcByModifiedHu(dataCr_partwodc);

% CR2�����W
meanDc2 = getDcByMean(dataCr2_partwodc);
huDc2 = getDcByHu(dataCr2_partwodc);
modifiedHuDc2 = getDcByModifiedHu(dataCr2_partwodc);

% �l�����W
meanDcMain = getDcByMean(data_wodc(rangeMain,:,:,:));

%% �C��DC�������IQ���ʕ\��
% CR�����W
figure(120)
T = ['�C��DC�������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr),'�j'];
plotIqForT(dataCr_partwodc,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% CR2�����W
figure(220)
T = ['�C��DC�������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr2),'�j'];
plotIqForT(dataCr2_partwodc,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% �l�����W
figure(320)
T = ['�C��DC�������IQ���ʂƊe��DC','�iRange�F',num2str(rangeMain),'�j'];
plotIqForT(data_wodc(rangeMain,:,:,:),rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% �C��DC������̈ʑ��\��
% CR�����W
figure(130)
T = ['�C��DC������̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
plotAllPhaseForT(dataCr_partwodc,rangeCr,dtCoh,'Title',T);

% CR2�����W
figure(230)
T = ['�C��DC������̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
plotAllPhaseForT(dataCr2_partwodc,rangeCr2,dtCoh,'Title',T);

% �l�����W
figure(230)
T = ['�C��DC������̈ʑ�','�iRange�F',num2str(rangeMain),'�j'];
plotAllPhaseForT(data_wodc(rangeMain,:,:,:),rangeMain,dtCoh,'Title',T);

%% CR�����W����`�b�v�Ԉʑ��G���ϓ������߁C�����S�����W�ɂ��ď���
phaseNoise = getPhaseNoise(dataCr_partwodc);
data_reduced = reducePhaseNoise( data_wodc, phaseNoise );

%% �`�b�v�ԎG���������IQ���ʕ\��
% CR�����W
figure(1000)
T = ['�`�b�v�Ԉʑ��G���������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr),'�j'];
plotIqForT(data_reduced,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% CR2�����W
figure(2000)
T = ['�`�b�v�Ԉʑ��G���������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr2),'�j'];
plotIqForT(data_reduced,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% �l�����W
figure(3000)
T = ['�`�b�v�Ԉʑ��G���������IQ���ʂƊe��DC','�iRange�F',num2str(rangeMain),'�j'];
plotIqForT(data_reduced,rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% �`�b�v�ԎG��������̈ʑ��\��
% CR�����W
figure(1100)
T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
plotAllPhaseForT(data_reduced,rangeCr,dtCoh,'Title',T);

% CR2�����W
figure(2100)
T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
plotAllPhaseForT(data_reduced,rangeCr2,dtCoh,'Title',T);

% �l�����W
figure(3100)
T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(rangeMain),'�j'];
plotAllPhaseForT(data_reduced,rangeMain,dtCoh,'Title',T);

%% �f�[�^�č\��
data = data_reduced;

% %% CR�̈ʑ��Z�o
% angle_rCR = unwrap(angle(data(rangeCr,:,:,:)));
% hosei_tmp = reshapeTxrx2No(angle_rCR);
% hosei_tmp2 = (squeeze(hosei_tmp)).';
% hosei_mean = mean(hosei_tmp2,2);
% hosei = repmat(hosei_mean,1,size(angle_rCR,2));
% 
% %% �f�[�^
% X_tmp = reshapeTxrx2No(data(range,:,:,:));
% X = (squeeze(X_tmp)).';
% X_hosei = zeros(16,size(X,2));
% phi = -1i * hosei;
% for k = 1:16
%     X_hosei(k,:) = X(k,:) .* exp(phi(k,:));
% end

%% �L�����u���[�V����
% data_reshaped = reshapeTxrx2No(data);
% 
% X = squeeze(data_reshaped(range,:,:)).';
% X_Cr = squeeze(data_reshaped(rangeCr,:,:)).';
% 
% angle_Cr = unwrap(angle(X_Cr));
% phi = mean(angle_Cr,2);
% 
% X_calibrated = X .* exp( -1i * phi );
% 
% % X_calibrated_2 = zeros(size(X_calibrated));
% % for k = 1:16
% %     X_calibrated_2(k,:) = X(k,:) .* exp( -1i * phi(k,:) );
% % end
% % if isequal( X_calibrated , X_calibrated_2 )
% %     disp('OK');
% % end
% 
% output_reshaped = zeros(1,size(X_calibrated,2),size(X_calibrated,1));
% output_reshaped(1,:,:) = X_calibrated.';
% 
% output = reshapeNo2Txrx(output_reshaped);


data_calibrated = calibration(data,rangeCr);

%% �L�����u���[�V������̈ʑ��\��
% CR�����W
figure(1300)
T = ['�L�����u���[�V������̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
plotAllPhaseForT(data_calibrated,rangeCr,dtCoh,'Title',T);

% CR2�����W
figure(2300)
T = ['�L�����u���[�V������̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
plotAllPhaseForT(data_calibrated,rangeCr2,dtCoh,'Title',T);

% �l�����W
figure(3300)
T = ['�L�����u���[�V������̈ʑ�','�iRange�F',num2str(rangeMain),'�j'];
plotAllPhaseForT(data_calibrated,rangeMain,dtCoh,'Title',T);


%%
toc