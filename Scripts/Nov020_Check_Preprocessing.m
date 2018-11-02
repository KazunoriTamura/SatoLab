%% 開始
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% パラメータ
rangeCr = 7;
rangeCr2 = 8;
rangeMain = 16;
dt = 0.000238;
scalePlotIqForCr = 15e4;
scalePlotIq = 15e4;
c = 299792458;
f_c = 79e9;
% sigma = 1e2;

%% 送信素子
Tx = zeros(4,3);
Tx1 = [-12.90e-3,-4.05e-3,0]; % [x,y,z]座標
Tx2 = [-12.90e-3,-1.35e-3,0];
Tx3 = [-12.90e-3,+1.35e-3,0];
Tx4 = [-12.90e-3,+4.05e-3,0];
% Tx1 = [-21.95e-3,-4.05e-3,0]; % [x,y,z]座標
% Tx2 = [-21.95e-3,-1.35e-3,0];
% Tx3 = [-21.95e-3,+1.35e-3,0];
% Tx4 = [-21.95e-3,+4.05e-3,0];
Tx(1,:) = Tx1;
Tx(2,:) = Tx2;
Tx(3,:) = Tx3;
Tx(4,:) = Tx4;

%% 受信素子
Rx = zeros(4,3);
Rx1 = [+5.00e-3,0,0]; % [x,y,z]座標
Rx2 = [+7.70e-3,0,0];
Rx3 = [10.40e-3,0,0];
Rx4 = [13.10e-3,0,0];
% Rx1 = [-4.05e-3,0,0]; % [x,y,z]座標
% Rx2 = [-1.35e-3,0,0];
% Rx3 = [+1.35e-3,0,0];
% Rx4 = [+4.05e-3,0,0];
Rx(1,:) = Rx1;
Rx(2,:) = Rx2;
Rx(3,:) = Rx3;
Rx(4,:) = Rx4;

%% 素子配列計算
r = zeros(16,3);
for tx = 1:4
    for rx = 1:4
        No = txrx2no(tx,rx);
        r(No,:) = Rx(rx,:) - Tx(tx,:);
    end
end

%% データ
data = data_coh;

%% 前処理
dtCoh = dt * nCohPoints;
n_time_max = size(data,2);
n_time_coh = 0:(n_time_max-1);
m_time_coh = n_time_coh * dtCoh;
lambda = c/f_c;
rawFile = dir('rd79_*.mat');
radarId = rawFile.name(25);

%% １号機用位相反転操作
if radarId == '1'
    data = conj(data);
end

%% 生データIQ平面表示
% CRレンジ
figure(1)
T = ['Raw Data','（Range：',num2str(rangeCr),'）'];
plotIqForT(data,rangeCr,'MaxValue',scalePlotIqForCr,'Title',T);

% CR2レンジ
figure(2)
T = ['Raw Data','（Range：',num2str(rangeCr2),'）'];
plotIqForT(data,rangeCr2,'MaxValue',scalePlotIqForCr,'Title',T);

% 人レンジ
figure(3)
T = ['Raw Data','（Range：',num2str(rangeMain),'）'];
plotIqForT(data,rangeMain,'MaxValue',scalePlotIq,'Title',T);

%% 楕円補正（８月パラメータ）
if radarId == '1'
    load('iqparams_1_180803_v4.mat');
elseif radarId == '2'
    load('iqparams_2_180905.mat');
else
    error('エラー！');
end
data_balanced = compensateImbalance(data,q,phi);
% data_balanced = data;

%% レンジ抽出
dataMain = data_balanced(rangeMain,:,:,:);
dataCr = data_balanced(rangeCr,:,:,:);
dataCr2 = data_balanced(rangeCr2,:,:,:);

%% 各種DCの計算
% CRレンジ
meanDc = getDcByMean(dataCr);
huDc = getDcByHu(dataCr);
modifiedHuDc = getDcByModifiedHu(dataCr);

% CR2レンジ
meanDc2 = getDcByMean(dataCr2);
huDc2 = getDcByHu(dataCr2);
modifiedHuDc2 = getDcByModifiedHu(dataCr2);

% 人レンジ
meanDcMain = getDcByMean(dataMain);


%% 楕円補正後のIQ平面表示
% CRレンジ
figure(100)
T = ['楕円補正後のIQ平面と各種DC','（Range：',num2str(rangeCr),'）'];
plotIqForT(dataCr,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% CR2レンジ
figure(200)
T = ['楕円補正後のIQ平面と各種DC','（Range：',num2str(rangeCr2),'）'];
plotIqForT(dataCr2,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% 人レンジ
figure(300)
T = ['楕円補正後のIQ平面と各種DC','（Range：',num2str(rangeMain),'）'];
plotIqForT(dataMain,rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% 楕円補正後の位相表示
% CRレンジ
figure(110)
T = ['楕円補正後の位相','（Range：',num2str(rangeCr),'）'];
plotAllPhaseForT(dataCr,rangeCr,dtCoh,'Title',T);

% CR2レンジ
figure(210)
T = ['楕円補正後の位相','（Range：',num2str(rangeCr2),'）'];
plotAllPhaseForT(dataCr2,rangeCr2,dtCoh,'Title',T);

% 人レンジ
figure(310)
T = ['楕円補正後の位相','（Range：',num2str(rangeMain),'）'];
plotAllPhaseForT(dataMain,rangeMain,dtCoh,'Title',T);


%% 微調整



%% 修正DCを引く
dataCr_partwodc = reduceInputDc(dataCr,modifiedHuDc);
dataCr2_partwodc = reduceInputDc(dataCr2,modifiedHuDc2);

%% 平均DCを引く
% data_wodc = data_balanced - mean(data_balanced,2);
data_wodc = data_balanced;
data_wodc(rangeCr,:,:,:) = dataCr_partwodc;
data_wodc(rangeCr2,:,:,:) = dataCr2_partwodc;

%% 各種DCの計算（再）
% CRレンジ
meanDc = getDcByMean(dataCr_partwodc);
huDc = getDcByHu(dataCr_partwodc);
modifiedHuDc = getDcByModifiedHu(dataCr_partwodc);

% CR2レンジ
meanDc2 = getDcByMean(dataCr2_partwodc);
huDc2 = getDcByHu(dataCr2_partwodc);
modifiedHuDc2 = getDcByModifiedHu(dataCr2_partwodc);

% 人レンジ
meanDcMain = getDcByMean(data_wodc(rangeMain,:,:,:));

%% 修正DC除去後のIQ平面表示
% CRレンジ
figure(120)
T = ['修正DC除去後のIQ平面と各種DC','（Range：',num2str(rangeCr),'）'];
plotIqForT(dataCr_partwodc,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% CR2レンジ
figure(220)
T = ['修正DC除去後のIQ平面と各種DC','（Range：',num2str(rangeCr2),'）'];
plotIqForT(dataCr2_partwodc,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% 人レンジ
figure(320)
T = ['修正DC除去後のIQ平面と各種DC','（Range：',num2str(rangeMain),'）'];
plotIqForT(data_wodc(rangeMain,:,:,:),rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% 修正DC除去後の位相表示
% CRレンジ
figure(130)
T = ['修正DC除去後の位相','（Range：',num2str(rangeCr),'）'];
plotAllPhaseForT(dataCr_partwodc,rangeCr,dtCoh,'Title',T);

% CR2レンジ
figure(230)
T = ['修正DC除去後の位相','（Range：',num2str(rangeCr2),'）'];
plotAllPhaseForT(dataCr2_partwodc,rangeCr2,dtCoh,'Title',T);

% 人レンジ
figure(230)
T = ['修正DC除去後の位相','（Range：',num2str(rangeMain),'）'];
plotAllPhaseForT(data_wodc(rangeMain,:,:,:),rangeMain,dtCoh,'Title',T);

%% CRレンジからチップ間位相雑音変動を求め，それを全レンジについて除去
phaseNoise = getPhaseNoise(dataCr_partwodc);
data_reduced = reducePhaseNoise( data_wodc, phaseNoise );

%% チップ間雑音除去後のIQ平面表示
% CRレンジ
figure(1000)
T = ['チップ間位相雑音除去後のIQ平面と各種DC','（Range：',num2str(rangeCr),'）'];
plotIqForT(data_reduced,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% CR2レンジ
figure(2000)
T = ['チップ間位相雑音除去後のIQ平面と各種DC','（Range：',num2str(rangeCr2),'）'];
plotIqForT(data_reduced,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

% 人レンジ
figure(3000)
T = ['チップ間位相雑音除去後のIQ平面と各種DC','（Range：',num2str(rangeMain),'）'];
plotIqForT(data_reduced,rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% チップ間雑音除去後の位相表示
% CRレンジ
figure(1100)
T = ['チップ間位相雑音除去後の位相','（Range：',num2str(rangeCr),'）'];
plotAllPhaseForT(data_reduced,rangeCr,dtCoh,'Title',T);

% CR2レンジ
figure(2100)
T = ['チップ間位相雑音除去後の位相','（Range：',num2str(rangeCr2),'）'];
plotAllPhaseForT(data_reduced,rangeCr2,dtCoh,'Title',T);

% 人レンジ
figure(3100)
T = ['チップ間位相雑音除去後の位相','（Range：',num2str(rangeMain),'）'];
plotAllPhaseForT(data_reduced,rangeMain,dtCoh,'Title',T);

%% データ再構成
data = data_reduced;

% %% CRの位相算出
% angle_rCR = unwrap(angle(data(rangeCr,:,:,:)));
% hosei_tmp = reshapeTxrx2No(angle_rCR);
% hosei_tmp2 = (squeeze(hosei_tmp)).';
% hosei_mean = mean(hosei_tmp2,2);
% hosei = repmat(hosei_mean,1,size(angle_rCR,2));
% 
% %% データ
% X_tmp = reshapeTxrx2No(data(range,:,:,:));
% X = (squeeze(X_tmp)).';
% X_hosei = zeros(16,size(X,2));
% phi = -1i * hosei;
% for k = 1:16
%     X_hosei(k,:) = X(k,:) .* exp(phi(k,:));
% end

%% キャリブレーション
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

%% キャリブレーション後の位相表示
% CRレンジ
figure(1300)
T = ['キャリブレーション後の位相','（Range：',num2str(rangeCr),'）'];
plotAllPhaseForT(data_calibrated,rangeCr,dtCoh,'Title',T);

% CR2レンジ
figure(2300)
T = ['キャリブレーション後の位相','（Range：',num2str(rangeCr2),'）'];
plotAllPhaseForT(data_calibrated,rangeCr2,dtCoh,'Title',T);

% 人レンジ
figure(3300)
T = ['キャリブレーション後の位相','（Range：',num2str(rangeMain),'）'];
plotAllPhaseForT(data_calibrated,rangeMain,dtCoh,'Title',T);


%%
toc