%% 開始
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% パラメータ
rangeCr = 7;
rangeCr2 = 9;
rangeMain = 16;
dt = 0.241e-3;
scalePlotIqForCr = 15e4;
scalePlotIq = 5e4;
c = 299792458;
f_c = 79e9;
sigma = 1e4;
useOnly8Chs = 0;    % 1=on, 0=off

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
Tx(1,:) = Tx4;
Tx(2,:) = Tx3;
Tx(3,:) = Tx2;
Tx(4,:) = Tx1;

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
data = data_coh(:,1:size(data_coh,2)/5,:,:);

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
% % CRレンジ
% figure(1)
% T = ['Raw Data','（Range：',num2str(rangeCr),'）'];
% plotIqForT(data,rangeCr,'MaxValue',scalePlotIqForCr,'Title',T);
% 
% % CR2レンジ
% figure(2)
% T = ['Raw Data','（Range：',num2str(rangeCr2),'）'];
% plotIqForT(data,rangeCr2,'MaxValue',scalePlotIqForCr,'Title',T);
% 
% % 人レンジ
% figure(3)
% T = ['Raw Data','（Range：',num2str(rangeMain),'）'];
% plotIqForT(data,rangeMain,'MaxValue',scalePlotIq,'Title',T);

%% 楕円補正
if radarId == '1'
    load('iqparams_1_180803_v4.mat');
elseif radarId == '2'
    load('iqparams_2_180905.mat');
else
    error('エラー！');
end
data_balanced = compensateImbalance(data,q,phi);
data_balanced = data;

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
% % CRレンジ
% figure(100)
% T = ['楕円補正後のIQ平面と各種DC','（Range：',num2str(rangeCr),'）'];
% plotIqForT(dataCr,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % CR2レンジ
% figure(200)
% T = ['楕円補正後のIQ平面と各種DC','（Range：',num2str(rangeCr2),'）'];
% plotIqForT(dataCr2,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % 人レンジ
% figure(300)
% T = ['楕円補正後のIQ平面と各種DC','（Range：',num2str(rangeMain),'）'];
% plotIqForT(dataMain,rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);
% 
% %% 楕円補正後の位相表示
% % CRレンジ
% figure(110)
% T = ['楕円補正後の位相','（Range：',num2str(rangeCr),'）'];
% plotAllPhaseForT(dataCr,rangeCr,dtCoh,'Title',T);
% 
% % CR2レンジ
% figure(210)
% T = ['楕円補正後の位相','（Range：',num2str(rangeCr2),'）'];
% plotAllPhaseForT(dataCr2,rangeCr2,dtCoh,'Title',T);
% 
% % 人レンジ
% figure(310)
% T = ['楕円補正後の位相','（Range：',num2str(rangeMain),'）'];
% plotAllPhaseForT(dataMain,rangeMain,dtCoh,'Title',T);


%% 微調整



%% 修正DCを引く
dataCr_partwodc = reduceInputDc(dataCr,modifiedHuDc);
dataCr2_partwodc = reduceInputDc(dataCr2,modifiedHuDc2);

%% 平均DCを引く
data_wodc = data_balanced - mean(data_balanced,2);
% data_wodc = data_balanced;
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
% % CRレンジ
% figure(120)
% T = ['修正DC除去後のIQ平面と各種DC','（Range：',num2str(rangeCr),'）'];
% plotIqForT(dataCr_partwodc,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % CR2レンジ
% figure(220)
% T = ['修正DC除去後のIQ平面と各種DC','（Range：',num2str(rangeCr2),'）'];
% plotIqForT(dataCr2_partwodc,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % 人レンジ
% figure(320)
% T = ['修正DC除去後のIQ平面と各種DC','（Range：',num2str(rangeMain),'）'];
% plotIqForT(data_wodc(rangeMain,:,:,:),rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);
 
%% 修正DC除去後の位相表示
% % CRレンジ
% figure(130)
% T = ['修正DC除去後の位相','（Range：',num2str(rangeCr),'）'];
% plotAllPhaseForT(dataCr_partwodc,rangeCr,dtCoh,'Title',T);
% 
% % CR2レンジ
% figure(230)
% T = ['修正DC除去後の位相','（Range：',num2str(rangeCr2),'）'];
% plotAllPhaseForT(dataCr2_partwodc,rangeCr2,dtCoh,'Title',T);
% 
% % 人レンジ
% figure(230)
% T = ['修正DC除去後の位相','（Range：',num2str(rangeMain),'）'];
% plotAllPhaseForT(data_wodc(rangeMain,:,:,:),rangeMain,dtCoh,'Title',T);

%% CRレンジからチップ間位相雑音変動を求め，それを全レンジについて除去
phaseNoise = getPhaseNoise(dataCr_partwodc);
data_reduced = reducePhaseNoise( data_wodc, phaseNoise );

%% チップ間雑音除去後のIQ平面表示
% % CRレンジ
% figure(1000)
% T = ['チップ間位相雑音除去後のIQ平面と各種DC','（Range：',num2str(rangeCr),'）'];
% plotIqForT(data_reduced,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % CR2レンジ
% figure(2000)
% T = ['チップ間位相雑音除去後のIQ平面と各種DC','（Range：',num2str(rangeCr2),'）'];
% plotIqForT(data_reduced,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % 人レンジ
% figure(3000)
% T = ['チップ間位相雑音除去後のIQ平面と各種DC','（Range：',num2str(rangeMain),'）'];
% plotIqForT(data_reduced,rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% チップ間雑音除去後の位相表示
% % CRレンジ
% figure(1100)
% T = ['チップ間位相雑音除去後の位相','（Range：',num2str(rangeCr),'）'];
% plotAllPhaseForT(data_reduced,rangeCr,dtCoh,'Title',T);
% 
% % CR2レンジ
% figure(2100)
% T = ['チップ間位相雑音除去後の位相','（Range：',num2str(rangeCr2),'）'];
% plotAllPhaseForT(data_reduced,rangeCr2,dtCoh,'Title',T);
% 
% % 人レンジ
% figure(3100)
% T = ['チップ間位相雑音除去後の位相','（Range：',num2str(rangeMain),'）'];
% plotAllPhaseForT(data_reduced,rangeMain,dtCoh,'Title',T);

%% キャリブレーション
data_calibrated = calibration(data_reduced,rangeCr);

% %% キャリブレーション後の位相表示
% % CRレンジ
% figure(1300)
% T = ['キャリブレーション後の位相','（Range：',num2str(rangeCr),'）'];
% plotAllPhaseForT(data_calibrated,rangeCr,dtCoh,'Title',T);
% 
% % CR2レンジ
% figure(2300)
% T = ['キャリブレーション後の位相','（Range：',num2str(rangeCr2),'）'];
% plotAllPhaseForT(data_calibrated,rangeCr2,dtCoh,'Title',T);
% 
% % 人レンジ
% figure(3300)
% T = ['キャリブレーション後の位相','（Range：',num2str(rangeMain),'）'];
% plotAllPhaseForT(data_calibrated,rangeMain,dtCoh,'Title',T);

%% 相関行列計算
X = squeeze(reshapeTxrx2No(data_calibrated(rangeMain,:,:,:))).';
if useOnly8Chs == 1
    X = extract8Chs(X.').';
end
Rxx = X * X' / size(X,2);
R = Rxx + max(diag(Rxx)).*eye(size(Rxx)).*sigma;
disp(['max(diag(Rxx))=',num2str(max(diag(Rxx)),'%e')]);
disp(['max(diag(Loading))=',num2str(max(diag(max(diag(Rxx)).*eye(size(Rxx)).*sigma)),'%e')]);
invR = inv(R);

%% サーチ
range_theta = 0:0.5:90;
range_phi = -180:0.5:180;
D = zeros(length(range_theta),length(range_phi));
if useOnly8Chs == 1
    r = extract8Chs(r.').';
end
count_theta = 0;
for theta_want_deg = range_theta
    count_theta = count_theta + 1;
    count_phi = 0;
    for phi_want_deg = range_phi
        count_phi = count_phi + 1;
        
        theta_want_rad = deg2rad(theta_want_deg);
        phi_want_rad = deg2rad(phi_want_deg);
        
%         V = steeringVector( theta_want_rad, phi_want_rad, lambda, r.' );

        a = steeringVector(theta_want_rad,phi_want_rad,lambda,r.');
        W = (invR*a)/(a'*invR*a);

        y = W' * X;
        
        D(count_theta,count_phi) = mean(abs(y).^2,2);
    end
end

D_norm = D./max(max(D));

%% 最大電力方向を求める
[row,col] = Index_Max(D);

theta_max_deg = range_theta(row);
phi_max_deg = range_phi(col);

disp(theta_max_deg);
disp(phi_max_deg);

%% プロット（長方形）
figure(300)
fig = gcf;
fig.OuterPosition = [100,300,3000,1200]; % [開始位置左，開始位置下，縦の長さ，横の長さ]
imagesc(range_phi,range_theta,pow2db(D_norm),[-20,0]);
colorbar;
axis xy
grid on
grid minor
% Title = ['（\theta_f,\phi_f）=（',num2str(theta_f_deg),'°',',',num2str(phi_f_deg),'°','）','，','（\theta_g,\phi_g）=（',num2str(theta_g_deg),'°',',',num2str(phi_g_deg),'°','）','，（sigma=',num2str(sigma),'）'];
% title(Title);
ax = gca;
ax.FontSize = 50;
ax.XLabel.String = 'Phi [degree]';
ax.XLabel.FontSize = 60;
ax.YLabel.String = 'Theta [dB]';
ax.YLabel.FontSize = 60;

%% プロット（円形）
figure(3003)
[ax, pax, cb] = spherical_pcolor(deg2rad(range_theta), deg2rad(range_phi), ...
    10*log10(D_norm), "az0", "left", "azdir", "cw",...
    "clim", [-30, 0], "fontsize", 40, "cblabel", "Level [dB]");
grid on
grid minor
Title1 = ['Range = ',num2str(rangeMain)];
% Title2 = ['（\theta_g,\phi_g）=（',num2str(theta_g_deg),'°',',',num2str(phi_g_deg),'°','）'];
Title3 = ['sigma=',num2str(sigma)];
ax.Title.String = {Title1;'';'';'';Title3};
ax.Title.FontSize = 45;
ax.Title.Position = [-160,0];
fig = gcf;
fig.Position = [100,100,2500,1000];

figName = ['Result_Search_sigma',num2str(sigma)];
print(['Figures/', figName, '.png'],'-dpng','-r0');


%%
toc