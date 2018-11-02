%% 開始
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints range
addpath('../../../GitHub/rd79codes/Functions')
addpath('../FunctionsOct/funcs')
addpath('../FunctionsOct')
set(0,'DefaultAxesFontSize',18)
tic

%% パラメータ
rangeCr = 8;
rangeCr2 = 9;
range = 13;
dt = 0.000238;
dtCoh = dt * nCohPoints;
scalePlotIqForCr = 4e5;
scalePlotIq = 5e4;
c = 299792458;
f_c = 79e9;
sigma = 1e2;

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

%% データ
data = data_coh;

%% 前処理
n_time_max = size(data,2);
n_time_coh = 0:(n_time_max-1);
m_time_coh = n_time_coh * dtCoh;
lambda = c/f_c;

%% 素子配列計算
r = zeros(16,3);
for tx = 1:4
    for rx = 1:4
        No = txrx2no(tx,rx);
        r(No,:) = Rx(rx,:) - Tx(tx,:);
    end
end

%% １号機用位相反転操作
% data = conj(data);

%% CRレンジについて，生データIQ平面表示
figure(1)
T = ['Raw Data','（Range：',num2str(rangeCr),'）'];
plotIqForT(data,rangeCr,'MaxValue',scalePlotIqForCr,'Title',T);

%% CR2レンジについて，生データIQ平面表示
figure(2)
T = ['Raw Data','（Range：',num2str(rangeCr2),'）'];
plotIqForT(data,rangeCr2,'MaxValue',scalePlotIqForCr,'Title',T);

%% 人レンジについて，生データIQ平面表示
figure(3)
T = ['Raw Data','（Range：',num2str(range),'）'];
plotIqForT(data,range,'MaxValue',scalePlotIq,'Title',T);

%% 楕円補正（８月パラメータ）
% load('iqparams_1_180803_v4.mat');
% data_balanced = compensateImbalance(data,q,phi);
data_balanced = data;

%% CRレンジについて，楕円補正後IQ平面表示（８月パラメータ）
figure(18)
T = ['楕円補正後','（Range：',num2str(rangeCr),'）'];
plotIqForT(data_balanced,rangeCr,'MaxValue',scalePlotIqForCr,'Title',T);

%% CR2レンジについて，楕円補正後IQ平面表示（８月パラメータ）
figure(28)
T = ['楕円補正後','（Range：',num2str(rangeCr2),'）'];
plotIqForT(data_balanced,rangeCr2,'MaxValue',scalePlotIqForCr,'Title',T);

%% 人レンジについて，楕円補正後IQ平面表示（８月パラメータ）
figure(28)
T = ['楕円補正後','（Range：',num2str(range),'）'];
plotIqForT(data_balanced,range,'MaxValue',scalePlotIq,'Title',T);

%% CRレンジ抽出
dataCr = data_balanced(rangeCr,:,:,:);
dataCr2 = data_balanced(rangeCr2,:,:,:);

%% CRレンジについて，各種DCの計算
meanDc = getDcByMean(dataCr);
huDc = getDcByHu(dataCr);
modifiedHuDc = getDcByModifiedHu(dataCr);

%% CR2レンジについて，各種DCの計算
meanDc2 = getDcByMean(dataCr2);
huDc2 = getDcByHu(dataCr2);
modifiedHuDc2 = getDcByModifiedHu(dataCr2);

%% CRレンジについて，IQ平面表示
figure(100)
T = ['楕円補正後のIQ平面と各種DC','（Range：',num2str(rangeCr),'）'];
plotIqForT(dataCr,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

%% CR2レンジについて，IQ平面表示
figure(200)
T = ['楕円補正後のIQ平面と各種DC','（Range：',num2str(rangeCr2),'）'];
plotIqForT(dataCr2,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

%% CRレンジについて，楕円補正後の位相表示
% figure(110)
% T = ['楕円補正後の位相','（Range：',num2str(rangeCr),'）'];
% plotAllPhaseForT(data,rangeCr,dtCoh,'Title',T);

%% CR2レンジについて，楕円補正後の位相表示
% figure(210)
% T = ['楕円補正後の位相','（Range：',num2str(rangeCr2),'）'];
% plotAllPhaseForT(data,rangeCr2,dtCoh,'Title',T);

%% 修正DCを引く
dataCr_partwodc = reduceInputDc(dataCr,modifiedHuDc);
dataCr2_partwodc = reduceInputDc(dataCr2,modifiedHuDc2);

%% 平均DCを引く
data_wodc = data_balanced - mean(data_balanced,2);
data_wodc(rangeCr,:,:,:) = dataCr_partwodc;
data_wodc(rangeCr2,:,:,:) = dataCr2_partwodc;

%% CRレンジの調整
% data_wodc(rangeCr,:,:,:) = dataCr_partwodc;
% data_wodc(rangeCr2,:,:,:) = dataCr2_partwodc;

% %% 修正DC除去後の各種DCの計算
% meanDc = getDcByMean(dataCR_partwodc);
% huDc = getDcByHu(dataCR_partwodc);
% modifiedHuDc = getDcByModifiedHu(dataCR_partwodc);
% 
% %% 修正DC除去後の各種DCの計算
% meanDc2 = getDcByMean(dataCR2_partwodc);
% huDc2 = getDcByHu(dataCR2_partwodc);
% modifiedHuDc2 = getDcByModifiedHu(dataCR2_partwodc);
% 
% %% CRレンジについて，修正DC除去後のIQ平面表示
% figure(105)
% T = ['修正DC除去後のIQ平面と各種DC','（Range：',num2str(rangeCR),'）'];
% plotIqForT(data,rangeCR,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);
% 
% %% CR2レンジについて，修正DC除去後のIQ平面表示
% figure(205)
% T = ['修正DC除去後のIQ平面と各種DC','（Range：',num2str(rangeCR2),'）'];
% plotIqForT(data,rangeCR2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% CRレンジについて，修正DC除去後の位相表示
figure(115)
T = ['チップ間位相差除去前の位相（修正DC除去後）','（Range：',num2str(rangeCr),'）'];
plotAllPhaseForT(data_wodc,rangeCr,dtCoh,'Title',T);

%% CR2レンジについて，修正DC除去後の位相表示
figure(215)
T = ['チップ間位相差除去前の位相（修正DC除去後）','（Range：',num2str(rangeCr2),'）'];
plotAllPhaseForT(data_wodc,rangeCr2,dtCoh,'Title',T);

%% 人レンジについて，修正DC除去後の位相表示
figure(315)
T = ['チップ間位相差除去前の位相（平均DC除去後）','（Range：',num2str(range),'）'];
plotAllPhaseForT(data_wodc,range,dtCoh,'Title',T);

%% CRレンジからチップ間位相雑音変動を求め，それを全レンジについて除去
phaseNoise = getPhaseNoise(dataCr_partwodc);
data_reduced = reducePhaseNoise( data_wodc, phaseNoise );

%% CRレンジについて，チップ間雑音除去後の位相表示
figure(1000)
T = ['チップ間位相雑音除去後の位相','（Range：',num2str(rangeCr),'）'];
plotAllPhaseForT(data_reduced,rangeCr,dtCoh,'Title',T);

%% CR2レンジについて，チップ間雑音除去後の位相表示
figure(2000)
T = ['チップ間位相雑音除去後の位相','（Range：',num2str(rangeCr2),'）'];
plotAllPhaseForT(data_reduced,rangeCr2,dtCoh,'Title',T);

%% Mainレンジについて，チップ間位相雑音除去後の位相表示
figure(3000)
T = ['チップ間位相雑音除去後の位相','（Range：',num2str(range),'）'];
plotAllPhaseForT(data_reduced,range,dtCoh,'Title',T);

%% データ再構成
data = data_reduced;

%% 人レンジについて，IQ平面表示
figure(301)
T = ['全処理後のIQ平面','（Range：',num2str(range),'）'];
plotIqForT(data,range,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

% %% ノイズレベル算出
% power = abs(data).^2;
% noiseLevel = estimateNoiseLevel(power);
% disp(noiseLevel)
% noiseLevel = noiseLevel./noiseLevel(1,1);
% one = ones(4,4);
% noiseLeInv = one./noiseLevel;
% 
% %% ノイズレベルの正規化
% data_constNoise = zeros(size(data));
% for iTx = 1:4
%     for iRx = 1:4
%         data_constNoise(:,:,iTx,iRx) = data(:,:,iTx,iRx) .* noiseLeInv(iTx,iRx);
%     end
% end
% data = data_constNoise;
% 

%% CRの位相算出
angle_rCR = unwrap(angle(data(rangeCr,:,:,:)));
hosei_tmp = reshapeTxrx2No(angle_rCR);
hosei_tmp2 = (squeeze(hosei_tmp)).';
hosei_mean = mean(hosei_tmp2,2);
hosei = repmat(hosei_mean,1,size(angle_rCR,2));

%% データ
X_tmp = reshapeTxrx2No(data(range,:,:,:));
X = (squeeze(X_tmp)).';
X_hosei = zeros(16,size(X,2));
phi = -1i * hosei;
for k = 1:16
    X_hosei(k,:) = X(k,:) .* exp(phi(k,:));
end

% X_hosei_check_reshaped = zeros(1,size(X_hosei,2),size(X_hosei,1));
% X_hosei_check_reshaped(1,:,:) = X_hosei.';
% X_hosei_check = reshapeNo2Txrx(X_hosei_check_reshaped);
% 
% figure(800)
% set(0,'DefaultAxesFontSize',40);
% plotAllPhaseForT(X_hosei_check,range,dtCoh,'Title','ビーム基準補正後の被験者レンジの位相');
% ax = gca;
% ax.TitleFontSizeMultiplier = 1.5;

%% 相関行列計算
Rxx = X_hosei * X_hosei' / size(X_hosei,2);
R = Rxx + max(diag(Rxx)).*eye(size(Rxx)).*sigma;
disp(max(diag(Rxx)));
% R = Rxx;
invR = inv(R);

%% ウェイトベクトルとそのループ
range_theta = 0:0.5:90;
range_phi = -180:0.5:180;
D = zeros(length(range_theta),length(range_phi));

% a = steeringVector(theta_f,phi_f,lambda,r.');
% W = (invR*a)/(a'*invR*a);

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
% Title1 = ['（\theta_f,\phi_f）=（',num2str(theta_f_deg),'°',',',num2str(phi_f_deg),'°','）'];
% Title2 = ['（\theta_g,\phi_g）=（',num2str(theta_g_deg),'°',',',num2str(phi_g_deg),'°','）'];
Title3 = ['sigma=',num2str(sigma)];
ax.Title.String = {'';'';'';'';Title3};
ax.Title.FontSize = 45;
ax.Title.Position = [-160,0];
fig = gcf;
fig.Position = [100,100,2500,1000];

% figName = ['Capon2d_(',num2str(theta_f_deg),'_',num2str(phi_f_deg),')_(',num2str(theta_g_deg),'_',num2str(phi_g_deg),')'];
% print(['Figures/', figName, '.png'],'-dpng','-r0');


%% ウェイトベクトルとそのループ
% range_theta = 0:0.4:50;
% range_phi = -180:0.4:180;
% D = zeros(length(range_theta),length(range_phi));
% 
% count_theta = 0;
% for theta_want_deg = range_theta
%     count_theta = count_theta + 1;
%     count_phi = 0;
%     for phi_want_deg = range_phi
%         count_phi = count_phi + 1;
%         
%         theta_want_rad = deg2rad(theta_want_deg);
%         phi_want_rad = deg2rad(phi_want_deg);
%         
%         a = steeringVector( theta_want_rad, phi_want_rad, lambda, r.' );
%         
%         W = (invR*a)/(a'*invR*a);
%         
%         y = W' * X_hosei;
%         
%         D(count_theta,count_phi) = mean(abs(y).^2);
%     end
% end

%% プロット(デシベル)
% D_max = max(max(D));
% % [Row,Col] = Index_Max(abs(D));
% % Row_dash = Row - 1;
% % Col_dash = Col - 1;
% figure(10000)
% imagesc(-180:0.4:180,0:0.4:50,10*log10(D/D_max),[-20 0])
% % hold on
% % plot(Row_dash,Col_dash,'o')
% axis xy
% T = strcat('(Azimuth,Nadir)=','(',num2str(phi_want_deg),'°',',',num2str(theta_want_deg),'°',')');
% title(T)
% c = colorbar;
% c.Label.String = 'Power [dB]';
% c.Label.FontSize = 65;
% ax = gca;
% ax.FontSize = 45;
% ax.XLabel.String = 'Azimuth angle [deg]';
% ax.XLabel.FontSize = 65;
% ax.YLabel.String = 'Nadir angle [deg]';
% ax.YLabel.FontSize = 65;


%% [0:360]から[-180:180]への変換
% D_0center = zeros(size(D,1),size(D,2));
% D_0center(:,181:361) = D(:,1:181);
% D_0center(:,1:180) = D(:,181:360);

%% プロット(デシベル)
% D_power = abs(D_0center).^2;
% D_power_max = max(max(D_power));
% [Row,Col] = Index_Max(abs(D));
% Row_dash = Row - 1;
% Col_dash = Col - 1;
% figure(10000)
% imagesc(-180:180,0:90,10*log10(D_power/D_power_max),[-20 0])
% % hold on
% % plot(Row_dash,Col_dash,'o')
% axis xy
% T = strcat('(Azimuth,Nadir)=','(',num2str(phi_want_deg),'°',',',num2str(theta_want_deg),'°',')');
% title(T)
% c = colorbar;
% c.Label.String = 'Power [dB]';
% c.Label.FontSize = 65;
% ax = gca;
% ax.FontSize = 45;
% ax.XLabel.String = 'Azimuth angle [deg]';
% ax.XLabel.FontSize = 65;
% ax.YLabel.String = 'Nadir angle [deg]';
% ax.YLabel.FontSize = 65;

%%
% figure(30000)
% clf
% [ax, pax, cb] = spherical_pcolor(deg2rad(0:0.1:90), deg2rad(-180:0.1:180), ...
%     10*log10(D / max(D(:))), "az0", "left", "azdir", "cw",...
%     "clim", [-20, 0], "fontsize", 40, "cblabel", "Level [dB]");
% fig = gcf;
% fig.OuterPosition = [100,100,1600,900];
% 
% figName = ['r',num2str(range),'_Circle_CP_NNL'];
% print(['Figures/', figName, '.png'],'-dpng','-r0');
% % Example: add a beam direction in the polar coordinates on the surface.
% % hold on
% % polarplot(pax, deg2rad(45), 40, "ro")
% % polarplot(pax,deg2rad(45), 40, "ro")
% % ax = gca;
% % ax.FontSize = 45;
% % ax.XLabel.String = 'Azimuth angle [deg]';
% % ax.XLabel.FontSize = 65;
% % ax.YLabel.String = 'Nadir angle [deg]';
% % ax.YLabel.FontSize = 65;
% 
% disp(Row_dash)
% disp(Col_dash)

%% 終了
toc