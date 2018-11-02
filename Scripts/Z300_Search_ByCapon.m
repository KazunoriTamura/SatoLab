%% �J�n
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints range
addpath('../../../GitHub/rd79codes/Functions')
addpath('../FunctionsOct/funcs')
addpath('../FunctionsOct')
set(0,'DefaultAxesFontSize',18)
tic

%% �p�����[�^
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
Tx(1,:) = Tx4;
Tx(2,:) = Tx3;
Tx(3,:) = Tx2;
Tx(4,:) = Tx1;

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

%% �f�[�^
data = data_coh;

%% �O����
n_time_max = size(data,2);
n_time_coh = 0:(n_time_max-1);
m_time_coh = n_time_coh * dtCoh;
lambda = c/f_c;

%% �f�q�z��v�Z
r = zeros(16,3);
for tx = 1:4
    for rx = 1:4
        No = txrx2no(tx,rx);
        r(No,:) = Rx(rx,:) - Tx(tx,:);
    end
end

%% �P���@�p�ʑ����]����
% data = conj(data);

%% CR�����W�ɂ��āC���f�[�^IQ���ʕ\��
figure(1)
T = ['Raw Data','�iRange�F',num2str(rangeCr),'�j'];
plotIqForT(data,rangeCr,'MaxValue',scalePlotIqForCr,'Title',T);

%% CR2�����W�ɂ��āC���f�[�^IQ���ʕ\��
figure(2)
T = ['Raw Data','�iRange�F',num2str(rangeCr2),'�j'];
plotIqForT(data,rangeCr2,'MaxValue',scalePlotIqForCr,'Title',T);

%% �l�����W�ɂ��āC���f�[�^IQ���ʕ\��
figure(3)
T = ['Raw Data','�iRange�F',num2str(range),'�j'];
plotIqForT(data,range,'MaxValue',scalePlotIq,'Title',T);

%% �ȉ~�␳�i�W���p�����[�^�j
% load('iqparams_1_180803_v4.mat');
% data_balanced = compensateImbalance(data,q,phi);
data_balanced = data;

%% CR�����W�ɂ��āC�ȉ~�␳��IQ���ʕ\���i�W���p�����[�^�j
figure(18)
T = ['�ȉ~�␳��','�iRange�F',num2str(rangeCr),'�j'];
plotIqForT(data_balanced,rangeCr,'MaxValue',scalePlotIqForCr,'Title',T);

%% CR2�����W�ɂ��āC�ȉ~�␳��IQ���ʕ\���i�W���p�����[�^�j
figure(28)
T = ['�ȉ~�␳��','�iRange�F',num2str(rangeCr2),'�j'];
plotIqForT(data_balanced,rangeCr2,'MaxValue',scalePlotIqForCr,'Title',T);

%% �l�����W�ɂ��āC�ȉ~�␳��IQ���ʕ\���i�W���p�����[�^�j
figure(28)
T = ['�ȉ~�␳��','�iRange�F',num2str(range),'�j'];
plotIqForT(data_balanced,range,'MaxValue',scalePlotIq,'Title',T);

%% CR�����W���o
dataCr = data_balanced(rangeCr,:,:,:);
dataCr2 = data_balanced(rangeCr2,:,:,:);

%% CR�����W�ɂ��āC�e��DC�̌v�Z
meanDc = getDcByMean(dataCr);
huDc = getDcByHu(dataCr);
modifiedHuDc = getDcByModifiedHu(dataCr);

%% CR2�����W�ɂ��āC�e��DC�̌v�Z
meanDc2 = getDcByMean(dataCr2);
huDc2 = getDcByHu(dataCr2);
modifiedHuDc2 = getDcByModifiedHu(dataCr2);

%% CR�����W�ɂ��āCIQ���ʕ\��
figure(100)
T = ['�ȉ~�␳���IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr),'�j'];
plotIqForT(dataCr,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

%% CR2�����W�ɂ��āCIQ���ʕ\��
figure(200)
T = ['�ȉ~�␳���IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr2),'�j'];
plotIqForT(dataCr2,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);

%% CR�����W�ɂ��āC�ȉ~�␳��̈ʑ��\��
% figure(110)
% T = ['�ȉ~�␳��̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
% plotAllPhaseForT(data,rangeCr,dtCoh,'Title',T);

%% CR2�����W�ɂ��āC�ȉ~�␳��̈ʑ��\��
% figure(210)
% T = ['�ȉ~�␳��̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
% plotAllPhaseForT(data,rangeCr2,dtCoh,'Title',T);

%% �C��DC������
dataCr_partwodc = reduceInputDc(dataCr,modifiedHuDc);
dataCr2_partwodc = reduceInputDc(dataCr2,modifiedHuDc2);

%% ����DC������
data_wodc = data_balanced - mean(data_balanced,2);
data_wodc(rangeCr,:,:,:) = dataCr_partwodc;
data_wodc(rangeCr2,:,:,:) = dataCr2_partwodc;

%% CR�����W�̒���
% data_wodc(rangeCr,:,:,:) = dataCr_partwodc;
% data_wodc(rangeCr2,:,:,:) = dataCr2_partwodc;

% %% �C��DC������̊e��DC�̌v�Z
% meanDc = getDcByMean(dataCR_partwodc);
% huDc = getDcByHu(dataCR_partwodc);
% modifiedHuDc = getDcByModifiedHu(dataCR_partwodc);
% 
% %% �C��DC������̊e��DC�̌v�Z
% meanDc2 = getDcByMean(dataCR2_partwodc);
% huDc2 = getDcByHu(dataCR2_partwodc);
% modifiedHuDc2 = getDcByModifiedHu(dataCR2_partwodc);
% 
% %% CR�����W�ɂ��āC�C��DC�������IQ���ʕ\��
% figure(105)
% T = ['�C��DC�������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCR),'�j'];
% plotIqForT(data,rangeCR,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);
% 
% %% CR2�����W�ɂ��āC�C��DC�������IQ���ʕ\��
% figure(205)
% T = ['�C��DC�������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCR2),'�j'];
% plotIqForT(data,rangeCR2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% CR�����W�ɂ��āC�C��DC������̈ʑ��\��
figure(115)
T = ['�`�b�v�Ԉʑ��������O�̈ʑ��i�C��DC������j','�iRange�F',num2str(rangeCr),'�j'];
plotAllPhaseForT(data_wodc,rangeCr,dtCoh,'Title',T);

%% CR2�����W�ɂ��āC�C��DC������̈ʑ��\��
figure(215)
T = ['�`�b�v�Ԉʑ��������O�̈ʑ��i�C��DC������j','�iRange�F',num2str(rangeCr2),'�j'];
plotAllPhaseForT(data_wodc,rangeCr2,dtCoh,'Title',T);

%% �l�����W�ɂ��āC�C��DC������̈ʑ��\��
figure(315)
T = ['�`�b�v�Ԉʑ��������O�̈ʑ��i����DC������j','�iRange�F',num2str(range),'�j'];
plotAllPhaseForT(data_wodc,range,dtCoh,'Title',T);

%% CR�����W����`�b�v�Ԉʑ��G���ϓ������߁C�����S�����W�ɂ��ď���
phaseNoise = getPhaseNoise(dataCr_partwodc);
data_reduced = reducePhaseNoise( data_wodc, phaseNoise );

%% CR�����W�ɂ��āC�`�b�v�ԎG��������̈ʑ��\��
figure(1000)
T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
plotAllPhaseForT(data_reduced,rangeCr,dtCoh,'Title',T);

%% CR2�����W�ɂ��āC�`�b�v�ԎG��������̈ʑ��\��
figure(2000)
T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
plotAllPhaseForT(data_reduced,rangeCr2,dtCoh,'Title',T);

%% Main�����W�ɂ��āC�`�b�v�Ԉʑ��G��������̈ʑ��\��
figure(3000)
T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(range),'�j'];
plotAllPhaseForT(data_reduced,range,dtCoh,'Title',T);

%% �f�[�^�č\��
data = data_reduced;

%% �l�����W�ɂ��āCIQ���ʕ\��
figure(301)
T = ['�S�������IQ����','�iRange�F',num2str(range),'�j'];
plotIqForT(data,range,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

% %% �m�C�Y���x���Z�o
% power = abs(data).^2;
% noiseLevel = estimateNoiseLevel(power);
% disp(noiseLevel)
% noiseLevel = noiseLevel./noiseLevel(1,1);
% one = ones(4,4);
% noiseLeInv = one./noiseLevel;
% 
% %% �m�C�Y���x���̐��K��
% data_constNoise = zeros(size(data));
% for iTx = 1:4
%     for iRx = 1:4
%         data_constNoise(:,:,iTx,iRx) = data(:,:,iTx,iRx) .* noiseLeInv(iTx,iRx);
%     end
% end
% data = data_constNoise;
% 

%% CR�̈ʑ��Z�o
angle_rCR = unwrap(angle(data(rangeCr,:,:,:)));
hosei_tmp = reshapeTxrx2No(angle_rCR);
hosei_tmp2 = (squeeze(hosei_tmp)).';
hosei_mean = mean(hosei_tmp2,2);
hosei = repmat(hosei_mean,1,size(angle_rCR,2));

%% �f�[�^
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
% plotAllPhaseForT(X_hosei_check,range,dtCoh,'Title','�r�[����␳��̔팱�҃����W�̈ʑ�');
% ax = gca;
% ax.TitleFontSizeMultiplier = 1.5;

%% ���֍s��v�Z
Rxx = X_hosei * X_hosei' / size(X_hosei,2);
R = Rxx + max(diag(Rxx)).*eye(size(Rxx)).*sigma;
disp(max(diag(Rxx)));
% R = Rxx;
invR = inv(R);

%% �E�F�C�g�x�N�g���Ƃ��̃��[�v
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

%% �v���b�g�i�����`�j
figure(300)
fig = gcf;
fig.OuterPosition = [100,300,3000,1200]; % [�J�n�ʒu���C�J�n�ʒu���C�c�̒����C���̒���]
imagesc(range_phi,range_theta,pow2db(D_norm),[-20,0]);
colorbar;
axis xy
grid on
grid minor
% Title = ['�i\theta_f,\phi_f�j=�i',num2str(theta_f_deg),'��',',',num2str(phi_f_deg),'��','�j','�C','�i\theta_g,\phi_g�j=�i',num2str(theta_g_deg),'��',',',num2str(phi_g_deg),'��','�j','�C�isigma=',num2str(sigma),'�j'];
% title(Title);
ax = gca;
ax.FontSize = 50;
ax.XLabel.String = 'Phi [degree]';
ax.XLabel.FontSize = 60;
ax.YLabel.String = 'Theta [dB]';
ax.YLabel.FontSize = 60;

%% �v���b�g�i�~�`�j
figure(3003)
[ax, pax, cb] = spherical_pcolor(deg2rad(range_theta), deg2rad(range_phi), ...
    10*log10(D_norm), "az0", "left", "azdir", "cw",...
    "clim", [-30, 0], "fontsize", 40, "cblabel", "Level [dB]");
grid on
grid minor
% Title1 = ['�i\theta_f,\phi_f�j=�i',num2str(theta_f_deg),'��',',',num2str(phi_f_deg),'��','�j'];
% Title2 = ['�i\theta_g,\phi_g�j=�i',num2str(theta_g_deg),'��',',',num2str(phi_g_deg),'��','�j'];
Title3 = ['sigma=',num2str(sigma)];
ax.Title.String = {'';'';'';'';Title3};
ax.Title.FontSize = 45;
ax.Title.Position = [-160,0];
fig = gcf;
fig.Position = [100,100,2500,1000];

% figName = ['Capon2d_(',num2str(theta_f_deg),'_',num2str(phi_f_deg),')_(',num2str(theta_g_deg),'_',num2str(phi_g_deg),')'];
% print(['Figures/', figName, '.png'],'-dpng','-r0');


%% �E�F�C�g�x�N�g���Ƃ��̃��[�v
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

%% �v���b�g(�f�V�x��)
% D_max = max(max(D));
% % [Row,Col] = Index_Max(abs(D));
% % Row_dash = Row - 1;
% % Col_dash = Col - 1;
% figure(10000)
% imagesc(-180:0.4:180,0:0.4:50,10*log10(D/D_max),[-20 0])
% % hold on
% % plot(Row_dash,Col_dash,'o')
% axis xy
% T = strcat('(Azimuth,Nadir)=','(',num2str(phi_want_deg),'��',',',num2str(theta_want_deg),'��',')');
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


%% [0:360]����[-180:180]�ւ̕ϊ�
% D_0center = zeros(size(D,1),size(D,2));
% D_0center(:,181:361) = D(:,1:181);
% D_0center(:,1:180) = D(:,181:360);

%% �v���b�g(�f�V�x��)
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
% T = strcat('(Azimuth,Nadir)=','(',num2str(phi_want_deg),'��',',',num2str(theta_want_deg),'��',')');
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

%% �I��
toc