%% �J�n
close all
clearvars -except datas_all rdrdate ECG data_coh nCohPoints 
addpath('../../../GitHub/rd79codes/Functions')
addpath('../Functions/funcs');
addpath('../Functions')
tic

%% �p�����[�^
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

%% �f�q�z��v�Z
r = zeros(16,3);
for tx = 1:4
    for rx = 1:4
        No = txrx2no(tx,rx);
        r(No,:) = Rx(rx,:) - Tx(tx,:);
    end
end

%% �f�[�^
data = data_coh(:,1:size(data_coh,2)/5,:,:);

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
% % CR�����W
% figure(1)
% T = ['Raw Data','�iRange�F',num2str(rangeCr),'�j'];
% plotIqForT(data,rangeCr,'MaxValue',scalePlotIqForCr,'Title',T);
% 
% % CR2�����W
% figure(2)
% T = ['Raw Data','�iRange�F',num2str(rangeCr2),'�j'];
% plotIqForT(data,rangeCr2,'MaxValue',scalePlotIqForCr,'Title',T);
% 
% % �l�����W
% figure(3)
% T = ['Raw Data','�iRange�F',num2str(rangeMain),'�j'];
% plotIqForT(data,rangeMain,'MaxValue',scalePlotIq,'Title',T);

%% �ȉ~�␳
if radarId == '1'
    load('iqparams_1_180803_v4.mat');
elseif radarId == '2'
    load('iqparams_2_180905.mat');
else
    error('�G���[�I');
end
data_balanced = compensateImbalance(data,q,phi);
data_balanced = data;

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
% % CR�����W
% figure(100)
% T = ['�ȉ~�␳���IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr),'�j'];
% plotIqForT(dataCr,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % CR2�����W
% figure(200)
% T = ['�ȉ~�␳���IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr2),'�j'];
% plotIqForT(dataCr2,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % �l�����W
% figure(300)
% T = ['�ȉ~�␳���IQ���ʂƊe��DC','�iRange�F',num2str(rangeMain),'�j'];
% plotIqForT(dataMain,rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);
% 
% %% �ȉ~�␳��̈ʑ��\��
% % CR�����W
% figure(110)
% T = ['�ȉ~�␳��̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
% plotAllPhaseForT(dataCr,rangeCr,dtCoh,'Title',T);
% 
% % CR2�����W
% figure(210)
% T = ['�ȉ~�␳��̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
% plotAllPhaseForT(dataCr2,rangeCr2,dtCoh,'Title',T);
% 
% % �l�����W
% figure(310)
% T = ['�ȉ~�␳��̈ʑ�','�iRange�F',num2str(rangeMain),'�j'];
% plotAllPhaseForT(dataMain,rangeMain,dtCoh,'Title',T);


%% ������



%% �C��DC������
dataCr_partwodc = reduceInputDc(dataCr,modifiedHuDc);
dataCr2_partwodc = reduceInputDc(dataCr2,modifiedHuDc2);

%% ����DC������
data_wodc = data_balanced - mean(data_balanced,2);
% data_wodc = data_balanced;
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
% % CR�����W
% figure(120)
% T = ['�C��DC�������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr),'�j'];
% plotIqForT(dataCr_partwodc,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % CR2�����W
% figure(220)
% T = ['�C��DC�������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr2),'�j'];
% plotIqForT(dataCr2_partwodc,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % �l�����W
% figure(320)
% T = ['�C��DC�������IQ���ʂƊe��DC','�iRange�F',num2str(rangeMain),'�j'];
% plotIqForT(data_wodc(rangeMain,:,:,:),rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);
 
%% �C��DC������̈ʑ��\��
% % CR�����W
% figure(130)
% T = ['�C��DC������̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
% plotAllPhaseForT(dataCr_partwodc,rangeCr,dtCoh,'Title',T);
% 
% % CR2�����W
% figure(230)
% T = ['�C��DC������̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
% plotAllPhaseForT(dataCr2_partwodc,rangeCr2,dtCoh,'Title',T);
% 
% % �l�����W
% figure(230)
% T = ['�C��DC������̈ʑ�','�iRange�F',num2str(rangeMain),'�j'];
% plotAllPhaseForT(data_wodc(rangeMain,:,:,:),rangeMain,dtCoh,'Title',T);

%% CR�����W����`�b�v�Ԉʑ��G���ϓ������߁C�����S�����W�ɂ��ď���
phaseNoise = getPhaseNoise(dataCr_partwodc);
data_reduced = reducePhaseNoise( data_wodc, phaseNoise );

%% �`�b�v�ԎG���������IQ���ʕ\��
% % CR�����W
% figure(1000)
% T = ['�`�b�v�Ԉʑ��G���������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr),'�j'];
% plotIqForT(data_reduced,rangeCr,'meanDc',meanDc,'huDc',huDc,'modifiedHuDc',modifiedHuDc,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % CR2�����W
% figure(2000)
% T = ['�`�b�v�Ԉʑ��G���������IQ���ʂƊe��DC','�iRange�F',num2str(rangeCr2),'�j'];
% plotIqForT(data_reduced,rangeCr2,'meanDc',meanDc2,'huDc',huDc2,'modifiedHuDc',modifiedHuDc2,'MaxValue',scalePlotIqForCr,'Title',T,'MarkerSize',5);
% 
% % �l�����W
% figure(3000)
% T = ['�`�b�v�Ԉʑ��G���������IQ���ʂƊe��DC','�iRange�F',num2str(rangeMain),'�j'];
% plotIqForT(data_reduced,rangeMain,'meanDc',meanDcMain,'MaxValue',scalePlotIq,'Title',T,'MarkerSize',5);

%% �`�b�v�ԎG��������̈ʑ��\��
% % CR�����W
% figure(1100)
% T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
% plotAllPhaseForT(data_reduced,rangeCr,dtCoh,'Title',T);
% 
% % CR2�����W
% figure(2100)
% T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
% plotAllPhaseForT(data_reduced,rangeCr2,dtCoh,'Title',T);
% 
% % �l�����W
% figure(3100)
% T = ['�`�b�v�Ԉʑ��G��������̈ʑ�','�iRange�F',num2str(rangeMain),'�j'];
% plotAllPhaseForT(data_reduced,rangeMain,dtCoh,'Title',T);

%% �L�����u���[�V����
data_calibrated = calibration(data_reduced,rangeCr);

% %% �L�����u���[�V������̈ʑ��\��
% % CR�����W
% figure(1300)
% T = ['�L�����u���[�V������̈ʑ�','�iRange�F',num2str(rangeCr),'�j'];
% plotAllPhaseForT(data_calibrated,rangeCr,dtCoh,'Title',T);
% 
% % CR2�����W
% figure(2300)
% T = ['�L�����u���[�V������̈ʑ�','�iRange�F',num2str(rangeCr2),'�j'];
% plotAllPhaseForT(data_calibrated,rangeCr2,dtCoh,'Title',T);
% 
% % �l�����W
% figure(3300)
% T = ['�L�����u���[�V������̈ʑ�','�iRange�F',num2str(rangeMain),'�j'];
% plotAllPhaseForT(data_calibrated,rangeMain,dtCoh,'Title',T);

%% ���֍s��v�Z
X = squeeze(reshapeTxrx2No(data_calibrated(rangeMain,:,:,:))).';
if useOnly8Chs == 1
    X = extract8Chs(X.').';
end
Rxx = X * X' / size(X,2);
R = Rxx + max(diag(Rxx)).*eye(size(Rxx)).*sigma;
disp(['max(diag(Rxx))=',num2str(max(diag(Rxx)),'%e')]);
disp(['max(diag(Loading))=',num2str(max(diag(max(diag(Rxx)).*eye(size(Rxx)).*sigma)),'%e')]);
invR = inv(R);

%% �T�[�`
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

%% �ő�d�͕��������߂�
[row,col] = Index_Max(D);

theta_max_deg = range_theta(row);
phi_max_deg = range_phi(col);

disp(theta_max_deg);
disp(phi_max_deg);

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
Title1 = ['Range = ',num2str(rangeMain)];
% Title2 = ['�i\theta_g,\phi_g�j=�i',num2str(theta_g_deg),'��',',',num2str(phi_g_deg),'��','�j'];
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