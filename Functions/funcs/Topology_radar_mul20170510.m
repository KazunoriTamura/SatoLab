function [ outputs ] = Topology_radar_mul20170510( sradar10, params)
%Topology_radar トポロジー法
% params はパラメータを格納した構造体。 以下の12個が必要。
% params.ncoh, params.Tr, params.cTr, params.HBR1, params.HBR2 
% params.Thcor, params.Thweight, params.Tc, params.NforF, params.CN 
% params.filt_length_low, params.filt_length_high
% radarはデータのアンラップ済み位相を格納したベクトル。(1xdatalength)

%% ハイパス
% filt_length_high = (params.filt_length_high);
% % % % % % filt_high = [0; hanning(int8(params.filt_length_high)); 0];
% % % % % % filt_high = filt_high/sum(filt_high);
% % % % % % % filt_length_high = params.filt_length_high+2;
% filt_length_high=length(filt.filt_high);
% radar_filt_high = filter(filt.filt_high, 1, radar);
% 
% dradar=radar_filt_high(round(filt_length_high/2):end);
% radar_sel=radar(1:length(dradar));
% radar_wot=radar_sel-dradar;
% % mradar=max(radar);
% % % % % % % figure(1)
% % % % % % % plot((1:length(radar))/params.cTr,radar/mradar,'-x',(1:length(dradar))/params.cTr,dradar/mradar,'-o')
% % % % % % % xlabel('sample number')
% % % % % % % legend('処理前', 'トレンド')
% % % % % % % ylim([-1 1])
% % % % % % % drawnow
% %% ローパス
% % filt_length_low = (params.filt_length_low);
% % filt_low = hanning(filt_length_low);
% % filt_low = [0; hanning(filt_length_low); 0];
% % filt_low = filt_low/sum(filt_low);
% filt_length_low = length(filt.filt_low);
% outputs.t_delay_filt=round(filt_length_low/2)*params.cTr;
% % filt_low = filt_low/sum(filt_low);
% 
% sradar10 = filter(filt.filt_low, 1, radar_wot);
% figure(10)
% plot((1:length(radar_wot))/params.cTr,radar_wot,(1:length(sradar10))/params.cTr, sradar10)
% xlabel('sample number')
% legend('トレンド除去後', 'トレンド除去平滑化後')
% 
% % % % % ylim([-1 1])
% % % % % drawnow
M=size(sradar10,2);
times=(0:M-1)*params.cTr;
%% Differentiating radar signal
dsradar=zeros(size(sradar10));
dsradar(1:end-1) = diff(sradar10)/params.cTr;

%% Detecting features
sradar = sradar10;
K=1;
K0=1;
K1=1;
K2=1;
K3=1;
%% For 1st diff
dif1_p=abs(sradar(2:M-1))>0.0;
dif2_1p=(sradar(2:M-1)-sradar(1:M-2))>0; 
dif2_1m=(sradar(2:M-1)-sradar(1:M-2))<0; 
dif2_2p=(sradar(3:M)-sradar(2:M-1))>0; 
dif2_2m=(sradar(3:M)-sradar(2:M-1))<0; 

forc1p=ones(1,M);
forc1m=ones(1,M)*-1;
forci1p=ones(1,M)*1i;
forci1m=ones(1,M)*-1i;
forci12p=ones(1,M)*1i/2;
forci12m=ones(1,M)*-1i/2;
%% PK
ipeak0=find(dif1_p.*dif2_1p.*dif2_2m)+1;
tpeak0=times(ipeak0);
apeak0=sradar(ipeak0);
cpeak0=forc1m(ipeak0);
%% VL
ipeak10=find(dif1_p.*dif2_1m.*dif2_2p)+1;
tpeak10=times(ipeak10);
apeak10=sradar(ipeak10);
cpeak10=forc1p(ipeak10);
%% For 2nd diff
ds1_p=abs(dsradar(2:M-1))>0.0;
dsdif_1p= (dsradar(2:M-1)-dsradar(1:M-2))>0;
dsdif_1m= (dsradar(2:M-1)-dsradar(1:M-2))<0;
dsdif_2p=(dsradar(3:M)-dsradar(2:M-1))>0;
dsdif_2m=(dsradar(3:M)-dsradar(2:M-1))<0;
%% Rising/Inflection
ipeak20=find(ds1_p.*dsdif_1p.*dsdif_2m)+1;
tpeak20=(times(ipeak20)+times(ipeak20+1))/2;
apeak20=(sradar(ipeak20)+sradar(ipeak20+1))/2;
pkk=dsradar(2:M-1)>0;
ipkk=find(ds1_p.*dsdif_1p.*dsdif_2m.*pkk)+1;
forcp20=forci12m;
forcp20(ipkk)=forci1p(ipkk);
cpeak20=forcp20(ipeak20);
ipeak20=ipeak20+0.5;
%% Falling/Inflection
ipeak30=find(ds1_p.*dsdif_1m.*dsdif_2p)+1;
tpeak30=(times(ipeak30)+times(ipeak30+1))/2;
apeak30=(sradar(ipeak30)+sradar(ipeak30+1))/2;
pkk=dsradar(2:M-1)<0;
ipkk=find(ds1_p.*dsdif_1m.*dsdif_2p.*pkk)+1;
forcp30=forci12p;
forcp30(ipkk)=forci1m(ipkk);
cpeak30=forcp30(ipeak30);
ipeak30=ipeak30+0.5;
tpeak000=[tpeak0 tpeak10 tpeak20 tpeak30];
ipeak000=[ipeak0 ipeak10 ipeak20 ipeak30];
apeak000=[apeak0 apeak10 apeak20 apeak30];
cpeak000=[cpeak0 cpeak10 cpeak20 cpeak30];
alls=sortrows([tpeak000' ipeak000' transpose(apeak000) transpose(cpeak000)]);
tpeak=alls(:,1)';
ipeak=alls(:,2)';
apeak=transpose(alls(:,3));
cpeak=transpose(alls(:,4));


%%
PEAKN = size(tpeak,2); %Number of feature points
PEAKN0 = size(tpeak0,2); %Number of PK
CORN = round(PEAKN0/times(end));
acor=zeros(PEAKN,CORN);
ibi = zeros(1,M);
ibim = zeros(1,M);
ibi0 = zeros(1,size(tpeak,2));

%% Feature sequence 1400 is determined by 1.8s/1.285ms
% 特徴点を中心に1400点の信号を抽出
% 特徴点が700点目以前にある場合は、それ以前は0としている。
gnas = zeros(PEAKN,params.NforF);

for I=1:PEAKN
    for J = -params.NforF/2:params.NforF/2-1
        if round(ipeak(I))+J >= 1 && round(ipeak(I))+J <= M
            gnas(I,J+params.NforF/2+1) = sradar(round(ipeak(I))+J);
        end
    end
end

%% checking topology

cv = complex(zeros(PEAKN,params.CN)); % complex vector
pcv = 0;
ncv = 0;
for I=1:PEAKN % the number of peaks
    cv1 = complex(zeros(1,params.CN));
    for K=-6:6 % uses 12 feature points around
        if I+K >= 1 && I+K <= PEAKN % the feature point must be within a valid range 例えば、一つ目の特徴点の場合、K>=0以降しか計算しない
            if abs(ipeak(I+K)-ipeak(I)) <= (params.CN-1)/2 % vector is 251 dimensional
                if cpeak(I+K) == 1i/2 % change the sign of imagenary part % Fig. 8 to Fig. 9
                    localcpeak = -1i/2;
                else 
                    if cpeak(I+K)== -1i/2 % change the sign of imagenary part % Fig. 8 to Fig. 9
                        localcpeak = 1i/2;
                    else
                        localcpeak = cpeak(I+K);
                    end
                end
                cv1(round(ipeak(I+K)-ipeak(I))+(params.CN-1)/2+1) = localcpeak;
            end
        end
    end
    cv(I,:) = cv1; % store the vector to matrix cv
    inds=find(cv1~=0);
    if(~isempty(inds))
        cv_sel=[0, cv1(inds), 0];
        inds_sel=[min(inds)-1, inds, max(inds)+1];
        cv2_temp=interp1(inds_sel,cv_sel,min(inds)-1:1:max(inds)+1);
        cv(I,min(inds):max(inds))=cv2_temp(2:end-1);
    end    
%     %  interpolation of phase
%     for L=2:params.CN-1
%         if cv1(L) == 0
%             L0=1;
%             while cv1(L-L0)==0 && L-L0 > 1
%                 L0 = L0 + 1;
%             end
%             PL = L-L0;
%             if cv1(PL)~=0
%                 pcv = cv1(PL);
%             else
%                 pcv = 0;
%             end
%             
%             L0=1;
%             while cv1(L+L0)==0 && L+L0 < params.CN
%                 L0 = L0 + 1;
%             end
%             NL = L+L0;
%             if cv1(NL)~=0
%                 ncv = cv1(NL);
%             else
%                 ncv = 0;
%             end
%             
%             if pcv ~= 0 && ncv ~= 0
%                 cv(I,L) = ((pcv)*(NL-L) + (ncv)*(L-PL))/(NL-PL);
%             end
%         end
%     end
end
cv_temp=cv;
cv = exp(1i*angle(cv));

%%
weightpeak = complex(zeros(PEAKN,PEAKN));
for I=1:PEAKN
    for J=1:PEAKN
        if cpeak(I) == cpeak(J)
            weightpeak(I,J) = cv(I,:)*cv(J,:)';
        end
    end
end

%%
SEARCHN = round(PEAKN/times(end))*2; %Calcurate length for correlation 2秒の中に何個特徴点があるか？
correlation = zeros(PEAKN,SEARCHN);
maxcor = zeros(1,PEAKN);
ibi3 = complex(zeros(1,PEAKN));
chosenpair = zeros(1,PEAKN);
avemaxcor = 1;
for I=1:PEAKN-SEARCHN
    for J=1:SEARCHN
        vec1 = gnas(I,:);
        vec2 = gnas(I+J,:);
        vec1 = vec1 - real(mean(vec1));
        vec2 = vec2 - real(mean(vec2));
        correlation(I,J) = vec1*vec2'/sqrt(vec1*vec1')/sqrt(vec2*vec2');%Correlation between 特徴点Iに対する信号,特徴点I+Jに対する信号
        if tpeak(I+J)-tpeak(I) < params.HBR1 || tpeak(I+J)-tpeak(I) > params.HBR2 || cpeak(I)~=cpeak(I+J) %Threshold for correlation
            correlation(I,J) = 0;
        end
    end
    
    % weighting
    mcv = 1:SEARCHN;
    if I<11
        weightcor = zeros(1,SEARCHN)+1;
    else
        bias0 = 0;
        weightsum = 0;
        for J=I-1:-1:I-10
            if maxcor(J)>0
                bias0 = bias0 + maxcor(J)*exp((J-I)/3); %忘却係数 3
                weightsum = weightsum + exp((J-1)/3);
            end
        end
        if weightsum ~= 0
            bias0 = bias0 / weightsum;
            weightcor = exp(-(mcv-bias0).^2/2/20^2); % 120 大きくなれば過去の影響を受けにくくなる。
        else
            weightcor = zeros(1,SEARCHN)+1;
        end
    end
    [maxval,maxcor(I)] = max(correlation(I,:).*weightcor);
    if max(correlation(I,:)) == 0
        maxcor(I) = 0;
    end
    if maxval < params.Thcor
        maxcor(I) = 0;
    end
    
    ibi3(I) = tpeak(maxcor(I)+I)-tpeak(I);
    chosenpair(I) = maxcor(I)+I; %Select the Feature point that has highest correlation
end

%%
localweight = zeros(1,PEAKN);
for I=1:PEAKN
    if chosenpair(I) ~= 0
        localweight(I) = max(0,real(weightpeak(I,chosenpair(I))));
    else
        localweight(I) = 0;
    end
end
localweight = localweight/max(abs(localweight));

ibi11 = complex(zeros(1,PEAKN));
for I=1:PEAKN
    if localweight(I) > params.Thweight
        ibi11(I) = ibi3(I);
    else
        ibi11(I) = 0;
    end
end
outputs.ibi11=ibi11;
outputs.ipeak=ipeak;
end

