%
% 20180513改訂 by 田村
%
%% 開始
clear
fld = dir('2018*');
dirname = strcat(fld.name,'/'); % ソースディレクトリ

%% ファイル名の取得
k = dir([dirname, 'rd79_*.bin']);
n_file = size(k,1);

%% ヘッダ読込
filename = k(1).name;
fd = fopen([dirname, filename]);
ident_str  = fread(fd,[1,8],'*char');
if(strcmp(ident_str,'MEIRDB02')==0)
    error('error: MEIRDBIN_IDENT_MISMATCH');
end
numofbyte  = fread(fd,1,'bit16');
numoffield = fread(fd,1,'bit16');
flag = fread(fd,1,'bit16');
sysid = int16(flag)/256;
modeid = rem(flag, 256);

if(sysid ~= 2)
    error('error: NOT 79GHz');
end

% 先頭ファイルの開始時刻
rdrdate = datetime(fread(fd,1,'bit16'), fread(fd,1,'bit16'), fread(fd,1,'bit16'), fread(fd,1,'bit16'),fread(fd,1,'bit16'), fread(fd,1,'bit16'), fread(fd,1,'bit16'));

nrange = fread(fd,1,'bit16');
max_ch = fread(fd,1,'bit16');
samples = fread(fd,1,'bit16');
errors = fread(fd,1,'bit16');
bitn_str = ['bit' num2str(numofbyte*8)];

fseek(fd,57,'bof');
radar_id  = fread(fd,[1,7],'*char');

% 出力ファイルの名前
if length(dirname) == 16
    filename_out = ['rd79_', datestr(rdrdate,'yymmdd_HHMMSS'), '_', radar_id, '.mat'];
else
    disp('フォルダ名が違います')
end

fclose(fd);

%% 格納用配列
datas_all = zeros(nrange, samples*n_file-1, sqrt(max_ch), sqrt(max_ch));

tic
for ik = 1:n_file
    %% ファイル読み込み(RAWDATA)
    filename = k(ik).name;
    fd = fopen([dirname, filename]);
    
    ident_str  = fread(fd,[1,8],'*char');
    if(strcmp(ident_str,'MEIRDB02')==0)
        error('error: MEIRDBIN_IDENT_MISMATCH');
    end
    numofbyte  = fread(fd,1,'bit16');
    numoffield = fread(fd,1,'bit16');
    flag = fread(fd,1,'bit16');
    sysid = int16(flag)/256;
    modeid = rem(flag, 256);
    
    if(sysid ~= 2)
        error('error: NOT 79GHz');
    end
    
    date = datetime(fread(fd,1,'bit16'), fread(fd,1,'bit16'), fread(fd,1,'bit16'), fread(fd,1,'bit16'),fread(fd,1,'bit16'), fread(fd,1,'bit16'), fread(fd,1,'bit16'));
    
    nrange = fread(fd,1,'bit16');
    max_ch = fread(fd,1,'bit16');
    samples = fread(fd,1,'bit16');
    errors(ik) = fread(fd,1,'bit16');
    
    fseek(fd,64,'bof');
    Rawdatas=fread(fd,[numoffield,Inf],bitn_str)';
    
    %% IQ-Range分割　I-Qに分かれて保存されているものを加算して複素に変換。その後、配列の形を変えて、Range-Time-Tx,Rxの組み合わせに変換します。Tx-Rxは別れていないのでRange-Time-4の大きさになります。
    dataiq=zeros(size(Rawdatas,1),size(Rawdatas,2)/2);
    for id=1:size(Rawdatas,2)/2
        dataiq(:,id)=Rawdatas(:,id*2-1)+Rawdatas(:,id*2)*1i; % Convert double to complex
    end
    
    datas=reshape((dataiq),nrange,size(Rawdatas,1)/nrange,size(Rawdatas,2)/2);
    
    %% Tx-Rx分割
    datas_2d=zeros(size(datas,1),size(datas,2),2,2);
    if(max_ch==4)
        datas_2d(:,:,1,1)=datas(:,:,1);%Tx1-Rx1
        datas_2d(:,:,1,2)=datas(:,:,2);%Tx1-Rx2
        datas_2d(:,:,2,1)=datas(:,:,3);%Tx2-Rx1
        datas_2d(:,:,2,2)=datas(:,:,4);%Tx2-Rx2
    elseif(max_ch==16)
        datas_2d(:,:,1,1)=datas(:,:,1);%Tx1-Rx1
        datas_2d(:,:,1,2)=datas(:,:,2);%Tx1-Rx2
        datas_2d(:,:,1,3)=datas(:,:,3);%Tx2-Rx1
        datas_2d(:,:,1,4)=datas(:,:,4);%Tx2-Rx2
        datas_2d(:,:,2,1)=datas(:,:,5);%Tx1-Rx1
        datas_2d(:,:,2,2)=datas(:,:,6);%Tx1-Rx2
        datas_2d(:,:,2,3)=datas(:,:,7);%Tx2-Rx1
        datas_2d(:,:,2,4)=datas(:,:,8);%Tx2-Rx2
        datas_2d(:,:,3,1)=datas(:,:,9);%Tx1-Rx1
        datas_2d(:,:,3,2)=datas(:,:,10);%Tx1-Rx2
        datas_2d(:,:,3,3)=datas(:,:,11);%Tx2-Rx1
        datas_2d(:,:,3,4)=datas(:,:,12);%Tx2-Rx2
        datas_2d(:,:,4,1)=datas(:,:,13);%Tx1-Rx1
        datas_2d(:,:,4,2)=datas(:,:,14);%Tx1-Rx2
        datas_2d(:,:,4,3)=datas(:,:,15);%Tx2-Rx1
        datas_2d(:,:,4,4)=datas(:,:,16);%Tx2-Rx2
    end
    %% アダマールキャリブレーション
    datas_cal=zeros(size(datas_2d));
    if(max_ch==4)
        datas_cal(:,:,1,:)=datas_2d(:,:,1,:)+datas_2d(:,:,2,:);
        datas_cal(:,:,2,:)=datas_2d(:,:,1,:)-datas_2d(:,:,2,:);
    elseif(max_ch==16)
        datas_cal(:,:,1,:)=datas_2d(:,:,1,:)+datas_2d(:,:,2,:)+datas_2d(:,:,3,:)+datas_2d(:,:,4,:);
        datas_cal(:,:,2,:)=datas_2d(:,:,1,:)-datas_2d(:,:,2,:)+datas_2d(:,:,3,:)-datas_2d(:,:,4,:);
        datas_cal(:,:,3,:)=datas_2d(:,:,1,:)+datas_2d(:,:,2,:)-datas_2d(:,:,3,:)-datas_2d(:,:,4,:);
        datas_cal(:,:,4,:)=datas_2d(:,:,1,:)-datas_2d(:,:,2,:)-datas_2d(:,:,3,:)+datas_2d(:,:,4,:);
    end
    
    %% 1ファイルに統合
    datas_all(:, 1+(ik-1)*samples:ik*samples, :, :) = datas_cal(:, :, :, :);
    
    % save(strrep(filename, 'bin', 'mat'),'datas_cal', 'date');
    fclose(fd);
end
toc

%% Save
tic
save(filename_out, 'datas_all', 'rdrdate', '-v7.3');
toc

%% 終了