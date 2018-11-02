% 20171026 Release
% H5�t�@�C�����g�p
% �����Z���T�[�ɂ��Ή��@���Ԃ�@�����ƃ`�F�b�N�͂��ĂȂ�
%
% 20180513���� by �c��
% 20180919���� by �c�� because Group����������ꍇ�̏������}�j���A���I�ɒǉ�
%
%% �J�n
tic
clear
fld = dir('opensignals_*.h5');
ECGfilename = fld.name;

%% �p�����[�^
group = 1;

%% Get ECG parameters
ECGinfo = h5info(ECGfilename);
for ik = 1:size(ECGinfo.Groups(group).Attributes, 1)
    if strcmp(ECGinfo.Groups(group).Attributes(ik).Name, 'device name')
        ECG.name = char(ECGinfo.Groups(group).Attributes(ik).Value);
        
    elseif strcmp(ECGinfo.Groups(group).Attributes(ik).Name, 'sampling rate')
        ECG.sample = double(ECGinfo.Groups(group).Attributes(ik).Value);
        
    elseif strcmp(ECGinfo.Groups(group).Attributes(ik).Name, 'resolution')
        ECG.resolution = double(ECGinfo.Groups(group).Attributes(ik).Value);
        
    elseif strcmp(ECGinfo.Groups(group).Attributes(ik).Name, 'date')
        ECG.date = datetime(strcat(ECGinfo.Groups(group).Attributes(ik).Value,'_',ECGinfo.Groups(group).Attributes(ik+1).Value), 'InputFormat', 'yyyy-MM-dd_HH:m:ss.SSS');
    
    elseif strcmp(ECGinfo.Groups(group).Attributes(ik).Name, 'channels')
        usechannel = ECGinfo.Groups(group).Attributes(ik).Value;
        
    elseif strcmp(ECGinfo.Groups(group).Attributes(ik).Name, 'nsamples')
        ECG.nsamples = double(ECGinfo.Groups(group).Attributes(ik).Value);
    end
end

ECG.time = (1:ECG.nsamples)/ECG.sample;

%% Get ECG data
% tmp = 1;
sensor_list = cell(1, size(usechannel,1));
for ik = 1 : size(usechannel,1)
    % �Z���T�[���擾
    sensor = h5readatt(ECGfilename, strcat('/', char(ECG.name), '/raw/channel_', num2str(usechannel(ik))), 'sensor');
    
    % �Z���T�[���̏d�����
    if ik >= 2
        counter = 1;
        for il = 1 : ik - 1
            if contains(char(sensor_list(il)), char(sensor))
                counter = counter + 1;
            end
        end
        
        if counter > 1
            sensor = cellstr(strcat(sensor, num2str(counter)));
        end
    end
    sensor_list(ik) = sensor;

    % �f�[�^�擾
    data = h5read(ECGfilename, strcat('/', char(ECG.name), '/raw/channel_', num2str(usechannel(ik))));
    data = double(data) / 2^ECG.resolution(ik); %���K��
    sensor_list_tmp = char(sensor_list(ik));
    sensor_list_tmp = strrep(sensor_list_tmp,'/','');
    sensor_list_tmp = strrep(sensor_list_tmp,'.','');
    ECG.datas.(sensor_list_tmp) = detrend(data); % �g�����h����
end

save(strcat('ECG_Group',num2str(group),'_',datestr(ECG.date, 'yyyymmdd_HHMMSS'), '.mat'), 'ECG');

%% IBI�ɕϊ�
% [ecgtime, ecgplot] = ecg2ibi(ECG.datas.ECG, ECG.sample);
% 
% save(strcat('ECG_', datestr(ECG.date, 'yyyymmdd_HHMMSS'), '_topo.mat'), 'ecgplot', 'ecgtime');

%% �I��
toc