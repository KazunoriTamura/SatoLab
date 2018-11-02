function [ radartime_correct, radarforplot_correct, radartime_wrong, radarforplot_wrong, Error, Error_correct ] = Judgement( m_time, radartime_all, radarforplot_all, ecgs )
% ���[�_�[�ɂ�鐄��l��ECG�ɂ��^�l�Ƃ̔�r
% CAR = Correct Answer Rate = ������ = �c���擾�� = ������/�S����_�� = �덷50msec�ȓ��̐���_��/�S����_��

%% �p�����[�^
allowance = 0.050;
decimal = 6;

%% �O����
stride = 1/10^(decimal);

%% �v�Z
N_all = length(radarforplot_all);

N_correct = 0;
N_wrong = 0;

newtime = 0:stride:m_time(end);
ecg_interp = interp1(ecgs.ecgtime, ecgs.ecgplot, newtime, 'liner', 'extrap'); % newtime�ɑ΂���ECG�̒l

radartime_correct = zeros(1,size(radartime_all,2));
radarforplot_correct = zeros(1,size(radartime_all,2));
radartime_wrong = zeros(1,size(radartime_all,2));
radarforplot_wrong = zeros(1,size(radartime_all,2));
Error = zeros(1,size(radartime_all,2));
Error_correct = zeros(1,size(radartime_all,2));
for i = 1:size(radartime_all,2)
    index = int64( radartime_all(i) * 10^(decimal) );
    true_value = ecg_interp(index);
    if isempty(true_value) == true
        disp('Judgement�֐����ɃG���[����')
    end
    estimation_value = radarforplot_all(i);
    error = estimation_value - true_value;
    Error(i) = error;
    if abs(error) < true_value * allowance 
        N_correct = N_correct + 1;
        radarforplot_correct(i) = estimation_value;
        radartime_correct(i) = radartime_all(i);
        Error_correct(i) = error;
    else
        N_wrong = N_wrong + 1;
        radarforplot_wrong(i) = estimation_value;
        radartime_wrong(i) = radartime_all(i);
    end
end
k_correct = find(radarforplot_correct == 0);
radarforplot_correct(k_correct) = [];
radartime_correct(k_correct) = [];
Error_correct(k_correct) = [];
k_wrong = find(radarforplot_wrong == 0);
radarforplot_wrong(k_wrong) = [];
radartime_wrong(k_wrong) = [];

%% �G���[�`�F�b�N
N_sum = N_correct + N_wrong;
if N_sum ~= N_all
   disp('Judgement�֐����ɃG���[����')
end

end

