function [output1, output2, time_start] = sync2(input1, input1_start, T1, input2, input2_start, T2)
%   2�M���𓯊�
%   sync(input1, input1_start, T1, input2, input2_start, T2)
%   input1 / input2                     : 1�����܂���4����(RANGE, SAMPLE, TX, RX)
%   input1_start / input2_start  : �e�f�[�^�n��̊J�n����, datetime�^
%   T1 / T2                                 : �e�f�[�^�n��T���v�����O�Ԋu
%   sync�̉��ǔŁ@�����Ɠ��쌟�؂ł��ĂȂ��@�N�������

%% Duration
if(ndims(input1) == 4 && ndims(input2) == 4)
    %����(Range, :, Tx, RX)�`���̏ꍇ
    input1_duration = seconds(T1 * size(input1, 2));
    input2_duration = seconds(T2 * size(input2, 2));
    flag = 1;
elseif(ndims(input1) == 4 && ismatrix(input2))
    % input1��4����, input2��1����
    input1_duration = seconds(T1 * size(input1, 2));
    input2_duration = seconds(T2 * numel(input2));
    flag = 2;
elseif(ismatrix(input1) && ndims(input2) == 4)
    % input1��1����, input2��4����
    input1_duration = seconds(T1 * numel(input1));
    input2_duration = seconds(T2 * size(input2, 2));
    flag = 3;
elseif(ismatrix(input1) && ismatrix(input2))
    %����1�����̏ꍇ
    input1_duration = seconds(T1 * numel(input1));
    input2_duration = seconds(T2 * numel(input2));
    flag = 4;
else
    error('���m�̓��̓t�H�[�}�b�g�ł��@����');
end

%% �I������
input1_end = input1_start + input1_duration;
input2_end = input2_start + input2_duration;

%% ���莞�Ԃ��d�Ȃ��Ă��Ȃ��ꍇ
if(input2_end < input1_start || input1_end < input2_start)
    error('���莞�Ԃ��d�Ȃ��Ă��܂���@����');
end

%% ���ԃI�t�Z�b�g
% input1�̎�������Ƃ���
offset_sec_before = seconds(input1_start - input2_start);
offset_index_before_input1 = floor(abs(offset_sec_before) / T1);
offset_index_before_input2 = floor(abs(offset_sec_before) / T2);

offset_sec_after = seconds(input1_end - input2_end);
offset_index_after_input1 = floor(abs(offset_sec_after) / T1);
offset_index_after_input2 = floor(abs(offset_sec_after) / T2);

%% ����
switch flag
    case 1
        % ����4�����z��
        if offset_sec_before >= 0
            % input2 -> input1�̏��ɊJ�n
            time_start = input1_start;
            
            if offset_sec_after <= 0
                % input1 -> input2�̏��ɏI��
                output1 = input1;
                output2 = input2(:, offset_index_before_input2 + 1 : end - offset_index_after_input2, :, :);
                
            else
                % input2 -> input1�̏��ɏI��
                output1 = input1(:, 1 : end - offset_index_after_input1, :, :);
                output2 = input2(:, offset_index_before_input2 : end, :, :);
            end
            
        else
            % input1 -> input2�̏��ɊJ�n
            time_start = input2_start;
            
            if offset_sec_after <= 0
                % input1 -> input2�̏��ɏI��
                output1 = input1(:, offset_index_before_input1 + 1 : end, :, :);
                output2 = input2(:, 1 : end - offset_index_after_input2, :, :);
                
            else
                % input2 -> input1�̏��ɏI��
                output1 = input1(:, offset_index_before_input1 + 1 : end - offset_index_after_input1, :, :);
                output2 = input2;
            end
        end
        
        
    case 2
        % input1��4����, input2��1����
        if offset_sec_before >= 0
            % input2 -> input1�̏��ɊJ�n
            time_start = input1_start;
            
            if offset_sec_after <= 0
                % input1 -> input2�̏��ɏI��
                output1 = input1;
                output2 = input2(offset_index_before_input2 + 1 : end - offset_index_after_input2);
                
            else
                % input2 -> input1�̏��ɏI��
                output1 = input1(:, 1 : end - offset_index_after_input1, :, :);
                output2 = input2(offset_index_before_input2 + 1 : end);
            end
            
        else
            % input1 -> input2�̏��ɊJ�n
            time_start = input2_start;
            
            if offset_sec_after <= 0
                % input1 -> input2�̏��ɏI��
                output1 = input1(:, offset_index_before_input1 + 1 : end, :, :);
                output2 = input2(1 : end - offset_index_after_input2);
                
            else
                % input2 -> input1�̏��ɏI��
                output1 = input1(:, offset_index_before_input1 + 1 : end - offset_index_after_input1, :, :);
                output2 = input2;
            end
        end
        
        
    case 3
        % input1��1����, input2��4����
        if offset_sec_before >= 0
            % input2 -> input1�̏��ɊJ�n
            time_start = input1_start;
            
            if offset_sec_after <= 0
                % input1 -> input2�̏��ɏI��
                output1 = input1;
                output2 = input2(:, offset_index_before_input2 + 1 : end - offset_index_after_input2, :, :);
                
            else
                % input2 -> input1�̏��ɏI��
                output1 = input1(1 : end - offset_index_after_input1);
                output2 = input2(:, offset_index_before_input2 + 1 : end, :, :);
            end
            
        else
            % input1 -> input2�̏��ɊJ�n
            time_start = input2_start;
            
            if offset_sec_after <= 0
                % input1 -> input2�̏��ɏI��
                output1 = input1(offset_index_before_input1 + 1 : end);
                output2 = input2(:, 1 : end - offset_index_after_input2, :, :);
                
            else
                % input2 -> input1�̏��ɏI��
                output1 = input1(offset_index_before_input1 + 1 : end - offset_index_after_input1);
                output2 = input2;
            end
        end

        
    case 4
        % ����1�����z��
        if offset_sec_before >= 0
            % input2 -> input1�̏��ɊJ�n
            time_start = input1_start;
            
            if offset_sec_after <= 0
                % input1 -> input2�̏��ɏI��
                output1 = input1;
                output2 = input2(offset_index_before_input2 + 1 : end - offset_index_after_input2);
                
            else
                % input2 -> input1�̏��ɏI��
                output1 = input1(1 : end - offset_index_after_input1);
                output2 = input2(offset_index_before_input2 + 1 : end);
            end
            
        else
            % input1 -> input2�̏��ɊJ�n
            time_start = input2_start;
            
            if offset_sec_after <= 0
                % input1 -> input2�̏��ɏI��
                output1 = input1(offset_index_before_input1 + 1 : end);
                output2 = input2(1 : end - offset_index_after_input2);
                
            else
                % input2 -> input1�̏��ɏI��
                output1 = input1(offset_index_before_input1 + 1 : end - offset_index_after_input1);
                output2 = input2;
            end
        end
        
end