function [output1, output2, time_start] = sync2(input1, input1_start, T1, input2, input2_start, T2)
%   2信号を同期
%   sync(input1, input1_start, T1, input2, input2_start, T2)
%   input1 / input2                     : 1次元または4次元(RANGE, SAMPLE, TX, RX)
%   input1_start / input2_start  : 各データ系列の開始時刻, datetime型
%   T1 / T2                                 : 各データ系列サンプリング間隔
%   syncの改良版　ちゃんと動作検証できてない　誰かやって

%% Duration
if(ndims(input1) == 4 && ndims(input2) == 4)
    %共に(Range, :, Tx, RX)形式の場合
    input1_duration = seconds(T1 * size(input1, 2));
    input2_duration = seconds(T2 * size(input2, 2));
    flag = 1;
elseif(ndims(input1) == 4 && ismatrix(input2))
    % input1が4次元, input2が1次元
    input1_duration = seconds(T1 * size(input1, 2));
    input2_duration = seconds(T2 * numel(input2));
    flag = 2;
elseif(ismatrix(input1) && ndims(input2) == 4)
    % input1が1次元, input2が4次元
    input1_duration = seconds(T1 * numel(input1));
    input2_duration = seconds(T2 * size(input2, 2));
    flag = 3;
elseif(ismatrix(input1) && ismatrix(input2))
    %共に1次元の場合
    input1_duration = seconds(T1 * numel(input1));
    input2_duration = seconds(T2 * numel(input2));
    flag = 4;
else
    error('未知の入力フォーマットです　死ね');
end

%% 終了時刻
input1_end = input1_start + input1_duration;
input2_end = input2_start + input2_duration;

%% 測定時間が重なっていない場合
if(input2_end < input1_start || input1_end < input2_start)
    error('測定時間が重なっていません　死ね');
end

%% 時間オフセット
% input1の時刻を基準とする
offset_sec_before = seconds(input1_start - input2_start);
offset_index_before_input1 = floor(abs(offset_sec_before) / T1);
offset_index_before_input2 = floor(abs(offset_sec_before) / T2);

offset_sec_after = seconds(input1_end - input2_end);
offset_index_after_input1 = floor(abs(offset_sec_after) / T1);
offset_index_after_input2 = floor(abs(offset_sec_after) / T2);

%% 同期
switch flag
    case 1
        % 共に4次元配列
        if offset_sec_before >= 0
            % input2 -> input1の順に開始
            time_start = input1_start;
            
            if offset_sec_after <= 0
                % input1 -> input2の順に終了
                output1 = input1;
                output2 = input2(:, offset_index_before_input2 + 1 : end - offset_index_after_input2, :, :);
                
            else
                % input2 -> input1の順に終了
                output1 = input1(:, 1 : end - offset_index_after_input1, :, :);
                output2 = input2(:, offset_index_before_input2 : end, :, :);
            end
            
        else
            % input1 -> input2の順に開始
            time_start = input2_start;
            
            if offset_sec_after <= 0
                % input1 -> input2の順に終了
                output1 = input1(:, offset_index_before_input1 + 1 : end, :, :);
                output2 = input2(:, 1 : end - offset_index_after_input2, :, :);
                
            else
                % input2 -> input1の順に終了
                output1 = input1(:, offset_index_before_input1 + 1 : end - offset_index_after_input1, :, :);
                output2 = input2;
            end
        end
        
        
    case 2
        % input1が4次元, input2が1次元
        if offset_sec_before >= 0
            % input2 -> input1の順に開始
            time_start = input1_start;
            
            if offset_sec_after <= 0
                % input1 -> input2の順に終了
                output1 = input1;
                output2 = input2(offset_index_before_input2 + 1 : end - offset_index_after_input2);
                
            else
                % input2 -> input1の順に終了
                output1 = input1(:, 1 : end - offset_index_after_input1, :, :);
                output2 = input2(offset_index_before_input2 + 1 : end);
            end
            
        else
            % input1 -> input2の順に開始
            time_start = input2_start;
            
            if offset_sec_after <= 0
                % input1 -> input2の順に終了
                output1 = input1(:, offset_index_before_input1 + 1 : end, :, :);
                output2 = input2(1 : end - offset_index_after_input2);
                
            else
                % input2 -> input1の順に終了
                output1 = input1(:, offset_index_before_input1 + 1 : end - offset_index_after_input1, :, :);
                output2 = input2;
            end
        end
        
        
    case 3
        % input1が1次元, input2が4次元
        if offset_sec_before >= 0
            % input2 -> input1の順に開始
            time_start = input1_start;
            
            if offset_sec_after <= 0
                % input1 -> input2の順に終了
                output1 = input1;
                output2 = input2(:, offset_index_before_input2 + 1 : end - offset_index_after_input2, :, :);
                
            else
                % input2 -> input1の順に終了
                output1 = input1(1 : end - offset_index_after_input1);
                output2 = input2(:, offset_index_before_input2 + 1 : end, :, :);
            end
            
        else
            % input1 -> input2の順に開始
            time_start = input2_start;
            
            if offset_sec_after <= 0
                % input1 -> input2の順に終了
                output1 = input1(offset_index_before_input1 + 1 : end);
                output2 = input2(:, 1 : end - offset_index_after_input2, :, :);
                
            else
                % input2 -> input1の順に終了
                output1 = input1(offset_index_before_input1 + 1 : end - offset_index_after_input1);
                output2 = input2;
            end
        end

        
    case 4
        % 共に1次元配列
        if offset_sec_before >= 0
            % input2 -> input1の順に開始
            time_start = input1_start;
            
            if offset_sec_after <= 0
                % input1 -> input2の順に終了
                output1 = input1;
                output2 = input2(offset_index_before_input2 + 1 : end - offset_index_after_input2);
                
            else
                % input2 -> input1の順に終了
                output1 = input1(1 : end - offset_index_after_input1);
                output2 = input2(offset_index_before_input2 + 1 : end);
            end
            
        else
            % input1 -> input2の順に開始
            time_start = input2_start;
            
            if offset_sec_after <= 0
                % input1 -> input2の順に終了
                output1 = input1(offset_index_before_input1 + 1 : end);
                output2 = input2(1 : end - offset_index_after_input2);
                
            else
                % input2 -> input1の順に終了
                output1 = input1(offset_index_before_input1 + 1 : end - offset_index_after_input1);
                output2 = input2;
            end
        end
        
end