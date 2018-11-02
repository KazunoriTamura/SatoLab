function [ rowVector ] = box2rowvector( Box )
% 4x4の行列を1x16の行ベクトルへ変換する
% ただし，変換の際に順番が「素子番号セットとチャンネルナンバーの対応関係」に準拠したものとなる

D = ndims(Box);
R = size(Box,1);

if D ~= 2
    error('入力は4x4で！');
end

if R ~= 4
    error('入力は4x4で！');
end

rowVector = zeros(1,16);
for iTx = 1:4
    for iRx = 1:4
        No = txrx2no(iTx,iRx);
        rowVector(No) = Box(iTx,iRx);
    end
end


end

