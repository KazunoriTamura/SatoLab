function [ Box ] = rowvector2box( rowVector )
% 1x16の行ベクトルを4x4の行列へ変換する
% ただし，変換の際に順番が「素子番号セットとチャンネルナンバーの対応関係」に準拠したものとなる

D = ndims(rowVector);
R = size(rowVector,1);

if D ~= 2
    error('入力は1x16で！');
end

if R ~= 1
    error('入力は1x16で！');
end

Box = zeros(4,4);
for iNo = 1:16
    
    [Tx,Rx] = no2txrx(iNo);
    Box(Tx,Rx) = rowVector(iNo);
    
end


end