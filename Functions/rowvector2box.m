function [ Box ] = rowvector2box( rowVector )
% 1x16�̍s�x�N�g����4x4�̍s��֕ϊ�����
% �������C�ϊ��̍ۂɏ��Ԃ��u�f�q�ԍ��Z�b�g�ƃ`�����l���i���o�[�̑Ή��֌W�v�ɏ����������̂ƂȂ�

D = ndims(rowVector);
R = size(rowVector,1);

if D ~= 2
    error('���͂�1x16�ŁI');
end

if R ~= 1
    error('���͂�1x16�ŁI');
end

Box = zeros(4,4);
for iNo = 1:16
    
    [Tx,Rx] = no2txrx(iNo);
    Box(Tx,Rx) = rowVector(iNo);
    
end


end