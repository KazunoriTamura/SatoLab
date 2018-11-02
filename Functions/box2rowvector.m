function [ rowVector ] = box2rowvector( Box )
% 4x4�̍s���1x16�̍s�x�N�g���֕ϊ�����
% �������C�ϊ��̍ۂɏ��Ԃ��u�f�q�ԍ��Z�b�g�ƃ`�����l���i���o�[�̑Ή��֌W�v�ɏ����������̂ƂȂ�

D = ndims(Box);
R = size(Box,1);

if D ~= 2
    error('���͂�4x4�ŁI');
end

if R ~= 4
    error('���͂�4x4�ŁI');
end

rowVector = zeros(1,16);
for iTx = 1:4
    for iRx = 1:4
        No = txrx2no(iTx,iRx);
        rowVector(No) = Box(iTx,iRx);
    end
end


end

