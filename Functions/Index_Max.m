function [ I_row, I_col ] = Index_Max( A )
% �s��̍ő�v�f��^����C���f�b�N�X�����߂�
[~,I]=max(A(:));
[I_row,I_col] = ind2sub(size(A),I);
end

