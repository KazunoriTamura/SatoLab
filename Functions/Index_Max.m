function [ I_row, I_col ] = Index_Max( A )
% 行列の最大要素を与えるインデックスを求める
[~,I]=max(A(:));
[I_row,I_col] = ind2sub(size(A),I);
end

