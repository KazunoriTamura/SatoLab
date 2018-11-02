function [ data_synthesized,W ] = Syn( data )
% Time x 16 �܂��� Time x 8 �܂��� Time x N �̓��̓f�[�^�ɑ΂���
% �ŗL�l�W�J����э������s��������̃f�[�^���o��

D = ndims(data);

if D ~= 2
    error('���͂͂Q�����ɂ��ĂˁI')
else
    Xt = data.';   % size(Xt) = 16 x Time
    Rxx = Xt * Xt' / size(Xt,2);
    [V,D] = eig(Rxx);
    [~,lcs] = max(max(D));
    data_synthesized = V(:,lcs)' * Xt;
    W = V(:,lcs);
end


end

