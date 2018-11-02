function [ data_synthesized,W ] = Syn( data )
% Time x 16 または Time x 8 または Time x N の入力データに対して
% 固有値展開および合成を行い統合後のデータを出力

D = ndims(data);

if D ~= 2
    error('入力は２次元にしてね！')
else
    Xt = data.';   % size(Xt) = 16 x Time
    Rxx = Xt * Xt' / size(Xt,2);
    [V,D] = eig(Rxx);
    [~,lcs] = max(max(D));
    data_synthesized = V(:,lcs)' * Xt;
    W = V(:,lcs);
end


end

