function [ outputs ] = extract8Chs( input )
% range x time x 16 または anything x 16 のデータから
% クロスチップチャネルを除いて
% range x time x 8 または anything x 8 のデータを取り出す

D = ndims(input);

if D ~= 2 && D ~= 3
    error('入力は２次元または３次元で！');
end

if D == 3
    U = size(input,3);
    if U ~= 16
        error('エラー！');
    end
    
    output = zeros(size(input,1),size(input,2),8);
    
    output(:,:,1) = input(:,:,1);
    output(:,:,2) = input(:,:,2);
    output(:,:,3) = input(:,:,5);
    output(:,:,4) = input(:,:,6);
    
    output(:,:,5) = input(:,:,11);
    output(:,:,6) = input(:,:,12);
    output(:,:,7) = input(:,:,15);
    output(:,:,8) = input(:,:,16);
end
    
if D == 2
    U = size(input,2);
    if U ~= 16
        error('エラー！');
    end
    
    output = zeros(size(input,1),8);
    
    output(:,1) = input(:,1);
    output(:,2) = input(:,2);
    output(:,3) = input(:,5);
    output(:,4) = input(:,6);
    
    output(:,5) = input(:,11);
    output(:,6) = input(:,12);
    output(:,7) = input(:,15);
    output(:,8) = input(:,16);
end

outputs = output;

end

