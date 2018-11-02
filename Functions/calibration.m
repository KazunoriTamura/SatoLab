function [ outputs ] = calibration( data, rangeCr )
% キャリブレーション
% CRレンジの位相を方位基準にする

D = ndims(data);
R = size(data,1);

if D ~= 4
    error('入力は４次元で！');
end

if R == 1
    error('入力は複数レンジ分で！');
end


data_reshaped = reshapeTxrx2No(data);

outputs = zeros(size(data));
for iRange = 1:R
    
    X = squeeze(data_reshaped(iRange,:,:)).';
    X_Cr = squeeze(data_reshaped(rangeCr,:,:)).';
    
    angle_Cr = unwrap(angle(X_Cr));
    phi = mean(angle_Cr,2);
    
    X_calibrated = X .* exp( -1i .* phi );
    
    output_reshaped = zeros(1,size(X_calibrated,2),size(X_calibrated,1));
    output_reshaped(1,:,:) = X_calibrated.';
    
    output = reshapeNo2Txrx(output_reshaped);
    
    outputs(iRange,:,:,:) = output;
end


end

