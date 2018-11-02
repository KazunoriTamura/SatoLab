function [ v ] = steeringVector( theta, phi, lambda, r )
% ステアリングベクトルの計算
L = lookDirectionVector(theta, phi);
s = r.' * L;
v = exp(1i*2*pi*s/lambda);
end