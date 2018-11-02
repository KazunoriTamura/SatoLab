function [ L ] = lookDirectionVector( theta, phi )
% 視線方向ベクトル

x = sin(theta) * cos(phi);
y = sin(theta) * sin(phi);
z = cos(theta);

L = [x,y,z].';

end

