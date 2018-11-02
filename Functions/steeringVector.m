function [ v ] = steeringVector( theta, phi, lambda, r )
% �X�e�A�����O�x�N�g���̌v�Z
L = lookDirectionVector(theta, phi);
s = r.' * L;
v = exp(1i*2*pi*s/lambda);
end