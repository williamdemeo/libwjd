% read in the data from file data.m
data
p = size(X,2);
[Q,R,E] = qr(X);
Qy = Q'*Y;
B = E * inv(R(1:p,1:p)) * Qy(1:p);
