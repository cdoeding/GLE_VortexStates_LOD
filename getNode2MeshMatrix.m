function T = getNode2MeshMatrix(T_H,l)
%GETNODE2MESHMATRIX Summary of this function goes here
%   Detailed explanation goes here

N_H = size(T_H.p,1); % number of coarse nodes
T_i = zeros(3,1);
T_j = zeros(3,1);
for k = 1:3
   T_i(k) = k;
   T_j(k) = T_H.t(l,k);
end
T = sparse(T_i,T_j,[1,1,1],3,N_H);

end

