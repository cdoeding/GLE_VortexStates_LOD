function M = assembleLocalMassMatrix(T_h,P0,l)
%ASSEMBLELOCALMASMATRIX Summary of this function goes here
%   Detailed explanation goes here

active_fine_triangles = find(logical(P0(:,l))); % indicies of fine triangles in current coarse triangle
T = T_h.t(active_fine_triangles,:); % fine trianlges in current coarse triangle
Nd = T_h.p; % fine nodes of fine trianlges
N_h = size(T_h.p,1); % number of fine nodes

m_i = zeros(9*size(T,1),1);
m_j = zeros(9*size(T,1),1);
m = zeros(9*size(T,1),1);
ind = 1;

M_loc = (1/24)*[2 1 1; 1 2 1; 1 1 2];

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node

    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

    for i = 1:3
        for j = 1:3
            m_i(ind) = tri(i);
            m_j(ind) = tri(j);
            m(ind) = detBT*M_loc(i,j);
            ind = ind + 1;
        end
    end
end

M = sparse(m_i,m_j,m,N_h,N_h);

end

