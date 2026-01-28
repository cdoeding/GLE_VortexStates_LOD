function A = assemble_stiffness_matrix(T,Nd,nodes2mesh)
%ASSEMBLE_STIFFNESS_MATRIX Summary of this function goes here
%   Detailed explanation goes here

a_i = zeros(9*size(T,1),1);
a_j = zeros(9*size(T,1),1);
a = zeros(9*size(T,1),1);
ind = 1;

grad = [-1, 1, 0; -1, 0, 1];

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node

    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    BTinv = inv(BT);
    D = BTinv*(BTinv');
    A_loc = 0.5*detBT*grad'*D*grad;

    %% assemble
    for i = 1:3
        for j = 1:3
            a_i(ind) = nodes2mesh(tri(i),1);
            a_j(ind) = nodes2mesh(tri(j),1);
            a(ind) = A_loc(i,j);
            ind = ind + 1;
        end
    end
end

del_index = find(a_i==0,1);
if isempty(del_index)
    del_index = length(a_i)+1;
end
A = sparse(a_i(1:del_index-1),a_j(1:del_index-1),a(1:del_index-1));
end