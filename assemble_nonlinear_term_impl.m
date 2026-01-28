function F = assemble_nonlinear_term_impl(u,T,Nd,nodes2mesh,parallel)
%ASSEMBLE_NONLINEAR_TERM_IMPL Summary of this function goes here
%   Detailed explanation goes here

if parallel
    F = assemble_nonlinear_term_impl_parallel(u,T,Nd,nodes2mesh);
else
    F = assemble_nonlinear_term_impl_sequential(u,T,Nd,nodes2mesh);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = assemble_nonlinear_term_impl_sequential(u,T,Nd,nodes2mesh)
% sequential implementation
F_i = zeros(9*size(T,1),1);
F_j = zeros(9*size(T,1),1);
FF = zeros(9*size(T,1),1);
ind = 1;

phi = { @(x) -x(1) - x(2) + 1;
    @(x) x(1);
    @(x) x(2)};

[quad,w] = getQuadrature(7);
no_of_basis = size(phi,1);
no_of_quad_points = length(w);

phi_in_quad = zeros(no_of_quad_points,no_of_basis);
for j = 1:no_of_quad_points
    for i = 1:no_of_basis
        phi_in_quad(j,i) = phi{i}(quad(:,j));
    end
end

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node

    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];

    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

    %% construct f = |u|^2 on triangle
    f_in_quad = zeros(no_of_quad_points,1);

    for j = 1:no_of_basis
        index = nodes2mesh(tri(j));
        if index ~= 0
            f_in_quad = f_in_quad + u(index)*phi_in_quad(:,j);
        end
    end

    f_in_quad = abs(f_in_quad).^2;

    %% assemble
    for i = 1:3
        for j = i:3
            e = detBT*evaluateQuadrature(transpose(f_in_quad.*phi_in_quad(:,i).*phi_in_quad(:,j)),w);

            if i == j
                F_i(ind) = nodes2mesh(tri(i));
                F_j(ind) = nodes2mesh(tri(j));
                FF(ind) = e;
                ind = ind + 1;
            else
                F_i(ind) = nodes2mesh(tri(i));
                F_j(ind) = nodes2mesh(tri(j));
                FF(ind) = e;
                ind = ind + 1;

                F_i(ind) = nodes2mesh(tri(j));
                F_j(ind) = nodes2mesh(tri(i));
                FF(ind) = e;
                ind = ind + 1;
            end
        end
    end
end

del_index = find(F_i==0,1);
if isempty(del_index)
    del_index = length(F_i)+1;
end
F = sparse(F_i(1:del_index-1),F_j(1:del_index-1),FF(1:del_index-1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = assemble_nonlinear_term_impl_parallel(u,T,Nd,nodes2mesh)
% parallel implementation
dim = size(Nd,1);

phi = { @(x) -x(1) - x(2) + 1;
    @(x) x(1);
    @(x) x(2)};

[quad,w] = getQuadrature(7);
no_of_basis = size(phi,1);
no_of_quad_points = length(w);

phi_in_quad = zeros(no_of_quad_points,no_of_basis);
for j = 1:no_of_quad_points
    for i = 1:no_of_basis
        phi_in_quad(j,i) = phi{i}(quad(:,j));
    end
end

spmd
    my_index_start = floor(size(T,1)*(spmdIndex - 1)/spmdSize + 1);
    my_index_end = floor(size(T,1)*(spmdIndex)/spmdSize);

    F_i = zeros(9*size(T,1),1);
    F_j = zeros(9*size(T,1),1);
    FF = zeros(9*size(T,1),1);
    ind = 1;

    for k = my_index_start:my_index_end
        tri = T(k,:); %node index of triangle
        z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
        z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
        z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node

        %% transformation to refenrence triangle
        BT = [z2(1)-z1(1), z3(1)-z1(1); ...
            z2(2)-z1(2), z3(2)-z1(2)];

        detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

        %% construct f = |u|^2 on triangle
        f_in_quad = zeros(no_of_quad_points,1);

        for j = 1:no_of_basis
            index = nodes2mesh(tri(j));
            if index ~= 0
                f_in_quad = f_in_quad + u(index)*phi_in_quad(:,j);
            end
        end

        f_in_quad = abs(f_in_quad).^2;

        %% assemble
        for i = 1:3
            for j = i:3
                e = detBT*evaluateQuadrature(transpose(f_in_quad.*phi_in_quad(:,i).*phi_in_quad(:,j)),w);

                if i == j
                    F_i(ind) = nodes2mesh(tri(i));
                    F_j(ind) = nodes2mesh(tri(j));
                    FF(ind) = e;
                    ind = ind + 1;
                else
                    F_i(ind) = nodes2mesh(tri(i));
                    F_j(ind) = nodes2mesh(tri(j));
                    FF(ind) = e;
                    ind = ind + 1;

                    F_i(ind) = nodes2mesh(tri(j));
                    F_j(ind) = nodes2mesh(tri(i));
                    FF(ind) = e;
                    ind = ind + 1;
                end
            end
        end
    end

    del_index = find(F_i==0,1);
    if isempty(del_index)
        del_index = length(F_i)+1;
    end

    F_worker = sparse(F_i(1:del_index-1),F_j(1:del_index-1),FF(1:del_index-1),dim,dim);

end

% collect results from workers
F = sparse(dim,dim);
for j = 1:length(my_index_start)
    F = F + F_worker{j};
end

end


