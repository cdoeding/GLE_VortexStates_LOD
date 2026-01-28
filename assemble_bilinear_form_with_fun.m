function S = assemble_bilinear_form_with_fun(A,kappa,T,Nd,parallel)
%ASSEMBLEGLOBALBILINEARFORM Summary of this function goes here
%   Detailed explanation goes here

if parallel
    S = assemble_bilinear_form_with_fun_parallel(A,kappa,T,Nd);
else
    S = assemble_bilinear_form_with_fun_sequential(A,kappa,T,Nd);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = assemble_bilinear_form_with_fun_sequential(A,kappa,T,Nd)
% sequential implementation

s_i = zeros(9*size(T,1),1);
s_j = zeros(9*size(T,1),1);
s_val = zeros(9*size(T,1),1);
ind = 1;

grad = {[-1; -1]; [1; 0]; [0; 1]};

phi = { @(x) -x(1) - x(2) + 1;
    @(x) x(1);
    @(x) x(2)};

[quad,w] = getQuadrature(7);
no_of_basis = size(grad,1);
no_of_quad_points = length(w);

grad_in_quad = zeros(2,no_of_quad_points,no_of_basis);
phi_in_quad = zeros(no_of_quad_points,no_of_basis);

for j = 1:no_of_quad_points
    for i = 1:no_of_basis
        grad_in_quad(:,j,i) = grad{i};
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

    b = [z1(1); z1(2)];

    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

    BTinv = inv(BT)';
    grad_in_BT = pagemtimes(BTinv,grad_in_quad);

    %% construct A on triangle
    A_in_quad = zeros(2,no_of_quad_points);

    for i = 1:no_of_quad_points
        A_in_quad(:,i) = A(BT*quad(:,i) + b);
    end

    %% assemble
    for i = 1:3
        for j = i:3
            int = diag((1i/kappa * grad_in_BT(:,:,i)+ phi_in_quad(:,i)'.*A_in_quad)'*(1i/kappa * grad_in_BT(:,:,j)+ phi_in_quad(:,j)'.*A_in_quad));
            e = detBT*evaluateQuadrature(transpose(int),w);
            if i == j
                s_i(ind) = tri(i);
                s_j(ind) = tri(i);
                s_val(ind) = e;
                ind = ind + 1;
            else
                s_i(ind) = tri(i);
                s_j(ind) = tri(j);
                s_val(ind) = conj(e);
                ind = ind + 1;

                s_i(ind) = tri(j);
                s_j(ind) = tri(i);
                s_val(ind) = e;
                ind = ind + 1;
            end
        end
    end
end

S = sparse(s_i,s_j,s_val);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = assemble_bilinear_form_with_fun_parallel(A,kappa,T,Nd)
% parallel implementation
dim = size(Nd,1);

grad = {[-1; -1]; [1; 0]; [0; 1]};

phi = { @(x) -x(1) - x(2) + 1;
    @(x) x(1);
    @(x) x(2)};

[quad,w] = getQuadrature(7);
no_of_basis = size(grad,1);
no_of_quad_points = length(w);

grad_in_quad = zeros(2,no_of_quad_points,no_of_basis);
phi_in_quad = zeros(no_of_quad_points,no_of_basis);

for j = 1:no_of_quad_points
    for i = 1:no_of_basis
        grad_in_quad(:,j,i) = grad{i};
        phi_in_quad(j,i) = phi{i}(quad(:,j));
    end
end

spmd
    my_index_start = floor(size(T,1)*(spmdIndex - 1)/spmdSize + 1);
    my_index_end = floor(size(T,1)*(spmdIndex)/spmdSize);

    s_i = zeros(9*size(T,1),1);
    s_j = zeros(9*size(T,1),1);
    s_val = zeros(9*size(T,1),1);
    ind = 1;

    for k = my_index_start:my_index_end
        tri = T(k,:); %node index of triangle
        z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
        z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
        z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node

        %% transformation to refenrence triangle
        BT = [z2(1)-z1(1), z3(1)-z1(1); ...
            z2(2)-z1(2), z3(2)-z1(2)];

        b = [z1(1); z1(2)];

        detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);

        BTinv = inv(BT)';
        grad_in_BT = pagemtimes(BTinv,grad_in_quad);

        %% construct A on triangle
        A_in_quad = zeros(2,no_of_quad_points);

        for i = 1:no_of_quad_points
            A_in_quad(:,i) = A(BT*quad(:,i) + b);
        end

        %% assemble
        for i = 1:3
            for j = i:3
                int = diag((1i/kappa * grad_in_BT(:,:,i)+ phi_in_quad(:,i)'.*A_in_quad)'*(1i/kappa * grad_in_BT(:,:,j)+ phi_in_quad(:,j)'.*A_in_quad));
                e = detBT*evaluateQuadrature(transpose(int),w);
                if i == j
                    s_i(ind) = tri(i);
                    s_j(ind) = tri(i);
                    s_val(ind) = e;
                    ind = ind + 1;
                else
                    s_i(ind) = tri(i);
                    s_j(ind) = tri(j);
                    s_val(ind) = conj(e);
                    ind = ind + 1;

                    s_i(ind) = tri(j);
                    s_j(ind) = tri(i);
                    s_val(ind) = e;
                    ind = ind + 1;
                end
            end
        end
    end

    del_index = find(s_i==0,1);
    if isempty(del_index)
        del_index = length(s_i)+1;
    end
    S_worker = sparse(s_i(1:del_index-1),s_j(1:del_index-1),s_val(1:del_index-1),dim,dim);

end

% collect results from workers
S = sparse(dim,dim);
for j = 1:length(my_index_start)
    S = S + S_worker{j};
end

end

