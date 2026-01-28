function [F_un_dn,F_dn_dn] = assembleNonlinearMatrices_un_dn_P1(u,d,T,Nd,parallel)
%   Summary of this function goes here
%   Detailed explanation goes here

% assemble the matrices F_un_dn and F_dn_dn with entries of the form
% ( Re( u^n (d^n)* ) phi_i , phi_j )_L2   and   ( |d^n|^2 phi_i , phi_j )_L2
% respectively

% the local part of the "nonlinearity matrix" is a 3x3 matrix. 
% we will store the values in a vector of length 3x3=9.

if parallel
    [F_un_dn,F_dn_dn] = assembleNonlinearMatrices_un_dn_P1_parallel(u,d,T,Nd);
else
    [F_un_dn,F_dn_dn] = assembleNonlinearMatrices_un_dn_P1_sequential(u,d,T,Nd);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F_un_dn,F_dn_dn] = assembleNonlinearMatrices_un_dn_P1_sequential(u,d,T,Nd)
F_un_dn_i = zeros(9*size(T,1),1); % stores the global i-index of the matrix
F_un_dn_j = zeros(9*size(T,1),1); % stores the global j-index of the matrix
F_un_dn_val = zeros(9*size(T,1),1); % stores the value of the element contribution

F_dn_dn_i = zeros(9*size(T,1),1); % stores the global i-index of the matrix
F_dn_dn_j = zeros(9*size(T,1),1); % stores the global j-index of the matrix
F_dn_dn_val = zeros(9*size(T,1),1); % stores the value of the element contribution

ind = 1;

% basis functions P1
phi = {@(x) -x(1) - x(2) + 1;
          @(x) x(1);
          @(x) x(2) };

no_of_basis = size(phi,1);

% quadrature weights (for a 7-point Gauss quadrature)
w = [ 0.225; 0.132394152788506; 0.125939180544827 ];
% components for quadrature points
alpha = [0.059715871789770; 0.797426985353087 ];
beta  = [0.470142064105115; 0.101286507323456 ];
% quadrature points
x_quad = [1/3, 1/3;
                 beta(1), beta(1);
                 alpha(1),beta(1);
                 beta(1),alpha(1);
                 beta(2), beta(2);
                 alpha(2),beta(2);
                 beta(2),alpha(2)];

% evaluate basis functions in quadrature points
% (3 basis functions, 7 quadrature points)
phi_in_quad_point = zeros(3,7);
for i=1:no_of_basis % index basis function
   for q=1:7 % index quad point
       phi_in_quad_point(i,q) = phi{i}(x_quad(q,:));
   end
end
   
for k = 1:size(T,1)
    tri = T(k,:); % Lagrange node indices of the triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
   
    % construct u^n and d^n on reference element;
    % store u^n and d^n in quadrature points 
    % (sum up contributions of each basis function)
   un_in_quad_point = zeros(7);
   dn_in_quad_point = zeros(7);
   for q=1:7
      un_in_quad_point(q) = u(tri(1))*phi_in_quad_point(1,q) + u(tri(2))*phi_in_quad_point(2,q) + ... 
             u(tri(3))*phi_in_quad_point(3,q);
      dn_in_quad_point(q) = d(tri(1))*phi_in_quad_point(1,q) + d(tri(2))*phi_in_quad_point(2,q) + ... 
             d(tri(3))*phi_in_quad_point(3,q);
   end
%    
    %% assemble
   for i = 1:no_of_basis
       for j = 1:no_of_basis

           % Exploit Symmetry f(x)*phi{i}(x)*phi{j}(x)=f(x)*phi{j}(x)*phi{i}(x)
           if ( i <= j )

              % integrate on reference element
              % order 6 quadrature is sufficient for computing optimal tau

              % integrate   Re(u^n (d^n)*)phi_i phi_j 
              e_1 = w(1) * un_in_quad_point(1) * conj( dn_in_quad_point(1) ) * phi_in_quad_point(i,1) * phi_in_quad_point(j,1) + ... 
                        w(2) * un_in_quad_point(2) * conj( dn_in_quad_point(2) ) * phi_in_quad_point(i,2) * phi_in_quad_point(j,2) + ... 
                        w(2) * un_in_quad_point(3) * conj( dn_in_quad_point(3) ) * phi_in_quad_point(i,3) * phi_in_quad_point(j,3) + ...
                        w(2) * un_in_quad_point(4) * conj( dn_in_quad_point(4) ) * phi_in_quad_point(i,4) * phi_in_quad_point(j,4) + ... 
                        w(3) * un_in_quad_point(5) * conj( dn_in_quad_point(5) ) * phi_in_quad_point(i,5) * phi_in_quad_point(j,5) + ...
                        w(3) * un_in_quad_point(6) * conj( dn_in_quad_point(6) ) * phi_in_quad_point(i,6) * phi_in_quad_point(j,6) + ...
                        w(3) * un_in_quad_point(7) * conj( dn_in_quad_point(7) ) * phi_in_quad_point(i,7) * phi_in_quad_point(j,7);
              e_1 = 0.5 * detBT * real(e_1);

              % integrate   |d^n|^2 phi_i phi_j 
              e_2 = w(1) * dn_in_quad_point(1) * conj( dn_in_quad_point(1) ) * phi_in_quad_point(i,1) * phi_in_quad_point(j,1) + ... 
                        w(2) * dn_in_quad_point(2) * conj( dn_in_quad_point(2) ) * phi_in_quad_point(i,2) * phi_in_quad_point(j,2) + ...
                        w(2) * dn_in_quad_point(3) * conj( dn_in_quad_point(3) ) * phi_in_quad_point(i,3) * phi_in_quad_point(j,3) + ...  
                        w(2) * dn_in_quad_point(4) * conj( dn_in_quad_point(4) ) * phi_in_quad_point(i,4) * phi_in_quad_point(j,4) + ... 
                        w(3) * dn_in_quad_point(5) * conj( dn_in_quad_point(5) ) * phi_in_quad_point(i,5) * phi_in_quad_point(j,5) + ...
                        w(3) * dn_in_quad_point(6) * conj( dn_in_quad_point(6) ) * phi_in_quad_point(i,6) * phi_in_quad_point(j,6) + ...
                        w(3) * dn_in_quad_point(7) * conj( dn_in_quad_point(7) ) * phi_in_quad_point(i,7) * phi_in_quad_point(j,7);
              e_2 = 0.5 * detBT * real(e_2);

              F_un_dn_i(ind) = tri(i);  % get global index of local node i
              F_un_dn_j(ind) = tri(j);  % get global index of local node j
              F_dn_dn_i(ind) = tri(i);
              F_dn_dn_j(ind) = tri(j);

              % integrate on reference element
              F_un_dn_val(ind) = e_1;
              F_dn_dn_val(ind) = e_2;

              ind = ind + 1;

              if ( i  <  j )

                  F_un_dn_i(ind) = tri(j);  % get global index of local node i
                  F_un_dn_j(ind) = tri(i);  % get global index of local node j
                  F_dn_dn_i(ind) = tri(j);
                  F_dn_dn_j(ind) = tri(i);

                  % integrate on reference element
                  F_un_dn_val(ind) = e_1;
                  F_dn_dn_val(ind) = e_2;

                  ind = ind + 1;

              end

           end
       end
   end


    % 
end

F_un_dn = sparse(F_un_dn_i,F_un_dn_j,F_un_dn_val);
F_dn_dn = sparse(F_dn_dn_i,F_dn_dn_j,F_dn_dn_val);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F_un_dn,F_dn_dn] = assembleNonlinearMatrices_un_dn_P1_parallel(u,d,T,Nd)
% parallel implementation
dim = size(Nd,1);

% basis functions P1
phi = {@(x) -x(1) - x(2) + 1;
          @(x) x(1);
          @(x) x(2) };

no_of_basis = size(phi,1);

% quadrature weights (for a 7-point Gauss quadrature)
w = [ 0.225; 0.132394152788506; 0.125939180544827 ];
% components for quadrature points
alpha = [0.059715871789770; 0.797426985353087 ];
beta  = [0.470142064105115; 0.101286507323456 ];
% quadrature points
x_quad = [1/3, 1/3;
                 beta(1), beta(1);
                 alpha(1),beta(1);
                 beta(1),alpha(1);
                 beta(2), beta(2);
                 alpha(2),beta(2);
                 beta(2),alpha(2)];

% evaluate basis functions in quadrature points
% (3 basis functions, 7 quadrature points)
phi_in_quad_point = zeros(3,7);
for i=1:no_of_basis % index basis function
   for q=1:7 % index quad point
       phi_in_quad_point(i,q) = phi{i}(x_quad(q,:));
   end
end

spmd
    my_index_start = floor(size(T,1)*(spmdIndex - 1)/spmdSize + 1);
    my_index_end = floor(size(T,1)*(spmdIndex)/spmdSize);

F_un_dn_i = zeros(9*size(T,1),1); % stores the global i-index of the matrix
F_un_dn_j = zeros(9*size(T,1),1); % stores the global j-index of the matrix
F_un_dn_val = zeros(9*size(T,1),1); % stores the value of the element contribution

F_dn_dn_i = zeros(9*size(T,1),1); % stores the global i-index of the matrix
F_dn_dn_j = zeros(9*size(T,1),1); % stores the global j-index of the matrix
F_dn_dn_val = zeros(9*size(T,1),1); % stores the value of the element contribution

ind = 1;
   
for k = my_index_start:my_index_end
    tri = T(k,:); % Lagrange node indices of the triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
   
    % construct u^n and d^n on reference element;
    % store u^n and d^n in quadrature points 
    % (sum up contributions of each basis function)
   un_in_quad_point = zeros(7);
   dn_in_quad_point = zeros(7);
   for q=1:7
      un_in_quad_point(q) = u(tri(1))*phi_in_quad_point(1,q) + u(tri(2))*phi_in_quad_point(2,q) + ... 
             u(tri(3))*phi_in_quad_point(3,q);
      dn_in_quad_point(q) = d(tri(1))*phi_in_quad_point(1,q) + d(tri(2))*phi_in_quad_point(2,q) + ... 
             d(tri(3))*phi_in_quad_point(3,q);
   end
%    
    %% assemble
   for i = 1:no_of_basis
       for j = 1:no_of_basis

           % Exploit Symmetry f(x)*phi{i}(x)*phi{j}(x)=f(x)*phi{j}(x)*phi{i}(x)
           if ( i <= j )

              % integrate on reference element
              % order 6 quadrature is sufficient for computing optimal tau

              % integrate   Re(u^n (d^n)*)phi_i phi_j 
              e_1 = w(1) * un_in_quad_point(1) * conj( dn_in_quad_point(1) ) * phi_in_quad_point(i,1) * phi_in_quad_point(j,1) + ... 
                        w(2) * un_in_quad_point(2) * conj( dn_in_quad_point(2) ) * phi_in_quad_point(i,2) * phi_in_quad_point(j,2) + ... 
                        w(2) * un_in_quad_point(3) * conj( dn_in_quad_point(3) ) * phi_in_quad_point(i,3) * phi_in_quad_point(j,3) + ...
                        w(2) * un_in_quad_point(4) * conj( dn_in_quad_point(4) ) * phi_in_quad_point(i,4) * phi_in_quad_point(j,4) + ... 
                        w(3) * un_in_quad_point(5) * conj( dn_in_quad_point(5) ) * phi_in_quad_point(i,5) * phi_in_quad_point(j,5) + ...
                        w(3) * un_in_quad_point(6) * conj( dn_in_quad_point(6) ) * phi_in_quad_point(i,6) * phi_in_quad_point(j,6) + ...
                        w(3) * un_in_quad_point(7) * conj( dn_in_quad_point(7) ) * phi_in_quad_point(i,7) * phi_in_quad_point(j,7);
              e_1 = 0.5 * detBT * real(e_1);

              % integrate   |d^n|^2 phi_i phi_j 
              e_2 = w(1) * dn_in_quad_point(1) * conj( dn_in_quad_point(1) ) * phi_in_quad_point(i,1) * phi_in_quad_point(j,1) + ... 
                        w(2) * dn_in_quad_point(2) * conj( dn_in_quad_point(2) ) * phi_in_quad_point(i,2) * phi_in_quad_point(j,2) + ...
                        w(2) * dn_in_quad_point(3) * conj( dn_in_quad_point(3) ) * phi_in_quad_point(i,3) * phi_in_quad_point(j,3) + ...  
                        w(2) * dn_in_quad_point(4) * conj( dn_in_quad_point(4) ) * phi_in_quad_point(i,4) * phi_in_quad_point(j,4) + ... 
                        w(3) * dn_in_quad_point(5) * conj( dn_in_quad_point(5) ) * phi_in_quad_point(i,5) * phi_in_quad_point(j,5) + ...
                        w(3) * dn_in_quad_point(6) * conj( dn_in_quad_point(6) ) * phi_in_quad_point(i,6) * phi_in_quad_point(j,6) + ...
                        w(3) * dn_in_quad_point(7) * conj( dn_in_quad_point(7) ) * phi_in_quad_point(i,7) * phi_in_quad_point(j,7);
              e_2 = 0.5 * detBT * real(e_2);

              F_un_dn_i(ind) = tri(i);  % get global index of local node i
              F_un_dn_j(ind) = tri(j);  % get global index of local node j
              F_dn_dn_i(ind) = tri(i);
              F_dn_dn_j(ind) = tri(j);

              % integrate on reference element
              F_un_dn_val(ind) = e_1;
              F_dn_dn_val(ind) = e_2;

              ind = ind + 1;

              if ( i  <  j )

                  F_un_dn_i(ind) = tri(j);  % get global index of local node i
                  F_un_dn_j(ind) = tri(i);  % get global index of local node j
                  F_dn_dn_i(ind) = tri(j);
                  F_dn_dn_j(ind) = tri(i);

                  % integrate on reference element
                  F_un_dn_val(ind) = e_1;
                  F_dn_dn_val(ind) = e_2;

                  ind = ind + 1;

              end

           end
       end
   end


    % 
end

del_index = find(F_un_dn_i==0,1);
    if isempty(del_index)
        del_index = length(s_i)+1;
    end

F_un_dn_worker = sparse(F_un_dn_i(1:del_index-1),F_un_dn_j(1:del_index-1),F_un_dn_val(1:del_index-1),dim,dim);
F_dn_dn_worker = sparse(F_dn_dn_i(1:del_index-1),F_dn_dn_j(1:del_index-1),F_dn_dn_val(1:del_index-1),dim,dim);

end

% collect results from workers
F_un_dn = sparse(dim,dim);
F_dn_dn = sparse(dim,dim);
for j = 1:length(my_index_start)
    F_un_dn = F_un_dn + F_un_dn_worker{j};
    F_dn_dn = F_dn_dn + F_dn_dn_worker{j};
end

end

