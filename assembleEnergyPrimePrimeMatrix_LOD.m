function Epp = assembleEnergyPrimePrimeMatrix_LOD(A,kappa,u,T,Nd,Phi)
%ASSEMBLEGLOBALNONLINEARMATRIX Summary of this function goes here
%   Detailed explanation goes here

%RePhi = real(Phi);
%ImPhi = imag(Phi);

% ARR Matrix with entries  ARR(j,i)
%  = \int_{\Omega} (1/kappa^2) \nabla \phi_i * \nabla \phi_j + (|A|^2 - 1 + 3 Re(u_h)^2 + Im(u_h)^2 ) \phi_i \phi_j
Arr = Phi*assembleARRMatrix(A,kappa,u,T,Nd)*Phi';
% AIR Matrix with entries AIR(j,i)
%  = \int_{\Omega} (1/kappa) \phi_j ( A * \nabla \phi_i) -  (1/kappa) \phi_i ( A * \nabla \phi_j)  + 2 * Re(u_h) Im(u_h) \phi_i \phi_j
Air = Phi*assembleAIRMatrix(A,kappa,u,T,Nd)*Phi';
% AII Matrix with entries AII(j,i)
%  = \int_{\Omega} (1/kappa^2) \nabla \phi_i * \nabla \phi_j + (|A|^2 - 1 + Re(u_h)^2 + 3 Im(u_h)^2 ) \phi_i \phi_j
Aii = Phi*assembleAIIMatrix(A,kappa,u,T,Nd)*Phi';

N = size(Arr,1);

Epp(1:N,1:N) = Arr;
Epp(1:N,N+1:2*N) = Air;
Epp(N+1:2*N,1:N) = transpose(Air);
Epp(N+1:2*N,N+1:2*N) = Aii;

end
 
% assemble entries  
% \int_{\Omega} (1/kappa^2) \nabla \phi_i * \nabla \phi_j + (|A|^2 - 1 + 3 Re(u_h)^2 + Im(u_h)^2 ) \phi_i \phi_j
function Arr = assembleARRMatrix(A,kappa,u,T,Nd)
% Summary of this function goes here
%   Detailed explanation goes here

disp("start assemble global ARR matrix")

Arr_i = zeros(9*size(T,1),1);
Arr_j = zeros(9*size(T,1),1);
Arr_val = zeros(9*size(T,1),1);
ind = 1;

grad1 = [-1; -1];
grad2 = [1; 0];
grad3 = [0; 1];

phi1 = @(x) -x(1) - x(2) + 1;
phi2 = @(x) x(1);
phi3 = @(x) x(2);

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node
    
    %% transformation to reference triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    b = [z1(1); z1(2)];
    
    BTinv = inv(BT)';
    D = (BTinv')*BTinv; % musste reihenfolge tauschen
    
    %% construct f on triangle  hh h  h b
    f1r = @(x) real(u(tri(1)))*phi1(x); % real contribution of phi1
    f2r = @(x) real(u(tri(2)))*phi2(x); % real contribution of phi2
    f3r = @(x) real(u(tri(3)))*phi3(x); % real contribution of phi3
    f1i = @(x) imag(u(tri(1)))*phi1(x); % imag contribution of phi1
    f2i = @(x) imag(u(tri(2)))*phi2(x); % imag contribution of phi2
    f3i = @(x) imag(u(tri(3)))*phi3(x); % imag contribution of phi3
    
    f = @(x) 3*(f1r(x) + f2r(x) + f3r(x))^2 + (f1i(x) + f2i(x) + f3i(x))^2 - 1 + norm(A(BT*x+b),2)^2;
    
    %% assemble
    %grad phi1 grad phi2
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad2'*D*grad1 + f(x)*phi2(x)*phi1(x),6); 

    Arr_i(ind) = tri(1);
    Arr_j(ind) = tri(2);
    Arr_val(ind) = e;
    ind = ind + 1;

    Arr_i(ind) = tri(2);
    Arr_j(ind) = tri(1);
    Arr_val(ind) = e;
    ind = ind + 1;
    
    %grad phi2 grad phi3
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad3'*D*grad2 + f(x)*phi3(x)*phi2(x),6);
    
    Arr_i(ind) = tri(2);
    Arr_j(ind) = tri(3);
    Arr_val(ind) = e;
    ind = ind + 1;

    Arr_i(ind) = tri(3);
    Arr_j(ind) = tri(2);
    Arr_val(ind) = e;
    ind = ind + 1;
    
    %grad phi3 grad phi1
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad1'*D*grad3 + f(x)*phi1(x)*phi3(x),6);
    
    Arr_i(ind) = tri(3);
    Arr_j(ind) = tri(1);
    Arr_val(ind) = e;
    ind = ind + 1;

    Arr_i(ind) = tri(1);
    Arr_j(ind) = tri(3);
    Arr_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi1 grad phi1
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad1'*D*grad1 + f(x)*phi1(x)*phi1(x),6); 
    
    Arr_i(ind) = tri(1);
    Arr_j(ind) = tri(1);
    Arr_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi2 grad phi2
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad2'*D*grad2 + f(x)*phi2(x)*phi2(x),6); 
    
    Arr_i(ind) = tri(2);
    Arr_j(ind) = tri(2);
    Arr_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi3 grad phi3
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad3'*D*grad3 + f(x)*phi3(x)*phi3(x),6); 
    
    Arr_i(ind) = tri(3);
    Arr_j(ind) = tri(3);
    Arr_val(ind) = e;
    ind = ind + 1;
    
end

Arr = sparse(Arr_i,Arr_j,Arr_val);

end

% assemble entries  
% \int_{\Omega} (1/kappa) \phi_j ( A * \nabla \phi_i) -  (1/kappa) \phi_k ( A * \nabla \phi_j)  + 2 * Re(u_h) Im(u_h) \phi_i \phi_j
function Air = assembleAIRMatrix(A,kappa,u,T,Nd)
% Summary of this function goes here
%   Detailed explanation goes here

disp("start assemble global AIR matrix")

Air_i = zeros(9*size(T,1),1);
Air_j = zeros(9*size(T,1),1);
Air_val = zeros(9*size(T,1),1);
ind = 1;

grad1 = [-1; -1];
grad2 = [1; 0];
grad3 = [0; 1];

phi1 = @(x) -x(1) - x(2) + 1;
phi2 = @(x) x(1);
phi3 = @(x) x(2);

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    b = [z1(1); z1(2)];
    
    BTinv = inv(BT)';
    D = BTinv*(BTinv');
    
    %% construct f on triangle
    f1r = @(x) real(u(tri(1)))*phi1(x); % real contribution of phi1
    f2r = @(x) real(u(tri(2)))*phi2(x); % real contribution of phi2
    f3r = @(x) real(u(tri(3)))*phi3(x); % real contribution of phi3
    f1i = @(x) imag(u(tri(1)))*phi1(x); % imag contribution of phi1
    f2i = @(x) imag(u(tri(2)))*phi2(x); % imag contribution of phi2
    f3i = @(x) imag(u(tri(3)))*phi3(x); % imag contribution of phi3
    
    f = @(x) 2*(f1r(x) + f2r(x) + f3r(x)) * (f1i(x) + f2i(x) + f3i(x));
    
    %% assemble
    %grad phi1 grad phi2
    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi1(x)*A(BT*x+b))'*BTinv*grad2 - (phi2(x)*A(BT*x+b))'*BTinv*grad1) + f(x)*phi2(x)*phi1(x),6); 
    Air_i(ind) = tri(1);
    Air_j(ind) = tri(2);
    Air_val(ind) = e;
    ind = ind + 1;

    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi2(x)*A(BT*x+b))'*BTinv*grad1 - (phi1(x)*A(BT*x+b))'*BTinv*grad2) + f(x)*phi1(x)*phi2(x),6); 
    Air_i(ind) = tri(2);
    Air_j(ind) = tri(1);
    Air_val(ind) = e;
    ind = ind + 1;
    
    %grad phi2 grad phi3
    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi2(x)*A(BT*x+b))'*BTinv*grad3 - (phi3(x)*A(BT*x+b))'*BTinv*grad2) + f(x)*phi3(x)*phi2(x),6); 
    Air_i(ind) = tri(2);
    Air_j(ind) = tri(3);
    Air_val(ind) = e;
    ind = ind + 1;
    
    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi3(x)*A(BT*x+b))'*BTinv*grad2 - (phi2(x)*A(BT*x+b))'*BTinv*grad3) + f(x)*phi2(x)*phi3(x),6); 
    Air_i(ind) = tri(3);
    Air_j(ind) = tri(2);
    Air_val(ind) = e;
    ind = ind + 1;
    
    %grad phi3 grad phi1
    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi3(x)*A(BT*x+b))'*BTinv*grad1 - (phi1(x)*A(BT*x+b))'*BTinv*grad3) + f(x)*phi1(x)*phi3(x),6); 
    Air_i(ind) = tri(3);
    Air_j(ind) = tri(1);
    Air_val(ind) = e;
    ind = ind + 1;
    
    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi1(x)*A(BT*x+b))'*BTinv*grad3 - (phi3(x)*A(BT*x+b))'*BTinv*grad1) + f(x)*phi3(x)*phi1(x),6); 
    Air_i(ind) = tri(1);
    Air_j(ind) = tri(3);
    Air_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi1 grad phi1
    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi1(x)*A(BT*x+b))'*BTinv*grad1 - (phi1(x)*A(BT*x+b))'*BTinv*grad1) + f(x)*phi1(x)*phi1(x),6); 
    Air_i(ind) = tri(1);
    Air_j(ind) = tri(1);
    Air_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi2 grad phi2
    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi2(x)*A(BT*x+b))'*BTinv*grad2 - (phi2(x)*A(BT*x+b))'*BTinv*grad2) + f(x)*phi2(x)*phi2(x),6); 
    Air_i(ind) = tri(2);
    Air_j(ind) = tri(2);
    Air_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi3 grad phi3
    e = detBT*integrate_unit_triangle(@(x) 1/kappa*((phi3(x)*A(BT*x+b))'*BTinv*grad3 - (phi3(x)*A(BT*x+b))'*BTinv*grad3) + f(x)*phi3(x)*phi3(x),6); 
    Air_i(ind) = tri(3);
    Air_j(ind) = tri(3);
    Air_val(ind) = e;
    ind = ind + 1;
    
end

Air = sparse(Air_i,Air_j,Air_val);

end

% assemble entries  
% \int_{\Omega} (1/kappa^2) \nabla \phi_i * \nabla \phi_j + (|A|^2 - 1 + Re(u_h)^2 + 3 Im(u_h)^2 ) \phi_i \phi_j
function Aii = assembleAIIMatrix(A,kappa,u,T,Nd)
% Summary of this function goes here
%   Detailed explanation goes here

disp("start assemble global AII matrix")

Aii_i = zeros(9*size(T,1),1);
Aii_j = zeros(9*size(T,1),1);
Aii_val = zeros(9*size(T,1),1);
ind = 1;

grad1 = [-1; -1];
grad2 = [1; 0];
grad3 = [0; 1];

phi1 = @(x) -x(1) - x(2) + 1;
phi2 = @(x) x(1);
phi3 = @(x) x(2);

for k = 1:size(T,1)
    tri = T(k,:); %node index of triangle
    z1 = Nd(T(k,1),:); %coordinates of 1st triangle node
    z2 = Nd(T(k,2),:); %coordinates of 2nd triangle node
    z3 = Nd(T(k,3),:); %coordinates of 3rd triangle node
    
    %% transformation to refenrence triangle
    BT = [z2(1)-z1(1), z3(1)-z1(1); ...
        z2(2)-z1(2), z3(2)-z1(2)];
    
    detBT = BT(1,1)*BT(2,2)-BT(1,2)*BT(2,1);
    b = [z1(1); z1(2)];
    
    BTinv = inv(BT)';
    D = (BTinv')*BTinv; % musste reihenfolge tauschen
    
    %% construct f on triangle
    f1r = @(x) real(u(tri(1)))*phi1(x); % real contribution of phi1
    f2r = @(x) real(u(tri(2)))*phi2(x); % real contribution of phi2
    f3r = @(x) real(u(tri(3)))*phi3(x); % real contribution of phi3
    f1i = @(x) imag(u(tri(1)))*phi1(x); % imag contribution of phi1
    f2i = @(x) imag(u(tri(2)))*phi2(x); % imag contribution of phi2
    f3i = @(x) imag(u(tri(3)))*phi3(x); % imag contribution of phi3
    
    f = @(x) (f1r(x) + f2r(x) + f3r(x))^2 + 3*(f1i(x) + f2i(x) + f3i(x))^2 - 1 + norm(A(BT*x+b),2)^2;
    
    %% assemble
    %grad phi1 grad phi2
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad1'*D*grad2 + f(x)*phi1(x)*phi2(x),6);
    
    Aii_i(ind) = tri(1);
    Aii_j(ind) = tri(2);
    Aii_val(ind) = e;
    ind = ind + 1;

    Aii_i(ind) = tri(2);
    Aii_j(ind) = tri(1);
    Aii_val(ind) = e;
    ind = ind + 1;
    
    %grad phi2 grad phi3
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad2'*D*grad3 + f(x)*phi2(x)*phi3(x),6);
    
    Aii_i(ind) = tri(2);
    Aii_j(ind) = tri(3);
    Aii_val(ind) = e;
    ind = ind + 1;
    
    Aii_i(ind) = tri(3);
    Aii_j(ind) = tri(2);
    Aii_val(ind) = e;
    ind = ind + 1;
    
    %grad phi3 grad phi1
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad3'*D*grad1 + f(x)*phi3(x)*phi1(x),6); 
    
    Aii_i(ind) = tri(3);
    Aii_j(ind) = tri(1);
    Aii_val(ind) = e;
    ind = ind + 1;
    
    Aii_i(ind) = tri(1);
    Aii_j(ind) = tri(3);
    Aii_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi1 grad phi1
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad1'*D*grad1 + f(x)*phi1(x)*phi1(x),6); 
    
    Aii_i(ind) = tri(1);
    Aii_j(ind) = tri(1);
    Aii_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi2 grad phi2
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad2'*D*grad2 + f(x)*phi2(x)*phi2(x),6);
    
    Aii_i(ind) = tri(2);
    Aii_j(ind) = tri(2);
    Aii_val(ind) = e;
    ind = ind + 1;
    
    
    %grad phi3 grad phi3
    e = detBT*integrate_unit_triangle(@(x) 1/(kappa^2)*grad3'*D*grad3 + f(x)*phi3(x)*phi3(x),6);
    
    Aii_i(ind) = tri(3);
    Aii_j(ind) = tri(3);
    Aii_val(ind) = e;
    ind = ind + 1;
    
end

Aii = sparse(Aii_i,Aii_j,Aii_val);

end