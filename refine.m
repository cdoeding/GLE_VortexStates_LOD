function [T,P,P0,P1dg] = refine(T)
%REFINE Summary of this function goes here
%   Detailed explanation goes here

% Construct data structure
[np] = size(T.p,1);
nt = size(T.t,1);
[e,nmbe] = computeEdges(T);
ne = size(e,1);
d2p = sparse(e,e(:,[2,1]),(1:ne)'*[1 1],np,np);
% New nodes from the mid points of each edge
newnode = 0.5.*(T.p(e(:,1),:)+T.p(e(:,2),:));
P = sparse((1:ne)'*[1 1],e,.5,ne,np);
T.p = [T.p; newnode];
P = [speye(np);P];
emarker = (np+1:np+ne)';
p = [T.t,emarker(d2p(T.t(:,1)+np*(T.t(:,2)-1))),...
    emarker(d2p(T.t(:,2)+np*(T.t(:,3)-1))),...
    emarker(d2p(T.t(:,3)+np*(T.t(:,1)-1)))];
T.t = [p(:,[1,4,6]);p(:,[4,2,5]);p(:,[6,5,3]);p(:,[4,5,6])];
tnew2t = repmat(1:nt,1,4).';
E = speye(nt);
P1dg = [kron(E,[1 0 0;.5 .5 0;.5 0 .5]);kron(E,[.5 .5 0;0 1 0;0 .5 .5]);
    kron(E,[.5 0 .5;0 .5 .5;0 0 1]);kron(E,[.5 .5 0;0 .5 .5;.5 0 .5])];
Dnodes = false(np,1);
Dnodes(T.Dnodes) = true;
T.Dnodes = find([Dnodes;Dnodes(e(:,1)) & Dnodes(e(:,2)) & nmbe==1]);
P0 = sparse((1:length(tnew2t))',tnew2t,1,length(tnew2t),nt);
end

