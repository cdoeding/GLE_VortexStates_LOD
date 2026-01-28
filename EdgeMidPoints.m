function [xmid,ymid] = EdgeMidPoints(p,t2e,t)
%EDGEMIDPOINTS Summary of this function goes here
%   Detailed explanation goes here

i=t(1,:); 
j=t(2,:); 
k=t(3,:); % triangle vertices 
t2e=t2e(:); % all edges in a long row
start=[j i i]; % start vertices of all edges
stop =[k k j]; % stop
xmid=(p(1,start)+p(1,stop))/2; % mid point x-coordinates 
ymid=(p(2,start)+p(2,stop))/2; % y- 
[~,idx]=unique(t2e); % remove duplicate edges 
xmid=xmid(idx); % unique edge x-coordinates 
ymid=ymid(idx); % unique edge y-coordinates 

end

