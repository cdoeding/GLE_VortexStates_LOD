function [x,w] = getQuadrature(p)
%GETQUADRATURE Summary of this function goes here
%   Detailed explanation goes here

if p == 1
    x = [1/3, 1/3]';
    w = 0.5;

elseif p == 3

    x = [0 0.5 0.5; 0.5 0 0.5];
    w = [1/6 1/6 1/6];

elseif p == 7
    w = 0.5*[0.225;
        0.132394152788506;
        0.132394152788506;
        0.132394152788506;
        0.125939180544827;
        0.125939180544827;
        0.125939180544827 ];

    alpha = [0.059715871789770;
        0.797426985353087 ];
    beta  = [0.470142064105115;
        0.101286507323456 ];

    x = [1/3, 1/3;
        beta(1), beta(1);
        alpha(1),beta(1);
        beta(1),alpha(1);
        beta(2), beta(2);
        alpha(2),beta(2);
        beta(2),alpha(2)]';

end
end

