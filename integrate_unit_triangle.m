function out = integrate_unit_triangle(f,p)

if p == 2
    out = 0.5*( f([1/3, 1/3]') );

elseif p == 3

    out = (0.5/3)*( f([0,0.5]') + f([0.5,0]') + f([0.5,0.5]') );

elseif p == 6
    w = [0.225;
        0.132394152788506;
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
        beta(2),alpha(2)];

    out = w(1)*f(x(1,:)') ...
        + w(2)*f(x(2,:)') ...
        + w(2)*f(x(3,:)') ...
        + w(2)*f(x(4,:)') ...
        + w(3)*f(x(5,:)') ...
        + w(3)*f(x(6,:)') ...
        + w(3)*f(x(7,:)');

    out = 0.5*out;

elseif p == 9
    out = 0;
    [x,wx] = get_gauss_nodes(5,0,1);
    [t,wt] = get_gauss_nodes(3,0,1);

    for i = 1:size(wx,2)
        for j = 1:size(wt,2)
            out = out + wx(i)*wt(j)*(1-x(i))*f([x(i), (1-x(i))*t(j)]');
        end
    end

end
end