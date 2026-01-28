function e = evaluateQuadrature(x,w)
%EVALUATEQUADRATURE Summary of this function goes here
%   Detailed explanation goes here

e = sum(x*diag(w),2);

end

