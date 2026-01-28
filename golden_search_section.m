%%%%%%----------Golden search section--------%%%%%%
function [TAU,F_TAU] =golden_search_section(a,b,f)

 % a=-1;                            % start of interval
 % b=4;                            % end of interval
 % f = @(tau) tau.^2 -2*tau +3;      % define function

t0=(1+sqrt(5))/2;  
epsilon=10^(-10);                  
      
tau1=a+(b-a)*t0/(2*t0+1);             
tau2=a+(b-a)*(1+t0)/(2*t0+1);

f_a = f(a);  % computing values in tau1, tau2 and boundary points
f_b = f(b);
f_tau1=f(tau1);                    
f_tau2=f(tau2);

%plot(tau1,f_tau1,'ro') ; 
%hold on
%plot(tau2,f_tau2,'ro');

imax= 1000; 
k=0; 
while ((abs(b-a)>epsilon) && (k<imax))
    k=k+1;
    A = [f_tau1 f_tau2 f_a f_b];
    if (min(A) == f_a) || (min(A) == f_tau1)
       
        b=tau2;
        tau2=tau1;
        tau1=a+(b-a)*t0/(2*t0+1);
        
        f_b = f(b);
        f_tau1=f(tau1);
        f_tau2=f(tau2);
        
        %plot(tau2,f_tau2,'bx');
    else
        a=tau1;
        tau1=tau2;
        tau2=a+(b-a)*(1+t0)/(2*t0+1);
        
        f_a = f(a);
        f_tau1=f(tau1);
        f_tau2=f(tau2);
        
        %plot(tau2,f_tau2,'bx');
    end
    
    k=k+1;
end
% chooses minimum point
if min(A) == f_tau1
    %sprintf('tau_min=%f', tau1)
    %sprintf('f(tau_min)=%f ', f_tau1)
    %plot(tau1,f_tau1,'rx')
    TAU = tau1;
    F_TAU = f_tau1;
elseif min(A) == f_tau2
    %sprintf('tau_min=%f', tau2)
    %sprintf('f(tau_min)=%f ', f_tau2)
    %plot(tau2,f_tau2,'rx')
    TAU = tau2;
    F_TAU = f_tau2;
elseif min(A) == f_a
    TAU = a;
    F_TAU = f_a;
else
    TAU = b;
    F_TAU = f_b;
end


%TAU
%F_TAU


