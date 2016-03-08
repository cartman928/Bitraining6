P=1;
lambda = 5;
x = 3;
for n = 1:10^6
stepsize = 10^(-3);
            %Primal Upadate
            for k = 1:10^5
            L_g = 2*x+lambda;
            x = x - stepsize*L_g;
            end
            %Dual Update
            D_g = x+P;
            lambda = lambda + stepsize*D_g       
end
%plot(lambda)