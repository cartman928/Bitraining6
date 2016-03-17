P=1;
lambda = -5;
x = -3;
for n = 1:10^5
stepsize = 10^(-1);
            %Primal Upadate
            L_g = 2*x+lambda;
            x = x - stepsize*L_g;
            %end
            L_p = x^2;
            %Dual Update
            D_g = x+P;
            lambda = lambda + stepsize*D_g;
            L_d = x^2+lambda*(x+P);
            [L_p-L_d x]
end
%plot(lambda)