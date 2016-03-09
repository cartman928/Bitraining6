P=1;
lambda = 0;
A = 1;
k = 1+4i;
x = 0;
for n = 1:10^6
stepsize = 10^(-3);
            %Primal Upadate
            for o = 1:10^5
            L_g = (x')*A-k+lambda*(x');
            x = x - stepsize*L_g'
            end
            L_p = (x')*A*x-k*x-(k*x)'+lambda*((norm(x))^2-P);
            %Dual Update
            D_g = (norm(x))^2-P;
            lambda = lambda + stepsize*D_g;
            L_d = (x')*A*x-k*x-(k*x)'+lambda*((norm(x))^2-P);
            [L_p L_d x]
end
%plot(lambda)