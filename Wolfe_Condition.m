P=1;
lambda = -5;
x = -5;
c1 = 10^(-4);
c2 = 0.1;
for n = 1:10^6

            %Primal Upadate
            for z = 1:10^5
                
                    %Armijo Rule(Sufficient Decrease)
                    for k = 1:10^5
                        alpha_k = k*10^(-3);
                        if (    (x+alpha_k*((-2*x)/norm(2*x)))^2+lambda*(x+alpha_k*((-2*x)/norm(2*x))+P) <= (x)^2+lambda*(x+P) + c1*alpha_k*((-2*x)/norm(2*x))*(2*x+lambda) )
                            alpha = alpha_k;
                            %Curvature Condition
                            [(x+alpha*((-2*x)/norm(2*x)))^2+lambda*(x+alpha*((-2*x)/norm(2*x))+P) c2*(2*x+lambda)];
                            
                            if( (x+alpha*((-2*x)/norm(2*x)))^2+lambda*(x+alpha*((-2*x)/norm(2*x))+P) >= c2*(2*x+lambda))
                            x = x+alpha*(-2*x)/norm(2*x)
                            disp('converges!');
                            break
                            end
                            
                            
                        end
                    end
            
            end
                   
                        
       a= 1;             
                  
                    for k = 1:10^5
                    L_g = 2*x+lambda;
                    x = x - alpha*L_g;
                    end
            L_p = x^2+lambda*(x+P);
            %Dual Update
            stepsize = 10^(-3);
            D_g = x+P;
            lambda = lambda + stepsize*D_g;
            L_d = x^2+lambda*(x+P);
            [L_p-L_d x]

end
%plot(lambda)