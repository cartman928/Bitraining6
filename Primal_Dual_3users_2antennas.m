function [v11, v12, v13, v21, v22, v23, v31, v32, v33] = Primal_Dual(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, P, n0)
%Duality

%Constants
lambda = [0;0;0];%Staring Point
v = zeros(18,1);%Staring Point
Lam = diag([lambda(1) lambda(1) lambda(2) lambda(2) lambda(3) lambda(3)]);
k = [(g1')*H11 (g1')*H12 (g1')*H13 (g2')*H21 (g2')*H22 (g2')*H23 (g3')*H31 (g3')*H32 (g3')*H33];
A = [(H11')*g1*(g1')*H11 (H11')*g1*(g1')*H12 (H11')*g1*(g1')*H13;(H12')*g1*(g1')*H11 (H12')*g1*(g1')*H12 (H12')*g1*(g1')*H13;(H13')*g1*(g1')*H11 (H13')*g1*(g1')*H12 (H13')*g1*(g1')*H13];
B = [(H21')*g2*(g2')*H21 (H21')*g2*(g2')*H22 (H21')*g2*(g2')*H23;(H22')*g2*(g2')*H21 (H22')*g2*(g2')*H22 (H22')*g2*(g2')*H23;(H23')*g2*(g2')*H21 (H23')*g2*(g2')*H22 (H23')*g2*(g2')*H23];
C = [(H31')*g3*(g3')*H31 (H31')*g3*(g3')*H32 (H31')*g3*(g3')*H33;(H32')*g3*(g3')*H31 (H32')*g3*(g3')*H32 (H32')*g3*(g3')*H33;(H33')*g3*(g3')*H31 (H33')*g3*(g3')*H32 (H33')*g3*(g3')*H33];
ABC =  [A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C];
%ABC =  [A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C]+10^(-2)*eye(18);
ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam];
%ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam]+10^(-2)*eye(18);
L_g = zeros(1,18);
%P/D Update
for n = 1:10^(3)
stepsize = 5*10^(-2);    
            
            %Primal Upadate
            L_g = -2*k+2*(v')*ABCLam;
            v = v - stepsize*L_g';
            Lp = 3-k*v-(k*v)'+v'*ABC*v+n0*(g1'*g1+g2'*g2+g3'*g3);
            
            %Dual Update
          
            D_g = [norm(v(1:2))^2+norm(v(7:8))^2+norm(v(13:14))^2-P;norm(v(3:4))^2+norm(v(9:10))^2+norm(v(15:16))^2-P;norm(v(5:6))^2+norm(v(11:12))^2+norm(v(17:18))^2-P];
            lambda = lambda + stepsize*D_g;
            
            Lam = diag([lambda(1) lambda(1) lambda(2) lambda(2) lambda(3) lambda(3)]);
            ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam];
            %ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam]+10^(-2)*eye(18);
            Ld = 3-k*v-(k*v)'+v'*ABCLam*v+n0*(g1'*g1+g2'*g2+g3'*g3)-P*(lambda(1)+lambda(2)+lambda(3));
      
            [real(Lp) real(Ld) norm(v(1:2))^2+norm(v(7:8))^2+norm(v(13:14))^2 norm(v(3:4))^2+norm(v(9:10))^2+norm(v(15:16))^2 norm(v(5:6))^2+norm(v(11:12))^2+norm(v(17:18))^2];
          
            
            %Lambda > 0
            
            if(lambda(1) < 10^(-8))
                lambda(1) = 0;
            end
            if(lambda(2) < 10^(-8))
                lambda(2) = 0;
            end
            if(lambda(3) < 10^(-8))
                lambda(3) = 0;
            end
            %}
            
            
             %Lambda < 0
            %{
            if(lambda(1) > -10^(-6))
                lambda(1) = 0;
            end
            if(lambda(2) > -10^(-6))
                lambda(2) = 0;
            end
            if(lambda(3) > -10^(-6))
                lambda(3) = 0;
            end
            %}
             
            lambda;
end

%V=(k*[inv([A+B+C+Lam]) 0*eye(6) 0*eye(6);0*eye(6)   inv([A+B+C+Lam])   0*eye(6);0*eye(6) 0*eye(6) inv([A+B+C+Lam])])';

v11 = v(1:2);
v21 = v(3:4);
v31 = v(5:6);

v12 = v(7:8);
v22 = v(9:10);
v32 = v(11:12);

v13 = v(13:14);
v23 = v(15:16);
v33 = v(17:18);

(norm(v11)^2+norm(v12)^2+norm(v13)^2);


%Verification
%{
for o = -10^(5):10^5
v = [10^(-5)*o;zeros(17,1)];
L(o+10^5+1) = 3-k*v-(k*v)'+v'*ABCLam*v+n0*(g1'*g1+g2'*g2+g3'*g3)-P*(lambda(1)+lambda(2)+lambda(3));
end
o = 1:2*10^5+1;
plot(o,L(o))
%axis([1 2*10^5+1 -5 5])
%}

end

    
    
