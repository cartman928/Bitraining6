function [v11, v12, v21, v22] = Primal_Dual_2users_2antennas(H11, H12, H21, H22, g1, g2, P, n0)
%Duality

%Constants
lambda = [0;0];%Staring Point
v = zeros(16,1);%Staring Point
Lam = diag([lambda(1) lambda(1) lambda(1) lambda(1) lambda(2) lambda(2) lambda(2) lambda(2)]);
k = [(g1')*H11 (g1')*H12 (g2')*H21 (g2')*H22];
A = [(H11')*g1*(g1')*H11 (H11')*g1*(g1')*H12;(H12')*g1*(g1')*H11 (H12')*g1*(g1')*H12];
B = [(H21')*g2*(g2')*H21 (H21')*g2*(g2')*H22;(H22')*g2*(g2')*H21 (H22')*g2*(g2')*H22];
AB =  [A+B 0*eye(8);0*eye(8) A+B];
%ABC =  [A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C]+10^(-2)*eye(18);
ABLam =  [A+B+Lam 0*eye(8);0*eye(8) A+B+Lam];
%ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam]+10^(-2)*eye(18);
L_g = zeros(1,16);
%P/D Update
for n = 1:10^(3)
stepsize = 5*10^(-2);    
            
            %Primal Upadate
            L_g = -2*k+2*(v')*ABLam;
            v = v - stepsize*L_g';
            Lp = 2-k*v-(k*v)'+v'*AB*v+n0*(g1'*g1+g2'*g2);
            
            %Dual Update
          
            D_g = [norm(v(1:4))^2+norm(v(9:12))^2-P;norm(v(5:8))^2+norm(v(13:16))^2];
            lambda = lambda + stepsize*D_g;
            
            Lam = diag([lambda(1) lambda(1) lambda(1) lambda(1) lambda(2) lambda(2) lambda(2) lambda(2)]);
            ABLam =  [A+B+Lam 0*eye(8);0*eye(8) A+B+Lam];
            %ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam]+10^(-2)*eye(18);
            Ld = 2-k*v-(k*v)'+v'*ABLam*v+n0*(g1'*g1+g2'*g2)-P*(lambda(1)+lambda(2));
      
            [real(Lp) real(Ld) norm(v(1:4))^2+norm(v(9:12))^2 norm(v(5:8))^2+norm(v(13:16))^2];
          
            
            %Lambda > 0
            
            if(lambda(1) < 10^(-8))
                lambda(1) = 0;
            end
            if(lambda(2) < 10^(-8))
                lambda(2) = 0;
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

v11 = v(1:4);
v21 = v(5:8);

v12 = v(9:12);
v22 = v(13:16);



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

    
    
