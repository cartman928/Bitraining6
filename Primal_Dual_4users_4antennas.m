function [v11, v12, v13, v14, v21, v22, v23, v24, v31, v32, v33, v34, v41, v42, v43, v44] = Primal_Dual_4users(H11, H12, H13, H14, H21, H22, H23, H24, H31, H32, H33, H34, H41, H42, H43, H44, g1, g2, g3, g4, P, n0)
%Duality

%Constants
lambda = [0;0;0;0];%Staring Point
v = zeros(64,1);%Staring Point
Lam = diag([lambda(1) lambda(1) lambda(1) lambda(1) lambda(2) lambda(2) lambda(2) lambda(2) lambda(3) lambda(3) lambda(3) lambda(3) lambda(4) lambda(4) lambda(4) lambda(4)]);
k = [(g1')*H11 (g1')*H12 (g1')*H13 (g1')*H14 (g2')*H21 (g2')*H22 (g2')*H23 (g2')*H24 (g3')*H31 (g3')*H32 (g3')*H33 (g3')*H34 (g4')*H41 (g4')*H42 (g4')*H43 (g4')*H44];

A = [(H11')*g1*(g1')*H11 (H11')*g1*(g1')*H12 (H11')*g1*(g1')*H13 (H11')*g1*(g1')*H14;(H12')*g1*(g1')*H11 (H12')*g1*(g1')*H12 (H12')*g1*(g1')*H13 (H12')*g1*(g1')*H14;(H13')*g1*(g1')*H11 (H13')*g1*(g1')*H12 (H13')*g1*(g1')*H13 (H13')*g1*(g1')*H14;(H14')*g1*(g1')*H11 (H14')*g1*(g1')*H12 (H14')*g1*(g1')*H13 (H14')*g1*(g1')*H14];

B = [(H21')*g2*(g2')*H21 (H21')*g2*(g2')*H22 (H21')*g2*(g2')*H23 (H21')*g2*(g2')*H24;(H22')*g2*(g2')*H21 (H22')*g2*(g2')*H22 (H22')*g2*(g2')*H23 (H22')*g2*(g2')*H24;(H23')*g2*(g2')*H21 (H23')*g2*(g2')*H22 (H23')*g2*(g2')*H23 (H23')*g2*(g2')*H24;(H24')*g2*(g2')*H21 (H24')*g2*(g2')*H22 (H24')*g2*(g2')*H23 (H24')*g2*(g2')*H24];

C = [(H31')*g3*(g3')*H31 (H31')*g3*(g3')*H32 (H31')*g3*(g3')*H33 (H31')*g3*(g3')*H34;(H32')*g3*(g3')*H31 (H32')*g3*(g3')*H32 (H32')*g3*(g3')*H33 (H32')*g3*(g3')*H34;(H33')*g3*(g3')*H31 (H33')*g3*(g3')*H32 (H33')*g3*(g3')*H33 (H33')*g3*(g3')*H34;(H34')*g3*(g3')*H31 (H34')*g3*(g3')*H32 (H34')*g3*(g3')*H33 (H34')*g3*(g3')*H34];

D = [(H41')*g4*(g4')*H41 (H41')*g4*(g4')*H42 (H41')*g4*(g4')*H43 (H41')*g4*(g4')*H44;(H42')*g4*(g4')*H41 (H42')*g4*(g4')*H42 (H42')*g4*(g4')*H43 (H42')*g4*(g4')*H44;(H43')*g4*(g4')*H41 (H43')*g4*(g4')*H42 (H43')*g4*(g4')*H43 (H43')*g4*(g4')*H44;(H44')*g4*(g4')*H41 (H44')*g4*(g4')*H42 (H44')*g4*(g4')*H43 (H44')*g4*(g4')*H44];

ABCD =  [A+B+C+D 0*eye(16) 0*eye(16) 0*eye(16);0*eye(16) A+B+C+D 0*eye(16) 0*eye(16);0*eye(16) 0*eye(16) A+B+C+D 0*eye(16);0*eye(16) 0*eye(16) 0*eye(16) A+B+C+D];
%ABC =  [A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C]+10^(-2)*eye(18);
ABCDLam =  [A+B+C+D+Lam 0*eye(16) 0*eye(16) 0*eye(16);0*eye(16) A+B+C+D+Lam 0*eye(16) 0*eye(16);0*eye(16) 0*eye(16)  A+B+C+D+Lam 0*eye(16);0*eye(16) 0*eye(16) 0*eye(16) A+B+C+D+Lam];
%ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam]+10^(-2)*eye(18);
L_g = zeros(1,64);

%P/D Update
for n = 1:5*10^(3)
stepsize = 5*10^(-2);    
            
            %Primal Upadate
            L_g = -2*k+2*(v')*ABCDLam;
            v = v - stepsize*L_g';
            Lp = 4-k*v-(k*v)'+v'*ABCD*v+n0*(g1'*g1+g2'*g2+g3'*g3+g4'*g4);
            
            %Dual Update
          
            D_g = [norm(v(1:4))^2+norm(v(17:20))^2+norm(v(33:36))^2+norm(v(49:52))^2-P;norm(v(5:8))^2+norm(v(21:24))^2+norm(v(37:40))^2+norm(v(53:56))^2-P;norm(v(9:12))^2+norm(v(25:28))^2+norm(v(41:44))^2+norm(v(57:60))^2-P;norm(v(13:16))^2+norm(v(29:32))^2+norm(v(45:48))^2+norm(v(61:64))^2-P];
            lambda = lambda + stepsize*D_g;
            
            Lam = diag([lambda(1) lambda(1) lambda(1) lambda(1) lambda(2) lambda(2) lambda(2) lambda(2) lambda(3) lambda(3) lambda(3) lambda(3) lambda(4) lambda(4) lambda(4) lambda(4)]);
            ABCDLam =  [A+B+C+D+Lam 0*eye(16) 0*eye(16) 0*eye(16);0*eye(16) A+B+C+D+Lam 0*eye(16) 0*eye(16);0*eye(16) 0*eye(16)  A+B+C+D+Lam 0*eye(16);0*eye(16) 0*eye(16) 0*eye(16) A+B+C+D+Lam];
            %ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam]+10^(-2)*eye(18);
            Ld = 4-k*v-(k*v)'+v'*ABCDLam*v+n0*(g1'*g1+g2'*g2+g3'*g3+g4'*g4)-P*(lambda(1)+lambda(2)+lambda(3)+lambda(4));
      
            [real(Lp) real(Ld) norm(v(1:4))^2+norm(v(17:20))^2+norm(v(33:36))^2+norm(v(49:52))^2 norm(v(5:8))^2+norm(v(21:24))^2+norm(v(37:40))^2+norm(v(53:56))^2 norm(v(9:12))^2+norm(v(25:28))^2+norm(v(41:44))^2+norm(v(57:60))^2 norm(v(13:16))^2+norm(v(29:32))^2+norm(v(45:48))^2+norm(v(61:64))^2]
          
            
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
            
            if(lambda(4) < 10^(-8))
                lambda(4) = 0;
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
v31 = v(9:12);
v41 = v(13:16);

v12 = v(17:20);
v22 = v(21:24);
v32 = v(25:28);
v42 = v(29:32);

v13 = v(33:36);
v23 = v(37:40);
v33 = v(41:44);
v43 = v(45:48);

v14 = v(49:52);
v24 = v(53:56);
v34 = v(57:60);
v44 = v(61:64);

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

    
    
