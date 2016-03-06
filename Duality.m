function [v11, v12, v13, v21, v22, v23, v31, v32, v33] = Duality(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, P, n0)
%Duality

%Constants
k = [g1'*H11 g1'*H12 g1'*H13 g2'*H21 g2'*H22 g2'*H23 g3'*H31 g3'*H32 g3'*H33];
A = [H11'*g1*g1'*H11 H11'*g1*g1'*H12 H11'*g1*g1'*H13;H12'*g1*g1'*H11 H12'*g1*g1'*H12 H12'*g1*g1'*H13;H13'*g1*g1'*H11 H13'*g1*g1'*H12 H13'*g1*g1'*H13];
B = [H21'*g2*g2'*H21 H21'*g2*g2'*H22 H21'*g2*g2'*H23;H22'*g2*g2'*H21 H22'*g2*g2'*H22 H22'*g2*g2'*H23;H23'*g2*g2'*H21 H23'*g2*g2'*H22 H23'*g2*g2'*H23];
C = [H31'*g3*g3'*H31 H31'*g3*g3'*H32 H31'*g3*g3'*H33;H32'*g3*g3'*H31 H32'*g3*g3'*H32 H32'*g3*g3'*H33;H33'*g3*g3'*H31 H33'*g3*g3'*H32 H33'*g3*g3'*H33];

%Dual Problem
%{
for o = -10^(4):10^4
lambda1 = o*10^(-3);
lambda2 = 2;
lambda3 = 2;
Lam = diag([lambda1 lambda1 lambda2 lambda2 lambda3 lambda3]);

g(o+10^(4)+1) = 3 - k*inv([A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam])*k'-P*(lambda1+lambda2+lambda3)+n0*(g1'*g1+g2'*g2+g3'*g3);
end
o = 1+0.5*10^4:(3/2)*10^4+1;
plot(o,g(o))
axis([1 2*10^4+1 -10^1 10^1])
%}


%Eigenvalues of ABC Matrix
E = sort((eig([A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C])));
%Least Eigenvalues
E(10);

%Gradients
%G = -P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6)));

%Search for Optimal Point
%stepsize = 1/real(E(18));
stepsize = 10^(-3);
lambda1 = real(-E(10)/2);
lambda2 = real(-E(10)/2);
lambda3 = real(-E(10)/2);

for n = 1:1*10^(5)
%{
lambda = real(lambda + stepsize*(-P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6)))))
real(-P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6))))
    %}
Lam = diag([lambda1 lambda1 lambda2 lambda2 lambda3 lambda3]);
gradient = [inv([A+B+C+Lam]) 0*eye(6) 0*eye(6);0*eye(6)   inv([A+B+C+Lam])   0*eye(6);0*eye(6) 0*eye(6) inv([A+B+C+Lam])]...
             *k'*k...
             *[inv([A+B+C+Lam]) 0*eye(6) 0*eye(6);0*eye(6)   inv([A+B+C+Lam])   0*eye(6);0*eye(6) 0*eye(6) inv([A+B+C+Lam])];
gradient1 = -P+gradient(1,1)+gradient(2,2)+gradient(7,7)+gradient(8,8)+gradient(13,13)+gradient(14,14);
lambda1 = real(lambda1 + stepsize*gradient1);

gradient2 = -P+gradient(3,3)+gradient(4,4)+gradient(9,9)+gradient(10,10)+gradient(15,15)+gradient(16,16);
lambda2 = real(lambda2 + stepsize*gradient2);

gradient3 = -P+gradient(5,5)+gradient(6,6)+gradient(11,11)+gradient(12,12)+gradient(17,17)+gradient(18,18);
lambda3 = real(lambda3 + stepsize*gradient3);

GG=norm(gradient1)+norm(gradient2)+norm(gradient3)

if((norm(gradient1)+norm(gradient2)+norm(gradient3))< 10^(-6)  )
            disp('converges!');
            break;%break the for loop if it's true the condition
end




end
GG=norm(gradient1)+norm(gradient2)+norm(gradient3)
%{
lambda1 = 10^(-8);
lambda2 = 10^(-8);
lambda3 = 10^(-8);
%}
Lam = diag([lambda1 lambda1 lambda2 lambda2 lambda3 lambda3]);

V=(k*[inv([A+B+C+Lam]) 0*eye(6) 0*eye(6);0*eye(6)   inv([A+B+C+Lam])   0*eye(6);0*eye(6) 0*eye(6) inv([A+B+C+Lam])])';

v11 = V(1:2);
v21 = V(3:4);
v31 = V(5:6);

v12 = V(1:2);
v22 = V(3:4);
v32 = V(5:6);

v13 = V(1:2);
v23 = V(3:4);
v33 = V(5:6);

(norm(v11)^2+norm(v12)^2+norm(v13)^2)




%Duality Gap
v = [v11;v12;v13;v21;v22;v23;v31;v32;v33];
g = 3 - k*inv([A+B+C+lambda1*eye(6) 0*eye(6) 0*eye(6);0*eye(6) A+B+C+lambda2*eye(6) 0*eye(6);0*eye(6) 0*eye(6) A+B+C+lambda3*eye(6)])*k'-P*(lambda1+lambda2+lambda3)+norm([g1;g2;g3])^2*n0;
L = 3-k*v-(k*v)'+v'*[A+B+C+lambda1*eye(6) 0*eye(6) 0*eye(6);0*eye(6) A+B+C+lambda2*eye(6) 0*eye(6);0*eye(6) 0*eye(6) A+B+C+lambda3*eye(6)]*v...
    -P*(lambda1+lambda2+lambda3)+norm([g1;g2;g3])^2*n0;


a=1;
%Normalize
%{
M1 = sqrt(norm(v11)^2+norm(v12)^2+norm(v13)^2);
v11 = v11/M1;
v12 = v12/M1;
v13 = v13/M1;

M2 = sqrt(norm(v21)^2+norm(v22)^2+norm(v23)^2);
v21 = v21/M2;
v22 = v22/M2;
v23 = v23/M2;

M3 = sqrt(norm(v31)^2+norm(v32)^2+norm(v33)^2);
v31 = v31/M3;
v32 = v32/M3;
v33 = v33/M3;
%}
end

    
    
