function [v11, v12, v13, v21, v22, v23, v31, v32, v33] = Primal_Dual(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, P, n0)
%Duality

%Constants
lambda = [0;0;0];
v = zeros(18,1);
Lam = diag([lambda(1) lambda(1) lambda(2) lambda(2) lambda(3) lambda(3)]);
k = [g1'*H11 g1'*H12 g1'*H13 g2'*H21 g2'*H22 g2'*H23 g3'*H31 g3'*H32 g3'*H33];
A = [H11'*g1*g1'*H11 H11'*g1*g1'*H12 H11'*g1*g1'*H13;H12'*g1*g1'*H11 H12'*g1*g1'*H12 H12'*g1*g1'*H13;H13'*g1*g1'*H11 H13'*g1*g1'*H12 H13'*g1*g1'*H13];
B = [H21'*g2*g2'*H21 H21'*g2*g2'*H22 H21'*g2*g2'*H23;H22'*g2*g2'*H21 H22'*g2*g2'*H22 H22'*g2*g2'*H23;H23'*g2*g2'*H21 H23'*g2*g2'*H22 H23'*g2*g2'*H23];
C = [H31'*g3*g3'*H31 H31'*g3*g3'*H32 H31'*g3*g3'*H33;H32'*g3*g3'*H31 H32'*g3*g3'*H32 H32'*g3*g3'*H33;H33'*g3*g3'*H31 H33'*g3*g3'*H32 H33'*g3*g3'*H33];
ABCLam =  [A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam];

%Priaml Problem
%L = 3-k*v-(k*v)'+v'*ABCLam*v+n0*(g1'*g1+g2'*g2+g3'*g3)-P*(lambda(1)+lambda(2)+lambda(3));
%L_g = -2*k+2*v'*ABCLam;

%Dual Problem
%D_g = [norm([v11;v12;v13])^2;norm([v21;v22;v23])^2;norm([v31;v32;v33])^2];

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

v = zeros(18,1);
stepsize = 10^(-5);
%Primal Upadate

for n = 1:10^8
L_g = -2*k+2*v'*ABCLam;
%norm([L_g(1:2) L_g(3:4) L_g(5:6) L_g(7:8) L_g(9:10) L_g(11:12) L_g(13:14) L_g(15:16) L_g(17:18)])
[v(1) norm(L_g(1))]
%v = v - stepsize*transpose(L_g);
v = v - stepsize*[(L_g(1))';zeros(17,1)];
end
%Dual Update
for n = 1:10^5
D_g = [norm([v(1:2);v(7:8);v(13:14)])^2;norm([v(3:4);v(9:10);v(15:16)])^2;norm([v(5:6);v(11:12);v(17:18)])^2];
lambda = lambda + stepsize*D_g;
end

%{
for o = -10^(5):10^5
lambda1 = o*10^(-5);
lambda2 = 5;
lambda3 = 5;
Lam = diag([lambda1 lambda1 lambda2 lambda2 lambda3 lambda3]);

g(o+10^5+1) = 3 - k*inv([A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam])*k'-P*(lambda1+lambda2+lambda3)+n0*(g1'*g1+g2'*g2+g3'*g3);
end
o = 1:2*10^5+1;
plot(o,g(o))
%axis([1 2*10^5+1 -5 5])
[M,I] = max(real(g(10^5+1-200:10^5+1+200)))

real(g(10^5+1-200:10^5+1+200));
%}


%Eigenvalues of ABC Matrix
E = sort((eig([A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C])));
%Least Eigenvalues
E(10);

%Gradients
%G = -P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6)));

%Search for Optimal Point
%stepsize = 1/real(E(18));
stepsize = 10^(-4);
%{
lambda1 = real(E(10)/2);
lambda2 = real(E(10)/2);
lambda3 = real(E(10)/2);
%}



lambda1 = 5;
lambda2 = 5;
lambda3 = 5;

%{
lambda1 = 0;
lambda2 = 0;
lambda3 = 0;
%}


for n = 1:10^(5)
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




a = real(3 - k*inv([A+B+C+Lam 0*eye(6) 0*eye(6);0*eye(6) A+B+C+Lam 0*eye(6);0*eye(6) 0*eye(6) A+B+C+Lam])*k'-P*(lambda1+lambda2+lambda3)+n0*(g1'*g1+g2'*g2+g3'*g3));
[a norm([gradient1 gradient2 gradient3]) lambda1 lambda2 lambda3]

if(norm([gradient1 gradient2 gradient3])< 10^(-4)  )
            disp('converges!');
            break;%break the for loop if it's true the condition
end




end

%{
lambda1 = n0;
lambda2 = n0;
lambda3 = n0;
%}
Lam = diag([lambda1 lambda1 lambda2 lambda2 lambda3 lambda3])

V=(k*[inv([A+B+C+Lam]) 0*eye(6) 0*eye(6);0*eye(6)   inv([A+B+C+Lam])   0*eye(6);0*eye(6) 0*eye(6) inv([A+B+C+Lam])])';

v11 = V(1:2);
v21 = V(3:4);
v31 = V(5:6);

v12 = V(7:8);
v22 = V(9:10);
v32 = V(11:12);

v13 = V(13:14);
v23 = V(15:16);
v33 = V(17:18);

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

    
    
