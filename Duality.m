function [v11, v12, v13, v21, v22, v23, v31, v32, v33] = Duality(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, P, n0)
%Duality
stepsize = 10^(-5);

%Constants
k = [g1'*H11 g1'*H12 g1'*H13 g2'*H21 g2'*H22 g2'*H23 g3'*H31 g3'*H32 g3'*H33];
A = [H11'*g1*g1'*H11 H11'*g1*g1'*H12 H11'*g1*g1'*H13;H12'*g1*g1'*H11 H12'*g1*g1'*H12 H12'*g1*g1'*H13;H13'*g1*g1'*H11 H13'*g1*g1'*H12 H13'*g1*g1'*H13];
B = [H21'*g2*g2'*H21 H21'*g2*g2'*H22 H21'*g2*g2'*H23;H22'*g2*g2'*H21 H22'*g2*g2'*H22 H22'*g2*g2'*H23;H23'*g2*g2'*H21 H23'*g2*g2'*H22 H23'*g2*g2'*H23];
C = [H31'*g3*g3'*H31 H31'*g3*g3'*H32 H31'*g3*g3'*H33;H32'*g3*g3'*H31 H32'*g3*g3'*H32 H32'*g3*g3'*H33;H33'*g3*g3'*H31 H33'*g3*g3'*H32 H33'*g3*g3'*H33];

%Dual Problem
for o = -10^(4):10^4
lambda = o*10^(-3);
g(o+10^(4)+1) = 3 - k*inv([A+B+C+lambda*eye(6) 0*eye(6) 0*eye(6);0*eye(6) A+B+C+lambda*eye(6) 0*eye(6);0*eye(6) 0*eye(6) A+B+C+lambda*eye(6)])*k'-3*lambda*P+n0*(g1'*g1+g2'*g2+g3'*g3);
end
o = 1+0.5*10^4:(3/2)*10^4+1;
plot(o,g(o))
axis([1 2*10^4+1 -10^1 10^1])

%Eigenvalues of ABC Matrix
E = sort((eig([A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C])));
%Least Eigenvalues
E(10);


k*[A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C]*k';
det([A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C]);

%Gradients
G = -P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6)));

%Search for Optimal Point
stepsize = 10^(-3);
lambda = real(-E(10)/2)
for n = 1:10^(4)
%{
lambda = real(lambda + stepsize*(-P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6)))))
real(-P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6))))
    %}
lambda = real(lambda + stepsize*(-P+trace(inv(A+B+C+lambda*eye(6))*k(1:6)'*k(1:6)*inv(A+B+C+lambda*eye(6)))))
real(-P+trace(inv(A+B+C+lambda*eye(6))*k(1:6)'*k(1:6)*inv(A+B+C+lambda*eye(6))))
end

-E(10)
-E(10)/2
P11 = Z11*g11+Z12*g21+Z13*g31;
P12 = Z11*g12+Z12*g22+Z13*g32;
P13 = Z11*g13+Z12*g23+Z13*g33;

P21 = Z21*g11+Z22*g21+Z23*g31;
P22 = Z21*g12+Z22*g22+Z23*g32;
P23 = Z21*g13+Z22*g23+Z23*g33;

P31 = Z31*g11+Z32*g21+Z33*g31;
P32 = Z31*g12+Z32*g22+Z33*g32;
P33 = Z31*g13+Z32*g23+Z33*g33;

v11 = R1\P11;
v12 = R1\P12;
v13 = R1\P13;

v21 = R2\P21;
v22 = R2\P22;
v23 = R2\P23;

v31 = R3\P31;
v32 = R3\P32;
v33 = R3\P33;

%Normalize

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
end

    
    