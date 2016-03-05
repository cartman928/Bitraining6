function [v11, v12, v13, v21, v22, v23, v31, v32, v33] = Duality(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, P, n0)
%Duality
stepsize = 10^(-5);

%Constants
k = [g1'*H11 g1'*H12 g1'*H13 g2'*H21 g2'*H22 g2'*H23 g3'*H31 g3'*H32 g3'*H33];
A = [H11'*g1*g1'*H11 H11'*g1*g1'*H12 H11'*g1*g1'*H13;H12'*g1*g1'*H11 H12'*g1*g1'*H12 H12'*g1*g1'*H13;H13'*g1*g1'*H11 H13'*g1*g1'*H12 H13'*g1*g1'*H13];
B = [H21'*g2*g2'*H21 H21'*g2*g2'*H22 H21'*g2*g2'*H23;H22'*g2*g2'*H21 H22'*g2*g2'*H22 H22'*g2*g2'*H23;H23'*g2*g2'*H21 H23'*g2*g2'*H22 H23'*g2*g2'*H23];
C = [H31'*g3*g3'*H31 H31'*g3*g3'*H32 H31'*g3*g3'*H33;H32'*g3*g3'*H31 H32'*g3*g3'*H32 H32'*g3*g3'*H33;H33'*g3*g3'*H31 H33'*g3*g3'*H32 H33'*g3*g3'*H33];

Aa = [H11'*g1*g1'*H11 H11'*g1*g1'*H12 ;H12'*g1*g1'*H11 H12'*g1*g1'*H12];
%Dual Problem
%{
for o = -10^(4):10^4
lambda1 = o*10^(-3);
lambda2 = 2;
lambda3 = 2;

g(o+10^(4)+1) = 3 - k*inv([A+B+C+lambda1*eye(6) 0*eye(6) 0*eye(6);0*eye(6) A+B+C+lambda2*eye(6) 0*eye(6);0*eye(6) 0*eye(6) A+B+C+lambda3*eye(6)])*k'-P*(lambda1+lambda2+lambda3)+n0*(g1'*g1+g2'*g2+g3'*g3);
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
stepsize = 10^(-2);
lambda1 = real(-E(10)/2);
lambda2 = real(-E(10)/2);
lambda3 = real(-E(10)/2);

for n = 1:10^(3)
%{
lambda = real(lambda + stepsize*(-P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6)))))
real(-P+norm(k)^2*trace(inv(A+B+C+lambda*eye(6))*inv(A+B+C+lambda*eye(6))))
    %}
gradient1 = real(-P+trace(inv(A+B+C+lambda1*eye(6))*k(1:6)'*k(1:6)*inv(A+B+C+lambda1*eye(6))));
lambda1 = real(lambda1 + stepsize*gradient1);

gradient2 = real(-P+trace(inv(A+B+C+lambda2*eye(6))*k(7:12)'*k(7:12)*inv(A+B+C+lambda2*eye(6))));
lambda2 = real(lambda2 + stepsize*gradient2);

gradient3 = real(-P+trace(inv(A+B+C+lambda3*eye(6))*k(13:18)'*k(13:18)*inv(A+B+C+lambda3*eye(6))));
lambda3 = real(lambda3 + stepsize*gradient3);



if((norm(gradient1)+norm(gradient2)+norm(gradient3))< 10^(-6)  )
            %disp('converges!');
            break;%break the for loop if it's true the condition
end




end

%{
lambda1 = 10^(-8);
lambda2 = 10^(-8);
lambda3 = 10^(-8);
%}

V1=(k(1:6)*inv(A+B+C+lambda1*eye(6)))';
V2=(k(7:12)*inv(A+B+C+lambda2*eye(6)))';
V3=(k(13:18)*inv(A+B+C+lambda3*eye(6)))';
v11 = V1(1:2);
v12 = V1(3:4);
v13 = V1(5:6);

v21 = V2(1:2);
v22 = V2(3:4);
v23 = V2(5:6);

v31 = V3(1:2);
v32 = V3(3:4);
v33 = V3(5:6);

(norm(v11)^2+norm(v12)^2+norm(v13)^2);




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

    
    