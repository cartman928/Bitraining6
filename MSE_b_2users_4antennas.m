function [v11, v12, v21, v22] = MSE_b(Z11, Z12, Z21, Z22, g1, g2, n0)
%update filters by MSE criterion

%MSE
R1 = (Z11*g1)*(Z11*g1)'+(Z12*g2)*(Z12*g2)'+n0*eye(4);
R2 = (Z21*g1)*(Z21*g1)'+(Z22*g2)*(Z22*g2)'+n0*eye(4);

P11 = Z11*g1;
P12 = Z12*g2;

P21 = Z21*g1;
P22 = Z22*g2;


v11 = R1\P11;
v12 = R1\P12;

v21 = R2\P21;
v22 = R2\P22;



%Normalize

M1 = sqrt(norm(v11)^2+norm(v12)^2);
v11 = v11/M1;
v12 = v12/M1;

M2 = sqrt(norm(v21)^2+norm(v22)^2);
v21 = v21/M2;
v22 = v22/M2;


end

    
    