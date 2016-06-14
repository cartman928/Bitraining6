function [v11, v12, v13, v21, v22, v23, v31, v32, v33] = MSE_b_3users_4antennas(Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33, g1, g2, g3, n0)
%update filters by MSE criterion

%MSE
R1 = (Z11*g1)*(Z11*g1)'+(Z12*g2)*(Z12*g2)'+(Z13*g3)*(Z13*g3)'+n0*eye(4);
R2 = (Z21*g1)*(Z21*g1)'+(Z22*g2)*(Z22*g2)'+(Z23*g3)*(Z23*g3)'+n0*eye(4);
R3 = (Z31*g1)*(Z31*g1)'+(Z32*g2)*(Z32*g2)'+(Z33*g3)*(Z33*g3)'+n0*eye(4);

P11 = Z11*g1;
P12 = Z12*g2;
P13 = Z13*g3;

P21 = Z21*g1;
P22 = Z22*g2;
P23 = Z23*g3;

P31 = Z31*g1;
P32 = Z32*g2;
P33 = Z33*g3;

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

    
    