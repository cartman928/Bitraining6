function [v11, v12, v13, v14, v21, v22, v23, v24, v31, v32, v33, v34, v41, v42, v43, v44] = MSE_b_4users(Z11, Z12, Z13, Z14, Z21, Z22, Z23, Z24, Z31, Z32, Z33, Z34, Z41, Z42, Z43, Z44, g1, g2, g3, g4, n0)
%update filters by MSE criterion

%MSE
R1 = (Z11*g1)*(Z11*g1)'+(Z12*g2)*(Z12*g2)'+(Z13*g3)*(Z13*g3)'+(Z14*g4)*(Z14*g4)'+n0*eye(2);
R2 = (Z21*g1)*(Z21*g1)'+(Z22*g2)*(Z22*g2)'+(Z23*g3)*(Z23*g3)'+(Z24*g4)*(Z24*g4)'+n0*eye(2);
R3 = (Z31*g1)*(Z31*g1)'+(Z32*g2)*(Z32*g2)'+(Z33*g3)*(Z33*g3)'+(Z34*g4)*(Z34*g4)'+n0*eye(2);
R4 = (Z41*g1)*(Z41*g1)'+(Z42*g2)*(Z42*g2)'+(Z43*g3)*(Z43*g3)'+(Z44*g4)*(Z44*g4)'+n0*eye(2);

P11 = Z11*g1;
P12 = Z12*g2;
P13 = Z13*g3;
P14 = Z14*g4;

P21 = Z21*g1;
P22 = Z22*g2;
P23 = Z23*g3;
P24 = Z24*g4;

P31 = Z31*g1;
P32 = Z32*g2;
P33 = Z33*g3;
P34 = Z34*g4;

P41 = Z41*g1;
P42 = Z42*g2;
P43 = Z43*g3;
P44 = Z44*g4;

v11 = R1\P11;
v12 = R1\P12;
v13 = R1\P13;
v14 = R1\P14;

v21 = R2\P21;
v22 = R2\P22;
v23 = R2\P23;
v24 = R2\P24;

v31 = R3\P31;
v32 = R3\P32;
v33 = R3\P33;
v34 = R3\P34;

v41 = R4\P41;
v42 = R4\P42;
v43 = R4\P43;
v44 = R4\P44;





%Normalize

M1 = sqrt(norm(v11)^2+norm(v12)^2+norm(v13)^2+norm(v14)^2);
v11 = v11/M1;
v12 = v12/M1;
v13 = v13/M1;
v14 = v14/M1;

M2 = sqrt(norm(v21)^2+norm(v22)^2+norm(v23)^2+norm(v24)^2);
v21 = v21/M2;
v22 = v22/M2;
v23 = v23/M2;
v24 = v24/M2;

M3 = sqrt(norm(v31)^2+norm(v32)^2+norm(v33)^2+norm(v34)^2);
v31 = v31/M3;
v32 = v32/M3;
v33 = v33/M3;
v34 = v34/M3;

M4 = sqrt(norm(v41)^2+norm(v42)^2+norm(v43)^2+norm(v44)^2);
v41 = v41/M4;
v42 = v42/M4;
v43 = v43/M4;
v44 = v44/M4;

end

    
    