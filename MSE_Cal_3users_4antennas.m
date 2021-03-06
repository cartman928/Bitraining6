function [Lp] = MSE_Cal(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11, v12, v13,v21,v22,v23,v31,v32,v33,g1, g2, g3, P, n0)

k = [(g1')*H11 (g1')*H12 (g1')*H13 (g2')*H21 (g2')*H22 (g2')*H23 (g3')*H31 (g3')*H32 (g3')*H33];
A = [(H11')*g1*(g1')*H11 (H11')*g1*(g1')*H12 (H11')*g1*(g1')*H13;(H12')*g1*(g1')*H11 (H12')*g1*(g1')*H12 (H12')*g1*(g1')*H13;(H13')*g1*(g1')*H11 (H13')*g1*(g1')*H12 (H13')*g1*(g1')*H13];
B = [(H21')*g2*(g2')*H21 (H21')*g2*(g2')*H22 (H21')*g2*(g2')*H23;(H22')*g2*(g2')*H21 (H22')*g2*(g2')*H22 (H22')*g2*(g2')*H23;(H23')*g2*(g2')*H21 (H23')*g2*(g2')*H22 (H23')*g2*(g2')*H23];
C = [(H31')*g3*(g3')*H31 (H31')*g3*(g3')*H32 (H31')*g3*(g3')*H33;(H32')*g3*(g3')*H31 (H32')*g3*(g3')*H32 (H32')*g3*(g3')*H33;(H33')*g3*(g3')*H31 (H33')*g3*(g3')*H32 (H33')*g3*(g3')*H33];
ABC =  [A+B+C 0*eye(12) 0*eye(12);0*eye(12) A+B+C 0*eye(12);0*eye(12) 0*eye(12) A+B+C];
v = [v11;v21;v31;v12;v22;v32;v13;v23;v33];
Lp = 3-k*v-(k*v)'+v'*ABC*v+n0*(g1'*g1+g2'*g2+g3'*g3);

end