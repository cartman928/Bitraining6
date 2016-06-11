function [g1, g2, g3] = MSE_f_2users(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11, v12, v13, v21, v22, v23, v31, v32, v33, n0)
%update filters by MSE criterion

%MSE
g1h = (H11*v11+H12*v21+H13*v31)'/...
     ( (H11*v11+H12*v21+H13*v31)*(H11*v11+H12*v21+H13*v31)'...
      +(H11*v12+H12*v22+H13*v32)*(H11*v12+H12*v22+H13*v32)'...
      +(H11*v13+H12*v23+H13*v33)*(H11*v13+H12*v23+H13*v33)'...
      +n0*eye(2)...
     );

g2h = (H21*v12+H22*v22+H23*v32)'/...
     ( (H21*v11+H22*v21+H23*v31)*(H21*v11+H22*v21+H23*v31)'...
      +(H21*v12+H22*v22+H23*v32)*(H21*v12+H22*v22+H23*v32)'...
      +(H21*v13+H22*v23+H23*v33)*(H21*v13+H22*v23+H23*v33)'...
      +n0*eye(2)...
     );
 

g1 = g1h';
g2 = g2h';
g3 = [0;0];

%Normalize
%{
if norm(g1)^2 > 1
g1 = g1/norm(g1);
end
if norm(g2)^2 > 1
g2 = g2/norm(g2);
end
if norm(g3)^2 > 1
g3 = g3/norm(g3);
end
%}
end

    
    