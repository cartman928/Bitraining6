function [g1, g2] = MSE_f_2users_2antennas(H11, H12, H21, H22, v11, v12, v21, v22, n0)
%update filters by MSE criterion

%MSE
g1h = (H11*v11+H12*v21)'/...
     ( (H11*v11+H12*v21)*(H11*v11+H12*v21)'...
      +(H11*v12+H12*v22)*(H11*v12+H12*v22)'...
      +n0*eye(4)...
     );

g2h = (H21*v12+H22*v22)'/...
     ( (H21*v11+H22*v21)*(H21*v11+H22*v21)'...
      +(H21*v12+H22*v22)*(H21*v12+H22*v22)'...
      +n0*eye(4)...
     );
 

g1 = g1h';
g2 = g2h';

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

    
    