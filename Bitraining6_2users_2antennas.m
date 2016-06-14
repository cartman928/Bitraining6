 %% Duality Method
%% Initialize Parameters
clc
clear

alpha = 0;    %coefficient for block fading model
beta = 0.8^2;  %attenuation loss from non-direct antennas
n0 = 10^(-3);    %noise variance
P = 1; %power constraint

iternums = 1:10; % number of iterations
N_Realizations = 10;

C1 = zeros(N_Realizations, length(iternums));
C2 = zeros(N_Realizations, length(iternums));
C3 = zeros(N_Realizations, length(iternums));
MSE = zeros(N_Realizations, length(iternums));
MSEb = zeros(N_Realizations, length(iternums));
MSEc = zeros(N_Realizations, length(iternums));

%% Start Loop
for Realization = 1 : N_Realizations
    Realization
        
    %% Random Channels
    H11 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    H22 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    H33 = (randn(2,2)+1i*randn(2,2))/sqrt(2);
    
    H12 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H13 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H21 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta);
    H23 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H31 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta); 
    H32 = (randn(2,2)+1i*randn(2,2))/sqrt(2/beta);
    
    %Backward Channel
    
    Z11 = H11';
    Z22 = H22';
    Z33 = H33';
    
    Z12 = H21';
    Z13 = H31';
    Z21 = H12';
    Z23 = H32';
    Z31 = H13';
    Z32 = H23';
  
    
    %% one iteration per block
    g1 = rand(2, 1) + 1i*rand(2, 1); 
    g2 = rand(2, 1) + 1i*rand(2, 1);
    g3 = rand(2, 1) + 1i*rand(2, 1);
    g1 = g1/norm(g1);
    g2 = g2/norm(g2);
    g3 = g3/norm(g3);
    
    g1b = rand(2, 1) + 1i*rand(2, 1); 
    g2b = rand(2, 1) + 1i*rand(2, 1);
    g3b = rand(2, 1) + 1i*rand(2, 1);
    g1b = g1b/norm(g1b);
    g2b = g2b/norm(g2b);
    g3b = g3b/norm(g3b);
    
    g1c = rand(2, 1) + 1i*rand(2, 1); 
    g2c = rand(2, 1) + 1i*rand(2, 1);
    g3c = rand(2, 1) + 1i*rand(2, 1);
    g1c = g1c/norm(g1c);
    g2c = g2c/norm(g2c);
    g3c = g3c/norm(g3c);
    
    %2 users
    g3 = [0;0];
    g3b = [0;0];
    g3c = [0;0];
    
 
    %{
    v11 = zeros(2, 1); 
    v12 = zeros(2, 1);
    v13 = zeros(2, 1); 
    v21 = zeros(2, 1); 
    v22 = zeros(2, 1);
    v23 = zeros(2, 1); 
    v31 = zeros(2, 1); 
    v32 = zeros(2, 1); 
    v33 = zeros(2, 1); 
    %}
    
    for numiters = 1:length(iternums)

        %% bi-directional training
            
            %%Backward Training: sudo-LS Algorithm
            %abs(error)<2*10^(-4)
            %[v11a, v12a, v13a, v21a, v22a, v23a, v31a, v32a, v33a] = Duality(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, P, n0);
            [v11b, v12b, v13b, v21b, v22b, v23b, v31b, v32b, v33b] = MSE_b_2users_2antennas(Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33, g1b, g2b, g3b, n0);
            [v11c, v12c, v13c, v21c, v22c, v23c, v31c, v32c, v33c] = MSE_b_2users_2antennas(Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33, g1c, g2c, g3c, n0);
            v12c = [0;0];
            v13c = [0;0];
            v21c = [0;0];
            v23c = [0;0];
            v31c = [0;0];
            v32c = [0;0];
            
            v11c = v11c/norm(v11c);
            v22c = v22c/norm(v22c);
            
            %[v11, v12, v13]
            %[v21, v22, v23]
            %[v31, v32, v33]
            %Power = [norm(v11)^2+norm(v12)^2+norm(v13)^2 norm(v21)^2+norm(v22)^2+norm(v23)^2 norm(v31)^2+norm(v32)^2+norm(v33)^2]
            %[v11, v12, v13, v21, v22, v23, v31, v32, v33] = Duality(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, P, n0);
            [v11, v12, v13, v21, v22, v23, v31, v32, v33] = Primal_Dual_2users_2antennas(H11, H12, H13, H21, H22, H23, H31, H32, H33, g1, g2, g3, P, n0);
            
            %{
            [norm(v11)^2+norm(v12)^2+norm(v13)^2 norm(v21)^2+norm(v22)^2+norm(v23)^2 norm(v31)^2+norm(v32)^2+norm(v33)^2]
            [norm(v11b)^2+norm(v12b)^2+norm(v13b)^2 norm(v21b)^2+norm(v22b)^2+norm(v23b)^2 norm(v31b)^2+norm(v32b)^2+norm(v33b)^2]
            [norm(v11c)^2+norm(v12c)^2+norm(v13c)^2 norm(v21c)^2+norm(v22c)^2+norm(v23c)^2 norm(v31c)^2+norm(v32c)^2+norm(v33c)^2]
            %}
            
            
            %%%Verification
            %{
            k = [g1'*H11 g1'*H12 g1'*H13 g2'*H21 g2'*H22 g2'*H23 g3'*H31 g3'*H32 g3'*H33];
            A = [H11'*g1*g1'*H11 H11'*g1*g1'*H12 H11'*g1*g1'*H13;H12'*g1*g1'*H11 H12'*g1*g1'*H12 H12'*g1*g1'*H13;H13'*g1*g1'*H11 H13'*g1*g1'*H12 H13'*g1*g1'*H13];
            B = [H21'*g2*g2'*H21 H21'*g2*g2'*H22 H21'*g2*g2'*H23;H22'*g2*g2'*H21 H22'*g2*g2'*H22 H22'*g2*g2'*H23;H23'*g2*g2'*H21 H23'*g2*g2'*H22 H23'*g2*g2'*H23];
            C = [H31'*g3*g3'*H31 H31'*g3*g3'*H32 H31'*g3*g3'*H33;H32'*g3*g3'*H31 H32'*g3*g3'*H32 H32'*g3*g3'*H33;H33'*g3*g3'*H31 H33'*g3*g3'*H32 H33'*g3*g3'*H33];
            
            va = [v11a;v21a;v31a;v12a;v22a;v32a;v13a;v23a;v33a];
            La = 3-k*va-(k*va)'+va'*[A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C]*va...
            +norm([g1;g2;g3])^2*n0
        
            vb = [v11b;v21b;v31b;v12b;v22b;v32b;v13b;v23b;v33b];
            Lb = 3-k*vb-(k*vb)'+vb'*[A+B+C 0*eye(6) 0*eye(6);0*eye(6) A+B+C 0*eye(6);0*eye(6) 0*eye(6) A+B+C]*vb...
            +norm([g1;g2;g3])^2*n0
            %maximize kv
            %minimize v'ABCv
            %}
            
            
            %%Forward Training: LS Algorithm
            [g1, g2, g3] = MSE_f_2users_2antennas(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11, v12, v13, v21, v22, v23, v31, v32, v33, n0);
            [g1b, g2b, g3b] = MSE_f_2users_2antennas(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11b, v12b, v13b, v21b, v22b, v23b, v31b, v32b, v33b, n0);
            [g1c, g2c, g3c] = MSE_f_2users_2antennas(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11c, v12c, v13c, v21c, v22c, v23c, v31c, v32c, v33c, n0);
            %[g1, g12, g13, g21, g2, g23, g31, g32, g3] = Duality(Z11, Z12, Z13, Z21, Z22, Z23, Z31, Z32, Z33, v11, v22, v33, P, n0);
            %[norm(g1)^2 norm(g2)^2 norm(g3)^2 norm(g1b)^2 norm(g2b)^2 norm(g3b)^2 norm(g1c)^2 norm(g2c)^2 norm(g3c)^2]
        
            %Power = [norm(g11)^2+norm(g12)^2+norm(g13)^2 norm(g21)^2+norm(g22)^2+norm(g23)^2 norm(g31)^2+norm(g32)^2+norm(g33)^2]

            %{
            g12 = zeros(2, 1);
            g13 = zeros(2, 1); 
            g21 = zeros(2, 1); 
            g23 = zeros(2, 1); 
            g31 = zeros(2, 1); 
            g32 = zeros(2, 1);
            %}
            %{
            MSE(Realization, numiters)=MSE_Cal(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11, v12, v13,v21,v22,v23,v31,v32,v33,g1, g2, g3, P, n0);
            MSEb(Realization, numiters)=MSE_Cal(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11b, v12b, v13b,v21b,v22b,v23b,v31b,v32b,v33b,g1b, g2b, g3b, P, n0);
            MSEc(Realization, numiters)=MSE_Cal(H11, H12, H13, H21, H22, H23, H31, H32, H33, v11c, v12c, v13c,v21c,v22c,v23c,v31c,v32c,v33c,g1c, g2c, g3c, P, n0);
            %}
            SINR1 = norm(g1'*(H11*v11+H12*v21+H13*v31))^2/((norm(g1'*(H11*v12+H12*v22+H13*v32)))^2+(norm(g1'*(H11*v13+H12*v23+H13*v33)))^2+n0*g1'*g1);
            SINR2 = norm(g2'*(H21*v12+H22*v22+H23*v32))^2/((norm(g2'*(H21*v11+H22*v21+H23*v31)))^2+(norm(g2'*(H21*v13+H22*v23+H23*v33)))^2+n0*g2'*g2);
            SINR3 = norm(g3'*(H31*v13+H32*v23+H33*v33))^2/((norm(g3'*(H31*v11+H32*v21+H33*v31)))^2+(norm(g3'*(H31*v12+H32*v22+H33*v32)))^2+n0*g3'*g3);
            C1(Realization, numiters) = abs(log2(1+SINR1));
            C2(Realization, numiters) = abs(log2(1+SINR2));
            C3(Realization, numiters) = abs(log2(1+SINR3));
            
            
            SINR1b = norm(g1b'*(H11*v11b+H12*v21b+H13*v31b))^2/((norm(g1b'*(H11*v12b+H12*v22b+H13*v32b)))^2+(norm(g1b'*(H11*v13b+H12*v23b+H13*v33b)))^2+n0*g1b'*g1b);
            SINR2b = norm(g2b'*(H21*v12b+H22*v22b+H23*v32b))^2/((norm(g2b'*(H21*v11b+H22*v21b+H23*v31b)))^2+(norm(g2b'*(H21*v13b+H22*v23b+H23*v33b)))^2+n0*g2b'*g2b);
            SINR3b = norm(g3b'*(H31*v13b+H32*v23b+H33*v33b))^2/((norm(g3b'*(H31*v11b+H32*v21b+H33*v31b)))^2+(norm(g3b'*(H31*v12b+H32*v22b+H33*v32b)))^2+n0*g3b'*g3b);
            C1b(Realization, numiters) = abs(log2(1+SINR1b));
            C2b(Realization, numiters) = abs(log2(1+SINR2b));
            C3b(Realization, numiters) = abs(log2(1+SINR3b));
            
            SINR1c = norm(g1c'*(H11*v11c+H12*v21c+H13*v31c))^2/((norm(g1c'*(H11*v12c+H12*v22c+H13*v32c)))^2+(norm(g1c'*(H11*v13c+H12*v23c+H13*v33c)))^2+n0*g1c'*g1c);
            SINR2c = norm(g2c'*(H21*v12c+H22*v22c+H23*v32c))^2/((norm(g2c'*(H21*v11c+H22*v21c+H23*v31c)))^2+(norm(g2c'*(H21*v13c+H22*v23c+H23*v33c)))^2+n0*g2c'*g2c);
            SINR3c = norm(g3c'*(H31*v13c+H32*v23c+H33*v33c))^2/((norm(g3c'*(H31*v11c+H32*v21c+H33*v31c)))^2+(norm(g3c'*(H31*v12c+H32*v22c+H33*v32c)))^2+n0*g3c'*g3c);
            C1c(Realization, numiters) = abs(log2(1+SINR1c));
            C2c(Realization, numiters) = abs(log2(1+SINR2c));
            C3c(Realization, numiters) = abs(log2(1+SINR3c));
            
            %{
            SINR1 = norm(g1'*(H11*v11+H12*v21+H13*v31))^2/((norm(g1'*(H11*v12+H12*v22+H13*v32)))^2+(norm(g1'*(H11*v13+H12*v23+H13*v33)))^2);
            SINR2 = norm(g2'*(H21*v12+H22*v22+H23*v32))^2/((norm(g2'*(H21*v11+H22*v21+H23*v31)))^2+(norm(g2'*(H21*v13+H22*v23+H23*v33)))^2);
            SINR3 = norm(g3'*(H31*v13+H32*v23+H33*v33))^2/((norm(g3'*(H31*v12+H32*v22+H33*v32)))^2+(norm(g3'*(H31*v12+H32*v22+H33*v32)))^2);
            C1(Realization, numiters) = abs(log2(1+SINR1));
            C2(Realization, numiters) = abs(log2(1+SINR2));
            C3(Realization, numiters) = abs(log2(1+SINR3));
            
            
            SINR1b = norm(g1b'*(H11*v11b+H12*v21b+H13*v31b))^2/((norm(g1b'*(H11*v12b+H12*v22b+H13*v32b)))^2+(norm(g1b'*(H11*v13b+H12*v23b+H13*v33b)))^2);
            SINR2b = norm(g2b'*(H21*v12b+H22*v22b+H23*v32b))^2/((norm(g2b'*(H21*v11b+H22*v21b+H23*v31b)))^2+(norm(g2b'*(H21*v13b+H22*v23b+H23*v33b)))^2);
            SINR3b = norm(g3b'*(H31*v13b+H32*v23b+H33*v33b))^2/((norm(g3b'*(H31*v12b+H32*v22b+H33*v32b)))^2+(norm(g3b'*(H31*v12b+H32*v22b+H33*v32b)))^2);
            C1b(Realization, numiters) = abs(log2(1+SINR1b));
            C2b(Realization, numiters) = abs(log2(1+SINR2b));
            C3b(Realization, numiters) = abs(log2(1+SINR3b));
            
            SINR1c = norm(g1c'*(H11*v11c+H12*v21c+H13*v31c))^2/((norm(g1c'*(H11*v12c+H12*v22c+H13*v32c)))^2+(norm(g1c'*(H11*v13c+H12*v23c+H13*v33c)))^2);
            SINR2c = norm(g2c'*(H21*v12c+H22*v22c+H23*v32c))^2/((norm(g2c'*(H21*v11c+H22*v21c+H23*v31c)))^2+(norm(g2c'*(H21*v13c+H22*v23c+H23*v33c)))^2);
            SINR3c = norm(g3c'*(H31*v13c+H32*v23c+H33*v33c))^2/((norm(g3c'*(H31*v12c+H32*v22c+H33*v32c)))^2+(norm(g3c'*(H31*v12c+H32*v22c+H33*v32c)))^2);
            C1c(Realization, numiters) = abs(log2(1+SINR1c));
            C2c(Realization, numiters) = abs(log2(1+SINR2c));
            C3c(Realization, numiters) = abs(log2(1+SINR3c));
            %}
            
            if norm(g1)^2 > 1
            g1 = g1/norm(g1);
            end
            if norm(g2)^2 > 1
            g2 = g2/norm(g2);
            end
            %{
            if norm(g3)^2 > 1
            g3 = g3/norm(g3);
            end
            %}
            
            %if norm(g1b)^2 > 1
            g1b = g1b/norm(g1b);
            %end
            %if norm(g2b)^2 > 1
            g2b = g2b/norm(g2b);
            %end
            %if norm(g3b)^2 > 1
            %g3b = g3b/norm(g3b);
            %end
            
            %if norm(g1c)^2 > 1
            g1c = g1c/norm(g1c);
            %end
            %if norm(g2c)^2 > 1
            g2c = g2c/norm(g2c);
            %end
            %if norm(g3c)^2 > 1
            %g3c = g3c/norm(g3c);
            %end
            
            end           
    
end




%% Plot C(bits/channel)
figure
hold on

p1=plot(iternums, mean(C1)+mean(C2),'o');
p2=plot(iternums, mean(C1b)+mean(C2b),'*');
p3=plot(iternums, mean(C1c)+mean(C2c),'s');

axis([1 numiters 0 40])

xlabel('Number of iterations')
ylabel('C(bits/channel)')
title('Capacity vs. Iterations: 3 users;P=1;n0=10^{-3};100 realizations')
legend([p1,p2,p3],'Cooperative Transmitters(Duality Method)','Cooperative Transmitters(R^{-1}p)','Simple Transmitters')

%{
p1=plot(iternums, mean(MSE),'o');
p2=plot(iternums, mean(MSEb),'*');
p3=plot(iternums, mean(MSEc),'s');

xlabel('Number of iterations')
ylabel('MSE')
title('MSE vs. Iterations: 3 users;P=1;n0=10^{-2};100 realizations')
legend([p1,p2,p3],'Cooperative Transmitters(Duality Method)','Cooperative Transmitters(R^{-1}p)','Simple Transmitters')
%}

%}
