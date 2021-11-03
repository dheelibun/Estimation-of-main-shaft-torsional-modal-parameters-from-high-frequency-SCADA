% Tikhonov regularisation to obtain displacement from velocity
clearvars
clc
%% 
t=0:0.02:600-0.02;  % time vector
%% inputs 
% rs - rotor speed
% gs - generator speed
% lambda - regularisation parameter

%% This section of the code needs to be executed only once for generating the Tikhonoc matrices
% load tik_matrices_mes
% load tik_matrix_inverse_lambda_1e-4.mat
n=numel(t)-2;
La=eye(n);
La(1,1)=sqrt(1/2); 
La(end,end)=sqrt(1/2); 
%precompute:
nn=n+2;
subdiagidx = 2:nn+1:nn^2;
diagidx = 1:nn+1:nn^2;
superdiagidx = nn+1:nn+1:nn^2;
subsuperidx = [subdiagidx, superdiagidx];
%gemeration
td = zeros(nn);
td(subdiagidx) = -1;
td(superdiagidx) = 1;
% 
% %%
Lc=td(2:end-1,:);
L=La*Lc;

Ltl=0.25*(L'*L);
Lta=L'*La;
I=eye(size(L,2));
save('tik_matrices.mat','-v7.3','La','Ltl','Lta','L')

%% 
load tik_matrices.mat


tic
% E1 is number of 10-simulations


for i=1:numel(E1)%:size(rs1,2)

    a=rs(loc:l1,ind(jj))-gs(loc:l1,ind(jj))/N; % rs- rotor speed, gs- generator speed, N - gear ratio

    sd=(Ltl+lambda(i)^2*I);  % lambda is the regularisation parameter

    x_reg=sd\(Lta*a*h*0.5);

    th(:,i)=x_reg(2:end-1,:);  % regularised torsional displacement

end
% save('tik_theta_lambda_1e-4','-v7.3','th')
toc
% The same code can be modified for finding the optimal lambda for each
% mean wind speed