%% Stiffness estimation using Collage method
clear vars
clc
tic
 
%% DTU 10MW ( parameters can be found from its report)
 
% t=t;
N=1;%61.82;  % gear ratio
h=0.02;  % time step
t=0:h:600-0.02; % time vector
n=numel(t(loc:loc1));
syms J C K D
%% DTU 10MW Rigid parameters
% kt=2317025352; % main shaft stiffness alone
% kflex=4.7522e8;  % main shaft + blade edgewise flexibility

%% inputs
% rs- rotor speed time series
% gs - generator speed time series
% gp - generator power time series



E1=1:216; % For DLC 1.2, there are 216 simulations
j=zeroes(numel(E1),1);
k=zeroes(numel(E1),1);
c=zeroes(numel(E1),1);


for ii=1:numel(E1)

y=gs(loc:loc1,E1(ii));
E0=y-y(1);
v=rs(loc:loc1,E1(ii))-gs(loc:loc1,E1(ii));%gradient(HT1,0.02);
gt=gp(loc:loc1,E1(ii))./y;  % generator torque calculation
Eaa=(h*cumtrapz(v)); % integration of damping term of the generator equation 
Eb1=(h*cumtrapz(th(:,E1(ii))));  % integration of stiffness term of the generator equation where 'th-torsional displacement' obtained using Tikhonov regularisation
   
Ef1=((h*cumtrapz(gt-mean(gt))));% integration of force term of the generator equation 

Ek=(J*E0-(C)*Eaa-K*Eb1+Ef1);  % generator equation, J- generator inertia, C- damping of the drivetrain, K- stiffness of the drivetrain
% 
e1=diff(E,J);
e2=diff(E,C);
e3=diff(E,K);
% e4=diff(E,D);
% 
sol=solve( e1==0,e2==0,e3==0, J,C, K);
% 
j(ii)=vpa(sol.J,10);
c(ii)=vpa(sol.C,10);
kk=vpa(sol.K,10);  % estimated stiffness
% d=vpa(sol.D,10);
    if isempty(kk)==1 || isnan(kk)==1
        k(ii)=0;
    else
        k(ii)=kk;
    end
end
toc