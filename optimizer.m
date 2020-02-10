% global parameters
D = 10;
d = 1;
global conditionnement = 0.2;
global n = 10;
global taille = 10;

% initial conditions
x0 =[D; d]; % D d
lb =[0; 0];
ub =[15; 5];

options = optimset('Display','iter', 'maxiter',400,'MaxFunEvals',600,'TolFun',1e-8, 'TolX',1e-8)

x = fmincon('objectiveFun', x0, [], [], [], [], lb, ub , 'nonlcong', options)