% global parameters
D = 10;
d = 1;
global conditionnement n taille;
conditionnement = 0.2;
n = 10;
taille = 5;

% initial conditions
x0 =[D; d]; % D d
lb =[1; 0.1];
ub =[15; 5];

options = optimset('Display','iter', 'maxiter',400,'MaxFunEvals',600,'TolFun',1e-8, 'TolX',1e-8)

x = fmincon('objectiveFun', x0, [], [], [], [], lb, ub , 'nonlcong', options)