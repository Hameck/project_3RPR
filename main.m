D = 10;
d = 3;

O = [0 0]';
A_1 = [0 0]';
A_2 = [D 0]';
A_3 = [D*1/2 D*sqrt(3)/2]';

rho_1 = 3*[cos(pi/3) sin(pi/3)]';
rho_2 = 3*[cos(2*pi/3) sin(2*pi/3)]';
rho_3 = 3*[cos(-pi/2) sin(-pi/2)]';

B_1 = A_1 + rho_1;
B_2 = A_2 + rho_2;
B_3 = A_3 + rho_3;

figure(1);
plot([A_1(1) B_1(1)],[A_1(2) B_1(2)], "-o", "linewidth", 3); hold on;
plot([A_2(1) B_2(1)],[A_2(2) B_2(2)], "-o", "linewidth", 3); hold on;
plot([A_3(1) B_3(1)],[A_3(2) B_3(2)], "-o", "linewidth", 3); hold on;
axis equal;