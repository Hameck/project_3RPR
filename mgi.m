% MGI

% global parameters

D = 10;
d = 1;

O = [0 0]';
A_1 = [0 0]';
A_2 = [D 0]';
A_3 = [D*1/2 D*sqrt(3)/2]';

% input = [x, y, phi]

x = 7;
y = 4;
phi = 3*pi/4;

% inverse kinematics

B_1 = [x-d*sin(phi+pi/2) y+d*cos(phi+pi/2)]';
B_2 = [x-d*sin(phi-5*pi/6) y+d*cos(phi-5*pi/6)]';
B_3 = [x-d*sin(phi-pi/6) y+d*cos(phi-pi/6)]';

P = [x y]';

% output = [rho_1 rho_2 rho_3]

rho_1 = norm(B_1 - A_1)
rho_2 = norm(B_2 - A_2)
rho_3 = norm(B_3 - A_3)

% draw

figure(1);
% plot([A_1(1) A_2(1) A_3(1) A_1(1)],[A_1(2) A_2(2) A_3(2) A_1(2)], "-o", "linewidth", 3); hold on;
plot([A_1(1) B_1(1)],[A_1(2) B_1(2)], "-o", "linewidth", 3); hold on;
plot([A_2(1) B_2(1)],[A_2(2) B_2(2)], "-o", "linewidth", 3); hold on;
plot([A_3(1) B_3(1)],[A_3(2) B_3(2)], "-o", "linewidth", 3); hold on;
plot([B_1(1) B_2(1) B_3(1) B_1(1)],[B_1(2) B_2(2) B_3(2) B_1(2)], "-o", "linewidth", 3);
plot([P(1)],[P(2)], "-o", "linewidth", 3); hold on;
axis equal;
grid on;