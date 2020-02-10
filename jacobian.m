% Jacobian

% parameters

E = [0, -1; 1, 0];

% inputs = [A_1 A_2 A_3 B_1 B_2 B_3 P]

rho_1 = norm(B_1 - A_1)
rho_2 = norm(B_2 - A_2)
rho_3 = norm(B_3 - A_3)

v_1 = (B_1-A_1)/norm(B_1-A_1);
v_2 = (B_2-A_2)/norm(B_2-A_2);
v_3 = (B_3-A_3)/norm(B_3-A_3);

% paralel Jacobian

A = [(E*v_1)',-(E*v_1)'*E*(P-B_1); ...
     (E*v_2)',-(E*v_2)'*E*(P-B_2); ...
     (E*v_3)',-(E*v_3)'*E*(P-B_3)];

% serial Jacobian 
B = diag([rho_1, rho_2, rho_3]);