function [G Geq]=nonlcong(x)

    global conditionnement pointsList;
       
    % conditionement
    condList = [];
    
    D = x(1)
    d = x(2)
    
    O = [0 0]';
    A_1 = [0 0]';
    A_2 = [D 0]';
    A_3 = [D*1/2 D*sqrt(3)/2]';
            
    for i=1:length(pointsList)
      for phi=-pi/6:pi/36:pi/6   
        % end-effector position
        x = pointsList(1,i);
        y = pointsList(2,i);
        
        % inverse kinematics

        B_1 = [x+d*cos(phi-5*pi/6) y+d*sin(phi-5*pi/6)]';
        B_2 = [x+d*cos(phi-pi/6) y+d*sin(phi-pi/6)]';
        B_3 = [x+d*cos(phi+pi/2) y+d*sin(phi+pi/2)]';

        P = [x y]';

        % output = [rho_1 rho_2 rho_3]

        rho_1 = norm(B_1 - A_1);
        rho_2 = norm(B_2 - A_2);
        rho_3 = norm(B_3 - A_3);

        % jacobian

        E = [0, -1; 1, 0];

        v_1 = (B_1-A_1)/norm(B_1-A_1);
        v_2 = (B_2-A_2)/norm(B_2-A_2);
        v_3 = (B_3-A_3)/norm(B_3-A_3);

        % paralel Jacobian

        A = [(v_1)',-(v_1)'*E*(P-B_1); ...
             (v_2)',-(v_2)'*E*(P-B_2); ...
             (v_3)',-(v_3)'*E*(P-B_3)]
        
        M = isnan(A)
         
        if(det(A) == 0 || any(M(:)))
          disp('Paralel singularity');
          figure(1);
          plot([A_1(1) B_1(1)],[A_1(2) B_1(2)], '-o', 'linewidth', 3); hold on;
          plot([A_2(1) B_2(1)],[A_2(2) B_2(2)], '-o', 'linewidth', 3); hold on;
          plot([A_3(1) B_3(1)],[A_3(2) B_3(2)], '-o', 'linewidth', 3); hold on;
          plot([B_1(1) B_2(1) B_3(1) B_1(1)],[B_1(2) B_2(2) B_3(2) B_1(2)], '-o', 'linewidth', 3);
          plot([P(1)],[P(2)], '-+', 'linewidth', 2); hold on;
          axis equal;
          grid on;
          pause()
          condList = [condList 1e-16];
          break; break;
        else
          % serial Jacobian 
          B = eye(3);

          % condition
          J = pinv(A)*B;
          cond(J)
          condList = [condList (conditionnement*cond(J)-1)];
        end
      end
    end
    
    G = [max(condList)];

    Geq = [];
end