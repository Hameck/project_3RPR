function [G Geq]=nonlcong(x)

    global n taille conditionnement;
    
    % Points d'interet
    pointsList = [];
    for i=-(taille/2):(taille/n):(taille/2)
        for j=-(taille/2):(taille/n):+(taille/2)
            pointsList = [pointsList, [i+x(1)*1/2;j+x(1)*sqrt(3)/2]];
        end
    end
    
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

        B_1 = [x-d*sin(phi+pi/2) y+d*cos(phi+pi/2)]';
        B_2 = [x-d*sin(phi-5*pi/6) y+d*cos(phi-5*pi/6)]';
        B_3 = [x-d*sin(phi-pi/6) y+d*cos(phi-pi/6)]';

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
             (v_3)',-(v_3)'*E*(P-B_3)];
             
        if(det(A) == 0)
          disp('Paralel singularity');
          condList = [condList 1e-16];
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