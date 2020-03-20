function [G Geq]=nonlcong(x)

    global conditionnement pointsList minPhi maxPhi maxJointLimit minJointLimit maxJointLimitEnabled minJointLimitEnabled condEnabled;
       
    % conditionement
    condList = [];
    % minArmLengths
    min_rho = [];
    % maxArmLengths
    max_rho = [];
    
    D = x(1);
    d = x(2);
    
    O = [0 0]';
    A_1 = [0 0]';
    A_2 = [D 0]';
    A_3 = [D*1/2 D*sqrt(3)/2]';
            
    for i=1:length(pointsList)
      for phi=deg2rad(minPhi):pi/36:deg2rad(maxPhi)

        % end-effector position
        x = pointsList(1,i)+D/2;
        y = pointsList(2,i)+ D*sqrt(3)/4;
        
        % inverse kinematics

        B_1 = [x+d*cos(phi-5*pi/6) y+d*sin(phi-5*pi/6)]';
        B_2 = [x+d*cos(phi-pi/6) y+d*sin(phi-pi/6)]';
        B_3 = [x+d*cos(phi+pi/2) y+d*sin(phi+pi/2)]';

        P = [x y]';
        
        rho_1 = norm(B_1 - A_1);
        rho_2 = norm(B_2 - A_2);
        rho_3 = norm(B_3 - A_3);
        
        if minJointLimitEnabled
            min_rho = [min_rho, min([rho_1, rho_2, rho_3])];
        end
        if maxJointLimitEnabled
            max_rho = [max_rho, max([rho_1, rho_2, rho_3])];
        end
        
        if condEnabled
            % jacobian
            E = [0, -1; 1, 0];

            v_1 = (B_1-A_1)/norm(B_1-A_1);
            v_2 = (B_2-A_2)/norm(B_2-A_2);
            v_3 = (B_3-A_3)/norm(B_3-A_3);

            % paralel Jacobian

            A = [(v_1)',-(v_1)'*E*(P-B_1); ...
                 (v_2)',-(v_2)'*E*(P-B_2); ...
                 (v_3)',-(v_3)'*E*(P-B_3)];

            M = isnan(A);

            if(det(A) == 0 || any(M(:)))
              disp('Paralel singularity');
              condList = [condList 1e-16];
            else
              % serial Jacobian 
              B = eye(3);

              % condition
              J = pinv(A)*B;
              condList = [condList (conditionnement*cond(J)-1)];
            end
        end
      end
    end
    if condEnabled
        G = [max(condList)];
    else
        G = [];
    end
    if minJointLimitEnabled
        minRho = min(min_rho);
        if minRho < minJointLimit
            G = [G 1e-16];
        end
    end
    if maxJointLimitEnabled
        maxRho = max(max_rho);
        if maxRho > maxJointLimit
            G = [G 1e-16];
        end
    end
    Geq = [];
end