clear all

% load in data for the pressure and vorticity at r = 1 calculated in the
% air
% also load in parameters that need to carry across from the two solutions:
% theta, dtheta, t, dt, Re
load(---,'pgr1','zeta_b','theta','dtheta','t','dt','Re')

% solution parameters
mu1 = 1;
We = 1;

% simulation parameters
dr = 0.005; % spatial discretisation step
tol = 1e-5; % convergence tolerance
rP = 1.9;   % P relaxation parameter
rV = 1.4;   % V relaxation parameter
rU = 1.6;   % U relaxation parameter

%solution spaces
r = 0:dr:1;
U = zeros(length(r),length(theta),length(t));
V = zeros(length(r),length(theta),length(t));
P = zeros(length(r),length(theta),length(t));
f1 = zeros(length(theta),length(t));           

for m = 2:length(t)
    
    % find P
    P(:,:,m) = P(:,:,m-1);
    P(length(r),:,m) = 1/We + pgr1(:,m);
    Per = 10;
    while Per > tol
        P_old = P(:,:,m);
        for i = 2:length(r)-1
            for j = 2:length(theta)-1
                P(i,j,m) = rP*((-2*r(i)^2/dr^2+1-2/dtheta^2)^(-1) * (-r(i)^2*(P(i+1,j,m)+P(i-1,j,m))/dr^2+r(i)*(P(i+1,j,m)-P(i-1,j,m))/2/dr-(P(i,j+1,m)+P(i,j-1,m))/dtheta^2)) + (1-rP)*P(i,j,m);
            end
        end
        for i = 2:length(r)
            P(i,1,m) = (4*P(i,2,m)-P(i,3,m))/3;
            P(i,length(theta),m) = (4*P(i,length(theta)-1,m)-P(i,length(theta)-2,m))/3;
        end
       Per = max(max(abs(P(:,:,m)-P_old)))/(sum(sum(abs(P(:,:,m))))/(length(r)*length(theta)));
    end
   
    % find velocities
    V(:,:,m) = V(:,:,m-1);
    U(:,:,m) = U(:,:,m-1);
    
    % begin iteration
    Uer = 10; Ver = 10;
    while Ver > tol || Uer > tol
        
        U_old = U(:,:,m);
        V_old = V(:,:,m);
        
        % calculate U
        for j = 2:length(theta)-1
            U(length(r),j,m) = rU*((4*U(length(r)-1,j,m)-U(length(r)-2,j,m)-dr*(V(length(r),j+1,m)-V(length(r),j-1,m))/dtheta)/3) + (1-rU)*U(length(r),j,m);
        end
        for i = length(r)-1:-1:2
            for j = 2:length(theta)-1
                U(i,j,m) = rU*((r(i)^2/dt+(2*r(i)^2/dr^2+2/dtheta^2)/mu1/Re)^(-1) * (r(i)*P(i,j,m)-r(i)^2*(P(i+1,j,m)-P(i-1,j,m))/2/dr+(r(i)^2*(U(i+1,j,m)+U(i-1,j,m))/dr^2-r(i)*(U(i+1,j,m)-U(i-1,j,m))/2/dr+(U(i,j+1,m)+U(i,j-1,m))/dtheta^2-(V(i,j+1,m)-V(i,j-1,m))/dtheta)/mu1/Re+r(i)^2*U(i,j,m-1)/dt)) + (1-rU)*U(i,j,m);
            end
        end
        for i = 2:length(r)
            U(i,1,m) = (4*U(i,2,m)-U(i,3,m))/3;
            U(i,length(theta),m) = (4*U(i,length(theta)-1,m)-U(i,length(theta)-2,m))/3;
        end    
        
        % calculate V
        for j = 2:length(theta)-1
            V(length(r),j,m) = rV*((3/2/dr-1)^(-1) * ((4*V(length(r)-1,j,m)-V(length(r)-2,j,m))/2/dr+mu1*zeta_b(j,m)-(U(length(r),j+1,m)-U(length(r),j-1,m))/2/dtheta)) + (1-rV)*V(length(r),j,m);
        end
        for i = length(r)-1:-1:2
            for j = 2:length(theta)-1
                V(i,j,m) = rV*((r(i)^2/dt+(2*r(i)^2/dr^2+2/dtheta^2)/mu1/Re)^(-1) * (-r(i)*(P(i,j+1,m)-P(i,j-1,m))/2/dtheta + (r(i)^2*(V(i+1,j,m)+V(i-1,j,m))/dr^2-r(i)*(V(i+1,j,m)-V(i-1,j,m))/2/dr+(V(i,j+1,m)+V(i,j-1,m))/dtheta^2+(U(i,j+1,m)-U(i,j-1,m))/dtheta)/mu1/Re+r(i)^2*V(i,j,m-1)/dt)) + (1-rV)*V(i,j,m);
            end
        end
        
        % calculate errors
        Uer = max(max(abs(U(:,:,m)-U_old)))/(sum(sum(abs(U(:,:,m))))/(length(r)*length(theta)));
        Ver = max(max(abs(V(:,:,m)-V_old)))/(sum(sum(abs(V(:,:,m))))/(length(r)*length(theta)));
   
    end   
    
    % calculate f_1
    f1(:,m) = f1(:,m-1) + dt*U(length(r),:,m)';
    
    % terminate code when f_1 gets large
    if max(abs(f1(:,m))) >= 100
        break
    end

end