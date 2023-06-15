clear all

% solution parameters
G = 1;  % slip length

% simulation parameters
x1 = -10; x2 = 10;          % spatial domain
t1 = -50; t2 = 30;          % start and end times
dt = 0.005;                 % time step
Nx = 401;                   % number of x points
Con_tol = 1e-6;             % global convergence tolerance
touchdown = 0.2;            % termination point
HilbertCoefficients = 5*Nx; % number of Hilbert Coefficients for Hilbert transform

% solution spaces
dx = (x2-x1)/(Nx-1);    
x = x1:dx:x2;
t = t1:dt:t2;
F = zeros(Nx,1);
P = zeros(Nx,1);
J = zeros(Nx,Nx);       % lubrication eqn: P coefficients
HilbPx = zeros(Nx,1);   % Hilbert transform of dP/dX
Px = zeros(Nx,1);       % dP/dX
f = zeros(Nx,1);        % lubrication eqn: JP = f

Nt = length(t); 

% initial conditions
F(:,1) = x.^2/2-t(1);

% start time advancement
for m = 1:Nt-1 
    
    % guess from previous time step
    P(:,m+1) = P(:,m);
    F(:,m+1) = F(:,m);
    F(1,m+1) = F(1,m) - dt; F(Nx,m+1) = F(Nx,m)-dt;
    
    % start iteration 
    Perror = 10; Ferror = 10;
    while Perror > Con_tol || Ferror > Con_tol
    
        Pprevious = P(:,m+1);
        Fprevious = F(:,m+1);
            
        % solve for P  
        for i = 2:Nx-1
            J(i,i-1) = -3*(F(i,m+1))^2*(F(i,m)+2*G)^2*(F(i+1,m+1)-F(i-1,m+1))/((F(i,m)+G)^2*4*dx^2)+(F(i,m+1))^3*(F(i,m)+4*G)/((F(i,m)+G)*dx^2);
        end
        for i = 1:Nx
            J(i,i) = -2*(F(i,m+1))^3*(F(i,m)+4*G)/((F(i,m)+G)*dx^2);
        end
        for i = 2:Nx-1
            J(i,i+1) =  3*(F(i,m+1))^2*(F(i,m)+2*G)^2*(F(i+1,m+1)-F(i-1,m+1))/((F(i,m)+G)^2*4*dx^2)+(F(i,m+1))^3*(F(i,m)+4*G)/((F(i,m)+G)*dx^2);
        end        
        for i = 2:Nx-1
            f(i) = 12*(F(i,m+1)-F(i,m))/dt;
        end
        
        P(2:Nx-1,m+1) = J(2:Nx-1,2:Nx-1) \ f(2:Nx-1);
        
        % relax pressure with previous time step
        P(2:Nx-1,m+1) = 0.5*(P(2:Nx-1,m+1)+P(2:Nx-1,m));
        
        % solve for F, with F(:,0) = F(:,1) + dt
        if m == 1 
            w = F(:,m)+dt;
        else
            w = F(:,m-1);
        end
        for i = 2:Nx-1
            Px(i) = (P(i+1,m+1)-P(i-1,m+1))/(2*dx);
        end
        HilbPx = imag(hilbert(Px,HilbertCoefficients));
        for i = 1:Nx          
            F(i,m+1) = dt^2*HilbPx(i) + 2*F(i,m) - w(i);
        end
  
        % calculate errors
        Perror = rme(P(:,m+1),Pprevious,Nx);
        Ferror = rme(F(:,m+1),Fprevious,Nx);
        
    end
    
    % terminate code if min(F) < touchdown
    if min(F(:,m+1)) < touchdown
        break
    end
     
end

% relative max error function
function [ r ] = rme(x,y,N)

    r = max(abs(x-y))/(sum(abs(x))/N);

end