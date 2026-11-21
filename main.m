%  //////////////////////////////////////////////////////////
%  /  Matlab's pseudo-code example for programming the DVM  /
%  /  220024 Aerodinàmica                                   /
%  /  ESEIAAT_UPC                                           /
%  //////////////////////////////////////////////////////////

% 
% Input parameters (all the distances are dimensionless)
% 

M = 64;    % number of panels

alfa = 4*2*pi/360; % angle of attack of analysis

cf_c = 0.25;
xh = 1 - cf_c;
eta_f = 10*2*pi/360; % flap deflection angle

ct_c = 0;
xt = 1 - ct_c;           
eta_t = 0;     


designacio = 230;
L = floor(designacio/100);
P = floor(mod(designacio, 100)/10);
Q = mod(designacio, 10);

xcm = 0.05*P;

if Q == 0
    % Standard series
    m_old = 0.1;
    diff_m = 1;
    while diff_m > 0.0001
        m = xcm + m_old * sqrt(m_old / 3);
        diff_m = abs(m - m_old);
        m_old = m;
    end

    qm = (3*m - 7*m^2 + 8*m^3 - 4*m^4)/sqrt(m - m^2) - 1.5*(1 - 2*m)*((pi/2) - asin(1 - 2*m));
    k1 = (6 * (0.15 * L)) / qm;

elseif Q == 1
    % Réflex series
    switch P
        case 2 
            m = 0.1300;
            k1 = 51.640;
        case 3 
            m = 0.2170;
            k1 = 15.793;
        case 4 
            m = 0.3180;
            k1 = 6.520;
        case 5 
            m = 0.4310;
            k1 = 3.191;
        otherwise
            error('Només es disposa de valors tabulats per a P=2, 3, 4 o 5 en perfils reflex.');
    end
    k21 = (3*(m - xcm)^2 - m^3) / (1 - m)^3;
end



% 
% (1) Geometry discretization
% 

N = M + 1 ;  % number of points along the chord

coord = zeros(2,N) ;   % allocate storage space for (x,z) coordinates   
beta = linspace(0, pi, N);
                                                                                                             
for i = 1:N   % calculate coordinates of the mean line
    
    coord(1,i) = 0.5 * (1 - cos(beta(i)));   % calculate x position using cosine distribution
    x = coord(1,i);
    
    if ( coord(1,i) < m )  % z coordinate according to NACA5 mean line functions (split point 'm')
        if Q == 0
            coord(2,i) = (k1/6) * (x^3 - 3*m*x^2 + m^2*(3-m)*x);
        else
            coord(2,i) = (k1/6) * ((x-m)^3 - k21*(1-m)^3*x - m^3*x + m^3);
        end
    else
        if Q == 0
            coord(2,i) = (k1*m^3 / 6) * (1 - x);
        else
            coord(2,i) = (k1/6) * (k21*(x-m)^3 - k21*(1-m)^3*x - m^3*x + m^3);
        end
    end
    
end


 % calculate z position of the hinges (flap and tab)
if ( xh < m )  % hinge z position according to NACA5 mean line functions
    if Q == 0
        zh = (k1/6) * (xh^3 - 3*m*xh^2 + m^2*(3-m)*xh);
    else
        zh = (k1/6) * ((xh-m)^3 - k21*(1-m)^3*xh - m^3*xh + m^3);
    end
else
    if Q == 0
        zh = (k1*m^3 / 6) * (1 - xh);
    else
        zh = (k1/6) * (k21*(xh-m)^3 - k21*(1-m)^3*xh - m^3*xh + m^3);
    end
end

if ( xt < m )  % tab hinge z position
    if Q == 0
        zt = (k1/6) * (xt^3 - 3*m*xt^2 + m^2*(3-m)*xt);
    else
        zt = (k1/6) * ((xt-m)^3 - k21*(1-m)^3*xt - m^3*xt + m^3);
    end
else
    if Q == 0
        zt = (k1*m^3 / 6) * (1 - xt);
    else
        zt = (k1/6) * (k21*(xt-m)^3 - k21*(1-m)^3*xt - m^3*xt + m^3);
    end
end

xt_new = xt; 
zt_new = zt;

if ( abs(eta_f) > 0.0 )  %  modify camber line coordinates to simulate flap deflection
    
    rot_f = [cos(-eta_f) -sin(-eta_f); sin(-eta_f) cos(-eta_f)];
    
    for i = 1:N
        if ( coord(1,i) > xh )  % rotate point an angle eta_f about (xh,zh)
            dx = coord(1,i) - xh;
            dz = coord(2,i) - zh;
            new_vec = rot_f * [dx; dz];
            coord(1,i) = xh + new_vec(1);
            coord(2,i) = zh + new_vec(2);
        end
    end
    
    % Update tab hinge position since it rotated with the flap
    dtx = xt - xh;
    dtz = zt - zh;
    new_t = rot_f * [dtx; dtz];
    xt_new = xh + new_t(1);
    zt_new = zh + new_t(2);
end

if ( abs(eta_t) > 0.0 )  % modify camber line coordinates to simulate tab deflection
    
    rot_t = [cos(-eta_t) -sin(-eta_t); sin(-eta_t) cos(-eta_t)];
    
    for i = 1:N
        x_orig = 0.5 * (1 - cos(beta(i))); % check original point position
        if ( x_orig > xt )  % rotate point an angle eta_t about (xt_new,zt_new)
            dx = coord(1,i) - xt_new;
            dz = coord(2,i) - zt_new;
            new_vec = rot_t * [dx; dz];
            coord(1,i) = xt_new + new_vec(1);
            coord(2,i) = zt_new + new_vec(2);
        end
    end
end


  
pchord = zeros(M) ;   % allocate storage space
pnorm = zeros(2,M) ;
ptang = zeros(2,M) ;
xvort = zeros(2,M) ;
xcont = zeros(2,M) ;

for i = 1:M  % calculate panels geometry
    
    x1 = coord(:,i)  ;  % panel nodes 
    x2 = coord(:,i+1) ;
    
    % calculate panel's chord, normal and tangent vectors
    pchord(i) = norm(x2 - x1);
    tx = (x2(1) - x1(1)) / pchord(i);
    tz = (x2(2) - x1(2)) / pchord(i);
    
    ptang(1,i) = tx ; % components of the tangent vector
    ptang(2,i) = tz ; 
    pnorm(1,i) = -tz ; % components of the normal vector
    pnorm(2,i) = tx  ;
    
    % calculate position of the vortex and control point
    xvort(:,i) = x1(:) + 0.25*pchord(i)*ptang(:,i) ;
    xcont(:,i) = x1(:) + 0.75*pchord(i)*ptang(:,i) ;
    
end

% plot the camber line with the corresponding points and normal vectors to
% check everything is ok. You can use the Matlab's function quiver(...) to plot normal vectors


% 
% (2) Assembly and solution of the equations system
% 

A = zeros(M,M) ;  % allocate variables
RHS = zeros(M,1) ;  

for i = 1:M  % scan control points
    
    xi = xcont(:,i)   ;  % control point's position
    ni = pnorm(:,i)   ;  % normal vector
    
    for j = 1:M  % scan vortices
        
        xv = xvort(:,j)  ;  % vortex position
        rx = xi(1) - xv(1);
        rz = xi(2) - xv(2);
        r2 = rx^2 + rz^2;
        
        u = rz / (2 * pi * r2)  ; % induced velocities
        w = -rx / (2 * pi * r2) ;
            
        A(i,j) = u * ni(1) + w * ni(2) ;  % influence coefficient
            
    end
    
    RHS(i) = -(cos(alfa) * ni(1) + sin(alfa) * ni(2))  ;  % projection of Uinf onto ni(:,i)
        
end

G = A \ RHS  ;  % circulation

% 
% (3) Loads calculation
% 

CL = 0.0 ; 
CMLE = 0.0  ;  
CH = 0.0 ;

for i = 1:M   % scan panels' vortices
    
    delta_CL = 2 * G(i)   ; % panel's lift using K-J
    delta_CM = - delta_CL * xvort(1,i)  ; % moment
        
    CL = CL + delta_CL ;   % accumulate panels' contributions
    CMLE = CMLE + delta_CM  ;
    
    % Flap moment about the hinge
    x1_orig = 0.5 * (1 - cos(beta(i)));
    if x1_orig >= xh - 1e-5
        dx_h = xvort(1,i) - xh; 
        delta_CH = - delta_CL * dx_h / (cf_c^2);
        CH = CH + delta_CH;
    end
    
end
fprintf('CL = %.6f\nCMLE = %.6f\nCH = %.6f\n', CL, CMLE, CH);