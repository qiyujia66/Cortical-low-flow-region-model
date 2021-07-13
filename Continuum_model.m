%%% The continuum model of computing low flow regions in cerebral
%%% cortex. 

tic
%% Setting computation domain geometry

set(0, 'DefaultLineLineWidth', 2); % set the default line width to 2
L = 600; % box size
H = 300; % box height
Length = 300; % length of the vessels
Z1 = 0; % top and bottom positions of the computational domain
Z2 = H;
z1 = (H - Length)/2; % top and bottom positions for the vessels
z2 = z1 + Length;


delta_y = 1/50*Length; % spacing of points along the vessels
v = z1+delta_y : delta_y : z2; % location of points along the vessels

w = numel(v); % number of discretization points needed per vessel

% The available locations of penetrating vessel arrangments
X1 = [5; 25; 35; 55; 10; 30; 40; 60; 5; 25; 35; 55;  10; 30; 40; 60;  5; 25; 35; 55;  10; 30; 40; 60; ]/60*L;
Y1 =  [5; 5; 5; 5; 15; 15; 15; 15; 25; 25; 25; 25;  35; 35; 35; 35;  45; 45; 45; 45;  55; 55; 55; 55; ]/60*L;
X2 = [15; 45; 20; 50; 15; 45; 20; 50; 15; 45; 20; 50; ]/60*L;
Y2 = [5; 5; 15; 15; 25; 25; 35; 35; 45; 45; 55; 55;]/60*L;
XXX = [X2; X1];
YYY = [Y2; Y1;];
figure;
plot(XXX, YYY, 'g*', 'linewidth', 8)
text(XXX, YYY, num2str([1:36]'),'FontSize',14)
axis([0 L 0 L])
ax = gca;
ax.FontSize = 16;
title('available locations of penetrating vessels','FontSize',14)

% Specifying which penetrating vessels are arterioles and the rest are venules 
arteriole_index = [13, 14, 2, 16, 17, 3, 19, 4, 20, 21, 5, 23, 6, 7, 26, 27, 8, 29, 30, 31, 10, 32, ...
    33, 12]; % One random arrangement
% arteriole_index = [9 5 1 34 26 18 10 6 2 36 28 20]; % An ordered honeycomb arrangment 
venule_index = setdiff(1:36, arteriole_index);

Xaaa = XXX(arteriole_index);
Yaaa = YYY(arteriole_index);
Xv = XXX(venule_index);
Yv = YYY(venule_index);

figure;
plot(Xaaa, Yaaa, 'r*', 'linewidth', 8)
text(Xaaa, Yaaa, num2str([1:numel(Xaaa)]'),'FontSize',14)
hold on
plot(Xv, Yv, 'b*', 'linewidth', 8)
% text(Xv, Yv, num2str([1:12]'),'FontSize',14)
axis([0 L 0 L])
ax = gca;
ax.FontSize = 16;
title('the penetrating vessel arrangement','FontSize',14)


ww = size([Xv; Xaaa], 1); %ww is the number of art-ven
line_x = reshape(([Xv; Xaaa]*ones(1, w))', size([Xv; Xaaa]*ones(1, w),1)*size([Xv; Xaaa]*ones(1, w),2), 1);
line_y = reshape(([Yv; Yaaa]*ones(1, w))', size([Yv; Yaaa]*ones(1, w),1)*size([Yv; Yaaa]*ones(1, w),2), 1);
line_z = reshape((v'*ones(1, ww)), size(v'*ones(1, ww),1)*size(v'*ones(1, ww),2), 1);
line = zeros(size(line_x,1)*1, 3); % xyz coordinates of art-ven
line(:,1) = line_x;
line(:,2) = line_y;
line(:,3) = line_z;

% Plotting a 3D sketch of the computation domain geometry
wwv=  size([Xv], 1); %ww is the number of art-ven
line_xv = reshape((Xv*ones(1, w))', size(Xv*ones(1, w),1)*size(Xv*ones(1, w),2), 1);
line_yv = reshape((Yv*ones(1, w))', size(Yv*ones(1, w),1)*size(Yv*ones(1, w),2), 1);
line_zv = reshape((v'*ones(1, wwv)), size(v'*ones(1, wwv),1)*size(v'*ones(1, wwv),2), 1);
wwa=  size([Xaaa], 1); %ww is the number of art-ven
line_xa = reshape((Xaaa*ones(1, w))', size(Xaaa*ones(1, w),1)*size(Xaaa*ones(1, w),2), 1);
line_ya = reshape((Yaaa*ones(1, w))', size(Yaaa*ones(1, w),1)*size(Yaaa*ones(1, w),2), 1);
line_za = reshape((v'*ones(1, wwa)), size(v'*ones(1, wwa),1)*size(v'*ones(1, wwa),2), 1);
figure;
plot3(line_xv, line_yv, line_zv,'*')
hold on
plot3(line_xa, line_ya, line_za,'r*')
axis([0 L 0 L 0 H])
title('a 3D sketch of the computation domain geometry','FontSize',14)

%% GG0 is Green's function of every source evaluated at other source locations

tic
n = 0 : 1 : 50;
m = n;

k = 2*pi*n/L; % Fourier modes
j = 2*pi*m/L;

S = sqrt(k.^2 + j'.^2); % sum of squares of wave numbers


GG0 = zeros( size(line,1),  size(line,1)); % Greens function for all numerical points

for I = 1 : size(line,1) % loop over every discretization point
    
    z0 = line(I, 3); % z-coordinate of source at this point
    
    a = sinh((z0 - H)*S);
    b = sinh(z0*S);
    c = cosh((z0 -H)*S);
    d = cosh(z0*S);
    
    e = -(a - c./d.*b).*S;
    
    B = 1./e/L^2;
    A = B.*c./d;
        
    B(1,1) = 0;
    A(1,1) = 0;
    
    
    for V = 1 : I
        if line(I, 3) > line(V, 3)
            GG0(I, V) = sum(2*(cos(k*(line(V, 1)-line(I, 1)) + j'*(line(V, 2)-line(I, 2))) + cos(k*(line(V, 1)-line(I, 1)) - j'*(line(V, 2)-line(I, 2))) ).*A.*cosh(S*line(V, 3)), 'all');
        else
            GG0(I, V) = sum(2*(cos(k*(line(V, 1)-line(I, 1)) + j'*(line(V, 2)-line(I, 2))) + cos(k*(line(V, 1)-line(I, 1)) - j'*(line(V, 2)-line(I, 2))) ).*B.*cosh(S*(line(V, 3) - H)), 'all');
        end
    end
    
end

GG0 = GG0+GG0' - diag(diag(GG0));


%% G is the Green's function of every source evaluated throughout the domain

n = [0 : 1 : 100, -99 : 1 : -1];
m = n;

k = 2*pi*n/L; % Fourier modes
j = 2*pi*m/L;

S = sqrt(k.^2 + j'.^2); % sum of squares of wave numbers

zrange = linspace(0, H, 12);
f_hat = zeros(size(m,2),size(m,2),size(zrange, 2));
G_hat = cell(size(line,1),1);
G = G_hat;
for I = 1 : size(line, 1)
    
    z0 = line(I, 3);
    
    a = sinh((z0 - H)*S);
    b = sinh(z0*S);
    c = cosh((z0 - H)*S);
    d = cosh(z0*S);
    
    e = -(a - c./d.*b).*S;
    
    B = 1./e/L^2;
    A = B.*c./d;
    B(1,1) = 0;
    A(1,1) = 0;
    
    for V = 1 : size(zrange, 2)
        
        if line(I, 3) < zrange(V)
            f_hat(: ,: ,V) =  B.*cosh(S*(zrange(V) - H));
        else
            f_hat(: ,: ,V) = A.*cosh(S*zrange(V));
        end
        G_hat{I}(: ,: ,V) = exp(-1i*k*line(I, 1) - 1i*j'*line(I, 2)).*f_hat(: ,: ,V);
        
        G{I}(: ,: ,V)  = ifft2(G_hat{I}(: ,: ,V))*size(n,2)^2;
        
    end
    
end

%% Compute the pressure and flow of the domain providing boundary conditions: pressure on top of arterioles 100, venules 50.

for beta = 16 % which number vessel will change its conductance
    pressure = zeros(size(f_hat));
    
    Xa = Xaaa([1:end]);
    Ya = Yaaa([1:end]);
    
    ww = size(Xv, 1)  + size(Xa, 1);
    G0 = GG0(1:ww*w, 1:ww*w);
    Laplacian = diag(-2*ones(w, 1)) + diag(ones(w-1, 1), 1) + diag(ones(w-1, 1), -1);
    Laplacian(w, w) = -1;
    Laplaciancell = repmat({Laplacian}, ww, 1);
    Laplacian = blkdiag(Laplaciancell{:});
    
    K = 0.8; % "Conductivity" of capillary bed
    K1 = 9000; % conductivity of penetrating arteriole
    K2 = 9000; % conductivity of penetrating venule
    p1 = 100; % pressure on top of arterioles
    p2 = 50;% pressure on top of venules
    
    l = size(G0, 1);
    
    p0 = zeros(l, 1);
    i = 1 : w : size(Xv, 1)*w+1-w;
    p0(i) = p2;
    
    p0(size(Xv, 1)*w+1 : w : size([Xa; Xv], 1)*w+1-w) = p1;
    
    z = repmat(ones(w, 1), ww, 1);
    
    % Change the conductance of vessel beta
    Kappa = [K1*ones(size(Xv, 1)*w, 1);K2*ones(size(Xa, 1)*w, 1)];
    % if blocking vessel beta then Kappa((beta-1+size([Xv], 1))*w+1:(beta+size([Xv], 1))*w) = 0; if dilating vessle beta then Kappa((beta-1+size([Xv], 1))*w+1:(beta+size([Xv], 1))*w) = 10* Kappa((beta-1+size([Xv], 1))*w+1:(beta+size([Xv], 1))*w);
    
%     Blocking arteriole beta
%     Kappa((beta-1+size(Xv, 1))*w+1:(beta+size(Xv, 1))*w) = 0;  %10* Kappa((beta-1+size(Xv, 1))*w+1:(beta+size(Xv, 1))*w);

%     Blocking venule beta
%     Kappa((beta-1)*w+1:beta*w) = 0;  %10* Kappa((beta-1+size(Xv, 1))*w+1:(beta+size(Xv, 1))*w);
  
    A = diag(Kappa)/K/delta_y* G0 * Laplacian - diag(z);
    b = -  diag(Kappa)/K/delta_y* G0 * p0;
    
    A = [A, ones(size(A,1), 1)];
    A = [A; Kappa'.*ones(1, size(A,1))*Laplacian,0];
    
    b = [b; -sum(Kappa.*p0)];
    
    A(w:w:ww*w, :) = 0;
    for i = w:w:ww*w
        A(i, i) = 1;
        A(i, i-1) = -1;
    end
    b(w:w:ww*w) = 0;
    
    pp = A\b;
    p = pp(1:end-1);
    
    % alpha is the strength of sources and sinks
    alpha = diag(Kappa)*(Laplacian*p + p0)/delta_y^2;
    
    for I = 1 : l
        for V = 1 : size(zrange, 2)
            pressure(: ,: ,V) = pressure(: ,: ,V) + 1/K*G{I}(: ,:, V)*alpha(I)*delta_y;
        end
    end
    pressure =  real(pressure)  + pp(end); % pp(end) is the zero wavelenth solution; pressure is the pressure solution.
   
    % Take gradient of pressure to obtain flux field
    [fy,fx,fz] = gradient(pressure, L/(size(pressure,2)-1), L/(size(pressure,1)-1), L/(size(pressure,3)-1));
    fx1 = K*fx;
    fy1 = K*fy;
    fz = K*fz;
    for V = 1 : size(zrange, 2)
        fx(: ,: ,V) = fy1(: ,: ,V)';
        fy(: ,: ,V) = fx1(: ,: ,V)';
    end
    
    flux_magnitude = sqrt (fx.^2 + fy.^2 + fz.^2);
    
    % Calculate the ratio of low flow regions and median flow along depth
    for Z = 1:1:10
        flux_magnitude1 = flux_magnitude(:, :, Z);
        [a3, a4] = find(flux_magnitude1 < 0.4*median(flux_magnitude1,'all'));
        median_flux(Z) =median(flux_magnitude1,'all');
        AA = meshgrid(1:200,1:200);
        B = AA';
        aa = [AA(1:end)', B(1:end)'];
        a1 = aa(:,1);
        a2 = aa(:,2);
        bb = [a3, a4];
        i = find(ismember(aa,bb,'rows'));
        ii = setdiff(1:numel(a1), i);        
        bad_ratio(Z) = numel(a3)/40000; % evaluate the low flow region ratio along the depth.
    end
    
    flux_magnitude1 = flux_magnitude(:, :, 6);
    [a3, a4] = find(flux_magnitude1 < 0.4*median(flux_magnitude1,'all'));
    
    % Plot low flow ratios
    figure; scatter(a3/200*L, a4/200*L)
    axis off ;
    hold on
    plot(Xa, Ya, 'r*', 'linewidth', 8)
    text(Xaaa, Yaaa, num2str([1:numel(Xaaa)]'),'FontSize',14)
    hold on
    plot(Xv, Yv, 'b*', 'linewidth', 8)
    axis([0 L 0 L])
    axis equal
    title('Low flow regions','FontSize',14)
    
end


toc