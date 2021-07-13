%%% The synthetic network model of computing low flow regions in cerebral
%%% cortex. 
tic
%% Set computation domain geometry

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
 
i = [13, 14, 2, 16, 17, 3, 19, 4, 20, 21, 5, 23, 6, 7, 26, 27, 8, 29, 30, 31, 10, 32, ...
    33, 12]; % One random arrangement
% i = [9 5 1 34 26 18 10 6 2 36 28 20]; % An ordered honeycomb arrangment 
o = setdiff(1:36, i);

figure;
plot(Xaaa, Yaaa, 'r*', 'linewidth', 8)
text(Xaaa, Yaaa, num2str([1:numel(Xaaa)]'),'FontSize',14)
hold on
plot(Xv, Yv, 'b*', 'linewidth', 8)
% text(Xv, Yv, num2str([1:12]'),'FontSize',14)
axis([0 L 0 L])
ax = gca;
ax.FontSize = 16;

% figure;
% voronoi(XXX, YYY)
% axis equal
% axis([0 L 0 L])


xs = rand(1, 5000)*L + 0;% rand*10 + 15;
ys = rand(1, 5000)*L + 0; %rand*10 + 25;
zs = rand(1, 5000)*Length + zmin;
XS = [xs', ys', zs'];

% All sample connecting points
figure
plot3(xs, ys, zs, '*')

distance = zeros(numel(xs), numel(XXX));
for i = 1:numel(xs)
    for j = 1:numel(XXX)
        distance(i,j) = norm([XS(i,1)-XXX(j), XS(i,2)-YYY(j)]);
    end
end

[B,I] = sort(distance,2);
nearst_2trees = I(:,1:2);
penetrating_vessels_position = [XXX, YYY];

%% Build the synthetic network

p0 = 100;
pt = 50;
delta_p = p0 - pt;
f = 1; % outflow flux
l = 1;
% 1 is the top boundary point
adj_tree = [];
position_tree = [];
Radius_tree = [];
for tree = 1 : size(penetrating_vessels_position, 1)
    
    if tree <= numel(Xaaa)
        f = 1;
    else
        f = 1;%numel(Xaaa)/numel(Xv);
    end

    
    position = [XXX(tree), YYY(tree), H; XXX(tree), YYY(tree), 0];
    adj = [1;2];
    [i1,~]=find(nearst_2trees==tree);
    xr = xs(i1);
    yr = ys(i1);
    zr = zs(i1);
    
    for iii = 1 : numel(xr)
        
        % new random vertex point
        Xr = [xr(iii), yr(iii), zr(iii)];
        % distance is the length from new point to other points
        distance = vecnorm(Xr-position,2,2);
        
        position = [position; xr(iii),yr(iii),zr(iii)];
        
        min_distance = zeros(size(adj,2), 1);
        
        for AddOnEdge = 1:size(adj,2)
            
            x1 = position(adj(1, AddOnEdge), 1);
            x2 = position(adj(2, AddOnEdge), 1);
            y1 = position(adj(1, AddOnEdge), 2);
            y2 = position(adj(2, AddOnEdge), 2);
            z1 = position(adj(1, AddOnEdge), 3);
            z2 = position(adj(2, AddOnEdge), 3);
            length = norm([x1-x2, y1-y2, z1-z2]);
            X1 = [x1, y1, z1];
            X2 = [x2, y2, z2];
            t = (- dot(X1-Xr, X2-X1))/length^2;
            
            % min_distance is the minimum distance from the new point to an
            % edge
            if t>0 && t<1
                min_distance(AddOnEdge) = norm(cross(Xr-X1, Xr-X2))/norm(X2-X1);
            else
                min_distance(AddOnEdge) = min([norm(Xr-X1), norm(Xr-X2)]);
            end
        end
        
        % if the new point is within 1 from another point or is closer with
        % another point than any edge, connect the new point to the existing
        % point
        if min(distance)< 10 || min(distance)< min(min_distance)
            [~, I] = min(distance);
            adj = [adj, [max(adj,[],'all')+1; I]];
            tot_adjac = [adj flipud(adj)];
            adjacency = sparse(tot_adjac(1,:), tot_adjac(2,:), ones(size(tot_adjac, 2), 1));
            G = graph(adjacency);
            
            all_nodes = unique( reshape(G.Edges.EndNodes, numel(G.Edges.EndNodes), 1) );
            
            [C,~,ic] = unique( reshape(G.Edges.EndNodes, numel(G.Edges.EndNodes), 1) );
            a_counts = accumarray(ic,1);
            value_counts = [C, a_counts];
            terminal_nodes = setdiff(value_counts((value_counts(:,2) == 1), 1), 1);
            terminal_edges= zeros(1, numel(terminal_nodes));
            for i = 1:numel(terminal_nodes)
                terminal_edges(i) = unique([find(adj(1,:)==terminal_nodes(i)); find(adj(2,:)==terminal_nodes(i))]);
            end
            
            % EP is the edge path form terminal nodes to inlet node
            EP = cell(1, numel(terminal_nodes));
            for i = 1 : numel(terminal_nodes)
                [P,d,edgepath] = shortestpath(G, terminal_nodes(i), 1);
                EP{i} = edgepath;
            end
            
            % 'flux(:,2)' is the flux of the edge 'flux(:,1)'
            [C,~,ic] = unique(cell2mat(EP));
            a_counts2 = accumarray(ic,1);
            flux = [C', a_counts2];
            
            r_list = flux(:,2).^(1/3);
            
            R = r_list;
            
        else
            
            
            [~,I] = min(min_distance);
            AddOnEdge = I;
            
            x1 = position(adj(1, AddOnEdge), 1);
            x2 = position(adj(2, AddOnEdge), 1);
            y1 = position(adj(1, AddOnEdge), 2);
            y2 = position(adj(2, AddOnEdge), 2);
            z1 = position(adj(1, AddOnEdge), 3);
            z2 = position(adj(2, AddOnEdge), 3);
            llength= norm([x1-x2, y1-y2, z1-z2]);
            
            
            
            if llength < 5
                N1 = 5;
            else
                N1 = floor(llength);
            end
            
            xb = linspace(x1,x2,N1);    xb = xb(2:end-1);
            yb = linspace(y1,y2,N1);    yb = yb(2:end-1);
            zb = linspace(z1,z2,N1);    zb = zb(2:end-1);
            
            xd = linspace(x1, Xr(1), N1);     xd = xd(2:end-1);
            yd = linspace(y1, Xr(2), N1);     yd = yd(2:end-1);
            zd = linspace(z1, Xr(3), N1);     zd = zd(2:end-1);
            
            xc = zeros(numel(xd) + (numel(xd)+1)*numel(xd)/2, 1);
            for i = 1 : numel(xd)
                xc((i+1)*i/2 : (i+1)*i/2+i) = linspace(xb(i), xd(i), 2+i-1);
            end
            yc = zeros(numel(xd) + (numel(xd)+1)*numel(xd)/2, 1);
            for i = 1 : numel(xd)
                yc((i+1)*i/2 : (i+1)*i/2+i) = linspace(yb(i), yd(i), 2+i-1);
            end
            zc = zeros(numel(xd) + (numel(xd)+1)*numel(xd)/2, 1);
            for i = 1 : numel(xd)
                zc((i+1)*i/2 : (i+1)*i/2+i) = linspace(zb(i), zd(i), 2+i-1);
            end
            
%             [~, ~, ib1] = intersect(xb, penetrating_vessels_position(:,1));
            ib1 = find(penetrating_vessels_position(:,1) == xb(1));
            ib2 = find(penetrating_vessels_position(:,2) == yb(1));
            
            if isempty(intersect(ib1, ib2)) == 1 % if the original vessel is not a penetrating vessel
                poslist = [xc, yc, zc]; % Samples points from triangular area for braching points
            else
                poslist = [xb', yb', zb']; % Samples points from original vessel for braching points if original vessel is a penetrating vessel
            end
            
            
            % when adding an edge, change 1 edge into 3
            adj = [adj(:, 1:AddOnEdge-1), [adj(1,AddOnEdge) adj(2,AddOnEdge) max(adj,[],'all')+1; ...
                max(adj,[],'all')+2   max(adj,[],'all')+2   max(adj,[],'all')+2],   adj(:, AddOnEdge+1 : end)];
            
            tot_adjac = [adj flipud(adj)];
            adjacency = sparse(tot_adjac(1,:), tot_adjac(2,:), ones(size(tot_adjac, 2), 1));
            G = graph(adjacency);
            
            all_nodes = unique( reshape(G.Edges.EndNodes, numel(G.Edges.EndNodes), 1) );
            
            [C,~,ic] = unique( reshape(G.Edges.EndNodes, numel(G.Edges.EndNodes), 1) );
            a_counts = accumarray(ic,1);
            value_counts = [C, a_counts];
            terminal_nodes = setdiff(value_counts((value_counts(:,2) == 1), 1), 1);
            terminal_edges= zeros(1, numel(terminal_nodes));
            for i = 1:numel(terminal_nodes)
                terminal_edges(i) = unique([find(adj(1,:)==terminal_nodes(i)); find(adj(2,:)==terminal_nodes(i))]);
            end
            middle_nodes = setdiff(all_nodes, [terminal_nodes;1]);
            
            % EP is the edge path form terminal nodes to inlet node
            EP = cell(1, numel(terminal_nodes));
            for i = 1 : numel(terminal_nodes)
                [P,d,edgepath] = shortestpath(G, terminal_nodes(i), 1);
                for j = 1 : numel(P)-1
                    EP{i}(j) = [find(adj(1,:) == P(j) & adj(2,:) == P(j+1)), find(adj(2,:) == P(j) & adj(1,:) == P(j+1))];
                end
            end
            
            % 'flux(:,2)' is the flux of the edge 'flux(:,1)'
            [C,~,ic] = unique(cell2mat(EP));
            a_counts2 = accumarray(ic,1);
            flux = [C', f*a_counts2];
            
            
            total_volume = zeros(numel(xb),1);
            flux_CV = zeros(numel(xb),1);
            kappa_list = zeros(numel(xb), size(G.Edges.EndNodes, 1));
            options = optimoptions('fsolve','SpecifyObjectiveGradient',false);
            r_list = flux(:,2).^(1/3);
            
            for j = 1 : numel(xb)
                
                new_position = [position; poslist(j,:)];
                l = vecnorm(new_position(G.Edges.EndNodes(:,1),:) - new_position(G.Edges.EndNodes(:,2),:), 2, 2);
                
                kappa_list(j,:) = r_list.^4./l;
                kappaij = sparse(tot_adjac(1,:), tot_adjac(2,:),[kappa_list(j,:) kappa_list(j,:)]);
                N = size(kappaij,2);
                Laplacian = spdiags(sum(kappaij,2),0,N,N)-kappaij;
                p_pts = [1; terminal_nodes];
                p_num = size(p_pts, 1);
                p_magnitude = [p0; pt*ones(size(terminal_nodes, 1), 1)];
                
                Laplacian(p_pts,:) = sparse(1:p_num,p_pts,ones(p_num,1),p_num,N);
                F1 = sparse(N,1);
                F1(p_pts) = p_magnitude;
                
                pN = Laplacian\F1;
                Flux = (pN(adj(1,:)) - pN(adj(2,:))).*kappa_list(j,:)';
                
                terminal_flux = abs(Flux(terminal_edges));
                flux_CV(j) = std(terminal_flux)/mean(terminal_flux);
                total_volume(j) = sum(l.*r_list.^2);
                
            end
            
            [~, I] = min(total_volume);
            new_position(end, :)  = poslist(I,:);
            %         Kappa = kappa_list(I,:);
            R = r_list;
            position = new_position;
        end
        G.Edges.Weight = R;
        x = position(:,1);
        y = position(:,2);
        z = position(:,3);
        LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
        
        end
    
    adj_tree{tree} = adj;
    position_tree{tree} = position;
    Radius_tree{tree} = R;
    
    
end

%% Put trees together to a connected network
adj_tree1= adj_tree;
edge_number = size(adj_tree1{1}, 2);
for i = 2 : numel(adj_tree1)
    adj_tree1{i} = adj_tree1{i} + max(adj_tree1{i-1},[],'all');
    edge_number = edge_number + size(adj_tree1{i}, 2);
end

full_adj = cell2mat(adj_tree1);
full_pos = cell2mat(position_tree');
full_radius = cell2mat(Radius_tree');
full_length = vecnorm(full_pos(full_adj(1,:),:) - full_pos(full_adj(2,:),:),2,2);
edge_center = (full_pos(full_adj(1,:),:) + full_pos(full_adj(2,:),:))/2;
% kappa_list = kappa_list(randperm(numel(kappa_list)));
U = unique(full_pos,'rows');
R = cell(size(full_pos,1),1);
for k = 1:size(U,1)
     R{k}=find(ismember(full_pos ,U(k,:), 'rows'));
end
Rows=R(cell2mat(cellfun(@(x)numel(x)>=2,R,'un',0))); %Duplicate vertex position rows index
for i = 1 : numel(Rows)
%     [r,c]=find(full_adj == Rows{i}(2));
    full_adj(full_adj == Rows{i}(2)) = Rows{i}(1);
end
 
%% Calculate pressure and flux field of network

[full_adj1, ia, ic] = unique(sort(full_adj)','rows','stable');
full_radius1 = full_radius(ia);
full_adj1 = full_adj1';
full_length1 =  full_length(ia);
edge_center1 = edge_center(ia,:);
resistance_list1 = 8*0.0091/pi*full_length1./(full_radius1/2).^4*4.*(1-0.863*exp(-full_radius1/14.3)+27.5*exp(-full_radius1/0.351));
kappa_list1 = 1./resistance_list1;

% 
[~, I] = sort(full_adj1(1,:)*10^6 + full_adj1(2,:));
full_adj1 = full_adj1(:,I);
full_radius1 = full_radius1(I);
full_length1 = full_length1(I);
edge_center1 = edge_center1(I,:);
kappa_list1 = kappa_list1(I);


tot_adjac = [full_adj1 flipud(full_adj1)];
adjacency = sparse(tot_adjac(1,:), tot_adjac(2,:), ones(size(tot_adjac, 2), 1));
G = graph(adjacency);
G.Edges.Weight = full_radius1;
nodes = unique([G.Edges.EndNodes(:,1); G.Edges.EndNodes(:,2)]);
x = full_pos(1:max(nodes),1);
y = full_pos(1:max(nodes),2);
z = full_pos(1:max(nodes),3);
LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
figure
plot(G,'XData',x,'YData',y,'ZData',z, 'LineWidth', LWidths,'NodeLabel',{},'EdgeColor', 'b', 'Marker','none')%,'edgealpha',1
axis([0 L 0 L])
axis equal

all_nodes = [1:max(full_adj1, [], 'all')]';
arteriole_top_nodes = all_nodes ((ismember(full_pos(all_nodes,:) ,[Xaaa,Yaaa,H*ones(numel(Xaaa),1)], 'rows')));
venule_top_nodes = all_nodes ((ismember(full_pos(all_nodes,:) ,[Xv, Yv,H*ones(numel(Xv),1)], 'rows')));
arteriole_nodes = all_nodes ((ismember(full_pos(all_nodes,1:2) ,[Xaaa,Yaaa], 'rows')));
venule_nodes = all_nodes ((ismember(full_pos(all_nodes,1:2) ,[Xv, Yv], 'rows')));
ia1 = find(ismember(full_adj1(1,:), [arteriole_nodes;venule_nodes]));
ib1 = find(ismember(full_adj1(2,:), [arteriole_nodes;venule_nodes]));
artvein_edges = intersect(ia1,ib1)';
ia2 = find(ismember(full_adj1(1,:), [arteriole_nodes]));
ib2 = find(ismember(full_adj1(2,:), [arteriole_nodes]));
art_edges = intersect(ia2,ib2)';
ia3 = find(ismember(full_adj1(1,:), [venule_nodes]));
ib3 = find(ismember(full_adj1(2,:), [venule_nodes]));
ven_edges = intersect(ia3,ib3)';

capillary = setdiff(1:numel(kappa_list1), artvein_edges);


all_nodes_count = unique( reshape(G.Edges.EndNodes, numel(G.Edges.EndNodes), 1) );

[C,~,ic] = unique( full_adj1(1:end) );
a_counts = accumarray(ic,1);
value_counts = [C', a_counts];
terminal_nodes = setdiff(value_counts((value_counts(:,2) == 1), 1), 1);

zeroflux_nodes = setdiff(terminal_nodes, terminal_nodes(( ismember(full_pos(terminal_nodes,1), XXX) + ismember(full_pos(terminal_nodes,2), YYY) + ismember(full_pos(terminal_nodes,3), 45) == 3)));
connect_to_zeroflux_nodes = zeros(1, numel(zeroflux_nodes));
for i = 1:numel(zeroflux_nodes)
[r,c] = find(full_adj1 == zeroflux_nodes(i));
connect_to_zeroflux_nodes(i) = full_adj1(setdiff([1,2],r), c);
end


arteriole_nodes1 = all_nodes ((ismember(full_pos(all_nodes,1:2) ,[Xaaa(2),Yaaa(2)], 'rows')));
ia = find(ismember(full_adj1(1,:), arteriole_nodes1));
ib = find(ismember(full_adj1(2,:), arteriole_nodes1));
art_edges1 = intersect(ia,ib)';
sort(full_radius1(art_edges1))
z = sort(full_pos(arteriole_nodes1,3), 'descend');
mean(kappa_list1(art_edges).*full_length1(art_edges))


kappaij = sparse(tot_adjac(1,:), tot_adjac(2,:),[kappa_list1 kappa_list1]);
N = size(kappaij,2);
Laplacian = spdiags(sum(kappaij,2),0,N,N)-kappaij;

p_pts = [arteriole_top_nodes; venule_top_nodes];
top_nodes = [arteriole_top_nodes; venule_top_nodes];
pos_art_ven = full_pos(top_nodes,:);

p_num = size(p_pts, 1);
p_magnitude = p0*ones(numel(p_pts), 1);
p_magnitude(1:numel(arteriole_top_nodes)) = [p0*ones(size(arteriole_top_nodes, 1), 1)];
p_magnitude(numel(arteriole_top_nodes)+1:numel(venule_top_nodes)+numel(arteriole_top_nodes)) = [pt*ones(size(venule_top_nodes, 1), 1)];

for i = 1 : size(zeroflux_nodes, 1)
    Laplacian(zeroflux_nodes(i),:) = sparse([1 1], [connect_to_zeroflux_nodes(i) zeroflux_nodes(i)],[-1 1], 1, N);
end
Laplacian(p_pts,:) = sparse(1:p_num,p_pts,ones(p_num,1),p_num,N);


F1 = sparse(N,1);
F1(p_pts) = p_magnitude;

pN = Laplacian\F1;
Flux = (pN(full_adj1(1,:)) - pN(full_adj1(2,:))).*kappa_list1;
% pN = full(pN(~isnan(pN)));
for i = 1 : numel(Rows)
%     [r,c]=find(full_adj == Rows{i}(2));
    pN(Rows{i}(2)) = pN(Rows{i}(1));
end

flux = full(Flux);

zz=setdiff(find(edge_center1(:,3)>0 & edge_center1(:,3)<=H), artvein_edges);
a = find(abs(flux(zz)) < 0.4*median(abs(flux(zz)),'all'));
bad_edge = zz(a);

% Plot low flow regions
figure;hold on
startv = [full_pos(full_adj1(1,:),1), full_pos(full_adj1(1,:),2), full_pos(full_adj1(1,:),3)];
endv = [full_pos(full_adj1(2,:),1), full_pos(full_adj1(2,:),2), full_pos(full_adj1(2,:),3)];
for i = 1:numel(bad_edge)
h = plot3([startv(bad_edge(i),1)';endv(bad_edge(i),1)'],...
    [startv(bad_edge(i),2)';endv(bad_edge(i),2)'],...
    [startv(bad_edge(i),3)';endv(bad_edge(i),3)']);
set(h,'color', '[0 0 1 0.5]', 'LineWidth', full_radius1(bad_edge(i)));
end
axis([0 L 0 L])
axis equal
hold on
plot(Xaaa, Yaaa, 'r*', 'linewidth', 8)
plot(Xv, Yv, 'k*', 'linewidth', 8)
view(3)
title('Low flow regions','FontSize',14)

%% Plot flux heatmap of network
map =jet(numel(Flux));
pressure_edge = pN(full_adj1(1,:)) + pN(full_adj1(2,:))/2;
% [~,~,index] = unique(pressure_edge);
[~,~,index] = unique(abs(Flux));
map = map(index,:);

G.Edges.Weight = full_radius1;
nodes = unique([G.Edges.EndNodes(:,1); G.Edges.EndNodes(:,2)]);
x = full_pos(1:max(nodes),1);
y = full_pos(1:max(nodes),2);
z = full_pos(1:max(nodes),3);

LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
figure
plot(G,'XData',x,'YData',y,'ZData',z, 'LineWidth', LWidths,'NodeLabel',{},'EdgeColor', map, 'Marker','none')%,'edgealpha',1
axis([0 L 0 L])
axis equal
hold on
plot(Xaaa, Yaaa, 'r*', 'linewidth', 8)
text(Xaaa, Yaaa, num2str([1:numel(Xaaa)]'),'FontSize',14)
hold on
plot(Xv, Yv, 'b*', 'linewidth', 8)
% text(Xv, Yv, num2str([1:12]'),'FontSize',14)
axis([0 L 0 L])
ax = gca;
ax.FontSize = 16;
title('flux heatmap of network','FontSize',14)

toc