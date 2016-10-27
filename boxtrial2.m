dbstop if error

%number of balls - to generalise with multiple sizes
balls = 625;
%size of balls, see above
size = 0.0125;
%number of dimensions - 2 for testing
nd = 2;

%distance between 2 balls - needs generalising for multiple balls
dist = size*2;
%precompute dist^2
dist2 = dist^2;

%generate random positions and velocities
pos = rand(balls,2).*(1-dist)+size;
v = (rand(balls,2)-0.5)./100.0;
pos = repmat(([1:25]-0.5)./25,25,1);
pos_temp = pos';
pos = [pos(:),pos_temp(:)];

%max # of steps for animation
n_step = 1000;
%time per step
t_step = 0.1;

%intialise step time, last calculation time and step index
t = 0;
t_sim = 0;
step = 0;
%initialise frame array for movie
F(n_step) = struct('cdata',[],'colormap',[]);
%initialise for last collision colours
ii=1;
jj=1;
kk=1;
ll=1;

%save initial momentum and energy
momentum = sum(v(:));
energy = sum(v(:,1).^2+v(:,2).^2);

%precompute triagnular numbers for collisions
tri = sort(linspace(1,balls-1,balls-1),'descend');
tri = fix(tri);
i_arr = NaN(sum(tri),1);
j_arr = i_arr;
tri = cumsum(tri);
tri = [0,tri];
for i=1:balls-1
    i_arr(tri(i)+1:tri(i+1)) = i;
    j_arr(tri(i)+1:tri(i+1)) = linspace(1,tri(i+1)-tri(i),tri(i+1)-tri(i))+i;
end

i_arr2=i_arr;
j_arr2=j_arr;

for i=1:nd-1
    i_arr2 = [i_arr2;i_arr+balls*i];
    j_arr2 = [j_arr2;j_arr+balls*i];
end
    
%plot initial set-up
figure(1);
hold off
scatter(pos(:,1),pos(:,2));
axis([0,1,0,1]);
drawnow

while step < n_step
    %calculate time to wall collision
    t_wall = -[(pos-size)./v,(pos-(1-size))./v];
    t_wall = reshape(t_wall,balls*nd,2);
    t_wall = max(t_wall,[],2);
    t_wall_x = t_wall(1:balls);
    min_x = min(t_wall_x(:));
    t_wall_y = t_wall((balls+1):(2*balls));
    min_y = min(t_wall_y(:));
    
    %calculate time to balls colliding
    %calculate quadratic coefficients
    a = sum((v(i_arr,:)-v(j_arr,:)).^2,2);
    b = 2*sum((v(i_arr,:)-v(j_arr,:)).*(pos(i_arr,:)-pos(j_arr,:)),2);
    c = sum((pos(i_arr,:)-pos(j_arr,:)).^2,2)-dist^2;
    %check if valid solutions
    chk = b.^2-4.*a.*c;
    chk(chk<0) = NaN;
    %calculate both quadratic solutions
    temp = [(-b+chk.^0.5)./(2.*a),(-b-chk.^0.5)./(2.*a)];
    t_col = min(temp,[],2);
    t_col(t_col<0)=NaN;
    [min_t_col,i_mtc] = min(t_col(:),[],'omitnan');
    
    j_mtc = j_arr(i_mtc);
    i_mtc = i_arr(i_mtc);
    
    %find soonest collision overall
    min_t = min([min_t_col,min_x,min_y],[],'omitnan');
        
    %plot at each time step up to collision
    if (t_sim+min_t) > (t+t_step)
        while (t_sim+min_t) > (t+t_step)
            t = t+t_step;
            t_diff = t-t_sim;
            pos_temp = pos + v.*t_diff;
            figure(1);
            hold off
            plot([0,1,1,0,0],[0,0,1,1,0])
            hold on
            scatter(pos_temp(:,1),pos_temp(:,2));
            scatter(pos_temp(ii,1),pos_temp(ii,2),'filled');
            scatter(pos_temp(jj,1),pos_temp(jj,2),'filled');
            scatter(pos_temp(kk,1),pos_temp(kk,2),'filled');
            scatter(pos_temp(ll,1),pos_temp(ll,2),'filled');
            axis([-0.1,1.1,-0.1,1.1]);
            drawnow
            step = step+1;
            F(step) = getframe;
            if numel(pos(pos<0))>0 || numel(pos(pos>1))>0
                %keyboard
            end
        end
    end
    %update sim time
    t_sim = t_sim + min_t;
    t_diff = min_t;
    
    %move positions
    pos = pos+v*t_diff;
    
    if numel(t_wall_x == min_t) > 0
        %adjust speeds for wall collision
        wh = find(t_wall_x == min_t);
        whl = find(pos(wh,1) < 0.5);
        whr = find(pos(wh,1) > 0.5);
        %pos(wh(whr),1) = 1-pos(wh(whr),1);
        %v(wh(whl),1) = -v(wh(whl),1);
        v(wh,1) = -v(wh,1);
        ii = wh;
    end
    
    if numel(t_wall_y == min_t) > 0
        %adjust speeds for wall collision
        wh = find(t_wall_y == min_t);
        v(wh,2) = -v(wh,2);
        jj = wh;
    end
    
    if numel(t_col == min_t) > 0
        %adjust speeds for interparticle collision
        wh = find(t_col == min_t);
        kk = i_arr(wh);
        ll = j_arr(wh);
        for k = 1:numel(wh)
            i_mtc = kk(k);
            j_mtc = ll(k);
            dij = pos(i_mtc,:)-pos(j_mtc,:);
            dij = dij./sum(dij.^2)^0.5;
            dv = sum((v(i_mtc,:)-v(j_mtc,:)).*dij).*dij;
            v(i_mtc,:) = v(i_mtc,:)-dv;
            v(j_mtc,:) = v(j_mtc,:)+dv;
        end
    end
    
    %add momentum and energy to array after each collision to check
    %consistency
    momentum = [momentum,sum(v(:))];
    energy = [energy,sum(v(:,1).^2+v(:,2).^2)];
    
    %keyboard
end

figure
x = linspace(0,numel(momentum)-1,numel(momentum));
plot(x, momentum)
hold on
plot(x, energy)

figure
movie(F,1,20)


