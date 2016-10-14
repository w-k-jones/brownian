dbstop if error

balls=40;

pos = rand(balls,2).*(1-dist)+dist/2;
vs = rand(balls,2)-0.5;
dist=0.02;

steps = 5;
tstp = 0;
lstp=0.01;
ttot = 0;
istp = 0;
nstp = steps/lstp;
F(nstp) = struct('cdata',[],'colormap',[]);
ii=1;
jj=1;

momentum = sum(vs(:))
energy = sum(vs(:,1).^2+vs(:,2).^2)

%continue until gt simulation time
while tstp<=steps
    
    %add momentum and energy to array after each collision to check
    %consistency
    momentum = [momentum,sum(vs(:))];
    energy = [energy,sum(vs(:,1).^2+vs(:,2).^2)];
    
    %create NaN arrays for wall/collision times at start of each loop to
    %overwrite all previous values.
    t_wall = NaN(balls,4);
    t_col = NaN(balls,balls);
    
    %loop over balls to compare to all other balls
    for i = 1:balls
        
        %calculate collision times with horizontal walls
        temp_wall_x = max([-(pos(i,1)-dist/2)/vs(i,1), -(pos(i,1)-(1-dist/2))/vs(i,1)]);
        temp_wall_y = max([-(pos(i,2)-(1-dist/2))/vs(i,2), -(pos(i,2)-dist/2)/vs(i,2)]);
        
        %write lt 0 values as NaN
        temp_wall_x(temp_wall_x<0) = NaN;
        temp_wall_y(temp_wall_y<0) = NaN;
        
        %all balls in bounds should have 2 t >0, 2t<0, if not get rid of
        %checks.
        t_wall(i,:) = min([temp_wall_x,temp_wall_y],[],'omitnan');
        
        if i==balls
            break
        end
        
        %loop over all higher index balls to check for collisions
        for j = i+1:balls
            %calculate quadratic coefficients
            a = sum((vs(i,:)-vs(j,:)).^2);
            b = 2*sum((vs(i,:)-vs(j,:)).*(pos(i,:)-pos(j,:)));
            c = sum((pos(i,:)-pos(j,:)).^2)-dist^2;
            %check if valid solutions
            chk = b^2-4*a*c;
            if chk>=0
                
                temp = [(-b+chk^0.5)/(2*a),(-b-chk^0.5)/(2*a)];
                %both solutions should be gt 0 if valid collision,
                %otherwise balls havepassed through each other!
                if numel(temp(temp>=0))==2
                    t_col(i,j) = min(temp);
                    time = min(temp);
                    %keyboard
                elseif numel(temp(temp>=0))==1
                    %keyboard
                end
                
            end
        end
    end
    
    t_wall(t_wall<0) = NaN;
    
    if (numel(t_col(isfinite(t_col)))>0) | (numel(t_wall(isfinite(t_wall)))>0)
        if numel(t_col(isfinite(t_col)))>0
            [min_t_col,ind] = min(t_col(:),[],'omitnan');
        else
            min_t_col = inf;
        end
        if numel(t_col(isfinite(t_col)))>0
            [min_t_wall,ind_wall] = min(t_wall(:),[],'omitnan');
        else
            min_t_wall = inf;
        end
        
        if min_t_col < min_t_wall
            while (ttot+min_t_col)>(tstp+lstp)
                if abs((tstp+lstp)-(ttot+min_t_col)) < lstp
                    pos = pos+(vs*abs((tstp+lstp)-(ttot+min_t_col)));
                else
                    pos = pos+vs*lstp;
                end
                figure(1);
                hold off
                scatter(pos(:,1),pos(:,2));
                hold on
                scatter(pos(ii,1),pos(ii,2),'filled');
                scatter(pos(jj,1),pos(jj,2),'filled');
                axis([0,1,0,1]);
                drawnow
                tstp = tstp+lstp;
                istp = istp+1;
                F(istp) = getframe;
                if numel(pos(pos<0))>0 | numel(pos(pos>1))>0
                    keyboard
                end
            end
            %keyboard
            ttot = ttot+min_t_col;
            [ii,jj] = ind2sub(size(t_col),ind);
            pos = pos+(vs*(ttot-tstp));
            
            dij = pos(jj,:)-pos(ii,:);
            
            sum(dij.^2);
            
            dij = dij./sum(dij.^2)^0.5;
            dv = sum((vs(ii,:)-vs(jj,:)).*dij).*dij;
            %keyboard
            if sum(dv)==0
                dv = [0.1,0.1]
            end
            
            %if sum([(vs(ii,:)-dv),(vs(jj,:)+dv)].^2) ~= sum([vs(ii,:),vs(jj,:)].^2)
            %keyboard
            %end
            
            vs(ii,:) = vs(ii,:)-dv;
            vs(jj,:) = vs(jj,:)+dv;
            
        else
            while (ttot+min_t_wall)>(tstp+lstp)
                if abs((tstp+lstp)-(ttot+min_t_wall)) < lstp
                    pos = pos+(vs*abs((tstp+lstp)-(ttot+min_t_wall)));
                else
                    pos = pos+vs*lstp;
                end
                figure(1);
                hold off
                scatter(pos(:,1),pos(:,2));
                hold on
                scatter(pos(ii,1),pos(ii,2),'filled');
                axis([0,1,0,1]);
                drawnow
                tstp = tstp+lstp;
                istp = istp+1;
                F(istp) = getframe;
                if numel(pos(pos<0))>0 | numel(pos(pos>1))>0
                    keyboard
                end
            end
            %keyboard
            ttot = ttot+min_t_wall;
            [ii,jj] = ind2sub(size(t_wall),ind_wall);
            %keyboard
            pos = pos+(vs*(ttot-tstp));
            if jj==2 | jj==4
                vs(ii,2) = -vs(ii,2);
            else
                vs(ii,1) = -vs(ii,1);
            end
            %keyboard
        end
    else
        break
    end
    t_col_save = t_col;
    t_wall_save = t_wall;
    
    
end

figure
x = linspace(0,numel(momentum)-1,numel(momentum))
plot(x, momentum)
hold on
plot(x, energy)

figure
movie(F,1,20)








