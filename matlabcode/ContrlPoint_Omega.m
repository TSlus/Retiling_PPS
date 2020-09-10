%% 每个面片上的控制顶点
ctrPs{numP}=[];
doctrP(numP) = 0; nullctrP(numP) = 0;
eps = 1e-10;

% load('n_choose_k.mat');
for i=1:numP
    mu=v_valence(i);
    L=abs(cos(pi/mu));  
    %在正方形内部产生4*(mu+1)^2个点
    len = 2 * ( mu + 1 );
    t = linspace(0,1,len); t = (t-0.5) * 2 * L;
    A = repmat(t,len,1);
    square=[]; square(:,2)=A(:);  %y坐标
    A = A'; square(:,1) = A(:);       %x坐标
    
    %将Q选出来，Loop网格上的交点
    Q = zeros(size(square,1), 3); ti = 1;
    Quv=zeros(size(square,1), 2); %留下领域内点q的uv坐标
    for q =1:size(square,1)
        x = square(q,1); y = square(q,2);
        theta = atan2(y,x);
        if theta < 0
           theta = theta + 2 * pi;
        end
        fidx = fix(theta/(2*pi/mu))+1;
        %正多边形脚点
        pu=[]; pu(:,1) = cos( (0:mu-1 )*( 2 * pi / mu)); 
        pu(:,2) = sin( (0:mu-1) * (2 * pi / mu) );
        a1=pu(fidx,:); fidxplus = fidx + 1;
        if fidx == mu 
            fidxplus = 1; 
        end
        a2 = pu(fidxplus,:);
        signs = cross([a1,0],[a2,0]); s = norm(signs);
        xy = [x,y];
        signs1 = cross([a1,0],[xy,0]);
        signs2 = cross([xy-a2,0],[a1-a2,0]);
        signs3 = cross([xy,0],[a2,0]);
        c1 = dot(signs1,signs); c2 = dot(signs2,signs); c3 = dot(signs3,signs);
        
        %排除空心点
        if c1 >= -eps && c2 >= -eps && c3 >= -eps
            lambda = norm(signs2)/s; view = norm(signs3)/s;%yita=norm(signs1)/s;
            %不直接计算，是为了保证重心坐标计算的点，在三角形内部
            lambda(abs(lambda) > 1) = 1;
            yita=1-lambda-view;
            yita(yita < 0) = 0;
            view=1 - lambda - yita;
            %将多边形第一个点映射到v,邻域顶点的第一个点
            %通过fidx
            p1_idx = oneRingPs{i}(fidx); p2_idx = oneRingPs{i}(fidxplus);
            
            %计算b(p)
            b1=vertices(i,:);b2=vertices(p1_idx,:);b3=vertices(p2_idx,:);
            qmesh = lambda * b1 + view * b2 + yita * b3;
            
            %找到[b1,b2,b3]或者[b1,b3,b2]对应的faces位置fi
            %在fi上计算与qmesh最近的点
            fi_u = vf_sparse(p1_idx, p2_idx);
            [~, pointidx] = min(sum(abs(loop_point{fi_u} - qmesh).^2,2));
            bq = loop_point{fi_u}(pointidx, :);
            Q(ti,:) = bq;
            Quv(ti, :)=[x,y];
            ti = ti+1;
        end   
    end
    Q = Q(1:ti-1,:);
    Quv = Quv(1:ti-1,:);
    
    m = mu+1; n = mu+1;
    if size(Q,1) <(m+1) * (n+1)   %保证(M'*M)非奇异
        disp('第i和邻域Q的数量不够，i=')
        disp(i);
        continue;
    end
    
    %将Q中的点取出来，当成控制点，做Omega(u)上的Bezier曲面
    m=mu+1; n = m;
    cB = zeros(1, m+1);
    %nchoosek(m,0:m);
    for ic=0:m
        %cB(ic+1) = nchoosek(m,ic);
        cB(ic+1) = n_choose_k(m+1, ic+1);
    end
    M = zeros(size(Q,1),(m+1)*(m+1)); % 计算最小二乘系数矩阵M
    for im=1:size(Q,1)
        u = Quv(im,1:2); u=u'; u = u / (2 * L) + 0.5;%参数范围在(0,1)*(0,1)
        Bu = cB .*(u .^(0:m)).*((1-u).^(m-(0:m)));
        Buv = Bu(1,:)'.* Bu(2,:);
        M(im,:) = Buv(:)';
    end
    
    bctrl=(M'*M)\(M'*Q);
    ctrPs{i} = bctrl;
end
