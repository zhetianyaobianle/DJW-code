clc;
clear;

n = 10; 
f = 10; 
r_t = 1; 
r_u = 1; 
x1=7;
y1 = 8;
x2=3;
y2=2;


[x,y]= meshgrid(1:1:f, 1: 1:n);

zong2 = 0.0 ;
        for j = 2:n 
            hou2 = Computefunction(n,f,x,y,x1,y1,x2,y2,r_t,r_u,j);
            zong2 =zong2 + hou2;
        end
        U=(y.*(y1-y2)*r_u/(n*f))+zong2;
        U

p = surf(x,y,U);

start=(y1-1)*n+x1;
target=(y2-1)*n+x2;


x=zeros(n*f,1);
for i=0:n-1
    for j=1:f
        x(i*n+j)=U(i+1,j);
    end
end
% x

obs=12;
obs1_x=2;obs1_y=1;
obs2_x=2;obs2_y=2;
obs3_x=2;obs3_y=4;
obs4_x=2;obs4_y=5;
obs5_x=1;obs5_y=6;
%obs6_x=5;obs6_y=7;
obs7_x=5;obs7_y=6;
obs8_x=4;obs8_y=3;
obs9_x=4;obs9_y=4;
obs10_x=7;obs10_y=4;
obs11_x=6;obs11_y=9;
obs12_x=5;obs12_y=9;


OBS=zeros(1, 12);
%for循环取值
for i = 1:obs
    if i==6
        continue
    end
    % 构造变量名
    var_name_x = ['obs', num2str(i), '_x'];
    var_name_y = ['obs', num2str(i), '_y'];
    %obs = ['obs', num2str(i), '_y'];

    % 使用 eval() 函数获取变量的值
    value_x = eval(var_name_x);
    value_y = eval(var_name_y);

    % 输出变量的值

    OBS(i)=(value_y-1)*n+value_x;
    x(OBS(i))=x(OBS(i))+0.2;
end



V=zeros(n,f);
for i = 0:n-1
    for j = 1:f
        fprintf('%d    ',x(i*4+j));
        V(i+1,j)=x(i*n+j);
    end
    disp(' ');
end


% V
surf(1:n,1:f,V);
xlabel(' {\it x}- axis');  ylabel(' {\it{y}}- axis');  zlabel(' U(x,y)/J ');
set(gca,'fontweight','bold','FontSize',27);
size_data = 70;
hold on;
scatter3(x1, y1, V(y1,x1),size_data, 'filled', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r'); % 添加一个红色点
hold on;
scatter3(x2, y2, V(y2,x2),size_data, 'filled', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g'); % 添加一个红色点




hold off;



 h=0.2;






while start~=target
    legend('', 'Starting point','Ending point', 'Location', 'northeast'); 
    nn=Nnext(x,start,f,n)
    %x(nn);

    d1=[fix((start-1)/n)+1,fix((nn-1)/n)+1];
    d2=[start-n*(fix((start-1)/n)),nn-n*(fix((nn-1)/n))];
    %Wri(d2,d1);
    d11=[d2(2),d1(2),x(nn)];
    d22=[d2(1),d1(1),x(start)];
    Wrii(d11,d22);

    start=nn;
end

% legend on;
legend('', 'Starting point','Ending point', 'Location', 'northeast'); 


function next=Nnext(x,start,n,m)
    %global m;global n;
    if start<=n
        if start==1
            minIndex = Min2(x,Min2(x,start+1,start+n),start+n+1);
        elseif start==n
            minIndex = Min2(x,Min2(x,start-1,start+n),start+n-1);
        else
            minIndex = Min3(x,Min3(x,start-1,start+1,start+n),start+n-1,start+n+1);
        end
    elseif start+n>m*n
        if mod((start-1),n)==0
            minIndex = Min2(x,Min2(x,start-n,start+1),start-n+1);
        elseif start==m*n
            minIndex = Min2(x,Min2(x,start-1,start-n),start-n-1);
        else
            minIndex = Min3(x,Min3(x,start-1,start-n,start+1),start-n-1,start-n+1);
        end
    else
       if mod((start-1),n)==0
           minIndex = Min3(x,Min3(x,start-n,start+1,start+n),start-n+1,start+n+1);
       elseif mod((start),n)==0
          minIndex = Min3(x,Min3(x,start-n,start+1,start+n),start-n-1,start+n-1);
       else
          minIndex = Min2(x,Min4(x,start-n,start-1,start+1,start+n),Min4(x,start-n-1,start+n-1,start-n+1,start+n+1));
       end
    end
    next=minIndex;

end

function wrii=Wrii(x1,y1)

    hold on;
    plot3([x1(1),y1(1)],[x1(2),y1(2)], [x1(3),y1(3)],'r-','LineWidth', 4); % 使用红色虚线连接两点
    hold off;
    wrii=0;

end
function wri=Wri(x1,y1)

    plot(x1,y1, 'r--','LineWidth', 2); % 使用红色虚线连接两点
    wri=0;
end
function mmin=Min2(x,a1,a2)
    if x(a1,1)<=x(a2,1)
        mmin=a1;
    else
        mmin=a2;
    end
end
function mmin=Min3(x,a1,a2,a3)

    if x(a1,1)<=x(a2,1)
        mmin=a1;
    else
        mmin=a2;
    end
    if x(mmin,1)>x(a3,1)
        mmin=a3;
    end
end
function mmin=Min4(x,a1,a2,a3,a4)
    if x(a1,1)<=x(a2,1)
        mmin=a1;
    else
        mmin=a2;
    end
    if x(mmin,1)>x(a3,1)
        mmin=a3;
    end
    if x(mmin,1)>x(a4,1)
        mmin=a4;
    end
end










