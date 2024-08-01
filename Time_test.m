clc;
clear;
tic;
m=200;
n=100;
r=1;
r0=1;
r1=1;
h=r/r0;
h1=r1/r0;
x1=0;
y1=0;

        for x2=0:1:n
            for y2=0:1:m
              
                     Z = recy9(x1,y1,x2,y2,m,n,r,h1,h);


            end
      
        end
        t=toc
       


function res=recy1(x1,y1,x2,y2,m,n,r,r0,r1,h1,h)
sum2=0;
for i=2:m
    sita=(i-1)*pi/(m);
    lm=1+h-h*cos(sita)+sqrt((1+h-h*cos(sita))^2-1);
    lmba=1/lm;
    A=lm-lmba;
    t=2+2*h-2*h*cos(sita);
    F1=(lm^x1-lmba^x1)/A;
    F2=(lm^(n-x1)-lmba^(n-x1))/A;
    Fn=(lm^n-lmba^n)/A;
    F0=0;
    g_11=(t-2)*(1-h1)*F1*F2+h1*(Fn+F0);

    F3=(lm^(n-x2)-lmba^(n-x2))/A;
    F4=(lm^(n-abs(x1-x2))-lmba^(n-abs(x1-x2)))/A;
    F5=(lm^(abs(x2-x1))-lmba^(abs(x2-x1)))/A;
    g_12=(t-2)*(1-h1)*F1*F3+h1*(F4+F5);

    F6=(lm^x2-lmba^x2)/A;
    g_22=(t-2)*(1-h1)*F6*F3+h1*(Fn+F0);

    s1=sin(y1*sita);
    s2=sin(y2*sita);

    F7=(lm^(n+1)-lmba^(n+1))/A;
    F8=(lm^(n-1)-lmba^(n-1))/A;

    DF_n=F7-Fn;
    DF_n_1=Fn-F8;

    sum1=(g_11*s1^2-2*g_12*s1*s2+g_22*s2^2)/(DF_n+(2*h1-1)*DF_n_1-2*h1);
    sum2=sum2+sum1;
end
res=(r0*r1*(y2-y1)^2)/(m*((n-1)*r1+r0))+2*r*sum2/m;
end



function res=recy7(x1,y1,x2,y2,m,n,r,r0,r1,h1,h)
sum2=0;
for i=2:m
    sita=(i-1)*pi/(m);
    t=2+2*h-2*h*cos(sita);
    chun=acosh(t/2);
    A=sinh(chun);
    U1=sinh(x1*chun)/A;
    U2=sinh((n-x1)*chun)/A;
    U3=sinh(n*chun)/A;
    U4=0;
    l_11=(t-2)*(1-h1)*U1*U2+h1*(U3+U4);

    U5 = sinh((n - x2) * chun) / A;
    U6 = sinh((n -abs( x1 - x2)) * chun) / A;
    U7 = sinh((abs( x1 - x2)) * chun) / A;
    l_12 = (t - 2) * (1 - h1) * U1 * U5 + h1 * (U6 + U7);

    U8=sinh(x2*chun)/A;
    l_22=(t-2)*(1-h1)*U8*U5+h1*(U3+U4);

    s1=sin(y1*(2*i-1)*pi/(2*m+1));
    s2=sin(y2*(2*i-1)*pi/(2*m+1));

    U9=sinh((n+1)*chun)/A;
    U10=sinh((n-1)*chun)/A;

    DU_n_1=U9-U3;
    DU_n_2=U3-U10;

    sum1=(l_11*s1^2-2*l_12*s1*s2+l_22*s2^2)/(DU_n_1+(2*h1-1)*DU_n_2-2*h1);
    sum2=sum2+sum1;
end
res=(r0*r1*(y2-y1)^2)/(m*((n-1)*r1+r0))+2*r*sum2/m;
end

function res=recy2(x1,y1,x2,y2,m,n,r,h1,h)
sum2=0;
for i=1:m
    sita=(2*i-1)*pi/(2*m+1);
    lm=1+h-h*cos(sita)+sqrt((1+h-h*cos(sita))^2-1);
    lmba=1/lm;
    A=lm-lmba;
    t=2+2*h-2*h*cos(sita);
    F1=(lm^x1-lmba^x1)/A;
    F2=(lm^(n-x1)-lmba^(n-x1))/A;
    Fn=(lm^n-lmba^n)/A;
    F0=0;
    g_11=(t-2)*(1-h1)*F1*F2+h1*(Fn+F0);

    F3=(lm^(n-x2)-lmba^(n-x2))/A;
    F4=(lm^(n-abs(x1-x2))-lmba^(n-abs(x1-x2)))/A;
    F5=(lm^(abs(x2-x1))-lmba^(abs(x2-x1)))/A;
    g_12=(t-2)*(1-h1)*F1*F3+h1*(F4+F5);

    F6=(lm^x2-lmba^x2)/A;
    g_22=(t-2)*(1-h1)*F6*F3+h1*(Fn+F0);

    s1=sin(y1*sita);
    s2=sin(y2*sita);

    F7=(lm^(n+1)-lmba^(n+1))/A;
    F8=(lm^(n-1)-lmba^(n-1))/A;

    DF_n=F7-Fn;
    DF_n_1=Fn-F8;

    sum1=(g_11*s1^2-2*g_12*s1*s2+g_22*s2^2)/(DF_n+(2*h1-1)*DF_n_1-2*h1);
    sum2=sum2+sum1;
end
res=4*r*sum2/(2*m+1);
end
function res=recy8(x1,y1,x2,y2,m,n,r,h1,h)
sum2=0;
for i=1:m
    sita=(2*i-1)*pi/(2*m+1);
    t=2+2*h-2*h*cos(sita);
    chun=acosh(t/2);
    A=sinh(chun);
    U1=sinh(x1*chun)/A;
    U2=sinh((n-x1)*chun)/A;
    U3=sinh(n*chun)/A;
    U4=0;
    l_11=(t-2)*(1-h1)*U1*U2+h1*(U3+U4);

    U5 = sinh((n - x2) * chun) / A;
    U6 = sinh((n -abs( x1 - x2)) * chun) / A;
    U7 = sinh((abs( x1 - x2)) * chun) / A;
    l_12 = (t - 2) * (1 - h1) * U1 * U5 + h1 * (U6 + U7);

    U8=sinh(x2*chun)/A;
    l_22=(t-2)*(1-h1)*U8*U5+h1*(U3+U4);

    s1=sin(y1*(2*i-1)*pi/(2*m+1));
    s2=sin(y2*(2*i-1)*pi/(2*m+1));

    U9=sinh((n+1)*chun)/A;
    U10=sinh((n-1)*chun)/A;

    DU_n_1=U9-U3;
    DU_n_2=U3-U10;

    sum1=(l_11*s1^2-2*l_12*s1*s2+l_22*s2^2)/(DU_n_1+(2*h1-1)*DU_n_2-2*h1);
    sum2=sum2+sum1;
end
res=4*r*sum2/(2*m+1);
end


function res=recy3(x1,y1,x2,y2,m,n,r,h1,h)
sum2=0;
for i=1:m
    sita=((2*i-1)*pi)/(2*m);
    lm=1+h-h*cos(sita)+sqrt((1+h-h*cos(sita))^2-1);
    lmba=1/lm;
    A=lm-lmba;
    t=2+2*h-2*h*cos(sita);
    F1=(lm^x1-lmba^x1)/A;
    F2=(lm^(n-x1)-lmba^(n-x1))/A;
    Fn=(lm^n-lmba^n)/A;
    F0=0;
    g_11=(t-2)*(1-h1)*F1*F2+h1*(Fn+F0);

    F3=(lm^(n-x2)-lmba^(n-x2))/A;
    F4=(lm^(n-abs(x1-x2))-lmba^(n-abs(x1-x2)))/A;
    F5=(lm^(abs(x2-x1))-lmba^(abs(x2-x1)))/A;
    g_12=(t-2)*(1-h1)*F1*F3+h1*(F4+F5);

    F6=(lm^x2-lmba^x2)/A;
    g_22=(t-2)*(1-h1)*F6*F3+h1*(Fn+F0);

    s1=sin(y1*sita);
    s2=sin(y2*sita);

    F7=(lm^(n+1)-lmba^(n+1))/A;
    F8=(lm^(n-1)-lmba^(n-1))/A;

    DF_n=F7-Fn;
    DF_n_1=Fn-F8;

    sum1=(g_11*s1^2-2*g_12*s1*s2+g_22*s2^2)/(DF_n+(2*h1-1)*DF_n_1-2*h1);
    sum2=sum2+sum1;
end
res=2*r*sum2/m;
end


function res=recy9(x1,y1,x2,y2,m,n,r,h1,h)
sum2=0;
for i=1:m
    sita=((2*i-1)*pi)/(2*m);
    t=2+2*h-2*h*cos(sita);
    chun=acosh(t/2);
    A=sinh(chun);
    U1=sinh(x1*chun)/A;
    U2=sinh((n-x1)*chun)/A;
    U3=sinh(n*chun)/A;
    U4=0;
    l_11=(t-2)*(1-h1)*U1*U2+h1*(U3+U4);

    % U5=sinh((n-x2)*chun)/A;
    % U6=sinh((n+x1-x2)*chun)/A;
    % U7=sinh((x2-x1)*chun)/A;
    % l_12=(t-2)*(1-h1)*U1*U5+h1*(U6+U7);
    U5 = sinh((n - x2) * chun) / A;
    U6 = sinh((n -abs( x1 - x2)) * chun) / A;
    U7 = sinh((abs( x1 - x2)) * chun) / A;
    l_12 = (t - 2) * (1 - h1) * U1 * U5 + h1 * (U6 + U7);

    U8=sinh(x2*chun)/A;
    l_22=(t-2)*(1-h1)*U8*U5+h1*(U3+U4);

    s1=sin(y1*sita);
    s2=sin(y2*sita);

    U9=sinh((n+1)*chun)/A;
    U10=sinh((n-1)*chun)/A;

    DU_n_1=U9-U3;
    DU_n_2=U3-U10;

    sum1=(l_11*s1^2-2*l_12*s1*s2+l_22*s2^2)/(DU_n_1+(2*h1-1)*DU_n_2-2*h1);
    sum2=sum2+sum1;
end
res=2*r*sum2/m;
end