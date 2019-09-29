% Linear shooting
% Approximate the solution of y"=(-2/x)y'+(2/x^2)y+ sin(lnx)/x^2
% for 1<=x<=2 with y(1)=1 and y(2)=2.

%Consider the boundary value problems (BVPs) for the second order differential equation of the form
%          (*) y?? = f(x, y, y?), a ? x ? b, y(a) = alpha and y(b) = beta.

%Existence and Uniqueness of Solutions of BVPs:
%Theorem: Suppose that f is continuous on the set:
%           D = {(x, y, y?); a ? x ? b, ?inf < y < inf,?inf < y' < inf}
%           and the partial derivatives fy and fy? are continuous on D. 
%If fy(x, y, y?)>  0, for all (x, y, y?) in D, and there exists a constant M such that |fy?(x, y, y?)|? M for all (x, y, y?) in D

%For the linear boundary value problem, f(x, y, y?)=p(x)y?+q(x)y + r(x). Observe that
%fy(x, y, y?)=q(x), fy?(x, y, y?)=p(x). 
%fy(x, y, y?) and fy?(x, y, y?) are continuous on D if and only if p(x), q(x) and r(x) are continuous for a ? x ? b. Now we check conditions in (i) and (ii).
%(i) fy(x, y, y?)=q(x) > 0 for a ? x ? b
%(ii) Since fy? is continuous on [a, b] , fy? is bounded
%So, if p(x), q(x) and r(x) are continuous for a ? x ? b, and q(x)> 0 for a ? x ? b, then this linear boundary value problem has a unique solution.


%Consider the solutions of the following two initial-value problems:
%       (**) (i) y?? = p(x)y? + q(x)y + r(x), a ? x ? b, y(a)=alpha, y?(a)= 0
%           (ii) y?? = p(x)y? + q(x)y, a ? x ? b, y(a)=0, y?(a)= 1
%say, y1(x) and y2(x). Let y(x) be the following linear combination of y1(x) and y2(x) :
%           y(x)=y1(x)+(beta-y1(b))/y2(b)*y2(x)



% Linear shooting
% Approximate the solution of y"=(-4/x)y'+(-2/(x^2))y+ (2*log(x)/(x^2))
% for 1<=x<=2 with y(1)=1/2 and y(2)=log(2).
function shooting(a,b,alpha,beta,n)
    p = @(x) (-4/x);  
    q = @(x) (-2/(x^2));
    r = @(x) (2*log(x)/(x^2)); 

    a=1; b=2; alpha =1/2; beta = log(2);  n = 20;
    h = (b-a)/n; u1 = alpha;    u2 = 0;    v1 = 0;    v2 = 1;
    u = zeros(2,n); v = zeros(2,n);
    
    fprintf('Cheking cond - (1) \n');
    xes=zeros(n,1);
    for i = 1 : n 
        xes(i)=a+(i)*h;
    end
    count=0;
    for i = 1 : n 
        if q(xes(i))<=0
            count=count+1;
            break
        end            
    end
       
    if count>0
        fprintf('we cannot say if this boundary value problem has a unique solution \n');
    else
        fprintf('This boundary value problem has a unique solution \n');
    end
    
     for i = 1 : n 
       x = a+(i-1)*h;
       t = x+0.5*h;
       k11 = h*u2;
       k12 = h*(p(x)*u2+q(x)*u1+r(x));
       k21 = h*(u2+0.5*k12);
       k22 = h*(p(t)*(u2+0.5*k12)+q(t)*(u1+0.5*k11)+r(t));
       k31 = h*(u2+0.5*k22);
       k32 = h*(p(t)*(u2+0.5*k22)+q(t)*(u1+0.5*k21)+r(t));
       t = x+h;
       k41 = h*(u2+k32);
       k42 = h*(p(t)*(u2+k32)+q(t)*(u1+k31)+r(t));
       u1 = u1+(k11+2*(k21+k31)+k41)/6;
       u2 = u2+(k12+2*(k22+k32)+k42)/6;
       k11 = h*v2;
       k12 = h*(p(x)*v2+q(x)*v1);
       t = x+0.5*h;
       k21 = h*(v2+0.5*k12);
       k22 = h*(p(t)*(v2+0.5*k12)+q(t)*(v1+0.5*k11));
       k31 = h*(v2+0.5*k22);
       k32 = h*(p(t)*(v2+0.5*k22)+q(t)*(v1+0.5*k21));
       t = x+h;
       k41 = h*(v2+k32);
       k42 = h*(p(t)*(v2+k32)+q(t)*(v1+k31));
       v1 = v1+(k11+2*(k21+k31)+k41)/6;
       v2 = v2+(k12+2*(k22+k32)+k42)/6;
       u(1,i) = u1;
       u(2,i) = u2;
       v(1,i) = v1;
       v(2,i) = v2;
     end

    res=zeros(n,1);
    res(1)=alpha; xes(1)=a;
    coff = (beta-u(1,n))/v(1,n);
     for i = 2 : n 
       res(i) = u(1,i)+coff*v(1,i);
     end

    plot(xes(:,1),res(:,1),'b-',xes(:,1),u(1,:),'r-.',xes(:,1),v(1,:),'m-.');
    title('Shooting Method: y"=Px*y+Qx*y+Rx, y(1)=alpha, y(2)=beta')
    text(1.2,0.4,'IVP: y_1"=Px*y+Qx*y+Rx, y_1(1)=alpha, y_1(1)=0')
    text(1.3,0.1,'IVP: y_2"=Px*y+Qx*y, y_2(1)=0, y_2(1)=1')
    text(1.2,0.6,'- y=y_1+(beta-y_1(2))/y_2(2)*y_2')
    axis([a b 0 1.2*res(n,1)])
end
