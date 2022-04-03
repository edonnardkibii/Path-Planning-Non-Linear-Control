%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%

global d0;
global d1;
global T;
global x0;
global x1;
global offset;
global f;
global xshift;

d0 = 2;
d1 = 2;
T = 10;


%Input Vektor for initial position
x0 = zeros(8,1);
x0(1) = 0;          %x
x0(2) = 0;          %y
x0(3) = 0;  %theta0
x0(4) = 0;  %theta1
x0(5) = 0;  %phi

x0(6) = 1;          %xi1
x0(7) = 0;          %xi2
x0(8) = 0;          %xi3

x1 = zeros(8,1);
x1(1) = 20;          %x
x1(2) = 20;          %y
x1(3) = -0.3;  %theta0
x1(4) = -0.3;   %theta1
x1(5) = -0.3;  %phi

x1(6) = 0;          %xi1
x1(7) = 0;          %xi2
x1(8) = 0;          %xi3

offset = zeros(5,1);
offset(1) = 0;
offset(2) = 0;
offset(3) = 0;
offset(4) = 0;
offset(5) = 0;

xshift = 1;
f = figure();
Recalculate();
updatePlot(0);

fig = uifigure('Resize','off');
uilabel(fig,'Position',[70,380,50,20],'Text','start');
uilabel(fig,'Position',[20,350,50,20],'Text','X');
u0(1) = uieditfield(fig,'numeric','Position',[70,350,50,20],'Value',x0(1));
uilabel(fig,'Position',[20,320,50,20],'Text','Y');
u0(2) = uieditfield(fig,'numeric','Position',[70,320,50,20],'Value',x0(2));
uilabel(fig,'Position',[20,290,50,20],'Text','Theta1');
u0(3) = uieditfield(fig,'numeric','Position',[70,290,50,20],'Value',x0(3));
uilabel(fig,'Position',[20,260,50,20],'Text','Theta2');
u0(4) = uieditfield(fig,'numeric','Position',[70,260,50,20],'Value',x0(4));
uilabel(fig,'Position',[20,230,50,20],'Text','Phi1');
u0(5) = uieditfield(fig,'numeric','Position',[70,230,50,20],'Value',x0(5));

uilabel(fig,'Position',[250,380,50,20],'Text','offset');
uilabel(fig,'Position',[200,350,50,20],'Text','X');
uof(1) = uieditfield(fig,'numeric','Position',[250,350,50,20],'Value',offset(1));
uilabel(fig,'Position',[200,320,50,20],'Text','Y');
uof(2) = uieditfield(fig,'numeric','Position',[250,320,50,20],'Value',offset(2));
uilabel(fig,'Position',[200,290,50,20],'Text','Theta1');
uof(3) = uieditfield(fig,'numeric','Position',[250,290,50,20],'Value',offset(3));
uilabel(fig,'Position',[200,260,50,20],'Text','Theta2');
uof(4) = uieditfield(fig,'numeric','Position',[250,260,50,20],'Value',offset(4));
uilabel(fig,'Position',[200,230,50,20],'Text','Phi1');
uof(5) = uieditfield(fig,'numeric','Position',[250,230,50,20],'Value',offset(5));

uilabel(fig,'Position',[450,380,50,20],'Text','end');
uilabel(fig,'Position',[400,350,50,20],'Text','X');
u1(1) = uieditfield(fig,'numeric','Position',[450,350,50,20],'Value',x1(1));
uilabel(fig,'Position',[400,320,50,20],'Text','Y');
u1(2) = uieditfield(fig,'numeric','Position',[450,320,50,20],'Value',x1(2));
uilabel(fig,'Position',[400,290,50,20],'Text','Theta1');
u1(3) = uieditfield(fig,'numeric','Position',[450,290,50,20],'Value',x1(3));
uilabel(fig,'Position',[400,260,50,20],'Text','Theta2');
u1(4) = uieditfield(fig,'numeric','Position',[450,260,50,20],'Value',x1(4));
uilabel(fig,'Position',[400,230,50,20],'Text','Phi1');
u1(5) = uieditfield(fig,'numeric','Position',[450,230,50,20],'Value',x1(5));

btn = uibutton(fig,...
                'Text','calculate',...
                'Position',[230 100 100 20],...
                'ButtonPushedFcn', @(btn,event) ButtonPushed(u0,u1,uof));

sld = uislider(fig,...
               'Position',[75 55 400 3],...
               'Limits',[0 T],...
               'ValueChangingFcn',@(sld,event) sliderMoving(event));


function ButtonPushed(u0,u1,uof)
    global x0
    global offset
    global x1
    global xshift
    xshift = u0(1).Value;
    x0(1) = 0;
    x0(2) = u0(2).Value;
    x0(3) = u0(3).Value;
    x0(4) = u0(4).Value;
    x0(5) = u0(5).Value;
    
    offset(1) = uof(1).Value;
    offset(2) = uof(2).Value;
    offset(3) = uof(3).Value;
    offset(4) = uof(4).Value;
    offset(5) = uof(5).Value;
    
    x1(1) = u1(1).Value - xshift;
    x1(2) = u1(2).Value;
    x1(3) = u1(3).Value;
    x1(4) = u1(4).Value;
    x1(5) = u1(5).Value;
    
    Recalculate();
    updatePlot(0);
end

function sliderMoving(event)
    i = event.Value();
    updatePlot(i)
end

function updatePlot(i)
    global d1
    global ist
    global state
    global soll
    global f
    global T
    global xshift
    clf(f)
    index = uint32(length(state(:,1))*i/T);
    plot(ist(1,:)+xshift,ist(2,:))
    hold on
    plot(soll(1,:)+xshift,soll(2,:))
    if index==0 
        index=1; 
    end
    DrawTruck([state(index,1)+xshift;state(index,2)],state(index,3),state(index,5))
    DrawTrailor([state(index,1)-d1*cos(state(index,4))+xshift;state(index,2)-d1*sin(state(index,4))],state(index,4))
    axis equal
end

function Recalculate()
    global d0
    global d1
    global T
    global x1
    global x0
    global offset
    global ist 
    global soll
    global t
    global state
    
    opts = odeset('RelTol',1e-12);

    coef = PathPlanner(x0,x1,d0,d1);
    initialState = [x0(1)+d1*cos(x0(4))+offset(1);x0(2)+d1*sin(x0(4))+offset(2);x0(3)+offset(3);x0(4)+offset(4);x0(5)+offset(5);x0(6:8)];
    [t,state] = ode45(@myfun,[0,T],initialState,opts,coef,x0(1),x1(1));
    ist = [state(:,1)-d1*cos(state(:,4)) , state(:,2)-d1*sin(state(:,4))]';
    soll = [[x0(1):0.1:x1(1)] ; polyval(coef,[x0(1):0.1:x1(1)])];
end


function dx_dt=myfun(t,state,coef,x_start,x_end)

    global d1;
    global d0;
    global T;
    
    tau = t/T;
    s = 3*tau^2-2*tau^3;
    
    x        = x_start+(x_end-x_start)*s;
    x_dot    = 1/T*(x_end-x_start)*(6*tau-6*tau^2);


    y_ref = polyval(coef,x);
    yx = polyval(polyder(coef),x);
    yxx = polyval(polyder(polyder(coef)),x);
    yxxx = polyval(polyder(polyder(polyder(coef))),x);
    yxxxx = polyval(polyder(polyder(polyder(polyder(coef)))),x);
    
    LD = Truck_1T_LieDeriv(state(1:5), state(6:8), d0, d1);
    
    ysigma = LD.Lf_h2;
    ysigma2 = LD.Lf2_h2;
    ysigma3 = LD.Lf3_h2;
    
    xsigma  = LD.Lf_h1;
    xsigma2 = LD.Lf2_h1;
    xsigma3 = LD.Lf3_h1;

    [x_ref_prime_4,y_ref_prime_4] =  Fourth_Derivatives_of_References(yx,yxx,yxxx,yxxxx);
    
    eta = sqrt(1+(yx)^2);
        
    x_ref_prime = 1/eta;
    y_ref_prime = yx/eta;
    
    x_ref_prime_2 = -yx*yxx/eta^4;
    y_ref_prime_2 = yxx/eta^2 - yx^2*yxx/eta^4;
    
    x_ref_prime_3 = -(yxx^2+yx*yxxx)/eta^5+4*yx^2*yxx^2/eta^7;
    y_ref_prime_3 = yxxx/eta^3-4*yxx^2*yx/eta^5-yx^2*yxxx/eta^5+4*yx^3*yxx^2/eta^7;
    
    e_x = x - (state(1) - d1 * cos(state(4)));
    e_y = y_ref - (state(2) - d1 * sin(state(4)));

    coef_k = poly(-3/d0*ones(1,4));
    k3 = coef_k(2);
    k2 = coef_k(3);
    k1 = coef_k(4);
    k0 = coef_k(5);
    
 
    nu_1 = x_ref_prime_4  +  k3 * (x_ref_prime_3-xsigma3) + k2 * (x_ref_prime_2-xsigma2) + k1 * (x_ref_prime-xsigma) + k0 *e_x; 
    nu_2 = y_ref_prime_4  +  k3 * (y_ref_prime_3-ysigma3) + k2 * (y_ref_prime_2-ysigma2) + k1 * (y_ref_prime-ysigma) + k0 *e_y;
    
    w = [LD.L_g1_Lf3_h1 LD.L_g2_Lf3_h1; LD.L_g1_Lf3_h2, LD.L_g2_Lf3_h2]^-1*([nu_1;nu_2]-[LD.Lf4_h1;LD.Lf4_h2]);
    w_1 = w(1);
    w_2 = w(2);

    dx_dt = zeros(8,1);
    v0 = state(6)/(cos(state(3)-state(4))) * eta * x_dot;
    dx_dt(1) = v0*cos(state(3));
    dx_dt(2) = v0*sin(state(3));
    dx_dt(3) = v0/d0*tan(state(5));
    dx_dt(4) = v0/d1*sin(state(3)-state(4));
    dx_dt(5) = w_2 * eta *x_dot;
    dx_dt(6) = state(7) * eta * x_dot;
    dx_dt(7) = state(8) * eta * x_dot;
    dx_dt(8) = w_1 * eta * x_dot;
end


