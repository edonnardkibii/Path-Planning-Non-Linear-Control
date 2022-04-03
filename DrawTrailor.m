function DrawTrailor(centerLocation,theta)

%%
%Setup
x = centerLocation(1);
y = centerLocation(2);

x_shift = x-x;
y_shift = y-y;



R = [cos(theta), -sin(theta);
     sin(theta),  cos(theta)];


%%
%Draw Car
x_coord = [x_shift-1, x_shift+1, x_shift+1, x_shift-1, x_shift-1];
y_coord = [y_shift-0.6, y_shift-0.6, y_shift+0.6, y_shift+0.6, y_shift-0.6];

x_rcoord = zeros(length(x_coord));
y_rcoord = zeros(length(y_coord));
for n=1:length(x_coord)
    x_rcoord(n) = R(1,1)*x_coord(n)+R(1,2)*y_coord(n);
    y_rcoord(n) = R(2,1)*x_coord(n)+R(2,2)*y_coord(n);
end
line(x_rcoord+x,y_rcoord+y)


%%
%Draw weheel 1
x_coord = [x_shift-0.2,x_shift+0.2,x_shift+0.2,x_shift-0.2,x_shift-0.2];
y_coord = [y_shift-0.45, y_shift-0.45, y_shift-0.55, y_shift-0.55, y_shift-0.45];
x_rcoord = zeros(length(x_coord));
y_rcoord = zeros(length(y_coord));
for n=1:length(x_coord)
    x_rcoord(n) = R(1,1)*x_coord(n)+R(1,2)*y_coord(n);
    y_rcoord(n) = R(2,1)*x_coord(n)+R(2,2)*y_coord(n);
end
line(x_rcoord+x,y_rcoord+y)

%%
%Draw weheel 2
x_coord = [x_shift-0.2,x_shift+0.2,x_shift+0.2,x_shift-0.2,x_shift-0.2];
y_coord = [y_shift+0.45, y_shift+0.45, y_shift+0.55, y_shift+0.55, y_shift+0.45];
x_rcoord = zeros(length(x_coord));
y_rcoord = zeros(length(y_coord));
for n=1:length(x_coord)
    x_rcoord(n) = R(1,1)*x_coord(n)+R(1,2)*y_coord(n);
    y_rcoord(n) = R(2,1)*x_coord(n)+R(2,2)*y_coord(n);
end

line(x_rcoord+x,y_rcoord+y)

