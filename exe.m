disp('TP 10');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.

function [ x, y ] = NURBS( interval, m, w, xp, yp, t )
  n = length(xp);
  x = zeros(length(interval), 1);
  y = zeros(length(interval), 1);
  for j=1:length(interval)
    denum = 0;
    for i=1:n
      denum = denum + w(i) * fonctionN(t, interval(j), m, i);
    end
    if denum == 0
      denum = 1;
    end
    for i=1:n
      x(j) = x(j) + xp(i) * w(i) * fonctionN(t, interval(j), m, i);
      y(j) = y(j) + yp(i) * w(i) * fonctionN(t, interval(j), m, i);
    end
    x(j) = x(j) / denum;
    y(j) = y(j) / denum;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.

m = 2;
k = 3;

xp = [ 0 -1 sqrt(3)/2 0 sqrt(3)/2 1 0];
yp = [ 0 0 -1/2 1 1/2 0 0];
t = [0, 0, 0,1/3,1/3,2/3,2/3, 1, 1, 1];
w = [1,1/2, 1,1/2, 1,1/2, 1];

n = length(xp);
interval = -1:0.01:1;

[x, y] = NURBS( interval, m, w, xp, yp, t);

figure;
plot(x, y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.

m = 2;
k = 3;

xp = [ 0 -1 -1 -1 0 1 1 1 0];
yp = [ 0 0 1/2 1 1 1 1/2 0 0];
t = [0, 0, 0,1/4,1/4,1/2,1/2, 3/4, 3/4, 1, 1, 1];
w = [1,sqrt(2)/2, 1,sqrt(2)/2, 1,sqrt(2)/2, 1,sqrt(2)/2, 1];

interval = 0:0.01:sqrt(3)/2+0.5;

[x, y] = NURBS( interval, m, w, xp, yp, t);

% figure;
% plot(x, y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4.

m = 2;
k = 3;

xp = [ 0 -1 sqrt(3)/2 0 sqrt(3)/2 1 0];
yp = [ 0 0 -1/2 1 1/2 0 0 ];
t = [0,0,0,1/3,1/3,2/3,2/3,1,1,1];
w = [1,1/2, 1,1/2, 1,1/2, 1];

interval = 0:0.01:sqrt(3)/2+0.5;

[x, y] = NURBS( interval, m, w, xp, yp, t);

% figure;
% plot(x, y);

pause;
