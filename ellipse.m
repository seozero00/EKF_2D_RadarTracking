function ellipse(x,y,a,b)
% x and y are the coordinates of the center of the circle
% r is the radius of the circle
% 0.1 is the angle step, bigger values will draw the circle faster but
% you might notice imperfections (not very smooth)
ang = 0:0.01:2*pi; 
xp = a*cos(ang);
yp = b*sin(ang);
plot(x+xp,y+yp);