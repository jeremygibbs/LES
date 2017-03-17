d = 4;
r = -2:2;
[rx,ry,rz] = meshgrid(r,r,r);
rv = sqrt(rx.^2 + ry.^2 + rz.^2);
gk = sqrt(6.0/pi/d/d)*exp(-6.0.*rv.^2/d/d);
(sum(sum(sum(rv))))/125
