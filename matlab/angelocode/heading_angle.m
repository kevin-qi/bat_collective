function theta = heading_angle(vx,vy)
A = vx<0;
theta = atan(vy./vx) + A.*pi;
theta = wrapTo2Pi(theta);
end

