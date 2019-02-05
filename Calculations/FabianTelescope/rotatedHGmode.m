function [ rotHGmode ] = rotatedHGmode( m,n,x,y,zx,zRx,zy,zRy,lambda,angle)
%Creates a rotated HG mode

HGmode = getHG(m,n,x,y,zx,zRx,zy,zRy,lambda);
rotated = imrotate(HGmode,angle);

if length(x(:,1))==length(rotated(:,1))&&length(y(1,:))==length(rotated(:,1))
    rotHGmode = rotated;
else
    x1crop = (length(rotated(1,:))-length(x))/2;
    x2crop = length(x)-1;
    y1crop = (length(rotated(:,1))-length(y))/2;
    y2crop = length(y)-1;
    rotHGmode = imcrop(rotated, [x1crop y1crop x2crop y2crop]);
end

end

