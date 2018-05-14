function [p,e] = mesh (ri, ro, NEr, NEth)
%this function generates the mesh for the FEM
% with NEr nodes  in r direction and 
% NEth nodes in the theta direction
% p will be list of points in x and y direction 
% e will be node numbers in elemnts

dr = (ro-ri)/NEr;
dth = 2*pi/NEth;
[r,th] = meshgrid(ri:dr:ro,0:dth:2*pi);
x = r.*cos(th);
y = r.*sin(th);
p = zeros((NEr+1)*NEth,2);
e = zeros(NEr*NEth,4);
for j=1:NEth
    for i=1:NEr+1
        p((j-1)*(NEr+1)+ i,:) = [x(j,i),y(j,i)];    
    end
end

l = 1;
cof = 0;

for j=1:NEth
    for i=1:NEr
        cof = (j)*(NEr+1)+ i;
        if( j == NEth)
            cof = i;
        end
        e(l,:) = [(j-1)*(NEr+1)+ i,(j-1)*(NEr+1)+ i+1,cof+1, cof];    
        l=l+1;
    end
end
