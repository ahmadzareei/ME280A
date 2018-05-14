function out= my_multiply(Ke_table,in,e)

out = zeros(size(e,1)+1,1);

out_m = zeros(2,1);
in_m = zeros(2,1);
Ke_m = zeros(1,3);

for element = 1:size(e,1)
    Ke_m = Ke_table(element,:);
    in_m = in(e(element,:));
    out_m(1) = Ke_m(1)*in_m(1) + Ke_m(2) * in_m(2);
    out_m(2) = Ke_m(2)*in_m(1) + Ke_m(3) * in_m(2);
    out(e(element,:)) = out(e(element,:)) + out_m;
end

