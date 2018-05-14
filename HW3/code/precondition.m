function Ke_table_pre = precondition(Ke_table,T,e);

Ke_table_pre = zeros(size(Ke_table));

for i=1:size(e,1)
    Ke_table_pre(i,:) = [Ke_table(i,1)*T(i).^2,...
                         Ke_table(i,2)*T(i)*T(i+1),...
                         Ke_table(i,3)*T(i+1)^2];
end

