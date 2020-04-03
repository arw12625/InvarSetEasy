function [poly_prod] = polyhedronProduct(poly_list)

Acell = cell(length(poly_list),1);
bcell = cell(length(poly_list),1);
Aecell = cell(length(poly_list),1);
becell = cell(length(poly_list),1);
for i = 1:length(poly_list)
    Acell{i} = poly_list{i}.A;
    bcell{i} = poly_list{i}.b;
    Aecell{i} = poly_list{i}.Ae;
    becell{i} = poly_list{i}.be;
end
A = blkdiag(Acell{:});
b = cell2mat(bcell);
Ae = blkdiag(Aecell{:});
be = cell2mat(becell);
poly_prod = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', be, 'irredundantHRep', true);

end

