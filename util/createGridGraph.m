function A = createGridGraph(width, height)

A = zeros(width * height, width * height);

c2i = @(x,y) height*(x-1)+y;
for i = 1:width
    for j = 1:height
        if j < height
            A(c2i(i,j),c2i(i,j+1)) = 1;
        end
        if j > 1
            A(c2i(i,j),c2i(i,j-1)) = 1;
        end
        if i < width
            A(c2i(i,j),c2i(i+1,j)) = 1;
        end
        if i > 1
            A(c2i(i,j),c2i(i-1,j)) = 1;
        end
    end
end

