function vecFunc = vectorizefunction(func, nobj)

    vecFunc = @newFunc;
    function yArray = newFunc(xArray)
        row = size(xArray,1);
        yArray = zeros(row, nobj);
        for i=1:row
            yArray(i,:) = func(xArray(i,:));
        end
    end
end