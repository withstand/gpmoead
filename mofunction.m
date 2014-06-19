function ret = mofunction(funcs)

ret = @newFunction;


    function [yArray, stdArray] = newFunction(xArray)
        n = numel(funcs);
        yArray = zeros(size(xArray,1), n);
        stdArray = yArray;
        for iFuncs = 1:n
            yStar = funcs{iFuncs}(xArray);
            yArray(:, iFuncs) = yStar(:,1);
            stdArray(:, iFuncs) = yStar(:,2);
        end
        
        
    end

end