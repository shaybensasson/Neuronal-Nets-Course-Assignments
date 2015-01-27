function main()

    function res = fun2(getMeACallback)
            res = getMeACallback(10)
    end

    callback = @(x) x+5;


    parfor i=1:10
        fun2(callback)
    end
           

end