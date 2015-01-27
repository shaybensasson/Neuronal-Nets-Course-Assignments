function I = feedExternalCurrent(p, percentile)
    I=0;    
    switch p.question
        case 1
           if(percentile<=1)
                I=35;
           end    
           
       case 2
           %if(percentile<=1)
                I=3;
           %end   
           
       case 3
           I=35;
       
        case 4
            if (percentile>=5 && percentile<=30) 
                I=20;
            elseif (percentile>=31 && percentile<=60)
                I=35;
            elseif (percentile>=61)
                I=55;
            end
        case 5
            %a: 
            %Mu=7; Sigma=100;
            %b: 
            %Mu=7; Sigma=10;
            %c
            %Mu=20; Sigma=100;
            %d
            %Mu=35; Sigma=150;
            %e
            %Mu=50; Sigma=10;
            %f
            %Mu=50; Sigma=100;
            global lastNormRndCurrent;
            global nextNormRndCalcPercentile;
    
            %I = normrnd(p.Mu,p.Sigma);
            
            if (percentile > nextNormRndCalcPercentile)
                lastNormRndCurrent = normrnd(p.Mu,p.Sigma);
                nextNormRndCalcPercentile = percentile + 1;

                
                %if (mod(percentile,10) == 0)
                    fprintf('normrnd calc at perc %d.\n', percentile);
                %end
            end
            I = lastNormRndCurrent;
        case 6
            %I=10;
            %I=30;
            %I=50;
            %I=100;
            I=250;
        case 7
            error('run multi compartment sim')
        case 8
            if(percentile>=3 && percentile<=30)
                I=-30;
            end
            
        otherwise
           I=p.NoQuestionConstantCurrent; 
           
   end
    
           
end