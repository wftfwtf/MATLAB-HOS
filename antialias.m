function [aliasfilter]=antialias(NROW,NCOL,mdim,strict)
aliasfilter=ones(NROW,NCOL);

if NROW==1
    [aliasfilter]=antialias1d(NROW,NCOL,mdim,strict,'x');
elseif NCOL==1
    [aliasfilter]=antialias1d(NROW,NCOL,mdim,strict,'y');
else
    for j=1:NCOL
        for i=1:NROW
            if NCOL>2 && abs(j-NCOL/2)> floor(NCOL/(mdim+1))+strict
                aliasfilter(i,j)=0;
            end
            
            if NROW>2 && abs(i-NROW/2)> floor(NROW/(mdim+1))+strict
                aliasfilter(i,j)=0;
            end
        end
    end
    aliasfilter=fftshift(aliasfilter);
    aliasfilter(1,1)=0;
end

end
