function [aliasfilter]=antialias1d(NROW,NCOL,mdim,strict,xory)
aliasfilter=ones(NROW,NCOL);
if xory=='x'
    for j=1:NCOL
        for i=1:NROW
            if abs(j-NCOL/2)> floor(NCOL/(mdim+1))+strict
                aliasfilter(i,j)=0;
            end
        end
    end
elseif xory=='y'
    for j=1:NCOL
        for i=1:NROW
            if abs(i-NROW/2)> floor(NROW/(mdim+1))+strict
                aliasfilter(i,j)=0;
            end
        end
    end
else
    error(message('assign filter direction'))
end
aliasfilter=fftshift(aliasfilter);
aliasfilter(1)=0;