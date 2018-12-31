function [xf]=fft1d2d(x,dim,direction)
    [~,c]=size(x);
    if dim==1
        if direction==1
            xf=fft(x,c,2);
        else
            xf=ifft(x,c,2);
        end
    else
        if direction==1
            xf=fft2(x);
        else
            xf=ifft2(x);
        end
    end
        
    
end