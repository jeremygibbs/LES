%%% Created by Rob Stoll for ME EN 7960 008 LES of Turb Flow Spring 2009 
%%% ----------------- 3D box filter in spectral domaim-------------- %%%
%%% note: f should be the ratio L/d (where d is the filter scale)    %%%
%%% f should correspond to a vector of integers with the same number %%%
%%% of elements as requested outputs
function [varargout] = box_filt_3D_real(x,f)

N=size(x);
x1 = N(1);
x2 = N(2);
x3 = N(3);
xf = x;

for III=1:nargout
    d  = x1/f(III);
    dh = d/2;
    for k=1:x3
        lk = k-dh;
        if (lk < 1); lk = lk + x3; end 
        
        bk = lk:lk+d;
        kk=find(mod(bk,x3+1)==0);
        bk(kk:end) = mod(bk(kk:end),x3+1)+1;
        
        for j=1:x2
            lj = j-dh;
            if (lj < 1); lj = lj + x2; end 
            
            bj = lj:lj+d;
            jj=find(mod(bj,x2+1)==0);
            bj(jj:end) = mod(bj(jj:end),x2+1)+1;
            
            for i=1:x1
                li = i-dh;
                if (li < 1); li = li + x1; end 
                
                bi = li:li+d;
                ii=find(mod(bi,x1+1)==0);
                bi(ii:end) = mod(bi(ii:end),x1+1)+1;
                
                xf(i,j,k) = mean(mean(mean(x(bi,bj,bk))));
            end
        end
    end
    varargout(III) = {xf};
end