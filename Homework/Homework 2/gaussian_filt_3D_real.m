%%% ----------------- 3D gaussian filter in real domaim------------- %%%
%%% note: f should be the ratio L/d (where d is the filter scale)    %%%
%%% f should correspond to a vector of integers with the same number %%%
%%% of elements as requested outputs

function [varargout] = gaussian_filt_3D_real(U,f)
    Nx = size(U,1);
    Ny = size(U,2);
    Nz = size(U,3);
    %width = (delta-1)/2;
    width = Nx/2;
    for III=1:nargout
        delta = Nx/f(III);
    count = 1;
    for r = -width:width
        Ggauss(count) = sqrt(6/(pi*delta^2))*exp(-6*abs(r)^2/delta^2);
        count = count+1;
    end
    Gtotal = sum(Ggauss);

    x = [1+Nx/2:Nx 1:Nx 1:Nx/2];
    y = [1+Ny/2:Ny 1:Ny 1:Ny/2];
    z = [1+Nz/2:Nz 1:Nz 1:Nz/2];
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                sum1 = 0;
                for r = -width:width
                    ii = i-r+Nx/2;
                    sum1 = sum1 + Ggauss(r+width+1)*(U(x(ii),j,k));
                end
                Utmp(i,j,k) = sum1;
            end
        end
    end
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                sum1 = 0;
                for r = -width:width
                    jj = j-r+Ny/2;
                    sum1 = sum1 + Ggauss(r+width+1)*(Utmp(i,y(jj),k));
                end
            Utmp2(i,j,k) = sum1;
            end
        end
    end
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                sum1 = 0;
                for r = -width:width
                    kk = k-r+Nz/2;
                    sum1 = sum1 + Ggauss(r+width+1)*(Utmp2(i,j,z(kk)));
                end
                Ubox(i,j,k) = sum1;
            end
        end
    end
    varargout(III) = {Ubox};
    end
end
