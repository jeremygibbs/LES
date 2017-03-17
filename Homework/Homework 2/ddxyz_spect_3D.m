%%% Created by Rob Stoll for ME EN 7960 008 LES of Turb Flow Spring 2009
%%% ----------------- spectral velocity gradients ------------------- %%%
function [dudx,dudy,dudz] = ddxyz_spect_3D(u)

N=size(u);
u_h = fftn(u); %take the 3D spectrum
dudx_h = u_h+sqrt(-1);
dudy_h = u_h+sqrt(-1);
dudz_h = u_h+sqrt(-1);

for k=1:N(3)
    kk=k-1;
    if(kk > N(3)/2); kk=kk-N(3); end

    for j=1:N(2)
        jj=j-1;
        if(jj > N(2)/2); jj=jj-N(2); end

        for i=1:N(1)
            ii=i-1;
            if(ii > N(3)/2); ii=ii-N(1); end

            %multiply by ki (wavenumber times imaginary number)
            if(kk == N(3)/2 | jj == N(2)/2 | ii == N(1)/2)
                dudx_h(i,j,k) = 0;
                dudy_h(i,j,k) = 0;
                dudz_h(i,j,k) = 0;
            else
                dudx_h(i,j,k) = u_h(i,j,k)*ii*sqrt(-1);
                dudy_h(i,j,k) = u_h(i,j,k)*jj*sqrt(-1);
                dudz_h(i,j,k) = u_h(i,j,k)*kk*sqrt(-1);
            end
        end
    end
end
dudx = ifftn(dudx_h,'symmetric');
dudy = ifftn(dudy_h,'symmetric');
dudz = ifftn(dudz_h,'symmetric');