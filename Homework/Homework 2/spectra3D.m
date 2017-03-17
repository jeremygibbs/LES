%%%% Created by Rob Stoll for ME EN 7960 008 LES of Turb Flow Spring 2009 
%%%% mfile to calculate 3D spectra from a periodic 3D velocity field 
%%%% bins is the number of wavenumber bins (shells), E is the returned
%%%% 3D Energy spectrum function [tke between wavenumbers k and k+dk] 
%%%% (vector) and Bn is the 3D wavenumber magnitude 
%%%% [sqrt(k_1^2+k_2^2+k_3^2)] bin centers (vector).

function [E,Bn]=spectra3D(U,V,W,bins)

N = size(U); %define the size of the field

dk = sqrt(3)*N(1)/2/bins; %space between two wavenumbers
Bn = (1:bins)*dk;        %define the wavenumber bins

U_hat=fftn(U); %3D Fourier transform of u
V_hat=fftn(V); %3D Fourier transform of v
W_hat=fftn(W); %3D Fourier transform of w

% Find the total squared magnitude
KE_hat = real(U_hat).^2 + imag(U_hat).^2 + ...
    real(V_hat).^2 + imag(V_hat).^2 + ...
    real(W_hat).^2 + imag(W_hat).^2;
KE_hat=1/prod(N)^2*KE_hat;

E = zeros(bins,1); %initialize the energy spectrum E

for k=1:N(1)
    kk=k-1;
    if(kk > N(1)/2); kk=kk-N(3); end
    
    for j=1:N(2)
        jj=j-1;
        if(jj > N(1)/2); jj=jj-N(2); end
        
        for i=1:N(3)
            ii=i-1;
            if(ii > N(1)/2); ii=ii-N(1); end
            
            %for each bin add up the total energy
            wv = sqrt(ii^2+jj^2+kk^2);
            cnt=0;
            for Pc=1:bins
                if(abs(Bn(Pc)-wv) < dk/2)
                    E(Pc)=E(Pc)+KE_hat(i,j,k);
                else
                    E(Pc)=E(Pc);
                end
            end
            
        end
    end
end
