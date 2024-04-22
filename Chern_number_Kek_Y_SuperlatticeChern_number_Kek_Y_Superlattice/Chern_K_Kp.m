clear
clc
%wave_fn1k
% wave_fn2k
% wave_fn3k
% wave_fn4k

Del=0.1;
A0=1.2;v0=1;vt=v0.*Del;At=A0.*Del.^2;
% here A0=eta_0, At=eta_tao
 A00=linspace(0,1,20);

% A02=linspace(0.2,1,20);
% 
% A00=[A00,A02];

dk=0.000001;

hb=1;%hbar
zi=-1;%polarisation parameter

dth=0.000001;%change in theta


 K=linspace(0.00001,30,20000);%momentun

 %we have calculated the Berry curvature for both the conduction band
 %the notation has been used throughout is psi_lamda,+

th=pi./4;
for jj=1:length(A00)
    jj
    A0=A00(jj);
    At=A0.*Del.^2;
 parfor ki=1:length(K)
     k=K(ki);
     
     psith= wave_fn2k(zi,A0,At,v0,vt,k,th,hb);
       psidkth= wave_fn2k(zi,A0,At,v0,vt,k+dk,th,hb);
       psithdth=wave_fn2k(zi,A0,At,v0,vt,k,th+dth,hb);
       
       t1=sin(th)*(psidkth'- psith')./dk+cos(th)*(psithdth'-psith')./dth./k;
        t2=cos(th)*( psidkth-psith)./dk-sin(th)*(psithdth-psith)./dth./k;
        Ak1=t1*t2;
        
        t1=sin(th)*(psidkth- psith)./dk+cos(th)*(  psithdth- psith)./dth./k;
        t2=cos(th)*(psidkth'- psith')./dk-sin(th)*(psithdth'- psith')./dth./k;
        Ak2=t2*t1; %Berry Curvature
       
       
        
    Ak=k.*(Ak1-Ak2);
    
%      end
%     
AK(ki)=2.*pi.*Ak;
     
end
  C2(jj)=imag(Simpson38_rule(K,AK)./(2.*pi));%Chern number
end
for jj=1:length(A00)
    A0=A00(jj);
    At=A0.*Del.^2;  
    jj
  parfor ki=1:length(K)
     k=K(ki);
      
    
       psith= wave_fn4k(zi,A0,At,v0,vt,k,th,hb);
       psidkth= wave_fn4k(zi,A0,At,v0,vt,k+dk,th,hb);
       psithdth=wave_fn4k(zi,A0,At,v0,vt,k,th+dth,hb);
       
       t1=sin(th)*(psidkth'- psith')./dk+cos(th)*(psithdth'-psith')./dth./k;
        t2=cos(th)*( psidkth-psith)./dk-sin(th)*(psithdth-psith)./dth./k;
        Ak1=t1*t2;
        
        t1=sin(th)*(psidkth- psith)./dk+cos(th)*(  psithdth- psith)./dth./k;
        t2=cos(th)*(psidkth'- psith')./dk-sin(th)*(psithdth'- psith')./dth./k;
        Ak2=t2*t1; % Berry Curvature
       
       
        
    Ak=k.*(Ak1-Ak2);
    
%     end

AK(ki)=2.*pi.*Ak;
     
    

end
  C4(jj)=imag(Simpson38_rule(K,AK)./(2.*pi));%Chern number
end
 plot(K,imag(AK))
 plot(A00,C2,'blue')
 hold
plot(A00,C4,'black')
% %         
% %         