c  Calculate Mittag-Leffler function by Direct Integration Method and Pade Approximation
CCCC   Written by Guoxing Lin, Clark University, Worcester, MA 01610 USA
c  In Direct Integration Method (DIM),    0 < alpha <= 2
C  In Pade Approximation, 0 < alpha < 1
ccccc  Reference for Direct Integration Method:
ccccc  Guoxing Lin, Analyze PFG anomalous diffusion via real space and phase space approaches,  Arxiv Nov 2017
ccccc  Reference for Direct Integration Method:
ccccc  C. Zeng, Y. Chen, Global Pade approximations of the generalized Mittag-Leffler function and its
ccccc  inverse, Fract. Calc. Appl. Anal. 18 (2015) 1492–1506.
cccc
cc    MLFGL(i)--- Mittag-Leffler Function calculated from Direct Integration Method
cc    MLFPD ----- Mittag-Leffler Function calculated from Pade Approximation
cc    alpha
cc    alpha1-----alpha
cc    bdelta----- total time
cc    t(j)  ----- indivadual time
cc    cg38 ------t(j)^alpha
cccc  exp(-cg38/gamma(1+alpha))----stretched exponential function

       program main
       integer*4 i,j,k,jtot,k0
       real :: alpha,alpha1,t(2402),a(2402),agal,cg38
       real ::  MLFGL(2402),mlfpd

       open(10,file='MLF_output.dat',status='new')

       write(*,*) 'Input alpha and time--bdelta'
       read(*,*) alpha,bdelta
       write(10,*) "alpha, total time=",alpha,bdelta

cccc    jtot equal 10*n   
       jtot=240*10

       if (alpha.lt.1) then
 15   write(10,258) 'j','time','time^alpha', 'MLFGL(-t^alpha)',
     $'MLFPD(-t^alpha)','exp(-t^alpha/gamma(1+alpha))'
       else
 16   write(10,259) 'j','time','time^alpha', 'MLFGL(-t^alpha)',
     $'exp(-t^alpha/gamma(1+alpha))'


       end if


       do 510 j=1, jtot
         t(j)=j*bdelta/jtot
         a(j)=-1./gamma(1.+alpha)
 510   continue
            
cccc       MLFGL(1,1)=1.
       MLFGL(1)=1.+a(1)*abs(t(1)-0.)**alpha
       cg38=(1.*bdelta/jtot)**alpha
       if (alpha.lt.1) then
       alpha1=alpha
       call amlfpd(alpha1,cg38,MLFPD,k0)
 155   write(10,261) 1, t(1),cg38, MLFGL(1),MLFPD,
     $exp(-cg38/gamma(1+alpha))
       else
 156   write(10,262) 1, t(1),cg38, MLFGL(1),
     $exp(-cg38/gamma(1+alpha))
       end if
       Do 220  j=2,jtot  
         write(*,*) j
          tga1=abs(((t(j)-0.))**alpha-(t(j)-t(1))**alpha)
          MLFGL(j)=1.+a(1)*tga1

       Do 231 k=2,j,1
           tga1=abs(((t(j)-t(k-1)))**alpha-(t(j)-t(k))**alpha)
           MLFGL(j)=MLFGL(j)+a(k)*tga1*MLFGL(k-1)
 231          continue
       cg38=(j*bdelta/jtot)**alpha
       if (alpha.lt.1) then
       alpha1=alpha
       call amlfpd(alpha1,cg38,MLFPD,k0)
 157   write(10,261) j, t(j),cg38, MLFGL(j),MLFPD,
     $exp(-cg38/gamma(1+alpha))
       else
 158   write(10,262) j, t(j),cg38, MLFGL(j),
     $exp(-cg38/gamma(1+alpha))
       end if
 220       continue
 258   format(A5,A11,A21,A18,A16,A29)
 259   format(A5,A11,A21,A16,A29)
 261   format(I5,5F16.8)
 262   format(I5,4F16.8)          
          end

            subroutine Amlfpd (alpha,z0,eab1,k10)
ccccc  based on the article C. Zeng, Y. Chen, Global Pade approximations of the generalized Mittag-Leffler function and its
ccccc  inverse, Fract. Calc. Appl. Anal. 18 (2015) 1492–1506.
            real :: alpha,eab1,z0
            real  ::  agp,agm,agm2,q0,q1
            integer*4 k10 
            agp=gamma(1.+alpha)
            agm=gamma(1.-alpha)
            agm2=gamma(1.-2*alpha)
            q0=(agp/agm-agp*agm/agm2)/(agp*agm-1.)
            q1=(agp-agm/agm2)/(agp*agm-1.)
            eab1=(1+z0/agm/q0)/(1.+q1*z0/q0+z0*z0/q0)
            end




