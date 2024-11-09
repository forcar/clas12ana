package org.clas.tools;

public class ElasLibNew {
 
    public ElasLibNew() {
		
    }
	
    public float elas(float eb, float theta, int param) {

//  Returns the e-p differential cross section (ub/sr)
//  as a function of the electron scattering angle theta.

//  eb: beam energy (GeV)
//  theta: scattering angle (deg)
//  param: 1=dipole ff 2=Bosted ff 3=Brash ff 7=A1-Mainz

    float e_prime;		  // energy of the scattered electron
    float gep;            // electric form factor of proton
    float  gmp;           // magnetic form factor of proton
    float q2,q1;	      // momentum transfer squared
    float alpha;          // Mott scattering constant, = e^2/4pi
    float tau;            // q^2/4M^2
    float eb;			  // beam energy
    float theta;          // electron scattering angle
    float rm_p;		      // proton mass
    float s2,c2,t2,theta2;
    float sig1,sig2;

    rm_p=0.93827;

    theta2 = Math.toRadians(theta)/2.;
    s2 = sin(theta2);
    c2 = cos(theta2);
    t2 = tan(theta2);
    
//  energy of the scattered electron from recoil of the target: 
    recoil = 1.+(2.*eb/rm_p)*s2^2;
    e_prime = eb/recoil;

//  Q squared, neglecting electron mass:
    q2 = 4.*eb*e_prime*s2^2;
    
//    call ff(q2,param,gmp,gep);

//  Mott (microbarns/sr)
    rmott = 389*(1./137.)^2*c2*2/(4.*eb^2*s2^4);
    
//  recoil correction 
    rmott = rmott/recoil;
    
//  epsilon
    tau = q2/(4.*rm_p^2);
    eps = 1+2*(1+tau)*t2^2;
    eps = 1./eps;;
    
//  form factor
    sig1 = (gep^2+tau*gmp^2)/(1.+tau)+2.*tau*gmp^2*t2^2;
    sig2 = (tau*gmp^2+eps*gep^2)/eps/(1+tau);

//  cross section
    elas = rmott*sig2;
     
    return elas;
    
    }
    }
    
    public float recoil(float eb, floattheta) {
    
    real eb,theta,rm_p,rec,e_prime;
    
    data rm_p/0.93827/
    
    theta2 = theta*3.1415926/180./2.
    s2 = sin(theta2)
    rec = 1.+(2.*eb/rm_p)*s2**2
    e_prime=eb/rec
    recoil=e_prime
    
    return
    end
    
    real function elasq(eb,q2,param)
    real q2,jac  
    data pi/3.141592/
    
    cx=1-q2/(2*eb*(eb-q2/2/0.938))
    ep=eb-q2/2/0.938
c      jac = pi*(1+q2/(2*ep*0.938))/(ep*eb)
    epr=eb/(1+eb*(1-cx))
    jac=epr**2/pi
    th=acos(cx)*180./pi
    elasq=elas(eb,th,param)*jac

    end
    
    real function elasradq(eb,q2,t,wcut,param)
    real q2,jac
    data pi/3.141592/
    
    cx=1-q2/(2*eb*(eb-q2/2/0.938))
    ep=eb-q2/2/0.938
c      jac = pi*(1+q2/(2*ep*0.938))/(ep*eb)
    epr=eb/(1+eb*(1-cx))
    jac=epr**2/pi
    th=acos(cx)*180./pi
    elasradq=elasrad(eb,th,1,t,0.,wcut,param,0.)*jac
    
    end
    
    real function asym(e0,q2,param,out)
    
    real mp,nu,e_pr,tau,v_l,snth2,csth2,tnth2,theta_e,
   &v_t,v_tpr,v_tlpr,ge,gm,snth_gam,csth_gam,theta_gam,pi
    real e0,q2,a1,a2,aep,aout(3)
    real x, y,r,ge2,gm2,den,qvec,phys
    data mp/.938/
    data pi/3.14159/

    nu 	= q2/2/mp
    e_pr	= e0-nu
    if (e_pr .le. 0.) return
    tau	= q2/4/mp**2
    v_l	= 1/(1+tau)**2
    snth2	= q2/4/e0/e_pr
    csth2	= 1-snth2
    if (csth2 .le. 0. .or. snth2 .le. 0.) return
    snth2	= sqrt(snth2)
    csth2	= sqrt(csth2)
    tnth2	= snth2/csth2
    theta_e	= 360*asin(snth2)/pi
    x		= tnth2**2
    y		= 1./(1.+tau)
    v_t	= y/2 + x
    v_tpr	=  tnth2*sqrt(x+y)
    v_tlpr	= -tnth2*y/sqrt(2.)
    
    call ff(q2,param,gm,ge)
    
    ge2	= ge**2
    gm2	= gm**2
    den	= v_l*(1.+tau)*ge2/gm2 + 2.*tau*v_t
    qvec	= (1.+tau)*q2
    qvec	= sqrt(qvec)
    snth_gam	= 2*snth2*csth2*e_pr/qvec
    csth_gam	= sqrt(1. - snth_gam**2)
    theta_gam	= 180.*asin(snth_gam)/pi
    aout(1)	= 2*tau*v_tpr/den
    aout(2)	= -2*sqrt(2.*tau*(1.+tau))*ge/gm*v_tlpr/den
    aout(3)	= aout(1)*csth_gam+aout(2)*snth_gam
    asym 	= aout(out)
    
    return
    end
    
    real function getff2(q2,opt,nf)
    integer nf
    real q2,opt
    
    call ff(q2,opt,gmp,gep)
    if (nf.eq.1) getff2=gep
    if (nf.eq.2) getff2=gmp 
    
    return
    end
    
    subroutine ff(q2,opt,gmp,gep)
    real opt
    real q2,gmp,gep,gmpp,gepp
    real q1,corr1,corr2
c      vector spix
c      vector spige
c      vector spigm
    
    q1 = sqrt(q2)
    
    if (opt.eq.1) then
    
c     Gep and Gmn dipole:
    gep	= 1./(1.+q2/0.71)**2
    gmp	= 2.7928*gep
    
    elseif (opt.eq.2) then
    
c     Bosted parameterization of form factors
c     Phys. Rev. C 51, 409

    corr1 = 1+0.62*q1+0.68*q2+2.8*q1**3+0.83*q1**4
    corr2 = 1+0.35*q1+2.44*q2+0.5*q1**3+1.04*q1**4+0.34*q1**5
    gep   = 1./corr1
    gmp   = 2.7928/corr2
    
    elseif (opt.eq.3) then
    
c     Brash parameterization based on Hall A GeP/GmP measurement
c     hep-ex\0111038 PRC 65, 051001

    corr1 = 1+0.116*q1+2.874*q2+0.241*q1**3+1.006*q1**4+0.345*q1**5
    corr2 = 1-0.13*(q2-0.04)
    gep   = corr2/corr1
    gmp   = 2.7928/corr1
    
    elseif (opt.eq.4) then
    
    corr2 = 1+0.35*q1+2.44*q2+0.5*q1**3+1.04*q1**4+0.34*q1**5
    gmp   = 2.7928/corr2
    gep   = 0.0
    
    elseif (opt.eq.5) then
    
    gep = 1.041/(1+q2/0.765)**2-0.041/(1+q2/6.2)**2
    gmp = 1.002/(1+q2/0.749)**2-0.002/(1+q2/6.0)**2
    gmp = 2.7928*gmp
    
    elseif (opt.eq.6) then
    
    gep = 1.041/(1+q2/0.765)**2-0.041/(1+q2/6.2)**2
    gmp = 1.002/(1+q2/0.749)**2-0.002/(1+q2/6.0)**2
    
    gepb = exp(-0.5*((q1-0.07)/0.27)**2)+exp(-0.5*((q1+0.07)/0.27)**2)
    gep = gep - 0.3*q2*gepb
    
    gmpb = exp(-0.5*((q1-0.35)/0.21)**2)+exp(-0.5*((q1+0.35)/0.21)**2)
    gmp = 2.7928*(gmp - 0.16*q2*gmpb)
    
    elseif (opt.eq.7) then
    
c Direct extraction of Ge and Gm from A1-Mainz data (Bernauer, arXiv:1307.6227v1      
    
c      gep = divdif(spige,spix,1000,q2,2)
c      gmp = 2.7928*divdif(spigm,spix,1000,q2,2)

    endif
    
    end
    
    double precision function spence(x)
    
    double precision x
    
    if (abs(x).lt.0.1) then
      spence = x+x**2/4.
    elseif (x.gt.0.99.and.x.lt.1.01) then
      spence = pi*pi/6.
    elseif (x.gt.-1.01.and.x.lt.-0.99) then
      spence = -pi*pi/12.
    elseif (x.gt.0) then
      spence = 0.1025+sintp(x)
    else
      spence = -0.0975+sintn(x)
    endif
    
    end
    
    double precision function sintp(x)
    implicit none
    double precision x
    double precision xstep,sum,y,arg
    integer i

    xstep	= (x-.1)/100.
    sum	= 0.
    y		= 0.1-xstep/2.
    do i=1,100
      y   	= y+xstep
      arg 	= abs(1.-y)
      sum 	= sum-dlog(arg)/y
    enddo
    sintp	= sum*xstep
    
    end

    double precision function sintn(x)
    implicit none
    double precision x,xa,ystep,y,sum
    integer i

    xa	= abs(x)
    ystep	= (xa-0.1)/100.
    sum	= 0.
    y		= 0.1-ystep/2.
    do i=1,100
      y	= y+ystep
      sum	= sum-dlog(1.+y)/y
    enddo
    sintn	= sum*ystep
    
    end
    
    real function bfunc(z)
    implicit none
    real z,xi1,xi2,xi

    xi1 = alog(1440.)-2.*alog(z)/3.
    xi2 = alog(183.)-alog(z)/3.
    xi  = xi1/xi2
    bfunc = (4./3.)*(1.+(z+1.)/(z+xi)/xi2/9.)
    
    end
    
    real function elasrad(es,theta_d,izn,t1,t2,wcut,param,rc)
    
c     Returns radiated cross section using equation II.6 of
c     Mo and Tsai, Rev. Mod. Phys. 41, 205-235 (1969).

c     es: incident electron energy (GeV)
c     theta_d: scattered electron angle (deg)
c     izn: index for target nucleus (1,2,3=H,C,Ta)
c     t1: target thickness (radiation lengths)
c     t2: outgoing e- pathlength (r.l.)
c     wcut: upper limit of W (should be > 0.938+W resolution - e.g. 1.0) (GeV)

c     This program calls function elas defined above.

    double precision spence,arg
    real me,mp,pi,alpha,z,t,es,theta_d,theta,delta
    real b,cst1,eta,bt,sigunp,sigel,eel,qs
    real mpi
    real znuc,deltac(27),del_mo,delta_t,idel
    real arg11,arg15,arg19,arg23
    real epr,e1,e3,e4,gamma4,beta4
    real radcor
    real snth
    real bfunc
    real tmass(3),znucc(3)

    data pi/3.1415926/
    data mpi/.13957/
    data me/0.000511/
    data alpha/7.2993e-3/
    data amu/0.931494061/
    data tmass/1.007276470,12.00,180.94788/
    data znucc/1.,6.,73./
    
    theta	= theta_d*pi/180.

    znuc	= znucc(izn)
    mp        = tmass(izn)*amu
    z		= 1.
    wc        = mp+wcut
    
    snth	= sin(theta)
    cst1	= 1.-cos(theta)
    eel	= es/(1.+es/mp*cst1)
    epr	= es+mp-eel
    
    epcut	= (es+(mp*mp-wc*wc)/2/mp)*eel/es
    delta	= eel-epcut
    
    gamma4	= epr/mp
    beta4	= sqrt(1.-1/gamma4**2)
    e1	= es
    e3	= eel
    e4	= epr
    eta	= es/eel
    qs	= 2.*es*eel*cst1
    
    deltac(1)=28./9.-13./6.*alog(qs/me**2)
    deltac(2)= (alog(qs/me**2) - 1. + 2.*alog(eta) )
   *     * ( 2.*alog(e1/delta) -      3.*alog(eta) )
    arg=(e3-e1)/e3
    deltac(3)=-spence(arg)
    deltac(4)=znuc**2*alog(e4/mp)
    deltac(5)=znuc**2*alog(mp/eta/delta)*
   * (alog((1.+beta4)/(1.-beta4))/beta4-2.)
    deltac(6)=znuc**2*0.5/beta4*(alog((1+beta4)/(1.-beta4))*
   * alog((e4+mp)/2/mp))
    arg=sqrt((e4-mp)*(1.+beta4)/(e4+mp)/(1.-beta4))
    deltac(7)=-znuc**2*spence(-arg)/beta4
    arg=(e3-mp)/e1
    deltac(8)=znuc*spence(arg)
    arg=(mp-e3)*mp/(2.*e3*e4-mp*e1)
    deltac(9)=-znuc*spence(arg)
    arg=2.*e3*(mp-e3)/(2*e3*e4-mp*e1)
    deltac(10)=znuc*spence(arg)
    arg11=(2*e3*e4-mp*e1)/e1/(mp-2*e3)
    arg11=abs(arg11)
    deltac(11)=znuc*alog(arg11)*alog(mp/2/e3)
    arg=(e3-e4)/e3
    deltac(12)=-znuc*spence(arg)
    arg=(e4-e3)*mp/(2.*e1*e4-mp*e3)
    deltac(13)=znuc*spence(arg)
    arg=2.*e1*(e4-e3)/(2*e1*e4-mp*e3)
    deltac(14)=-znuc*spence(arg)
    arg15=(2*e1*e4-mp*e3)/e3/(mp-2*e1)
    arg15=abs(arg15)
    deltac(15)=-znuc*alog(arg11)*alog(mp/2/e1)
    arg=(e1-mp)/e1
    deltac(16)=-znuc*spence(arg)
    arg=(mp-e1)/e1
    deltac(17)=znuc*spence(arg)
    arg=2.*(mp-e1)/mp
    deltac(18)=-znuc*spence(arg)
    arg19=abs(mp/(2*e1-mp))
    deltac(19)=-znuc*alog(arg19)*alog(mp/2/e1)
    arg=(e3-mp)/e3
    deltac(20)=znuc*spence(arg)
    arg=(mp-e3)/e3
    deltac(21)=-znuc*spence(arg)
    arg=2*(mp-e3)/mp
    deltac(22)=znuc*spence(arg)
    arg23=abs(mp/(2*e3-mp))
    deltac(23)=znuc*alog(arg23)*alog(mp/2/e3)
    arg=(e1-e3)/e1
    deltac(24)=-spence(arg)
    arg=(e4-mp)*(1-beta4)/(e4+mp)/(1+beta4)
    arg=sqrt(arg)
    deltac(25)=znuc**2*spence(arg)/beta4
    arg=(e4-mp)/(e4+mp)
    arg=sqrt(arg)
    deltac(26)=-znuc**2*spence(arg)/beta4
    deltac(27)=znuc**2*spence(-arg)/beta4

    del_mo  = 0.
    do idel=1,27
       del_mo = del_mo + deltac(idel)
    enddo
    del_mo  = -alpha*del_mo/pi

    sigunp = 1.0
    if (izn.eq.1.and.rc.eq.0) sigunp  = elas(es,theta_d,param)

c     Straggling correction delta_t
    
    b		= bfunc(znuc)
    delta_s  = -0.5*b*t1*(alog(es/eta**2/delta))
    delta_p1 = -0.5*b*t1*(alog(eel/delta))
    delta_p2 = -    b*t2*(alog(eel/delta))
    delta_t  = delta_s+delta_p1+delta_p2
    
    radcor1  = 1.+del_mo+delta_t
    radcor   = exp(del_mo+delta_t)
    print *, del_mo,delta_t
    print *, es,eel,del_mo,delta_s,delta_p1,delta_p2
    print *, 1/radcor,1/radcor1

    elasrad = radcor*sigunp
    
    end


    real function elasradnew(es,theta_d,izn,t1,t2,t3,wcut,param,rc)
    
c     Returns radiated cross section using equation II.6 of
c     Mo and Tsai, Rev. Mod. Phys. 41, 205-235 (1969).

c     es: incident electron energy (GeV)
c     theta_d: scattered electron angle (deg)
c     izn: index for target nucleus (1,2,3=H,C,Ta)
c     t1: target thickness (radiation lengths)
c     t2: outgoing e- pathlength (r.l.)
c     t3: windows (has b factor included)
c     wcut: upper limit of W (should be > 0.938+W resolution - e.g. 1.0) (GeV)

c     This program calls function elas defined above.

    double precision spence,arg
    real me,mp,pi,alpha,z,t,es,theta_d,theta,delta
    real b,cst1,eta,bt,sigunp,sigel,eel,qs
    real mpi
    real znuc,deltac(27),del_mo,delta_t,idel
    real arg11,arg15,arg19,arg23
    real epr,e1,e3,e4,gamma4,beta4
    real radcor
    real snth
    real bfunc
    real tmass(3),znucc(3)

    data pi/3.1415926/
    data mpi/.13957/
    data me/0.000511/
    data alpha/7.2993e-3/
    data amu/0.931494061/
    data tmass/1.007276470,12.00,180.94788/
    data znucc/1.,6.,73./
    
    theta	= theta_d*pi/180.

    znuc	= znucc(izn)
    mp        = tmass(izn)*amu
    z		= 1.
    wc        = mp+wcut
    
    snth	= sin(theta)
    cst1	= 1.-cos(theta)
    eel	= es/(1.+es/mp*cst1)
    epr	= es+mp-eel
    
    epcut	= (es+(mp*mp-wc*wc)/2/mp)*eel/es
    delta	= eel-epcut
    
    gamma4	= epr/mp
    beta4	= sqrt(1.-1/gamma4**2)
    e1	= es
    e3	= eel
    e4	= epr
    eta	= es/eel
    qs	= 2.*es*eel*cst1
    
    deltac(1)=28./9.-13./6.*alog(qs/me**2)
    deltac(2)= (alog(qs/me**2) - 1. + 2.*znuc*alog(eta) )
   *     * ( 2.*alog(e1/delta) -      3.*alog(eta) )
    arg=(e3-e1)/e3
    deltac(3)=-spence(arg)
    deltac(4)=-znuc**2*alog(e4/mp)
    deltac(5)=znuc**2*alog(mp/eta/delta)*
   * (alog((1.+beta4)/(1.-beta4))/beta4-2.)
    deltac(6)=znuc**2*0.5/beta4*(alog((1+beta4)/(1.-beta4))*
   * alog((e4+mp)/2/mp))
    arg=sqrt((e4-mp)*(1.+beta4)/(e4+mp)/(1.-beta4))
    deltac(7)=-znuc**2*spence(-arg)/beta4
    arg=(mp-e3)/e1
    deltac(8)=znuc*spence(-arg)
    arg=(mp-e3)*mp/(2.*e3*e4-mp*e1)
    deltac(9)=-znuc*spence(arg)
    arg=2.*e3*(mp-e3)/(2*e3*e4-mp*e1)
    deltac(10)=znuc*spence(arg)
    arg11=(2*e3*e4-mp*e1)/e1/(mp-2*e3)
    arg11=abs(arg11)
    deltac(11)=znuc*alog(arg11)*alog(mp/2/e3)
    arg=(e4-e3)/e3
    deltac(12)=-znuc*spence(-arg)
    arg=(e4-e3)*mp/(2.*e1*e4-mp*e3)
    deltac(13)=znuc*spence(arg)
    arg=2.*e1*(e4-e3)/(2*e1*e4-mp*e3)
    deltac(14)=-znuc*spence(arg)
    arg15=(2*e1*e4-mp*e3)/e3/(mp-2*e1)
    arg15=abs(arg15)
    deltac(15)=-znuc*alog(arg11)*alog(mp/2/e1)
    arg=(mp-e1)/e1
    deltac(16)=-znuc*spence(-arg)
    arg=(mp-e1)/e1
    deltac(17)=znuc*spence(arg)
    arg=2.*(mp-e1)/mp
    deltac(18)=-znuc*spence(arg)
    arg19=abs(mp/(2*e1-mp))
    deltac(19)=-znuc*alog(arg19)*alog(mp/2/e1)
    arg=(mp-e3)/e3
    deltac(20)=znuc*spence(-arg)
    arg=(mp-e3)/e3
    deltac(21)=-znuc*spence(arg)
    arg=2*(mp-e3)/mp
    deltac(22)=znuc*spence(arg)
    arg23=abs(mp/(2*e3-mp))
    deltac(23)=znuc*alog(arg23)*alog(mp/2/e3)
    arg=(e1-e3)/e1
    deltac(24)=-spence(arg)
    arg=(e4-mp)*(1-beta4)/(e4+mp)/(1+beta4)
    arg=sqrt(arg)
    deltac(25)=znuc**2*spence(arg)/beta4
    arg=(e4-mp)/(e4+mp)
    arg=sqrt(arg)
    deltac(26)=-znuc**2*spence(arg)/beta4
    deltac(27)=znuc**2*spence(-arg)/beta4

    del_mo  = 0.
    do idel=7,7
       del_mo = del_mo + deltac(idel)
    enddo
    del_mo  = -alpha*del_mo/pi

    sigunp = 1.0
    if (izn.eq.1.and.rc.eq.0) sigunp  = elas(es,theta_d,param)

c     Straggling correction delta_t
    
    b	       = bfunc(znuc)
    delta_t1 = -0.5*b*t1*alog(es/eta**2/delta)
    delta_t2 = -0.5*b*t2*alog(eel/delta)
    delta_t3 = -      t3*alog(eel/delta)
    delta_t  = delta_t1+delta_t2+delta_t3
    
    radcor1  = 1.+del_mo+delta_t
    radcor   = exp(del_mo+delta_t)
c      print *, del_mo,delta_t
c      print *, es,eel,del_mo,delta_t1,delta_t2,delta_t3
c      print *, 1/radcor,1/radcor1

    elasradnew = radcor*sigunp

    end

}
