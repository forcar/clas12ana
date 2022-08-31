package org.clas.analysis;

import org.clas.tools.KinLib;
import org.clas.viewer.DetectorMonitor;

import org.jlab.clas.physics.Particle;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

import java.util.ArrayList;
import java.util.HashMap;

public class ECelas extends DetectorMonitor {
	
	KinLib kin = new KinLib();
	double kinqw[] = null;
	boolean doWA = true, doEV = true;
	double EB, beamEnergy;
	double mp = kin.mp;
	DataBank RecPart = null, RecCal = null;
	HashMap<Integer,ArrayList<Integer>> part2calo  = null;
	
    public ECelas(String name) {
        super(name);
        dgmActive=true; 
        this.setDetectorTabNames("WAGON",
        		                 "EVENT",
        		                 "KCOR",
        		                 "RCOR");
        this.use123Buttons(true);
        this.useSliderPane(true);
        this.init();
        this.localinit();
        this.localclear();       
    }
    
    public void localinit() {
        System.out.println("ECelas.localinit()"); 
        tl.setFitData(Fits);
    } 
    
    public void localclear() {
    	System.out.println("ECelas.localclear()");
    	isAnalyzeDone = false;
    	getDataGroup().clear();
    	runlist.clear();
    	Fits.clear();
    	FitSummary.clear();
    	tl.Timeline.clear();
    	runslider.setValue(0);
    }
    
    public void init(int run) {
	    setRunNumber(run);
	    runlist.add(run);  
	    EB = beamEnergy = getBeamEnergy(run);
	    kinqw = kin.getQW(beamEnergy, 0.2*beamEnergy, 6.0);   	
    }
    
    @Override  
    public void createHistos(int run) {
	    System.out.println("ECelas:createHistos("+run+")");
	    init(run);
    	createWAGON(0);
    	createWAGON(1);
    	createWAGON(2);
    	createWAGON(3);
    	createEVENT(0);
    	createEVENT(1);
    	createEVENT(2);
    	createKCOR(0);
    	createKCOR(1);
    	histosExist = true;
    }
    
    public void createWAGON(int st) {
    	switch (st) {        
        case 0: 
        	dgm.add("WAGON",6,3,0,st,getRunNumber());         	
        	for(int is=1; is<7; is++) dgm.makeH1("wa0"+is,180, 160,200,-1,"SECTOR "+is,"#Delta#phi (DEG)",  "",1,0,"11001100");
         	for(int is=1; is<7; is++) dgm.makeH1("wa1"+is,100, -1,  1, -1," ",         "#Delta EBEAM (GEV)","",1,0,"11001100");
        	for(int is=1; is<7; is++) dgm.makeH1("wa2"+is,100,-20, 20, -1," ",         "#Delta#theta (DEG)","",1,0,"11001100");
        	break;
        case 1:
        	dgm.add("WAGON",6,3,0,st,getRunNumber());
        	double xlo=kin.th2vsth1(33, beamEnergy), xhi=kin.th2vsth1(5, beamEnergy), yhi=kin.th1vsth2(xlo+2, beamEnergy); 
        	F1D f1 = new F1D("k1","th1vsth2(x,"+beamEnergy+")",xlo+5,Math.min(80,xhi)); f1.setLineColor(1);f1.setLineWidth(2);f1.setLineStyle(3);
        	for(int is=1; is<7; is++) {dgm.makeH2("wb0"+is,80, xlo, 80, 80, 5, yhi, -1,"SECTOR "+is,"#theta p+ (DEG)", "#theta e- (DEG)"); dgm.addDataSet(f1,-2);}
        	
        	xlo=kin.p2vsp1(beamEnergy*0.985, beamEnergy); xhi=kin.p2vsp1(beamEnergy*0.8,beamEnergy);       	
        	F1D f2 = new F1D("k2","p1vsp2(x,"+beamEnergy+")",xlo,xhi); f2.setLineColor(1);f2.setLineWidth(2);f2.setLineStyle(3);
        	for(int is=1; is<7; is++) {dgm.makeH2("wb1"+is, 50, xlo, xhi, 50, beamEnergy*0.8, beamEnergy,-1," ","MOM p+ (GeV)", "MOM e- (GeV)");  dgm.addDataSet(f2,-2);}   
        	
        	xhi=kin.th2vsth1(5, beamEnergy); double ylo=kin.p2vsp1(beamEnergy*0.985, beamEnergy); yhi=kin.p2vsp1(beamEnergy*0.8,beamEnergy); xlo=kin.th2vsp2(yhi, beamEnergy); 
        	F1D f3 = new F1D("k3","p2vsth2(x,"+beamEnergy+")",xlo+5,Math.min(80,xhi)); f3.setLineColor(1);f3.setLineWidth(2);f3.setLineStyle(3);
        	for(int is=1; is<7; is++) {dgm.makeH2("wb2"+is, 80, xlo, 80, 50, ylo, yhi,-1," ","#theta p+ (DEG)", "MOM p+ (GeV)");  dgm.addDataSet(f3,-2);}
        	break;
        case 2:
        	dgm.add("WAGON",6,3,0,st,getRunNumber());
        	double xlo1=kin.th2vsth1(33, beamEnergy);
        	for(int is=1; is<7; is++) {dgm.makeH2("wc0"+is, 80, xlo1, 80, 40, -2, 2,   0,"SECTOR "+is,"#theta p+ (DEG)", "#Delta #theta e- (DEG)");}
        	xlo=kin.p2vsp1(beamEnergy*0.985, beamEnergy); xhi=kin.p2vsp1(beamEnergy*0.8,beamEnergy);  
        	for(int is=1; is<7; is++) {dgm.makeH2("wc1"+is, 50, xlo, xhi, 40, 0.8, 1.1,1," ","MOM p+ (GeV)", "MOM / KIN e-"); }   
        	for(int is=1; is<7; is++) {dgm.makeH2("wc2"+is, 80, xlo1, 80, 40, 0.6, 1.4,1," ","#theta p+ (DEG)", "MOM / KIN p+");}
            break;
        case 3:
    		double thmin=5.2, thmax=11, wmax=1.3;
    		double pmax=kin.ep_from_w(beamEnergy, thmin, 0)*1.02, pmin=kin.ep_from_w(beamEnergy, thmax, wmax);
    		System.out.println(pmin);
        	dgm.add("WAGON",6,4,0,st,getRunNumber());        
        	for(int is=1; is<7; is++) dgm.makeH2("wd0"+is, 80,thmin,thmax,80,pmin,pmax,-1,"SECTOR "+is,"#theta neg (deg)","p (GeV)");
            for(int is=1; is<7; is++) dgm.makeH2("wd1"+is, 80,thmin,thmax,80,pmin,pmax,-1," ",          "#theta e- (deg)","p (GeV)");
            for(int is=1; is<7; is++) dgm.makeH2("wd2"+is, 80,thmin,thmax,80,pmin,pmax,-1," ",         "#theta pi- (deg)","p (GeV)");
            for(int is=1; is<7; is++) dgm.makeH2("wd3"+is, 80,thmin,thmax,80,pmin,pmax,-1," ",        "#theta null (deg)","p (GeV)");
    	}
    }
    
    public void createEVENT(int st) {
    	int nb = 100, is;
    	switch (st) {        
        case 0: 
        	dgm.add("EVENT",6,3,0,st,getRunNumber());
        	for(is=1; is<7; is++)  dgm.makeH2("ev0"+is, nb,0.7,kinqw[1]*1.1,60,4,40,-1,"SECTOR "+is,"W (GEV)","#theta (DEG)"); 
        	for(is=1; is<7; is++)  dgm.makeH2("ev1"+is, nb,0.7,kinqw[1]*1.1,60,4,40,-1,"",          "W (GEV)","#theta (DEG)"); 
        	for(is=1; is<7; is++) {dgm.makeH1("ev2a"+is,nb,0.7,kinqw[1]*1.1,-1,"","W (GEV)","",1,0,"1000000"); 
        	                       dgm.makeH1("ev2b"+is,nb,0.7,kinqw[1]*1.1,-2,"","W (GEV)","",1,1,"1000000");}
        	break;
        case 1:
        	dgm.add("EVENT",3,2,0,st,getRunNumber());         	
        	for(is=1; is<7; is++) {dgm.makeH1("wth1"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=6-8",  "W (GEV)","",is==6?9:is,0,"1000000");}
        	for(is=1; is<7; is++) {dgm.makeH1("wth2"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=8-10", "W (GEV)","",is==6?9:is,0,"1000000");}
        	for(is=1; is<7; is++) {dgm.makeH1("wth3"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=10-15","W (GEV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("wth4"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=15-20","W (GEV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("wth5"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=20-25","W (GEV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("wth6"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=25-30","W (GEV)","",is==6?9:is,0,"1000000");} 
        	break;
        case 2:
        	dgm.add("EVENT",3,2,0,st,getRunNumber()); 
        	double be=beamEnergy;
        	for(is=1; is<7; is++) {dgm.makeH1("pth1"+is,nb,0.45,kin.p1vsth1(6, be)*1.1,-(is==1?1:2),"#theta=6-8",  "P (GEV)","",is==6?9:is,0,"1000000");}
        	for(is=1; is<7; is++) {dgm.makeH1("pth2"+is,nb,0.45,kin.p1vsth1(8, be)*1.1,-(is==1?1:2),"#theta=8-10", "P (GEV)","",is==6?9:is,0,"1000000");}
        	for(is=1; is<7; is++) {dgm.makeH1("pth3"+is,nb,0.45,kin.p1vsth1(10,be)*1.1,-(is==1?1:2),"#theta=10-15","P (GEV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("pth4"+is,nb,0.45,kin.p1vsth1(15,be)*1.1,-(is==1?1:2),"#theta=15-20","P (GEV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("pth5"+is,nb,0.45,kin.p1vsth1(20,be)*1.1,-(is==1?1:2),"#theta=20-25","P (GEV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("pth6"+is,nb,0.45,kin.p1vsth1(25,be)*1.1,-(is==1?1:2),"#theta=25-30","P (GEV)","",is==6?9:is,0,"1000000");}         }
    }
    
    public void createKCOR(int st) {
    	switch (st) {
    	case 0:
    		dgm.add("KCOR",6,2,0,st,getRunNumber());
    		for(int is=1; is<7; is++) dgm.makeH2("kc0"+is,100,40,80,50,-5,5, -0,"SECTOR "+is,"#theta p+ (DEG)","#Delta#phi (DEG)"); 
    		for(int is=1; is<7; is++) dgm.makeH2("kc1"+is,100, 8,31,50,-5,5, -0," ",         "#theta e- (DEG)","#Delta#phi (DEG)"); 
    		break;
    	case 1:
    		dgm.add("KCOR",6,2,0,st,getRunNumber());
    		for(int is=1; is<7; is++) dgm.makeH2("epth0"+is,50, 5,31,50,-0.1,0.1, -0,"SECTOR "+is,"#theta e- (DEG)","#Delta P (GEV)");    		
    		for(int is=1; is<7; is++) dgm.makeH2("epph0"+is,50,-8,37,70,-0.1,0.1, -0," ",           "#phi e- (DEG)","#Delta P (GEV)");     		
    	}
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plot("WAGON");
    	plot("EVENT");    
    	plot("KCOR");
    }
    
    @Override    
    public void processEvent(DataEvent event) {   	
    	if(!processFilter(event)) return;
    	if(doWA) processWA(event);
    	if(doEV) processEV(event);
    }
    
    public boolean processFilter(DataEvent event) {
        RecCal   = event.hasBank("REC::Calorimeter") ? event.getBank("REC::Calorimeter"):null;
        RecPart  = event.hasBank("REC::Particle")    ? event.getBank("REC::Particle"):null;       
    	boolean test1 = RecPart!=null && RecPart.rows()!=0;
    	boolean test2 = RecCal !=null && RecCal.rows()!=0;
    	if(test2) part2calo = mapByIndex(RecCal,"pindex"); 
    	return test1 && test2;
    }
    
    public Particle getPart(int ipart) {
    	int    pid  = RecPart.getInt("pid",ipart);
    	double epx  = RecPart.getFloat("px",ipart);
    	double epy  = RecPart.getFloat("py",ipart);
    	double epz  = RecPart.getFloat("pz",ipart);
    	return new Particle(pid>0?pid:-11,epx,epy,epz);    	    	
    }
     
    public boolean processWA(DataEvent event) {
    	    	
    	Particle neg=null,eneg=null,pineg=null,noneg=null;
         
        ArrayList<Integer> eleCandi = new ArrayList<>();
        ArrayList<Integer>    eleEB = new ArrayList<>();
        ArrayList<Integer>    proEB = new ArrayList<>();
        ArrayList<Integer>    pimEB = new ArrayList<>();
        ArrayList<Integer>    noPID = new ArrayList<>();

        for (int ipart=0; ipart<RecPart.rows(); ipart++) {
        	
            int pid    = RecPart.getInt("pid",ipart);
            int charge = RecPart.getInt("charge",ipart);
            int status = RecPart.getInt("status",ipart);       

            final boolean isFD = (int)(Math.abs(status)/1000) == 2;

            if (isFD && charge < 0 && pid==0)   noPID.add(ipart);
            if (isFD && charge < 0 && Math.abs(pid)>0) eleCandi.add(ipart);
            if (isFD && pid==11)       eleEB.add(ipart);
            if (isFD && pid==-211)     pimEB.add(ipart);
            if (pid==2212)             proEB.add(ipart);
        }
        
        for (Integer ipart : eleCandi) {
        	neg = getPart(ipart);
    		int s = part2calo.containsKey(ipart) ? RecCal.getInt("sector", part2calo.get(ipart).get(0)):0; 
    		if(s>0) dgm.fill("wd0"+s,neg.theta()/kin.d2r,neg.p());	
        }
        
        for (Integer ipart : eleEB) {
            eneg = getPart(ipart);
    		int s = part2calo.containsKey(ipart) ? RecCal.getInt("sector", part2calo.get(ipart).get(0)):0; 
    		if(s>0) dgm.fill("wd1"+s,eneg.theta()/kin.d2r,eneg.p());	
        }
        
        for (Integer ipart : pimEB) {
            pineg = getPart(ipart);
    		int s = part2calo.containsKey(ipart) ? RecCal.getInt("sector", part2calo.get(ipart).get(0)):0; 
    		if(s>0) dgm.fill("wd2"+s,pineg.theta()/kin.d2r,pineg.p());	
        }
        
        for (Integer ipart : noPID) {
            noneg = getPart(ipart);
    		int s = part2calo.containsKey(ipart) ? RecCal.getInt("sector", part2calo.get(ipart).get(0)):0;   		
    		if(s>0) dgm.fill("wd3"+s,noneg.theta()/kin.d2r,noneg.p());	
        }
        
        if (eleCandi.size()==1 && proEB.size()==1) {
        	
        	double cthe = neg.pz()/neg.p(), sthe=Math.sqrt(1-cthe*cthe), ep=neg.p(), ethe=neg.theta()/kin.d2r;
        	
        	int    pid  = RecPart.getInt("pid", eleCandi.get(0));
        	double ppx  = RecPart.getFloat("px",proEB.get(0));
        	double ppy  = RecPart.getFloat("py",proEB.get(0));
        	double ppz  = RecPart.getFloat("pz",proEB.get(0));
        	double pp   = Math.sqrt(ppx*ppx + ppy*ppy + ppz*ppz);
        	double cthp = ppz/pp, sthp=Math.sqrt(1-cthp*cthp);
        	Particle pos = new Particle(2212,ppx,ppy,ppz);

        	double ebeam  = (mp*(cthe + sthe*cthp/sthp)-mp)/(1-cthe); 
        	double delthe = Math.acos(cthe) - Math.atan(pp*sthp/(beamEnergy-pp*cthp));
        	double delphi = (neg.phi()-pos.phi())/kin.d2r;
        	if (delphi<0) delphi = 360+delphi;
        	
        	double pthe = Math.acos(cthp)/kin.d2r;
        	
        	if (RecCal!=null) {
        		int s = part2calo.containsKey(eleCandi.get(0)) ? RecCal.getInt("sector", part2calo.get(eleCandi.get(0)).get(0)):0;         		
        		if(s>0) {dgm.fill("wa0"+s, delphi);
        		         if(delphi>171 && delphi<194) {
        		        	 dgm.fill("wa1"+s, beamEnergy-ebeam);
        		         	 dgm.fill("wa2"+s, delthe/kin.d2r);	
        		         	 dgm.fill("wb0"+s, pthe, ethe);
        		         	 dgm.fill("wb1"+s, pp, ep);
        		         	 dgm.fill("wb2"+s, pthe, pp);
        		         	 dgm.fill("wc0"+s, pthe, ethe-kin.th1vsth2(pthe,beamEnergy));
        		         	 dgm.fill("wc1"+s, pp,ep/kin.p1vsp2(pp,beamEnergy));
        		         	 dgm.fill("wc2"+s, pthe, pp/kin.p2vsth2(pthe,beamEnergy));
        		         	 dgm.fill("kc0"+s, pthe,delphi-180);
        		         	 dgm.fill("kc1"+s, ethe,delphi-180);
        		         }
                         double nu = beamEnergy - ep;
                         double q2 = 2*beamEnergy*ep*(1-cthe);
                         double w2 = -q2 + mp*mp + 2*mp*nu;
                         if (w2>0 && pid==11) {
                         	double w = Math.sqrt(w2);
                         	dgm.fill("ev1"+s,w,Math.acos(cthe)/kin.d2r); dgm.fill("ev2b"+s,w);
                         }        		         
        		}
        	}

        	return delthe > -0.026 ? true : false;
        } 
        return false; 	
    }

    public void processEV(DataEvent event) { 
    	int ist = getElecTriggerSector();
        for(int loop = 0; loop < RecPart.rows(); loop++){
        	if(RecPart.getInt("pid",loop)!=0) {
        		Particle p = new Particle(RecPart.getInt("pid", loop),
        		RecPart.getFloat("px", loop),
        		RecPart.getFloat("py", loop),
        		RecPart.getFloat("pz", loop),
        		RecPart.getFloat("vx", loop),
        		RecPart.getFloat("vy", loop),
        		RecPart.getFloat("vz", loop));
                p.setProperty("beta", RecPart.getFloat("beta", loop));
                   
                if (p.pid()==11 && part2calo.containsKey(loop)) {
                    double nu = beamEnergy - p.e();
                    double q2 = 2*beamEnergy*p.p()*(1-Math.cos(p.theta()));
                    double w2 = -q2 + mp*mp + 2*mp*nu;
                    if (w2>0) {
                    	int    s = RecCal.getInt("sector",part2calo.get(loop).get(0));
                		float lv = RecCal.getFloat("lv",  part2calo.get(loop).get(0));          
                		float lw = RecCal.getFloat("lw",  part2calo.get(loop).get(0));   
                    	double w = Math.sqrt(w2), the=p.theta()/kin.d2r, phi=p.phi()/kin.d2r;
                    	phi = (phi<-10?360+phi:phi)-(s-1)*60;
                    	dgm.fill("ev0"+s,w,the); dgm.fill("ev2a"+s,w);
                    	boolean trig = s==ist; 
                    	if(true && lv>19 && lw>19) {
                    	if(the>= 6&&the <8) {dgm.fill("wth1"+s,w);dgm.fill("pth1"+s,p.e());}
                    	if(the>= 8&&the<10) {dgm.fill("wth2"+s,w);dgm.fill("pth2"+s,p.e());}
                    	if(the>=10&&the<15) {dgm.fill("wth3"+s,w);dgm.fill("pth3"+s,p.e());}
                    	if(the>=15&&the<20) {dgm.fill("wth4"+s,w);dgm.fill("pth4"+s,p.e());}
                    	if(the>=20&&the<25) {dgm.fill("wth5"+s,w);dgm.fill("pth5"+s,p.e());}
                    	if(the>=25&&the<30) {dgm.fill("wth6"+s,w);dgm.fill("pth6"+s,p.e());} 
   		         	    dgm.fill("epth0"+s, the, kin.p1vsth1(the, beamEnergy)-p.e());
   		         	    dgm.fill("epph0"+s, phi, kin.p1vsth1(the, beamEnergy)-p.e());                    	
                    	}
                    }                   
                }
        	}
        } 
    }
    
    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(), getRunNumber());
    }   
        
}
