package org.clas.analysis;

import org.clas.tools.EBMCEngine;
import org.clas.tools.Event;
import org.clas.tools.KinLib;
import org.clas.viewer.DetectorMonitor;

import org.jlab.clas.physics.Particle;
import org.jlab.geom.prim.Point3D;
import org.jlab.groot.math.F1D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ECelas extends DetectorMonitor {
	
	EBMCEngine  ebmce = new EBMCEngine();
    Event          ev = new Event();
	KinLib        kin = new KinLib();
	
	double kinqw[] = null, EB, beamEnergy, mp = kin.mp;
	boolean isMC = false;

	HashMap<Integer,ArrayList<Integer>> part2calo  = null;	
	
	IndexedList<List<Particle>> ecpart = new IndexedList<List<Particle>>(1);
	List<Particle>  epart  = new ArrayList<Particle>();
	List<Particle>  ppart2 = new ArrayList<Particle>();
	List<Particle>  ppart4 = new ArrayList<Particle>();
	IndexGenerator ig = new IndexGenerator();
	
    public ECelas(String name) {
        super(name);
        dgmActive=true; 
        setDetectorTabNames("WAGON",
        		            "EVENT",
        		            "KCOR",
        		            "RCOR",
        		            "XSEC");
        use123Buttons(true);
        useSliderPane(true);
        useECEnginePane(true);
        init();
        localinit();
        localclear();       
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
	    kinqw = kin.getQW(beamEnergy, 0.8*beamEnergy, 6.0);  
	    if (run<100) isMC=true;
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
    	createXSEC(0);
    	histosExist = true;
    }
    
    public void createWAGON(int st) {
    	double xlo,xhi,ylo,yhi;
    	int lc=7, lw=2, ls=1;
    	switch (st) {        
        case 0: 
        	dgm.add("WAGON",6,3,0,st,getRunNumber());         	
        	for(int is=1; is<7; is++) dgm.makeH1("wa0"+is,180, 160,200,-1,"SECTOR "+is,"#Delta#phi (deg)",  "",1,0,"11001100");
         	for(int is=1; is<7; is++) dgm.makeH1("wa1"+is,100, -1,  1, -1," ",         "#Delta EB (GeV)",   "",1,0,"11001100");
        	for(int is=1; is<7; is++) dgm.makeH1("wa2"+is,100,-20, 20, -1," ",         "#Delta#theta (deg)","",1,0,"11001100");
        	break;
        case 1:
        	dgm.add("WAGON",6,3,0,st,getRunNumber());
        	
        	xlo=kin.th2vsth1(33, beamEnergy); xhi=kin.th2vsth1(5, beamEnergy); yhi=kin.th1vsth2(xlo+2, beamEnergy); 
        	F1D f1 = new F1D("k1","th1vsth2(x,"+beamEnergy+")",xlo+3.5,Math.min(80,xhi)); f1.setLineColor(lc); f1.setLineWidth(lw); f1.setLineStyle(ls);
        	for(int is=1; is<7; is++) {dgm.makeH2("wb0"+is,80, xlo, 70, 80, 5, yhi, -1,"SECTOR "+is,"#theta p+ (deg)", "#theta e- (deg)"); 
        	dgm.addDataSet(f1,-2);}
        	
        	xlo=kin.th2vsth1(33, beamEnergy);             xhi=kin.p2vsp1(beamEnergy*0.1,beamEnergy); 
        	ylo=kin.p2vsp1(beamEnergy*0.985, beamEnergy); yhi=kin.p2vsp1(beamEnergy*0.5,beamEnergy);  
        	F1D f3 = new F1D("k3","p2vsth2(x,"+beamEnergy+")",xlo+4,67); f3.setLineColor(lc); f3.setLineWidth(lw); f3.setLineStyle(ls);
        	for(int is=1; is<7; is++) {dgm.makeH2("wb2"+is, 80, xlo, 70, 50, ylo, yhi,-1," ","#theta p+ (deg)", "MOM p+ (GeV)");  
        	dgm.addDataSet(f3,-2);}
        	
        	xlo=kin.p2vsp1(beamEnergy*0.985, beamEnergy); xhi=kin.p2vsp1(beamEnergy*0.6,beamEnergy);       	
        	F1D f2 = new F1D("k2","p1vsp2(x,"+beamEnergy+")",xlo,xhi); f2.setLineColor(lc); f2.setLineWidth(lw); f2.setLineStyle(ls);
        	for(int is=1; is<7; is++) {dgm.makeH2("wb1"+is, 50, xlo, xhi, 50, beamEnergy*0.6, beamEnergy,-1," ","MOM p+ (GeV)", "MOM e- (GeV)");  
        	dgm.addDataSet(f2,-2);}   
        	break;
        case 2:
        	dgm.add("WAGON",6,7,0,st,getRunNumber());
        	double xlo1=kin.th2vsth1(33, beamEnergy);
        	for(int is=1; is<7; is++) {dgm.makeH2("wc0"+is, 80,xlo1, 60, 40, -2.0, 2.0,0,"SECTOR "+is,"#theta p+ (DEG)", "#Delta #theta e- (DEG)");}
        	for(int is=1; is<7; is++) {dgm.makeH2("wc2"+is, 80,xlo1, 60, 40,  0.6, 1.4,1," ","#theta p+ (DEG)", "MOM / KIN p+");}
        	for(int is=1; is<7; is++) {dgm.makeH2("wc00"+is,16,  14, 30, 40, -0.5, 0.5,0," ","#theta e- (DEG)", "#Delta EB (GeV)");}
        	for(int is=1; is<7; is++) {dgm.makeH2("wc20"+is,18,  23, 42, 40, -0.5, 0.5,0," ","#theta p+ (DEG)", "#Delta EB (GeV)");}
        	xlo=kin.p2vsp1(beamEnergy*0.985, beamEnergy); xhi=kin.p2vsp1(beamEnergy*0.5,beamEnergy);  
        	for(int is=1; is<7; is++) {dgm.makeH2("wc1"+is, 50, xlo, xhi, 40, 0.8, 1.1,1," ","MOM p+ (GeV)", "MOM / KIN e-"); }   
        	for(int is=1; is<7; is++) {dgm.makeH2("wc3"+is, 18,  23,  42, 40, 0.8, 1.1,1," ","#theta p+ (deg)", "MOM / KIN #theta e-"); }   
        	for(int is=1; is<7; is++) {dgm.makeH2("wc4"+is, 40, 0.8, 1.1, 40,-0.5, 0.5,0," ","MOM / KIN #theta e-", "#Delta E (GeV)"); }   
           break;
        case 3:
        	dgm.add("WAGON",6,4,0,st,getRunNumber());        
        	for(int is=1; is<7; is++) dgm.makeH2("wd0"+is,100,10,60,4,1,5,-1,"SECTOR "+is,"#theta p+ (deg)","STATUS");
        	String str = "RED: CD=1    BLK: FD=1    GRN: FD=1 CD=0    BLU: FD=1 CD=1";
        	for(int is=1; is<7; is++) {
        		for(int ist=1; ist<5; ist++) {
        			dgm.makeH1("wd1"+is+ist,100,10,60,-(ist==1?1:2),str,"#theta p+ (deg)","",ist==4?9:ist,0,"1000000");
            		dgm.cc("wd1"+is+ist,false,false,0,400f,0,0);
        		}
        	}
    	}
    }
    
    public void createEVENT(int st) {
    	int nb = 100, is;
    	switch (st) {        
        case 0: 
        	dgm.add("EVENT",6,5,0,st,getRunNumber());
        	for(is=1; is<7; is++)  dgm.makeH2("ev0"+is, nb,0.7,kinqw[1]*1.1,60,4,40,-1,"SECTOR "+is,"W (GEV)","#theta (deg)"); 
        	for(is=1; is<7; is++)  dgm.makeH2("ev1"+is, nb,0.7,kinqw[1]*1.1,60,4,40,-1,"",          "W (GEV)","#theta (deg)"); 
        	for(is=1; is<7; is++) {dgm.makeH1("ev2a"+is,nb,0.7,kinqw[1]*1.1,-1,"","W (GeV)","",1,0,"1000000"); 
        	                       dgm.makeH1("ev2b"+is,nb,0.7,kinqw[1]*1.1,-2,"","W (GeV)","",1,1,"1000000");}
        	for(is=1; is<7; is++)  dgm.makeH2("ev3"+is, nb, 0.7,1.1,60,-8,0,-1,"","W (GeV)","Vertex (cm)"); 
        	for(is=1; is<7; is++)  dgm.makeH2("ev4"+is, 30,-0.5,0.5,30,-8,0,-1,"","#Delta EB (GeV)","Vertex (cm)"); 
        	break;
        case 1:
        	dgm.add("EVENT",3,2,0,st,getRunNumber());         	
        	for(is=1; is<7; is++) {dgm.makeH1("wth1"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=5-8",  "W (GeV)","",is==6?9:is,0,"1000000");}
        	for(is=1; is<7; is++) {dgm.makeH1("wth2"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=8-10", "W (GeV)","",is==6?9:is,0,"1000000");}
        	for(is=1; is<7; is++) {dgm.makeH1("wth3"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=10-15","W (GeV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("wth4"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=15-20","W (GeV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("wth5"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=20-25","W (GeV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("wth6"+is,nb,0.7,kinqw[1]*1.1,-(is==1?1:2),"#theta=25-30","W (GeV)","",is==6?9:is,0,"1000000");} 
        	break;
        case 2:
        	dgm.add("EVENT",3,2,0,st,getRunNumber()); 
        	double be=beamEnergy;
        	for(is=1; is<7; is++) {dgm.makeH1("pth1"+is,nb,0.45,kin.p1vsth1(6, be)*1.1,-(is==1?1:2),"#theta=6-8",  "P (GeV)","",is==6?9:is,0,"1000000");}
        	for(is=1; is<7; is++) {dgm.makeH1("pth2"+is,nb,0.45,kin.p1vsth1(8, be)*1.1,-(is==1?1:2),"#theta=8-10", "P (GeV)","",is==6?9:is,0,"1000000");}
        	for(is=1; is<7; is++) {dgm.makeH1("pth3"+is,nb,0.45,kin.p1vsth1(10,be)*1.1,-(is==1?1:2),"#theta=10-15","P (GeV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("pth4"+is,nb,0.45,kin.p1vsth1(15,be)*1.1,-(is==1?1:2),"#theta=15-20","P (GeV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("pth5"+is,nb,0.45,kin.p1vsth1(20,be)*1.1,-(is==1?1:2),"#theta=20-25","P (GeV)","",is==6?9:is,0,"1000000");} 
        	for(is=1; is<7; is++) {dgm.makeH1("pth6"+is,nb,0.45,kin.p1vsth1(25,be)*1.1,-(is==1?1:2),"#theta=25-30","P (GeV)","",is==6?9:is,0,"1000000");}         }
    }
    
    public void createKCOR(int st) {
    	switch (st) {
    	case 0:
    		dgm.add("KCOR",6,4,0,st,getRunNumber());
    		for(int is=1; is<7; is++) dgm.makeH2("kc0"+is,   50, 20,60,30,-2.0,4.0, -0,"SECTOR "+is,"#theta p+ (deg)","#Delta#phi (deg)"); 
    		for(int is=1; is<7; is++) dgm.makeH2("kc1"+is,   50,  8,31,30,-2.0,4.0, -0," ","#theta e- (deg)","#Delta#phi (deg)"); 
    		for(int is=1; is<7; is++) dgm.makeH2("kepth0"+is,25,  5,25,25,-0.2,0.2, -0," ","#theta e- (deg)","#Delta P (GeV)");    		
    		for(int is=1; is<7; is++) dgm.makeH2("kepph0"+is,25,-15,27,25,-0.2,0.2, -0," ","#phi e- (deg)",  "#Delta P (GeV)");     		
    	}
    }
    
    public void createXSEC(int st) {
    	switch (st) {
    	case 0:
    		dgm.add("XSEC",6,4,0,st,getRunNumber());
    		for(int is=1; is<7; is++) dgm.makeH2("xs1"+is,100,4,18,100,-32,32,   -0,"W < 1.065 ",    "S"+is+" #theta e- (deg)","#phi e- (deg)");     		
    		for(int is=1; is<7; is++) dgm.makeH2("xs2"+is,100,4,18,100,-32,32,   -0,"#Delta#phi < 3","S"+is+" #theta e- (deg)","#phi e- (deg)");     		
    		for(int is=1; is<7; is++) dgm.makeH2("xs3"+is,100,-220,-80,100,-100,100,-1," ",             "S"+is+" PCAL HX (cm)", "PCAL HY (cm)");     		
   		    for(int is=1; is<7; is++) dgm.makeH2("xs4"+is,100,-220,-80,100,-100,100,-1," ",             "S"+is+" PCAL X (cm)",  "PCAL Y (cm)");     		
//    		for(int is=1; is<7; is++) dgm.makeH2("xs4"+is,100,-220,-80,100,-100,100,-1," ",             "S"+is+" ECIN X (cm)",    "ECIN Y (cm)");     		
    	}    	
    }
    
    @Override       
    public void plotHistos(int run) {
    	if(!histosExist) return;
    	plot("WAGON");
    	plot("EVENT");    
    	plot("KCOR");
    	plot("XSEC");
    }
    
    @Override    
    public void processEvent(DataEvent event) {  
    	
    	ev.init(event);   	
    	ev.setHipoEvent(isHipo3Event);
    	ev.setEventNumber(getEventNumber());
    	ev.requireOneElectron(false);
    	ev.setElecTriggerSector(ev.getElecTriggerSector(shiftTrigBits(getRunNumber())));
	    
    	if(!ev.procEvent(event)) return;
    	
    	epart.clear(); ppart2.clear(); ppart4.clear();
    	
    	epart  = makePART(11,  2,0.1f);
    	ppart2 = makePART(2212,2,0.1f);
    	ppart4 = makePART(2212,4,0.1f);
    	
    	fillHists();

    }
    
    public List<Particle> makePART(int pid, int stat, float pmin) {
    	List<Particle> olist = new ArrayList<Particle>();        
        for (Particle p : ev.getPART(pmin,pid,stat)) olist.add(p);
        return olist;    	
    }

    public Point3D squeeze(Point3D xyz, int det, int e_sect) {        
        xyz.rotateZ(Math.toRadians(-60*(e_sect-1)));
        xyz.rotateY(Math.toRadians(-25.));
        xyz.translateXYZ(-(det==1?40:50),0,0);
        xyz.rotateY(Math.toRadians(25.));
//        xyz.rotateZ(Math.toRadians(60*(e_sect-1)));
        return xyz;
    }
     
    public void fillHists() {
    	
      	if (epart.size()!=1) return;
      	
      	Particle pcal = new Particle(), ecin = new Particle(), elec = new Particle();
     	
        List<Particle> e = ev.getECAL((int)epart.get(0).getProperty("pindex"));
        
        if (e.size()>0 && e.get(0).getProperty("layer")==1) {pcal = e.get(0); elec = e.get(0);}
        if (e.size()>1 && e.get(1).getProperty("layer")==4)  ecin = e.get(1);
        
        int s = (int) elec.getProperty("sector");
        
        double cthe = Math.cos(elec.theta()), sthe=Math.sqrt(1-cthe*cthe), ep=elec.p();
        double ethe = Math.toDegrees(elec.theta()), ephi = Math.toDegrees(elec.phi());
    	double   nu = beamEnergy - ep;
    	double   q2 = 2*beamEnergy*ep*(1-cthe);
    	double   w2 = -q2 + mp*mp + 2*mp*nu;
    	double  vze = elec.vz();
    	
    	if(w2<0) return;
    	
    	double w = Math.sqrt(w2);

        float lv = (float) elec.getProperty("lv");
        float lw = (float) elec.getProperty("lw");
        
        float x1=-1000,y1=-1000,x2=-1000,y2=-1000, hx1=-1000, hy1=-1000;
        
        Point3D pxyz = ev.getRotTiltPoint(new Point3D((float)pcal.getProperty("x"),(float)pcal.getProperty("y"),(float)pcal.getProperty("z")),s);
        x1 = (float) pxyz.x();
        y1 = (float) pxyz.y();
        
        Point3D hpxyz = ev.getRotTiltPoint(new Point3D((float)pcal.getProperty("hx"),(float)pcal.getProperty("hy"),(float)pcal.getProperty("hz")),s);
        hx1 = (float) hpxyz.x();
        hy1 = (float) hpxyz.y();

//        if (e.size()>1 && e.get(1).getProperty("layer")==4) {
//        	pxyz = ev.getRotTiltPoint(new Point3D((float)ecin.getProperty("x"),(float)ecin.getProperty("y"),(float)ecin.getProperty("z")),s);
//        	x2 = (float) pxyz.x();
//        	y2 = (float) pxyz.y();
//        }  
        
        boolean trig = s == getElecTriggerSector(shiftTrigBits(getRunNumber())); 
        
        double phie = (ephi<-20 ? 360+ephi : ephi)-(s-1)*60;
        
        if(true && lv>19 && lw>19) {
        	if(ethe>= 5&&ethe <8) {dgm.fill("wth1"+s,w);dgm.fill("pth1"+s,ep);}
            if(ethe>= 8&&ethe<10) {dgm.fill("wth2"+s,w);dgm.fill("pth2"+s,ep);}
            if(ethe>=10&&ethe<15) {dgm.fill("wth3"+s,w);dgm.fill("pth3"+s,ep);}
            if(ethe>=15&&ethe<20) {dgm.fill("wth4"+s,w);dgm.fill("pth4"+s,ep);}
            if(ethe>=20&&ethe<25) {dgm.fill("wth5"+s,w);dgm.fill("pth5"+s,ep);}
            if(ethe>=25&&ethe<30) {dgm.fill("wth6"+s,w);dgm.fill("pth6"+s,ep);} 
		    dgm.fill("kepth0"+s, ethe, kin.p1vsth1(ethe, beamEnergy)-ep);
		    dgm.fill("kepph0"+s, phie, kin.p1vsth1(ethe, beamEnergy)-ep); 
        }
        
        double delphi=0, delphi2=0, delphi4=0, pphi2=-1000, pthe2=-1000, pphi4=-1000, pthe4=-1000;
        
        if(ppart2.size()==1) {
        	pphi2 = Math.toDegrees(ppart2.get(0).phi());
        	pthe2 = Math.toDegrees(ppart2.get(0).theta());
        	delphi2 = ephi-pphi2;
        	if(delphi2<0) delphi2 = 360+delphi2;
        }
        
        if(ppart4.size()==1) {
        	pphi4 = Math.toDegrees(ppart4.get(0).phi());
        	pthe4 = Math.toDegrees(ppart4.get(0).theta());
        	delphi4 = ephi-pphi4;
        	if(delphi4<0) delphi4 = 360+delphi4;
        }
  
        boolean dphi = false;
        boolean dphi2 = delphi2>177 && delphi2<183;
        boolean dphi4 = delphi4>177 && delphi4<183;  
        
        if (w<1.065) {
        if(dphi4 &&                     ppart4.size()==1) {dgm.fill("wd1"+s+2,pthe4);dgm.fill("wd0"+s,pthe4,1);}
        if(dphi2 && ppart2.size()==1)                     {dgm.fill("wd1"+s+1,pthe2);dgm.fill("wd0"+s,pthe2,2);}
        if(dphi2 && ppart2.size()==1 && ppart4.size()==0) {dgm.fill("wd1"+s+3,pthe2);dgm.fill("wd0"+s,pthe2,3);}
        if(dphi2 && ppart2.size()==1 && ppart4.size()==1) {dgm.fill("wd1"+s+4,pthe2);dgm.fill("wd0"+s,pthe2,4);}
        }
        
        Particle prot = null;
        if (useFD && ppart2.size()==1) {dphi = dphi2; delphi = delphi2; prot = ppart2.get(0);} //use FD proton
        if (useCD && ppart4.size()==1) {dphi = dphi4; delphi = delphi4; prot = ppart4.get(0);} //use CD proton        
        
        dgm.fill("ev0"+s,w,ethe); 
        dgm.fill("ev2a"+s,w); 
        dgm.fill("ev3"+s,w,vze);
    	
	    if(w>1.065) return;
	    	
        double   cthp = Math.cos(prot.theta()), sthp=Math.sqrt(1-cthp*cthp), pp=prot.p(), pthe=Math.toDegrees(prot.theta());
        double   pphi = Math.toDegrees(prot.phi());      
        double    vzp = prot.vz();
        double delthe = Math.toDegrees(elec.theta() - Math.atan(pp*sthp/(beamEnergy-pp*cthp)));
        double  ebeam = (mp*(cthe + sthe*cthp/sthp)-mp)/(1-cthe); 
       	double   berr = ebeam-beamEnergy;

        dgm.fill("wa0"+s, delphi); 
        
        dgm.fill("xs1"+s, ethe, phie);
        dgm.fill("xs3"+s, hx1,hy1);
        dgm.fill("xs4"+s, x1,y1);
//        dgm.fill("xs4"+s, x2,y2);
        
        if(dphi) {    
            dgm.fill("xs2"+s, ethe, phie);
          	dgm.fill("ev1"+s,w,ethe); 
        	dgm.fill("ev2b"+s,w);
        	dgm.fill("ev4"+s,berr,vzp);
        	
        	dgm.fill("wa1"+s, berr); 
        	dgm.fill("wa2"+s, delthe);	
        	
        	dgm.fill("wb0"+s, pthe, ethe);
        	dgm.fill("wb1"+s, pp, ep);
        	dgm.fill("wb2"+s, pthe, pp);
        	
        	dgm.fill("wc0"+s, pthe, ethe-kin.th1vsth2(pthe,beamEnergy));
        	dgm.fill("wc2"+s, pthe, pp/kin.p2vsth2(pthe,beamEnergy));
        	dgm.fill("wc00"+s,ethe, berr);
        	dgm.fill("wc20"+s,pthe, berr);
         	dgm.fill("wc1"+s, pp,ep/kin.p1vsp2(pp,beamEnergy));
        	dgm.fill("wc3"+s, pthe, pp/kin.p2vsth1(ethe,beamEnergy));
        	dgm.fill("wc4"+s, pp/kin.p2vsth1(ethe,beamEnergy),berr);
        }
        
        dgm.fill("kc0"+s, pthe,delphi-180);
        dgm.fill("kc1"+s, ethe,delphi-180);

    }
    
    public void analyze() {
    	System.out.println(getDetectorName()+".Analyze() "); 
    }

    public void plot(String tabname) {     
    	dgm.drawGroup(tabname,0,getActive123(), getRunNumber());
    }   
        
}
