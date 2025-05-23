package org.clas.tools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.Vector3;
import org.jlab.detector.base.DetectorLayer;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;
import org.jlab.geom.base.Detector;
import org.jlab.geom.base.Layer;
import org.jlab.geom.component.ScintillatorPaddle;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Point3D;

public class Event {
     
    public int run=0,event=0,unixtime=0;
    public long timestamp=0,trigger=0;
    public float starttime=0, RF_TIME=0;
	
	private DataEvent      de = null;
	private DataBank partBank = null;
	private DataBank caloBank = null;
	private DataBank ftofBank = null;
	private DataBank clusBank = null;
	private DataBank caliBank = null;
	private DataBank caleBank = null;
	private DataBank htccBank = null;
	private DataBank peakBank = null;
	public  DataBank trajBank = null;
	public  DataBank   mcBank = null;
	
	public int TRpid=0, MCpid=0;
	private boolean isHipo3Event = false; 
	public boolean     isMC = false;	
	public boolean   isMuon = false;
	public boolean   isPhys = false;
	
	public IndexedList<List<Particle>> partmap = new IndexedList<List<Particle>>(1);   //getRECParticle(pid...))
	public List<Particle>              part    = new ArrayList<Particle>();	           //getPART()
	Map<Integer,List<Integer>>         partMap = new HashMap<Integer,List<Integer>>(); //loadMapByIndex(bank,pid)
	Map<Integer,List<Integer>>         caloMap = new HashMap<Integer,List<Integer>>();
	Map<Integer,List<Integer>>         htccMap = new HashMap<Integer,List<Integer>>();
	Map<Integer,List<Integer>>         ftofMap = new HashMap<Integer,List<Integer>>();
	public Map<Integer,List<Integer>>  trajMap = new HashMap<Integer,List<Integer>>();
	public Map<Integer,List<Integer>>  trajDet = new HashMap<Integer,List<Integer>>();
	
	private boolean hasRUNconfig = false;
	private boolean hasRECevent  = false;
	private boolean hasRECscintillator = false;
	private boolean hasRECcalorimeter = false;
	private boolean hasRECcaloextras = false;	
	private boolean hasRECcherenkov = false;
	private boolean hasECALclusters = false;
	private boolean hasECALcalib = false;
	private boolean hasECALpeaks = false;
	private boolean hasRECparticle = false;
	private boolean hasRECtrack = false;
	private boolean hasRECtraj = false;
	private boolean hasMCParticle = false;
	private boolean hasHTCC = false;
	private boolean hasHTCCrec = false;
	
	private boolean requireOneElectron = false;
	
	int[]        tb = new int[32];
	public int[] N1 = new int[6];
	public int[] N2 = new int[6];
	
	public  int trigger_sect=0;
	public  int eventNumber=0;
	private int nelec=0;
	public boolean debug = false;
	
	private float timeshift = 0f;
	private int  startTimeCut = -100;
	
	public float tpol = 0;
	public float spol = 0;
	
    public List<Particle> pmc = new ArrayList<>();
    public List<Vector3>  pmv = new ArrayList<>();
	
    Detector ecDetector =  null;
    
	public Event() {
		System.out.println("Event() instantiated");
	}
	
	public void init(DataEvent event) {	
		this.de = event;
		trigger = 0;
		starttime = -100; 
		part.clear();
		partmap.clear();
		partMap.clear();
		ftofMap.clear();
		caloMap.clear();
		htccMap.clear();
		trajMap.clear();
		trajDet.clear();
		partBank = null;
		caloBank = null;
		ftofBank = null;
		clusBank = null;
		caliBank = null;
		peakBank = null;		
		mcBank   = null;

		hasRUNconfig       = de.hasBank("RUN::config");
		hasRECcalorimeter  = de.hasBank("REC::Calorimeter");
		hasRECcaloextras   = de.hasBank("REC::CaloExtras");
		hasECALclusters    = de.hasBank("ECAL::clusters");	
		hasRECevent        = de.hasBank("REC::Event");
		hasRECscintillator = de.hasBank("REC::Scintillator");
		hasRECcherenkov    = de.hasBank("REC::Cherenkov");
		hasECALpeaks       = de.hasBank("ECAL::peaks");
		hasECALcalib       = de.hasBank("ECAL::calib");	
		hasRECparticle     = de.hasBank("REC::Particle");
		hasRECtrack        = de.hasBank("REC::Track");
		hasRECtraj         = de.hasBank("REC::Traj");
		hasHTCC            = de.hasBank("HTCC::adc");
		hasHTCCrec         = de.hasBank("HTCC::rec");
		hasMCParticle      = de.hasBank("MC::Particle");
		
		if(hasRUNconfig) processRUNconfig();						
		isMC = hasMCParticle;
	}
	
	private boolean isGoodEvent() {
		isMuon   =  isMuon() && hasECALpeaks && hasECALcalib;
		isPhys   = !isMuon && hasRUNconfig && hasRECevent && hasRECparticle;		
		return isPhys || isMuon;		
	}
	
	private boolean isMuon() { //MCpid (MC data) or cosmic run configuration
		return MCpid==-13 || (hasECALclusters && !hasRECcalorimeter);
	}
	
	private boolean isGoodRows() {
		return caloBank.rows()==clusBank.rows();
	}
	
	private void fail(int val) {
		if(val==0) return;
		System.out.println(" "); System.out.println("EVENT: "+eventNumber+" FAIL "+val);		
	}
	
	public boolean procEvent(DataEvent event) {
		
        if(hasMCParticle)          processMCparticle();
        
        if(!isGoodEvent()) return false;
	    
	    if(hasRECevent)            processRECevent();
	    if(requireOneElectron && !countElectronTriggers(false)) return false;
	    if(!isMuon  && starttime > startTimeCut) {
	    	if(hasRECscintillator) processRECscintillator();
	    	if(hasRECcalorimeter)  processRECcalorimeter();
	    	if(hasRECcherenkov)    processRECcherenkov();
	    	if(hasECALclusters)    processECALclusters();
	    	if(hasECALclusters && hasRECcalorimeter && !isGoodRows()) return false;
	    	if(hasECALpeaks)       processECALpeaks();
	    	if(hasECALcalib)       processECALcalib(); //pass1
	    	if(hasRECcaloextras)   processRECcaloextras(); //pass2
	    	if(hasRECtrack)        processRECtrack();
	    	if(hasRECtraj)         processRECtraj();
	    	if(hasRECparticle)     processRECparticle();
	    	
			return true;
		}
	    if(isMuon) {
			initpartmap(13); //initialize with cosmic muon (pid=13) sign ignored!
	    	processECALclusters();
	    	processECALpeaks();
	    	processECALcalib();
	    	return true;
	    }
	    return false;
	}

	private void initpartmap(int pid) {
		partmap.add(new ArrayList<Particle>(),pid);
		Particle p = new Particle(pid, 0,0,0,0,0,0);                         			
		p.setProperty("status", 2000);
		p.setProperty("pindex", 0);
		p.setProperty("beta",0);
		p.setProperty("chi2pid", 0);
		p.setProperty("ppid", pid);
		partmap.getItem(pid).add(p);
		part.add(p);
	}
	
	public void setGeometry(Detector ecDetector) {
		this.ecDetector = ecDetector; 
	}
	
	public void setMCpid(int val) {
		MCpid = val;
	}
	
	public void setHipoEvent(boolean val) {
		this.isHipo3Event = val;
	}
	
	public void setTrigger(int trig) {
		TRpid = trig;
	}
	
	public void setElecTriggerSector(int val) {
		trigger_sect = val;
	}
	
	public void setEventNumber(int val) {
		eventNumber = val;
	}
	
	public void requireOneElectron(boolean val) {
		requireOneElectron = val;
	}
	
	public void setTimeShift(float val) {
		timeshift = val;
	}
	
	public void setStartTimeCut(int val) {
		startTimeCut = val;
	}
	
	private void processRUNconfig() {		
		storeRUNconfig(de.getBank("RUN::config"));
	}
	
	private void processRECevent() {		
		storeRECevent(de.getBank("REC::Event"));  	
	}	
	
	private void processRECparticle() {		
	 	storeRECparticle(de.getBank("REC::Particle"));	 
	}
	
	private void processRECscintillator() {		
		storeRECscintillator(de.getBank("REC::Scintillator"));				
	}	
	
	private void processRECcalorimeter() {		
		storeRECcalorimeter(de.getBank("REC::Calorimeter"));  		
	}
	
	private void processRECcaloextras() {		
		storeRECcaloextras(de.getBank("REC::CaloExtras"));  		
	}
	
	private void processRECcherenkov() {		
		storeRECcherenkov(de.getBank("REC::Cherenkov"));  		
	}
	
	private void processRECtrack() {		
		storeRECtrack(de.getBank("REC::Track"));		
	}
	
	private void processRECtraj() {		
		storeRECtraj(de.getBank("REC::Traj"));		
	}	
	
	private void processECALclusters() {
		storeECALclusters(de.getBank("ECAL::clusters"));		
	}	
	
	private void processECALcalib() {		
		storeECALcalib(de.getBank("ECAL::calib"));		
	}
	
	private void processECALpeaks() {		
		storeECALpeaks(de.getBank("ECAL::peaks"));		
	}
	
	private void processMCparticle() {
		storeMCParticle(de.getBank("MC::Particle"));
	}
	
	private void storeRUNconfig(DataBank bank) {
		this.run       = bank.getInt("run",0);
		this.event     = bank.getInt("event",0);
		this.unixtime  = bank.getInt("unixtime",0);
		this.trigger   = bank.getLong("trigger",0);
		this.timestamp = bank.getLong("timestamp",0);	
		this.tpol      = -bank.getFloat("torus",0);
		this.spol      = -bank.getFloat("solenoid",0);
	}
	
	private void storeRECevent(DataBank bank) {
		this.starttime = isHipo3Event ? bank.getFloat("STTime", 0):
                                        bank.getFloat("startTime", 0);	
		this.RF_TIME   = bank.getFloat("RFTime",0);
	}
	
	private void storeRECparticle(DataBank bank) {
		partBank = bank;
		getPART();                                            // load List<Particle> part from REC::Particle
		partMap  = loadMapByIndex(partBank,"pid");	          // load partMap which maps pindex to pid from REC::Particle
    	getRECparticle(0,-11,11,22,2112,2212,-2212,211,-211); // load partmap for all PIDs from partMap
	}
	
	private void storeRECcalorimeter(DataBank bank) {
		caloBank = bank;
		caloMap  = loadMapByIndex(caloBank,"pindex"); // "pindex" mapped to REC::Calorimeter index
	}
		
	private void storeRECcherenkov(DataBank bank) {
		htccBank = bank;
		htccMap  = loadMapByIndex(htccBank,"pindex");
	}
	
	private void storeRECscintillator(DataBank bank) {	
		ftofBank = bank;
		ftofMap  = loadMapByIndex(ftofBank,"pindex");		
	}	
	
	private void storeRECtrack(DataBank bank) {	
		
	}
	
	private void storeRECtraj(DataBank bank) {	
		trajBank = bank;
		trajMap  = loadMapByIndex(trajBank,"pindex");
		trajDet  = loadMapByIndex(trajBank,"detector");
	}	
	
	private void storeECALclusters(DataBank bank) {	
		clusBank = bank;
	}
	
	private void storeRECcaloextras(DataBank bank) {
		caleBank = bank;
	}
	
	private void storeECALcalib(DataBank bank) {	
		caliBank = bank;		
	}
	
	private void storeECALpeaks(DataBank bank) {	
		peakBank = bank;		
	}
	
	private void storeMCParticle(DataBank bank) {	
		mcBank = bank;		
        pmc.clear(); pmv.clear();
        for (int i=0; i<bank.rows(); i++) {
        	double px = bank.getFloat("px",i);
            double py = bank.getFloat("py",i);
            double pz = bank.getFloat("pz",i);
            double vx = bank.getFloat("vx",i);
            double vy = bank.getFloat("vy",i);
            double vz = bank.getFloat("vz",i);
            int   pid = bank.getInt("pid",i);     
//            if(pid==MCpid) {pmc.add(new Particle(pid, px, py, pz, vx, vy, vz)); pmv.add(new Vector3(px,py,pz));}
            pmc.add(new Particle(pid, px, py, pz, vx, vy, vz)); 
            pmv.add(new Vector3(px,py,pz));
         }            
	}

	private void getPART() {
        for(int i = 0; i < partBank.rows(); i++){           	
            int      pid = partBank.getInt("pid", i);              
            float     px = partBank.getFloat("px", i);
            float     py = partBank.getFloat("py", i);
            float     pz = partBank.getFloat("pz", i);
            float     vx = partBank.getFloat("vx", i);
            float     vy = partBank.getFloat("vy", i);
            float     vz = partBank.getFloat("vz", i);
            float   beta = partBank.getFloat("beta", i);
            float   chi2 = Math.abs(partBank.getFloat("chi2pid", i));
            float   chrg = partBank.getByte("charge", i);
            short status = (short) Math.abs(partBank.getShort("status", i));
                       
            Particle p = new Particle();
            
            if(pid==0) p.initParticleWithMass(0, px, py, pz, vx, vy, vz);
            if(pid!=0) p.initParticle(      pid, px, py, pz, vx, vy, vz);
            p.setProperty("chrg", chrg);
            p.setProperty("ppid", pid);
            p.setProperty("status", status);
            p.setProperty("beta", beta);
            p.setProperty("index",i);                        
            p.setProperty("chi2pid",chi2);
            
            part.add(i,p);     //Lists do not support sparse indices !!!!!            
        }                		
	}
	
	//returns all entries with momentum > thr and requested pid, status, charge
    public List<Particle> getPART(int pid, double thr, int ist, int q) { 
    	List<Particle> olist = new ArrayList<Particle>(); 
    	for (Particle p : getParticle(pid)) {
    		int status = (int) (Math.abs(p.getProperty("status"))/1000);
    		int   chrg = (int) p.getProperty("chrg");
    		System.out.println("PART"+pid+" status "+status+" chrg "+chrg+" p "+(float)p.p()
    		+" th "+(float)Math.toDegrees(p.theta())+" ph "+(float)Math.toDegrees((p.phi())));
    		if(status==ist && p.p()>=thr && chrg==q) olist.add(p); 
    	}          	
       return olist;    	
    }
	
	//returns all entries with momentum > thr and requested pid, status
    public List<Particle> getPART(int pid, double thr, int ist) { 
    	List<Particle> olist = new ArrayList<Particle>();    
    	for (Particle p : getParticle(pid)) {
    		int status = (int) (Math.abs(p.getProperty("status"))/1000);
    		if(status==ist && p.p()>=thr) olist.add(p); 
    	}          	
       return olist;    	
    }
    
    //returns REC::Particle entries with momentum > thr and requested pid
    public List<Particle> getPART(int pid, double thr) { 
    	List<Particle> olist = new ArrayList<Particle>();    
    	for (Particle p : getParticle(pid)) {
    		int status = (int) (Math.abs(p.getProperty("status"))/1000);
    		if(status==2 && p.p()>=thr) olist.add(p); 
    	}          	
       return olist;    	
    }
    
    //returns REC::Particle entries with requested pid	
	public List<Particle> getParticle(int ipid) { 
		List<Particle> olist = new ArrayList<Particle>();
	    IndexGenerator ig = new IndexGenerator();                
	    for (Map.Entry<Long,List<Particle>>  entry : partmap.getMap().entrySet()){
	           int pid = ig.getIndex(entry.getKey(), 0);  
	           if(ipid==pid) for (Particle p : entry.getValue()) olist.add(p);
	    }	
	    return olist;
	}

	//Use partMap to convert REC::Particle to partmap with tpid index. 
	private void getRECparticle(int...ttpid) {
		for (int tpid : ttpid) {
		if(partMap.containsKey(tpid)) {
			for(int ipart : partMap.get(tpid)){  //retrieve tpid indexed pindex to REC::Particle
				int      pid = partBank.getInt("pid",  ipart);    
				float     px = partBank.getFloat("px", ipart);
				float     py = partBank.getFloat("py", ipart);
				float     pz = partBank.getFloat("pz", ipart);
				float     vx = partBank.getFloat("vx", ipart);
				float     vy = partBank.getFloat("vy", ipart);
				float     vz = partBank.getFloat("vz", ipart);
				float   beta = partBank.getFloat("beta", ipart);
	            float   chi2 = Math.abs(partBank.getFloat("chi2pid", ipart));
	            float   chrg = partBank.getByte("charge", ipart);
				short status = (short) Math.abs(partBank.getShort("status", ipart));
			
				Particle p = new Particle(); 
				
	            if(pid==0) p.initParticleWithMass(0, px, py, pz, vx, vy, vz);
	            if(pid!=0) p.initParticle(      pid, px, py, pz, vx, vy, vz);

	            p.setProperty("chrg", chrg);
				p.setProperty("status", status);
				p.setProperty("pindex", ipart);
				p.setProperty("beta", beta);
				p.setProperty("chi2pid", chi2);

				int ip = pid<0 ? Math.abs(pid)+1 : pid;	//index ip must be +			
				if(!partmap.hasItem(ip)) {partmap.add(new ArrayList<Particle>(),ip);} 
				    partmap.getItem(ip).add(p);
			}		
		}
		}
	}
	
	public void reportElectrons(String tag) {
		nelec++; System.out.println("Evnt "+eventNumber+" Nelec "+nelec+" "+tag);
	}
	
	private float getBeta(Particle p) {
		double path = p.getProperty("path");
		double time = p.getProperty("time");
		return (float) (path/(time-starttime-timeshift)/29.979f);
	}	
	
    public int     getFDTrigger()            {return (int)(trigger)&0x000000000ffffffff;}
    public int     getCDTrigger()            {return (int)(trigger>>32)&0x00000000ffffffff;}
    public boolean isGoodFD()                {return  getFDTrigger()>0;}    
    public boolean isTrigBitSet(int bit)     {int mask=0; mask |= 1<<bit; return isTrigMaskSet(mask);}
    public boolean isTrigMaskSet(int mask)   {return (getFDTrigger()&mask)!=0;}	
        
    public int getElecTriggerSector(Boolean shift) { //shift: true=outbending e- false=inbending e-
    	int[] tb = new int[32];
    	int tbsum=0, ts=0;
    	for (int i = 31; i >= 0; i--) {tb[i] = ((trigger & (1 << i))!=0)?1:0;}
    	for (int i=shift?8:1; i<(shift?14:7); i++) {
    		tbsum+=tb[i]; if(tb[i]>0) ts=shift?i-7:i;
    	}
    	return (tbsum==0||tbsum>1) ? 0:ts;
    }
    
    public int getElecTriggerSector(int f) { //shift: 0,1,2
    	int shift=1+f*7;
    	int[] tb = new int[32];
    	int tbsum=0, ts=0;
    	for (int i = 31; i >= 0; i--) {tb[i] = ((trigger & (1 << i))!=0)?1:0;}
    	for (int i=shift; i<shift+6; i++) {
    		tbsum+=tb[i]; if(tb[i]>0) ts=i-f*7;
    	}
    	return (tbsum==0||tbsum>1) ? 0:ts;
    }
	
	public boolean countElectronTriggers(boolean debug) {
		int tbsum=0;
		for (int i = 31; i >= 0; i--) {tb[i] = ((trigger & (1 << i))!=0)?1:0;}
//		for (int i = 31; i >= 0; i--) {tb[i] = (isTrigBitSet(i))?1:0;}	
		for (int i=1; i<7; i++) tbsum+=tb[i];
		if(debug) {
			System.out.println(" ");
			System.out.print("Run "+run+" Event "+event+" Trigger Sector "+trigger_sect);
			System.out.println(" Bits "+tb[0]+" "+tb[1]+tb[2]+tb[3]+tb[4]+tb[5]+tb[6]+" "+tbsum);
			if(tbsum==1) {for (int i=1; i<7; i++) {if(tb[i]==1) N1[i-1]++;}}
			if(tbsum>1)  {for (int i=1; i<7; i++) {if(tb[i]==1) N2[i-1]++;}}
			System.out.println(N1[0]+" "+N1[1]+" "+N1[2]+" "+N1[3]+" "+N1[4]+" "+N1[5]);
			System.out.println(N2[0]+" "+N2[1]+" "+N2[2]+" "+N2[3]+" "+N2[4]+" "+N2[5]);
		}		
		if(tbsum==0||tbsum>1) return false;
		return true;			
	}	

	public List<Particle> getFTOF(int ipart) {
		List<Particle> ftofpart = new ArrayList<Particle>();
		if(ftofMap.containsKey(ipart)) {
			for(int imap : ftofMap.get(ipart)) {				
				Particle p = new Particle();                         
				p.setProperty("sector", ftofBank.getByte("sector",imap)); 
				p.setProperty("layer",  ftofBank.getByte("layer",imap));
				p.setProperty("strip",  ftofBank.getShort("component",imap));
				p.setProperty("index",  ftofBank.getShort("index",imap));
				p.setProperty("energy", ftofBank.getFloat("energy",imap));
				p.setProperty("time",   ftofBank.getFloat("time",imap));
				p.setProperty("path",   ftofBank.getFloat("path",imap));                   		                    
				p.setProperty("x",      ftofBank.getFloat("x",imap));
				p.setProperty("y",      ftofBank.getFloat("y",imap)); 
				p.setProperty("z",      ftofBank.getFloat("z",imap)); 
				p.setProperty("hx",     ftofBank.getFloat("hx",imap));
				p.setProperty("hy",     ftofBank.getFloat("hy",imap));
				p.setProperty("hz",     ftofBank.getFloat("hz",imap));
				
				if(trajMap.containsKey(ipart)) {
					for(int tmap : trajMap.get(ipart)) {
					   if(trajBank.getInt("detector",tmap)==12&&trajBank.getInt("layer",tmap)==(p.getProperty("layer"))) {
							p.setProperty("tx",     trajBank.getFloat("x",tmap));
							p.setProperty("ty",     trajBank.getFloat("y",tmap));
							p.setProperty("tz",     trajBank.getFloat("z",tmap));
							p.setProperty("cz",     trajBank.getFloat("cz",tmap));
					   }
					}
				}				
				Point3D res = getResidual(p);
				p.setProperty("resx",res.x());p.setProperty("resy",res.y());p.setProperty("resz",res.z());
				ftofpart.add(p);
			}
		}
		return ftofpart;
	}
	
	public List<Particle> getHTCC(int ipart) {
		List<Particle> htccpart = new ArrayList<Particle>();
		if(htccMap.containsKey(ipart)) {
			for(int imap : htccMap.get(ipart)) {				
				Particle p = new Particle();                         
				p.setProperty("sector", htccBank.getByte("sector",imap)); 
				p.setProperty("index",  htccBank.getShort("index",imap));
				p.setProperty("det",    htccBank.getInt("detector",imap));
				p.setProperty("nphe",   htccBank.getFloat("nphe",imap));
				p.setProperty("time",   htccBank.getFloat("time",imap));
				p.setProperty("path",   htccBank.getFloat("path",imap));                   		                    
				p.setProperty("x",      htccBank.getFloat("x",imap));
				p.setProperty("y",      htccBank.getFloat("y",imap)); 
				p.setProperty("z",      htccBank.getFloat("z",imap)); 
			
			    if(trajMap.containsKey(ipart)) {
			    	for(int tmap : trajMap.get(ipart)) {
			    		if(trajBank.getInt("detector",tmap)==15) {
			    			p.setProperty("tx",     trajBank.getFloat("x",tmap));
			    			p.setProperty("ty",     trajBank.getFloat("y",tmap));
			    			p.setProperty("tz",     trajBank.getFloat("z",tmap));
			    		}
			    	}
			    }
				htccpart.add(p);
			}			    
		}
		return htccpart;
	}
	
	public int getECALMULT(int is, int layer) {
		int n=0;
		for(int i = 0; i < caloBank.rows(); i++){  
			 if(is==caloBank.getByte("sector",i) && layer==caloBank.getByte("layer", i)) n++;
		}
		return n;
	}
	
    public int getDet(int layer) {
	    int[] il = {0,0,0,1,1,1,2,2,2}; // layer 1-3: PCAL 4-6: ECinner 7-9: ECouter  
	    return il[layer-1];
	}

	public int[][][] getECALPID(int e_sect) {
		int[][][] out = new int[6][3][6];
		for(int i = 0; i < caloBank.rows(); i++){
			 int is  = caloBank.getByte("sector",i);
			 int il  = getDet(caloBank.getByte("layer",i));
			 int pid = partBank.getInt("pid", caloBank.getShort("pindex",i));
			 boolean good_n = pid==2112;
			 boolean good_g = pid==22;
			 boolean good_neut = good_n || good_g;
			 boolean good_char = pid!=0 && !good_neut;
		     int ipid=-1;
			 if(good_char&&is==e_sect) ipid=4;
			 if(good_n   &&is==e_sect) ipid=0;
			 if(good_g   &&is==e_sect) ipid=1;
			 if(good_char&&is!=e_sect) ipid=5;
			 if(good_n   &&is!=e_sect) ipid=2;
			 if(good_g   &&is!=e_sect) ipid=3;
			 if(ipid != -1) out[ipid][il][is-1]++;
		}
		return out;		
	}
    
    // returns list of ECAL particle clusters indexed by sector 
    public IndexedList<List<Particle>> getECALParticles(List<Particle> list) {

    	IndexedList<List<Particle>> olist = new IndexedList<List<Particle>>(1);  
    	
    	for (Particle p : list) {
    		int ip = (int)p.getProperty("pindex");    		
    		int is = (int) getECAL(ip).get(0).getProperty("sector"); //only first ECAL cluster mapped to particle is returned
    		if (!olist.hasItem(is)) olist.add(new ArrayList<Particle>(), is); 
			     olist.getItem(is).add(p);
    	}
    	return olist;
    }
    
    // returns list of ECAL particle clusters indexed by sector, layer 
    public IndexedList<List<Particle>> getECALClusters(List<Particle> list) { 
    	
        IndexedList<List<Particle>> olist = new IndexedList<List<Particle>>(2); //sector,layer 
        
    	for (Particle p : list) {
    		int ip = (int)p.getProperty("pindex");
    		for (Particle ec : getECAL(ip)) { // all ECAL clusters mapped to particle are returned
    			int   is = (int)   ec.getProperty("sector");
    			int   il = (int)   ec.getProperty("layer");
    			if (!olist.hasItem(is,il)) olist.add(new ArrayList<Particle>(), is,il); 
    			     olist.getItem(is,il).add(ec);
    		}
    	}
    	return olist;
    }
    
    // filter ECAL clusters by PID and multiplicity, indexed by sector
    public IndexedList<List<Particle>> filterECALClusters(IndexedList<List<Particle>> list, int n, int pid) {                                                                                                     	
        IndexedList<List<Particle>> olist = new IndexedList<List<Particle>>(1);       
		IndexGenerator ig = new IndexGenerator();
    	for (Map.Entry<Long,List<Particle>>  entry : list.getMap().entrySet()){
			int is = ig.getIndex(entry.getKey(), 0); 
			int il = ig.getIndex(entry.getKey(), 1); 
			int ip = (int) entry.getValue().get(0).getProperty("pindex");
			if(entry.getValue().size()>0 && entry.getValue().size()<n && Math.abs(part.get(ip).getProperty("ppid"))==pid) {
				if(!olist.hasItem(is)) olist.add(new ArrayList<Particle>(), is); 
				for (Particle p : entry.getValue()) olist.getItem(is).add(p);
			}
    	}
    	return olist;    	
    }
		
	public List<Particle> getECAL(int ipart) {
		return isMuon ? getECALMUON():getECALPHYS(ipart);		
	}
	
	public List<Particle> getECALMUON() { //for cosmic muon runs
		List<Particle> ecalpart = new ArrayList<Particle>();		
		for (int i = 0; i<clusBank.rows(); i++) {
			Particle p = new Particle(); 
			p.copy(part.get(0)); //initialized with pid=13 for isMuon=true
			p.setProperty("sector", clusBank.getByte("sector",i)); 
			p.setProperty("layer",  clusBank.getByte("layer",i));
			p.setProperty("pindex", 0);
			p.setProperty("index",  i);
			p.setProperty("energy", clusBank.getFloat("energy",i)*1e3);
			p.setProperty("time",   clusBank.getFloat("time",i));
			p.setProperty("path",   0);                   		                    
			p.setProperty("x",      clusBank.getFloat("x",i));
			p.setProperty("y",      clusBank.getFloat("y",i)); 
			p.setProperty("z",      clusBank.getFloat("z",i)); 
			p.setProperty("hx",     0);
			p.setProperty("hy",     0);
			p.setProperty("hz",     0);
			p.setProperty("lu",     0);
			p.setProperty("lv",     0);
			p.setProperty("lw",     0);
			p.setProperty("du",     clusBank.getFloat("widthU",i));
			p.setProperty("dv",     clusBank.getFloat("widthV",i));
			p.setProperty("dw",     clusBank.getFloat("widthW",i));
			p.setProperty("iu",    (clusBank.getInt("coordU", i)-4)/8+1);
			p.setProperty("iv",    (clusBank.getInt("coordV", i)-4)/8+1);
			p.setProperty("iw",    (clusBank.getInt("coordW", i)-4)/8+1);
			
			int ical = (int) p.getProperty("index");
			
			if(caliBank!=null) {
				p.setProperty("raweu", caliBank.getFloat("rawEU", ical)*1e3);
				p.setProperty("rawev", caliBank.getFloat("rawEV", ical)*1e3);
				p.setProperty("rawew", caliBank.getFloat("rawEW", ical)*1e3);
				p.setProperty("receu", caliBank.getFloat("recEU", ical)*1e3);
				p.setProperty("recev", caliBank.getFloat("recEV", ical)*1e3);
				p.setProperty("recew", caliBank.getFloat("recEW", ical)*1e3);
			}
			
			if(peakBank!=null && clusBank!=null) {
				int idu = clusBank.getInt("idU",ical)-1, idv = clusBank.getInt("idV",ical)-1, idw = clusBank.getInt("idW",ical)-1;
				p.setProperty("ustat", peakBank.getInt("status", idu));					
				p.setProperty("vstat", peakBank.getInt("status", idv));					
				p.setProperty("wstat", peakBank.getInt("status", idw));
				
				Point3D pc = new Point3D(p.getProperty("x"),p.getProperty("y"),p.getProperty("z"));
				p.setProperty("leffu", getLeff(pc,getPeakline(idu,pc,peakBank))); //readout distance U
				p.setProperty("leffv", getLeff(pc,getPeakline(idv,pc,peakBank))); //readout distance V
				p.setProperty("leffw", getLeff(pc,getPeakline(idw,pc,peakBank))); //readout distance W
			}
			
			ecalpart.add(p);
		}
		return ecalpart;
	}
	
	//returns List<Particle> for requested pindex "ipart" with ECAL relevant variables from REC::Calorimeter, REC::CaloExtras, REC::Traj, ECAL::clusters, ECAL::peaks
	public List<Particle> getECALPHYS(int ipart) { //for physics runs
		List<Particle> ecalpart = new ArrayList<Particle>();
		if(caloMap.containsKey(ipart)) {
			for(int imap : caloMap.get(ipart)) { //retrieve REC::Calorimeter index corresponding to REC::Particle index ipart			
				Particle p = new Particle(); 
				p.copy(part.get(        caloBank.getShort("pindex",imap)));
				p.setProperty("sector", caloBank.getByte("sector",imap)); 
				p.setProperty("layer",  caloBank.getByte("layer",imap)); //layer=1,4,7
				p.setProperty("pindex", caloBank.getShort("pindex",imap));
				p.setProperty("index",  caloBank.getShort("index",imap));
				p.setProperty("energy", caloBank.getFloat("energy",imap)*1e3);
				p.setProperty("time",   caloBank.getFloat("time",imap));
				p.setProperty("path",   caloBank.getFloat("path",imap));  
				p.setProperty("x",      caloBank.getFloat("x",imap));
				p.setProperty("y",      caloBank.getFloat("y",imap)); 
				p.setProperty("z",      caloBank.getFloat("z",imap)); 	
				p.setProperty("hx",     caloBank.getFloat("hx",imap));
				p.setProperty("hy",     caloBank.getFloat("hy",imap));
				p.setProperty("hz",     caloBank.getFloat("hz",imap));
				p.setProperty("lu",     caloBank.getFloat("lu",imap));
				p.setProperty("lv",     caloBank.getFloat("lv",imap));
				p.setProperty("lw",     caloBank.getFloat("lw",imap));
				p.setProperty("du",     caloBank.getFloat("du",imap));
				p.setProperty("dv",     caloBank.getFloat("dv",imap));
				p.setProperty("dw",     caloBank.getFloat("dw",imap));
				p.setProperty("m2u",    caloBank.getFloat("m2u",imap));
				p.setProperty("m2v",    caloBank.getFloat("m2v",imap));
				p.setProperty("m2w",    caloBank.getFloat("m2w",imap));
				p.setProperty("pid",    part.get((int)p.getProperty("pindex")).pid());
				
				int ical = (int) p.getProperty("index");				
				int  lay = (int) p.getProperty("layer"); 
				int  pid = (int) p.getProperty("pid");	
				
				float  bet = (float) part.get((int)p.getProperty("pindex")).getProperty("beta");
				float tim1 = (float) p.getProperty("time"), tim2=0;
				float  pat = (float) p.getProperty("path");
				
//				boolean goodMatch = clusBank!=null && is==clusBank.getByte("sector", ical) && lay==clusBank.getByte("layer",ical);
//				boolean goodMatch = clusBank!=null ;
				
//				if(!goodMatch) continue;
				
				if(clusBank==null && caleBank!=null) { //use REC::CaloExtras which is column mapped to REC::Calorimeter
					int  is = (int) p.getProperty("sector");
					p.setProperty("iu",   (caleBank.getInt("dbstU",   imap)/10));
					p.setProperty("iv",   (caleBank.getInt("dbstV",   imap)/10));
					p.setProperty("iw",   (caleBank.getInt("dbstW",   imap)/10));					
					p.setProperty("raweu", caleBank.getFloat("rawEU", imap)*1e3);
					p.setProperty("rawev", caleBank.getFloat("rawEV", imap)*1e3);
					p.setProperty("rawew", caleBank.getFloat("rawEW", imap)*1e3);
					p.setProperty("receu", caleBank.getFloat("recEU", imap)*1e3);
					p.setProperty("recev", caleBank.getFloat("recEV", imap)*1e3);
					p.setProperty("recew", caleBank.getFloat("recEW", imap)*1e3);	
					p.setProperty("leffu", getLeff(is,lay+0,(int)p.getProperty("iu"),p.getProperty("x"),p.getProperty("y"),p.getProperty("z")));//readout distance U
					p.setProperty("leffv", getLeff(is,lay+1,(int)p.getProperty("iv"),p.getProperty("x"),p.getProperty("y"),p.getProperty("z")));//readout distance V
					p.setProperty("leffw", getLeff(is,lay+2,(int)p.getProperty("iw"),p.getProperty("x"),p.getProperty("y"),p.getProperty("z")));//readout distance W
				}
				
				if(clusBank!=null) {
					p.setProperty("cstat",  clusBank.getInt("status", ical));					
					p.setProperty("iu",    (clusBank.getInt("coordU", ical)-4)/8+1);
					p.setProperty("iv",    (clusBank.getInt("coordV", ical)-4)/8+1);
					p.setProperty("iw",    (clusBank.getInt("coordW", ical)-4)/8+1);					
					p.setProperty("x",      clusBank.getFloat("x",ical)); //override caloBank with clusBank
					p.setProperty("y",      clusBank.getFloat("y",ical)); //override caloBank with clusBank
					p.setProperty("z",      clusBank.getFloat("z",ical)); //override caloBank with clusBank
					p.setProperty("energy", clusBank.getFloat("energy",ical)*1e3); //override caloBank with clusBank
					p.setProperty("time",   clusBank.getFloat("time",ical)); //override caloBank with clusBank
					tim2 = (float) p.getProperty("time");
				}
				
//				p.setVector(p.pid(),p.getProperty("x"),p.getProperty("y"),p.getProperty("z"),p.vx(),p.vy(),p.vz());	//why did i do this? lcs 10-9-24				
				
				if(clusBank!=null && caliBank!=null) {	//use ECAL::calib which is column mapped to ECAL::clusters				
					p.setProperty("raweu", caliBank.getFloat("rawEU", ical)*1e3);
					p.setProperty("rawev", caliBank.getFloat("rawEV", ical)*1e3);
					p.setProperty("rawew", caliBank.getFloat("rawEW", ical)*1e3);
					p.setProperty("receu", caliBank.getFloat("recEU", ical)*1e3);
					p.setProperty("recev", caliBank.getFloat("recEV", ical)*1e3);
					p.setProperty("recew", caliBank.getFloat("recEW", ical)*1e3);
				}
/*				
				if(true) {
					int  is = (int) p.getProperty("sector");					
					float leffu = getLeff(is,lay+0,(int)p.getProperty("iu"),p.getProperty("x"),p.getProperty("y"),p.getProperty("z"));
					float leffv = getLeff(is,lay+1,(int)p.getProperty("iv"),p.getProperty("x"),p.getProperty("y"),p.getProperty("z"));
					float leffw = getLeff(is,lay+2,(int)p.getProperty("iw"),p.getProperty("x"),p.getProperty("y"),p.getProperty("z"));
					System.out.println("sector "+is+" lay "+lay+" leff1: "+leffu+" "+leffv+" "+leffw);
					
					Point3D pc = new Point3D(p.getProperty("x"),p.getProperty("y"),p.getProperty("z"));
					int iidu = clusBank.getInt("idU",ical)-1, iidv = clusBank.getInt("idV",ical)-1, iidw = clusBank.getInt("idW",ical)-1;

					float lleffu = getLeff(pc,getPeakline(iidu,pc,peakBank));  
					float lleffv = getLeff(pc,getPeakline(iidv,pc,peakBank));  
					float lleffw = getLeff(pc,getPeakline(iidw,pc,peakBank));  
					System.out.println("sector "+is+" lay "+lay+" leff2: "+lleffu+" "+lleffv+" "+lleffw);					
				}
*/				
				if(peakBank!=null && clusBank!=null) {
					int idu = clusBank.getInt("idU",ical)-1, idv = clusBank.getInt("idV",ical)-1, idw = clusBank.getInt("idW",ical)-1;
					p.setProperty("ustat", peakBank.getInt("status", idu));					
					p.setProperty("vstat", peakBank.getInt("status", idv));					
					p.setProperty("wstat", peakBank.getInt("status", idw));
					
					Point3D pc = new Point3D(p.getProperty("x"),p.getProperty("y"),p.getProperty("z"));
					p.setProperty("leffu", getLeff(pc,getPeakline(idu,pc,peakBank))); //readout distance U
					p.setProperty("leffv", getLeff(pc,getPeakline(idv,pc,peakBank))); //readout distance V
					p.setProperty("leffw", getLeff(pc,getPeakline(idw,pc,peakBank))); //readout distance W
				}
				
				p.setProperty("beta", getBeta(p)); //override caloBank with clusBank if it exists
				
//				boolean good = debug && lay>0 && Math.abs(pid)==22;

//				if(good) System.out.println(eventNumber+" "+starttime+" "+pid+" "+lay+" "+tim1+" "+tim2+" "+pat+" "+bet+" "+(tim1-starttime)+" "+pat/bet/29.979f);

				if(trajMap.containsKey(ipart)) {
					int  is = (int) p.getProperty("sector");
					for(int tmap : trajMap.get(ipart)) {
					   if(trajBank.getInt("detector",tmap)==7 &&
					     (trajBank.getInt("layer",tmap)==p.getProperty("layer") ||
					      trajBank.getInt("layer",tmap)==p.getProperty("layer")+1) ) {
							p.setProperty("tx",     trajBank.getFloat("x",tmap));
							p.setProperty("ty",     trajBank.getFloat("y",tmap));
							p.setProperty("tz",     trajBank.getFloat("z",tmap));
							float cx = trajBank.getFloat("cx", tmap);
							float cy = trajBank.getFloat("cy", tmap);
							float cz = trajBank.getFloat("cz", tmap);
							Point3D xyz = new Point3D(cx,cy,cz);
							xyz.rotateZ(Math.toRadians(-60*(is-1)));
							xyz.rotateY(Math.toRadians(-25.));
							p.setProperty("cz", xyz.z());
					   }
					}
				}
				
				Point3D res = getResidual(p);
				p.setProperty("resx",res.x());p.setProperty("resy",res.y());p.setProperty("resz",res.z());
				
				ecalpart.add(p);
			}
		}
		return ecalpart;
	}
    
    public float getLeff(Point3D point, Line3D peakline) {
    	return (float) point.distance(peakline.end());
    }
	
    public float getLeff(int is, int il, int ic, double x, double y, double z) { // is: 1-6 il: 1,2,3=PCAL 4,5,6=ECIN 7,8,9=ECOU
        Line3D peakline = getPeakline(is,il,ic);
        Point3D po = new Point3D(peakline.origin().x(),peakline.origin().y(),peakline.origin().z());
        Point3D pe = new Point3D(peakline.end().x(),   peakline.end().y(),   peakline.end().z());
        return getLeff(new Point3D(x,y,z),new Line3D(po,pe));
    }
    
    public Line3D getPeakline(int is, int il, int ic) { //calculate peakline from geometry package
        int superlayer = (int) ((il-1)/3); //0=PCAL 1=ECIN 2=ECOU
        int localLayer = (il-1)%3;         //0=U 1=V 2=W
        int pcalz = DetectorLayer.PCAL_Z;
        int ecinz = DetectorLayer.EC_INNER_Z;
        int ecouz = DetectorLayer.EC_OUTER_Z;
        int off = superlayer==0 ? pcalz : (superlayer==1 ? ecinz : ecouz);
        Layer detLayer = ecDetector.getSector(is-1).getSuperlayer(superlayer).getLayer(localLayer+off); //localLayer+off=9,10,11 for U,V,W planes
        return ((ScintillatorPaddle) detLayer.getComponent(ic-1)).getLine();
    }
	    
    public Line3D getPeakline(int iid, Point3D point, DataBank bank) { //get peakline from ECAL::peaks
        Point3D  po = new Point3D(bank.getFloat("xo",iid),
                                  bank.getFloat("yo",iid),
                                  bank.getFloat("zo",iid));
        Point3D  pe = new Point3D(bank.getFloat("xe",iid),
                                  bank.getFloat("ye",iid),
                                  bank.getFloat("ze",iid));
        return  new Line3D(po,pe);
    }
	
	public float getVar(List<Particle> list, String property, int layer) {
		for (Particle p: list ) {
			if (p.getProperty("layer") == layer) {
				if(p.hasProperty(property)) return (float)p.getProperty(property);
			}
		}
		return 0f;
	}
	
	public Point3D getHxyz(Particle p) {
    	float hx = p.hasProperty("tx")?(float)p.getProperty("tx"):(float)p.getProperty("hx");
    	float hy = p.hasProperty("ty")?(float)p.getProperty("ty"):(float)p.getProperty("hy");
    	float hz = p.hasProperty("tz")?(float)p.getProperty("tz"):(float)p.getProperty("hz");
    	return new Point3D(hx,hy,hz);
	}
	
    public Point3D getResidual(Particle p) { 
    	Point3D hxyz = getHxyz(p);
        Point3D  xyz = new Point3D(hxyz.x()-(float)p.getProperty("x"),hxyz.y()-(float)p.getProperty("y"),hxyz.z()-(float)p.getProperty("z"));	
        return getRotTiltPoint(xyz,(int) p.getProperty("sector"));
    }
    
    public Point3D getRotTiltPoint(Point3D xyz, int is) {
        xyz.rotateZ(Math.toRadians(-60*(is-1)));
        xyz.rotateY(Math.toRadians(-25)); 	
        return xyz;    	
    }
/*	
	public void getRECparticleOLD(int tpid) {
//		if(tpid==11) System.out.println(" ");
		if (!partMap.containsKey(tpid)) return;
		for(int ipart : partMap.get(tpid)){  
			int      pid = partBank.getInt("pid",  ipart);              
			float     px = partBank.getFloat("px", ipart);
			float     py = partBank.getFloat("py", ipart);
			float     pz = partBank.getFloat("pz", ipart);
			float     vx = partBank.getFloat("vx", ipart);
			float     vy = partBank.getFloat("vy", ipart);
			float     vz = partBank.getFloat("vz", ipart);
			float   beta = partBank.getFloat("beta", ipart);
			short status = (short) Math.abs(partBank.getShort("status", ipart));
			
			if(ftofMap.containsKey(ipart)) {
			for(int imap : ftofMap.get(ipart)) {				
				Particle p = new Particle();                         
				p.initParticle(pid, px, py, pz, vx, vy, vz);                	   
				p.setProperty("ppid", pid);
				p.setProperty("status", status);
				p.setProperty("pindex",ipart); 
				p.setProperty("sector", ftofBank.getByte("sector",imap)); 
				p.setProperty("layer",  ftofBank.getByte("layer",imap));
				p.setProperty("index",  ftofBank.getShort("index",imap));
				p.setProperty("energy", ftofBank.getFloat("energy",imap));
				p.setProperty("time",   ftofBank.getFloat("time",imap));
				p.setProperty("path",   ftofBank.getFloat("path",imap));                   		                    
				p.setProperty("x",      ftofBank.getFloat("x",imap));
				p.setProperty("y",      ftofBank.getFloat("y",imap)); 
				p.setProperty("z",      ftofBank.getFloat("z",imap)); 
				p.setProperty("hx",     ftofBank.getFloat("hx",imap));
				p.setProperty("hy",     ftofBank.getFloat("hy",imap));
				p.setProperty("hz",     ftofBank.getFloat("hz",imap));
				
				int ip = pid<0?Math.abs(pid)+1:pid;				
				if(!part.hasItem(ip,1)) {part.add(new ArrayList<Particle>(),ip,1);}
				    part.getItem(ip,1).add(p);
			}
			}
			
			if(caloMap.containsKey(ipart)) {
			for(int imap : caloMap.get(ipart)) {				
				Particle p = new Particle();                         
				p.initParticle(pid, px, py, pz, vx, vy, vz);               	   
				p.setProperty("ppid",   pid);
				p.setProperty("status", status);
				p.setProperty("pindex", ipart); 
				p.setProperty("sector", caloBank.getByte("sector",imap)); 
				p.setProperty("layer",  caloBank.getByte("layer",imap));
				p.setProperty("index",  caloBank.getShort("index",imap));
				p.setProperty("energy", caloBank.getFloat("energy",imap)*1000);
				p.setProperty("time",   caloBank.getFloat("time",imap));
				p.setProperty("path",   caloBank.getFloat("path",imap));                   		                    
				p.setProperty("x",      caloBank.getFloat("x",imap));
				p.setProperty("y",      caloBank.getFloat("y",imap)); 
				p.setProperty("z",      caloBank.getFloat("z",imap)); 
				p.setProperty("hx",     caloBank.getFloat("hx",imap));
				p.setProperty("hy",     caloBank.getFloat("hy",imap));
				p.setProperty("hz",     caloBank.getFloat("hz",imap));				

				int ical = (int) p.getProperty("index");
				if(clusBank!=null) {
					p.setProperty("iu",    (clusBank.getInt("coordU", ical)-4)/8+1);
					p.setProperty("iv",    (clusBank.getInt("coordV", ical)-4)/8+1);
					p.setProperty("iw",    (clusBank.getInt("coordW", ical)-4)/8+1);
				}
				p.setProperty("beta", (pid==2112||pid==22)?newBeta(p):beta);
				
				if(trajMap.containsKey(ipart)) {
					for(int tmap : trajMap.get(ipart)) {
					   if(trajBank.getInt("detector",tmap)==7&&trajBank.getInt("layer",tmap)==(p.getProperty("layer")+1)) {
							p.setProperty("tx",     trajBank.getFloat("x",tmap));
							p.setProperty("ty",     trajBank.getFloat("y",tmap));
							p.setProperty("tz",     trajBank.getFloat("z",tmap));
							p.setProperty("cz",     trajBank.getFloat("cz",tmap));
//							System.out.println(p.getProperty("x")+" "+p.getProperty("y")+" "+p.getProperty("z"));
//							System.out.println(p.getProperty("hx")+" "+p.getProperty("hy")+" "+p.getProperty("hz"));
//							System.out.println(p.getProperty("tx")+" "+p.getProperty("ty")+" "+p.getProperty("tz"));
					   }
					}
				}				
				
				int ip = pid<0?Math.abs(pid)+1:pid;				
				if(!part.hasItem(ip,0)) {part.add(new ArrayList<Particle>(),ip,0);}
				    part.getItem(ip,0).add(p);  				
			}
			}			
		}		
	}	
*/	
    /**
     * @param fromBank the bank containing the index variable
     * @param idxVarName the name of the index variable
     * @return map with keys being the index in toBank and values the indices in fromBank
     */
    public static Map<Integer,List<Integer>> loadMapByIndex(DataBank fromBank, String idxVarName) {
        Map<Integer,List<Integer>> map=new HashMap<Integer,List<Integer>>();
        if (fromBank!=null) {
            for (int iFrom=0; iFrom<fromBank.rows(); iFrom++) {
                final int iTo = fromBank.getInt(idxVarName,iFrom);
                if (!map.containsKey(iTo)) map.put(iTo,new ArrayList<Integer>()); 
                map.get(iTo).add(iFrom);
            }
        }
        return map;
    }	
	
}
