package org.clas.tools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jlab.clas.physics.Particle;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedList.IndexGenerator;
import org.jlab.geom.prim.Point3D;

public class Event {
     
    public int run=0,event=0,unixtime=0;
    public long timestamp=0,trigger=0;
    public float starttime=0;
	
	private DataEvent      ev = null;
	private DataBank partBank = null;
	private DataBank caloBank = null;
	private DataBank clusBank = null;
	public int TRpid = 0;
	private boolean isHipo3Event; 
	
	public IndexedList<List<Particle>> part = new IndexedList<List<Particle>>(1);
	Map<Integer,List<Integer>>      caloMap = new HashMap<Integer,List<Integer>>();
	Map<Integer,List<Integer>>      partMap = new HashMap<Integer,List<Integer>>();
	
	private boolean hasRUNconfig = false;
	private boolean hasRECevent  = false;
	private boolean hasRECscintillator = false;
	private boolean hasRECcalorimeter = false;
	private boolean hasECALclusters = false;
	private boolean hasRECparticle = false;
	private boolean hasRECtrack = false;
	private boolean hasECALcalib = false;
	
	int[]        tb = new int[32];
	public int[] N1 = new int[6];
	public int[] N2 = new int[6];
	
	private int trigger_sect=0;
	private int eventNumber=0;
	private int nelec=0;
	
	public Event() {
		
	}
	
	public void init() {	 
		trigger=0;
		starttime = -100;
		part.clear();
		partMap.clear();
		caloMap.clear();
		partBank = null;
		caloBank = null;
		clusBank = null;
	}
	
	public boolean filter() {
		hasRUNconfig       = ev.hasBank("RUN::config");
		hasRECevent        = ev.hasBank("REC::Event");
		hasRECscintillator = ev.hasBank("REC::Scintillator");
		hasRECcalorimeter  = ev.hasBank("REC::Calorimeter");
		hasECALclusters    = ev.hasBank("ECAL::clusters");
		hasRECparticle     = ev.hasBank("REC::Particle");
		hasRECtrack        = ev.hasBank("REC::Track");
		hasECALcalib       = ev.hasBank("ECAL::calib");	
//		return hasRUNconfig&&hasRECevent&&hasRECparticle&&hasRECcalorimeter&&hasECALclusters;
		return hasRUNconfig&&hasRECevent&&hasRECparticle&&hasRECcalorimeter;
	}
	
	public boolean procEvent(DataEvent event) {
		this.ev = event;
        init(); 
        if(!filter()) return false;
	    if(hasRUNconfig)      processRUNconfig();
	    if(!countElectronTriggers(false)) return false;
	    if(hasRECevent)       processRECevent();
	    if(hasRECcalorimeter) processRECcalorimeter();
	    if(hasECALclusters)   processECALclusters();
	    if(hasRECcalorimeter) processRECscintillator();
	    if(hasRECtrack)       processRECtrack();
	    if(hasECALcalib)      processECALcalib();		
	    if(hasRECparticle)    processRECparticle();
		if(starttime > -100) {
			getRECparticle(11);
			getRECparticle(22);
			getRECparticle(2112);
			getRECparticle(211);
			getRECparticle(-211);
		}
	    return true;
	}
	
	public List<Particle> getParticle(int ipid) {		
		List<Particle> pout = new ArrayList<Particle>();
	    IndexGenerator ig = new IndexGenerator();                
	    for (Map.Entry<Long,List<Particle>>  entry : part.getMap().entrySet()){
	           int pid = ig.getIndex(entry.getKey(), 0);   
	           if(ipid==pid) {for (Particle pp : entry.getValue()) pout.add(pp);} 
	    }	
	    return pout;
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
	
	public void processRUNconfig() {		
		storeRUNconfig(ev.getBank("RUN::config"));
	}
	
	public void processRECevent() {		
		storeRECevent(ev.getBank("REC::Event"));  	
	}	
	
	public void processRECparticle() {		
	 	storeRECparticle(ev.getBank("REC::Particle"));	 
	}
	
	public void processRECscintillator() {		
		storeRECscintillator(ev.getBank("REC::Scintillator"));				
	}	
	
	public void processRECcalorimeter() {		
		storeRECcalorimeter(ev.getBank("REC::Calorimeter"));  		
	}
	
	public void processRECtrack() {		
		storeRECtrack(ev.getBank("REC::Track"));		
	}
	
	public void processECALclusters() {
		storeECALclusters(ev.getBank("ECAL::clusters"));		
	}	
	
	public void processECALcalib() {		
		storeECALcalib(ev.getBank("ECAL::calib"));		
	}
	
	public void storeRUNconfig(DataBank bank) {
		this.run       = bank.getInt("run",0);
		this.event     = bank.getInt("event",0);
		this.unixtime  = bank.getInt("unixtime",0);
		this.trigger   = bank.getLong("trigger",0);
		this.timestamp = bank.getLong("timestamp",0);		
	}
	
	public void storeRECevent(DataBank bank) {
		this.starttime = isHipo3Event ? bank.getFloat("STTime", 0):
                                        bank.getFloat("startTime", 0);		
	}
	
	public void storeRECparticle(DataBank bank) {
		partBank = bank;
		partMap  = loadMapByIndex(partBank,"pid");	
	}

	
	public void storeRECcalorimeter(DataBank bank) {
		caloBank = bank;
		caloMap  = loadMapByIndex(caloBank,"pindex");
	}
	
	public void storeRECscintillator(DataBank bank) {	
		
	}	
	
	public void storeRECtrack(DataBank bank) {	
		
	}
	
	public void storeECALclusters(DataBank bank) {	
		clusBank = bank;
	}
	
	public void storeECALcalib(DataBank bank) {	
		
	}		
	
	public void reportElectrons(String tag) {
		nelec++;System.out.println("Evnt "+eventNumber+" Nelec "+nelec+" "+tag);
	}
	
	public float newBeta(Particle p) {
		double path= p.getProperty("path");
		double time= p.getProperty("time");
		return (float) (path/(time-starttime)/29.97f);
	}	
	
    public int     getFDTrigger()            {return (int)(trigger)&0x000000000ffffffff;}
    public int     getCDTrigger()            {return (int)(trigger>>32)&0x00000000ffffffff;}
    public boolean isGoodFD()                {return  getFDTrigger()>0;}    
    public boolean isTrigBitSet(int bit)     {int mask=0; mask |= 1<<bit; return isTrigMaskSet(mask);}
    public boolean isTrigMaskSet(int mask)   {return (getFDTrigger()&mask)!=0;}	
	
	public boolean countElectronTriggers(boolean debug) {
		int tbsum=0;
		for (int i = 31; i >= 0; i--) {tb[i] = ((trigger & (1 << i))!=0)?1:0;}
//		for (int i = 31; i >= 0; i--) {tb[i] = (isTrigBitSet(i))?1:0;}	
		for (int i=1; i<7; i++) tbsum+=tb[i];
		if(debug) {
			System.out.print("Run "+run+" Event "+event+" Trigger Sector "+trigger_sect);
			System.out.println("Bits "+tb[0]+" "+tb[1]+tb[2]+tb[3]+tb[4]+tb[5]+tb[6]+" "+tbsum);
			if(tbsum==1) {for (int i=1; i<7; i++) {if(tb[i]==1) N1[i-1]++;}}
			if(tbsum>1)  {for (int i=1; i<7; i++) {if(tb[i]==1) N2[i-1]++;}}
			System.out.println(N1[0]+" "+N1[1]+" "+N1[2]+" "+N1[3]+" "+N1[4]+" "+N1[5]);
			System.out.println(N2[0]+" "+N2[1]+" "+N2[2]+" "+N2[3]+" "+N2[4]+" "+N2[5]);
		}
		if(tbsum==0||tbsum>1) return false;
		return true;			
	}
	
	
	public void getRECparticle(int tpid) {
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
			if(caloMap.containsKey(ipart)) {
			for(int icalo : caloMap.get(ipart)) {				
				Particle p = new Particle();                         
				p.initParticle(pid, px, py, pz, vx, vy, vz);                	   
				p.setProperty("ppid", pid);
				p.setProperty("status", status);
				p.setProperty("pindex",ipart); 
				p.setProperty("sector", caloBank.getByte("sector",icalo)); 
				p.setProperty("layer",  caloBank.getByte("layer",icalo));
				p.setProperty("index",  caloBank.getShort("index",icalo));
				p.setProperty("energy", caloBank.getFloat("energy",icalo)*1000);
				p.setProperty("time",   caloBank.getFloat("time",icalo));
				p.setProperty("path",   caloBank.getFloat("path",icalo));
				Point3D xyz = new Point3D(caloBank.getFloat("x",icalo),
										  caloBank.getFloat("y",icalo),
                    		              caloBank.getFloat("z",icalo));                    		                    
				p.setProperty("x",      xyz.x());
				p.setProperty("y",      xyz.y());
				p.setProperty("z",      xyz.z());
				p.setProperty("hx",     caloBank.getFloat("hx",icalo));
				p.setProperty("hy",     caloBank.getFloat("hy",icalo));
				p.setProperty("hz",     caloBank.getFloat("hz",icalo));
				int ical = (int) p.getProperty("index");
				if(clusBank!=null) {
				p.setProperty("iu",    (clusBank.getInt("coordU", ical)-4)/8+1);
				p.setProperty("iv",    (clusBank.getInt("coordV", ical)-4)/8+1);
				p.setProperty("iw",    (clusBank.getInt("coordW", ical)-4)/8+1);
				}
				p.setProperty("beta", (pid==2112||pid==22)?newBeta(p):beta);
				int ip = pid<0?Math.abs(pid)+1:pid;
				
				if(!part.hasItem(ip)) {part.add(new ArrayList<Particle>(),ip);}
				    part.getItem(ip).add(p);  
				
			}
			}			
		}		
	}	
	
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