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

public class Event {
     
    public int run=0,event=0,unixtime=0;
    public long timestamp=0,trigger=0;
    public float starttime=0;
	
	private DataEvent ev = null;
	private DataBank caloBank = null;
	public int TRpid = 0;
	
	public IndexedList<List<Particle>> part = new IndexedList<List<Particle>>(2);
	Map<Integer,List<Integer>> caloMap = null;
	
	int[]        tb = new int[32];
	public int[] N1 = new int[6];
	public int[] N2 = new int[6];
	
	public Event() {
		
	}
	
	public void init() {		
		trigger=0;
		starttime = -100;	
		part.clear();
	}
	
	public void procEvent(DataEvent event) {
		this.ev = event;
        init(); 
	    if(ev.hasBank("RUN::config"))        processRUNconfig();
	    if(ev.hasBank("REC::Event"))         processRECevent();
	    if(ev.hasBank("REC::Calorimeter"))   processRECcalorimeter();
	    if(ev.hasBank("REC::Particle"))      processRECparticle();
	    if(ev.hasBank("REC::Scintillator"))  processRECscintillator();
	    if(ev.hasBank("REC::Track"))         processRECtrack();
	    if(ev.hasBank("ECAL::clusters"))     processECALclusters();
	    if(ev.hasBank("ECAL::calib"))        processECALcalib();		
	}
	
	public List<Particle> getParticle(int ipid) {		
		List<Particle> pout = new ArrayList<Particle>();
	    IndexGenerator ig = new IndexGenerator();
	    for (Map.Entry<Long,List<Particle>>  entry : part.getMap().entrySet()){
	           long hash = entry.getKey();
	           int pid = ig.getIndex(hash, 0);
	           int sec = ig.getIndex(hash, 1);
	           if(ipid==pid) {for (Particle pp : part.getItem(pid,sec)) pout.add(pp);}
	    }	    
	    return pout;
	}
	
	public void setTrigger(int trig) {
		TRpid = trig;
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
		this.starttime = bank.getFloat("startTime",0);
	}
	
	public void storeRECparticle(DataBank bank) {
        if(starttime > -100) {        	
            for(int i = 0; i < bank.rows(); i++){           	
                int      pid = bank.getInt("pid", i);              
                float     px = bank.getFloat("px", i);
                float     py = bank.getFloat("py", i);
                float     pz = bank.getFloat("pz", i);
                float     vx = bank.getFloat("vx", i);
                float     vy = bank.getFloat("vy", i);
                float     vz = bank.getFloat("vz", i);
                float   beta = bank.getFloat("beta", i);
                short status = (short) Math.abs(bank.getShort("status", i));
                for (int ic : caloMap.get(i)) {
                	if (pid!=0) {
                      Particle p = new Particle();                         
                	  p.initParticle(pid, px, py, pz, vx, vy, vz);                	   
                	  p.setProperty("ppid", pid);
                      p.setProperty("status", status);
                      p.setProperty("pindex",i); 
         	          p.setProperty("sector", caloBank.getByte("sector",ic)); 
                      p.setProperty("layer",  caloBank.getByte("layer",ic));
                      p.setProperty("index",  caloBank.getShort("index",ic));
                      p.setProperty("energy", caloBank.getFloat("energy",ic)*1000);
         	          p.setProperty("time",   caloBank.getFloat("time",ic));
         	          p.setProperty("path",   caloBank.getFloat("path",ic));
         	          p.setProperty("x",      caloBank.getFloat("x",ic));
         	          p.setProperty("y",      caloBank.getFloat("y",ic));
         	          p.setProperty("z",      caloBank.getFloat("z",ic));
                      p.setProperty("beta", (pid==2112)?newBeta(p):beta);
         	          int ip = pid<0?Math.abs(pid)+1:pid;
         	          int is = caloBank.getByte("sector",ic); 
         	          if(!part.hasItem(ip,is)) part.add(new ArrayList<Particle>(),ip,is);
         	              part.getItem(ip,is).add(p);
                	}
                }                         
            }            
        }				
	}
	
	public float newBeta(Particle p) {
		double path= p.getProperty("path");
		double time= p.getProperty("time");
		return (float) (path/(time-starttime-2.0f)/29.97f);
	}
	
	public void storeRECscintillator(DataBank bank) {	
		
	}
	
	public void storeRECcalorimeter(DataBank bank) {
		caloBank = bank;
		caloMap = loadMapByIndex(caloBank,"pindex");
	}	
	
	public void storeRECtrack(DataBank bank) {	
		
	}
	
	public void storeECALclusters(DataBank bank) {	
		
	}
	
	public void storeECALcalib(DataBank bank) {	
		
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
			System.out.print("Run "+run+" Event "+event+" ");
			System.out.println("Bits "+tb[0]+" "+tb[1]+tb[2]+tb[3]+tb[4]+tb[5]+tb[6]+" "+tbsum);
			if(tbsum==1) {for (int i=1; i<7; i++) {if(tb[i]==1) N1[i-1]++;}}
			if(tbsum>1)  {for (int i=1; i<7; i++) {if(tb[i]==1) N2[i-1]++;}}
			System.out.println(N1[0]+" "+N1[1]+" "+N1[2]+" "+N1[3]+" "+N1[4]+" "+N1[5]);
			System.out.println(N2[0]+" "+N2[1]+" "+N2[2]+" "+N2[3]+" "+N2[4]+" "+N2[5]);
		}
		if(tbsum==0) return false;
		return true;			
	}
	
    /**
     * @param fromBank the bank containing the index variable
     * @param idxVarName the name of the index variable
     * @return map with keys being the index in toBank and values the indices in fromBank
     */
     public static Map<Integer,List<Integer>> loadMapByIndex( 
             DataBank fromBank,
             String idxVarName) {
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
