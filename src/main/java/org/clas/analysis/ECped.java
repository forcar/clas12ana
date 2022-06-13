package org.clas.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.swing.filechooser.FileSystemView;

import org.clas.viewer.DetectorMonitor;
import org.jlab.utils.groups.IndexedList;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.utils.groups.IndexedList.IndexGenerator;

public class ECped extends DetectorMonitor {
	
	String path;
	
	int crate,slot;
	int[] ec_crate = {1,7,13,19,25,31};
    int[] pc_crate = {3,9,15,21,27,33};
    int[] ecal_crate = {1,3,7,9,13,15,19,21,25,27,31,33};
    int[] slotmax = {19,17}; int slm;
    
	IndexedList<Float>  ped = new IndexedList<Float>(3);
	IndexedTable dbped;
		
    public ECped(String name) {
    	super(name);		
    	path = FileSystemView.getFileSystemView().getHomeDirectory().toString()+"/ECMON/PED/";
    }
    
    public void initCCDB(int runno) {
        dbped  = cm.getConstants(runno, "/daq/fadc/ec");
    }
    
    public List<String> getList(String path) { 
    	System.out.println(path);
    	List<String> out = new ArrayList();
    	String strCurrentLine;
       
    	try (BufferedReader br = new BufferedReader(new FileReader(path))) {
    	while ((strCurrentLine = br.readLine()) != null) out.add(strCurrentLine);
    	} catch (IOException e) {
    		e.printStackTrace();    		
    	}    	
    	return out;
    }
    
    public void parseList(List<String> list) {    	    	    	
    	for (String item : list) {
    		String parts1[] = item.split("\\s+"); 
    		if(parts1[0].equals("FADC250_SLOT")) slot = Integer.parseInt(parts1[1]);
    		if(slot<slm && parts1[0].equals("FADC250_ALLCH_PED")) for (int i=0; i<16; i++) {ped.add(Float.parseFloat(parts1[i+1]),crate,slot,i);}
    	}     	
    }
    
    public void processFile() {
		for (int detid=0; detid<2; detid++) { slm=slotmax[detid];
			for (int is=1; is<7; is++) {
				crate = detid==0 ? ec_crate[is-1]:pc_crate[is-1];
				parseList(getList(path+"adc"+(detid==0?"ecal":"pcal")+is+"_ped.cnf"));
			}
		}
    }
    
	public void writeNewTable(int run) {
		
		String line = new String();
		
		try { 
			File outputFile = new File(path+"fadc_"+run);
			FileWriter outputFw = new FileWriter(outputFile.getAbsoluteFile());
			BufferedWriter outputBw = new BufferedWriter(outputFw);
			
			System.out.println("ECped.writefile(fadc_"+run+")");
			
			IndexGenerator ig = new IndexGenerator();
			for (Map.Entry<Long,Float>  entry : ped.getMap().entrySet()){
				long hash = entry.getKey();
				int cr = ig.getIndex(hash, 0);
				int sl = ig.getIndex(hash, 1);
				int ch = ig.getIndex(hash, 2);
				line = cr+" "+sl+" "+ch+" "+
				ped.getItem(cr,sl,ch)+" "+
				dbped.getIntValue("nsb",cr,sl,ch)+" "+		
				dbped.getIntValue("nsa",cr,sl,ch)+" "+		
				dbped.getIntValue("tet",cr,sl,ch)+" "+	
				dbped.getIntValue("window_offset",cr,sl,ch)+" "+		
				dbped.getIntValue("window_size",cr,sl,ch);	

				if (line!=null) {
				      System.out.println(line);
				      outputBw.write(line);
				      outputBw.newLine();
				}
			}			
			outputBw.close();
			outputFw.close();			
		}
		
		catch(IOException ex) {
			System.out.println("Error writing file '" );                   
			ex.printStackTrace();
		}
		
	}
    
	public static void main(String[] args) {		
		ECped reader = new ECped("ECped");
		int run = 15928;	 	
		if(args.length==0) reader.initCCDB(run); 
		if(args.length==1) reader.initCCDB(Integer.parseInt(args[0])); 
	    reader.processFile();
	    reader.writeNewTable(run);
	}    
}
