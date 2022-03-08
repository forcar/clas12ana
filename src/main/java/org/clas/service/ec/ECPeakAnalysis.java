package org.clas.service.ec;

import java.util.List;

/**
 *
 * @author gavalian
 * modified lcsmith 2/22/22
 */

public class ECPeakAnalysis {
	
	static int splitStrip;
    
    public static int getPeakSplitIndex(List<ECPeak> peaks){        
        for(int i = 0; i < peaks.size(); i++){
            splitStrip = peaks.get(i).getSplitStrip(); //index of strip used to split peak
            if(splitStrip>=0) return i; //index of peak tagged to be split
        }
        return -1;
    }
    
    public static void splitPeaks(List<ECPeak> peaks){
        
        while(true){ //repeat processing all peaks until no split found
            int index = getPeakSplitIndex(peaks);
        	if(ECCommon.debugSplit) System.out.println("New Iteration "+splitStrip);
            if(index<0){
                return; // no split was found in any peak.  Exit.
            } else {
                ECPeak  peak = peaks.get(index); //retrieve tagged peak with split candidate
                peaks.remove(index); //tagged peak removed from list 
                peaks.addAll(peak.splitPeak(splitStrip)); //two split peaks returned to list
            }
        }
    }
    
}
