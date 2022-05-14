package org.clas.tools;

import javax.swing.JFrame;

import org.jlab.clas.pdg.PDGDatabase;
import org.jlab.clas.pdg.PDGParticle;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.ParticleGenerator;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.clas.physics.Vector3;
import org.jlab.clas.reactions.IDecay;
import org.jlab.groot.data.H2F;
import org.jlab.groot.graphics.EmbeddedCanvas;

/**
 *
 * @author gavalian (modifications by lcsmith)
 */
public class TwoBodyDecay {
    
    int decayParticleID1;
    int decayParticleID2;
    int parentParticleID;
    
    LorentzVector decayProd1;
    LorentzVector decayProd2;
    
    float cos_min=-1f, cos_max=+1f;
    
    public TwoBodyDecay() {        
        this.setDecayParticle(111);
        this.setDecayProducts(22,22);
    }
    
    public TwoBodyDecay(int parentID, int childid1, int childid2) {
        this.setDecayParticle(parentID);
        this.setDecayProducts(childid1,childid2);
    }

    public void setDecayParticle(int id) {
        this.parentParticleID = id;
    }

    public void setDecayProducts(int pid1, int pid2) {
        decayProd1 = new LorentzVector();
        decayProd2 = new LorentzVector();
        this.decayParticleID1 = pid1;
        this.decayParticleID2 = pid2;
    }
    
    public void setCosRange(float min, float max) {
    	cos_min = min;
    	cos_max = max;
    }
    
    private double getPhiRandom() {
        return Math.random()*2.0*Math.PI-Math.PI;
    }
    
    public float getCosThetaRandom(){
        return (float) (cos_min + Math.random()*(cos_max-cos_min));
    }
    
    public void decayParticles(PhysicsEvent event) {
        double cosTheta = this.getCosThetaRandom();
        double phi      = this.getPhiRandom();

        PDGParticle p1 = PDGDatabase.getParticleById(decayParticleID1);
        PDGParticle p2 = PDGDatabase.getParticleById(decayParticleID2);
    
        Particle mother = event.getParticleByPid(parentParticleID, 0);
        LorentzVector  vector = new LorentzVector();
        vector.copy(mother.vector());
        LorentzVector[] vec = DecayKinematics.getDecayParticlesLab(vector,
                p1.mass(), p2.mass(), Math.acos(cosTheta), phi);
        decayProd1.copy(vec[0]);
        decayProd2.copy(vec[1]);
        
        int index = event.getParticleIndex(parentParticleID, 0);
        event.removeParticle(index);
        
        event.addGeneratedParticle(new Particle(decayParticleID1,
                vec[0].px(),  vec[0].py(), vec[0].pz(),
                mother.vertex().x(),mother.vertex().y(),mother.vertex().z()
        ));
        
        event.addGeneratedParticle(new Particle(decayParticleID2,
                vec[1].px(),  vec[1].py(), vec[1].pz(),
                mother.vertex().x(),mother.vertex().y(),mother.vertex().z()
        ));                
    }
    
    public void pi0Demo(Particle p1, Particle p2) {
    	
    }
    
    public static void main(String[] args) {
    	
        JFrame frame = new JFrame("Pizero");
        frame.setSize(800,800);    	
        
        EmbeddedCanvas c = new EmbeddedCanvas(); 
    	
    	H2F h = new H2F("decay",50, -10, 10, 50,-10,10);
    	
    	TwoBodyDecay decay = new TwoBodyDecay();
    	decay.setCosRange(-1,1); 
    	
    	ParticleGenerator pg = new ParticleGenerator(111);    	
    	pg.setRange(1,8,24,26,55,65); 
    	
    	for (int i=0; i<10000; i++) {
        	PhysicsEvent event = new PhysicsEvent(); event.clear();
        	event.addParticle(pg.getParticle());
        	decay.decayParticles(event);
    	
        	Particle p1 = event.getGeneratedParticle(0); 
        	Particle p2 = event.getGeneratedParticle(1);
    	
        	float copa = (float) p1.cosTheta(p2);    
        	float  opa = (float) Math.toDegrees(Math.acos(copa));
        	float    x = (float) Math.abs(((p1.e()-p2.e())/(p1.e()+p2.e())));
        	float  ivm = (float) Math.sqrt(2*p1.e()*p2.e()*(1-copa));
    	
        	float dthe = (float)(Math.toDegrees(p1.theta())-Math.toDegrees(p2.theta()));
        	float dphi = (float)(Math.toDegrees(p1.phi())  -Math.toDegrees(p2.phi()));
    	
        	h.fill(dphi, dthe); 
    	}
 	
    	c.cd(0); c.draw(h);     
        frame.add(c);
        frame.setLocationRelativeTo(null);
        frame.setVisible(true);
    	
//    	System.out.println("OPA = "+opa+" X = "+x+" IVM= "+ivm);
//    	System.out.println(event.toLundStringGenerated());    	
    }

}
