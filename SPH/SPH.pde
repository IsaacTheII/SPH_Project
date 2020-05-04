//import Particle;

Particle p;
Cell root;
ArrayList<Particle> particles;


void setup() {
  p = new Particle(new PVector(0, 0, 0), new PVector(0.001, 0, 0), 1, 1, 0, 0);
  p.periodic_boundary();
  root = new Cell(new PVector(0,0,0), new PVector(1,1,1), 8, null);
  particles = new ArrayList<Particle>();
  particles.add(p);
  particles.add(null);
  particles.add(p);
  particles.add(p);
  particles.add(p);
  particles.add(p);
}

void draw(){
  p.update_pos(10);
  //particles.add(p);
  for(Particle p : particles){
    println(p.r);
  }
  //p.periodic_boundary();
}
