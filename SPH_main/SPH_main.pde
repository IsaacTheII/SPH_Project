import peasy.*;

float a = 0.;
int leaf_size = 8;
int speed_up = 1;
int param = 1; // 1 = random, 0 = uniform
int iter = 21;
float e_ini = 1;
int nn = 16;
boolean dim = false;
float courant = 0.15;
int size = 1000;

Simulation sim;

PeasyCam cam;

void setup() {
  size(800, 800, P3D);
  colorMode(HSB, 1);
  smooth();
  
  cam = new PeasyCam(this, width/2, height/2,0,1000);


  sim = new Simulation(leaf_size, param, iter, e_ini, nn, dim, courant, size);
  sim.calc_forces();
}

void draw() {
  if(frameCount%10==0) println(frameRate);
  background(0);
  translate(width/2, height/2);
  if (dim) {
    rotateY(a);
    a += 0.002;
  }

  // implement constant time per frame here instead of calculations per frame
  for (int i = 0; i < speed_up; i++) {
    sim.update();
  }
  sim.show_particles();
}
