float a = 0.;
int leaf_size = 8;
int speed_up = 1;
int param = 1; // 1 = random, 0 = uniform
int iter = 50;
float e_ini = 10000000.;
int nn = 32;
boolean dim = false;
float courant = 0.05;
int size = 700;

Simulation sim;

void setup() {
  size(1920, 1080, P3D);
  colorMode(HSB, 1);
  smooth();


  sim = new Simulation(leaf_size, param, iter, e_ini, nn, dim, courant, size);
  sim.calc_forces();
}

void draw() {

  background(0);
  translate(width/2, height/2);
  if (dim) {
    rotateY(a);
    a += 0.002;
  }
  //if (frameCount % 50 == 0) {
  //  sim.reset_e();
  //}

  // implement constant time per frame here instead of calculations per frame
  for (int i = 0; i < speed_up; i++) {
    sim.update();
  }
  println(frameRate);
  sim.show_particles();
}
