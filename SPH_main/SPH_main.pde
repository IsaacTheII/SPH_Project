import peasy.*;

float a = 0.;
int leaf_size = 8;
int speed_up = 1;
int param = 0; // 2 = testing environment (low amount of particles), 1 = random, 0 = uniform
int iter = 27;
float e_ini = 1;
int nn = 16;
boolean dim = false;
float courant = 0.05;
int size = 1000;
float v_ini = 2;
int btype = 2;
float frame_dt = 0.;
float dt_thresh = 0.0025;
int total_framecount = 0;


Simulation sim;

PeasyCam cam;

void setup() {
  size(1200, 1200, P3D);
  colorMode(HSB, 1);
  smooth();

  cam = new PeasyCam(this, width/2, height/2, 0, 1000);


  sim = new Simulation(leaf_size, param, iter, e_ini, nn, dim, courant, size, v_ini, btype);
  sim.calc_forces();
}

void draw() {
  if (frameCount%100==0) {
    println(frameRate);
    //sim.reset_e();
  }
  background(0);

  translate(width/2, height/2);
  if (dim) {
    rotateY(a);
    a += 0.002;
  }

  // implement constant time per frame here instead of calculations per frame
  //int counter = 0;
  //while (frame_dt < dt_thresh) {
  //  sim.update();
  //  frame_dt += sim.dt;
  //  counter++;
  //  total_framecount++;
  //  if (total_framecount % 1000 == 0) sim.reset_e();
  //}
  //println(counter);
  sim.update();

  sim.show_particles();
  frame_dt = frame_dt%dt_thresh;
}
