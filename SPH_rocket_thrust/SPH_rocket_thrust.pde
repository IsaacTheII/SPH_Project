import peasy.*;

float a = 0.;
int leaf_size = 8;
int speed_up = 1;
int iter = 27;
float e_ini = 1;
int nn = 8;
boolean dim = false;
float courant = 0.05;
int size = 400;
float v_ini = 10;
int btype = 1;
int chambertype = 1;
int atmospheretype = 1;
float frame_dt = 0.;
float dt_thresh = 0.0006;
int total_framecount = 0;
int w = 4;
int h = 1;

Simulation sim;

void setup() {
  size(1600, 800, P2D);
  colorMode(HSB, 1);
  smooth();

  sim = new Simulation(leaf_size, iter, e_ini, nn, dim, courant, size, v_ini, btype, atmospheretype, chambertype, w, h);
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

   //implement constant time per frame here instead of calculations per frame
  int counter = 0;
  while (frame_dt < dt_thresh) {
    sim.update();
    frame_dt += sim.dt;
    counter++;
    total_framecount++;
    //if (total_framecount % 1000 == 0) sim.reset_e();
  }
  //println(counter);
  //total_framecount++;
  //if (total_framecount % 100 == 0) sim.reset_e();
  //sim.update();

  sim.show_particles();
  frame_dt = frame_dt%dt_thresh;
  saveFrame("movie/SPH_####.jpg"); //<>//
}
