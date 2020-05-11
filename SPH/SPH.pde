float a = 0;
float b = 0;
float c = 0;


Node root;
int leaf_size = 8;
int speed_up = 1;
int param = 1; // 1 = random, 0 = uniform
int iter = 21;
float e_ini = 100;
int nn = 32;
boolean dim = false;
float courant = 0.05;
int size = 700;

Simulation sim;

void setup() {
  size(1920, 1080, P3D);
  sim = new Simulation(leaf_size, speed_up, param, iter, e_ini, nn, dim, courant, size);
  colorMode(HSB, 1);

  sim.calc_forces();
}

void draw() {

  background(0);
  translate(width/2, height/2);
  //sim.nn_density(sim.root);



  //rotateX(a);
  //rotateY(b);
  //rotateZ(c);
  //a += 0.001;
  //b += 0.005;
  //c += 0.003;

  for (int i = 0; i < speed_up; i++) {
    sim.update();
  }
  sim.show_particles();
}
