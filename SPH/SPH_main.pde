float a = 0;
float b = 0;
float c = 0;


Node root;
int leaf_size = 8;
int speed_up = 1;
int num_steps = 1;
int param = 1; // 1 = random, 0 = uniform
int iter = 10;
float e_ini = 100;
int nn = 32;
boolean dim = true;
float courant = 0.1;
int size = 700;

Simulation sim;

void setup() {
  size(1920, 1080, P3D);
  sim = new Simulation(leaf_size, speed_up, num_steps, param, iter, e_ini, nn, dim, courant, size);
  sim.treeBuild(sim.root, 0);
  sim.nn_density(sim.root);
  colorMode(HSB, 1);

  //sim.calc_forces()
  //noSmooth();
}

void draw() {
  background(0);
  translate(width/2, height/2);

  if (dim) {
    rotateX(a);
    rotateY(b);
    rotateZ(c);
    a += 0.001;
    b += 0.005;
    c += 0.003;
  }

  //sim.update()
  sim.show_particles();
  println(frameRate);
}
