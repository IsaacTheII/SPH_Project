class Simulation { //<>// //<>// //<>// //<>// //<>//
  int leaf_size;
  int num_particles;
  int iter;
  float e_ini;
  int nn;
  boolean dim;
  float courant;
  float gamma;
  float sigma;
  Node root;
  ArrayList<Particle> particles;
  int size;
  float t_global = 0;
  float dt = 0;
  float max_val = 10;
  float v_ini;
  float rotationSpeed = 0.002;
  int btype;
  ArrayList<Boundary> boundaries;
  PVector rlow;
  PVector rhigh;
  ArrayList<Particle> garbage_colleciton = new ArrayList<Particle>();

  Simulation(int leaf_size_, int param_, int iter_, float e_ini_, int nn_, boolean dim_, float courant_, int size_, float v_ini_, int btype_) {
    leaf_size = leaf_size_;
    iter = iter_;
    e_ini = e_ini_;
    nn = nn_;
    dim = dim_;
    courant = courant_;
    size = size_;
    v_ini = v_ini_;
    btype = btype_;
    boundaries = new ArrayList<Boundary>();

    if (dim) {
      gamma = 5./3.; // gamma = (f+2)/f , f = number of degrees of freedom
      sigma = 8./PI;
      num_particles = int(pow(iter, 3));
    } else {
      gamma = 2.;
      sigma = 40./(7 * PI);
      num_particles = int(pow(iter, 2));
    }
    particles = new ArrayList<Particle>();
    rlow = getRlow();
    rhigh = getRhigh();
    root = new Node(0, 10, rlow, rhigh, dim); // hardcoded 10 particles to start with btype = 5
    read_data(param_);
    createBoundaries();
  }

  PVector getRlow() {
    if (dim) {
      return new PVector(0, 0, 0);
    } else {
      return new PVector(0, 0);
    }
  }

  PVector getRhigh() {
    if (dim) {
      return new PVector(1, 1, 1);
    } else {
      return new PVector(1, 1);
    }
  }


  void read_data(int param) {
    println("Creating ", num_particles, " particles");

    randomSeed(0);
    if (dim) {
      if (param == 1) {
        for (int i = 0; i < num_particles-1; i++) {
          float x = random(1);
          float y = random(1);
          float z = random(1);
          Particle particle = new Particle(new PVector(x, y, z), new PVector(0, 0, 0), new PVector(0, 0, 0), 1./num_particles, 1.);
          particles.add(particle);
        }
        Particle particle = new Particle(new PVector(0.5, 0.5, 0.5), new PVector(0, 0, 0), new PVector(0, 0, 0), 1./num_particles, e_ini);
        particles.add(particle);
      } else {
        float spacing = 1. / iter;
        for (int x = 0; x < iter; x++) {
          for (int y = 0; y < iter; y++) {
            for (int z = 0; z < iter; z++) {
              if ((x == floor(iter/2)) && (y == floor(iter/2)) && (z == floor(iter/2))) {
                Particle particle = new Particle(new PVector(spacing * x + 0.5 * spacing, spacing * y + 0.5 * spacing, spacing * z + 0.5 * spacing), new PVector(0, 0, 0), new PVector(0, 0, 0), 1./num_particles, e_ini);
                particles.add(particle);
              } else {
                PVector pos = new PVector(spacing * x + 0.5 * spacing + randomGaussian() / 100, spacing * y + 0.5 * spacing + randomGaussian() / 100, spacing * z + 0.5 * spacing + randomGaussian() / 100);
                Particle particle = new Particle(pos, new PVector(0, 0, 0), new PVector(0, 0, 0), 1./num_particles, 1.);
                particles.add(particle);
              }
            }
          }
        }
      }
    } else {
      if (param == 1 && btype < 5) {
        for (int i = 0; i < num_particles-1; i++) {
          float x = random(1);
          float y = random(1);
          Particle particle = new Particle(new PVector(x, y), new PVector(v_ini, 0), new PVector(0, 0), 1./num_particles, 1);
          particles.add(particle);
        }
        Particle particle = new Particle(new PVector(0.5, 0.5), new PVector(v_ini, 0), new PVector(0, 0), 1./num_particles, e_ini);
        particles.add(particle);
      } else if (param == 2 && btype < 5) {
        for (int i = 0; i < num_particles-1; i++) {
          float x = random(1);
          float y = random(1);
          Particle particle = new Particle(new PVector(x, y), new PVector(v_ini, 0), new PVector(0.0, 0.0), 1./num_particles, 1);
          particles.add(particle);
        }
        Particle particle = new Particle(new PVector(0.01, 0.5), new PVector(v_ini, -v_ini), new PVector(0.0, 0.0), 1./num_particles, 1);
        particles.add(particle);
      } else if (param == 0 && btype < 5) {
        float spacing = 1. / iter;
        for (int x = 0; x < iter; x++) {
          for (int y = 0; y < iter; y++) {
            if ((x == floor(iter/3)) && (y == floor(iter/2))) {
              Particle particle = new Particle(new PVector(spacing * x + 0.5 * spacing, spacing * y + 0.5 * spacing), new PVector(v_ini, 0), new PVector(0, 0), 1./num_particles, e_ini);
              particles.add(particle);
            } else {
              PVector pos = new PVector(spacing * x + 0.5 * spacing + randomGaussian() / 1000, spacing * y + 0.5 * spacing + randomGaussian() / 1000);
              Particle particle = new Particle(pos, new PVector(v_ini, 0), new PVector(0, 0), 1./num_particles, 1.);
              particles.add(particle);
            }
          }
        }
      } else if (btype == 5) {
        println("Generating just a few particles.");
        for (int i = 0; i < 10; i++) {
          float x = random(1);
          float y = random(1);
          Particle particle = new Particle(new PVector(x, y), new PVector(0.0, 0.0), new PVector(0.0, 0.0), 1./num_particles, 1);
          particles.add(particle);
        }
      } else if (btype >= 6) {
        println("Generating wing particles.");
        for (int i = 0; i < pow(iter, 2); i++) {
          float x = random(0.1);
          float y = random(1);
          Particle particle = new Particle(new PVector(x, y), new PVector(v_ini, 0.0), new PVector(0.0, 0.0), 1./num_particles, 1);
          particles.add(particle);
        }
      }
    }
    println("Done creating particles!\n");
  }

  void createBoundaries() {

    if (btype == 0) { 
      // normal upper and lower boundaries (base)
      boundaries.add(new Boundary(new PVector(2.0, 0.0), new PVector(-1.0, 0.0)));
      boundaries.add(new Boundary(new PVector(-1.0, 1.0), new PVector(2.0, 1.0)));
      //boundaries.add(new Boundary(new PVector(-1.0, 0.0), new PVector(2.0, 0.0)));
      //boundaries.add(new Boundary(new PVector(2.0, 1.0), new PVector(-1.0, 1.0)));
    } else if (btype == 1) {
      // base combined with a line segment angled at 45 degrees
      boundaries.add(new Boundary(new PVector(2.0, 0.0), new PVector(-1.0, 0.0)));
      boundaries.add(new Boundary(new PVector(-1.0, 1.0), new PVector(2.0, 1.0)));
      //boundaries.add(new Boundary(new PVector(-1.0, 0.0), new PVector(2.0, 0.0)));
      //boundaries.add(new Boundary(new PVector(2.0, 1.0), new PVector(-1.0, 1.0)));
      //boundaries.add(new Boundary(new PVector(0.4, 0.4), new PVector(0.6, 0.6)));
      boundaries.add(new Boundary(new PVector(0.6, 0.6), new PVector(0.4, 0.4)));
    } else if (btype == 2) {
      // base combined with two line segment angled at 45 degrees each
      boundaries.add(new Boundary(new PVector(2.0, 0.0), new PVector(-1.0, 0.0)));
      boundaries.add(new Boundary(new PVector(-1.0, 1.0), new PVector(2.0, 1.0)));
      //boundaries.add(new Boundary(new PVector(-1.0, 0.0), new PVector(2.0, 0.0)));
      //boundaries.add(new Boundary(new PVector(2.0, 1.0), new PVector(-1.0, 1.0)));
      boundaries.add(new Boundary(new PVector(0.3, 0.6), new PVector(0.2, 0.5)));
      boundaries.add(new Boundary(new PVector(0.2, 0.5), new PVector(0.3, 0.4)));
      //boundaries.add(new Boundary(new PVector(0.2, 0.5), new PVector(0.3, 0.6)));
      //boundaries.add(new Boundary(new PVector(0.3, 0.4), new PVector(0.2, 0.5)));
    } else if (btype == 3) {
      // base combined with a straight line element
      boundaries.add(new Boundary(new PVector(2.0, 0.0), new PVector(-1.0, 0.0)));
      boundaries.add(new Boundary(new PVector(-1.0, 1.0), new PVector(2.0, 1.0)));
      boundaries.add(new Boundary(new PVector(0.3, 0.6), new PVector(0.3, 0.4)));
    } else if (btype == 4) {
      // base combined with simplified "rocket thruster" with iterative particle placement
      boundaries.add(new Boundary(new PVector(2.0, 0.0), new PVector(-1.0, 0.0)));
      boundaries.add(new Boundary(new PVector(-1.0, 1.0), new PVector(2.0, 1.0)));
      boundaries.add(new Boundary(new PVector(0.0, 0.0), new PVector(0.4, 0.45)));
      boundaries.add(new Boundary(new PVector(0.0, 1.0), new PVector(0.4, 0.55)));
      boundaries.add(new Boundary(new PVector(0.4, 0.45), new PVector(0.6, 0.4)));
      boundaries.add(new Boundary(new PVector(0.4, 0.55), new PVector(0.6, 0.6)));

      // Reversely orientated boundaries. Not necessary?
      boundaries.add(new Boundary(new PVector(0.4, 0.45), new PVector(0.0, 0.0)));
      boundaries.add(new Boundary(new PVector(0.4, 0.55), new PVector(0.0, 1.0)));
      boundaries.add(new Boundary(new PVector(0.6, 0.4), new PVector(0.4, 0.45)));
      boundaries.add(new Boundary(new PVector(0.6, 0.6), new PVector(0.4, 0.55)));
    } else if (btype == 5) {
      // same as btype = 4
      //combustion champer
      boundaries.add(new Boundary(new PVector(0.3, 0.29), new PVector(0.4, 0.45)));
      boundaries.add(new Boundary(new PVector(0.3, 0.71), new PVector(0.4, 0.55)));
      boundaries.add(new Boundary(new PVector(0.09, 0.3), new PVector(0.31, 0.3)));
      boundaries.add(new Boundary(new PVector(0.09, 0.7), new PVector(0.31, 0.7)));
      boundaries.add(new Boundary(new PVector(0.1, 0.3), new PVector(0.1, 0.7)));
      //boundaries.add(new Boundary(new PVector(0.1, 0.45), new PVector(0.1, 0.6)));
      //boundaries.add(new Boundary(new PVector(0.1, 0.65), new PVector(0.1, 0.7)));

      //methan oxygen pipes
      //boundaries.add(new Boundary(new PVector(0.4, 0.45), new PVector(0.6, 0.4)));
      //boundaries.add(new Boundary(new PVector(0.4, 0.55), new PVector(0.6, 0.6)));

      //nozzle





      // Reversely orientated boundaries. Not necessary?
      //boundaries.add(new Boundary(new PVector(0.4, 0.45), new PVector(0.0, 0.0)));
      //boundaries.add(new Boundary(new PVector(0.4, 0.55), new PVector(0.0, 1.0)));
      boundaries.add(new Boundary(new PVector(0.6, 0.4), new PVector(0.4, 0.45)));
      boundaries.add(new Boundary(new PVector(0.6, 0.6), new PVector(0.4, 0.55)));
    } else if (btype == 6) {
      //Wing
      boundaries.add(new Boundary(new PVector(0.32, 0.52), new PVector(0.375, 0.5)));
      boundaries.add(new Boundary(new PVector(0.32, 0.52), new PVector(0.315, 0.55)));
      boundaries.add(new Boundary(new PVector(0.375, 0.5), new PVector(0.425, 0.5)));
      boundaries.add(new Boundary(new PVector(0.425, 0.5), new PVector(0.7, 0.6)));
      boundaries.add(new Boundary(new PVector(0.315, 0.55), new PVector(0.7, 0.6)));
      boundaries.add(new Boundary(new PVector(-0.1, 0.0), new PVector(1.1, 0.0)));
      boundaries.add(new Boundary(new PVector(-0.1, 1.0), new PVector(1.1, 1.0)));
    } else if (btype == 7) {
      ArrayList<PVector> points = new ArrayList<PVector>();
      //Wing
      PVector p1 = new PVector(0.76, 0.6);
      PVector p2 = new PVector(0.4, 0.515);
      PVector p3 = new PVector(0.32, 0.515);
      PVector p4 = new PVector(0.3, 0.5125);
      PVector p5 = new PVector(0.285, 0.5);
      PVector p6 = new PVector(0.2963, 0.4875);
      PVector p7 = new PVector(0.325, 0.47);
      PVector p8 = new PVector(0.375, 0.46);
      PVector p9 = new PVector(0.45, 0.46);
      PVector p10 = new PVector(0.55, 0.48);
      points.add(p1);
      points.add(p2);
      points.add(p3);
      points.add(p4);
      points.add(p5);
      points.add(p6);
      points.add(p7);
      points.add(p8);
      points.add(p9);
      points.add(p10);
      float factor = 0.7;
      PVector center = new PVector(0.5, 0.5);
      for (PVector p : points) {
        PVector shifted = PVector.sub(p, center);
        boundaries.add(new Boundary(p1, p2));
        p.set(PVector.add(PVector.mult(shifted, factor), center));
      }
      
      for (int i = 0; i < 10; i++){
        boundaries.add(new Boundary(points.get(i), points.get((i+1)%10)));
      }

      boundaries.add(new Boundary(new PVector(-0.1, 0.0), new PVector(1.1, 0.0)));
      boundaries.add(new Boundary(new PVector(-0.1, 1.0), new PVector(1.1, 1.0)));
      boundaries.add(new Boundary(new PVector(0.0, 0.0), new PVector(0.0, 1.0)));
    }
  }

  int getNexdim(int dimension) {
    if (dim) {
      return (dimension + 1)%3;
    } else {
      return (dimension + 1)%2;
    }
  }


  int treeBuild(Node node, int dimension) {
    if (node.isleaf(leaf_size)) {
      //for (int i = node.start; i < node.end; i++) {
      //  Particle p = particles.get(i);
      //  node.m += p.m;
      //  node.cm = PVector.add(PVector.mult(p.pos, p.m), node.cm);
      //}
      //node.cm = PVector.mult(node.cm, 1/node.m);
      node.leaf = true;
      return 0;
    }
    int nexdim = getNexdim(dimension);
    float v;
    int s;
    if (dimension == 0) {
      v = node.center.x;
      s = partition(node.start, node.end, v, dimension);
      if ((s > node.start) && (s < node.end)) {
        PVector rhigh = node.rhigh.copy();
        PVector rlow = node.rlow.copy();
        rhigh.x = v;
        rlow.x = v;
        Node left = new Node(node.start, s, node.rlow, rhigh, dim);
        Node right = new Node(s, node.end, rlow, node.rhigh, dim);
        node.left = left;
        node.right = right;
        treeBuild(node.left, nexdim);
        treeBuild(node.right, nexdim);
        //node.m += node.left.m + node.right.m;
        //PVector leftcm = PVector.mult(node.left.cm, node.left.m);
        //PVector rightcm = PVector.mult(node.right.cm, node.right.m);
        //node.cm = PVector.div(PVector.add(leftcm, rightcm), node.m);
      } else if (!(s > node.start)) {
        node.rlow.x = v;
        node.center = PVector.div(PVector.add(node.rlow, node.rhigh), 2);
        treeBuild(node, nexdim);
      } else if (!(s < node.end)) {
        node.rhigh.x = v;
        node.center = PVector.div(PVector.add(node.rlow, node.rhigh), 2);
        treeBuild(node, nexdim);
      }
    } else if (dimension == 1) {
      v = node.center.y;
      s = partition(node.start, node.end, v, dimension);
      if ((s > node.start) && (s < node.end)) {
        PVector rhigh = node.rhigh.copy();
        PVector rlow = node.rlow.copy();
        rhigh.y = v;
        rlow.y = v;
        Node left = new Node(node.start, s, node.rlow, rhigh, dim);
        Node right = new Node(s, node.end, rlow, node.rhigh, dim);
        node.left = left;
        node.right = right;
        treeBuild(node.left, nexdim);
        treeBuild(node.right, nexdim);
        //node.m += node.left.m + node.right.m;
        //PVector leftcm = PVector.mult(node.left.cm, node.left.m);
        //PVector rightcm = PVector.mult(node.right.cm, node.right.m);
        //node.cm = PVector.div(PVector.add(leftcm, rightcm), node.m);
      } else if (!(s > node.start)) {
        node.rlow.y = v;
        node.center = PVector.div(PVector.add(node.rlow, node.rhigh), 2);
        treeBuild(node, nexdim);
      } else if (!(s < node.end)) {
        node.rhigh.y = v;
        node.center = PVector.div(PVector.add(node.rlow, node.rhigh), 2);
        treeBuild(node, nexdim);
      }
    } else {
      println("3rd dimension cut?");
      exit();
      v = node.center.z;
      s = partition(node.start, node.end, v, dimension);
      if ((s > node.start) && (s < node.end)) {
        PVector rhigh = node.rhigh.copy();
        PVector rlow = node.rlow.copy();
        rhigh.z = v;
        rlow.z = v;
        Node left = new Node(node.start, s, node.rlow, rhigh, dim);
        Node right = new Node(s, node.end, rlow, node.rhigh, dim);
        node.left = left;
        node.right = right;
        treeBuild(node.left, nexdim);
        treeBuild(node.right, nexdim);
        //node.m += node.left.m + node.right.m;
        //PVector leftcm = PVector.mult(node.left.cm, node.left.m);
        //PVector rightcm = PVector.mult(node.right.cm, node.right.m);
        //node.cm = PVector.div(PVector.add(leftcm, rightcm), node.m);
      } else if (!(s > node.start)) {
        node.rlow.z = v;
        node.center = PVector.div(PVector.add(node.rlow, node.rhigh), 2);
        treeBuild(node, nexdim);
      } else if (!(s < node.end)) {
        node.rhigh.z = v;
        node.center = PVector.div(PVector.add(node.rlow, node.rhigh), 2);
        treeBuild(node, nexdim);
      }
    }
    return 0;
  }

  int partition(int start, int end, float v, int dim) {
    int i = start;
    int j = end - 1;

    if (dim == 0) {
      while ((i <= j) && (particles.get(i).pos.x < v)) {
        i++;
      }
      while ((i <= j) && (particles.get(j).pos.x >= v)) {
        j--;
      }
      while (i < j) {
        Particle t = particles.get(i);
        particles.set(i, particles.get(j));
        particles.set(j, t);
        i++;
        j--;
        while (particles.get(i).pos.x < v) {
          i++;
        }
        while (particles.get(j).pos.x >= v) {
          j--;
        }
      }
      return i;
    } else if (dim == 1) {
      while ((i <= j) && (particles.get(i).pos.y < v)) {
        i++;
      }
      while ((i <= j) && (particles.get(j).pos.y >= v)) {
        j--;
      }
      while (i < j) {
        Particle t = particles.get(i);
        particles.set(i, particles.get(j));
        particles.set(j, t);
        i++;
        j--;
        while (particles.get(i).pos.y < v) {
          i++;
        }
        while (particles.get(j).pos.y >= v) {
          j--;
        }
      }
      return i;
    } else {
      while ((i <= j) && (particles.get(i).pos.z < v)) {
        i++;
      }
      while ((i <= j) && (particles.get(j).pos.z >= v)) {
        j--;
      }
      while (i < j) {
        Particle t = particles.get(i);
        particles.set(i, particles.get(j));
        particles.set(j, t);
        i++;
        j--;
        while (particles.get(i).pos.z < v) {
          i++;
        }
        while (particles.get(j).pos.z >= v) {
          j--;
        }
      }
      return i;
    }
  }


  void nn_density(Node root) {
    if (dim) {
      for (Particle p : particles) {
        p.nn_search_3d(root, nn, particles);
      }
    } else {
      for (Particle p : particles) {
        p.nn_search_2d(root, nn, particles, boundaries);
      }
    }
    if (dim) {
      for (Particle p : particles) {
        p.rho = 0;
        for (ParticleTuple tup : p.n_closest) {
          p.rho += tup.p.m * p.monoghan_kernel_3d(tup.dist, sigma);
        }
        if (p.rho > max_val) {
          max_val = p.rho;
        }
      }
    } else {
      max_val = 0;
      for (Particle p : particles) {
        p.rho = 0;
        for (ParticleTuple tup : p.n_closest) {
          if (!tup.p.isReflectedParticle) {
            p.rho += tup.p.m * p.monoghan_kernel_2d(tup.dist, sigma);
          }
        }
        if (p.rho > max_val) {
          max_val = p.rho;
        }
        // add small positive number to rho to avoid zero density
      }
    }
  }

  void calcsound() {
    for (Particle p : particles) {
      p.c_sound = pow(gamma * (gamma - 1) * p.e_pred, 0.5);
    }
  }

  void calc_dt() {
    float c_max = 0.;
    float h_min = MAX_FLOAT;
    for (Particle p : particles) {
      c_max = max(p.c_sound, c_max);
      h_min = min(p.h, h_min);
    }
    dt = h_min/c_max * courant;
    t_global += dt;
  }

  void reset_ae() {
    if (dim) {
      for (Particle p : particles) {
        p.a.set(0., 0., 0.);
        p.e_dot = 0.;
      }
    } else {
      for (Particle p : particles) {
        p.a.set(0., 0.);
        p.e_dot = 0.;
      }
    }
  }

  void reset_e() {
    for (Particle p : particles) {
      p.e = 1.;
    }
  }


  void nn_sphforce() {
    reset_ae();
    for (Particle p : particles) {
      p.calc_ae(gamma, sigma, dim);
    }
  }

  void calc_forces() {
    treeBuild(root, 0);
    nn_density(root);
    calcsound();
    calc_dt();
    nn_sphforce();
  }

  void boundary_condtions(Particle p) {

    if (dim) {
      p.pos.set((p.pos.x + 1) % 1., (p.pos.y + 1) % 1., (p.pos.z + 1) % 1.);
    } else if (!dim && btype >= 5) {
      if (p.pos.x > 1.0 || p.pos.x < 0.0) {
        garbage_colleciton.add(p);
      }
      if (p.pos.y > 1.0 || p.pos.y < 0.0) {
        garbage_colleciton.add(p);
      }
    } else {
      if (p.pos.x >= 1.0) {
        p.pos.set(0.0, random(1));
        p.vel.set(v_ini, 0);
        reset_e();
      } else if (p.pos.x <= 0.0) {
        p.pos.set(0.0, random(1));
        p.vel.set(v_ini, 0);
        reset_e();
      }
      //if (p.pos.y > 1.0) {
      //  p.pos.set(p.pos.x, 2 - p.pos.y);
      //  p.vel.set(p.vel.x, -p.vel.y);
      //} else if (p.pos.y < 0.0) {
      //  p.pos.set(p.pos.x, - p.pos.y);
      //  p.vel.set(p.vel.x, -p.vel.y);
      //}
      //p.pos.set((p.pos.x + 1) % 1., (p.pos.y + 1) % 1.);
      //p.pos.set(p.pos.x, (p.pos.y + 1) % 1.);
    }
  }

  void ignite() {
    if (btype == 5) {      
      for (Particle p : particles) {
        if ( p.pos.x >= .38 && p.pos.x <= .42) {
          if ( p.pos.y >= 0.48 && p.pos.y <= 0.52) {
            p.e *= 1.01;
          }
        }
      }
    }
  }

  void drift1() {
    for (Particle p : particles) {

      // Calculate the next position
      p.temp_pos = PVector.add(p.pos, PVector.mult(p.vel, dt * 0.5));

      // Check if crossed a boundary
      for (Boundary b : boundaries) {
        b.checkBoundary(p);
      }
      p.pos = p.temp_pos;
      boundary_condtions(p);
      p.v_pred = PVector.add(p.vel, PVector.mult(p.a, dt * 0.5));
      p.e_pred = p.e + p.e_dot * 0.5 * dt;
    }
  }

  void drift2() {
    for (Particle p : particles) {

      // Calculate the next position
      p.temp_pos = PVector.add(p.pos, PVector.mult(p.vel, dt * 0.5));

      // Checked if crossed a boundary
      for (Boundary b : boundaries) {
        b.checkBoundary(p);
      }
      p.pos = p.temp_pos;
      boundary_condtions(p);
    }
  }

  void kick() {
    for (Particle p : particles) {
      p.vel.add(PVector.mult(p.a, dt));
      p.e += p.e_dot * dt;

      //p.vel.add(0, v_ini * dt);
    }
  }

  void update() {
    particles.removeAll(garbage_colleciton);
    garbage_colleciton.clear();
    if (btype == 5 && max_val <= 20) {
      //if (btype == 5 && particles.size() < pow(iter, 2)) {
      float x = random(0.15, 0.2);
      float y = random(0.4, 0.6);
      particles.add(new Particle(new PVector(x, y), new PVector(v_ini, 0.0), new PVector(0.0, 0.0), 1./num_particles, 50));
      //ignite();
    } else if (btype >= 6) {
      while (particles.size() < pow(iter, 2)) {
        float x = 0.001;
        float y = random(1);
        particles.add(new Particle(new PVector(x, y), new PVector(v_ini, 0.0), new PVector(0.0, 0.0), 1./num_particles, 1));
        //ignite();
      }
    }
    this.root = new Node(0, particles.size(), new PVector(0, 0), new PVector(1, 1), dim);
    drift1();
    ignite();
    calc_forces();
    kick();
    drift2();
  }

  void show_particles() {
    root.show(size);
    for (Boundary b : boundaries) {
      b.drawBoundary(size);
    }
    if (dim) {
      for (Particle p : particles) {
        p.show_3d(size, max_val);
      }
    } else {  
      for (Particle p : particles) {
        p.show_2d(size, max_val);
      }
    }
  }
}
