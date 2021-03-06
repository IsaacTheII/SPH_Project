class Simulation {
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
  float max_val = 0;

  Simulation(int leaf_size_, int param_, int iter_, float e_ini_, int nn_, boolean dim_, float courant_, int size_) {
    leaf_size = leaf_size_;
    iter = iter_;
    e_ini = e_ini_;
    nn = nn_;
    dim = dim_;
    courant = courant_;
    size = size_;
    if (dim) {
      gamma = 5./3.;
      sigma = 8./PI;
      num_particles = int(pow(iter, 3));
    } else {
      gamma = 2.;
      sigma = 40./(7 * PI);
      num_particles = int(pow(iter, 2));
    }
    particles = new ArrayList<Particle>();
    PVector rlow = getRlow();
    PVector rhigh = getRhigh();
    root = new Node(0, num_particles, rlow, rhigh, dim);
    read_data(param_);
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
              if ((spacing * x + 0.5 * spacing == 0.5) && (spacing * y + 0.5 * spacing == 0.5) && (spacing * z + 0.5 * spacing == 0.5)) {
                Particle particle = new Particle(new PVector(spacing * x + 0.5 * spacing, spacing * y + 0.5 * spacing, spacing * z + 0.5 * spacing), new PVector(0, 0, 0), new PVector(0, 0, 0), 1./num_particles, e_ini);
                particles.add(particle);
              } else {
                Particle particle = new Particle(new PVector(spacing * x + 0.5 * spacing, spacing * y + 0.5 * spacing, spacing * z + 0.5 * spacing), new PVector(0, 0, 0), new PVector(0, 0, 0), 1./num_particles, 1.);
                particles.add(particle);
              }
            }
          }
        }
      }
    } else {
      if (param == 1) {
        for (int i = 0; i < num_particles-1; i++) {
          float x = random(1);
          float y = random(1);
          Particle particle = new Particle(new PVector(x, y), new PVector(0, 0), new PVector(0, 0), 1./num_particles, 1);
          particles.add(particle);
        }
        Particle particle = new Particle(new PVector(0.5, 0.5), new PVector(0, 0), new PVector(0, 0), 1./num_particles, e_ini);
        particles.add(particle);
      } else {
        float spacing = 1. / iter;
        for (int x = 0; x < iter; x++) {
          for (int y = 0; y < iter; y++) {
            if ((spacing * x + 0.5 * spacing == 0.5) && (spacing * y + 0.5 * spacing == 0.5)) {
              Particle particle = new Particle(new PVector(spacing * x + 0.5 * spacing, spacing * y + 0.5 * spacing), new PVector(0, 0), new PVector(0, 0), 1./num_particles, e_ini);
              particles.add(particle);
            } else {
              Particle particle = new Particle(new PVector(spacing * x + 0.5 * spacing, spacing * y + 0.5 * spacing), new PVector(0, 0), new PVector(0, 0), 1./num_particles, 1.);
              particles.add(particle);
            }
          }
        }
      }
    }
    println("Done creating particles!\n");
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
        treeBuild(node, nexdim);
      } else if (!(s < node.end)) {
        node.rhigh.x = v;
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
        treeBuild(node, nexdim);
      } else if (!(s < node.end)) {
        node.rhigh.y = v;
        treeBuild(node, nexdim);
      }
    } else {
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
        treeBuild(node, nexdim);
      } else if (!(s < node.end)) {
        node.rhigh.z = v;
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
        p.nn_search_2d(root, nn, particles);
      }
    }
    max_val = 0.;
    if (dim) {
      for (Particle p : particles) {
        p.rho = 0;
        for (ParticleTuple tup : p.n_closest) {
          p.rho += tup.p.m * p.monoghan_kernel_3d(tup.dist, sigma);
        }
        if (p.pa > max_val) {
          max_val = p.pa;
        }
      }
    } else {
      for (Particle p : particles) {
        p.rho = 0;
        for (ParticleTuple tup : p.n_closest) {
          p.rho += tup.p.m * p.monoghan_kernel_2d(tup.dist, sigma);
        }
        if (p.pa > max_val) {
          max_val = p.pa;
        }
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

  void drift1() {
    for (Particle p : particles) {
      p.pos.add(PVector.mult(p.vel, dt * 0.5));
      if (dim) {
        p.pos.set((p.pos.x + 1) % 1., (p.pos.y + 1) % 1., (p.pos.z + 1) % 1.);
      } else {
        p.pos.set((p.pos.x + 1) % 1., (p.pos.y + 1) % 1.);
      }
      p.v_pred = PVector.add(p.vel, PVector.mult(p.a, dt * 0.5));
      p.e_pred = p.e + p.e_dot * 0.5 * dt;
    }
  }

  void drift2() {
    for (Particle p : particles) {
      p.pos.add(PVector.mult(p.vel, dt * 0.5));
      if (dim) {
        p.pos.set((p.pos.x + 1) % 1., (p.pos.y + 1) % 1., (p.pos.z + 1) % 1.);
      } else {
        p.pos.set((p.pos.x + 1) % 1., (p.pos.y + 1) % 1.);
      }
    }
  }

  void kick() {
    for (Particle p : particles) {
      p.vel.add(PVector.mult(p.a, dt));
      p.e += p.e_dot * dt;
    }
  }

  void update() {
    drift1();
    calc_forces();
    kick();
    drift2();
  }



  void show_particles() {
    //root.show(size);
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
