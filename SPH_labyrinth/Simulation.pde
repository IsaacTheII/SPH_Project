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
  float t_global = 0.;
  float dt = 0.;
  float max_val = 10;
  float v_ini;
  float rotationSpeed = 0.002;
  int btype;
  ArrayList<Boundary> boundaries;
  PVector rlow;
  PVector rhigh;
  ArrayList<Particle> garbage_colleciton = new ArrayList<Particle>();
  int w;
  int h;
  float scale = 3.5;

  Simulation(int leaf_size_, int param_, int iter_, float e_ini_, int nn_, boolean dim_, float courant_, int size_, float v_ini_, int btype_, int w_, int h_) {
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
    w = w_;
    h = h_;
    gamma = 1.1;
    sigma = 40./(7 * PI);
    num_particles = int(pow(iter, 2));

    particles = new ArrayList<Particle>();
    rlow = new PVector(0, 0);
    rhigh = new PVector(w, h);
    ;
    root = new Node(0, 10, rlow, rhigh, dim); // hardcoded 10 particles to start with btype = 5
    read_data();
    createBoundaries();
  }



  void read_data() {
    println("Creating ", num_particles, " particles");

    randomSeed(0);
    for (int i = 0; i <= nn + 1; i++) particles.add(new Particle(new PVector(.34 + random(0.05), 6/scale - random(.9)/scale), new PVector(v_ini, 0.0), new PVector(0.0, 0.0), 1./num_particles, 1));
    println("Done creating particles!\n");
  }

  void createBoundaries() {
    boundaries.add(new Boundary(new PVector(1, 1).div(scale), new PVector(12, 1).div(scale)));
    boundaries.add(new Boundary(new PVector(13, 1).div(scale), new PVector(13, 6).div(scale)));
    boundaries.add(new Boundary(new PVector(13, 6).div(scale), new PVector(1, 6).div(scale)));
    boundaries.add(new Boundary(new PVector(1, 6).div(scale), new PVector(1, 1).div(scale)));

    boundaries.add(new Boundary(new PVector(1, 5).div(scale), new PVector(2, 5).div(scale)));
    boundaries.add(new Boundary(new PVector(1, 3).div(scale), new PVector(2, 3).div(scale)));
    boundaries.add(new Boundary(new PVector(2, 3).div(scale), new PVector(2, 4).div(scale)));
    boundaries.add(new Boundary(new PVector(2, 2).div(scale), new PVector(4, 2).div(scale)));
    boundaries.add(new Boundary(new PVector(3, 2).div(scale), new PVector(3, 4).div(scale)));
    boundaries.add(new Boundary(new PVector(3, 4).div(scale), new PVector(5, 4).div(scale)));
    boundaries.add(new Boundary(new PVector(4, 2).div(scale), new PVector(4, 3).div(scale)));
    boundaries.add(new Boundary(new PVector(4, 4).div(scale), new PVector(4, 5).div(scale)));
    boundaries.add(new Boundary(new PVector(4, 5).div(scale), new PVector(6, 5).div(scale)));
    boundaries.add(new Boundary(new PVector(5, 5).div(scale), new PVector(5, 6).div(scale)));
    boundaries.add(new Boundary(new PVector(5, 4).div(scale), new PVector(5, 3).div(scale)));
    boundaries.add(new Boundary(new PVector(5, 3).div(scale), new PVector(6, 3).div(scale)));
    boundaries.add(new Boundary(new PVector(6, 3).div(scale), new PVector(6, 2).div(scale)));
    boundaries.add(new Boundary(new PVector(6, 2).div(scale), new PVector(7, 2).div(scale)));
    boundaries.add(new Boundary(new PVector(5, 1).div(scale), new PVector(5, 2).div(scale)));
    boundaries.add(new Boundary(new PVector(8, 1).div(scale), new PVector(8, 3).div(scale)));

    boundaries.add(new Boundary(new PVector(7, 3).div(scale), new PVector(7, 4).div(scale)));
    boundaries.add(new Boundary(new PVector(7, 4).div(scale), new PVector(8, 4).div(scale)));
    boundaries.add(new Boundary(new PVector(8, 4).div(scale), new PVector(8, 5).div(scale)));
    boundaries.add(new Boundary(new PVector(5, 5).div(scale), new PVector(5, 6).div(scale)));
    boundaries.add(new Boundary(new PVector(8, 5).div(scale), new PVector(12, 5).div(scale)));

    boundaries.add(new Boundary(new PVector(9, 3).div(scale), new PVector(9, 4).div(scale)));
    boundaries.add(new Boundary(new PVector(6, 5).div(scale), new PVector(6, 4).div(scale)));
    boundaries.add(new Boundary(new PVector(7, 5).div(scale), new PVector(7, 6).div(scale)));
    boundaries.add(new Boundary(new PVector(12, 5).div(scale), new PVector(8, 5).div(scale)));
    boundaries.add(new Boundary(new PVector(10, 4).div(scale), new PVector(13, 4).div(scale)));
    boundaries.add(new Boundary(new PVector(7, 3).div(scale), new PVector(12, 3).div(scale)));
    boundaries.add(new Boundary(new PVector(12, 3).div(scale), new PVector(12, 2).div(scale)));
    boundaries.add(new Boundary(new PVector(9, 2).div(scale), new PVector(12, 2).div(scale)));
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
      if (p.pos.x >= w || p.pos.x <= 0.0 || p.pos.y >= h || p.pos.y <= 0.0) {
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

    if (max_val <= 100) {
      println(max_val);
      for (int i = 0; i <= 2; i++) particles.add(new Particle(new PVector(.34 + random(0.05), 6/scale - random(.9)/scale), new PVector(v_ini, 0.0), new PVector(0.0, 0.0), 1./num_particles, 1));
    }
    
    this.root = new Node(0, particles.size(), new PVector(0, 0), new PVector(w, h), dim);
    drift1();
    calc_forces();
    kick();
    drift2();
  }

  void show_particles() {
    //root.show(size, w);
    stroke(0, 0, 1);
    for (Boundary b : boundaries) {
      b.drawBoundary(size, w, h);
    } 
    for (Particle p : particles) {
      p.show_2d(size, max_val, w, h);
    }
  }
}
