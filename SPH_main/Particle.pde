class Particle {
  PVector pos;
  PVector temp_pos;
  PVector vel;
  PVector v_pred;
  PVector a;
  float m;
  float e;
  float e_pred;
  float e_dot = 0;
  float h;
  float rho;
  float c_sound;
  float pa;
  ArrayList<ParticleTuple> n_closest;
  boolean reflected = false;

  Particle(PVector pos_, PVector vel_, PVector a_, float m_, float e_) {
    pos = pos_;
    temp_pos = pos;
    vel = vel_;
    a = a_;
    v_pred = vel_;
    m = m_;
    e = e_;
    e_pred = e_;
  }


  // heap sort
  //ArrayList<NodeTuple> heapify(ArrayList<NodeTuple> list_tup, int i) {
  //  int m = i;
  //  int l = 2 * i + 1;
  //  int r = 2 * (i + 1);
  //  if ((l < list_tup.size()) && (list_tup.get(l).dist < list_tup.get(m).dist)) {
  //    m = l;
  //  }
  //  if ((r < list_tup.size()) && (list_tup.get(r).dist < list_tup.get(m).dist)) {
  //    m = l;
  //  }
  //  if (i != m) {
  //    NodeTuple temp = list_tup.get(m);
  //    list_tup.set(m, list_tup.get(i));
  //    list_tup.set(i, temp);
  //  }
  //  return list_tup;
  //}
  //
  //ArrayList<NodeTuple> insert_into_sort_list(ArrayList<NodeTuple> list_tup, NodeTuple tup) {
  //  list_tup.add(tup);
  //  for (int i = floor(list_tup.size()/2)-1; i >= 0; i--) {
  //    list_tup = heapify(list_tup, i);
  //  }
  //  return list_tup;
  //}

  ArrayList<NodeTuple> insert_into_sort_list(ArrayList<NodeTuple> list_tup, NodeTuple tup) {
    for (int i = 0; i < list_tup.size(); i++) {
      NodeTuple temp = list_tup.get(i);
      if (temp.dist > tup.dist) {
        list_tup.add(i, tup);
        return list_tup;
      }
    }
    list_tup.add(tup);
    return list_tup;
  }

  ArrayList<ParticleTuple> insert_into_sort_list(ArrayList<ParticleTuple> list_tup, ParticleTuple tup) {
    for (int i = 0; i < list_tup.size(); i++) {
      if (list_tup.get(i).dist > tup.dist) {
        list_tup.add(i, tup);
        return list_tup;
      }
    }
    list_tup.add(tup);
    return list_tup;
  }

  void insert_in_q(PVector distance, Particle p) {
    if (n_closest.get(n_closest.size() - 1).dist > distance.mag()) {
      n_closest.remove(n_closest.size() - 1);
      insert_into_sort_list(n_closest, new ParticleTuple(p, distance.mag(), distance));
    }
  }






  void nn_search_2d(Node node, int nn, ArrayList<Particle> particles, ArrayList<Boundary> boundaries) {
    ArrayList<NodeTuple> node_q_sorted = new ArrayList<NodeTuple>();
    //for (int x_off = -1; x_off < 2; x_off++) {
    //  for (int y_off = -1; y_off < 2; y_off++) {
    //    PVector offset = new PVector(x_off, y_off);
    //    node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node, node.min_node_dist_pb(pos, offset), offset));
    //  }
    //}

    for (int x_off = -1; x_off < 2; x_off++) {
      PVector offset = new PVector(x_off, 0);
      node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node, node.min_node_dist_pb(pos, offset), offset));
    }


    n_closest = new ArrayList<ParticleTuple>();
    for (int i = 0; i <= nn; i++) {
      n_closest.add(new ParticleTuple(this, MAX_FLOAT, new PVector(0, 0)));
    }
    NodeTuple temp = node_q_sorted.get(0);
    node_q_sorted.remove(0);
    while (n_closest.get(n_closest.size() - 1).dist > temp.dist) {
      if (temp.node.leaf) {
        for (int i = temp.node.start; i < temp.node.end; i++) {
          PVector distance = PVector.sub(pos, PVector.add(particles.get(i).pos, temp.offset));
          insert_in_q(distance, particles.get(i));
        }
      } else {
        NodeTuple left = new NodeTuple(temp.node.left, temp.node.left.min_node_dist_pb(pos, temp.offset), temp.offset);
        NodeTuple right = new NodeTuple(temp.node.right, temp.node.right.min_node_dist_pb(pos, temp.offset), temp.offset);
        node_q_sorted = insert_into_sort_list(node_q_sorted, left);
        node_q_sorted = insert_into_sort_list(node_q_sorted, right);
      }
      if (node_q_sorted.size() <= 0) {
        break;
      } else {
        temp = node_q_sorted.get(0);
        node_q_sorted.remove(0);
      }
    }
    h = n_closest.get(n_closest.size() - 1).dist;

    for (Boundary b : boundaries) {
      //PVector P1 = b.startPoint;
      //PVector P2 = b.endPoint;
      //float nominator = (P2.y - P1.y) * pos.x - (P2.x - P1.x) * pos.y + P2.x * P1.y - P2.y * P1.x;
      //float denominator = pow(pow(P2.y - P1.y, 2) + pow(P2.x - P1.x, 2), 0.5);
      //float d_new = abs(nominator)/ denominator;
      float d_new = b.distanceToBoundary(this);
      if (d_new < h) {
        ArrayList<Particle> reflected = b.reflectParticles(n_closest);
        for (Particle p : reflected) {
          PVector distance = PVector.sub(pos, p.pos);
          insert_in_q(distance, p);
        }

        if (h < n_closest.get(n_closest.size() - 1).dist) {        
          println("h before: ", h);
          println("h after:  ", n_closest.get(n_closest.size() - 1).dist);
        }
        h = n_closest.get(n_closest.size() - 1).dist;
      }
    }
  }



  void nn_search_3d(Node node, int nn, ArrayList<Particle> particles) {
    ArrayList<NodeTuple> node_q_sorted = new ArrayList<NodeTuple>();
    for (int x_off = -1; x_off < 2; x_off++) {
      for (int y_off = -1; y_off < 2; y_off++) {
        for (int z_off = -1; z_off < 2; z_off++) {
          PVector offset = new PVector(x_off, y_off, z_off);
          node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node, node.min_node_dist_pb(pos, offset), offset));
        }
      }
    }
    n_closest = new ArrayList<ParticleTuple>();
    for (int i = 0; i <= nn; i++) {
      n_closest.add(new ParticleTuple(this, MAX_FLOAT, new PVector(0, 0, 0)));
    }
    NodeTuple temp = node_q_sorted.get(0);
    node_q_sorted.remove(0);
    while (n_closest.get(n_closest.size() - 1).dist > temp.dist) {
      if (temp.node.leaf) {
        for (int i = temp.node.start; i < temp.node.end; i++) {
          PVector distance = PVector.sub(pos, PVector.add(particles.get(i).pos, temp.offset));
          if (n_closest.get(n_closest.size() - 1).dist > distance.mag()) {
            n_closest.remove(n_closest.size() - 1);
            n_closest = insert_into_sort_list(n_closest, new ParticleTuple(particles.get(i), distance.mag(), distance));
          }
        }
      } else {
        NodeTuple left = new NodeTuple(temp.node.left, temp.node.left.min_node_dist_pb(pos, temp.offset), temp.offset);
        NodeTuple right = new NodeTuple(temp.node.right, temp.node.right.min_node_dist_pb(pos, temp.offset), temp.offset);
        node_q_sorted = insert_into_sort_list(node_q_sorted, left);
        node_q_sorted = insert_into_sort_list(node_q_sorted, right);
      }
      if (node_q_sorted.size() <= 0) {
        break;
      } else {
        temp = node_q_sorted.get(0);
        node_q_sorted.remove(0);
      }
    }
    h = n_closest.get(n_closest.size() - 1).dist;
  }

  float monoghan_kernel_3d(float dist, float sigma) {
    float dist_h = dist/h;
    if (dist_h < 0.5) {
      return pow(sigma/h, 3) * (6 * pow(dist_h, 3) - 6 * pow(dist_h, 2) + 1);
    } else if (dist_h <= 1) {
      return pow(sigma/h, 3) * 2 * pow(1 - dist_h, 3);
    } else {
      return 0;
    }
  }

  float monoghan_kernel_2d(float dist, float sigma) {
    float dist_h = dist/h;
    if (dist_h < 0.5) {
      return pow(sigma/h, 2) * (6 * pow(dist_h, 3) - 6 * pow(dist_h, 2) + 1);
    } else if (dist_h <= 1) {
      return pow(sigma/h, 2) * 2 * pow(1 - dist_h, 3);
    } else {
      return 0;
    }
  }

  PVector grad_monoghan_kernel_3d(float dist, PVector dist_v, float sigma) {
    float dist_h = dist/h;
    PVector result = new PVector(0, 0, 0);
    if (dist_h < 0.5) {
      dist_v.normalize(result);
      return result.mult(6 * pow(sigma/h, 4) * (3 * pow(dist_h, 2)  - 2 * dist_h));
    } else if (dist_h <= 1.) {
      dist_v.normalize(result);
      return result.mult(6 * pow(sigma/h, 4) * -1 * pow(1 - dist_h, 2));
    } else {
      return result;
    }
  }

  PVector grad_monoghan_kernel_2d(float dist, PVector dist_v, float sigma) {
    float dist_h = dist/h;
    PVector result = new PVector(0, 0);
    if (dist_h < 0.5) {
      dist_v.normalize(result);
      return result.mult(6 * pow(sigma/h, 3) * (3 * pow(dist_h, 2)  - 2 * dist_h));
    } else if (dist_h <= 1.) {
      dist_v.normalize(result);
      return result.mult(6 * pow(sigma/h, 3) * -1. * pow(1. - dist_h, 2));
    } else {
      return result;
    }
  }

  float calc_viscosity(Particle other, PVector dist) {
    float alpha = 1.;
    float beta = 2.;
    float eta2 = 0.001;
    PVector vab = PVector.sub(vel, other.vel);
    if (vab.dot(dist) < 0) {
      float muab = ((h + other.h) * 0.5 * vab.dot(dist)) / (dist.dot(dist) + eta2);
      return -alpha * (c_sound + other.c_sound) * 0.5 * muab + beta * pow(muab, 2);
    } else {
      return 0.;
    }
  }


  void calc_ae(float gamma, float sigma, boolean dim) {
    pa = pow(c_sound, 2)/(gamma * rho);
    for (ParticleTuple closest : n_closest) {
      float pb = pow(closest.p.c_sound, 2)/(gamma * closest.p.rho);
      float pi_ab = calc_viscosity(closest.p, closest.offset);
      PVector vab = PVector.sub(vel, closest.p.vel);
      float Fab = pa + pb + pi_ab;
      PVector contribution_a;
      if (dim) {
        contribution_a = grad_monoghan_kernel_3d(closest.dist, closest.offset, sigma);
      } else {
        contribution_a = grad_monoghan_kernel_2d(closest.dist, closest.offset, sigma);
      }
      a.add(PVector.mult(contribution_a, -0.5 * closest.p.m * Fab));
      closest.p.a.add(PVector.mult(contribution_a, 0.5 * m * Fab));

      e_dot += closest.p.m * (pa + pi_ab) * vab.dot(PVector.mult(contribution_a, 0.5));
      closest.p.e_dot += m * (pb + pi_ab) * vab.dot(PVector.mult(contribution_a, 0.5));
    }
  }

  void show_2d(int size, float max_val) {
    float col = map(pow(rho, 1), 0, pow(max_val, 1), 0, 1);
    stroke(col, 1, 1);
    strokeWeight(map(pow(rho, 1), 0, pow(max_val, 1), 5, 20));
    float x = map(pos.x, 0, 1, -size/2, size/2);
    float y = map(pos.y, 0, 1, -size/2, size/2);
    point(x, y);
  }

  void show_3d(int size, float max_val) {
    float col = map(pow(rho, 1), 0, pow(max_val, 1), 0, 1);
    stroke(col, 1, 1);
    strokeWeight(map(pow(rho, 1), 0, pow(max_val, 1), 5, 20));
    float x = map(pos.x, 0, 1, -size/2, size/2);
    float y = map(pos.y, 0, 1, -size/2, size/2);
    float z = map(pos.z, 0, 1, -size/2, size/2);
    point(x, y, z);
  }
}
