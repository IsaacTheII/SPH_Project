class Particle {
  PVector pos;
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
  ArrayList<ParticleTuple> n_closest;

  Particle(PVector pos_, PVector vel_, PVector a_, float m_, float e_) {
    pos = pos_;
    vel = vel_;
    a = a_;
    v_pred = vel_;
    m = m_;
    e = e_;
    e_pred = e_;
  }

  ArrayList<NodeTuple> insert_into_sort_list(ArrayList<NodeTuple> list_tup, NodeTuple tup) {
    int range = list_tup.size();
    for (int i = 0; i < range; i++) {
      NodeTuple temp = list_tup.get(i);
      if (temp.dist > tup.dist) {
        println("before", list_tup.get(i).dist);
        list_tup.add(i, tup);
        println("after", list_tup.get(i).dist);
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


  void nn_search_2d(Node node, int nn, ArrayList<Particle> particles) {
    ArrayList<NodeTuple> node_q_sorted = new ArrayList<NodeTuple>();
    for (int x_off = -1; x_off < 2; x_off++) {
      for (int y_off = -1; y_off < 2; y_off++) {
        PVector offset = PVector.sub(node.rhigh, node.rhigh);
        offset.x = offset.x * x_off;
        offset.y = offset.y * y_off;
        node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node, node.min_node_dist_pb(pos, offset), offset));
      }
    }
    n_closest = new ArrayList<ParticleTuple>();
    for (int i = 0; i <= nn; i++) {
      n_closest.add(new ParticleTuple(this, MAX_FLOAT, new PVector(0, 0)));
    }
    //ParticleTuple furthest = n_closest.get(n_closest.size() - 1);
    while (n_closest.get(n_closest.size() - 1).dist > node_q_sorted.get(0).dist) {
      NodeTuple temp = node_q_sorted.get(0);
      node_q_sorted.remove(0);
      if (temp.node.leaf) {
        for (int i = node.start; i < node.end; i++) {
          PVector distance = PVector.sub(pos, PVector.add(particles.get(i).pos, temp.offset));
          if (n_closest.get(n_closest.size() - 1).dist > distance.mag()) {
            n_closest.remove(n_closest.size() - 1);
            n_closest = insert_into_sort_list(n_closest, new ParticleTuple(particles.get(i), distance.mag(), distance));
          }
        }
      } else {
        node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node.left, node.left.min_node_dist_pb(pos, temp.offset), temp.offset));
        node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node.right, node.right.min_node_dist_pb(pos, temp.offset), temp.offset));
      }
      if (node_q_sorted.size() <= 0) {
        break;
      }
    }
  }

  void nn_search_3d(Node node, int nn, ArrayList<Particle> particles) {
    ArrayList<NodeTuple> node_q_sorted = new ArrayList<NodeTuple>();

    for (int x_off = -1; x_off < 2; x_off++) {
      for (int y_off = -1; y_off < 2; y_off++) {
        for (int z_off = -1; z_off < 2; z_off++) {

          PVector offset = PVector.sub(node.rhigh, node.rlow);
          offset.x = offset.x * x_off;
          offset.y = offset.y * y_off;
          offset.z = offset.z * z_off;
          node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node, node.min_node_dist_pb(pos, offset), offset));
        }
      }
    }
    n_closest = new ArrayList<ParticleTuple>();
    for (int i = 0; i <= nn; i++) {
      n_closest.add(new ParticleTuple(this, MAX_FLOAT, new PVector(0, 0, 0)));
    }
    while (n_closest.get(n_closest.size() - 1).dist > node_q_sorted.get(0).dist) {
      NodeTuple temp = node_q_sorted.get(0);
      node_q_sorted.remove(0);
      println(temp.node.leaf);

      if (temp.node.leaf) {
        for (int i = node.start; i < node.end; i++) {
          PVector distance = PVector.sub(pos, PVector.add(particles.get(i).pos, temp.offset));
          if (n_closest.get(n_closest.size() - 1).dist > distance.mag()) {
            n_closest.remove(n_closest.size() - 1);
            n_closest = insert_into_sort_list(n_closest, new ParticleTuple(particles.get(i), distance.mag(), distance));
          }
        }
      } else {
        node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node.left, node.left.min_node_dist_pb(pos, temp.offset), temp.offset));
        node_q_sorted = insert_into_sort_list(node_q_sorted, new NodeTuple(node.right, node.right.min_node_dist_pb(pos, temp.offset), temp.offset));
      }

      if (node_q_sorted.size() <= 0) {

        break;
      }
    }
    println(n_closest.get(n_closest.size() - 1).dist);

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



  void show(int size, float max_rho) {
    float col = map(rho, 0, max_rho, 0, 1);
    stroke(col, 1, 1);
    strokeWeight(4);
    float x = map(pos.x, 0, 1, -size/2, size/2);
    float y = map(pos.y, 0, 1, -size/2, size/2);
    point(x, y);
  }

  void show_3d(int size, float max_rho) {
    float col = map(rho, 0, max_rho, 0, 1);
    stroke(col, 1, 1);
    strokeWeight(4);
    float x = map(pos.x, 0, 1, -size/2, size/2);
    float y = map(pos.y, 0, 1, -size/2, size/2);
    float z = map(pos.z, 0, 1, -size/2, size/2);
    point(x, y, z);
  }
}
