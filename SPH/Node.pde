class Node {
  Node left;
  Node right;
  int start;
  int end;
  PVector rlow;
  PVector rhigh;
  PVector center;
  float rmax;
  float m = 0;
  PVector cm;
  float size;
  boolean leaf = false;
  boolean dim;

  Node(int start_, int end_, final PVector rlow_, final PVector rhigh_, boolean dim_) {
    start = start_;
    end = end_;      
    rlow = rlow_;
    rhigh = rhigh_;
    center = PVector.add(rlow, rhigh);
    center = PVector.div(center, 2);
    PVector dist = PVector.sub(rhigh_, rlow_);
    rmax = dist.mag()/2;
    dim = dim_;
    if (dim) {
      size = max(dist.x, dist.y, dist.z);
      cm = new PVector(0, 0, 0);
    } else {
      size = max(dist.x, dist.y);
      cm = new PVector(0, 0);
    }
  }


  boolean isleaf(int leafsize) {
    return end - start <= leafsize;
  }


  float min_node_dist_pb(PVector pos, PVector offset) {
    PVector dist = PVector.add(center, offset);
    dist = PVector.sub(pos, dist);
    return dist.mag() - rmax;
  }

  void display(int size_) {
    stroke(0, 0, 1);
    strokeWeight(1);
    noFill();
    if (dim) {
      float w = map(rhigh.x - rlow.x, 0, 1, 0, size_);
      float h = map(rhigh.y - rlow.y, 0, 1, 0, size_);
      float d = map(rhigh.z - rlow.z, 0, 1, 0, size_);

      float x = map(center.x, 0, 1, -size_/2, size_/2);
      float y = map(center.y, 0, 1, -size_/2, size_/2);
      float z = map(center.z, 0, 1, -size_/2, size_/2);

      pushMatrix();
      translate(x, y, z);
      box(w, h, d);
      popMatrix();
      if (!leaf) {
        right.display(size_);
        left.display(size_);
      }
    } else {
      float xl = map(rlow.x, 0, 1, -size_/2, size_/2);
      float yl = map(rlow.y, 0, 1, -size_/2, size_/2);
      float xh = map(rhigh.x, 0, 1, -size_/2, size_/2);    
      float yh = map(rhigh.y, 0, 1, -size_/2, size_/2);
      rectMode(CORNERS);
      rect(xl, yl, xh, yh);
      if (!leaf) {
        right.display(size_);
        left.display(size_);
      }
    }
  }
}
