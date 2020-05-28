class ParticleTuple {
  Particle p;
  float dist;
  PVector offset;
  ParticleTuple(Particle p_, float dist_, PVector offset_) {
    p = p_;
    dist = dist_;
    offset = offset_;
  }
}

class NodeTuple {
  Node node;
  float dist;
  PVector offset;
  NodeTuple(Node node_, float dist_, PVector offset_) {
    node = node_;
    dist = dist_;
    offset = offset_;
  }
}
