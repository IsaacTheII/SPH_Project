import java.util.*;


class Particle {
  PVector r;
  PVector v;
  PVector v_pred;
  PVector a;
  float m;
  float radius;
  float e;
  float e_pred;
  float e_dot;
  float c;
  float rho;
  float h;
  ArrayList<Particle> n_clostest;


  Particle(PVector pos, PVector vel, float _m, float _radius, float _e, float _c) {
    this.r = pos;
    this.v = vel;
    this.v_pred = vel.copy();
    this.a = new PVector(0, 0, 0);
    this.m = _m;
    this.radius = _radius;
    this.e = _e;
    this.e_pred = _e;
    this.e_dot = 0.0;
    this.c = _c;
    this.rho = 0.0;
    this.h = 0.0;
    this.n_clostest = new ArrayList<Particle>();
  }

  void periodic_boundary() {
    if (this.r.x >= 1) {
      this.r.x -= 1;
    } else if (this.r.x < 0) {
      this.r.x+=1;
    }
    if (this.r.y >= 1) {
      this.r.y -= 1;
    } else if (this.r.y < 0) {
      this.r.y+=1;
    }
    if (this.r.z >= 1) {
      this.r.z -= 1;
    } else if (this.r.z < 0) {
      this.r.z+=1;
    }
  }

  float dist(Particle other) {
    return this.r.dist(other.r);
  }

  void update_pos(float dt) {
    this.r.add(PVector.mult(this.v, dt));
    this.periodic_boundary();
  }

  void update_vel(float dt) {
    this.v.add(PVector.mult(this.a, dt));
  }
}
