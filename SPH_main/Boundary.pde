class Boundary { //<>//

  PVector startPoint;
  PVector endPoint;
  PVector boundaryVector;
  PVector normal;
  PVector crossingPoint;  
  float det;


  Boundary(PVector startPoint_, PVector endPoint_) {
    startPoint = startPoint_;
    endPoint = endPoint_;
    boundaryVector = PVector.sub(endPoint, startPoint); // pointless?
    //normal = new PVector(startPoint.x - endPoint.x, startPoint.y - endPoint.y);
    //normal = new PVector(startPoint.y - endPoint.y, startPoint.x - endPoint.x);
    normal = new PVector(endPoint.y - startPoint.y, startPoint.x - endPoint.x);
    normal.mult(1/normal.mag());
    crossingPoint = new PVector();
  }

  boolean hasCrossed(Particle particle) {

    float t_nominator = (particle.pos.x - startPoint.x) * (startPoint.y - endPoint.y) - (particle.pos.y - startPoint.y)*(startPoint.x - endPoint.x);
    float t_denominator = (particle.pos.x - particle.temp_pos.x)*(startPoint.y - endPoint.y) - (particle.pos.y - particle.temp_pos.y)*(startPoint.x - endPoint.x);

    if (t_denominator == 0.0) {
      return false; // lines are parallel
    }

    float t = t_nominator/t_denominator;

    if (t < 0.0 || t > 1.0) {
      return false; // lines don't cross
    }

    crossingPoint.set(particle.pos.x + t*(particle.temp_pos.x - particle.pos.x), particle.pos.y + t*(particle.temp_pos.y - particle.pos.y));
    //println("Boundaries: ", startPoint, endPoint);
    //println("Has crossed at: ", crossingPoint);
    //println("Particle pos / Particle temp_pos", particle.pos, particle.temp_pos);
    //println("Normal to boundary: ", normal);
    //exit();

    return true;
  }

  void reflectAtBoundary(Particle particle) {
    // Reflect particle at the boundary
    PVector overShoot = PVector.sub(particle.temp_pos, crossingPoint);
    float factorNormal = -overShoot.dot(normal);
    //particle.pos = PVector.add(particle.temp_pos, PVector.mult(normal, 2*factorNormal));
    //particle.vel = PVector.sub(particle.pos, crossingPoint).setMag(particle.vel.mag());
    particle.temp_pos.add(PVector.mult(normal, 2*factorNormal));
    particle.vel = PVector.sub(particle.temp_pos, crossingPoint).setMag(particle.vel.mag());
  }

  void checkBoundary(Particle particle) {

    if (hasCrossed(particle)) {
      reflectAtBoundary(particle);
    }
  }

  void drawBoundary(int size) {
    float x1 = map(startPoint.x, 0, 1, -size/2, size/2);
    float y1 = map(startPoint.y, 0, 1, -size/2, size/2);
    float x2 = map(endPoint.x, 0, 1, -size/2, size/2);
    float y2 = map(endPoint.y, 0, 1, -size/2, size/2);

    // Draw boundaries
    strokeWeight(10);
    line(x1, y1, x2, y2);
  }
}
