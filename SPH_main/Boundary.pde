class Boundary { //<>//

  PVector startPoint;
  PVector endPoint;
  PVector boundaryVector;
  PVector normal;
  PVector crossingPoint;  
  PVector center;
  float det;
  boolean isStationary = true;

  Boundary(PVector startPoint_, PVector endPoint_) {
    startPoint = startPoint_;
    endPoint = endPoint_;
    boundaryVector = PVector.sub(endPoint, startPoint); // pointless?
    //normal = new PVector(startPoint.x - endPoint.x, startPoint.y - endPoint.y);
    //normal = new PVector(startPoint.y - endPoint.y, startPoint.x - endPoint.x);
    normal = new PVector(startPoint.y - endPoint.y, endPoint.x - startPoint.x);
    normal.div(normal.mag());
    crossingPoint = new PVector();
    center = PVector.div(PVector.add(startPoint, endPoint), 2);
  }

  boolean hasCrossed(Particle particle) {

    float t_nominator = (particle.pos.x - startPoint.x) * (startPoint.y - endPoint.y) - (particle.pos.y - startPoint.y)*(startPoint.x - endPoint.x);
    float t_denominator = (particle.pos.x - particle.temp_pos.x)*(startPoint.y - endPoint.y) - (particle.pos.y - particle.temp_pos.y)*(startPoint.x - endPoint.x);
    float t = t_nominator/t_denominator;

    float u_nominator =  (particle.pos.x - particle.temp_pos.x)*(particle.pos.y - startPoint.y) - (particle.pos.y - particle.temp_pos.y)*(particle.pos.x - startPoint.x);
    float u = -u_nominator/t_denominator;

    if (t_denominator == 0.0) {
      return false; // lines are parallel
    } else if (t <= 1.0 && t >= 0.0 && u <= 1.0 && u >= 0.0) {

      crossingPoint.set(particle.pos.x + t*(particle.temp_pos.x - particle.pos.x), particle.pos.y + t*(particle.temp_pos.y - particle.pos.y));
      return true;
    } else { 
      return false; // lines don't cross
    }
  }

  float distanceToBoundary(Particle particle) {

    PVector P1 = particle.pos;
    PVector P2 = PVector.sub(P1, normal);

    float t_nominator = (P1.x - startPoint.x) * (startPoint.y - endPoint.y) - (P1.y - startPoint.y)*(startPoint.x - endPoint.x);
    float t_denominator = (P1.x - P2.x)*(startPoint.y - endPoint.y) - (P1.y - P2.y)*(startPoint.x - endPoint.x);
    float t = t_nominator/t_denominator;

    float u_nominator =  (P1.x - P2.x)*(particle.pos.y - startPoint.y) - (P1.y - P2.y)*(P1.x - startPoint.x);
    float u = -u_nominator/t_denominator;

    if (t_denominator == 0.0) {
      return MAX_FLOAT; // lines are parallel
    } else if (t <= 1.0 && t >= 0.0 && u <= 1 && u >= 0.0) {
      crossingPoint.set(particle.pos.x + t*(particle.temp_pos.x - particle.pos.x), particle.pos.y + t*(particle.temp_pos.y - particle.pos.y));
      return PVector.sub(crossingPoint, P1).mag();
    } else { 
      return MAX_FLOAT; // lines don't cross
    }
  }

  void reflectAtBoundary(Particle particle) {
    // Reflect particle at the boundary
    PVector overShoot = PVector.sub(particle.temp_pos, crossingPoint);
    float factorNormal = overShoot.dot(normal);
    //particle.pos = PVector.add(particle.temp_pos, PVector.mult(normal, 2*factorNormal));
    //particle.vel = PVector.sub(particle.pos, crossingPoint).setMag(particle.vel.mag());
    particle.temp_pos.add(PVector.mult(normal, -2*factorNormal));
    particle.vel = PVector.sub(particle.temp_pos, crossingPoint).setMag(particle.vel.mag());
  }

  void checkBoundary(Particle particle) {

    if (hasCrossed(particle)) {
      reflectAtBoundary(particle);
    }
  }

  ArrayList<Particle> reflectParticles(ArrayList<ParticleTuple> neighbours) {
    ArrayList<Particle> reflected = new ArrayList<Particle>();
    for (ParticleTuple tup : neighbours) {
      Particle p = tup.p;

      //dist = PVector.sub(pos, PVector.add(particles.get(i).pos, temp.offset));

      PVector pos = PVector.sub(p.pos, tup.offset); // changed from "add" to "sub" -> seems to work (symmetrical behaviour)
      float nominator = (endPoint.y - startPoint.y) * pos.x - (endPoint.x - startPoint.x) * pos.y + endPoint.x * startPoint.y - endPoint.y * startPoint.x;
      float denominator = pow(pow(endPoint.y - startPoint.y, 2) + pow(endPoint.x - startPoint.x, 2), 0.5);
      float d = abs(nominator)/ denominator;
      PVector newPos = PVector.add(pos, PVector.mult(normal, -2 * d));

      PVector tempPos = PVector.add(pos, p.vel);
      nominator = (endPoint.y - startPoint.y) * tempPos.x - (endPoint.x - startPoint.x) * tempPos.y + endPoint.x * startPoint.y - endPoint.y * startPoint.x;
      denominator = pow(pow(endPoint.y - startPoint.y, 2) + pow(endPoint.x - startPoint.x, 2), 0.5);
      d = abs(nominator)/ denominator;
      tempPos = PVector.add(tempPos, PVector.mult(normal, -2 * d));
      PVector newVel = PVector.sub(tempPos, newPos);
      Particle particle = new Particle(newPos, newVel, new PVector(0, 0), p.m, p.e);
      particle.isReflectedParticle = true;
      reflected.add(particle);
    }
    return reflected;
  }

  void rotateBoundary(float alpha) { // rotate a (line segment) boundary around its center point (in 2D)

    float s = sin(alpha);
    float c = cos(alpha);

    // translate point back to origin:
    startPoint.x -= center.x;
    startPoint.y -= center.y;
    endPoint.x -= center.x;
    endPoint.y -= center.y;

    // rotate point around origin
    float newStartPointx = startPoint.x * c - startPoint.y * s;
    float newStartPointy = startPoint.x * s + startPoint.y * c;
    float newEndPointx = endPoint.x * c - endPoint.y * s;
    float newEndPointy = endPoint.x * s + endPoint.y * c;

    // translate point back:
    startPoint.x = newStartPointx + center.x;
    startPoint.y = newStartPointy + center.y;
    endPoint.x = newEndPointx + center.x;
    endPoint.y = newEndPointy + center.y;
    
    // update normal vector
    normal = new PVector(startPoint.y - endPoint.y, endPoint.x - startPoint.x);
    normal.div(normal.mag());
    center = PVector.div(PVector.add(startPoint, endPoint), 2);
  }

  void drawBoundary(int size) {
    float x1 = map(startPoint.x, 0, 1, -size/2, size/2);
    float y1 = map(startPoint.y, 0, 1, -size/2, size/2);
    float x2 = map(endPoint.x, 0, 1, -size/2, size/2);
    float y2 = map(endPoint.y, 0, 1, -size/2, size/2);

    // Draw boundaries
    strokeWeight(5);
    line(x1, y1, x2, y2);
  }
}
