


class Cell{
  PVector bot;
  PVector top;
  float mass;
  int max_size;
  PVector center_of_volume;
  float radius;
  Cell parent;
  boolean isLeaf;
  float pivot;
  int split_Dim;
  boolean has_bary_weight;
  PVector center_of_mass;
  Cell child_bot;
  Cell child_top;
  ArrayList<Particle> particles;
  
  Cell(PVector low, PVector high, int _max_size, Cell _parent){
    bot = low;
    top = high;
    mass = 0;
    max_size = _max_size;
    center_of_volume = PVector.sub(high, low).div(2);
    radius = PVector.sub(high, center_of_volume).mag();
    parent = _parent;
    isLeaf = true;
    has_bary_weight = false;
  }
  
   ArrayList<Particle> get_particles(){
     ArrayList<Particle> all_particles = new ArrayList<>();
     ArrayList<Cell> cell_queue = new ArrayList<Cell>();
     cell_queue.add(this);
     while(cell_queue.size() > 0){
       Cell cell = cell_queue.remove(0);
       if(cell.isLeaf){
         all_particles.addAll(cell.particles);
       }
       else{
         cell_queue.add(cell.child_bot);
         cell_queue.add(cell.child_top);
       }
     }
     return all_particles;    
   }
   
   
   
   void split(){
     isLeaf = false;
     
     
   }
  
  
  
}
