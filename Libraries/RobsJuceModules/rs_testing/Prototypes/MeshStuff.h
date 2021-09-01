
//=================================================================================================

// maybe rename to rsSurfaceMeshGenerator, use rsVector3D instead of 2D for the geometry to make 
// more sense
// maybe make other classes that can create more general patterns of connectivity

template<class T>
class rsMeshGenerator2D
{

public:


  //rsMeshGenerator2D() {}

  /** Enumeration of the available topologies. This determines, how the vertices are connected to
  other vertices in the mesh. Although the names are suggestive, the topology itself does not 
  imply a certain shape in 3D space...tbc... */
  enum class Topology
  {
    plane,        // edges of the parameter rectangle are not connected to anything
    cylinderV,    // vertical cyclinder: left edge is connected to right edge
    cylinderH,    // horizontal cyclinder: top edge is connected to bottom edge
    torus,        // connects left to right and top to bottom
    mobiusStrip,  // like cylinder but with right edge reversed
    kleinBottle   // like torus but with top and right edge reversed

    // cone         // top edge is connected to an additional tip vertex
    // doubleCone   // top and bottom edges are connected to additional tip vertices
    // closedCylinder // vertices of top and bottom edges are connected to other vertices on the 
                      // same edge (with offsets of Nu/2), forming a star when seen from above
  };
  // the doubleCone and closedCylinder topologies can also be used for a sphere - the actual 
  // interpretation as 3D shape is determined by the geometry, i.e. by the associated 3D mesh
  // i think, the kleinBottle is wrong: one edge must be connected in normal and the other in reverse 
  // mode

  /** Enumeration of the available geometries. This determines, how the (u,v)-coordinates of 
  vertices in the parameterMesh are mapped to (x,y,z)-coordinates in the spatialMesh. */
  /*
  enum class Geometry
  {
    plane,
    cylinder
  };
  */
  // or maybe let the user provide functions fx,fy,fz for the mapping - this is more flexible 
  // and/or maybe write another class rsMeshMapper or something that simplifies this
  


  // make a similar enum class for the geometry...maybe for the user, it would be more convenient
  // to just select a shape that determines both, topology and geometry


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setNumSamples(int Nx, int Ny)
  {
    this->Nx = Nx;
    this->Ny = Ny;
  }

  void setParameterRange(T x0, T x1, T y0, T y1)
  {
    this->x0 = x0;
    this->x1 = x1;
    this->y0 = y0;
    this->y1 = y1;
  }
  // maybe get rid or if not, rename to setRange

  void setTopology(Topology newTopology)
  {
    topology = newTopology;
  }

  /*
  void setVertexCoordinates(int i, int i, T x, T y)
  {
    int k = flatIndex(i, j);
    parameterMesh.setVertexData(k, rsVector2D<T>(x, y));
  }
  */
  // can be used to manually set the coordinates of vertex with given index pair


  // setTopology, setGeometry


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int flatIndex(int i, int j) const { return i  * Ny + j; }




  //-----------------------------------------------------------------------------------------------
  // \name Retrieval

  /** Should be called after setting up the desired mesh parameters. */
  void updateMeshes()
  {
    updateParameterMesh();
    updateSpatialMesh();
  }
  // get rid of this!

  /** Returns the parameter mesh in u,v space. This is what a PDE solver needs to pass to
  rsNumericDifferentiator::gradient2D to compute numerical approximations to the partial
  derivatives. They are taken with respect to our parameters u and v, which take the role of x
  and y in the notation of gradient2D. ...it's a bit confusing that there, the name u is used for
  the dependent variable: our parameter u here maps to the x there and the u there is the scalar
  field f(u,v) in a PDE...but these are the conventions: u is conventionally used as first
  parameter in surface theory and also as scalar field variable in PDE theory...tbc... */
  rsGraph<rsVector2D<T>, T> getParameterMesh() const { return parameterMesh; }
  // renamed to getMesh

  /** Returns the spatial mesh in (x,y,z)-space that corresponds to the parameter mesh in
  (u,v)-space. This is what a visualizer needs to display the results...tbc...  */
  rsGraph<rsVector3D<T>, T> getSpatialMesh()   const { return spatialMesh; }
  // the functionality of mapping a 2D mesh of (u,v) parameters should be moved to another class,
  // maybe rsMeshMapper_2D_to_3D
  // maybe they should return const references to the members






protected:

  //-----------------------------------------------------------------------------------------------
  // \name Internal

  void updateParameterMesh()
  {
    parameterMesh.clear();
    addParameterVertices();
    addConnections();
    updateParameterCoordinates();
  }
  // todo: do not create the mesh from scratch if not necessary - for example, if only the topology
  // has changed, we may keep the vertices and need only recompute the edges.

  void updateSpatialMesh()
  {

  }


  void addParameterVertices()
  {
    for(int i = 0; i < Nx; i++)
      for(int j = 0; j < Ny; j++)
        parameterMesh.addVertex(rsVector2D<T>(T(i), T(j)));
  }
  // -maybe addVertex should allow to pre-allocate memory for the edges - then we can pre-compute
  //  the number of vertices and pre-allocate - maybe each vertex could also pre-allocate space for
  //  the edges (use 4)


  /** Adds the connections between the vertices, i.e. the graph's edges - but we use the term 
  "edge" here to mean the edges of our parameter rectangle, so to avoid confusion, we call the
  graph-theoretic edges "connections" here. */
  void addConnections()
  {
    addCommonConnections();
    addTopologySpecificConnections();
  }

  /** Adds the connections that are always present, regardless of selected topology. */
  void addCommonConnections()
  {
    connectInner();             // rename to connectInnerStraight

    connectInnerDiagonal();
    // -make optional or maybe connectInner should do the switch 
    // -seems to lead to memory corruption - an assert gets triggered by rsImage - weird!

    connectTop();
    connectBottom();
    connectLeft();
    connectRight();
    connectCorners();
  }

  /** Adds the extra connections that depend on the selected topology. */
  void addTopologySpecificConnections()
  {
    using TP = Topology;
    switch(topology)  
    {
    case TP::cylinderV: {    connectLeftToRight();  } break;
    case TP::cylinderH: {    connectTopToBottom();  } break;
    case TP::torus: {        connectLeftToRight();
                             connectTopToBottom();  } break;
    case TP::mobiusStrip: {  connectLeftToRightReversed(); } break;
    case TP::kleinBottle: {  connectLeftToRightReversed();
                             connectTopToBottomReversed(); } break;
    default: { }   // topology is plane - nothing to do
    }
  }



  // maybe have functions: addInnerConnectionsDirect, addInnerConnectionsDiagonal, 
  // addTopConnections, addBottomConnections, addLeftConnections, addRightConnections, 
  // addTopLeftConnections ..or shorther: connectInner, connectTop, connectBottom

  /** Connects all inner vertices to their 4 direct neighbours with an edge. */
  void connectInner()
  {
    for(int i = 1; i < Nx-1; i++) {
      for(int j = 1; j < Ny-1; j++) {
        int k = flatIndex(i, j);
        parameterMesh.addEdge(k, west( i, j));
        parameterMesh.addEdge(k, east( i, j));
        parameterMesh.addEdge(k, north(i, j));
        parameterMesh.addEdge(k, south(i, j)); }}
  }
  // maybe make a version to connect to diagonal neighbors, too



  void connectInnerDiagonal()
  {
    for(int i = 1; i < Nx-1; i++) {
      for(int j = 2; j < Ny-1; j++) {
        int k = flatIndex(i, j);
        parameterMesh.addEdge(k, northEast(i, j));
        parameterMesh.addEdge(k, northWest(i, j));
        parameterMesh.addEdge(k, southEast(i, j));
        //int k2 = southWest(i, j);
        parameterMesh.addEdge(k, southWest(i, j)); 
        // This call creates memory corruption error later on, namely when rsImage::allocateMemory
        // is called. WTF?!
      }
    }

    // The j-loop starts at 2, because of weird memory corruptions in testTransportEquation when it 
    // starts at 1 as it actually should. I think, they come from the southWest-edge - for some 
    // value of i, with j==1, there must be something going wrong -> figure out and debug and write
    // a note about what it was. It's a kinda lucky circumstance that i found the location of the 
    // bug, when its effect occurs totally elsewhere. This should be documented for reference.
    // For debug:
    //int i = 2, j = 1;  // i = 1 or Nu-2 do not trigger it
    //int k = flatIndex(i, j);
    //parameterMesh.addEdge(k, southWest(i, j));

  }

  // maybe split into 2 functions for west and east diagonals


  /** Connects the vertices at the top, except for the (top-left and top-right) corners to their
  left, right and bottom neighbours. */
  void connectTop()
  {
    int j = Ny-1;
    for(int i = 1; i < Nx-1; i++) {
      int k = flatIndex(i, j);
      parameterMesh.addEdge(k, west( i, j));
      parameterMesh.addEdge(k, east( i, j));
      parameterMesh.addEdge(k, south(i, j)); }
    // todo: depending on topology, maybe have wrap-around-connections - but maybe all topology
    // dependent connections should be consolidated in one function
  }

  void connectBottom()
  {
    int j = 0;
    for(int i = 1; i < Nx-1; i++) {
      int k = flatIndex(i, j);
      parameterMesh.addEdge(k, west( i, j));
      parameterMesh.addEdge(k, east( i, j));
      parameterMesh.addEdge(k, north(i, j)); }
  }

  void connectLeft()
  {
    int i = 0;
    for(int j = 1; j < Ny-1; j++) {
      int k = flatIndex(i, j);
      parameterMesh.addEdge(k, north(i, j));
      parameterMesh.addEdge(k, south(i, j));
      parameterMesh.addEdge(k, east( i, j)); }
  }

  void connectRight()
  {
    int i = Nx-1;
    for(int j = 1; j < Ny-1; j++) {
      int k = flatIndex(i, j);
      parameterMesh.addEdge(k, north(i, j));
      parameterMesh.addEdge(k, south(i, j));
      parameterMesh.addEdge(k, west( i, j)); }
  }

  void connectCorners()
  {
    int i, j, k;

    i = 0; j = 0; k = flatIndex(i, j);        // bottom left
    parameterMesh.addEdge(k, east( i, j));
    parameterMesh.addEdge(k, north(i, j));

    i = Nx-1; j = 0;  k = flatIndex(i, j);    // bottom right
    parameterMesh.addEdge(k, west( i, j));
    parameterMesh.addEdge(k, north(i, j));

    i = 0; j = Ny-1;  k = flatIndex(i, j);    // top left
    parameterMesh.addEdge(k, east( i, j));
    parameterMesh.addEdge(k, south(i, j));

    i = Nx-1; j = Ny-1; k = flatIndex(i, j);  // top right
    parameterMesh.addEdge(k, west( i, j));
    parameterMesh.addEdge(k, south(i, j));
  }

  void connectLeftToRight()
  {
    for(int j = 0; j < Ny; j++) {
      int k1 = flatIndex(0,    j);
      int k2 = flatIndex(Nx-1, j);
      parameterMesh.addEdge(k1, k2, true); }
  }

  void connectLeftToRightReversed()
  {
    for(int j = 0; j < Ny; j++) {
      int k1 = flatIndex(0,    j     );
      int k2 = flatIndex(Nx-1, Ny-1-j);
      parameterMesh.addEdge(k1, k2, true); }
  }

  void connectTopToBottom()
  {
    for(int i = 0; i < Nx; i++) {
      int k1 = flatIndex(i, 0   );
      int k2 = flatIndex(i, Ny-1);
      parameterMesh.addEdge(k1, k2, true); }
  }

  void connectTopToBottomReversed()
  {
    for(int i = 0; i < Nx; i++) {
      int k1 = flatIndex(i,      0   );
      int k2 = flatIndex(Nx-1-i, Ny-1);
      parameterMesh.addEdge(k1, k2, true); }
  }


  void updateParameterCoordinates()
  {
    for(int i = 0; i < Nx; i++) {
      for(int j = 0; j < Ny; j++) {
        int k = flatIndex(i, j);
        T x = rsLinToLin(T(i), T(0), T(Nx-1), x0, x1);
        T y = rsLinToLin(T(j), T(0), T(Ny-1), y0, y1);
        parameterMesh.setVertexData(k, rsVector2D<T>(x, y)); }}
  }

  // Functions to compute vertex indices for neighbors:
  int east(     int i, int j) const { return flatIndex(i+1, j  ); }
  int west(     int i, int j) const { return flatIndex(i-1, j  ); }
  int north(    int i, int j) const { return flatIndex(i,   j+1); }
  int south(    int i, int j) const { return flatIndex(i,   j-1); }

  int northEast(int i, int j) const { return flatIndex(i+1, j+1); }
  int northWest(int i, int j) const { return flatIndex(i-1, j+1); }
  int southEast(int i, int j) const { return flatIndex(i+1, j-1); }

  int southWest(int i, int j) const { return flatIndex(i-1, j-1); }
  // this seems to cause weird problems in connectInnerDiagonal



  //-----------------------------------------------------------------------------------------------
  // \name Data

  int Nx = 0;   // number of vertices along x-coordinate
  int Ny = 0;   // number of vertices along y-coordinate

  T x0 = T(0), x1 = T(1);    // lower and upper limit for x-coordinate
  T y0 = T(0), y1 = T(1);    // lower and upper limit for y-coordinate
  // use x0,etc
  // maybe use 0,1 always - always operate on normalized coordinates..or maybe not?

  Topology topology = Topology::plane;

  rsGraph<rsVector2D<T>, T> parameterMesh; // rename to mesh
  rsGraph<rsVector3D<T>, T> spatialMesh;   // move elsewhere - don't make a god class!

};
// move implementations out of the class

// -todo: interpret the (x,y) values in the 2D grid and (u,v)-parameters that are later mapped to
//  (x,y,z)-coordinates, according to some geometry settings
// -let user select topology and geometry...but maybe the geometric aspects should be done by a 
//  separate class - maybe we need to keep two meshes: one in 2D parameter space and one in 3D 
//  actual space (maybe the latter can be also 4D for the Klein bottle)
// -the PDE solver uses the 2D mesh in parameter space, but the interpretation and visualization is 
//  done by means of the 3D mesh - also, the weights (distances) for the 2D parameter mesh have to
//  be computed by help of the 3D mesh in user space - but optionally, we may also use the default
//  weighting
// -geometries: plane, torus, cylinder, mobiusStrip, kleinBottle, cone, doubleCone, disc, sphere
// -use should set up stuff via setters and then retrieve 2 meshes via 
//  getParameterMesh, getSpatialMesh - the 1st is 2D and used in  the PDE solver, the 2nd is for
//  visualization
// -a double-cone or spherical topology can be achieved in 1 of two ways:
//  (1) close a cylinder on top and bottom by connecting top and bottom row vertices like: 
//      bottom: v(i,0) to v((i+Nx/2)%Nx, 0), top: v(i,Ny-1) to v((i+Nx/2)%Nx, Ny-1)   i=0..Nx-1
//      this makes only sense when Nx is even (right? what happens, when Nx is odd? try it!)
//  (2) add two additional vertice for the two cone cusps (top, bottom) and connecting top and 
//      bottom rows to these cusp vertices
//  -for a simple cone, do this only for the top row, a simple cone can also degenerate to a disc
//   when the height of the cone is zero (i.e. z-values are all equal)
// -for best flexibility, the user 3 functions fx(u,v), fy(u,v), fz(u,v) to compute coordinates

//=================================================================================================

template<class T>
void randomizeVertexPositions(rsGraph<rsVector2D<T>, T>& mesh, T dx, T dy, 
  int minNumNeighbors = 0, int seed = 0);