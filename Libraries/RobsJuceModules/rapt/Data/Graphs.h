#ifndef RAPT_GRAPHS_H
#define RAPT_GRAPHS_H

/** This file contains data structures to represent graphs in the sense of sets of vertices that
may be connected by edges. */


//=================================================================================================

/** Class to represent a graph that has data associated with each vertex but not with the edges. An
edge between any pair of vertices can just exist or don't exist. We keep the vertices in a 
std::vector, so they can be addressed by their indices and each vertex contains its data and a 
std::vector<int> with the indices of its connected ("neighbor") vertices, so this is basically an 
adjacency list representation. An example use is for irregular meshes of vertices for solving 
partial differential equations - in this case, the data type for the vertices could be 
rsVector2D<float> or similar. @see rsNumericDifferentiator::gradient2D for an example. */

template<class TVtx, class TEdg>
class rsGraph
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Data

  /** Structure to represent an edge between two vertices. It stores the index of the target 
  vertex, the source vertex is given implicitly by the fact that edges are stored as adjacency 
  lists inside the vertices themselves, so the source vertex is always the vertex that owns this 
  edge. It's probably easiest to visualize the vertices storing their outgoing edges, but this 
  interpretation as an edge being "outgoing" is up to client code - if it wants to, it may also 
  interpret the list of stored edges in each vertex as incoming. Each edge may also store 
  associated data, such as a weight or whatever - if you don't need any edge data, pass 
  (for example) rsEmptyType for the template parameter TEdg. */
  struct Edge
  {
    Edge(int newTarget, const TEdg& newData) : target(newTarget), data(newData) {}
    int target;   // index of target vertex of this edge
    TEdg data;    // data associated with this edge
  };

  /** Structure to hold one vertex containing its associated data and a list of edges emanating 
  from it. You can interpret the edges as outgoing or incoming - whatever is most convenient for 
  your particular use case. If you don't want any data stored at the vertices, you may pass (for 
  example) rsEmptyType for the template parameter TVtx. */
  struct Vertex  // maybe make it a class
  {
    Vertex(const TVtx& newData) : data(newData) {}
    TVtx data;                    // data associated with the vertex
    std::vector<Edge> edges;
    // new

    int getNumEdges()              const { return (int) edges.size(); }
    const TEdg& getEdgeData(int j) const { return edges[j].data;      }


    std::vector<int> neighbors;   // array of vertex indices that are neighbours of this vertex
    // old - get rid
  };


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Adds a new vertex with the given associated data. */
  void addVertex(const TVtx& data = TVtx()) { vertices.push_back(Vertex(data)); }
  // O(?)   numEdges? in the worst case, each vertex must be copied - together with its adjacency 
  // list

  /** Modifies the data stored at vertex i. */
  void setVertexData(int i, const TVtx& data) { vertices[i].data = data; }
  // O(1)


  void addEdge(int i, int j, const TEdg& data = TEdg(1), bool bothWays = false)
  {
    vertices[i].edges.push_back(j, data);
    if(bothWays)
      vertices[j].edges.push_back(i, data);
  }
  // O(vertices[i].numEdges + vertices[j].numEdges)



  /** Connects vertex i to vertex j by an edge. If the optional boolean parameter "bothWays" is 
  true, it also adds the edge from j to i. */
  void addEdgeOld(int i, int j, bool bothWays = false)
  { 
    vertices[i].neighbors.push_back(j);
    if(bothWays)
      vertices[j].neighbors.push_back(i);
  }
  // obsolete soon



  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of vertices. */
  int getNumVertices() const { return (int) vertices.size(); }
  // O(1)

  /** Returns a const reference to the data that is stored in vertex i. */
  const TVtx& getVertexData(int i) const { return vertices[i].data; }
  // O(1)

  /** Returns the number of edges emanating from vertex i. */
  int getNumEdges(int i) const { return (int) vertices[i].egdes.size(); }
  // O(1)

  /** Returns the number of edges in the whole graph. Note that if there's an edge from vertex i to
  vertex j and also one from vertex j to vertex i (the same edge backwards), both are counted. So
  if you use this class for an undirected graph by adding symmetric edges, you may want to divide 
  by 2. */
  int getNumEdges() const
  {
    int n = 0;
    for(int i = 0; i < getNumVertices(); i++)
      n += getNumEdges(i);
    return n;
  }
  // O(numVertices)

  /** Returns a const reference to the array of edges emanating from vertex i. */
  const std::vector<Edge>& getEdges(int i) const { return vertices[i].edges; }
  // O(1)

  /** Returns a const reference to the data that is stored in the j-th edge emanating from vertex 
  i. Note that j is (in general) not the index of the target vertex - it's the index at which the 
  edge occurs in the adjacency list of vertex i. */
  const TEdg& getEdgeData(int i, int j) const { return vertices[i].getEdgeData(j); }
  // O(1)




  // old - get rid:
  /** Returns the number of neighbors that are adjacent to vertex i. */
  int getNumNeighbors(int i) const { return (int) vertices[i].neighbors.size(); }
  // rename to getNumEdgesFrom

  /** Returns a const reference to the array of neighbors of vertex i. */
  const std::vector<int>& getNeighbors(int i) const { return vertices[i].neighbors; }
  // rename to getEdgesFrom






protected:

  std::vector<Vertex> vertices;

};

// ToDo:
// -implement functions like isConnected(int i, int j), containsDuplicateEdges(), 
//  containsDuplicateVertices, isSymmetric (i.e. non-directed)
// -maybe instead of sdt::vector, we could use rsSortedSet
// -instead of having each vertex maintain a list of adjacent vertices, we could have an explicit
//  array of edges - which data-structure is better may depend on the situation and maybe it makes
//  sense to have both variants
// -maybe allow (optionally) data to be associated with each edge
// -if we want a graph without vertex- or edge data, we can pass an empty struct for the respective
//  template parameter

#endif