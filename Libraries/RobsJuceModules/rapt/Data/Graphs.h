#ifndef RAPT_GRAPHS_H
#define RAPT_GRAPHS_H

/** This file contains data structures to represent graphs in the sense of sets of vertices that
may be connected by edges. */

//=================================================================================================

/** Class to represent a graph that may have data stored at the vertices and/or the edges. What 
kind of data is stored is determined by the template parameters TVtx and TEdg respectively. If you 
don't want any data stored at vertices and/or edges, you can pass an empty struct for the 
respective type (such as rsEmptyType) as dummy data type to satisfy the compiler without creating 
storage overhead.

The graph is represented as a std::vector of vertices (of type rsGraph::Vertex). Each vertex holds 
its data (if any) and a std::vector of edges (of type rsGraph::Edge), so it's an adjacency list 
based representation. The Edge objects themselves contain a field for the target vertex index, i.e.
the index of the vertex (in our std::vector of vertices), this edge points to and a data field 
(which can, for example, contain a weight, but may also be empty). Whether you want to interpret 
these "emanating" edges for each vertex as outgoing or incoming is up to you - the class itself is 
oblivious of such an interpretation, even though the terminology of "target", "points to" etc. may 
suggest an outgoing interpretation. The order of the vertices in the graph and of the emanating 
edges of each vertex is determined by the order they are added to the graph by client code. The 
class itself enforces no particular order and client code should better not assume any.

An example use is for irregular meshes of vertices for solving partial differential equations. In 
this case, the data type for the vertices TVtx could be rsVector2D<float> or similar (todo: use the 
edge data type TEdg for the weights later). ...tbc... @see rsNumericDifferentiator::gradient2D for 
an example. */

template<class TVtx, class TEdg> // types for the data stored at the vertices and edges
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
  class Edge
  {
  public:
    Edge(int newTarget, const TEdg& newData) : target(newTarget), data(newData) {}

    void setData(const TEdg& newData) { data = newData; }

    int         getTarget() const { return target; }
    const TEdg& getData()   const { return data;   }

  protected:
    int target;   // index of target vertex of this edge
    TEdg data;    // data stored at this edge
  };

  /** Class to hold one vertex containing its associated data and a list of edges emanating from 
  it. You can interpret the edges as outgoing or incoming - whatever is most convenient for your 
  particular use case. If you don't want any data stored at the vertices, you may pass (for 
  example) rsEmptyType for the template parameter TVtx. */
  class Vertex
  {
  public:
    Vertex(const TVtx& newData) : data(newData) {}

    void setData(const TVtx& newData)         { data = newData;         }
    void addEdge(const Edge& edge)            { edges.push_back(edge);  }
    void setEdgeData(int j, const TEdg& data) { edges[j].setData(data); }
    void removeEdge(int k)                    { rsRemove(edges, k);     }

    const TVtx& getData()            const { return data;                 }
    int         getNumEdges()        const { return (int) edges.size();   }
    int         getEdgeTarget(int j) const { return edges[j].getTarget(); }
    const TEdg& getEdgeData(int j)   const { return edges[j].getData();   }
    bool hasEdgeTo(int j) const 
    { 
      for(int k = 0; k < (int)edges.size(); k++) {
        if(edges[k].getTarget() == j)
          return true; }
      return false;
    } // needs test

  protected:
    TVtx data;               // data stored at the vertex
    std::vector<Edge> edges; // edges emanating from this vertex

    friend class rsGraph<TVtx, TEdg>;
  };

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Adds a new vertex with the given associated data. */
  void addVertex(const TVtx& data = TVtx()) { vertices.push_back(Vertex(data)); }
  // O(?)   numEdges? in the worst case, each vertex must be copied - together with its adjacency 
  // list

  /** Modifies the data stored at vertex i. */
  void setVertexData(int i, const TVtx& data) { vertices[i].setData(data); }
  // O(1)


  void addEdge(int i, int j) { addEdge(i, j, TEdg(1), false); }

  void addEdge(int i, int j, const TEdg& data, bool bothWays = false)
  {
    rsAssert(isVertexIndexValid(i) && isVertexIndexValid(j), "Invalid vertex index");
    vertices[i].addEdge(Edge(j, data));
    if(bothWays)
      vertices[j].addEdge(Edge(i, data));
  }
  // O(vertices[i].numEdges + vertices[j].numEdges)

  /** Convenience function to add an edge with a default value of 1, possibly symmetrically. */
  void addEdge(int i, int j, bool bothWays) { addEdge(i, j, TEdg(1), bothWays); }

  /** Removes an edge from vertex i to vertex j, if such an edge exists and returns true, iff 
  there was such an edge. If there is no edge from i to j, it does nothing and returns false.  */
  //bool removeEdgeAt(int i, int j) {}

  /** Removes the k-th edge emanating from vertex i. Note that this is *not* in general an edge 
  between vertex i and vertex k. If you want to remove and edge by giving the idices of the 
  vertices it connects, use removeEdge */
  void removeEdgeAtIndex(int i, int k)
  {
    rsAssert(isVertexIndexValid(i) && vertices[i].getNumEdges() > k);
    vertices[i].removeEdge(k);
  }


  void setEdgeData(int i, int k, const TEdg& data) 
  { 
    rsAssert(isVertexIndexValid(i) && vertices[i].getNumEdges() > k);
    vertices[i].setEdgeData(k, data); 
  }
  // rename to setEdgeDataByIndex

  /** Clears the edges of all vertices. */
  void clearEdges()
  {
    for(size_t i = 0; i < vertices.size(); i++)
      vertices[i].edges.clear();
  }

  /** Clears the array of vertices. Resets the graph into its initial, empty state. */
  void clear() { vertices.clear(); }


  // todo: setEdgeData(int i, int j, const TEdg& data), removeEdge(i, j), removeVertex(i)

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns true, iff i is a valid vertex index, i.e. an index within the proper range. */
  bool isVertexIndexValid(int i)
  {
    return i >= 0 && i < (int)vertices.size();
  }

  /** Returns the number of vertices. */
  int getNumVertices() const { return (int) vertices.size(); }
  // O(1)

  /** Returns a const reference to the data that is stored in vertex i. */
  const TVtx& getVertexData(int i) const { return vertices[i].getData(); }
  // O(1)

  /** Returns the number of edges emanating from vertex i. */
  int getNumEdges(int i) const { return vertices[i].getNumEdges(); }
  // O(1)

  /** Returns the number of edges in the whole graph. Note that if there's an edge from vertex i to
  vertex j and also one from vertex j to vertex i (the same edge backwards), both are counted. So
  if you use this class for representing undirected graph by adding symmetric edges, you may want 
  to divide by 2. */
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
  const TEdg& getEdgeData(int i, int k) const { return vertices[i].getEdgeData(k); }
  // O(1)
  // rename to getEdgeDataByIndex

  /** Returns the index of the target vertex of the j-th edge emanating from vertex i. */
  int getEdgeTarget(int i, int k) const { return vertices[i].getEdgeTarget(k); }
  // rename to getEdgeTargetByIndex

  /** Returns true, iff for each edge from vertex i to j there also exists an edge from j to i. In 
  this case, the graph can be seen as undirected. */
  bool isSymmetric() const
  {
    for(int i = 0; i < (int) vertices.size(); i++) {
      for(int j = 0; j < (int) vertices[i].edges.size(); j++) {
        int k = getEdgeTarget(i, j);
        if(!vertices[k].hasEdgeTo(i))
          return false;   }}
    return true;
  }
  // needs test
  // maybe reverse roles of j and k - we should generally use i,j for vertex indices and k for the
  // vertex-index in the edge-list of vertex i
  // there should not only exist the backward edge, but it should also have the same data!


protected:

  std::vector<Vertex> vertices;

};

// ToDo:
// -implement functions like isConnected(int i, int j), containsDuplicateEdges(), 
//  containsDuplicateVertices, isSymmetric (i.e. non-directed), etc.
// -maybe instead of std::vector, we could use rsSortedSet - but maybe that should be deferred to
//  a different implementation, with similar interface but different goals (like faster 
//  determination whether two vertices are connected)
//  ...but maybe that could be realized by a 3rd (and 4th) template parameter which for the storage
//  container(s) which defaults to std::vector (maybe via a partial specialization). this could 
//  then be replaced by something like rsSortedSet, if needed. It could be desirable to keep the 
//  vertices and/or edges within each vertex sorted - care has to be taken: when reordering 
//  vertices, edge targets must also be updated
// -instead of having each vertex maintain a list of adjacent vertices, we could have an explicit
//  array of edges - which data-structure is better may depend on the situation and maybe it makes
//  sense to have both variants
// -maybe numEdges could be cached - but maybe in a subclass
// -maybe make a class rsGraphAlgorithms that contain static functions that manipulate graphs, the
//  rsGraph class itself should provide only basic barebone functionality (if the algorithms were 
//  added as methods, the class may grow unreasonably large - the number of such algorithms is too 
//  vast)

#endif