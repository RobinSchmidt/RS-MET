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

template<class T>
class rsGraphWithVertexData
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Adds a new vertex with the given associated data. */
  void addVertex(const T& data) { vertices.push_back(Vertex(data)); }

  /** Connects vertex i to vertex j by an edge. If the optional boolean parameter "bothWays" is 
  true, it also adds the edge from j to i. */
  void addEdge(int i, int j, bool bothWays = false)  
  { 
    vertices[i].neighbors.push_back(j);
    if(bothWays)
      vertices[j].neighbors.push_back(i);
  }

  /** Modifies the data associated with vertex i. */
  void setVertexData(int i, const T& data) { vertices[i].data = data; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of vertices. */
  int getNumVertices() const { return (int) vertices.size(); }

  /** Returns the number of neighbors that are adjacent to vertex i. */
  int getNumNeighbors(int i) const { return (int) vertices[i].neighbors.size(); }

  /** Returns a const reference to the array of neighbors of vertex i. */
  const std::vector<int>& getNeighbors(int i) const { return vertices[i].neighbors; }

  /** Returns a const reference to the data that is associated with vertex i. */
  const T& getVertexData(int i) const { return vertices[i].data; }

  //-----------------------------------------------------------------------------------------------
  // \name Data

protected:

  /** Structure to hold one vertex containing its associated data and a list of indices of 
  vertices that are connected to this vertex. */
  struct Vertex
  {
    Vertex(const T& newData) : data(newData) {}
    T data;                       // data associated with the vertex
    std::vector<int> neighbors;   // array of vertex indices that are neighbours of this vertex
  };

  std::vector<Vertex> vertices;

};

// ToDo:
// -implement functions like isConnected(int i, int j), containsDuplicateEdges(), 
//  containsDuplicateVertices
// -maybe instead of sdt::vector, we could use rsSortedSet
// -instead of having each vertex maintain a list of adjacent vertices, we could have an explicit
//  array of edges - which data-structure is better may depend on the situation and maybe it makes
//  sense to have both variants
// -maybe allow (optionally) data to be associated with each edge

#endif