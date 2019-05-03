#ifndef RAPT_POLYGON2D_H
#define RAPT_POLYGON2D_H

/** This is a class for representing 2-dimensional polygons in the plane.

\todo: implement area computation */

template<class RealType>
class rsPolygon2D
{

public:

  //friend class rsAffineTransform2D<RealType>;  // why?

  /** \name Construction/Destruction */

  /** Constructor. Creates an empty polygon (i.e. it does not yet have any vertices). You can add
  vertices using addVertex */
  rsPolygon2D()
  {

  }

  /** Creates a regular polygon with the given number of vertices.
  \todo comment the other parameters, remove the rotation parameter (it is redundant with the
  intial angle?) */
  rsPolygon2D(int numVertices, rsPoint2D<RealType> center = rsPoint2D<RealType>(0, 0),
    RealType horizontalRadius = 1, RealType verticalRadius = 1,
    RealType initialAngle = 0, RealType rotation = 0)
  {
    vertices.ensureAllocatedSize(numVertices);
    int i;
    rsPoint2D<RealType> v;

    double scaler = 2.0*PI / numVertices;
    for(i = 0; i < numVertices; i++)
    {
      double a = initialAngle + scaler * i; // angle
      v.x = horizontalRadius * cos(a);
      v.y = verticalRadius   * sin(a);
      addVertex(v);
    }
    // \todo: optimize this using sinCos and/or trigonometric recursion


    // apply rotation, if necessary
    if(rotation != 0)
    {
      // \todo create an rsAffineTransform2D implementing the desired rotation and apply it to
      // all  vertices
    }

    // translate all vertices by the desired center:
    for(i = 0; i < numVertices; i++)
      vertices[i] += center;
  }


  /** Constructor that creates a regular polygon with the given number of vertices, given radius
  (defined as distance of the vertices from the center) that is centered at the given center. The
  first vertex will be at a point that is a distance of "radius" away from the center at an angle
  of "rotation" (given in radians). For example, to create unit square, you would call it like
  ....
  \todo: compare this to the more general constructor above - actually this one here is redundant
  -> get rid of it
  */
  rsPolygon2D(int numVertices,
    RealType radius = 1,
    rsPoint2D<RealType> center = rsPoint2D<RealType>(0, 0),
    RealType rotation = 0)
  {
    vertices.ensureAllocatedSize(numVertices);
    int i;
    rsPoint2D<RealType> v;

    // create vertices at equally spaced angles around the origin, beginning with (radius, 0):
    v.x = radius;
    v.y = 0;
    addVertex(v);
    for(i = 1; i < numVertices; i++)
    {
      double a = 2.0 * i * PI / numVertices; // angle
      v.x = radius * cos(a);
      v.y = radius * sin(a);
      addVertex(v);
    }
    // \todo: optimize this using sinCos and/or trigonometric recursion


    // apply rotation, if necessary
    if(rotation != 0)
    {
      // \todo create an rsAffineTransform2D implementing the deired rotation and apply it to all
      // vertices
    }

    // translate all vertices by the desired center:
    for(i = 0; i < numVertices; i++)
      vertices[i] += center;
  }


  /** \name Setup */

  /** Clears all the vertices such that the polygon is "empty" again. */
  void clearVertices()
  {
    vertices.clear();
  }

  /** Allocates memory large enough to hold at least "numVerticesToReserve" vertices. You can
  call this function before successive calls to addVertex when you know the number of vertices to
  be added in advance. This optimizes memory allocation (i.e. avoid re-allocations when the
  vertex array grows). */
  void reserveVertexMemory(int numVerticesToReserve)
  {
    vertices.ensureAllocatedSize(numVerticesToReserve);
  }

  /** Adds a vertex to the polygon. */
  void addVertex(rsPoint2D<RealType> vertexToAdd)
  {
    vertices.appendElement(vertexToAdd);
  }

  //void reserveMemoryForMoreVertices(int numberToReserve)
  // ...to accelerate performance of addVertex


  /** \name Inquiry */

  /** Returns one of the vertices. The caller must be sure that the index is in the valid
  range. */
  inline rsPoint2D<RealType> getVertex(int index) const
  {
    return vertices.getElement(index);
  }

  /** Returns the number of vertices in the polygon. */
  inline int getNumVertices() const
  {
    return vertices.getNumElements();
  }

  /** Returns true when this polygon contains the given point, false otherwise. */
  //bool containsPoint(const rsPoint2D<RealType> &p) const
  //{
  //  return isPointInsideTriangle(p, *this);
  //}

  //rsDblRectangle getBoundingBox();


  /** \name Operators */

  /** Compares two polygons for equality. */
  bool operator==(const rsPolygon2D& p) const
  {
    return vertices == p.vertices;
  }

  /** Compares two polygons for inequality. */
  bool operator!=(const rsPolygon2D& p) const
  {
    return !(*this == p);
  }

  // intersectWith(const rsPolygon2D& p);
  //   returns a polygon that represents the intersection of "this" and "p"

  // makeRegularPolygon(rsPoint2D<RealType> center, RealType radius, int numVertices);
  //    -> factory function to create a regular polygon

protected:

  /** \name Data */

  std::vector<rsPoint2D<RealType>> vertices;

};

// a typedef'd explicit instantiation for coordinates of type double:
//typedef rsPolygon2D<double> rsDblPolygon2D;


#endif
