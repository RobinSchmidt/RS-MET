#ifndef RS_RSAGG_H
#define RS_RSAGG_H

#include "../../Math/Geometry/Polygon2D.h"
#include "../Misc/Image.h"
using namespace RSLib;

namespace rsAGG
{

  /** Renders a filled polygon with given color into the given image using 
  agg::rasterizer_scanline_aa and agg::renderer_aa. */
  void drawFilledPolygon(const rsPolygon2D<double> &polygon, rsImageRegionRGBA &imageRegion, 
    const rsColorRGBA &color);

  void drawOutlinedPolygon(const rsPolygon2D<double> &polygon, rsImageRegionRGBA &imageRegion, 
    const rsColorRGBA &color, double strokeWidth);
    // add parameters: strokeStyle (solid, dashed, dotted etc.), joinStyle, capStyle
    // or maybe just pass a const reference to the rsGraphicsRenderer2DState object

  /** Creates a polygon that represents the outline of some input polygon according to the 
  selected line-thickness using agg::vcgen...? . */
  //rsPolygon2D<double> generatePolygonOutline(const rsPolygon2D<double> inputPolygon);

}

#endif
