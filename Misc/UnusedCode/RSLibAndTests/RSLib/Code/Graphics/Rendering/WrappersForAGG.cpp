#include "AGG_2.4/include/agg_basics.h"
#include "AGG_2.4/include/agg_rendering_buffer.h"
#include "AGG_2.4/include/agg_rasterizer_scanline_aa.h"
//#include "AGG_2.4/include/agg_rasterizer_outline.h"   // rasterizes outlined shapes?
#include "AGG_2.4/include/agg_scanline_p.h"             // packed scanlines
#include "AGG_2.4/include/agg_scanline_u.h"             // unpacked scanlines
//#include "AGG_2.4/include/agg_scanline_bin.h"         // binary scanlines (pixel on/off)
#include "AGG_2.4/include/agg_renderer_scanline.h"
//#include "AGG_2.4/include/agg_renderer_primitives.h"  // renders primitives efficiently?
#include "AGG_2.4/include/agg_conv_stroke.h"            // for rendering outlined polygons
//#include "AGG_2.4/include/agg_conv_dash.h"            // for rendering dashed outlines
//#include "AGG_2.4/include/agg_conv_curve.h"
//#include "AGG_2.4/include/agg_conv_contour.h"
//#include "AGG_2.4/include/agg_conv_smooth_poly1.h"
//#include "AGG_2.4/include/agg_conv_marker.h"
//#include "AGG_2.4/include/agg_arrowhead.h"
//#include "AGG_2.4/include/agg_vcgen_markers_term.h"
#include "AGG_2.4/include/agg_path_storage.h"           // for rendering outlined polygons

#include "WrappersForAGG.h"

// taken from pixel_formats.h in the examples folder of the AGG distribution package:
#if defined RS_SYSTEM_LINUX
  #include "AGG_2.4/include/agg_pixfmt_rgba.h"
  #define pix_format agg::pix_format_argb32
  typedef agg::pixfmt_argb32 pixfmt;
  typedef agg::pixfmt_argb32_pre pixfmt_pre;
  typedef agg::rgba8 color_type;
  typedef agg::order_argb component_order;
#else
  #include "AGG_2.4/include/agg_pixfmt_rgba.h"
  #define pix_format agg::pix_format_bgra32
  typedef agg::pixfmt_bgra32 pixfmt;
  typedef agg::pixfmt_bgra32_pre pixfmt_pre;
  typedef agg::rgba8 color_type;
  typedef agg::order_bgra component_order;
#endif

namespace rsAGG
{
  typedef agg::renderer_base<pixfmt> renderer_base;
  typedef agg::renderer_scanline_aa_solid<renderer_base> renderer_aa;

  inline agg::rgba8 toAggColor(const rsColorRGBA &c)
  {
    return agg::rgba8(c.r, c.g, c.b, c.a);
  }

  void renderRasterizedScanlines(agg::rasterizer_scanline_aa<> &rasterizer, 
    rsImageRegionRGBA &imageRegion, const rsColorRGBA &color)
  {
    agg::rendering_buffer rbuf((rsUint8*) imageRegion.getPointerToPixel(0, 0), 
      imageRegion.getWidth(), imageRegion.getHeight(), imageRegion.getLineStrideInBytes());
    pixfmt pixf(rbuf);
    renderer_base rb(pixf);
    renderer_aa ren_aa(rb);
    ren_aa.color(toAggColor(color));
    agg::scanline_u8 scanline;
    agg::render_scanlines(rasterizer, scanline, ren_aa);
  }

  void drawFilledPolygon(const rsPolygon2D<double> &polygon, rsImageRegionRGBA &imageRegion, 
    const rsColorRGBA &color)
  {
    agg::rasterizer_scanline_aa<> ras;
    ras.move_to_d(polygon.getVertex(0).x, polygon.getVertex(0).y);
    for(int i = 1; i < polygon.getNumVertices(); i++)
    {
      //double x = polygon.getVertex(i).x;
      //double y = polygon.getVertex(i).y;
      ras.line_to_d(polygon.getVertex(i).x, polygon.getVertex(i).y);
    }
    renderRasterizedScanlines(ras, imageRegion, color);
  }

  void drawOutlinedPolygon(const rsPolygon2D<double> &polygon, rsImageRegionRGBA &imageRegion, 
    const rsColorRGBA &color, double strokeWidth)
  {
    if( polygon.getNumVertices() < 1 )
      return;

    agg::rasterizer_scanline_aa<> ras;    
    agg::path_storage path;

    path.move_to(polygon.getVertex(0).x, polygon.getVertex(0).y);
    for(int i = 0; i < polygon.getNumVertices(); i++)
      path.line_to(polygon.getVertex(i).x, polygon.getVertex(i).y);
    path.close_polygon();

    agg::line_cap_e           cap = agg::butt_cap;
    //if(m_cap.cur_item() == 1) cap = agg::square_cap;
    //if(m_cap.cur_item() == 2) cap = agg::round_cap;
        
    agg::line_join_e           join = agg::miter_join;
    //if(m_join.cur_item() == 1) join = agg::miter_join_revert;
    //if(m_join.cur_item() == 2) join = agg::round_join;
    //if(m_join.cur_item() == 3) join = agg::bevel_join;

    agg::conv_stroke<agg::path_storage> stroke(path);
    stroke.line_join(join);
    stroke.line_cap(cap);
    stroke.miter_limit(8.0);        // later: use a function-parameter
    stroke.width(strokeWidth);
    ras.add_path(stroke);

    renderRasterizedScanlines(ras, imageRegion, color);
  }
 
  /*
  rsPolygon2D<double> generatePolygonOutline(const rsPolygon2D<double> pIn)
  {
    rsPolygon2D<double> pOut;  // output polygon
    pOut = pIn;                // preliminary
    return pOut;  
  }
  */

}
