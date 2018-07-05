

//-------------------------------------------------------------------------------------------------
// Triangles

// code based on:  https://www.youtube.com/watch?v=9A5TVh6kPLA
// v0 should be left to v1 and v2 below the line connecting v0 and v1
void drawTriangleFlatTop(rsImageDrawerFFF& drw, 
  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color)
{
  static const float d = 0.5f;                // offset of a pixel's center from its index/coord
  float m0 = (v2.x - v0.x) / (v2.y - v0.y);   // inverse of slope of line from v0 to v2
  float m1 = (v2.x - v1.x) / (v2.y - v1.y);   // inverse of slope of line from v1 to v2
  for(int y = (int)ceil(v0.y-d); y < (int)ceil(v2.y-d); y++) { // loop over scanlines
    float px0 = m0 * (float(y) - v0.y + d) + v0.x;             // start x-coord
    float px1 = m1 * (float(y) - v1.y + d) + v1.x;             // end x-coord
    for(int x = (int)ceil(px0-d); x < (int)ceil(px1-d); x++)   // loop over pixels in scanline
      drw.plot(x, y, color);
  }
  // the ceil-function together with the offset of 0.5 amounts to the top-left rule
}
// maybe compute px0, px1 incrementally, i.e. init to v0.x, v1.x and incerement by dx0, dx1 in each
// iteration, where dx0 =
// but not in the prototype

// v1 should be left to v2 and v0 above the line connecting v1 and v2
void drawTriangleFlatBottom(rsImageDrawerFFF& drw, 
  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color)
{
  float m0 = (v1.x - v0.x) / (v1.y - v0.y);   // inverse of slope of line from v0 to v1
  float m1 = (v2.x - v0.x) / (v2.y - v0.y);   // inverse of slope of line from v0 to v2
  int ys = (int) ceil(v0.y - 0.5f);           // y-coord of first scanline
  int ye = (int) ceil(v2.y - 0.5f);           // y-coord of scanline after the last line drawn
  for( int y = ys; y < ye; y++ ) {            // loop over scanlines
    float px0 = m0 * (float(y) + 0.5f - v0.y) + v0.x; // start x-coord
    float px1 = m1 * (float(y) + 0.5f - v0.y) + v0.x; // end x-coord
    int xs = (int) ceil(px0 - 0.5f);                  // start pixel
    int xe = (int) ceil(px1 - 0.5f);                  // end pixel (after the last pixel drawn)
    for(int x = xs; x < xe; x++)                      // loop over pixels in current scanline
      drw.plot(x, y, color);
  }
}
// compactify this further (get rid of ys,ye,xs,xe, use d=0.5)

// i think, the anti-aliased version should let loop indices start at = floor(...) and end at
// <= ceil(...)

// maybe wrap into class rsPolygonDrawer

void drawTriangle(rsImageDrawerFFF& drw, 
  const rsVector2DF& v0, const rsVector2DF& v1, const rsVector2DF& v2, float color)
{
  // use pointers so we can swap (for sorting purposes)
  typedef rsVector2DF Vec2; // for convenience
  const Vec2* pv0 = &v0;
  const Vec2* pv1 = &v1;
  const Vec2* pv2 = &v2;
  // todo: use pointers as arguments - maybe provide convenience function that takes const 
  // references

  // sort vertices by y:
  if(pv1->y < pv0->y) std::swap(pv0, pv1);
  if(pv2->y < pv1->y) std::swap(pv1, pv2);
  if(pv1->y < pv0->y) std::swap(pv0, pv1);

  if(pv0->y == pv1->y) {        // triangle is flat top
    if(pv1->x < pv0->x) std::swap(pv0, pv1); // sort top vertices by x
    drawTriangleFlatTop(drw, *pv0, *pv1, *pv2, color);
  }
  else if(pv1->y == pv2->y) {   // triangle is flat bottom
    if(pv2->x < pv1->x) std::swap( pv1,pv2 ); // sort bottom vertices by x
    drawTriangleFlatBottom(drw, *pv0, *pv1, *pv2, color);
  }
  else {
    // split general triangle into flat-top and flat-bottom:
    const float alpha = (pv1->y - pv0->y) / (pv2->y - pv0->y);
    const Vec2 vi = *pv0 + alpha * (*pv2 - *pv0);    // splitting vertex by linear interpolation between v0 and v2
    if(pv1->x < vi.x) { // long side is on the right (major right)
      drawTriangleFlatBottom(drw, *pv0, *pv1,   vi, color);
      drawTriangleFlatTop(   drw, *pv1,   vi, *pv2, color);
    }
    else {              // long side is on the left (major left)
      drawTriangleFlatBottom(drw, *pv0,   vi, *pv1, color);
      drawTriangleFlatTop(   drw,   vi, *pv1, *pv2, color);
    }
  }
}