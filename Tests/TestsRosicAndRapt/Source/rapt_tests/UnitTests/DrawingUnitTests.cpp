bool triangleRasterization()
{
  bool r = true;

  typedef rsVector2DF Vec2;    // for convenience
  typedef Vec2 V;              // even shorter
  float c = 1.0f;              // color (gray value)
  rsImageF img(13,10);         // image to draw on
  rsImageDrawerFFF drw(&img);  // drawer object
  drw.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);

  Vec2 A = V(0.5,0.5), B = V(4.5,2.5), C = V(2.5,0.5);
  float x = 0, y = 0;

  std::vector<Vec2> 
    triangle = { A, B, C },
    square   = { V(x, y), V(x, y+1), V(x+1, y+1), V(x+1, y) },
    polygon  = clipPolygon(triangle, square);

  //drawTriangleAntiAliasedProto(drw, A, B, C, c);
  //drawTriangleAntiAliased(drw, V(0.5,0.5), V(2.5,0.5), V(4.5,2.5), c);


  A = V(7.5f,1.5f), B = V(1.5f,5.5f), C = V(9.5f,7.5f);
  drawTriangleAntiAliasedSpanBased(drw, A, B, C, c);

  A = V(7,1), B = V(1,5), C = V(9,7);
  //drawTriangleAntiAliasedSpanBased(drw, A, B, C, c);


  // we should test at least with integer and half-integer vertex coordinates - see, if the
  // loop min/max values always get the correct values

  writeScaledImageToFilePPM(img, "TriangleTest.ppm", 16);
  return r;
}

bool triangleRasterization2()
{
  // We have 3 versions of the triangle drawing function - the naive, slow version, and the 
  // box-based and span-based optimizations - we check here, if the optimized versions produce the 
  // same results as the naive version.

  bool r = true;
  int numTriangles = 1; // it takes quite a lot of time to draw triangles, so we use a small number

                        // create and set up images and drawer objects:
  int w = 20;
  int h = 20;
  rsImageF img1(w, h), img2(w, h), img3(w, h);
  rsImageDrawerFFF drw1(&img1), drw2(&img2), drw3(&img3);
  drw1.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);
  drw2.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);
  drw3.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);

  // random number generators for the triangle vertices and color:
  rsNoiseGeneratorF coordinateGenerator;
  coordinateGenerator.setRange(-10, 30);
  rsNoiseGeneratorF colorGenerator;
  colorGenerator.setRange(0, 1);

  // draw triangles using the different functions and compare results:
  float color;
  rsVector2DF a, b, c;
  for(int i = 1; i <= numTriangles; i++)
  {
    // we clear because otherwise, it's likely that all pixel soon saturate to white:
    img1.clear();
    img2.clear();
    img3.clear();

    // produce a random triangle with random color:
    a.x = coordinateGenerator.getSample();
    a.y = coordinateGenerator.getSample();
    b.x = coordinateGenerator.getSample();
    b.y = coordinateGenerator.getSample();
    c.x = coordinateGenerator.getSample();
    c.y = coordinateGenerator.getSample();
    color = colorGenerator.getSample();

    // draw it with the different versions of the drawing
    // maybe draw only if the triangle is counterclockwise? hmm...
    drawTriangleAntiAliasedProto(   drw1, a, b, c, color);
    drawTriangleAntiAliasedBoxBased(drw2, a, b, c, color);

    drawTriangleAntiAliasedSpanBased(drw3, a, b, c, color);
    // bottom scanline looks wrong


    writeScaledImageToFilePPM(img1, "TriangleNaive.ppm", 16);
    writeScaledImageToFilePPM(img2, "TriangleBoxOptimization.ppm", 16);
    writeScaledImageToFilePPM(img3, "TriangleSpanOptimization.ppm", 16);



    // compare drawing results:
    r &= img2.areAllPixelsEqualTo(&img1);
    r &= img3.areAllPixelsEqualTo(&img1);
    int dummy = 0;
  }

  //writeImageToFilePPM(img1, "RandomTriangles.ppm");
  return r;
}