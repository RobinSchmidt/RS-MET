bool triangleRasterization()
{
  bool r = true;

  typedef rsVector2DF Vec2;    // for convenience
  typedef Vec2 V;              // even shorter
  float c = 1.0f;              // color (gray value)
  //rsImageF img(13,10);         // image to draw on
  rsImageF img(20,20);         // image to draw on
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


  // a somewhat left-pointing triangle with half-integer vertices:
  A = V(1.5f,1.5f), B = V(3.5f,8.5f), C = V(7.5f,4.5f);
  drawTriangleAntiAliasedSpanBased(drw, A, B, C, c);

  // a somewhat right-pointing triangle with half-integer vertices:
  A = V(7.5f,1.5f), B = V(1.5f,5.5f), C = V(9.5f,7.5f);
  //drawTriangleAntiAliasedSpanBased(drw, A, B, C, c);

  // with integer vertices:
  A = V(7,1), B = V(1,5), C = V(9,7);
  //drawTriangleAntiAliasedSpanBased(drw, A, B, C, c);


  // an example that turned out to be problematic in the random test:
  //A = V(3.7f,15.6f), B = V(10.3f,13.2f), C = V(21.8f,19.5f);
  A = V(5.3f,15.7f), B = V(10.9f,11.7f), C = V(6.0f,-0.1f);
  img.clear();
  drawTriangleAntiAliasedSpanBased(drw, A, B, C, c);
  writeScaledImageToFilePPM(img, "TestTriangleSpan.ppm", 16);
  img.clear();
  drawTriangleAntiAliasedProto(drw, A, B, C, c);
  writeScaledImageToFilePPM(img, "TestTriangleNaive.ppm", 16);
  // seems like the last (19) pixel in the second-to-last (18) line is drawn twice


  // we should test at least with integer and half-integer vertex coordinates - see, if the
  // loop min/max values always get the correct values

  //writeScaledImageToFilePPM(img, "TriangleTest.ppm", 16);
  return r;
}

bool triangleRasterization2()
{
  // We have 3 versions of the triangle drawing function - the naive, slow version, and the 
  // box-based and span-based optimizations - we check here, if the optimized versions produce the 
  // same results as the naive version.

  bool r = true;
  int numTriangles = 10; // it takes quite a lot of time to draw triangles, so we use a small number

                        // create and set up images and drawer objects:
  int w = 20;
  int h = 20;
  rsImageF img1(w, h), img2(w, h), img3(w, h);
  rsImageDrawerFFF drw1(&img1), drw2(&img2), drw3(&img3);
  drw1.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);
  drw2.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);
  drw3.setBlendMode(rsImageDrawerFFF::BLEND_ADD_CLIP);

  // random number generators for the triangle vertices and color:
  rsNoiseGenerator<float> coordinateGenerator;
  coordinateGenerator.setRange(-10, 30);
  rsNoiseGenerator<float> colorGenerator;
  colorGenerator.setRange(0, 1);

  // draw triangles using the different functions and compare results:
  float color;
  rsVector2D<float> a, b, c;
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

    // draw it with the different versions of the drawing function:
    drawTriangleAntiAliasedProto(   drw1, a, b, c, color);
    drawTriangleAntiAliasedBoxBased(drw2, a, b, c, color);
    drawTriangleAntiAliasedSpanBased(drw3, a, b, c, color);



    // compare drawing results:
    float tol = float(1.e-5);
    r &= img2.areAllPixelsEqualTo(&img1, tol);
    r &= img3.areAllPixelsEqualTo(&img1, tol);
      // remaining error seems to be of numerical nature (the naive and box-based versiona produce a 
      // nonzero value of the order of 10^-5 where the optimized version returns 0 bcs the pixel is 
      // considered already outside the nonzero range) - but 10^-5 seems a bit much, even for 
      // single precision (i'd expect something like 10^-7)...hmmm...

    // for debug:
    //writeScaledImageToFilePPM(img1, "TriangleNaive.ppm", 16);
    //writeScaledImageToFilePPM(img2, "TriangleBoxOptimization.ppm", 16);
    //writeScaledImageToFilePPM(img3, "TriangleSpanOptimization.ppm", 16);
    std::vector<float> v1, v2, v3, err;
    v1 = img1.toStdVector();
    v2 = img2.toStdVector();
    v3 = img3.toStdVector();
    err = v1-v3;
    int errMaxIdx = RAPT::rsArrayTools::maxDeviationIndex(&v1[0], &v3[0], w*h);
    float errMax = err[errMaxIdx];

    int dummy = 0;
  }

  //writeImageToFilePPM(img1, "RandomTriangles.ppm");
  return r;
}