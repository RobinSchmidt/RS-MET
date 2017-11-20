using namespace RSLib;

rsColorRGBA::rsColorRGBA(rsUint8 red, rsUint8 green, rsUint8 blue, rsUint8 alpha) 
: r(red)
, g(green)
, b(blue)
, a(alpha)
{

}

rsColorRGBA::rsColorRGBA(int index)
{
  switch( index )
  {
  case black:       *this = rsColorRGBA(  0,   0,   0, 255); break;
  case blue:        *this = rsColorRGBA(  0,   0, 255, 255); break;
  case cyan:        *this = rsColorRGBA(  0, 255, 255, 255); break;
  case green:       *this = rsColorRGBA(  0, 255,   0, 255); break;
  case magenta:     *this = rsColorRGBA(255,   0, 255, 255); break;
  case red:         *this = rsColorRGBA(255,   0,   0, 255); break;
  case transparent: *this = rsColorRGBA(  0,   0,   0,   0); break;
  case white:       *this = rsColorRGBA(255, 255, 255, 255); break;
  case yellow:      *this = rsColorRGBA(255, 255,   0, 255); break;
  default:          *this = rsColorRGBA(  0,   0,   0,   0); break;
  }
}

rsColorRGBA::rsColorRGBA(float grayValue)
{
  r = g = b = (unsigned char) (grayValue * 255.f);
  a = 255;
}


rsColorFloatRGBA::rsColorFloatRGBA(float red, float green, float blue, float alpha)
: r(red)
, g(green)
, b(blue)
, a(alpha)
{

}



