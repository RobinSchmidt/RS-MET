#ifndef RS_POLYGON2DTESTS_H
#define RS_POLYGON2DTESTS_H

#include "../../../Shared/Shared.h"

bool testPolygon2D();

bool testRegularPolygonCreation2D(std::string &reportString);
//bool testPointInsidePolygon2D(std::string &reportString);

bool convexPolygonClipping(std::string &reportString);

bool pixelCoverage(std::string &reportString);



#endif
