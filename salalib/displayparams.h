#pragma once

struct DisplayParams
{
   enum { AXMANESQUE = 0, GREYSCALE = 1, MONOCHROME = 2, DEPTHMAPCLASSIC = 3, PURPLEORANGE = 4, BLUERED = 5, HUEONLYAXMANESQUE = 6,
          VIRIDIS = 7, MAGMA = 8 }; // 7, 8: perceptually uniform, colour-blind friendly scales
   float blue;
   float red;
   int colorscale;
   DisplayParams()
   { blue = 0.0f; red = 1.0f; colorscale = 0; }
};
