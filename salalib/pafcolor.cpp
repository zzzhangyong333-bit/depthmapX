// sala - a component of the depthmapX - spatial network analysis platform
// Copyright (C) 2011-2012, Tasos Varoudis

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <math.h>
#include <float.h> // _finite support

#include <salalib/pafcolor.h>  // <- sala actually includes vertex.h for us

static unsigned int g_nicecolor[] = { 
   0x003333DD, // 0 blue
   0x003388DD, // 1
   0x0022CCDD, // 2
   0x0022CCBB, // 3
   0x0022DD88, // 4
   0x0088DD22, // 5
   0x00BBCC22, // 6
   0x00DDCC22, // 7
   0x00DD8833, // 8
   0x00DD3333, // 9 red
};

// Test a range designed to try to keep consistent saturation and brightness of g_nicecolor, and only move hue
static unsigned int g_nicecolorhsb[] = {
   0x003333DD, // 0 blue
   0x003377DD, // 1
   0x0033BBDD, // 2
   0x0033DDBB, // 3
   0x0033DD55, // 4
   0x0055DD33, // 5
   0x00BBDD33, // 6
   0x00DDBB33, // 7
   0x00DD7733, // 8
   0x00DD3333, // 9 red
};

static unsigned int g_greyscale[] = { 
   0x00000000, // 0 black
   0x00444444, // 1
   0x00777777, // 2
   0x00AAAAAA, // 3
   0x00CCCCCC, // 4
   0x00EEEEEE, // 5
   0x00FFFFFF, // 6 white
};

static unsigned int g_bluered[] = {
   0x004575B4, // 0 blue
   0x0091BFDB, 
   0x00E0F3F8, 
   0x00FFFFBF, 
   0x00FEE090, 
   0x00FC8D59, 
   0x00D73027  // 6 red 
};

static unsigned int g_purpleorange[] = {
   0x00542788, // 0 purple
   0x00998EC3, // 1
   0x00D8DAEB, // 2
   0x00F7F7F7, // 3
   0x00FEE0B6, // 4
   0x00F1A340, // 5
   0x00B35806  // 6 orange
};

// Perceptually uniform scales (from matplotlib, CC0): 11 evenly spaced stops, interpolated linearly.
// Unlike the banded scales above these are continuous, read correctly in greyscale print and are
// safe for the common forms of colour blindness.
static unsigned int g_viridis[] = {
   0x00440154, 0x00482475, 0x00414487, 0x00355F8D, 0x002A788E, 0x0021918C,
   0x0022A884, 0x0044BF70, 0x007AD151, 0x00BDDF26, 0x00FDE725
};

static unsigned int g_magma[] = {
   0x00000004, 0x00140E36, 0x003B0F70, 0x00641A80, 0x008C2981, 0x00B73779,
   0x00DE4968, 0x00F7705C, 0x00FE9F6D, 0x00FECF92, 0x00FCFDBF
};

static unsigned int interpolateStops(const unsigned int *stops, int count, double field)
{
   double pos = field * (count - 1);
   int lower = int(pos);
   if (lower >= count - 1) {
      return stops[count - 1];
   }
   if (lower < 0) {
      return stops[0];
   }
   double t = pos - lower;
   unsigned int a = stops[lower], b = stops[lower + 1];
   unsigned int result = 0;
   for (int shift = 0; shift <= 16; shift += 8) {
      double ca = (a >> shift) & 0xFF, cb = (b >> shift) & 0xFF;
      unsigned int c = (unsigned int)(ca + (cb - ca) * t + 0.5);
      result |= (c & 0xFF) << shift;
   }
   return result;
}

// htmlByte converts a normalised number to an HTML safe byte

unsigned char htmlByte(double colorByte)
{
   // Quick mod - TV
#if defined(_MSC_VER)   
   return (unsigned char((colorByte + 0.0333) * 15.0) * 0x11);
#else
   return ((unsigned char)((colorByte + 0.0333) * 15.0) * 0x11);
#endif   
}

PafColor& PafColor::makeColor(double field, DisplayParams dp)
{
   // Quick mod - TV
   if (field == -1.0 || std::isnan(field)) {
      // -1.0 is (currently) a nan value, set alpha channel to 0 (transparent)
      switch (dp.colorscale) {
         case DisplayParams::MONOCHROME:
         case DisplayParams::GREYSCALE:
            m_color = 0x00000000;   // <- monochrome and greyscale, simply hide
            break;
         default:
            // if in colour, then show greyed out:
            m_color = 0x007f7f7f;   // <- grey retained for visibility on certain values
            break;
      }
      return *this;
   }
   if (dp.blue > dp.red) {
      field = 1.0 - field;
      dp.blue = 1.0f - dp.blue;
      dp.red = 1.0f - dp.red;
   }
   if (dp.colorscale == DisplayParams::DEPTHMAPCLASSIC) {
      makeDepthmapClassic(field, dp.blue, dp.red);
   }
   else {
      field = (field - dp.blue) / (dp.red - dp.blue);
      // Quick mod - TV
      if (std::isnan(field)) {
         field = 0.5;
      }
      if (field > 1.0) {
         field = 1.0;
      }
      else if (field < 0.0) {
         field = 0.0;
      }
      switch(dp.colorscale) {
      case DisplayParams::AXMANESQUE:
         makeAxmanesque(field);
         break;
      case DisplayParams::HUEONLYAXMANESQUE:
         makeHueOnlyAxmanesque(field);
         break;
      case DisplayParams::PURPLEORANGE:
         makePurpleOrange(field);
         break;
      case DisplayParams::BLUERED:
         makeBlueRed(field);
         break;
      case DisplayParams::VIRIDIS:
         makeViridis(field);
         break;
      case DisplayParams::MAGMA:
         makeMagma(field);
         break;
      case DisplayParams::GREYSCALE:
      case DisplayParams::MONOCHROME:
         makeGreyScale(field);
         break;
      }
   }
   return *this;
}

// this makes an Axman-like colour range

PafColor& PafColor::makeAxmanesque( double field )
{
   m_color = 0xff000000 | g_nicecolor[int((field - 1e-9) * 10.0)];
   return *this;
}

PafColor& PafColor::makeHueOnlyAxmanesque( double field )
{
   m_color = 0xff000000 | g_nicecolorhsb[int((field - 1e-9) * 10.0)];
   return *this;
}

// this makes a purple-orange scheme that is red-green colour-blind safe

PafColor& PafColor::makePurpleOrange( double field )
{
   m_color = 0xff000000 | g_purpleorange[int((field - 1e-9) * 7.0)];
   return *this;
}

// this makes a blue-red scheme that is red-green colour-blind safe

PafColor& PafColor::makeBlueRed( double field )
{
   m_color = 0xff000000 | g_bluered[int((field - 1e-9) * 7.0)];
   return *this;
}

PafColor& PafColor::makeViridis( double field )
{
   m_color = 0xff000000 | interpolateStops(g_viridis, 11, field);
   return *this;
}

PafColor& PafColor::makeMagma( double field )
{
   m_color = 0xff000000 | interpolateStops(g_magma, 11, field);
   return *this;
}

// this makes a greyscale colour range

PafColor& PafColor::makeGreyScale( double field )
{
   m_color = 0xff000000 | g_greyscale[int((field - 1e-9) * 7.0)];
   return *this;
}

// note, makeDepthmapClassic converts to a safe HTML colour

PafColor& PafColor::makeDepthmapClassic( double field, double blue, double red )
{
   m_color = 0xff000000; // set alpha to 255, solid colour
   double green = blue + (red-blue) / 10.0;
   // NB previously included colour muting: the 1.0 was originally 0.9 to mute the colours slightly
   if (field >= 0.0 && field < blue) {
      setr(htmlByte(0.5 * (blue - field)/blue * 1.0));
      // Quick mod - TV
#if defined(_MSC_VER)
      setb(unsigned char(0xFF));
#else
      setb((unsigned char)(0xFF));
#endif
   }
   else if (field >= blue && field < (green+blue)/2.0) {
      // Quick mod - TV
#if defined(_MSC_VER)      
      setb(unsigned char(0xFF));
#else
      setb((unsigned char)(0xFF));
#endif      
      setg(htmlByte((2.0*(field - blue)/(green-blue)) * 1.0));
   }
   else if (field >= (green+blue)/2.0 && field < green) {
      setb(htmlByte((2.0*(green - field)/(green-blue)) * 1.0));
      // Quick mod - TV
#if defined(_MSC_VER)
      setg(unsigned char(0xFF));
#else
      setg((unsigned char)(0xFF));
#endif      
   }
   else if (field >= green && field < (green+red)/2.0 ) {
      // Quick mod - TV
#if defined(_MSC_VER)      
      setg(unsigned char(0xFF));
#else
      setg((unsigned char)(0xFF));
#endif      
      setr(htmlByte((2.0*(field - green)/(red-green)) * 1.0));
   }
   else if (field >= (green+red)/2.0 && field < red) {
      setg(htmlByte((2.0*(red - field)/(red-green)) * 1.0));
      // Quick mod - TV
#if defined(_MSC_VER)      
      setr(unsigned char(0xFF));
#else
      setr((unsigned char)(0xFF));
#endif      
   }
   else if (field >= red) {
      // Quick mod - TV
#if defined(_MSC_VER)      
      setr(unsigned char(0xFF));
#else
      setr((unsigned char)(0xFF));
#endif      
      setb(htmlByte(0.5 * (field - red)/(1.0 - red) * 1.0));
   }
   return *this;
}
