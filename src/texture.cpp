#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

  Color Texture::sample(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.


      if (sp.lsm == L_ZERO) {
          if (sp.psm == P_NEAREST) {
              return sample_nearest(sp.p_uv, 0);
          }
          return sample_bilinear(sp.p_uv, 0);
      }

      else if (sp.lsm == L_NEAREST) {

          float level = get_level(sp);
          
          level = round(level);
          
          //cout << level << "\n";
          if (level < 0) {
              level = 0;
          }
          if (level > mipmap.size()) {
              level = mipmap.size();
          }

          if (sp.psm == P_NEAREST) {
              return sample_nearest(sp.p_uv, level);
          }
          return sample_bilinear(sp.p_uv, level);

      }

      else if (sp.lsm == L_LINEAR) {

          float level = get_level(sp);
          float levela = floor(level);
          float levelb = ceil(level);
          float d = level - levela;

          if (levela < 0) {
              levela = 0;
          }
          if (levela >= mipmap.size()) {
              levela = mipmap.size() - 1;
          }
          if (levelb < 0) {
              levelb = 0;
          }
          if (levelb >= mipmap.size()) {
              levelb = mipmap.size() - 1;
          }

          if (sp.psm == P_NEAREST) {
              return sample_nearest(sp.p_uv, levela) * (1-d) + sample_nearest(sp.p_uv, levelb) * d;
          }
          return sample_bilinear(sp.p_uv, levela) * (1 - d) + sample_bilinear(sp.p_uv, levelb) * d;


      }




// return magenta for invalid level
    return Color(1, 0, 1);
  }

  float Texture::get_level(const SampleParams& sp) {
    // TODO: Task 6: Fill this in.

      Vector2D x = Vector2D(sp.p_dx_uv.x * (width - 1), sp.p_dx_uv.y * (height - 1));

      Vector2D y = Vector2D(sp.p_dy_uv.x * (width - 1), sp.p_dy_uv.y * (height - 1));


      float xnorm = x.norm();
      float ynorm = y.norm();

      float m = max(xnorm, ynorm);

      return log2(m);
  }

  Color MipLevel::get_texel(int tx, int ty) {
    return Color(&texels[tx * 3 + ty * width * 3]);
  }

  Color Texture::sample_nearest(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
    auto& mip = mipmap[level];
    
    Vector2D new_uv = Vector2D(uv.x * (mip.width - 1), uv.y * (mip.height - 1));
    Vector2D rounded = Vector2D(round(new_uv.x), round(new_uv.y));

    return mip.get_texel(rounded.x, rounded.y);
    //cout << "scaled uv: " << new_uv << "  ";

    uv = new_uv;
    Vector2D d = new_uv - Vector2D(0.5, 0.5);

    float ufloor = floor(d.x);
    float vfloor = floor(d.y);
    float uceil = ceil(d.x);
    float vceil = ceil(d.y);

    //cout << ufloor << " " << uceil << " " << vfloor << " " << vceil << "\n";
    
    

    Vector2D u00 = Vector2D(ufloor, vfloor);
    Vector2D u10 = Vector2D(uceil, vfloor);
    Vector2D u01 = Vector2D(ufloor, vceil);
    Vector2D u11 = Vector2D(uceil, vceil);

    float du00 = (d - u00).norm();
    float du01 = (d - u01).norm();
    float du10 = (d - u10).norm();
    float du11 = (d - u11).norm();

    float minNorm = min(du00, min(du01, min(du10, du11)));




    if (du00 == minNorm) { 

        //cout << mip.get_texel(u00.x, u00.y);
        return mip.get_texel(u00.x, u00.y );
    }
    else if (du01 == minNorm) {
        return mip.get_texel(u01.x, u01.y);

    }
    else if (du10 == minNorm) {
        return mip.get_texel(u10.x, u10.y);

    }
    else if (du11 == minNorm) {
        return mip.get_texel(u11.x, u11.y);

    }

    // return magenta for invalid level
    return Color(1, 0, 1);
  }

  Color Texture::sample_bilinear(Vector2D uv, int level) {
    // TODO: Task 5: Fill this in.
      //cout << level;
    auto& mip = mipmap[level];
    Vector2D new_uv = Vector2D(uv.x * (mip.width - 1), uv.y * (mip.height - 1));

    uv = new_uv;
    Vector2D d = new_uv; // -Vector2D(0.5, 0.5);
    
    float ufloor = floor(d.x);
    float vfloor = floor(d.y);
    float uceil = ceil(d.x);
    float vceil = ceil(d.y);

    float s = d.x - ufloor;
    float t = d.y - vfloor;


    //cout << ufloor << " " << uceil << " " << vfloor << " " << vceil << "\n";
    Color u00 = mip.get_texel(ufloor, vfloor);
    Color u10 = mip.get_texel(uceil, vfloor);
    
    Color u01 = mip.get_texel(ufloor, vceil);
    
    Color u11 = mip.get_texel(uceil, vceil);
    
    
    //cout << u00 << u10 << u01 << u11 << "\n";

    Color u0 = u00 * s + (1 - s) * u10;
    Color u1 = u01 * s + (1 - s) * u11;

    //cout << u0 << u1 << "\n";

    Color f = u0 * t + (1 - t) * u1;
    return f;
    


    // return magenta for invalid level
    return Color(1, 0, 1);
  }



  /****************************************************************************/

  // Helpers

  inline void uint8_to_float(float dst[3], unsigned char* src) {
    uint8_t* src_uint8 = (uint8_t*)src;
    dst[0] = src_uint8[0] / 255.f;
    dst[1] = src_uint8[1] / 255.f;
    dst[2] = src_uint8[2] / 255.f;
  }

  inline void float_to_uint8(unsigned char* dst, float src[3]) {
    uint8_t* dst_uint8 = (uint8_t*)dst;
    dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
    dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
    dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  }

  void Texture::generate_mips(int startLevel) {

    // make sure there's a valid texture
    if (startLevel >= mipmap.size()) {
      std::cerr << "Invalid start level";
    }

    // allocate sublevels
    int baseWidth = mipmap[startLevel].width;
    int baseHeight = mipmap[startLevel].height;
    int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

    numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
    mipmap.resize(startLevel + numSubLevels + 1);

    int width = baseWidth;
    int height = baseHeight;
    for (int i = 1; i <= numSubLevels; i++) {

      MipLevel& level = mipmap[startLevel + i];

      // handle odd size texture by rounding down
      width = max(1, width / 2);
      //assert (width > 0);
      height = max(1, height / 2);
      //assert (height > 0);

      level.width = width;
      level.height = height;
      level.texels = vector<unsigned char>(3 * width * height);
    }

    // create mips
    int subLevels = numSubLevels - (startLevel + 1);
    for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
      mipLevel++) {

      MipLevel& prevLevel = mipmap[mipLevel - 1];
      MipLevel& currLevel = mipmap[mipLevel];

      int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
      int currLevelPitch = currLevel.width * 3; // 32 bit RGB

      unsigned char* prevLevelMem;
      unsigned char* currLevelMem;

      currLevelMem = (unsigned char*)&currLevel.texels[0];
      prevLevelMem = (unsigned char*)&prevLevel.texels[0];

      float wDecimal, wNorm, wWeight[3];
      int wSupport;
      float hDecimal, hNorm, hWeight[3];
      int hSupport;

      float result[3];
      float input[3];

      // conditional differentiates no rounding case from round down case
      if (prevLevel.width & 1) {
        wSupport = 3;
        wDecimal = 1.0f / (float)currLevel.width;
      }
      else {
        wSupport = 2;
        wDecimal = 0.0f;
      }

      // conditional differentiates no rounding case from round down case
      if (prevLevel.height & 1) {
        hSupport = 3;
        hDecimal = 1.0f / (float)currLevel.height;
      }
      else {
        hSupport = 2;
        hDecimal = 0.0f;
      }

      wNorm = 1.0f / (2.0f + wDecimal);
      hNorm = 1.0f / (2.0f + hDecimal);

      // case 1: reduction only in horizontal size (vertical size is 1)
      if (currLevel.height == prevLevel.height) {
        //assert (currLevel.height == 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          for (int ii = 0; ii < wSupport; ii++) {
            uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
            result[0] += wWeight[ii] * input[0];
            result[1] += wWeight[ii] * input[1];
            result[2] += wWeight[ii] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (3 * i), result);
        }

        // case 2: reduction only in vertical size (horizontal size is 1)
      }
      else if (currLevel.width == prevLevel.width) {
        //assert (currLevel.width == 1);

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          result[0] = result[1] = result[2] = 0.0f;
          for (int jj = 0; jj < hSupport; jj++) {
            uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
            result[0] += hWeight[jj] * input[0];
            result[1] += hWeight[jj] * input[1];
            result[2] += hWeight[jj] * input[2];
          }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + (currLevelPitch * j), result);
        }

        // case 3: reduction in both horizontal and vertical size
      }
      else {

        for (int j = 0; j < currLevel.height; j++) {
          hWeight[0] = hNorm * (1.0f - hDecimal * j);
          hWeight[1] = hNorm;
          hWeight[2] = hNorm * hDecimal * (j + 1);

          for (int i = 0; i < currLevel.width; i++) {
            wWeight[0] = wNorm * (1.0f - wDecimal * i);
            wWeight[1] = wNorm * 1.0f;
            wWeight[2] = wNorm * wDecimal * (i + 1);

            result[0] = result[1] = result[2] = 0.0f;

            // convolve source image with a trapezoidal filter.
            // in the case of no rounding this is just a box filter of width 2.
            // in the general case, the support region is 3x3.
            for (int jj = 0; jj < hSupport; jj++)
              for (int ii = 0; ii < wSupport; ii++) {
                float weight = hWeight[jj] * wWeight[ii];
                uint8_to_float(input, prevLevelMem +
                  prevLevelPitch * (2 * j + jj) +
                  3 * (2 * i + ii));
                result[0] += weight * input[0];
                result[1] += weight * input[1];
                result[2] += weight * input[2];
              }

            // convert back to format of the texture
            float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
          }
        }
      }
    }
  }

}
