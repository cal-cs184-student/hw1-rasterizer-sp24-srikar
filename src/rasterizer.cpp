#include "rasterizer.h"
#include <vector>


using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
      // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
      // NOTE: You are not required to implement proper supersampling for points and lines
      // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

      if (sample_buffer[y * width + x] == Color::White) {
          sample_buffer[y * width + x] = c;
      }
      else {
          sample_buffer[y * width + x] = -1 * sample_buffer[y * width + x] + (c);
      }
      sample_buffer[y * width + x] = c;

  }

  void RasterizerImp::fill_pixel(size_t x, size_t y, size_t k, Color c) {
      // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
      // NOTE: You are not required to implement proper supersampling for points and lines
      // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)


      //cout << x << " " << y << " " << c << "\n";

      
      

      sample_buffer[y * width * sample_rate + x * sample_rate + k] = c;

  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    for (int k = 0; k < sample_rate; k++) {
        fill_pixel(sx, sy, k, color);
    }

    
    return;
  }

  void RasterizerImp::rasterize_point(float x, float y, float k, Color color) {
      // fill in the nearest pixel
      int sx = (int)floor(x);
      int sy = (int)floor(y);
      int sk = (int)floor(k);

      // check bounds
      if (sx < 0 || sx >= width) return;
      if (sy < 0 || sy >= height) return;
      /*if (sk < 0 || sk >= sample_rate) return;*/

      fill_pixel(sx, sy, sk, color);
      return;
  }


  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {

        for (int k = 0; k < sample_rate; k++) {
            rasterize_point(pt[0], pt[1], k, color);
        }
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0, float x1, float y1, float x2, float y2, Color color) {
      
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

      int n = sqrt(sample_rate) - 1;
      vector<vector<float>> v{ {0.5}, {0.25, 0.75}, {1.0/6.0, 3.0/6.0, 5.0/6.0}, {1.0/8.0, 3.0/8.0, 5.0/8.0, 7.0/8.0}};

      int xmin = min(x0, min(x1, x2));
      int xmax = max(x0, max(x1, x2));
      int ymin = min(y0, min(y1, y2));
      int ymax = max(y0, max(y1, y2));


      for (int x = xmin; x <= xmax; x++) {
          for (int y = ymin; y <= ymax; y++) {
              float count = 0;

              for (int di = 0; di <= n; di++) {
                  float dx = v[n][di];
                  /*dx = ((float)rand()) / RAND_MAX;*/

                  for (int dj = 0; dj <= n; dj++) {
                      float dy = v[n][dj];
                      /*dy = ((float)rand()) / RAND_MAX;*/

                      float sx = x + dx;
                      float sy = y + dy;

                      float d01 = (sx - x0) * (y1 - y0) - (sy - y0) * (x1 - x0);
                      float d12 = (sx - x1) * (y2 - y1) - (sy - y1) * (x2 - x1);
                      float d02 = (sx - x2) * (y0 - y2) - (sy - y2) * (x0 - x2);

                      int k = dj * sqrt(sample_rate) + di;



                      if (d01 >= 0 && d12 >= 0 && d02 >= 0) {
                          rasterize_point(x, y, k, color);
                      }
                      else if (d01 <= 0 && d12 <= 0 && d02 <= 0) {
                          rasterize_point(x, y, k, color);
                      }

                  }
              }
          }
      }

    // TODO: Task 2: Update to implement super-sampled rasterization



  }
  

  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0, float x1, float y1, Color c1, float x2, float y2, Color c2){

      int n = sqrt(sample_rate) - 1;
      vector<vector<float>> v{ {0.5}, {0.25, 0.75}, {1.0 / 6.0, 3.0 / 6.0, 5.0 / 6.0}, {1.0 / 8.0, 3.0 / 8.0, 5.0 / 8.0, 7.0 / 8.0} };

      int xmin = min(x0, min(x1, x2));
      int xmax = max(x0, max(x1, x2));
      int ymin = min(y0, min(y1, y2));
      int ymax = max(y0, max(y1, y2));

      for (int x = xmin; x <= xmax; x++) {
          for (int y = ymin; y <= ymax; y++) {

              for (int di = 0; di <= n; di++) {
                  float dx = v[n][di];
                  /*dx = ((float)rand()) / RAND_MAX;*/

                  for (int dj = 0; dj <= n; dj++) {
                      float dy = v[n][dj];

                      float sx = x + dx;
                      float sy = y + dy;

                      float a = (-1 * (sx - x1) * (y2 - y1) + (sy - y1) * (x2 - x1));
                      a = a / (-1 * (x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));

                      float b = (-1 * (sx - x2) * (y0 - y2) + (sy - y2) * (x0 - x2));
                      b = b / (-1 * (x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));

                      float c = (-1 * (sx - x0) * (y1 - y0) + (sy - y0) * (x1 - x0));
                      c = c / (-1 * (x2 - x0) * (y1 - y0) + (y2 - y0) * (x1 - x0));

                      int k = dj * sqrt(sample_rate) + di;

                      if (a >= 0 && b >= 0 && c >= 0) {
                          rasterize_point(x, y, k, c0 * a + c1 * b + c2 * c);
                      }
                      else if (a <= 0 && b <= 0 && c <= 0) {
                          rasterize_point(x, y, k, c0 * a + c1 * b + c2 * c);
                      }

                  }

              }

              


          }

      }


  }

  

  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0, float x1, float y1, float u1, float v1, float x2, float y2, float u2, float v2, Texture& tex){
      int xmin = min(x0, min(x1, x2));
      int xmax = max(x0, max(x1, x2));
      int ymin = min(y0, min(y1, y2));
      int ymax = max(y0, max(y1, y2));

      //cout << x0 << " " << y0 << " " << u0 << " " << v0 << " " << x1 << " " << y1 << " " << u1 << " " << v1 << " " << x2 << " " << y2 << " " << u2 << " " << v2 << "\n";

      int n = sqrt(sample_rate) - 1;
      vector<vector<float>> v{ {0.5}, {0.25, 0.75}, {1.0 / 6.0, 3.0 / 6.0, 5.0 / 6.0}, {1.0 / 8.0, 3.0 / 8.0, 5.0 / 8.0, 7.0 / 8.0} };

      for (int x = xmin; x <= xmax; x++) {
          for (int y = ymin; y <= ymax; y++) {

              for (int di = 0; di <= n; di++) {
                  float dx = v[n][di];


                  for (int dj = 0; dj <= n; dj++) {
                      float dy = v[n][dj];
                      
                      float sx = x + dx;
                      float sy = y + dy;

                      int k = dj * sqrt(sample_rate) + di;

                      float a = (-1 * (sx - x1) * (y2 - y1) + (sy - y1) * (x2 - x1));
                      a = a / (-1 * (x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));

                      float b = (-1 * (sx - x2) * (y0 - y2) + (sy - y2) * (x0 - x2));
                      b = b / (-1 * (x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));

                      float c = (-1 * (sx - x0) * (y1 - y0) + (sy - y0) * (x1 - x0));
                      c = c / (-1 * (x2 - x0) * (y1 - y0) + (y2 - y0) * (x1 - x0));
                      
                      float d01 = (sx - x0) * (y1 - y0) - (sy - y0) * (x1 - x0);
                      float d12 = (sx - x1) * (y2 - y1) - (sy - y1) * (x2 - x1);
                      float d02 = (sx - x2) * (y0 - y2) - (sy - y2) * (x0 - x2);

                      float u = u0 * a + u1 * b + u2 * c;
                      float v = v0 * a + v1 * b + v2 * c;
                      Vector2D uv = Vector2D(u, v);

                      if (!((d01 >= 0 && d12 >= 0 && d02 >= 0) || (d01 <= 0 && d12 <= 0 && d02 <= 0))) {
                          continue;
                      }

                      float ax = (-1 * (sx + 1 - x1) * (y2 - y1) + (sy - y1) * (x2 - x1)) / (-1 * (x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                      float bx = (-1 * (sx + 1 - x2) * (y0 - y2) + (sy - y2) * (x0 - x2)) / (-1 * (x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                      float cx = (-1 * (sx + 1 - x0) * (y1 - y0) + (sy - y0) * (x1 - x0)) / (-1 * (x2 - x0) * (y1 - y0) + (y2 - y0) * (x1 - x0));

                      float ay = (-1 * (sx - x1) * (y2 - y1) + (sy + 1 - y1) * (x2 - x1)) / (-1 * (x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                      float by = (-1 * (sx - x2) * (y0 - y2) + (sy + 1 - y2) * (x0 - x2)) / (-1 * (x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                      float cy = (-1 * (sx - x0) * (y1 - y0) + (sy + 1 - y0) * (x1 - x0)) / (-1 * (x2 - x0) * (y1 - y0) + (y2 - y0) * (x1 - x0));

                      float ux = u0 * ax + u1 * bx + u2 * cx;
                      float vx = v0 * ax + v1 * bx + v2 * cx;

                      float uy = u0 * ay + u1 * by + u2 * cy;
                      float vy = v0 * ay + v1 * by + v2 * cy;

                      SampleParams sp = SampleParams();
                      sp.p_uv = uv;
                      sp.p_dx_uv = Vector2D(ux, vx) - uv;
                      sp.p_dy_uv = Vector2D(uy, vy) - uv;
                      sp.lsm = lsm;
                      sp.psm = psm;

                      Color col = tex.sample(sp);
                      rasterize_point(x, y, k, col);


                      //if (lsm == L_ZERO) {
                      //    
                      //    if (psm == P_NEAREST) {
                      //        rasterize_point(x, y, k, tex.sample_nearest(uv, 0));
                      //    }
                      //    else if (psm == P_LINEAR) {
                      //        rasterize_point(x, y, k, tex.sample_bilinear(uv, 0));
                      //    }
                      //}


                      


                    
                        

                      
                      

                  }
              }


          }

      }


  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;

    
    this->sample_buffer.resize(width * height * rate, Color::White);
  }
  

  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support
      

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
          float r = 0;
          float g = 0;
          float b = 0;

          for (int k = 0; k < sample_rate; k += 1) {
              Color col = sample_buffer[y * width * sample_rate + x * sample_rate + k];
              
              r += col.r;
              g += col.g;
              b += col.b;
              
          }
          r /= sample_rate;
          g /= sample_rate;
          b /= sample_rate;
          Color new_col = Color(r, g, b);


        //Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&new_col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
