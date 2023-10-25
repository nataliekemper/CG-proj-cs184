#include "rasterizer.h"

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


    int root = (int) sqrt(sample_rate);

    for (int i = root*y; i < root*(y+1); i++)
      for (int j = root*x; j < root*(x+1); j++)
        sample_buffer[i * (root*width) + j] = c;
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

    fill_pixel(sx, sy, color);
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
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }

  // Checks whether point is within the half plane of a line defined by:
  // P0 = (x0, y0)
  // P1 = (x1, y1)
  // point = (px, py)
  // Returns true if the point within or on the line
  double RasterizerImp::lineEquation(float x0, float y0,
                                   float x1, float y1,
                                   float px, float py) {

      Vector2D P0 = Vector2D((double) x0, (double) y0);
      Vector2D P1 = Vector2D((double) x1, (double) y1);
      Vector2D target = Vector2D((double) px, (double) py);
      // Normal line of P0 -> P1
      Vector2D L = P1 - P0;
      Vector2D N = Vector2D(-(L.y), L.x);
      Vector2D V = target - P0;
      double inside = dot(N, V);
      return inside;
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
      int root = (int) sqrt(sample_rate);
      
      // imin = minimum starting x, rounded to the nearest x+0.5
      float imin = max(0.0, round(min(x0, min(x1, x2)) - 0.5)) + 0.5;
      // imax = maximum x of triangle or window width, whichever comes first
      float imax = min((float)root*width, max(x0, max(x1, x2)));

      float jmin = max(0.0, round(min(y0, min(y1, y2)) - 0.5)) + 0.5;
      float jmax = min((float)root*height, max(y0, max(y1, y2)));


      for (float i = imin; i < imax; i++) {
        for (float j = jmin; j < jmax; j++) {
              double A = lineEquation(x0*root, y0*root, x1*root, y1*root, i, j);
              double B = lineEquation(x1*root, y1*root, x2*root, y2*root, i, j);
              double C = lineEquation(x2*root, y2*root, x0*root, y0*root, i, j);
              if ((A >= 0 && B >= 0 && C >= 0) || (A <= 0 && B <= 0 && C <= 0)) {
                  int i_int = (int) floor(i);
                  int j_int = (int) floor(j);
                  if (i_int > 0 && j_int > 0 && i_int < root*width && j_int < root*height) {
                    sample_buffer[j_int * (root*width) + i_int] = color;
                  }
              }
          }
      }
    // TODO: Task 2: Update to implement super-sampled rasterization
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle

    int root = (int) sqrt(sample_rate);
    // Pre-scale all the coordinates to the high-resolution sample bufffer
    x0 *= root;
    y0 *= root;
    x1 *= root;
    y1 *= root;
    x2 *= root;
    y2 *= root;
      
    // imin = minimum starting x, rounded to the nearest x+0.5
    float imin = max(0.0, round(min(x0, min(x1, x2)) - 0.5)) + 0.5;
    // imax = maximum x of triangle or window width, whichever comes first
    float imax = min((float)root*width, max(x0, max(x1, x2)));

    float jmin = max(0.0, round(min(y0, min(y1, y2)) - 0.5)) + 0.5;
    float jmax = min((float)root*height, max(y0, max(y1, y2)));


    for (float i = imin; i < imax; i++) {
        for (float j = jmin; j < jmax; j++) {
            double A = lineEquation(x0, y0, x1, y1, i, j);
            double B = lineEquation(x1, y1, x2, y2, i, j);
            double C = lineEquation(x2, y2, x0, y0, i, j);
            if ((A >= 0 && B >= 0 && C >= 0) || (A <= 0 && B <= 0 && C <= 0)) {

                double alpha = ( -(i - x1)*(y2 - y1) + (j - y1)*(x2 - x1) )
                /( -(x0 - x1)*(y2 - y1) + (y0 - y1) * (x2 - x1));
                double beta = ( -(i - x2)*(y0 - y2) + (j - y2)*(x0 - x2) )
                /( -(x1 - x2)*(y0 - y2) + (y1 - y2) * (x0 - x2));
                double gamma = 1 - alpha - beta;

                Color color = alpha * c0 + beta * c1 + gamma * c2;

                int i_int = (int) floor(i);
                int j_int = (int) floor(j);
                if (i_int > 0 && j_int > 0 && i_int < root*width && j_int < root*height) {
                    sample_buffer[j_int * (root*width) + i_int] = color;
          }
        }
      }
    }
  }


  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
      int root = (int) sqrt(sample_rate);
      // Pre-scale all the coordinates to the high-resolution sample bufffer
      x0 *= root;
      y0 *= root;
      x1 *= root;
      y1 *= root;
      x2 *= root;
      y2 *= root;

      Vector2D T0 = Vector2D(u0, v0);
      Vector2D T1 = Vector2D(u1, v1);
      Vector2D T2 = Vector2D(u2, v2);
      Color color = Color(1, 0, 0);

      // imin = minimum starting x, rounded to the nearest x+0.5
      float imin = max(0.0, round(min(x0, min(x1, x2)) - 0.5)) + 0.5;
      // imax = maximum x of triangle or window width, whichever comes first
      float imax = min((float)root*width, max(x0, max(x1, x2)));

      float jmin = max(0.0, round(min(y0, min(y1, y2)) - 0.5)) + 0.5;
      float jmax = min((float)root*height, max(y0, max(y1, y2)));

      for (float i = imin; i < imax; i++) {
        for (float j = jmin; j < jmax; j++) {
          double A = lineEquation(x0, y0, x1, y1, i, j);
          double B = lineEquation(x1, y1, x2, y2, i, j);
          double C = lineEquation(x2, y2, x0, y0, i, j);
          if ((A >= 0 && B >= 0 && C >= 0) || (A <= 0 && B <= 0 && C <= 0)) {

            SampleParams params;
            params.p_uv = barycentric(x0, y0, x1, y1, x2, y2, i, j, T0, T1, T2);
            params.p_dx_uv = barycentric(x0, y0, x1, y1, x2, y2, i+1, j, T0, T1, T2);
            params.p_dy_uv = barycentric(x0, y0, x1, y1, x2, y2, i, j+1, T0, T1, T2);
            params.psm = this->psm;
            params.lsm = this->lsm;

            color = tex.sample(params);

            int i_int = (int) floor(i);
            int j_int = (int) floor(j);
            if (i_int > 0 && j_int > 0 && i_int < root*width && j_int < root*height) {
                sample_buffer[j_int * (root*width) + i_int] = color;
            }
          }
        }
      }
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
  }

  Vector2D RasterizerImp::barycentric(float x0, float y0, float x1, float y1, float x2, float y2,
                       float i, float j, Vector2D T0, Vector2D T1, Vector2D T2) {
    double alpha = ( -(i - x1)*(y2 - y1) + (j - y1)*(x2 - x1) )
      /( -(x0 - x1)*(y2 - y1) + (y0 - y1) * (x2 - x1));;
    double beta = ( -(i - x2)*(y0 - y2) + (j - y2)*(x0 - x2) )
      /( -(x1 - x2)*(y0 - y2) + (y1 - y2) * (x0 - x2));
    double gamma = 1 - alpha - beta;
    return alpha * T0 + beta * T1 + gamma * T2;
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;


    this->sample_buffer.resize(width * height * sample_rate, Color::White);
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

    int root = (int) sqrt(sample_rate);

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color avg = Color(0, 0, 0);

        for (int i = root*y; i < root*(y+1); i++) {
          for (int j = root*x; j < root*(x+1); j++) {
            avg += sample_buffer[i * (width*root) + j];
          }
        }

        Color col = avg * (1.0/sample_rate);

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
