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
  void RasterizerImp::fill_pixel(float x, float y, Color c, bool supersampled) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
      int sq_sr = sqrt(this->sample_rate);
      float sx = (x*sq_sr);
      float sy = (y*sq_sr);
      if (supersampled){
          for (int i=0;i<sq_sr;i++){
              for (int j=0;j<sq_sr;j++){
                  // check bounds
                  sample_buffer[(size_t)((sy+j)*width*sq_sr+(sx+i))] = c;
              }
          }
      } else {
          sample_buffer[(size_t)(sy * width*sq_sr + sx)] = c;
      }
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

    fill_pixel((double)sx, (double)sy, color, true);
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

//Perform line test for triangle
    bool line_test(float pt_x, float pt_y, float x0, float y0, float x1, float y1, float x2, float y2){
        float l1 = -(pt_x-x0)*(y1-y0) + (pt_y-y0)*(x1-x0);
        float l2 = -(pt_x-x1)*(y2-y1) + (pt_y-y1)*(x2-x1);
        float l3 = -(pt_x-x2)*(y0-y2) + (pt_y-y2)*(x0-x2);
        return (l1 >= 0 && l2 >= 0 && l3 >= 0) || (l1 <= 0 && l2 <= 0 && l3 <= 0);
    }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
//      int max_x = (int)max({x0,x1,x2});
//      int max_y = (int)max({y0,y1,y2});
//      int min_x = (int)min({x0,x1,x2});
//      int min_y = (int)min({y0,y1,y2});
//      for (int i = min_y; i <= max_y; i++) {
//          for (int j = min_x; j <= max_x; j++){
//              double pt[] = {(j+0.5), (i+0.5)};
//              if (line_test(pt[0], pt[1], x0, y0, x1, y1, x2, y2)) {
//                rasterize_point(pt[0], pt[1], color);
//              }
//          }
//      }
    // TODO: Task 2: Update to implement super-sampled rasterization
      int sq_sr = sqrt(this->sample_rate);
      int max_x = (int)max({x0,x1,x2})*sq_sr;
      int max_y = (int)max({y0,y1,y2})*sq_sr;
      int min_x = (int)min({x0,x1,x2})*sq_sr;
      int min_y = (int)min({y0,y1,y2})*sq_sr;
      for (int i = min_y; i <= max_y; i++) {
          for (int j = min_x; j <= max_x; j++){
              float pt[] = {(float)(j+0.5), (float)(i+0.5)};
              if (line_test(pt[0], pt[1], x0*sq_sr, y0*sq_sr, x1*sq_sr, y1*sq_sr, x2*sq_sr, y2*sq_sr)) {
                  // check bounds
                  if (j < 0 || j >= width*sq_sr) return;
                  if (i < 0 || i >= height*sq_sr) return;
                  fill_pixel((float)j/sq_sr, (float)i/sq_sr, color, false);
//                rasterize_point(pt[0]/sqrt(sample_rate), pt[1]/sqrt(sample_rate), color);
              }
          }
      }
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
//      float x = (x0 + x1 + x2)/3;
//      float y = (y0 + y1 + y2)/3;
      
      //rasterize
      int sq_sr = sqrt(this->sample_rate);
      int max_x = (int)max({x0,x1,x2})*sq_sr;
      int max_y = (int)max({y0,y1,y2})*sq_sr;
      int min_x = (int)min({x0,x1,x2})*sq_sr;
      int min_y = (int)min({y0,y1,y2})*sq_sr;
      for (int y = min_y; y <= max_y; y++) {
          for (int x = min_x; x <= max_x; x++){
              float pt[] = {(float)(x+0.5), (float)(y+0.5)};
              if (line_test(pt[0], pt[1], x0*sq_sr, y0*sq_sr, x1*sq_sr, y1*sq_sr, x2*sq_sr, y2*sq_sr)) {
                  float a = (-(x - x1)*(y2 - y1) + (y - y1)*(x2 - x1))/(-(x0 - x1)*(y2-y1) + (y0 - y1)*(x2 - x1));
                  float b = (-(x - x2)*(y0 - y2) + (y - y2)*(x0 - x2))/(-(x1 - x2)*(y0-y2) + (y1 - y2)*(x0 - x2));
                  float g = 1 - a - b;
                  Color col = a*c0 + b*c1 + g*c2;
                  // check bounds
                  if (x < 0 || x >= width*sq_sr) return;
                  if (y < 0 || y >= height*sq_sr) return;
                  fill_pixel((float)x/sq_sr, (float)y/sq_sr, col, false);
              }
          }
      }

  }




  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
      int sq_sr = sqrt(this->sample_rate);
      int max_x = (int)max({x0,x1,x2})*sq_sr;
      int max_y = (int)max({y0,y1,y2})*sq_sr;
      int min_x = (int)min({x0,x1,x2})*sq_sr;
      int min_y = (int)min({y0,y1,y2})*sq_sr;
      for (int y = min_y; y <= max_y; y++) {
          for (int x = min_x; x <= max_x; x++){
              float pt[] = {(float)(x+0.5), (float)(y+0.5)};
              if (line_test(pt[0], pt[1], x0*sq_sr, y0*sq_sr, x1*sq_sr, y1*sq_sr, x2*sq_sr, y2*sq_sr)) {
                  float a = (-(x - x1)*(y2 - y1) + (y - y1)*(x2 - x1))/(-(x0 - x1)*(y2-y1) + (y0 - y1)*(x2 - x1));
                  float b = (-(x - x2)*(y0 - y2) + (y - y2)*(x0 - x2))/(-(x1 - x2)*(y0-y2) + (y1 - y2)*(x0 - x2));
                  float g = 1 - a - b;
                  double u = a*u0 + b*u1 + g*u2;
                  double v = a*v0 + b*v1 + g*v2;
                  SampleParams params;
                  params.psm = this->psm;
                  params.lsm = this->lsm;
                  //change:
                  params.p_uv = {u*tex.width, v*tex.height}; // ?
                  Color col = Color(1,0,1);
                  if (params.psm == P_NEAREST){
                      col = tex.sample_nearest(params.p_uv);
                  } else {
                      col = tex.sample_bilinear(params.p_uv);
                  }
                  // check bounds
                  if (x < 0 || x >= width*sq_sr) return;
                  if (y < 0 || y >= height*sq_sr) return;
                  fill_pixel((float)x/sq_sr, (float)y/sq_sr, col, false);
              }
          }
      }
      
      
     
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle




  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;
      clear_buffers(); //added


//    this->sample_buffer.resize(width * height, Color::White); //changed
      this->sample_buffer.resize(width * height * this->sample_rate, Color::White);
      
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;
    clear_buffers(); //added


//    this->sample_buffer.resize(width * height, Color::White); //changed
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
      int sq_sr = sqrt(this->sample_rate);

    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
//          Color scol = sample_buffer[y * width + x];
          size_t sx = sq_sr * x;
          size_t sy = sq_sr * y;
          float sum_r = 0.0;
          float sum_g = 0.0;
          float sum_b = 0.0;
          //sum rgb values over kernel
          for (int i=0;i<sq_sr;i++){
              for (int j=0;j<sq_sr;j++){
                  Color c = sample_buffer[(sy+j)*width*sq_sr+(sx+i)];
                  sum_r += c.r;
                  sum_g += c.g;
                  sum_b += c.b;
              }
          }
          // set color to avg of rgb values (divide by sample_rate)
          Color col = Color(sum_r/(this->sample_rate),sum_g/(this->sample_rate),sum_b/(this->sample_rate));
//        Color col = sample_buffer[y*sq_sr * width + x*sq_sr];
        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
