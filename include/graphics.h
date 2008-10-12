#ifndef ARIADNE_GRAPHICS_H
#define ARIADNE_GRAPHICS_H

#define HAVE_GTK_H

namespace Ariadne {

class Box;
class Polytope;

class Graphic {
  public:
    ~Graphic();
    Graphic();
    void set_bounding_box(const Box&);
    void plot(const Box&);
    void plot(const Polytope&);
    void clear();
    void display();
    void write(const char* filename);
  public:
    class Impl;
  private:
    Impl* _impl;
};

Graphic& operator<<(Graphic& g, const Box& bx);

} // namespace Ariadne

#endif // ARIADNE_GRAPHICS_H
