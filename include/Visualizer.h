#ifndef VISUALIZER_H
#define VISUALIZER_H

enum FieldType{E,H};

class Visualizer:public mglDraw{
public:
  Visualizer(Solver* solver);
  Solver* solver;
  int framedelay;
  FieldType fieldtype;
  int Draw(mglGraph *gr);
  void makeGif(const char* file="Output.gif",const uint Nframes=20);
  void Test(const char* file="Output.gif",const uint Nframes=20);
  int run();
};

#endif // VISUALIZER_H
