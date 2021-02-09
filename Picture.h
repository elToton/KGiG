#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

using byte = unsigned char;

struct Point {
    float x;
    float y;
};

struct Pixel {
    unsigned char A;
    unsigned char B;
    unsigned char C;
};
enum ColorSpace {
    RGB,
    HSL,
    HSV,
    YCbCr_601,
    YCbCr_709,
    YCoCg,
    CMY
};

class Picture{
private:
    int Width, Height, ColourDepth, Type;
    byte *PictureData;
    const double OrderedMatrix[8][8] ={
            {0.0 / 64.0, 48.0 / 64.0, 12.0 / 64.0, 60.0 / 64.0, 3.0 / 64.0, 51.0 / 64.0, 15.0 / 64.0, 63.0 / 64.0},
            {32.0 / 64.0, 16.0 / 64.0, 44.0 / 64.0, 28.0 / 64.0, 35.0 / 64.0, 19.0 / 64.0, 47.0 / 64.0, 31.0 / 64.0},
            {8.0 / 64.0, 56.0 / 64.0, 4.0 / 64.0, 52.0 / 64.0, 11.0 / 64.0, 59.0 / 64.0, 7.0 / 64.0, 55.0 / 64.0},
            {40.0 / 64.0, 24.0 / 64.0, 36.0 / 64.0, 20.0 / 64.0, 43.0 / 64.0, 27.0 / 64.0, 39.0 / 64.0, 23.0 / 64.0},
            {2.0 / 64.0, 50.0 / 64.0, 14.0 / 64.0, 62.0 / 64.0, 1.0 / 64.0, 49.0 / 64.0, 13.0 / 64.0, 61.0 / 64.0},
            {34.0 / 64.0, 18.0 / 64.0, 46.0 / 64.0, 30.0 / 64.0, 33.0 / 64.0, 17.0 / 64.0, 45.0 / 64.0, 29.0 / 64.0},
            {10.0 / 64.0, 58.0 / 64.0, 6.0 / 64.0, 54.0 / 64.0, 9.0 / 64.0, 57.0 / 64.0, 5.0 / 64.0, 53.0 / 64.0},
            {42.0 / 64.0, 26.0 / 64.0, 38.0 / 64.0, 22.0 / 64.0, 41.0 / 64.0, 25.0 / 64.0, 37.0 / 64.0, 21.0 / 64.0}
    };
    const double JarvisJudiceNinkeMatrix[3][5] = {
            {0, 0, 0, 7, 5},
            {3, 5, 7, 5, 3},
            {1, 3, 5, 3, 1}
    };
    const double SierraMatrix[3][5] = {
            {0, 0, 0, 5, 3},
            {2, 4, 5, 4, 2},
            {0, 2, 3, 2, 0}
    };
    const double AtkinsonMatrix[3][5] = {
            {0, 0, 0, 1, 1},
            {0, 1, 1, 1, 0},
            {0, 0, 1, 0, 0}
    };
    const double HalftoneMatrix[4][4] = {
            {13.0 / 16.0, 11.0 / 16.0, 4.0 / 16.0, 8.0 / 16.0},
            {6.0 / 16.0, 0, 3.0 / 16.0, 15.0 / 16.0},
            {14.0 / 16.0, 1.0 / 16.0, 2.0 / 16.0, 7.0 / 16.0},
            {9.0 / 16.0, 5.0 / 16.0, 10.0 / 16.0, 12.0 / 16.0}
    };
    double RGBCorrection(double value);
    double revCorrection(double value);
    ColorSpace currentCS;
    void ConvertRGB();
    Pixel *image;
public:
    Picture();
    float correction;
    explicit Picture(char *path);
    static void WriteBinary(char* path, char *vector, int countbytes);
    bool Read(char* path);
    void Inverse();
    void Mirror(const char string[2]);
    void Rotation(const char string[2]);
    void Output(const char *string);
    void DrawLine(Point start, Point end, byte colourdepth, float thickness, float gamma);
    void DrawLine(float x0, float y0, float x1, float y1, byte brightness, float thickness, float gamma);
    void DrawPoint(int x, int y, double depth, byte color, double gamma);
    void WithOutDithering(unsigned int bits);
    void Ordered(unsigned int bits);
    void Random(unsigned int bits);
    void FloydSteinberg(unsigned int bits);
    void JarvisJudiceNinke(unsigned int bits);
    void Sierra(unsigned int bits);
    void Atkinson(unsigned int bits);
    void Halftone(unsigned int bits);
    void Grad();
    void SetCS(ColorSpace);
    void ReadOneFile(string filename);
    void ReadThreeFiles(string filename);
    void WriteToOneFile(string filename);
    void WriteToThreeFiles(string filename);
    void ConvertAny(ColorSpace convert);
};

