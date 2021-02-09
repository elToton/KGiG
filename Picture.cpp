#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include "Picture.h"
#include <fstream>
#include <algorithm>

#define MAX_VALUE 255

using namespace std;

using byte = unsigned char;



Picture::Picture(char *path){
    error_code ec{};
    if (ec != error_code{}) {
        cout << ec.message() << endl;
        exit(1);
    }
    Read(path);
    if (Type != 5 && Type != 6) {
        cout << "wrong type" << endl;
        exit(1);
    }
}
//запись метаданных
void Picture::WriteBinary(char* path, char *vector, int countbytes){
    FILE *outfile = fopen(path, "w");
    if (outfile){
        fwrite(vector, 1, countbytes, outfile);
        fclose(outfile);
    } else {
        cout << "error output file" << endl;
        exit(1);
    }
};

bool Picture::Read(char* path){
    FILE *infile = fopen(path, "rb");
    if (!infile){
        cout << "error input file" << endl;
        return -1;
    }
    //считывание хэдера
    fscanf(infile, "P%d\n%d %d\n%d\n", &this -> Type, &this -> Width, &this -> Height, &this -> ColourDepth);
    //считывание данных
    if (Type == 5) {
        PictureData = static_cast<byte *>(malloc(sizeof(byte) * Width * Height));
        fread(PictureData,1,Width*Height,infile);
    }else{
        //так как rgb
        PictureData = static_cast<byte *>(malloc(sizeof(byte) * Width * Height*3));
        fread(PictureData,1,Width*Height*3,infile);
    }
    return true;
}

//запись метаданных в файл
void Picture::Output(const char *path) {
    stringstream ss;
    ss << Width;
    string WidthS;
    ss >> WidthS;
    ss.clear();
    ss << Height;
    string HeightS;
    ss >> HeightS;
    ss.clear();
    ss << ColourDepth;
    string ColourDepthS;
    ss >> ColourDepthS;
    ss.clear();
    //выделение памяти в файл для хэдера и данных в байтах
    int BuffSize;
    if (Type == 5) {
        BuffSize = sizeof(Width) + sizeof(Height) + sizeof(ColourDepth) + Width * Height + 6;
    } else if (Type == 6){
        BuffSize = sizeof(Width) + sizeof(Height) + sizeof(ColourDepth) + Width * Height * 3 + 6;
    }
    char Buff[BuffSize];
    //запись метаданных
    Buff[0] = 'P';
    Buff[1] = char(Type + '0');
    Buff[2] = '\n';
    int i = 3;

    for (auto c : WidthS) {
        Buff[i] = c;
        i += 1;
    }
    Buff[i] = ' ';
    i += 1;
    for (auto c : HeightS) {
        Buff[i] = c;
        i += 1;
    }
    Buff[i] = '\n';
    i += 1;
    for (auto c : ColourDepthS) {
        Buff[i] = c;
        i += 1;
    }
    Buff[i] = '\n';
    i += 1;
    //запись данных
    if(Type == 5) {
        for (int j = 0; j < Width * Height; j++) {
            Buff[i] = PictureData[j];
            i += 1;
        }
    }else{
        for (int j = 0; j < Width * Height * 3; j++) {
            Buff[i] = PictureData[j];
            i += 1;
        }
    }
    WriteBinary(const_cast<char *>(path), Buff, sizeof(Buff));
}

//реверс цвета
void Picture::Inverse(){
    int countbytes;
    if (Type == 5){
        countbytes = Width*Height;
    }else{
        countbytes=Width*Height*3;
    }
    for (int i = 0; i < countbytes; i++) {
        PictureData[i] = ~PictureData[i];
    }
}
//отражение по гор/вер
void Picture::Mirror(const char *direction){
    //по горизонтали
    if (direction[0] == 'H'){
        //для моно
        if (Type == 5){
            uint_fast64_t i,j;
            for (i = 0; i < Height; i++) {
                for (j = 0; j < Width/2; j++) {
                    swap(PictureData[i * Width + j], PictureData[i * Width + Width - 1 - j]);
                }
            }
        } else{//для ргб
            uint_fast64_t i,j;
            for (i = 0; i < Height; i++){
                for (j = 0; j < Width*3/2; j+=3){
                    swap(PictureData[i * Width * 3 + j], PictureData[i * Width * 3 + Width * 3 - j - 3]);
                    swap(PictureData[i * Width * 3 + j + 1], PictureData[i * Width * 3 + Width * 3 - j - 2]);
                    swap(PictureData[i * Width * 3 + j + 2], PictureData[i * Width * 3 + Width * 3 - j - 1]);
                }
            }
        }
    }
    //по вертикали
    if (direction[0] == 'V'){
        //для моно
        if (Type == 5){
            uint_fast64_t i,j;
            for (i = 0; i < Width; i++){
                for (j = 0; j < Height/2; j++){
                    swap(PictureData[j * Width + i], PictureData[(Height - 1 - j) * Width + i]);
                }
            }
        } else{//для ргб
            uint_fast64_t i,j;
            for (i = 0; i < Width*3; i++){
                for (j = 0; j < Height/2; j++){
                    swap(PictureData[j * Width * 3 + i], PictureData[(Height - 1 - j) * Width * 3 + i]);
                }
            }
        }
    }
}
//разворот
void Picture::Rotation(const char *direction){
    //по часовой
    if (direction[0] == 'R'){
        //для моно
        if (Type == 5){
            byte* TurnedMetaData;
            uint_fast64_t TurnedWidth, TurnedHeight;
            TurnedHeight = Width;
            TurnedWidth = Height;
            TurnedMetaData = static_cast<byte *>(malloc(Width * Height));
            for (uint_fast64_t i = 0; i < Height; i++){
                for (uint_fast64_t j = 0; j < Width; j++) {
                    TurnedMetaData[(TurnedWidth - i - 1) + j * TurnedWidth] = PictureData[i * Width + j];
                }
            }
            Width = TurnedWidth;
            Height = TurnedHeight;
            PictureData = TurnedMetaData;
        } else{//для ргб
            byte* TurnedMetaData;
            uint64_t TurnedWidth, TurnedHeight;
            TurnedHeight = Width;
            TurnedWidth = Height;
            TurnedMetaData = static_cast<byte *>(malloc(Width * Height*9));
            for (uint64_t i = 0; i < Height; i++) {
                for (uint64_t j = 0; j < Width; j++) {
                    TurnedMetaData[j * 3 + i * TurnedWidth * 3] = PictureData[(Height - 1 - j) * Width * 3 + i * 3];
                    TurnedMetaData[j * 3 + 1 + i * TurnedWidth * 3] = PictureData[(Height - 1 - j) * Width * 3 + i * 3 + 1];
                    TurnedMetaData[j * 3 + 2 + i * TurnedWidth * 3] = PictureData[(Height - 1 - j) * Width * 3 + i * 3 + 2];
                }
            }
            Width = TurnedWidth;
            Height = TurnedHeight;
            PictureData = TurnedMetaData;
        }
    }
    //против часовой
    if (direction[0] == 'L'){
        //для моно
        if (Type == 5){
            byte* TurnedData;
            uint64_t TurnedWidth, TurnedHeight;
            TurnedHeight = Width;
            TurnedWidth = Height;
            TurnedData = static_cast<byte *>(malloc(Width * Height));
            for (uint64_t i = 0; i < Height; i++){
                for (uint64_t j = 0; j < Width; j++){
                    TurnedData[j * TurnedWidth + i] = PictureData[(i + 1) * Width - 1 - j];
                }
            }
            Width = TurnedWidth;
            Height = TurnedHeight;
            PictureData = TurnedData;
        } else{//для ргб
            byte* TurnedData;
            uint64_t TurnedWidth, TurnedHeight;
            TurnedHeight = Width;
            TurnedWidth = Height;
            TurnedData = static_cast<byte *>(malloc(Width * Height * 9));
            for (uint64_t i = 0; i < Height; i++) {
                for (uint64_t j = 0; j < Width; j++) {
                    TurnedData[(TurnedHeight - 1 - j) * TurnedWidth * 3 + i * 3] = PictureData[i * Width * 3 + j * 3];
                    TurnedData[(TurnedHeight - 1 - j) * TurnedWidth * 3 + i * 3 + 1] = PictureData[i * Width * 3 + j * 3 + 1];
                    TurnedData[(TurnedHeight - 1 - j) * TurnedWidth * 3 + i * 3 + 2] = PictureData[i * Width * 3 + j * 3 + 2];
                }
            }
            Width = TurnedWidth;
            Height = TurnedHeight;
            PictureData = TurnedData;
        }
    }
}

//..................................................................................................................

void Picture::DrawLine(Point start, Point end, byte color, float thickness, float gamma) {
    if (thickness <= 0) {
        cout << "wrong thickness" << endl;
        return;
    }
    int px;
    int py;
    if (start.x > end.x) swap(start.x, end.x);
    if (start.y > end.y) swap(start.y, end.y);
    bool steep = abs(end.y - start.y) > abs(end.x - start.x);
    if (steep) {
        swap(start.x, start.y);
        swap(end.x, end.y);
    }
    auto intPart = [](double x) -> int { return (int) x; };
    auto plot = [&](int x, int y, double intensity) -> void {
        if (steep) DrawPoint(y, x, 1.0 - intensity, color, gamma);
        else DrawPoint(x, y, 1.0 - intensity, color, gamma);
    };
    auto distance = [](Point a, Point b) -> double {
        return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
    };

    if (start.y == end.y){
        double y = end.y + thickness / 2;
        py = floor(y) - (floor(thickness)/2);
        int n = floor(thickness);
        for(int i = 0; i < n; i++){
            for (px = start.x; px < end.x; px++) {
                plot(px, py, y - start.y + thickness / 2);
            }
            py++;

        }
    }

    if (start.x == end.x){
        double x = end.x + thickness / 2;
        px = floor(x) - (floor(thickness)/2);
        int n = floor(thickness);
        for(int i = 0; i < n; i++){
            for (py = start.y; py < end.y; py++) {
                plot(px, py, x - start.x + thickness / 2);
            }
            px++;

        }
    }

    if ((start.x != end.x) and (start.y != end.y)) {
        double dx = end.x - start.x;
        double dy = end.y - start.y;
        double gradient = dy / dx;

        double y = start.y + gradient * (round(start.x) - start.x);

        for (int x = round(start.x); x <= round(end.x); x++) {
            for (int plotY = intPart(y - (thickness - 1) / 2);
                 plotY <= intPart(y - (thickness - 1) / 2 + thickness); plotY++) {
                plot(x, plotY, min(1.0, (thickness + 1.0) / 2.0 - fabs(y - plotY)));
            }
            y += gradient;
        }



        Point plotStart = {round(start.x), round(start.y)};
        for (int plotX = round(start.x) - thickness / 2; plotX < round(start.x); plotX++) {
            y = start.y + gradient * (plotX - start.x);
            for (int plotY = int(y - (thickness - 1) / 2.0);
                 plotY <= int(y - (thickness - 1) / 2.0 + thickness); plotY++) {
                plot(plotX, plotY, min(1.0, (thickness + 0.5) / 2.0 -
                                            distance({(float) plotX, (float) plotY}, {plotStart.x, plotStart.y})));
            }
        }

        Point plotEnd = {round(end.x), round(end.y)};
        for (int plotX = round(end.x) + 1; plotX <= round(end.x) + thickness / 2; plotX++) {
            y = start.y + gradient * (plotX - start.x);
            for (int plotY = int(y - (thickness - 1) / 2.0);
                 plotY <= int(y - (thickness - 1) / 2.0 + thickness); plotY++) {
                plot(plotX, plotY, min(1.0, (thickness + 0.5) / 2.0 -
                                            distance({(float) plotX, (float) plotY}, {plotEnd.x, plotEnd.y})));
            }
        }
        return;
    }
}

void Picture::DrawLine(float x0, float y0, float x1, float y1, byte color, float thickness, float gamma) {
    DrawLine({x0, y0}, {x1, y1}, color, thickness, gamma);
}

void Picture::DrawPoint(int x, int y, double depth, byte color, double gamma) {
    depth = max(min(depth, 1.0), 0.0);
    if (depth < 0 || depth > 1) {
        cout << "wrong depth" << endl;
        return;
    }
    if (y < 0 || y >= Height || x < 0 || x >= Width) {
        cout << "wrong coords" << endl;
        return;
    }
    double k = color / (double) PictureData[Width * y + x];
    PictureData[Width * y + x] *= pow((1 - k) * depth + k, 1.0 / gamma);
}

//......................................................................................................................

void Picture::WithOutDithering(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    vector<double> error;
    for (int i = 0; i < Width * Height; i++) {
        buff = revCorrection(PictureData[i]) / (double) 255;
        buff *= value - 1;
        buff = round(buff);
        PictureData[i] = round(RGBCorrection(buff * (255 / (value - 1))));
    }
}
/*
void Picture::Ordered(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    vector<double> error;
    for (int i = 0; i < Height*3; i+=3){
        for (int j = 0; j < Width*3; j += 3) {
            for (int k = 0; k < 3; k++){
                buff = (revCorrection(PictureData[i * Width + j + k]) +
                        (255 / (bits)) * (OrderedMatrix[(i/3) % 8][(j/3)  % 8] - 0.5)) /
                       255;
                if (buff < 0) buff = 0;
                buff *= value - 1;
                buff = round(buff);
                PictureData[i * Width + j + k] = round(RGBCorrection(buff * (255 / (value - 1))));
            }
        }
    }
}*/

void Picture::Ordered(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    vector<double> error;
    for (int i = 0; i < Height; i++){
        for (int j = 0; j < Width; j++){
            buff = (revCorrection(PictureData[i * Width + j]) + (255 / (bits)) * (OrderedMatrix[i % 8][j % 8] - 0.5)) / 255;
            if (buff < 0) buff = 0;
            buff *= value - 1;
            buff = round(buff);
            PictureData[i * Width + j] = round(RGBCorrection(buff * (255 / (value - 1))));
        }
    }
}


void Picture::Random(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    vector<double> error;
    srand(666);
    for (int i = 0; i < Height*Width; i++) {
        buff = (revCorrection(PictureData[i])) + 255/static_cast<double>(rand()/32767. - 0,75)/static_cast<double>(255);
        if (buff < 0) buff = 0;
        buff *= value;
        buff = round(buff);
        PictureData[i] = round(RGBCorrection(buff * (255 / static_cast<double>(value - 1))));


        /*for (int j = 0; j < Width; j++) {
            buff = (revCorrection(PictureData[i * Width + j]) + (255 / (bits)) * ((double) rand() / 32767.0 - 0.75)) / (double) 255;
            if (buff < 0) buff = 0;
            double noise =  (double) rand() / RAND_MAX - 0.75;
            value += noise / bits;
            value = min(max(value, 0), 0);
            buff *= value;
            buff = round(buff);
            PictureData[i * Width + j] = round(RGBCorrection(buff * (255 / (value - 1))));
        }*/
    }
}

  void Picture::FloydSteinberg(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    double kerror;
    vector<double> error;
    error.resize(Width*Height, 0);
    for (int i = 0; i < Height; i++){
        for (int j = 0; j < Width; j++){
            buff = (revCorrection(PictureData[i * Width + j]) + error[i * Width + j]) / (double)255;
            buff *= (value - 1);
            buff = round(buff);
            buff *= 255 / (value - 1);
            kerror = PictureData[i * Width + j] + error[i * Width + j] - buff;
            PictureData[i * Width + j] = buff;
            if (j + 1 < Width) error[i * Width + j + 1] += kerror * (7.0 / 16.0);
            if (i + 1 < Height){
                if (j + 1 < Width) error[(i + 1) * Width + j + 1] += kerror * (1.0 / 16.0);
                error[(i + 1) * Width + j] += kerror * (5.0 / 16.0);
                if ((i - 1 > 0)&&(j - 1 > 0)) error[(i - 1) * Width + j - 1] += kerror * (3.0 / 16.0);
            }
        }
    }
}


void Picture::JarvisJudiceNinke(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    double kerror;
    vector<double> error;
    error.resize(Width * Height, 0);
    for (int i = 0; i < Height; i++){
        for (int j = 0; j < Width; j++){
            buff = (revCorrection(PictureData[i * Width + j]) + error[i * Width + j]) / (double)255;
            buff *= (value - 1);
            buff = round(buff);
            buff *= 255 / (value - 1);
            kerror = PictureData[i * Width + j] + error[i * Width + j] - buff;
            PictureData[i * Width + j] = buff;
            for (int k = 0; k <= 2; k++){
                for (int v = -2; v <= 2; v++){
                    if (i + k < Height){
                        if ((k == 0) && (v > 0)){
                            if (j + v < Width) error[(i + k) * Width + j + v] += kerror * JarvisJudiceNinkeMatrix[k][2 + v] / 48.0;
                        }
                        else{
                            if ((j + v < Width) && (j + v > 0)) error[(i + k) * Width + j + v] += kerror * JarvisJudiceNinkeMatrix[k][2 + v] / 48.0;
                        }
                    }
                }
            }
        }
    }
}

void Picture::Sierra(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    double kerror;
    vector<double> error;
    error.resize(Width * Height, 0);
    for (int i = 0; i < Height; i++){
        for (int j = 0; j < Width; j++){
            buff = (revCorrection(PictureData[i * Width + j]) + error[i * Width + j]) / (double)255;
            buff *= (value - 1);
            buff = round(buff);
            buff *= 255 / (value - 1);
            kerror = PictureData[i * Width + j] + error[i * Width + j] - buff;
            PictureData[i * Width + j] = buff;
            for (int k = 0; k <= 2; k++){
                for (int v = -2; v <= 2; v++){
                    if (i + k < Height){
                        if ((k == 0) && (v > 0)){
                            if (j + v < Width) error[(i + k) * Width + j + v] += kerror * SierraMatrix[k][2 + v] / 32.0;
                        }
                        else{
                            if ((j + v < Width) && (j + v > 0)) error[(i + k) * Width + j + v] += kerror * SierraMatrix[k][2 + v] / 32.0;
                        }
                    }
                }
            }
        }
    }
}

void Picture::Atkinson(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    double kerror;
    vector<double> error;
    error.resize(Width * Height, 0);
    for (int i = 0; i < Height; i++){
        for (int j = 0; j < Width; j++){
            buff = (revCorrection(PictureData[i * Width + j]) + error[i * Width + j]) / (double)255;
            buff *= (value - 1);
            buff = round(buff);
            buff *= 255 / (value - 1);
            kerror = PictureData[i * Width + j] + error[i * Width + j] - buff;
            PictureData[i * Width + j] = buff;
            for (int k = 0; k <= 2; k++){
                for (int v = -2; v <= 2; v++){
                    if (i + k < Height){
                        if ((k == 0) && (v > 0)){
                            if (j + v < Width) error[(i + k) * Width + j + v] += kerror * AtkinsonMatrix[k][2 + v] / 8.0;
                        }
                        else{
                            if ((j + v < Width) && (j + v > 0)) error[(i + k) * Width + j + v] += kerror * AtkinsonMatrix[k][2 + v] / 8.0;
                        }
                    }
                }
            }
        }
    }
}

void Picture::Halftone(unsigned int bits){
    int value = pow(2, bits);
    double buff;
    vector<double> error;
    for (int i = 0; i < Height; i++){
        for (int j = 0; j < Width; j++){
            buff = (revCorrection(PictureData[i * Width + j]) + (255 / (bits)) * (HalftoneMatrix[i % 4][j % 4] - 0.75)) / 255;
            if (buff < 0) buff = 0;
            buff *= value;
            buff = round(buff);
            PictureData[i * Width + j] = round(RGBCorrection(buff * (255 / (value - 1))));
        }
    }
}

/*void Picture::Grad(){
    for (int i = 0; i < Height*3; i+=3){
        for (int j = 0; j < Width*3; j+=3) {
            PictureData[i * Width + j] = RGBCorrection(((double) j / Width) * 255.0);
            PictureData[i * Width + j + 1] = RGBCorrection(((double) j / Width) * 255.0);
            PictureData[i * Width + j + 2] = RGBCorrection(((double) j / Width) * 255.0);
        }
    }
}*/


void Picture::Grad(){
    for (int i = 0; i < Height; i++){
        for (int j = 0; j < Width; j++) PictureData[i * Width + j] = RGBCorrection(((double)j / Width) * 255.0);
    }
}

double Picture::RGBCorrection(double value){
    value = value / 255;
    if (value > 1) value = 1;
    if (correction == 0){
        if (value < 0.0031308) return  value * 12.92 * 255;
        else return 255 * ((211.0 * pow(value, 0.4166) - 11.0) / 200.0);
    }
    else return 255 * pow(value, correction);
}

double Picture::revCorrection(double value){
    value = value / 255;
    if (correction == 0){
        if (value < 0.04045) return  (255 * value) / 12.92;
        else return 255 * (pow((200.0 * value + 11.0) / 211.0, 2.4));
    }
    else return 255 * pow(value, 1 / correction);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Picture::Picture() : Height(0), Width(0), image(nullptr), currentCS(ColorSpace::RGB){}

void Picture::ReadOneFile(string filename){
    ifstream input;
    input.open(filename, ios_base::in | ios_base::binary);
    string header;
    input >> header;
    input >> Width >> Height;
    int maxvalue;
    input >> maxvalue;
    char partPixel;
    char* image_char = &partPixel;
    image = new Pixel[Width * Height];
    input.read(image_char, 1);
    for (int i = 0; i < Width * Height; i++){
        input.read(image_char, 1);
        image[i].A = *image_char;
        input.read(image_char, 1);
        image[i].B = *image_char;
        input.read(image_char, 1);
        image[i].C = *image_char;
    }
    input.close();
}

void Picture::ReadThreeFiles(string filename){
    ifstream inputA, inputB, inputC;
    filename.insert(filename.size() - 4, "_1");
    inputA.open(filename, ios_base::in | ios_base::binary);
    filename.replace(filename.size() - 5, 1, "2");
    inputB.open(filename, ios_base::in | ios_base::binary);
    filename.replace(filename.size() - 5, 1, "3");
    inputC.open(filename, ios_base::in | ios_base::binary);

    string headerA, headerB, headerC;
    inputA >> headerA;
    inputB >> headerB;
    inputC >> headerC;
    inputA >> Width >> Height;
    int widthB, widthC, heightB, heightC;
    inputB >> widthB >> heightB;
    inputC >> widthC >> heightC;
    int maxvalue;
    inputA >> maxvalue;
    inputB >> maxvalue;
    inputC >> maxvalue;
    char partPixel;
    char* image_char = &partPixel;
    image = new Pixel[Width * Height];
    inputA.read(image_char, 1);
    inputB.read(image_char, 1);
    inputC.read(image_char, 1);
    for (int i = 0; i < Width * Height; i++){
        inputA.read(image_char, 1);
        image[i].A = *image_char;
        inputB.read(image_char, 1);
        image[i].B = *image_char;
        inputC.read(image_char, 1);
        image[i].C = *image_char;
    }
}

void Picture::WriteToOneFile(string filename){
    ofstream output;
    output.open(filename, ios_base::out | ios_base::binary);
    output << "P6" << '\n';
    output << Width << ' ' << Height << '\n' << MAX_VALUE << '\n';
    char* image_char = (char*)image;
    output.write(image_char, Width * Height * 3);
    output.close();
}

void Picture::WriteToThreeFiles(string filename){
    ofstream outputA, outputB, outputC;
    filename.insert(filename.size() - 4, "_1");
    outputA.open(filename, ios_base::out | ios_base::binary);
    filename.replace(filename.size() - 5, 1, "2");
    outputB.open(filename, ios_base::out | ios_base::binary);
    filename.replace(filename.size() - 5, 1, "3");
    outputC.open(filename, ios_base::out | ios_base::binary);
    outputA << "P5" << '\n';
    outputA << Width << ' ' << Height << '\n' << MAX_VALUE << '\n';
    outputB << "P5" << '\n';
    outputB << Width << ' ' << Height << '\n' << MAX_VALUE << '\n';
    outputC << "P5" << '\n';
    outputC << Width << ' ' << Height << '\n' << MAX_VALUE << '\n';
    char* part_pixel;
    for (int i = 0; i < Width * Height; i++){
        outputA.write((char*)&image[i].A, 1);
        outputB.write((char*)&image[i].B, 1);
        outputC.write((char*)&image[i].C, 1);
    }
}

void Picture::ConvertAny(ColorSpace convert){
    if (currentCS == convert) return;
    ConvertRGB();
    if (convert == RGB) return;
    double Max, Min, H, S, L, C, X, V, M;
    double Y, Cb, Cr, Co, Cg, Kr, Kg, Kb, R, G, B;
    char Ho;
    switch (convert){
        case HSL:
        case HSV:
            for (int i = 0; i < Width * Height; i++) {
                R = image[i].A / 255.0;
                G = image[i].B / 255.0;
                B = image[i].C / 255.0;
                Max = max(R, max(G, B));
                Min = min(R, min(G, B));
                V = Max;
                C = Max - Min;
                L = V - C / 2.0;
                if (C == 0) H = 0;
                else{
                    if (V == R)
                        H = (60.0) * ((G - B) / C);
                    else if (V == G)
                        H = (60.0) * (2 + (B - R) / C);
                    else if (V == B)
                        H = (60.0) * (4 + (R - G) / C);
                    else
                        H = 0;
                }
                if (convert == HSV){
                    S = (V == 0) ? 0 : C / V;
                    image[i].C = V * 255.0;
                }
                if (convert == HSL) {
                    S = ((L == 0) || (L == 1)) ? 0 : ((V - L) / min(L, 1 - L));
                    image[i].C = L * 255.0;
                }
                image[i].B = S * 255.0;
                image[i].A = (H/360.0)*255.0;

            }
            currentCS = convert;
            break;
        case YCbCr_601:
        case YCbCr_709:
            if (convert == YCbCr_601) {
                Kr = 0.299;
                Kg = 0.587;
                Kb = 0.114;
                currentCS = YCbCr_601;
            }
            else{
                Kr = 0.0722;
                Kg = 0.2126;
                Kb = 0.7152;
                currentCS = YCbCr_709;
            }
            for (int i = 0; i < Width * Height; i++){
                R = image[i].A / 255.0;
                G = image[i].B / 255.0;
                B = image[i].C / 255.0;
                Y = Kr * R + Kg * G + Kb * B;
                Cb = 0.5 * ((B - Y) / (1.0 - Kb));
                Cr = 0.5 * ((R - Y) / (1.0 - Kr));
                image[i].A = Y * 255.0;
                image[i].B = (Cb + 0.5)*255.0;
                image[i].C = (Cr + 0.5)*255.0;
            }
            break;
        case YCoCg:
            for (int i = 0; i < Width * Height; i++) {
                R = image[i].A / 255.0;
                G = image[i].B / 255.0;
                B = image[i].C / 255.0;
                Y = R / 4 + G / 2 + B / 4;
                Co = R / 2 - B / 2;
                Cg = -R / 4 + G / 2 - B / 4;
                image[i].A = Y * 255.0;
                image[i].B = (Co + 0.5) * 255.0;
                image[i].C = (Cg + 0.5) * 255.0;
            }
            currentCS = YCoCg;
            break;
        case CMY:
            for (int i = 0; i < Width * Height; i++) {
                R = image[i].A / 255.0;
                G = image[i].B / 255.0;
                B = image[i].C / 255.0;
                C = 1 - R;
                M = 1 - G;
                Y = 1 - B;
                image[i].A = C * 255.0;
                image[i].B = M  * 255.0;
                image[i].C = Y * 255.0;
            }
            currentCS = CMY;
            break;
    }
}

void Picture::ConvertRGB(){
    double H, S, L, C, H_D, X, m, R, G, B, Y, Cb, Cr, Co, Cg, M, Kr, Kg, Kb;
    if (currentCS == RGB) return;
    switch (currentCS){
        case HSL:
        case HSV:
            for (int i = 0; i < Height * Width; i++){
                H = (image[i].A / 255.0) * 360.0;
                S = image[i].B / 255.0;
                L = image[i].C / 255.0;
                H_D = H/ 60;
                if (currentCS == HSL) {
                    C = (1 - abs(2 * L - 1)) * S;
                    X = C * (1 - abs(fmod(H_D, 2) - 1));
                    m = L - C / 2.0;
                }
                else{
                    C = S * L;
                    X = C * (1.0 - abs(fmod(H_D, 2) - 1.0));
                    m = L - C;
                }
                m *= 255.0;
                if ((H_D >= 0) && (H_D <= 1)){
                    image[i].A = C * 255.0 + m;
                    image[i].B = X * 255.0 + m;
                    image[i].C = 0 + m;
                }
                if ((H_D > 1) && (H_D <= 2)){
                    image[i].A = X * 255.0 + m;
                    image[i].B = C * 255.0 + m;
                    image[i].C = 0 + m;
                }
                if ((H_D > 2) && (H_D <= 3)){
                    image[i].A = 0 + m;
                    image[i].B = C * 255.0 + m;
                    image[i].C = X * 255.0 + m;
                }
                if ((H_D > 3) && (H_D <= 4)){
                    image[i].A = 0 + m;
                    image[i].B = X * 255.0 + m;
                    image[i].C = C * 255.0 + m;
                }
                if ((H_D > 4) && (H_D <= 5)){
                    image[i].A = X * 255.0 + m;
                    image[i].B = 0 + m;
                    image[i].C = C * 255.0 + m;
                }
                if ((H_D > 5) && (H_D <= 6)){
                    image[i].A = C * 255.0 + m;
                    image[i].B = 0 + m;
                    image[i].C = X * 255.0 + m;
                }
            }
            break;
        case YCbCr_601:
        case YCbCr_709:
            if (currentCS == YCbCr_601) {
                Kr = 0.299;
                Kg = 0.587;
                Kb = 0.114;
            }
            else{
                Kr = 0.0722;
                Kg = 0.2126;
                Kb = 0.7152;
            }
            for (int i = 0; i < Width * Height; i++) {
                Y = image[i].A / 255.0;
                Cb = (image[i].B / 255.0) - 0.5;
                Cr = (image[i].C / 255.0) - 0.5;
                R = (Y + Cr * (2.0 - 2.0 * Kr));
                G = (Y - (Kb / Kg) * (2.0 - 2.0 * Kb) * Cb - (Kr / Kg) * (2.0 - 2.0 * Kr) * Cr);
                B = (Y + (2.0 - 2.0 * Kb) * Cb);
                if (R < 0) R = 0;
                if (G < 0) G = 0;
                if (B < 0) B = 0;
                if (R > 1) R = 1;
                if (G > 1) G = 1;
                if (B > 1) B = 1;
                image[i].A = R * 255.0;
                image[i].B = G * 255.0;
                image[i].C = B * 255.0;
            }
            break;
        case YCoCg:
            for (int i = 0; i < Width * Height; i++) {
                Y = image[i].A / 255.0;
                Co = (image[i].B / 255.0) - 0.5;
                Cg = (image[i].C / 255.0) - 0.5;
                R = Y + Co - Cg;
                G = Y + Cg;
                B = Y - Co - Cg;
                if (R < 0) R = 0;
                if (G < 0) G = 0;
                if (B < 0) B = 0;
                if (R > 1) R = 1;
                if (G > 1) G = 1;
                if (B > 1) B = 1;
                image[i].A = R * 255.0;
                image[i].B = G * 255.0;
                image[i].C = B * 255.0;
            }
            break;
        case CMY:
            for (int i = 0; i < Width * Height; i++) {
                C = image[i].A / 255.0;
                M = image[i].B / 255.0;
                Y = image[i].C / 255.0;
                R = 1 - C;
                G = 1 - M;
                B = 1 - Y;
                image[i].A = R * 255.0;
                image[i].B = G * 255.0;
                image[i].C = B * 255.0;
            }
            break;
    }
    currentCS = RGB;
}

void Picture::SetCS(ColorSpace current){
    currentCS = current;
}