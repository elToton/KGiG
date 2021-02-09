#include "Picture.h"
#include <iostream>

using namespace std;

//enum {inv, mirh, mirv, rotr, rotl} num;

int main(int argc, char** argv){
    Picture picture;
    //..................................................................................................................
    /*int num = argv[3][0] - '0';
    switch (num){
        case inv:
            picture.Inverse();
            break;
        case mirh:
            picture.Mirror("H");
            break;
        case mirv:
            picture.Mirror("V");
            break;
        case rotr:
            picture.Rotation("R");
            break;
        case rotl:
            picture.Rotation("L");
            break;
        default:
            cout << "you wrote wrong number"<< endl;
    }*/
    //..................................................................................................................
    /*byte color = stoi(argv[3]);
    float thickness = stof(argv[4]);
    float x0 = stof(argv[5]);
    float y0 = stof(argv[6]);
    float x1 = stof(argv[7]);
    float y1 = stof(argv[8]);
    float gamma = 2.2;
    if (argc == 10) {
        gamma = stof(argv[9]);
    }
    picture.DrawLine(x0, y0, x1, y1, color, thickness, gamma);*/
    //..................................................................................................................
    /*bool gradient = atoi(argv[3]);
    int dithering = atoi(argv[4]);
    int bits = atoi(argv[5]);
    picture.correction = atof(argv[6]);
    if (gradient) picture.Grad();
    switch(dithering){
        case 0:
            picture.WithOutDithering(bits);
            break;
        case 1:
            picture.Ordered(bits);
            break;
        case 2:
            picture.Random(bits);
            break;
        case 3:
            picture.FloydSteinberg(bits);
            break;
        case 4:
            picture.JarvisJudiceNinke(bits);
            break;
        case 5:
            picture.Sierra(bits);
            break;
        case 6:
            picture.Atkinson(bits);
            break;
        case 7:
            picture.Halftone(bits);
            break;
    }*/
    //..................................................................................................................
    ColorSpace convertion;
    bool ThreeInput, ThreeOutput;
    string InFile, OutFile, buffer;
    for (int position = 1; position < argc; position++){
        buffer = argv[position];
        if (buffer == "-f"){
            position++;
            buffer = argv[position];
            if (buffer == "RGB"){
                picture.SetCS(RGB);
                continue;
            }
            if (buffer =="HSL"){
                picture.SetCS(HSL);
                continue;
            }
            if (buffer == "HSV"){
                picture.SetCS(HSV);
                continue;
            }
            if (buffer == "YCbCr.601"){
                picture.SetCS(YCbCr_601);
                continue;
            }
            if (buffer == "YCbCr.709"){
                picture.SetCS(YCbCr_709);
                continue;
            }
            if (buffer == "YCoCg"){
                picture.SetCS(YCoCg);
                continue;
            }
            if (buffer == "CMY"){
                picture.SetCS(CMY);
                continue;
            }
            return 1;
        }
        if (buffer =="-t"){
            position++;
            buffer = argv[position];
            if (buffer == "RGB"){
                convertion = RGB;
                continue;
            }
            if (buffer == "HSL"){
                convertion = HSL;
                continue;
            }
            if (buffer == "HSV"){
                convertion = HSV;
                continue;
            }
            if (buffer == "YCbCr.601"){
                convertion = YCbCr_601;
                continue;
            }
            if (buffer == "YCbCr.709"){
                convertion = YCbCr_709;
                continue;
            }
            if (buffer == "YCoCg"){
                convertion = YCoCg;
                continue;
            }
            if (buffer == "CMY"){
                convertion = CMY;
                continue;
            }
            return 1;
        }
        if (buffer =="-i"){
            position++;
            buffer = argv[position];
            if (buffer =="1") ThreeInput = false;
            else ThreeInput = true;
            position++;
            buffer = argv[position];
            InFile = buffer;
            continue;
        }
        if (buffer =="-o"){
            position++;
            buffer = argv[position];
            if (buffer =="1") ThreeOutput = false;
            else ThreeOutput = true;
            position++;
            buffer = argv[position];
            OutFile = buffer;
            continue;
        }
    }
    if (ThreeInput) picture.ReadThreeFiles(InFile);
    else picture.ReadOneFile(InFile);
    picture.ConvertAny(convertion);
    if (ThreeOutput) picture.WriteToThreeFiles(OutFile);
    else picture.WriteToOneFile(OutFile);
    return 0;
}