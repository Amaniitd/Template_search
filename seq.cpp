#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main(int argc, char **argv)
{
    if(argc != 3){
        cout<<"Usage: ./a.out dataImage queryImage\n";
        return 0; 
    }
    string data_image_path = argv[1];
    string query_image_path = argv[2];

    // reading a matrix from a file with first line m*n
    ifstream data_image_file(data_image_path);
    ifstream query_image_file(query_image_path);
    int data_image_m, data_image_n;
    int query_image_m, query_image_n;
    data_image_file >> data_image_m >> data_image_n;//m = cols (x) n = rows (y)
    query_image_file >> query_image_m >> query_image_n;

    int data_image[data_image_m*data_image_n*3];
    int query_image[query_image_m*query_image_n*3];

    for (int i = 0; i < data_image_m*data_image_n*3; i++)
    {
            data_image_file >> data_image[i];
    }

    for (int i = 0; i < query_image_m*query_image_n*3; i++)
    {
            query_image_file >> query_image[i];
    }

    //Limits of the bounding box
    float xmin[] = {0, 0, -0.707*query_image_n};//0 = -45, 1 = 0, 2 = 45
    float xmax[] = {(query_image_m + query_image_n)*707, query_image_m, query_image_m*0.707};
    float ymin[] = {-0.707* query_image_m, 0, 0 };
    float ymax[] = {query_image_n*0.707, query_image_n, (query_image_m+query_image_n)*0.707};
    
    double grey_sum=0;
    bool RMS = false;
    int row = data_image_n;
    int col = data_image_m;  
    
    double TH2 = 0.5;//for difference between grey of the q and d
    double TH1 = 0.1;//RMSD difference
    
    int topper[3];
    int bestscore = 255;

    double query_grey = 0;
    for (int i=0;i<query_image_m; i++){//calculating the grey average of the query
        for (int j=0;j<query_image_n; j++){
            for (int c=0;c<3;c++)
                query_grey += query_image[i*query_image_n*3 + j*3+c];
        }
    }
    query_grey = query_grey/(query_image_m*query_image_n*3);

    for (int theta=0; theta<3; theta++){
        for (int x=0; x< data_image_m ;x++){
            for (int y=0; y< data_image_n; y++){
                grey_sum=0;

                for (int X= x+ xmin[theta]; X< x + xmax[theta]; X++ ){//grey avg of image
                    for (int Y=y+ymin[theta]; Y< y + ymax[theta]; Y++){
                        for (int color=0; color<3; color)
                            grey_sum += data_image[(col-1 - Y)*(row)*3 + X*3 + color];
                    }
                }
                grey_sum = grey_sum/( (ymax[theta] - ymin[theta])*(xmax[theta] - xmin[theta])*3);


                if ( abs(grey_sum - query_grey)<TH2){//prelims passed //calculate the RMSD
                    double square_diff=0.0;
                    for(int qx=0; qx < query_image_m; qx++){
                        for (int qy=0; qy<query_image_n; qy++){
                            float rotated_x = 0.0;//rotated x
                            float rotated_y = 0.0;
                            int x0 = floor(rotated_x);
                            int x1 = ceil(rotated_x);
                            int y0 = floor(rotated_y);
                            int y1 = ceil(rotated_y);
                            float color_val=0.0;
                            float query_val=0.0;
                            for (int c=0;c<3;c++){
                                int z00=255,z01=255,z10=255,z11=255;
                                if(y0<col && y0>=0){
                                    if(x0<row && x0>=0)
                                        z00 = data_image[(col-y0-1)*row*3 + x0*3 + c];
                                    if(x1<row && x1>=0){
                                        z10 = data_image[(col-y0-1)*row*3 + x1*3 + c];
                                    }
                                }
                                
                                if(y1<col && y1>=0){
                                    if(x0<row && x0>=0)
                                        z01 = data_image[(col-y1-1)*row*3 + x0*3 + c];
                                    if(x1<row && x1>=0){
                                        z11 = data_image[(col-y1-1)*row*3 + x1*3 + c];
                                    }
                                }
                                
                                color_val = z00*(x1 - rotated_x)*(y1 - rotated_y) + z10*(x0 - rotated_x)*(1-rotated_y);
                                color_val+= z01*(x1-rotated_x)*(y0-rotated_y) + z11*(x0-rotated_x)*(y0-rotated_y); //outputs the 
                                query_val = query_image[(query_image_m - qy-1)*query_image_n*3 + qx*3 + c];
                                square_diff += pow(color_val-query_val, 2);
                            
                            }
                        }
                    }
                    double RMSD = sqrt(square_diff/(query_image_m*query_image_n*3));
                    if (RMSD < TH1 && RMSD < bestscore){
                        topper[0] = x;
                        topper[1] = y;
                        topper[2] = theta;
                        bestscore = RMSD;
                    }

                }
                    
            }
        }
    }
    cout << "data_image_m: " << data_image_m << endl;
    cout << "data_image_n: " << data_image_n << endl;
    cout << "query_image_m: " << query_image_m << endl;
    cout << "query_image_n: " << query_image_n << endl;

    return 0;
}