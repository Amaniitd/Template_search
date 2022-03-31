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
    cout<<"file opened\n";
    int data_image_m, data_image_n;
    int query_image_m, query_image_n;
    data_image_file >> data_image_m >> data_image_n;//m = cols (x) n = rows (y)
    query_image_file >> query_image_m >> query_image_n;
    
    cout<<data_image_m<<','<<data_image_n<<'\n';
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
    data_image_file.close();
    query_image_file.close();
    cout<<"file read\n";
    //Limits of the bounding box
    float xmin[] = {0, 0, -0.707*query_image_n};//0 = -45, 1 = 0, 2 = 45
    float xmax[] = {(query_image_m + query_image_n)*707, query_image_m, query_image_m*0.707};
    float ymin[] = {-0.707* query_image_m, 0, 0 };
    float ymax[] = {query_image_n*0.707, query_image_n, (query_image_m+query_image_n)*0.707};
    float cos[] = {0.707, 1, 0.707};
    float sin[] = {-0.707,0, 0.707};
    float grey_sum=0;
    bool RMS = false;
    int row = data_image_n;
    int col = data_image_m;
    
    float TH2 = 10.0;//for difference between grey of the q and d
    float TH1 = 0.1;//RMSD difference
    
    int topper[3];
    int bestscore = 255;

    float query_grey = 0;
    for (int i=0;i<query_image_m; i++){//calculating the grey average of the query
        for (int j=0;j<query_image_n; j++){
            for (int c=0;c<3;c++)
                query_grey += query_image[i*query_image_n*3 + j*3+c];
        }
    }
    query_grey = query_grey/(query_image_m*query_image_n*3);
    cout<<"query grey "<<query_grey<<'\n';
    int iter =0;
    for (int theta=0; theta<3; theta++){
        for (int x=0; x< data_image_m ;x++){
            for (int y=0; y< data_image_n; y++){
                grey_sum=0;
                for (int X= x+ xmin[theta]; X< x + xmax[theta]; X++ ){//grey avg of image
                    for (int Y=y+ymin[theta]; Y< y + ymax[theta]; Y++){
                        for (int color=0; color<3; color++){
                            if(X< data_image_m && Y<data_image_n){
                                grey_sum += data_image[(row-1 - Y)*(col)*3 + X*3 + color];
                            }else{
                                grey_sum +=255;
                            }
                        }
                    }
                }
                grey_sum = grey_sum/( (ymax[theta] - ymin[theta])*(xmax[theta] - xmin[theta])*3);
                
                /*if(abs(grey_sum - query_grey)<TH2){//for checking filter work
                    topper[0] = x;
                    topper[1] = y;
                    topper[2] = theta;
                    TH2 = abs(grey_sum - query_grey);
                    cout<<"grey diff = "<<TH2<<'\n';
                    cout<<"tuple "<<topper[0]<<','<<topper[1]<<','<<topper[2]<<'\n';
                }*/
                if ( abs(grey_sum - query_grey)<TH2){//filtered - calculate the RMSD
                    float square_diff=0.0;
                    for(int qx=0; qx < query_image_m; qx++){
                        for (int qy=0; qy<query_image_n; qy++){
                            float rotated_x = x + qx*cos[theta] - qy*sin[theta];//rotated x
                            float rotated_y = y + qx*sin[theta] + qy*cos[theta];
                            int x0 = floor(rotated_x);
                            int x1 = ceil(rotated_x);
                            int y0 = floor(rotated_y);
                            int y1 = ceil(rotated_y);
                            float color_val=0.0;
                            float query_val=0.0;

                            for (int c=0;c<3;c++){
                                color_val=255;
                                int z00=255,z01=255,z10=255,z11=255;
                                if(x0==x1 && y0==y1 && x0<col && y0<row){
                                    color_val = data_image[(row-1-y0)*col*3 + x0*3+c];
                                }else if(x0==x1 && x0<col){
                                    if(y0>=0 && y0<row){
                                        z00 = data_image[(row-y0-1)*col*3 + x0*3 + c];
                                    }
                                    if(y1>=0 && y1<row){
                                        z01 = data_image[(row-y1-1)*col*3 + x0*3 + c];
                                    }
                                    color_val = (y1-rotated_y)*z00+(y0-rotated_y)*z01;
                                }else if(y0==y1 && y0<row){
                                    if(x0>=0 && x0<col){
                                        z00 = data_image[(row-y0-1)*col*3 + x0*3 + c];
                                    }
                                    if(x1>=0 && x1<col){
                                        z10 = data_image[(row-y0-1)*col*3 + x1*3 + c];
                                    }
                                    color_val = (x1-rotated_y)*z00+(x0-rotated_y)*z10;
                                }
                                else if(x0!=x1 && y0!=y1){
                                    //bilinear over the 
                                    if(y0<col && y0>=0){
                                        if(x0<row && x0>=0)
                                            z00 = data_image[(row-y0-1)*col*3 + x0*3 + c];
                                        if(x1<row && x1>=0){
                                            z10 = data_image[(row-y0-1)*col*3 + x1*3 + c];
                                        }
                                    }
                                    
                                    if(y1<col && y1>=0){
                                        if(x0<row && x0>=0)
                                            z01 = data_image[(row-y1-1)*col*3 + x0*3 + c];
                                        if(x1<row && x1>=0){
                                            z11 = data_image[(row-y1-1)*col*3 + x1*3 + c];
                                        }
                                    }
                                    
                                    color_val = z00*(x1 - rotated_x)*(y1 - rotated_y) + z10*(x0 - rotated_x)*(y1-rotated_y);
                                    color_val+= z01*(x1-rotated_x)*(y0-rotated_y) + z11*(x0-rotated_x)*(y0-rotated_y);
                                } 
                                query_val = query_image[(query_image_n - qy-1)*query_image_m*3 + qx*3 + c];
                                square_diff += pow((color_val-query_val), 2);
                            
                            }
                        }
                    }
                    float RMSD = sqrt(square_diff/(query_image_m*query_image_n*3));
                    if (RMSD < TH1 && RMSD < bestscore){
                        topper[0] = x;
                        topper[1] = y;
                        topper[2] = theta;
                        bestscore = RMSD;
                        cout<<"grey sum "<<grey_sum<<", RMSD "<<RMSD<<'\n';
                        cout<<"tuple "<<topper[0]<<','<<topper[1]<<','<<topper[2]<<'\n';
                    }

                }
                    
            }
        }
    }
    cout<<"FINAL RESULT = "<<topper[0]<<','<<topper[1]<<','<<topper[2];
    /*cout << "data_image_m: " << data_image_m << endl;
    cout << "data_image_n: " << data_image_n << endl;
    cout << "query_image_m: " << query_image_m << endl;
    cout << "query_image_n: " << query_image_n << endl;*/

    return 0;
}