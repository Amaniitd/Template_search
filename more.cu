#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
using namespace std;

/*
index 0  - -45, 1 - 0, 2 - 45 
*/
//range = xmin + xmax + ymin + ymax
__global__ void filter(int * data_image, int * range, int query_grey,int row, int col, int TH2, char * filtered){//filer the candidates for calculting the RMSD
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if(x>=col || y>=row){
        return;
    }
    int theta = blockIdx.z;
    float grey_sum=0;
    for (int Y= y+ range[theta]; Y< y + range[3+theta]; Y++ ){//grey avg of image
        for (int X=x+range[6+theta]; X< x + range[9+theta]; X++){
            for (int color=0; color<3; color++){
            
                if(X>=0 && Y>= 0 && X< row && Y < col){
                    grey_sum += data_image[(row-1 - X)*(col)*3 + Y*3 + color];
                }else{
                    grey_sum +=255;
                }
            }
        }
    }
    grey_sum = grey_sum/( (range[9+theta] - range[6+theta])*(range[3+theta] - range[theta])*3);

    if ( (grey_sum - query_grey)<TH2 && (grey_sum - query_grey) > -1*TH2){
        //printf("%d,%d,%d,%0.6f\n", x,y,theta,grey_sum);
        filtered[(row-1-y)*col*3 + x*3 + theta] = 1;
    }else{
        filtered[(row-1-y)*col*3 + x*3 + theta] = 0;
    }
    
}

__global__ void RMSDselect(int * data_image, int * query_image,float* cossin,
int * coordinates,const float TH1, const int row, const int col, const int query_row, const int query_col, float * top, const int maxIndex){
    
    int index = blockIdx.x*blockDim.x + threadIdx.x;
    if(index >= maxIndex){
        return;
    }
    int x = coordinates[3*index];
    int y = coordinates[3*index+1];
    int theta = coordinates[3*index+2];
    
    top[index] = -1;
    float square_diff=0.0;

    for(int qx=0; qx < query_row; qx++){
        for (int qy=0; qy<query_col; qy++){
            float rotated_x = x + qx*cossin[theta] + qy*cossin[3+theta];//rotated x
            float rotated_y = y + qy*cossin[theta] - qx*cossin[3+theta];
            int x0 = floor(rotated_x);
            int x1 = ceil(rotated_x);
            int y0 = floor(rotated_y);
            int y1 = ceil(rotated_y);
            float color_val=0.0;
            float query_val=0.0;

            for (int c=0;c<3;c++){
                color_val=255;
                int z00=255,z01=255,z10=255,z11=255;
                if(x0==x1 && y0==y1){
                    if(y0>=0 && x0>=0 && y0<col && x0<row){
                        color_val = data_image[(row-1-x0)*col*3 + y0*3+c];
                    }
                }else if(y0==y1 && y0>=0 && y0<col){
                    if(x0>=0 && x0<row){
                        z00 = data_image[(row-x0-1)*col*3 + y0*3 + c];
                    }
                    if(x1>=0 && x1<row){
                        z10 = data_image[(row-x1-1)*col*3 + y0*3 + c];
                    }
                    color_val = (x1-rotated_x)*z00+(x0-rotated_x)*z10;
                }else if (x0==x1 && x0>=0 && x0<row){
                    if(y0>=0 && y0<col){
                        z00 = data_image[(row-x0-1)*col*3 + y0*3 + c];
                    }
                    if(y1>=0 && y1<col){
                        z01 = data_image[(row-x0-1)*col*3 + y1*3 + c];
                    }
                    color_val = (y1-rotated_y)*z00+(y0-rotated_y)*z01;
                }
                else if(x0!=x1 && y0!=y1){
                    //bilinear over the 
                    if(x0<row && x0>=0){
                        if(y0<col && y0>=0)
                            z00 = data_image[(row-x0-1)*col*3 + y0*3 + c];
                        if(y1<col && y1>=0){
                            z01 = data_image[(row-x0-1)*col*3 + y1*3 + c];
                        }
                    }
                    
                    if(x1<row && x1>=0){
                        if(y0<col && y0>=0)
                            z10 = data_image[(row-x1-1)*col*3 + y0*3 + c];
                        if(y1<col && y1>=0){
                            z11 = data_image[(row-x1-1)*col*3 + y1*3 + c];
                        }
                    }
                    
                    color_val = z00*(x1 - rotated_x)*(y1 - rotated_y) + z10*(x0 - rotated_x)*(y1-rotated_y);
                    color_val+= z01*(x1-rotated_x)*(y0-rotated_y) + z11*(x0-rotated_x)*(y0-rotated_y);
                } 
                query_val = query_image[(query_row - qx-1)*query_col*3 + qy*3 + c];
                square_diff += (color_val-query_val)*(color_val-query_val);
            
            }
        }
    }
    float RMSD = sqrt(square_diff/(query_row*query_col*3));//figure out a way to find the square root
    if(x==200 && y==200){
        printf("RMSD = %.1f, TH1 = %.1f\n",RMSD, TH1);
    }
    if(RMSD < TH1){
        top[index]=RMSD;
    }

}

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
    cout<<"size of data "<<sizeof(data_image)/sizeof(data_image[0])<<'\n';
    cout<<"file read\n";
    //0 = -45, 1 = 0, 2 = 45
    int range[] = {0, 0, floor(-0.707*query_image_n),
    ceil((query_image_m + query_image_n)*0.707), query_image_m, ceil(query_image_m*0.707),
    floor(-0.707* query_image_m), 0, 0,
    ceil(query_image_n*0.707), query_image_n, ceil((query_image_m+query_image_n)*0.707)};
    
    
    
    char filtered[data_image_n*data_image_m*3];
    int *dimage, *dqimage;
    char * dfiltered;
    int * drange;
    
    cudaMalloc(&drange, 12*4);
    cudaMalloc(&dimage,  4*data_image_n*data_image_m*3);
    cudaMalloc(&dfiltered, data_image_n*data_image_m*3);//allocate an array for to mark the filtered values
    cudaMalloc(&dqimage,  4*query_image_n*query_image_m*3);

    cudaMemcpy(drange, range, 12*4, cudaMemcpyHostToDevice);
    cudaMemcpy(dimage, data_image,  data_image_m*data_image_n*3*4, cudaMemcpyHostToDevice);
    cudaMemcpy(dqimage, query_image,  query_image_m*query_image_n*3*4, cudaMemcpyHostToDevice);
    
    
    float TH2 = 10;//for difference between grey of the q and d
    //float TH1 = 20;//RMSD difference
    
    int row = data_image_m;
    int col = data_image_n;
    int query_row = query_image_m;
    int query_col = query_image_n;

    float query_grey = 0;
    for (int i=0;i<query_row*query_col*3; i++){//calculating the grey average of the query
        query_grey += query_image[i];
    }
    query_grey = query_grey/(query_image_m*query_image_n*3);
    cout<<"query grey "<<query_grey<<'\n';
    
    
    int rootx = (col+31)/32;
    int rooty = (row+31)/32;
    
    
    dim3 blocks(rootx,rooty,3);
    dim3 threads(32,32,1);//many extra threads will be there
    
    filter<<<blocks, threads>>>(dimage, drange, query_grey, row, col, TH2, dfiltered);

    cudaDeviceSynchronize();
    cudaMemcpy(filtered, dfiltered, row*col*3, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    cout<<cudaGetErrorString(err)<<'\n';

    vector<int> coordinates;
    for (int x=0;x<row; x++){
        for (int y=0; y<col;y++){
            for (int theta=0;theta<3;theta++){
                if(filtered[(row-x-1)*col*3 + y*3 + theta]==1){
                    if(x==200 && y==200){
                        cout<<"yes!!!\n";
                    }
                    //cout<<j<<','<<i<<','<<theta<<'\n';
                    coordinates.push_back(x);
                    coordinates.push_back(y);
                    coordinates.push_back(theta);
                }
            }
        }
    }
    cudaFree(dfiltered);
    
    
    int no_filtered = coordinates.size()/3;
    cout<<"no_filtered = "<<no_filtered<<'\n';
    int * dcoordinates;
    float * dcossin;
    float * dtops;
    float tops[no_filtered];
    float TH1 = 9;
    float cossin[] = {0.707, 1, 0.707, -0.707, 0, 0.707};
    cudaMalloc(&dcossin, 6*sizeof(float));
    cudaMemcpy(dcossin, cossin, 6*sizeof(float), cudaMemcpyHostToDevice);
    cudaMalloc(&dtops,no_filtered*sizeof(float));

    cudaMalloc(&dcoordinates, coordinates.size()*4);
    cudaMemcpy(dcoordinates, coordinates.data(), coordinates.size()*4, cudaMemcpyHostToDevice);
    
    int nblocks = (no_filtered+1023)/1024;
    //int nthreads = no_filtered/nblocks + no_filtered%nblocks;
    RMSDselect<<<nblocks, 1024>>>(dimage, dqimage, dcossin, dcoordinates,TH1, row, col, query_row, query_col, dtops, no_filtered);
    
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    cout<<cudaGetErrorString(err)<<'\n';
    cudaMemcpy(tops, dtops, no_filtered*sizeof(float), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    for (int i=0;i<no_filtered;i++){
        if(tops[i]!=-1){
            //if(coordinates[3*i+2] == 0){
                cout<<coordinates[3*i]<<','<<coordinates[3*i+1]<<','<<coordinates[3*i+2]<<'\n';
            //}
        }
    }
    cout<<"completed\n";

    return 0;
}