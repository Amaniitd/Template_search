#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdio.h>
using namespace std;

/*


*/
//range = xmin + xmax + ymin + ymax
__global__ void filter(int * data_image, float * range, int query_grey,int row, int col, int TH2, char * filtered){//filer the candidates for calculting the RMSD
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if(x>=col || y>=row){
        return;
    }
    for (int theta=0; theta<3; theta++){
        float grey_sum=0;
        for (int X= x+ range[theta]; X< x + range[3+theta]; X++ ){//grey avg of image
            for (int Y=y+range[6+theta]; Y< y + range[9+theta]; Y++){
                for (int color=0; color<3; color++){
                
                    if(X>=0 && Y>= 0 && X< col && Y < row){
                        grey_sum += data_image[(row-1 - Y)*(col)*3 + X*3 + color];
                    }else{
                        grey_sum +=255;
                    }
                }
            }
        }
        grey_sum = grey_sum/( (range[9+theta] - range[6+theta])*(range[3+theta] - range[theta])*3);

        if ( (grey_sum - query_grey)<TH2 && (grey_sum - query_grey) > -1*TH2){
            filtered[(row-1-y)*col*3 + x*3 + theta] = 1;
        }else{
            filtered[(row-1-y)*col*3 + x*3 + theta] = 0;
        }

    }
}

__global__ void RMSDselect(int * dataImage, int * queryImage,float* cossin,
int * coordinates,int TH1, int row, int col, int n, int * top){
    
    int index = blockIdx.y*blockDim.y + blockIdx.x*blockDim.x + threadIdx.x + threadIdx.y;

    int x = coordinates[3*index];
    int y = coordinates[3*index+1];
    int theta = coordinates[3*index+2]

    float square_diff=0.0;
    for(int qx=0; qx < query_image_m; qx++){
        for (int qy=0; qy<query_image_n; qy++){
            float rotated_x = x + qx*cossin[theta] - qy*cossin[3+theta];//rotated x
            float rotated_y = y + qx*cossin[3+theta] + qy*cossin[theta];
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
                    if(x0>=0 && y0>=0 && x0<col && y0<row){
                        color_val = data_image[(row-1-y0)*col*3 + x0*3+c];
                    }
                }else if(x0==x1 && x0>=0 && x0<col){
                    if(y0>=0 && y0<row){
                        z00 = data_image[(row-y0-1)*col*3 + x0*3 + c];
                    }
                    if(y1>=0 && y1<row){
                        z01 = data_image[(row-y1-1)*col*3 + x0*3 + c];
                    }
                    color_val = (y1-rotated_y)*z00+(y0-rotated_y)*z01;
                }else if (y0==y1 && y0>=0 && y0<row){
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
                    if(y0<row && y0>=0){
                        if(x0<col && x0>=0)
                            z00 = data_image[(row-y0-1)*col*3 + x0*3 + c];
                        if(x1<col && x1>=0){
                            z10 = data_image[(row-y0-1)*col*3 + x1*3 + c];
                        }
                    }
                    
                    if(y1<row && y1>=0){
                        if(x0<col && x0>=0)
                            z01 = data_image[(row-y1-1)*col*3 + x0*3 + c];
                        if(x1<col && x1>=0){
                            z11 = data_image[(row-y1-1)*col*3 + x1*3 + c];
                        }
                    }
                    
                    color_val = z00*(x1 - rotated_x)*(y1 - rotated_y) + z10*(x0 - rotated_x)*(y1-rotated_y);
                    color_val+= z01*(x1-rotated_x)*(y0-rotated_y) + z11*(x0-rotated_x)*(y0-rotated_y);
                } 
                query_val = query_image[(query_image_n - qy-1)*query_image_m*3 + qx*3 + c];
                square_diff += (color_val-query_val)*(color_val-query_val);
            
            }
        }
    }
    float RMSD = sqrt(square_diff/(query_image_m*query_image_n*3));//figure out a way to find the square root
    if (RMSD < TH1 && RMSD < bestscore){
        top[index] = x;
        top[index+1] = y;
        top[index+2] = theta;
        //cout<<"grey sum "<<grey_sum<<", RMSD "<<RMSD<<'\n';
        //cout<<"tuple "<<topper[0]<<','<<topper[1]<<','<<topper[2]<<'\n';
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

    float range[] = {0, 0, -0.707*query_image_n,(query_image_m + query_image_n)*707, query_image_m, query_image_m*0.707,
    -0.707* query_image_m, 0, 0 , query_image_n*0.707, query_image_n, (query_image_m+query_image_n)*0.707};//0 = -45, 1 = 0, 2 = 45
    
    
    char filtered[data_image_n*data_image_m*3];
    for(int i=0;i<data_image_n*data_image_m*3; i++){
        filtered[i]=2;
    }
    int *dimage, *dqimage;
    char *dfiltered;
    float * drange;
    
    cudaMalloc(&drange, 12*sizeof(float));
    cudaMalloc(&dimage,  4*data_image_n*data_image_m*3);
    cudaMalloc(&dfiltered, data_image_n*data_image_m*3);//allocate an array for to mark the filtered values
    cudaMalloc(&dqimage,  4*query_image_n*query_image_m*3);

    cudaMemcpy(dfiltered, filtered, data_image_m*data_image_n*3, cudaMemcpyHostToDevice);
    cudaMemcpy(drange, range, 12*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dimage, data_image,  data_image_m*data_image_n*3*4, cudaMemcpyHostToDevice);
    cudaMemcpy(dqimage, query_image,  query_image_m*query_image_n*3*4, cudaMemcpyHostToDevice);
    
    


    int row = data_image_n;
    int col = data_image_m;
    
    float TH2 = (query_image_m*query_image_n)/5;//for difference between grey of the q and d
    //float TH1 = 20;//RMSD difference
    
    

    float query_grey = 0;
    for (int i=0;i<query_image_m; i++){//calculating the grey average of the query
        for (int j=0;j<query_image_n; j++){
            for (int c=0;c<3;c++)
                query_grey += query_image[i*query_image_n*3 + j*3+c];
        }
    }
    query_grey = query_grey/(query_image_m*query_image_n*3);
    cout<<"query grey "<<query_grey<<'\n';
    
    
    int rootx = ceil(sqrt(col));
    int rooty = ceil(sqrt(row));
    
    
    dim3 blocks(rootx,rooty,1);
    dim3 threads(rootx,rooty,1);//many extra threads will be there
    
    filter<<<blocks, threads>>>(dimage, dqimage, drange, query_grey, row, col, TH2, dfiltered);

    cudaDeviceSynchronize();
    cudaMemcpy(filtered, dfiltered, row*col*3, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    vector<int> coordinates;
    cudaError_t err = cudaGetLastError();
    for (int i=0;i<row; i++){
        for (int j=0; j<col;j++){
            for (int theta=0;theta<3;theta++){
                if(filtered[(row-i-1)*col*3 + j*3 + theta]==1){
                    cout<<j<<','<<i<<','<<theta<<'\n';
                    coordinates.push_back(j);
                    coordinates.push_back(i);
                    coordinates.push_back(theta);
                }
            }
        }
    }
    int n = 2;
    float *dcossin;//, *dRMSD;
    int * dtops;
    int topn[3*n];
    int TH1 = 10;
    float cossin[] = {0.707, 1, 0.707, -0.707, 0, 0.707};
    cudaMalloc(&dcossin, 6*sizeof(float));
    cudaMalloc(&dtops,4*coordinates.size());
    cudaMemcpy(&dcossin, cossin, 6*sizeof(float), cudaMemcpyHostToDevice);

    cudaMalloc(&dcoordinates, coordinates.size()*4);
    cudaMemcpy(dcoordinates, coordinates.data(), coordinates.size()*4, cudaMemcpyHostToDevice);
    
    int four = ceil(sqrt(sqrt(coordinates.size()/3)));
    dim3 lessblocks(four,four,1);
    dim3 lessthreads(four,four,1);
    RMSDselect<<<lessblocks, lessthreads>>>(dimage, dqimage, dcossin, row, col, n, dtops);

    cout<<cudaGetErrorString(err)<<'\n';

    return 0;
}
