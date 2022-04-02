#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
using namespace std;

/*


*/
//range = xmin + xmax + ymin + ymax
__global__ void filter(int * data_image, int * query_image, float * range, int query_grey,int row, int col, int TH2, char * filtered){//filer the candidates for calculting the RMSD
    int x = blockIdx.y * blockDim.y + threadIdx.y;
    int y = blockIdx.x * blockDim.x + threadIdx.x;
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
            printf("%d,%d,%d,%0.6f\n", x,y,theta,grey_sum);
            filtered[(row-1-y)*col*3 + x*3 + theta] = 1;
        }else{
            filtered[(row-1-y)*col*3 + x*3 + theta] = 0;
        }

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
    
    //Limits of the bounding box
    //RMSD exclusive data
    //bfloat16 * drange, *dcos, *dsin;//, *dRMSD;

    //cudaMalloc(&drange, 12*2);
    /*bfloat16 cos[] = {0.707, 1, 0.707};
    bfloat16 sin[] = {-0.707,0, 0.707};
    cudaMalloc(&dcos, 3*2);
    cudaMalloc(&dsin, 3*2);
    cudaMalloc(&dRMSD, 2);
    cudaMemcpy(dcos, range, 12*2, cudaMemcpyHostToDevice);
    cudaMemcpy(dsin, range, 12*2, cudaMemcpyHostToDevice);
    int topper[3];
    int bestscore = 255;*/


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
    
    
    int four = ceil(sqrt(sqrt(row*col)));
    
    dim3 blocks(four,four,1);
    dim3 threads(four,four,1);//many extra threads will be there
    
    filter<<<blocks, threads>>>(dimage, dqimage, drange, query_grey, row, col, TH2, dfiltered);

    cudaDeviceSynchronize();
    cudaMemcpy(filtered, dfiltered, row*col*3, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    for (int i=0;i<row; i++){
        for (int j=0; j<col;j++){
            for (int theta=0;theta<3;theta++){
                if(filtered[(row-i-1)*col*3 + j*3 + theta]==1){
                    cout<<j<<','<<i<<','<<theta<<'\n';
                }else if(filtered[(row-i-1)*col*3 + j*3 + theta] == 2){
                    cout<<"not computed "<<j<<','<<i<<','<<theta<<'\n';
                }
            }
        }
    }
    cout<<cudaGetErrorString(err)<<'\n';

    return 0;
}
