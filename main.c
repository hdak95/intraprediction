#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define WIDTH       512
#define HEIGHT      512
#define BLOCK_SIZE  4
#define NMODE       6

typedef unsigned char BYTE;

void MemFree_2D         (BYTE** arr, int height);           //2D memory free
void MemFree_2D_int     (int** arr,int height);
void FileRead           (char* filename, BYTE** img_in, int width, int height);
void FileWrite          (char* filename, BYTE** img_out, int width, int height);
void Encode             (BYTE** img_ori, BYTE** img_pred, int** img_resi, BYTE** img_recon);

BYTE** MemAlloc_2D      (int width, int height);
int**  MemAlloc_2D_int  (int width, int height);

int intra_prediction    (BYTE* ori, BYTE* ref, BYTE* planar_ref,BYTE (pred)[NMODE][BLOCK_SIZE*BLOCK_SIZE], int (resi)[NMODE][BLOCK_SIZE*BLOCK_SIZE], BYTE (recon)[NMODE][BLOCK_SIZE*BLOCK_SIZE]);
int intra_dc            (BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int intra_hor           (BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int intra_ver           (BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int intra_DL            (BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int intra_DR            (BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);
int intra_Planar            (BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon);

float GetPSNR           (BYTE** img_ori, BYTE** img_dist, int width, int height);


int main()
{
    BYTE **img_ori, **img_pred, **img_recon, ** img_in_R, **img_in_G, **img_in_B, **img_pred_R, **img_pred_G, **img_pred_B;
    int **img_resi_R, **img_resi_G, **img_resi_B;

    int i,j;
    
    img_ori     = MemAlloc_2D(WIDTH, HEIGHT*3);
    img_pred    = MemAlloc_2D(WIDTH, HEIGHT*3);
    img_recon   = MemAlloc_2D(WIDTH, HEIGHT*3);

    img_in_R    = MemAlloc_2D(WIDTH, HEIGHT);
    img_in_B    = MemAlloc_2D(WIDTH, HEIGHT);
    img_in_G    = MemAlloc_2D(WIDTH, HEIGHT);

    img_pred_R  = MemAlloc_2D(WIDTH, HEIGHT);
    img_pred_G  = MemAlloc_2D(WIDTH, HEIGHT);
    img_pred_B  = MemAlloc_2D(WIDTH, HEIGHT);

    img_resi_R  = MemAlloc_2D_int(WIDTH, HEIGHT);
    img_resi_G  = MemAlloc_2D_int(WIDTH, HEIGHT);
    img_resi_B  = MemAlloc_2D_int(WIDTH, HEIGHT);

    FileRead("./Lena(512x512).rgb",img_ori,WIDTH,HEIGHT*3);

    for(i = 0 ; i < HEIGHT ; i++){
        memcpy(img_in_R[i],img_ori[i],sizeof(BYTE) * WIDTH);
    }
    for(i = 0 ; i < HEIGHT ; i++){
        memcpy(img_in_G[i],img_ori[i + HEIGHT],sizeof(BYTE) * WIDTH);
    }
    for(i = 0 ; i < HEIGHT ; i++){
        memcpy(img_in_B[i],img_ori[i + (HEIGHT<<1) ],sizeof(BYTE) * WIDTH);
    }

    //////////////////////////////////////////////////////
    /*          Intra Prediction Processing             */
    //////////////////////////////////////////////////////

    Encode(img_in_R, img_pred_R, img_resi_R, &img_recon[0]);

    Encode(img_in_G, img_pred_G, img_resi_G, &img_recon[HEIGHT]);

    Encode(img_in_B, img_pred_B, img_resi_B, &img_recon[HEIGHT*2]);

    //////////////////////////////////////////////////////

    // merging result image
    for(i=0;i<HEIGHT;i++){
        memcpy(img_pred[i],img_pred_R[i],sizeof(BYTE) * WIDTH);
    }
    for(i=0;i<HEIGHT;i++){
        memcpy(img_pred[i+HEIGHT],img_pred_G[i],sizeof(BYTE) * WIDTH);
    }
    for(i=0;i<HEIGHT;i++){
        memcpy(img_pred[i+HEIGHT*2],img_pred_B[i],sizeof(BYTE) * WIDTH);
    }

    // get psnr

    printf("PREDICTION VS ORIGINAL PSNR : %.2f\n", GetPSNR(img_ori,img_pred,WIDTH,HEIGHT*3));

    printf("RECON VS ORIGINAL PSNR : %.2f\n", GetPSNR(img_ori,img_recon,WIDTH,HEIGHT*3));

    // get result

    FileWrite("[Intra]Lena(512x512).rgb",img_pred,WIDTH,HEIGHT*3);
    
    FileWrite("[Recon]Lena(512x512).rgb",img_recon,WIDTH,HEIGHT*3);

    // memory free
    MemFree_2D(img_in_R,HEIGHT);
    MemFree_2D(img_in_G,HEIGHT);
    MemFree_2D(img_in_B,HEIGHT);
    MemFree_2D(img_pred_R,HEIGHT);
    MemFree_2D(img_pred_G,HEIGHT);
    MemFree_2D(img_pred_B,HEIGHT);
    MemFree_2D(img_ori , HEIGHT*3);
    MemFree_2D(img_pred, HEIGHT*3);
    MemFree_2D_int(img_resi_R,HEIGHT);
    MemFree_2D_int(img_resi_G,HEIGHT);
    MemFree_2D_int(img_resi_B,HEIGHT);

    return 0;
}
BYTE** MemAlloc_2D(int width, int height){
    BYTE** arr;
    int i;

    arr = (BYTE**)malloc(sizeof(BYTE*) * height);
    for(i=0; i<height; i++)
        arr[i]= (BYTE*)malloc(sizeof(BYTE) * width);

    return arr;
}

int** MemAlloc_2D_int(int width, int height){
    int** arr;
    int i;

    arr = (int**)malloc(sizeof(int*) * height);
    for(i=0; i<height; i++)
        arr[i] = (int*)malloc(sizeof(int) * width);

    return arr;
}

void MemFree_2D(BYTE** arr, int height){
    int i;
    for(i=0; i<height; i++){
        free(arr[i]);
    }
    free(arr);
}

void MemFree_2D_int(int** arr, int height){
    int i;
    for(i=0; i<height; i++){
        free(arr[i]);
    }
    free(arr);
}

void FileRead(char* filename, BYTE** img_in, int width, int height){
    FILE* fp_in;
    int i;
    fp_in = fopen(filename, "rb");
    for(i = 0 ; i < height ; i++)
        fread(img_in[i], sizeof(BYTE), width, fp_in);
    fclose(fp_in);
}

void FileWrite(char* filename, BYTE** img_out, int width, int height){
    FILE* fp_out;
    int i;

    fp_out = fopen(filename, "wb");
    for(i = 0 ; i < height ; i++)
        fwrite(img_out[i],sizeof(BYTE), width, fp_out);
    fclose(fp_out);
}

float GetPSNR(BYTE** img_ori, BYTE** img_dist, int width, int height){      //  PSNR calculation
    float mse=0;
    int i,j;

    for(i=0 ; i < height ; i++){                        // MSE calculation
        for(j = 0 ; j< width ; j++){
            mse += ((img_ori[i][j] - img_dist[i][j]) * (img_ori[i][j] - img_dist[i][j])) / (float)(width*height);
        }
    }
    return 10*(float)log10((255*255)/mse);      //  PSNR
}

void Encode(BYTE** img_ori, BYTE** img_pred, int** img_resi, BYTE** img_recon){
    int i,j,m,n;
    int best_mode;
    int min_SAD,temp_SAD;

    static BYTE ori   [BLOCK_SIZE*BLOCK_SIZE];
    static BYTE ref   [BLOCK_SIZE*3+1];
    static BYTE planar_ref [BLOCK_SIZE*2+2];
    static BYTE pred  [NMODE][BLOCK_SIZE*BLOCK_SIZE];
    static BYTE recon [NMODE][BLOCK_SIZE*BLOCK_SIZE];
    static int resi   [NMODE][BLOCK_SIZE*BLOCK_SIZE];

    int counter[NMODE]={0,};


    BYTE** img_padding = MemAlloc_2D(WIDTH + 1, HEIGHT +1);

    for(i=0; i < HEIGHT; i++){
        for(j=0;j<WIDTH;j++){
            img_padding[i+1][j+1]=img_ori[i][j];
        }
    }

    for(i=0; i<HEIGHT; i++)
        img_padding[i+1][0] = 128;
    
    for(i=0; i<WIDTH + 1; i++)
        img_padding[0][i] = 128;

    //intra prediction loop
    for(i=0; i < HEIGHT; i += BLOCK_SIZE){
        for(j=0; j < WIDTH ; j += BLOCK_SIZE){
            // get original block
            for(m=0;m<BLOCK_SIZE;m++){
                for(n=0;n<BLOCK_SIZE;n++){
                    ori[m * BLOCK_SIZE + n] = img_ori[i + m][j + n];
                }
            }

            //get reference samples
            if(j != WIDTH - BLOCK_SIZE){
                for(m=0;m<2*BLOCK_SIZE+1;m++)
                    ref[m] = img_padding[i][j+m];
            }
            else{
                for(m=0;m<BLOCK_SIZE+1;m++)
                    ref[m] = img_padding[i][j+m];
                for(m=BLOCK_SIZE+1;m<2*BLOCK_SIZE+1;m++)
                    ref[m]=128;
            }
            if(j != 0){
                for(m=0;m<BLOCK_SIZE;m++)
                    ref[m + (2*BLOCK_SIZE) + 1] = img_padding[(i+1)+m][j];
            }
            else{
                for(m=0;m<BLOCK_SIZE;m++)
                    ref[m + (2*BLOCK_SIZE) + 1] = 128;
            }

            //get planar reference samples
            if(j != WIDTH - BLOCK_SIZE){
                for(m=0;m<BLOCK_SIZE+1;m++)
                    planar_ref[m]=img_padding[i][j+m+1];
            }
            else{
                for(m=0;m<BLOCK_SIZE;m++)
                    planar_ref[m]=img_padding[i][j+m+1];
                planar_ref[BLOCK_SIZE]=img_padding[i][j+m];
            }
            if(i!=HEIGHT-BLOCK_SIZE){
                for(m=0;m<BLOCK_SIZE+1;m++)
                    planar_ref[m+BLOCK_SIZE+1]=img_padding[i+m+1][j];
            }
            else{
                for(m=0;m<BLOCK_SIZE;m++)
                    planar_ref[m+BLOCK_SIZE+1]=img_padding[i+m+1][j];
                planar_ref[2*BLOCK_SIZE+1]=img_padding[i+m][j];
            }
            
            // search for best mode
            best_mode = intra_prediction(ori, ref, planar_ref,
                                         pred, resi, recon);
            counter[best_mode]++;
           
            // generate reconstructed image
            for(m=0;m<BLOCK_SIZE;m++){
                for(n=0;n<BLOCK_SIZE;n++){
                    img_pred [i+m][j+n]=pred [best_mode][m*BLOCK_SIZE + n];
                    img_resi [i+m][j+n]=resi [best_mode][m*BLOCK_SIZE + n];
                    img_recon[i+m][j+n]=recon[best_mode][m*BLOCK_SIZE + n];
                }
            }
        }
    }
    for(i=0;i<NMODE;i++){
        printf("%d mode : %d\n",i,counter[i]);
    }
    MemFree_2D(img_padding,HEIGHT + 1);
}

int intra_prediction(BYTE* ori, BYTE* ref, BYTE* planar_ref,
                     BYTE (pred)[NMODE][BLOCK_SIZE*BLOCK_SIZE]  ,
                     int (resi)[NMODE][BLOCK_SIZE*BLOCK_SIZE]   ,
                     BYTE (recon)[NMODE][BLOCK_SIZE*BLOCK_SIZE] ){
    static int SAD[NMODE];
    int min_SAD,best_mode,i,j,n,m;

    SAD[0] = intra_ver(ori, ref, pred[0], resi[0], recon[0]);
    
    SAD[1] = intra_hor(ori, ref, pred[1], resi[1], recon[1]);

    SAD[2] = intra_dc (ori, ref, pred[2], resi[2], recon[2]);

    SAD[3] = intra_DL (ori, ref, pred[3], resi[3], recon[3]);

    SAD[4] = intra_DR (ori, ref, pred[4], resi[4], recon[4]);

    SAD[5] = intra_Planar (ori, planar_ref, pred[5], resi[5], recon[5]);

    best_mode = 0;
    min_SAD = SAD[0];

    for(n = 1 ; n < NMODE ; n++){
        if(min_SAD > SAD[n]){
            min_SAD = SAD[n];
            best_mode = n;
        }
    }
    return best_mode;
}

int intra_dc(BYTE* ori, BYTE* ref,
             BYTE* pred, int* resi, BYTE* recon){
    int i,j;
    int dc_Val=BLOCK_SIZE;
    int SAD=0;

    for(i = 0 ; i < BLOCK_SIZE ; i++)
        dc_Val += ref[i + 1] + ref[BLOCK_SIZE * 2 + i + 1];
    dc_Val = (int)(dc_Val+(log2(BLOCK_SIZE)+1)/2) >> ((int)log2(BLOCK_SIZE)+1);

    for(i = 0 ; i < BLOCK_SIZE ; i++){
        for(j = 0 ; j < BLOCK_SIZE ; j++){
            if(i==0&&j==0){
                pred[i*BLOCK_SIZE+j]=((ref[1]+2*dc_Val+ref[2*BLOCK_SIZE+1]+2)+1)>>2;
            }else if(i==0){
                pred[i*BLOCK_SIZE+j]=((ref[j+1]+3*dc_Val+2)+1)>>2;
            }else if(j==0){
                pred[i*BLOCK_SIZE+j]=((ref[2*BLOCK_SIZE+i+1]+3*dc_Val+2)+1)>>2;
            }else{
                pred[i*BLOCK_SIZE+j]=dc_Val;
            }
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            resi[i*BLOCK_SIZE+j]=ori[i*BLOCK_SIZE+j]-pred[i*BLOCK_SIZE+j];
            SAD+=abs(resi[i*BLOCK_SIZE+j]);
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            recon[i*BLOCK_SIZE+j]=resi[i*BLOCK_SIZE+j]+pred[i*BLOCK_SIZE+j];
        }
    }

    return SAD;
}

int intra_hor(BYTE* ori, BYTE* ref,
              BYTE* pred, int* resi, BYTE* recon){
    int i,j;
    int SAD=0;

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            pred[i*BLOCK_SIZE+j]=ref[2*BLOCK_SIZE+i+1];
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            resi[i*BLOCK_SIZE+j]=ori[i*BLOCK_SIZE+j]-pred[i*BLOCK_SIZE+j];
            SAD+=abs(resi[i*BLOCK_SIZE+j]);
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            recon[i*BLOCK_SIZE+j]=resi[i*BLOCK_SIZE+j]+pred[i*BLOCK_SIZE+j];
        }
    }

    return SAD;
}

int intra_ver(BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon){
    int i,j;
    int SAD=0;

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            pred[i*BLOCK_SIZE+j]=ref[j+1];
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            resi[i*BLOCK_SIZE+j]=ori[i*BLOCK_SIZE+j]-pred[i*BLOCK_SIZE+j];
            SAD+=abs(resi[i*BLOCK_SIZE+j]);
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            recon[i*BLOCK_SIZE+j]=resi[i*BLOCK_SIZE+j]+pred[i*BLOCK_SIZE+j];
        }
    }

    return SAD;
}

int intra_DL(BYTE* ori, BYTE* ref, BYTE* pred,int* resi, BYTE* recon){
    int i,j;
    int SAD=0;

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            pred[i*BLOCK_SIZE+j]=ref[i+j+2];
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            resi[i*BLOCK_SIZE+j]=ori[i*BLOCK_SIZE+j]-pred[i*BLOCK_SIZE+j];
            SAD+=abs(resi[i*BLOCK_SIZE+j]);
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            recon[i*BLOCK_SIZE+j]=resi[i*BLOCK_SIZE+j]+pred[i*BLOCK_SIZE+j];
        }
    }

    return SAD;
}

int intra_DR(BYTE* ori, BYTE* ref, BYTE* pred,int* resi, BYTE* recon){
    int i,j;
    int SAD=0;

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            if(i<=j){
                pred[i*BLOCK_SIZE+j]=ref[j-i];
            }
            else{
                pred[i*BLOCK_SIZE+j]=ref[2*BLOCK_SIZE+i-j];
            }
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            resi[i*BLOCK_SIZE+j]=ori[i*BLOCK_SIZE+j]-pred[i*BLOCK_SIZE+j];
            SAD+=abs(resi[i*BLOCK_SIZE+j]);
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            recon[i*BLOCK_SIZE+j]=resi[i*BLOCK_SIZE+j]+pred[i*BLOCK_SIZE+j];
        }
    }

    return SAD;
}

int intra_Planar(BYTE* ori, BYTE* ref, BYTE* pred, int* resi, BYTE* recon){
    int i,j;
    int SAD=0;

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            pred[i*BLOCK_SIZE+j]=(int)((BLOCK_SIZE-1-j)*ref[BLOCK_SIZE+i+1]+(j+1)*ref[BLOCK_SIZE]+(BLOCK_SIZE-1-i)*ref[i]+(i+1)*ref[2*BLOCK_SIZE+1]+BLOCK_SIZE+(log2(BLOCK_SIZE)+1)/2)>>((int)log2(BLOCK_SIZE)+1);
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            resi[i*BLOCK_SIZE+j]=ori[i*BLOCK_SIZE+j]-pred[i*BLOCK_SIZE+j];
            SAD+=abs(resi[i*BLOCK_SIZE+j]);
        }
    }

    for(i=0;i<BLOCK_SIZE;i++){
        for(j=0;j<BLOCK_SIZE;j++){
            recon[i*BLOCK_SIZE+j]=resi[i*BLOCK_SIZE+j]+pred[i*BLOCK_SIZE+j];
        }
    }

    return SAD;
    //return 999999;
}